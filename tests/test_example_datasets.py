"""
Integration tests for PyBasin using the two bundled example datasets.

These tests run a full single-scenario model for each example dataset and
pin diagnostic error metrics to recorded reference values, so a silently
changed model result fails the test.

The vitrinite reflectance, apatite fission track and (U-Th)/He metrics are
mean absolute errors between observed and simulated values, which respond
continuously to a change in the model output. They replace the earlier
goodness-of-fit pins, which saturate (a goodness-of-fit value can sit at 0
or 1 across a wide range of simulated values) and so barely constrain the
model state. Temperature keeps its goodness-of-fit value because the data
are one-sided borehole temperature constraints; that makes the temperature
error collapse to 0 and carry no information here.

The reference values were recorded from the base-case scenario of each
bundled dataset. If an intentional change to the physics shifts them, run
the model once and update REFERENCE_METRICS below.

Run with:
    pytest tests/test_example_datasets.py

Skip slow tests:
    pytest tests/test_example_datasets.py -m "not slow"
"""

import importlib.util
import inspect
import os
import sys

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so that pybasin.py and lib/ are
# importable regardless of where pytest is invoked from.
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

# pybasin.py imports matplotlib.pyplot at module level; force the non-
# interactive Agg backend before that import so the tests work in headless
# environments (CI, SSH sessions, etc.).
import matplotlib  # noqa: E402
matplotlib.use("Agg")

import pybasin  # noqa: E402  (import after sys.path / matplotlib setup)
import lib.model_input_io as model_input_io  # noqa: E402


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _load_params(dataset_name):
    """Load the ModelParameters and ParameterRanges classes from a dataset."""
    params_path = os.path.join(
        PROJECT_ROOT, "input_data", dataset_name, "pybasin_params.py"
    )
    spec = importlib.util.spec_from_file_location("pybasin_params", params_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.ModelParameters, module.ParameterRanges


def _vr_mean_absolute_error(vr_model_data):
    """
    Mean absolute error between observed and simulated vitrinite reflectance.

    vr_model_data is the VR_model_data structure returned by the model; its
    last element is the vr_data_well dataframe holding the observed VR and
    the simulated_vr interpolated at each sample depth. Returns np.nan when
    VR is not simulated.
    """
    if vr_model_data is None:
        return np.nan
    vr_data_well = vr_model_data[-1]
    residual = vr_data_well["VR"] - vr_data_well["simulated_vr"]
    return float(np.nanmean(np.abs(residual)))


def _run_single_scenario(dataset_name):
    """
    Load a dataset, disable file/figure output, and run one model scenario.

    Returns a dict of diagnostic metrics
        {"T_gof", "vr_mae", "aft_age_error", "he_age_error"}
    with np.nan for metrics that the dataset does not simulate.

    The vitrinite reflectance, apatite fission track and (U-Th)/He metrics
    are mean absolute errors between observed and simulated values (in VR
    units and in My), which respond continuously to a change in the model
    output. The goodness-of-fit values these replace saturate (they can sit
    at 0 or 1 across a wide range of simulated values) and so are much less
    diagnostic. Temperature keeps its goodness-of-fit value because the data
    are one-sided borehole temperature constraints, which makes the
    temperature error collapse to 0 and carry no information here.
    """
    input_dir = os.path.join(PROJECT_ROOT, "input_data", dataset_name)

    ModelParameters, ParameterRanges = _load_params(dataset_name)
    params = ModelParameters()

    # Suppress all file and figure output so tests leave no side-effects.
    params.make_model_data_fig = False
    params.save_model_run_data = False

    # Validate that all required input CSVs are present.
    model_input_io.check_input_data_files(input_dir, params)

    # Read all input data.
    (well_strats, strat_info_mod, salinity_bnd_df,
     T_data_df, vr_data_df,
     aft_samples, aft_ages,
     he_samples, he_data,
     salinity_data, surface_temp, litho_props,
     pressure_data, porosity_data) = model_input_io.read_model_input_data(
        input_dir, params
    )

    # Build the list of model scenarios; only use the first (base-case) one.
    scenario_param_names, scenario_param_list = pybasin.setup_model_scenarios_new(
        ParameterRanges
    )
    if len(scenario_param_names) == 0:
        scenario_param_names = [None]
        scenario_param_list = [[None]]

    first_scenario_params = scenario_param_list[0]

    well = params.wells[0]
    well_number = 0

    well_strat, well_strat_orig = model_input_io.select_well_strat(well, well_strats)

    # Neither example dataset uses salinity simulation.
    surface_salinity_well = None

    # Pre-allocate a one-row results dataframe (mirrors what main() does).
    attributes = inspect.getmembers(
        params, lambda a: not inspect.isroutine(a)
    )
    param_names_list = [
        a[0] for a in attributes
        if not (a[0].startswith("__") and a[0].endswith("__"))
    ]
    columns = (
        ["model_run", "model_error", "well", "computational_time"]
        + param_names_list
        + [
            "well", "T_gof", "vr_gof", "aft_age_gof", "aft_age_error",
            "he_gof", "he_error", "mean_gof",
            "resetting_depth_model_min", "resetting_depth_model_max",
            "resetting_depth_data_min", "non-resetting_depth_data_max",
        ]
    )
    model_results_df = pd.DataFrame(index=[0], columns=columns)
    model_results_series = model_results_df.loc[0]

    # Output directories are required by the function signature but no files
    # will be written because save_model_run_data and save_burial_csv_files
    # are both False.
    output_dir = os.path.join(PROJECT_ROOT, "model_output", dataset_name)
    csv_output_dir = os.path.join(output_dir, "thermal_history_datafiles")

    result = pybasin.update_model_params_and_run_model_new(
        model_scenario_number=0,
        pybasin_params=params,
        param_names=scenario_param_names,
        param_set=first_scenario_params,
        well_number=well_number,
        well=well,
        model_results_series=model_results_series,
        well_strat=well_strat,
        well_strat_orig=well_strat_orig,
        strat_info_mod=strat_info_mod,
        surface_temp=surface_temp,
        surface_salinity_well=surface_salinity_well,
        litho_props=litho_props,
        T_data=T_data_df,
        vr_data_df=vr_data_df,
        aft_samples=aft_samples,
        aft_ages=aft_ages,
        he_samples=he_samples,
        he_data=he_data,
        salinity_data=salinity_data,
        csv_output_dir=csv_output_dir,
        output_dir=output_dir,
        log_screen_output=False,
        record_data=True,
        save_burial_csv_files=False,
    )

    (_, _, _,
     _model_run_data,
     _T_model_data, T_gof, _T_r2,
     _C_data,
     _vr_gof, _vr_r2, VR_model_data,
     _aft_age_gof, aft_age_error, _aft_age_r2, _AFT_data,
     _he_age_gof, he_age_error, _he_age_r2,
     _He_model_data,
     _pressure_gof, _pressure_r2, _Pressure_model_data,
     _model_results_series_updated) = result

    return {
        "T_gof": T_gof,
        "vr_mae": _vr_mean_absolute_error(VR_model_data),
        "aft_age_error": aft_age_error,
        "he_age_error": he_age_error,
    }


# ---------------------------------------------------------------------------
# Reference metrics recorded from the base-case scenario of each bundled
# dataset. The VR, AFT and He entries are mean absolute errors (see
# _run_single_scenario); T_gof is a goodness-of-fit value. np.nan marks a
# metric that the dataset does not simulate. Update these if an intentional
# physics change shifts the result.
# ---------------------------------------------------------------------------

REFERENCE_METRICS = {
    "example_dataset_1": {
        "T_gof": 1.0,
        "vr_mae": 0.09553769192011514,
        "aft_age_error": 103.35923645636582,
        "he_age_error": np.nan,
    },
    "example_dataset_2": {
        "T_gof": np.nan,
        "vr_mae": np.nan,
        "aft_age_error": np.nan,
        "he_age_error": 6.85405497468589,
    },
}


def _assert_metrics_match_reference(dataset_name, metrics):
    """Assert each metric matches its recorded reference (nan matches nan)."""
    reference = REFERENCE_METRICS[dataset_name]
    for name, expected in reference.items():
        got = metrics[name]
        if np.isnan(expected):
            assert np.isnan(got), (
                f"{dataset_name} {name}: expected nan, got {got!r}"
            )
        else:
            assert got == pytest.approx(expected, rel=1e-3), (
                f"{dataset_name} {name}: expected {expected!r}, got {got!r}"
            )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_example_dataset_1():
    """
    End-to-end test using example_dataset_1 (Roer Valley Graben, NDW-01).

    This dataset simulates temperature, vitrinite reflectance, and apatite
    fission track ages. The vitrinite reflectance and fission track mean
    absolute errors, and the temperature goodness-of-fit, must match their
    recorded reference values.
    """
    metrics = _run_single_scenario("example_dataset_1")
    _assert_metrics_match_reference("example_dataset_1", metrics)


@pytest.mark.slow
def test_example_dataset_2():
    """
    End-to-end test using example_dataset_2 (Molasse Basin, E40).

    This dataset simulates apatite (U-Th)/He ages; there is no temperature
    data for the E40 outcrop, so T_gof is nan. The He mean absolute age
    error must match its recorded reference value.
    """
    metrics = _run_single_scenario("example_dataset_2")
    _assert_metrics_match_reference("example_dataset_2", metrics)
