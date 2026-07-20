"""
Integration tests for PyBasin using the two bundled example datasets.

These tests run a full single-scenario model for each example dataset and
verify that the key goodness-of-fit metrics are finite values in [0, 1].

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


def _run_single_scenario(dataset_name):
    """
    Load a dataset, disable file/figure output, and run one model scenario.

    Returns the tuple
        (T_gof, vr_gof, aft_age_gof, he_age_gof)
    with np.nan for metrics that are not simulated.
    """
    input_dir = os.path.join(PROJECT_ROOT, "input_data", dataset_name)

    ModelParameters, ParameterRanges = _load_params(dataset_name)
    params = ModelParameters()

    # Suppress all file and figure output so tests leave no side-effects.
    params.make_model_data_fig = False
    params.save_model_run_data = False

    # Validate that all required input CSVs are present.
    pybasin.check_input_data_files(input_dir, params)

    # Read all input data.
    (well_strats, strat_info_mod, salinity_bnd_df,
     T_data_df, vr_data_df,
     aft_samples, aft_ages,
     he_samples, he_data,
     salinity_data, surface_temp, litho_props) = pybasin.read_model_input_data(
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

    well_strat, well_strat_orig = pybasin.select_well_strat(well, well_strats)

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
     vr_gof, _vr_r2, _VR_model_data,
     aft_age_gof, _aft_age_error, _aft_age_r2, _AFT_data,
     he_age_gof, _he_age_error, _he_age_r2,
     _He_model_data,
     _model_results_series_updated) = result

    return T_gof, vr_gof, aft_age_gof, he_age_gof


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_example_dataset_1():
    """
    End-to-end test using example_dataset_1 (Roer Valley Graben, NDW-01).

    This dataset simulates temperature, vitrinite reflectance, and apatite
    fission track ages. All three goodness-of-fit values must be finite
    numbers in [0, 1].
    """
    T_gof, vr_gof, aft_age_gof, he_age_gof = _run_single_scenario(
        "example_dataset_1"
    )

    assert np.isfinite(T_gof), f"T_gof is not finite: {T_gof}"
    assert 0.0 <= T_gof <= 1.0, f"T_gof out of range [0, 1]: {T_gof}"

    assert np.isfinite(vr_gof), f"vr_gof is not finite: {vr_gof}"
    assert 0.0 <= vr_gof <= 1.0, f"vr_gof out of range [0, 1]: {vr_gof}"

    assert np.isfinite(aft_age_gof), f"aft_age_gof is not finite: {aft_age_gof}"
    assert 0.0 <= aft_age_gof <= 1.0, f"aft_age_gof out of range [0, 1]: {aft_age_gof}"


@pytest.mark.slow
def test_example_dataset_2():
    """
    End-to-end test using example_dataset_2 (Molasse Basin, E40).

    This dataset simulates apatite (U-Th)/He ages; there is no temperature
    data for the E40 outcrop, so T_gof is expected to be nan.  The He
    goodness-of-fit value must be a finite number in [0, 1].
    """
    T_gof, vr_gof, aft_age_gof, he_age_gof = _run_single_scenario(
        "example_dataset_2"
    )

    # No temperature measurements exist for the E40 outcrop, so T_gof is nan.
    assert np.isnan(T_gof) or (
        np.isfinite(T_gof) and 0.0 <= T_gof <= 1.0
    ), f"T_gof is neither nan nor in [0, 1]: {T_gof}"

    assert np.isfinite(he_age_gof), f"he_age_gof is not finite: {he_age_gof}"
    assert 0.0 <= he_age_gof <= 1.0, f"he_age_gof out of range [0, 1]: {he_age_gof}"
