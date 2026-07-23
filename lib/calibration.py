"""
Simple, flexible model calibration helpers for pybasin.

Provides a function to run one pybasin model scenario with one or more
input parameters overridden to chosen values, and a grid-search calibration
routine that scans a range of values for one parameter to find the
best-fitting value against a chosen goodness-of-fit metric.

This does not use pybasin's own multiple-model-run machinery
(ParameterRanges / the "_s" suffix convention in pybasin_params.py), which
runs a fixed, user-specified list of scenarios. Instead, run_single_scenario()
lets calibration code choose parameter values programmatically, eg. from a
grid search or an optimizer, one model run at a time.
"""

import importlib.util
import inspect
import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

import pybasin
from . import model_input_io


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_params(dataset_name):

    """
    load a fresh copy of the ModelParameters and ParameterRanges classes
    for one input dataset

    a fresh copy is loaded every time (rather than reusing a cached
    module), so that repeated calls, eg. from a calibration loop, do not
    leak parameter changes from one run into the next
    """

    params_path = os.path.join(
        PROJECT_ROOT, 'input_data', dataset_name, 'pybasin_params.py')
    spec = importlib.util.spec_from_file_location('pybasin_params', params_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module.ModelParameters, module.ParameterRanges


def run_single_scenario(dataset_name, param_overrides=None, well=None):

    """
    run a single pybasin model scenario for one example dataset, with one
    or more input parameters optionally overridden

    :param dataset_name: name of the input_data subdirectory, eg.
        "example_dataset_2"
    :param param_overrides: dict, optional. maps a ModelParameters
        attribute name to the value it should be set to for this run, eg.
        {"exhumed_thicknesses": np.array([2500.0])}. any attribute
        normally set in pybasin_params.py can be overridden this way, not
        just exhumed_thicknesses.
    :param well: str, optional. well to run; defaults to the first well
        in ModelParameters.wells

    :return: dict with keys T_gof, T_r2, vr_gof, vr_r2, aft_age_gof,
        aft_age_error, aft_age_r2, he_age_gof, he_age_error, he_age_r2
        (np.nan for metrics that are not simulated for this dataset), plus
        model_results_series, the full row of results (all input
        parameters and goodness-of-fit statistics) that would also be
        written to a model_results_*.csv file by a normal model run
    """

    input_dir = os.path.join(PROJECT_ROOT, 'input_data', dataset_name)

    ModelParameters, _ParameterRanges = load_params(dataset_name)
    params = ModelParameters()

    # no need for figures or saved .pck files during calibration
    params.make_model_data_fig = False
    params.save_model_run_data = False

    if param_overrides:
        for name, value in param_overrides.items():
            if not hasattr(params, name):
                msg = ('"%s" is not a parameter in the ModelParameters class '
                      'of %s/pybasin_params.py' % (name, dataset_name))
                raise AttributeError(msg)
            setattr(params, name, value)

    model_input_io.check_input_data_files(input_dir, params)

    (well_strats, strat_info_mod, salinity_bnd_df,
     T_data_df, vr_data_df,
     aft_samples, aft_ages,
     he_samples, he_data,
     salinity_data, surface_temp, litho_props,
     pressure_data, porosity_data) = model_input_io.read_model_input_data(
        input_dir, params)

    # deliberately ignore the dataset's own ParameterRanges "_s" scenario
    # list here: calibration takes over the role of choosing parameter
    # values, so param_overrides above must be the only source of
    # parameter values for this run
    scenario_param_names = [None]
    first_scenario_params = [None]

    if well is None:
        well = params.wells[0]
    well_number = params.wells.index(well)

    well_strat, well_strat_orig = model_input_io.select_well_strat(well, well_strats)

    surface_salinity_well = None
    if getattr(params, 'well_specific_surface_salinity_bnd', False) is True:
        surface_salinity_well = model_input_io.select_well_salinity_bnd(
            well, salinity_bnd_df)

    # pre-allocate a one-row results dataframe, mirrors what main() does
    attributes = inspect.getmembers(
        params, lambda a: not inspect.isroutine(a))
    param_names_list = [a[0] for a in attributes
                        if not (a[0].startswith('__') and a[0].endswith('__'))]
    columns = (['model_run', 'model_error', 'well', 'computational_time']
              + param_names_list
              + ['well', 'T_gof', 'vr_gof', 'aft_age_gof', 'aft_age_error',
                 'he_gof', 'he_error', 'mean_gof',
                 'resetting_depth_model_min', 'resetting_depth_model_max',
                 'resetting_depth_data_min', 'non-resetting_depth_data_max'])
    model_results_df = pd.DataFrame(index=[0], columns=columns)
    model_results_series = model_results_df.loc[0]

    output_dir = os.path.join(PROJECT_ROOT, 'model_output', dataset_name)
    csv_output_dir = os.path.join(output_dir, 'thermal_history_datafiles')

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
        show_progress=False,
    )

    (_, _, _,
     _model_run_data,
     _T_model_data, T_gof, T_r2,
     _C_data,
     vr_gof, vr_r2, _VR_model_data,
     aft_age_gof, aft_age_error, aft_age_r2, _AFT_data,
     he_age_gof, he_age_error, he_age_r2,
     _He_model_data,
     _pressure_gof, _pressure_r2, _Pressure_model_data,
     model_results_series_updated) = result

    return {
        'T_gof': T_gof, 'T_r2': T_r2,
        'vr_gof': vr_gof, 'vr_r2': vr_r2,
        'aft_age_gof': aft_age_gof, 'aft_age_error': aft_age_error,
        'aft_age_r2': aft_age_r2,
        'he_age_gof': he_age_gof, 'he_age_error': he_age_error,
        'he_age_r2': he_age_r2,
        'model_results_series': model_results_series_updated,
    }


def grid_search(dataset_name, param_name, values, metric='he_age_gof',
                maximize=True, well=None, fixed_overrides=None):

    """
    scan a range of values for one input parameter and record the
    resulting goodness of fit, to find the best-fitting value

    this is a simple, robust calibration approach: since goodness-of-fit
    metrics from a burial and thermal history model can be noisy or
    non-smooth functions of a parameter (eg. because of the discrete
    timestepping of the burial history), a grid search is more reliable
    than eg. a gradient-based optimizer, at the cost of needing to choose
    a value range and resolution in advance. see refine_grid_search()
    below to zoom in once a promising range has been located.

    :param dataset_name: name of the input_data subdirectory, eg.
        "example_dataset_2"
    :param param_name: name of the ModelParameters attribute to calibrate,
        eg. "exhumed_thicknesses". the values in `values` are assigned to
        this attribute directly, so if the parameter is an array in
        pybasin_params.py (eg. exhumed_thicknesses = np.array([3000.0])),
        pass full arrays here too, eg. [np.array([v]) for v in [...]]
    :param values: sequence of parameter values to try
    :param metric: str, optional. which entry of run_single_scenario()'s
        return dict to use as the score, default "he_age_gof". use an
        "..._error" or "..._rmse" style metric with maximize=False
    :param maximize: bool, optional. True (default) if a higher metric
        value is a better fit, appropriate for all "_gof" and "_r2"
        metrics. False if a lower value is a better fit, appropriate for
        "_error" or "_rmse" metrics
    :param well: str, optional, see run_single_scenario()
    :param fixed_overrides: dict, optional. additional parameter
        overrides (see run_single_scenario()) to hold fixed for every run

    Not every parameter value tried is guaranteed to correspond to a valid
    burial history (eg. an exhumation timing that is inconsistent with the
    stratigraphy raises a ValueError); such values are recorded with a NaN
    score rather than aborting the whole search, and a warning is logged
    for each one.

    :return: pandas.DataFrame with one row per value tried, columns
        "value" and the chosen metric (NaN for values that raised an
        error). df.attrs["best_value"] and df.attrs["best_<metric>"] hold
        the best-fitting value found (ignoring NaN rows) and its score.
    """

    rows = []
    for value in values:
        overrides = dict(fixed_overrides) if fixed_overrides else {}
        overrides[param_name] = value
        try:
            result = run_single_scenario(dataset_name, param_overrides=overrides,
                                         well=well)
            score = result[metric]
        except Exception as e:
            logger.warning('model run failed for %s = %s, recording NaN: %s'
                           % (param_name, value, e))
            score = np.nan
        rows.append({'value': value, metric: score})

    df = pd.DataFrame(rows)

    if not df[metric].notna().any():
        raise RuntimeError('every model run in this grid search failed, '
                           'see the logged warnings above for the errors')

    best_i = df[metric].idxmax() if maximize is True else df[metric].idxmin()
    df.attrs['best_value'] = df.loc[best_i, 'value']
    df.attrs['best_' + metric] = df.loc[best_i, metric]

    return df
