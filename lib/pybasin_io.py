"""
functions for saving and loading pybasin model run output data

pybasin used to save model run output (.pck files) as a single plain,
positional list, which had to be unpacked in the exact same order
everywhere it was read. this module replaces that with a named
dict of xarray.Dataset / dict structures, and provides a loader that
transparently upgrades files saved with the old list based format.
"""

import pickle

import numpy as np
import xarray as xr

OUTPUT_FORMAT_VERSION = 2


def build_grid_dataset(time_array_bp, surface_temp_array, basal_hf_array,
                       z_nodes, active_nodes, T_nodes, node_strat, node_age):

    """
    package the burial and thermal history grid into a labelled xarray Dataset

    dimensions are time (model timesteps) and node (model nodes, ie.
    discretized layers of the stratigraphic column)
    """

    ds = xr.Dataset(
        data_vars={
            'z': (['time', 'node'], z_nodes),
            'T': (['time', 'node'], T_nodes),
            'active': (['time', 'node'], active_nodes),
            'surface_temperature': (['time'], surface_temp_array),
            'basal_heat_flow': (['time'], basal_hf_array),
        },
        coords={
            'time': time_array_bp,
            'node_strat': ('node', np.asarray(node_strat)),
            'node_age': ('node', node_age),
        },
        attrs={
            'time_units': 'years before present',
            'depth_units': 'm',
            'temperature_units': 'degrees C',
            'heat_flow_units': 'W m-2',
        },
    )

    return ds


def _pack_T_data(T_data):

    if T_data is None:
        return None

    (T_depth, T_obs, T_obs_sigma, T_data_type, T_gof, T_rmse) = T_data

    return {
        'depth': T_depth,
        'observed': T_obs,
        'observed_sigma': T_obs_sigma,
        'data_type': T_data_type,
        'gof': T_gof,
        'rmse': T_rmse,
    }


def _pack_C_data(C_data):

    if C_data is None:
        return None

    (C_nodes, surface_salinity_array, salinity_lwr_bnd, salinity_depth,
     salinity_data, salinity_data_unc, salinity_RMSE, q_solute_bottom,
     q_solute_top) = C_data

    return {
        'C_nodes': C_nodes,
        'surface_salinity': surface_salinity_array,
        'salinity_lower_bnd': salinity_lwr_bnd,
        'salinity_depth': salinity_depth,
        'salinity_observed': salinity_data,
        'salinity_observed_unc': salinity_data_unc,
        'salinity_rmse': salinity_RMSE,
        'q_solute_bottom': q_solute_bottom,
        'q_solute_top': q_solute_top,
    }


def _pack_VR_data(VR_model_data):

    if VR_model_data is None:
        return None

    (vr_nodes, vr_depth, vr_obs, vr_min, vr_max, vr_obs_sigma, vr_GOF,
     vr_rmse, vr_data_well) = VR_model_data

    return {
        'vr_nodes': vr_nodes,
        'depth': vr_depth,
        'observed': vr_obs,
        'observed_min': vr_min,
        'observed_max': vr_max,
        'observed_sigma': vr_obs_sigma,
        'gof': vr_GOF,
        'rmse': vr_rmse,
        'data': vr_data_well,
    }


def _pack_simulated_AFT_data(simulated_AFT_data):

    if simulated_AFT_data is None:
        return None

    (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
     aft_ln_mean_nodes, aft_ln_std_nodes, aft_node_times_burial,
     aft_node_zs, aft_node_times, aft_node_temps) = simulated_AFT_data

    return {
        'age_nodes': aft_age_nodes,
        'age_nodes_min': aft_age_nodes_min,
        'age_nodes_max': aft_age_nodes_max,
        'ln_mean_nodes': aft_ln_mean_nodes,
        'ln_std_nodes': aft_ln_std_nodes,
        'node_times_burial': aft_node_times_burial,
        'node_zs': aft_node_zs,
        'node_times': aft_node_times,
        'node_temps': aft_node_temps,
    }


def _pack_AFT_data(AFT_data):

    if AFT_data is None:
        return None

    (simulated_AFT_data, aft_sample_names, aft_age_depth, aft_age,
     aft_age_stderr_min, aft_age_stderr_plus, aft_length_mean,
     aft_length_std, aft_age_samples, single_grain_aft_ages,
     single_grain_aft_ages_se_min, single_grain_aft_ages_se_plus,
     aft_age_bins, aft_age_pdfs, aft_age_GOF, aft_age_error,
     aft_sample_times, aft_sample_temps, time_array_bp,
     z_aft_samples, T_samples, aft_data_samples) = AFT_data

    return {
        'simulated': _pack_simulated_AFT_data(simulated_AFT_data),
        'sample_names': aft_sample_names,
        'age_depth': aft_age_depth,
        'age': aft_age,
        'age_stderr_min': aft_age_stderr_min,
        'age_stderr_plus': aft_age_stderr_plus,
        'length_mean': aft_length_mean,
        'length_std': aft_length_std,
        'age_samples': aft_age_samples,
        'single_grain_ages': single_grain_aft_ages,
        'single_grain_ages_se_min': single_grain_aft_ages_se_min,
        'single_grain_ages_se_plus': single_grain_aft_ages_se_plus,
        'age_bins': aft_age_bins,
        'age_pdfs': aft_age_pdfs,
        'gof': aft_age_GOF,
        'age_error': aft_age_error,
        'sample_times': aft_sample_times,
        'sample_temps': aft_sample_temps,
        'time_array_bp': time_array_bp,
        'z_samples': z_aft_samples,
        'T_samples': T_samples,
        'data': aft_data_samples,
    }


def _pack_simulated_AHe_data(simulated_AHe_data):

    if simulated_AHe_data is None:
        return None

    (ahe_age_nodes, ahe_age_nodes_min, ahe_age_nodes_max,
     ahe_node_times_burial, ahe_node_zs) = simulated_AHe_data

    return {
        'age_nodes': ahe_age_nodes,
        'age_nodes_min': ahe_age_nodes_min,
        'age_nodes_max': ahe_age_nodes_max,
        'node_times_burial': ahe_node_times_burial,
        'node_zs': ahe_node_zs,
    }


def _pack_AHe_data(AHe_data):

    if AHe_data is None:
        return None

    (ahe_sample_depths, ahe_ages_all_samples, ahe_ages_all_samples_SE,
     ahe_age_bin, ahe_age_pdfs, modeled_ahe_age_samples,
     modeled_ahe_age_samples_min, modeled_ahe_age_samples_max,
     ahe_age_gof, ahe_age_error, simulated_AHe_data,
     ahe_data_samples) = AHe_data

    return {
        'sample_depths': ahe_sample_depths,
        'ages_all_samples': ahe_ages_all_samples,
        'ages_all_samples_se': ahe_ages_all_samples_SE,
        'age_bin': ahe_age_bin,
        'age_pdfs': ahe_age_pdfs,
        'modeled_age_samples': modeled_ahe_age_samples,
        'modeled_age_samples_min': modeled_ahe_age_samples_min,
        'modeled_age_samples_max': modeled_ahe_age_samples_max,
        'gof': ahe_age_gof,
        'age_error': ahe_age_error,
        'simulated': _pack_simulated_AHe_data(simulated_AHe_data),
        'data': ahe_data_samples,
    }


def build_model_run_data(time_array_bp, surface_temp_array, basal_hf_array,
                         z_nodes, active_nodes, T_nodes, node_strat, node_age,
                         T_data, C_data, VR_model_data, AFT_data, AHe_data):

    """
    package a single model run's output data into the current, named
    dict / xarray based output format

    this replaces the old plain, positional list based format that
    pybasin used to save .pck files with. use load_model_run_data() or
    normalize_model_run_data() to read data saved in either format.
    """

    return {
        'format_version': OUTPUT_FORMAT_VERSION,
        'grid': build_grid_dataset(time_array_bp, surface_temp_array,
                                   basal_hf_array, z_nodes, active_nodes,
                                   T_nodes, node_strat, node_age),
        'T_data': _pack_T_data(T_data),
        'C_data': _pack_C_data(C_data),
        'VR_model_data': _pack_VR_data(VR_model_data),
        'AFT_data': _pack_AFT_data(AFT_data),
        'AHe_data': _pack_AHe_data(AHe_data),
    }


def _normalize_old_format(data):

    """
    convert data saved with the old, plain positional list based output
    format (pybasin versions before the named dict / xarray format) into
    the current named format
    """

    (time_array_bp, surface_temp_array, basal_hf_array,
     z_nodes, active_nodes, T_nodes, node_strat, node_age,
     T_data, C_data, VR_model_data, AFT_data, AHe_data) = data

    return build_model_run_data(
        time_array_bp, surface_temp_array, basal_hf_array,
        z_nodes, active_nodes, T_nodes, node_strat, node_age,
        T_data, C_data, VR_model_data, AFT_data, AHe_data)


def normalize_model_run_data(data):

    """
    return model run output data in the current named dict / xarray
    format, converting it first if it was saved with the old plain
    list based format
    """

    if isinstance(data, dict) and 'format_version' in data:
        return data

    if isinstance(data, (list, tuple)):
        return _normalize_old_format(data)

    msg = 'unrecognized model run output data format: %s' % type(data)
    raise ValueError(msg)


def load_model_run_data(path):

    """
    load a pybasin model run output (.pck) file, saved with either the
    current named dict / xarray format or the older plain list format,
    and return it in the current named format
    """

    with open(path, 'rb') as fin:
        data = pickle.load(fin)

    return normalize_model_run_data(data)
