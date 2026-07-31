"""
PyBasin, version 0.1

recoded in 2014-2015 to simplify the code and to use .csv files
+ pandas for model input

Elco Luijendijk, Goettingen University

<elco.luijendijk@geo.uni-goettingen.de>

"""

import sys
import os
import argparse
import importlib.util
import ast
import logging
# from runpy import run_path


import datetime
import pickle
import itertools
import time
import inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from multiprocessing import Pool

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s',
                    stream=sys.stdout)
logger = logging.getLogger(__name__)

import lib.pybasin_lib as pybasin_lib
import lib.pybasin_figures as pybasin_figures
import lib.pybasin_io as pybasin_io
import lib.model_data_comparison as model_data_comparison
import lib.model_input_io as model_input_io

# helium diffusion algortihm by Meesters and Dunai (2003)
try:
    import lib.helium_diffusion_models as he
except ImportError:
    logger.warning('failed to import native U-Th/He module')


# make sure multi-threading for numpy is turned off (this slows down the heat
# flow solution a lot...)
os.environ['OPENBLAS_NUM_THREADS'] = '1'

# set to True to use vectorized thermochron functions (for debugging/benchmarking)
vectorize_thermochron = True


def run_model_and_compare_to_data(well_number, well, well_strat,
                                  strat_info_mod, pybasin_params,
                                  surface_temp, surface_salinity_well,
                                  litho_props,
                                  csv_output_dir,
                                  output_dir,
                                  model_scenario_number,
                                  model_results_series,
                                  T_data, vr_data_df,
                                  aft_samples, aft_ages,
                                  he_samples, he_data, salinity_data,
                                  pressure_data=None,
                                  vr_method='easyRo',
                                  save_csv_files=True,
                                  show_progress=True):
    """
    run basin & thermal history model and compare modeled  and observed temperature, salinity,
    vitrinite reflectance, fission track ages and/or (U-Th)/He ages

    """

    # run burial history model
    model_result_vars = \
        pybasin_lib.run_burial_hist_model(well_number, well, well_strat,
                                          strat_info_mod, pybasin_params,
                                          surface_temp, surface_salinity_well,
                                          litho_props,
                                          csv_output_dir,
                                          model_scenario_number,
                                          save_csv_files=save_csv_files,
                                          show_progress=show_progress)

    simulate_fluid_flow = getattr(pybasin_params, 'simulate_fluid_flow', False)

    if simulate_fluid_flow is True:
        [P_ex_nodes, permeability_cells, q_fluid_top, q_fluid_bottom] = \
            model_result_vars[-4:]
        model_result_vars = model_result_vars[:-4]

    if pybasin_params.simulate_salinity is False:
        [geohist_df, time_array, time_array_bp,
         surface_temp_array, basal_hf_array,
         z_nodes, T_nodes, active_nodes,
         n_nodes, n_cells,
         node_strat, node_age,
         prov_start_nodes, prov_end_nodes, porosity_nodes, k_nodes,
         rho_nodes] = \
            model_result_vars
    else:
        [geohist_df, time_array, time_array_bp,
         surface_temp_array, basal_hf_array,
         z_nodes, T_nodes, active_nodes,
         n_nodes, n_cells,
         node_strat, node_age,
         prov_start_nodes, prov_end_nodes, porosity_nodes, k_nodes,
         rho_nodes,
         C_nodes, surface_salinity_array,
         salinity_lwr_bnd, Dw, q_solute_top, q_solute_bottom] = \
            model_result_vars

    # find out if exhumation end has changed
    exhumed_units = [unit[0] == '-' for unit in geohist_df.index]
    unit_names = list(geohist_df.index)
    exhumed_unit_blocks = []
    start = 0
    while True in exhumed_units[start:]:
        start = start + exhumed_units[start:].index(True)
        end = start + exhumed_units[start:].index(False)
        unit_block = [unit_names[start], unit_names[end-1]]
        exhumed_unit_blocks.append(unit_block)
        start = end

    # calculate cooling during exhumation
    for i, exhumed_unit_block in enumerate(exhumed_unit_blocks):
        start = geohist_df.loc[exhumed_unit_block[-1], 'age_bottom']
        end = geohist_df.loc[exhumed_unit_block[0], 'age_top']

        if end < 0:
            end = 0

        # find temperatures at start and end
        if start * 1e6 in time_array_bp and end * 1e6 in time_array_bp:
            start_ind = np.where(time_array_bp / 1e6 == start)[0][0]
            end_ind = np.where(time_array_bp / 1e6 == end)[0][0]

        else:
            logger.info('could not find exact start and end of exhumation in time array')
            logger.info('using closest time instead')
            start_ind = np.argmin(np.abs(time_array_bp / 1e6 - start))
            end_ind = np.argmin(np.abs(time_array_bp / 1e6 - end))

        # calculate cooling
        T_cooling = T_nodes[start_ind] - T_nodes[end_ind]

        # filter out eroded formations
        T_cooling_preserved = T_cooling[active_nodes[end_ind]]

        # store results
        model_results_series['start_exhumation_phase_%i' % i] = start
        model_results_series['end_exhumation_phase_%i' % i] = end
        model_results_series['mean_cooling_exhumation_phase_%i' % i] = T_cooling_preserved.mean()
        model_results_series['min_cooling_exhumation_phase_%i' % i] = T_cooling_preserved.min()
        model_results_series['max_cooling_exhumation_phase_%i' % i] = T_cooling_preserved.max()

    # record max temperature and depth
    model_results_series['max_temperature'] = T_nodes.max()
    model_results_series['max_present_temperature'] = \
        T_nodes[-1, active_nodes[-1]].max()
    model_results_series['max_depth'] = z_nodes.max()

    cebs_input = 'cebs.py'
    if pybasin_params.use_strat_map_input is True \
            and os.path.isfile(cebs_input) is True:

        logger.info('reading model input from cebs.py')

        # from . import cebs
        # model_results_series = cebs.present_temp_in_given_depth(
        #    z_nodes, model_scenario_number, T_nodes, model_results_series)
        # model_results_df = cebs.thermal_conductivity(
        #    model_results_series, node_strat, geohist_df, model_scenario_number,
        #    k_nodes)
        # model_results_df = cebs.porosity(model_results_df, node_strat,
        #                                 geohist_df, model_scenario_number,
        #                                 porosity_nodes)

        msg = "the cebs.py module has been deprecated, please contact the developer of the code"
        raise ValueError(msg)

    vr_nodes = None

    ################################
    # simulate vitrinite reflectance
    ################################
    if pybasin_params.simulate_VR is True:

        # find if there are VR samples for this well
        ind = ((vr_data_df['well'] == well)
               & (vr_data_df['depth'] < z_nodes[-1].max()))
        vr_data_well = vr_data_df[ind].copy()

        # interpolate vitrinite reflectance data
        if True in ind.values \
                or pybasin_params.calculate_thermochron_for_all_nodes is True:
            logger.info('calculating vitrinite reflectance for n=%i nodes' % n_nodes)

            vr_nodes = pybasin_lib.calculate_vr(T_nodes,
                                                active_nodes,
                                                time_array,
                                                n_nodes, vr_method=vr_method,
                                                vectorize_thermochron=vectorize_thermochron,
                                                show_progress=show_progress)

            # store surface and bottom VR value
            model_results_series['vr_surface'] = vr_nodes[-1, active_nodes[-1]][0]
            model_results_series['vr_bottom'] = vr_nodes[-1, active_nodes[-1]][-1]
            
            if pybasin_params.use_strat_map_input is True \
                    and os.path.isfile(cebs_input) is True:
                # model_results_series = cebs.vr_top_bot(
                #      model_results_series, node_strat,
                #     vr_nodes, geohist_df, model_scenario_number)
                # model_results_series = cebs.vr_middle(model_results_series, node_strat,
                #     vr_nodes, geohist_df,
                #     model_scenario_number, z_nodes)
                msg = "the strat_map_input option has been discontinued"
                raise ValueError(msg)

    if pybasin_params.simulate_salinity is True:
        # store depth to 1 g/L aand 0.035 kg/kg salinity

        # calculate density from salinity and temperature

        density = pybasin_lib.equations_of_state_batzle1992(
            np.zeros_like(T_nodes[-1, active_nodes[-1]]),
            T_nodes[-1, active_nodes[-1]],
            C_nodes[-1, active_nodes[-1]])

        C_nodes_gl = C_nodes[-1, active_nodes[-1]] * density
        depth_to_1g_per_l = np.interp([1.0], C_nodes_gl,
                                      z_nodes[-1, active_nodes[-1]])
        depth_to_seawater_salinity = np.interp([pybasin_params.salinity_seawater],
                                               C_nodes[-1, active_nodes[-1]],
                                               z_nodes[-1, active_nodes[-1]])

        model_results_series['depth_to_C=1gL-1'] = depth_to_1g_per_l
        model_results_series['depth_to_C=0.035kg/kg'] = depth_to_seawater_salinity

    if (pybasin_params.simulate_AFT is True
            or pybasin_params.simulate_He is True):
        calculate_thermochron_for_all_nodes = \
            pybasin_params.calculate_thermochron_for_all_nodes
        if pybasin_params.make_model_data_fig is True:
            calculate_thermochron_for_all_nodes = True

    ##############################################################
    # simulate apatite fission track ages and length distributions
    ##############################################################
    simulated_AFT_data = None
    location_has_AFT = False

    if pybasin_params.simulate_AFT is True:

        # find if there is any aft data for this well:
        ind = ((aft_samples['well'] == well) & (aft_samples['depth'] <= z_nodes[-1].max() + 1.0))
        aft_data_well = aft_samples[ind].copy()

        if True in ind.values:
            location_has_AFT = True
        else:
            logger.info('no AFT data found for this location')

        (modeled_aft_age_samples,
         modeled_aft_age_samples_min,
         modeled_aft_age_samples_max,
         aft_ln_mean_samples,
         aft_ln_std_samples,
         aft_sample_times_burial,
         aft_sample_zs,
         aft_sample_times,
         aft_sample_temps,
         simulated_AFT_data,
         z_aft_samples,
         T_samples) = pybasin_lib.assemble_data_and_simulate_aft(
            pybasin_params.resample_timesteps,
            pybasin_params.provenance_time_nt,
            n_nodes, time_array_bp,
            z_nodes, T_nodes, active_nodes,
            prov_start_nodes, prov_end_nodes,
            pybasin_params.annealing_kinetics_values,
            pybasin_params.annealing_kinetic_param,
            surface_temp,
            aft_data_well,
            calculate_thermochron_for_all_nodes=calculate_thermochron_for_all_nodes,
            annealing_eq=pybasin_params.annealing_equation,
            C0=pybasin_params.C0,
            C1=pybasin_params.C1,
            C2=pybasin_params.C2,
            C3=pybasin_params.C3,
            alpha=pybasin_params.alpha,
            location_has_AFT=location_has_AFT,
            provenance_start_temp=pybasin_params.provenance_start_temp,
            vectorize_thermochron=vectorize_thermochron,
            show_progress=show_progress)

        # store surface and bottom VR value
        if simulated_AFT_data is not None:

            (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
             aft_ln_mean_nodes, aft_ln_std_nodes,
             aft_node_times_burial, aft_node_zs,
             aft_node_times, aft_node_temps) = simulated_AFT_data

            model_results_series['aft_age_surface_min'] = aft_age_nodes[active_nodes[-1]][0].min()
            model_results_series['aft_age_surface_max'] = aft_age_nodes[active_nodes[-1]][0].max()
            model_results_series['aft_age_bottom_min'] = aft_age_nodes[active_nodes[-1]][-1].min()
            model_results_series['aft_age_bottom_max'] = aft_age_nodes[ active_nodes[-1]][-1].max()

            calculate_resetting_depth = True
            nodata_val = -99999.9

            # age cutoff for sample/node to be considered fully reset
            # todo: add this as a parameter to pybasin_params.py
            full_resetting_age = 0.5

            if aft_age_nodes_min.min() < full_resetting_age:
                ind_min = aft_age_nodes_min < full_resetting_age
                model_results_series['full_resetting_depth_aft_model_min'] = z_nodes[-1][ind_min][0]
                model_results_series['full_resetting_depth_aft_model_min'] = T_nodes[-1][ind_min][0]
            else:
                model_results_series['full_resetting_depth_aft_model_min'] = nodata_val
                model_results_series['full_resetting_depth_aft_model_min'] = nodata_val

            if aft_age_nodes_max.min() < full_resetting_age:
                ind_max = aft_age_nodes_max < full_resetting_age
                model_results_series['full_resetting_depth_aft_model_max'] = z_nodes[-1][ind_max][0]
                model_results_series['full_resetting_T_aft_model_min'] = T_nodes[-1][ind_max][0]
            else:
                model_results_series['full_resetting_depth_aft_model_max'] = nodata_val
                model_results_series['full_resetting_T_aft_model_min'] = nodata_val

            if calculate_resetting_depth is True and pybasin_params.calculate_thermochron_for_all_nodes is True:

                # modeled resetting depth
                ind_reset_min = aft_age_nodes_min <= node_age
                ind_reset_max = aft_age_nodes_max <= node_age
                if True in ind_reset_min:
                    model_results_series['partial_resetting_depth_aft_model_min'] = \
                        z_nodes[-1][ind_reset_min].min()
                else:
                    model_results_series['partial_resetting_depth_aft_model_min'] = nodata_val
                if True in ind_reset_max:
                    model_results_series['partial_resetting_depth_aft_model_max'] = \
                        z_nodes[-1][ind_reset_max].min()
                else:
                    model_results_series['partial_resetting_depth_aft_model_max'] = nodata_val

    #################################
    # simulate apatite (U-Th)/He ages
    #################################
    location_has_He_data = False

    if pybasin_params.simulate_He is True:

        resample_t = pybasin_params.resample_timesteps
        nt_prov = pybasin_params.provenance_time_nt

        # find if there is any U-Th/He data for this well:
        location_has_He_data = False
        if he_samples is not None:
            ind = he_samples['location'] == well
            he_samples_well = he_samples[ind].copy()

            if True in ind.values:
                location_has_He_data = True

                # check if any samples are at or near the bottom of the model
                model_bottom = z_nodes[-1].max()
                for _, he_row in he_samples_well.iterrows():
                    sample_depth = he_row['depth']
                    depth_diff = model_bottom - sample_depth
                    if abs(depth_diff) <= 1.0:
                        logger.warning(f'''Warning: He sample "{he_row['sample']}" in well "{well}" is located at or beyond the bottom of the model (sample depth: {sample_depth:.1f} m, model bottom: {model_bottom:.1f} m). This sample may be excluded from the model-data comparison. Consider extending the well stratigraphy slightly.''')
                    elif abs(depth_diff) <= 5.0:
                        logger.info(f'''Note: He sample "{he_row['sample']}" in well "{well}" is close to the bottom of the model (sample depth: {sample_depth:.1f} m, model bottom: {model_bottom:.1f} m). Verify that the sample falls within the modeled depth range. Small discrepancies may result from slight inaccuracies in decompaction / backstripping.''')

        decay_constant_238U = pybasin_params.decay_constant_238U
        decay_constant_235U = pybasin_params.decay_constant_235U
        decay_constant_232Th = pybasin_params.decay_constant_232Th

        (modeled_he_age_samples,
         modeled_he_age_samples_min,
         modeled_he_age_samples_max,
         he_node_times_burial,
         he_node_zs,
         simulated_He_data) = pybasin_lib.assemble_data_and_simulate_he(
            he_samples_well,
            he_data,
            decay_constant_238U,
            decay_constant_235U,
            decay_constant_232Th,
            n_nodes,
            resample_t,
            nt_prov,
            time_array_bp,
            z_nodes,
            T_nodes,
            active_nodes,
            prov_start_nodes,
            prov_end_nodes,
            surface_temp,
            calculate_thermochron_for_all_nodes=calculate_thermochron_for_all_nodes,
            C0=pybasin_params.C0,
            C1=pybasin_params.C1,
            C2=pybasin_params.C2,
            C3=pybasin_params.C3,
            alpha=pybasin_params.alpha,
            ahe_method=pybasin_params.ahe_method,
            provenance_start_temp=pybasin_params.provenance_start_temp,
            log_tT_paths=pybasin_params.log_tT_paths, tT_path_filename=pybasin_params.datafile_output_dir,
            vectorize_thermochron=vectorize_thermochron,
            show_progress=show_progress)

        # store surface and bottom VR value
        if simulated_He_data is not None:

            (he_age_nodes, he_age_nodes_min, he_age_nodes_max,
             he_node_times_burial, he_node_zs) = simulated_He_data

            he_age_nodes_array = np.array(he_age_nodes)
            he_age_nodes_min_array = np.min(np.array(he_age_nodes_min), axis=1)
            he_age_nodes_max_array = np.max(np.array(he_age_nodes_max), axis=1)

            model_results_series['he_age_surface_min'] = he_age_nodes_array[active_nodes[-1]][0].min()
            model_results_series['he_age_surface_max'] = he_age_nodes_array[active_nodes[-1]][0].max()
            model_results_series['he_age_bottom_min'] = he_age_nodes_array[active_nodes[-1]][-1].min()
            model_results_series['he_age_bottom_max'] = he_age_nodes_array[active_nodes[-1]][-1].max()

            calculate_resetting_depth = True
            nodata_val = -99999.9

            # age cutoff for sample/node to be considered fully reset
            # todo: add this as a parameter to pybasin_params.py
            full_resetting_age = 0.5

            if he_age_nodes_min_array.min() < full_resetting_age:
                ind_min = he_age_nodes_min_array < full_resetting_age
                model_results_series['full_resetting_depth_he_model_min'] = z_nodes[-1][ind_min][0]
                model_results_series['full_resetting_T_he_model_min'] = T_nodes[-1][ind_min][0]
            else:
                model_results_series['full_resetting_depth_he_model_min'] = nodata_val
                model_results_series['full_resetting_T_he_model_min'] = nodata_val

            if he_age_nodes_max_array.min() < full_resetting_age:
                ind_max = he_age_nodes_max_array < full_resetting_age
                model_results_series['full_resetting_depth_he_model_max'] = z_nodes[-1][ind_max][0]
                model_results_series['full_resetting_T_he_model_max'] = T_nodes[-1][ind_max][0]
            else:
                model_results_series['full_resetting_depth_he_model_max'] = nodata_val
                model_results_series['full_resetting_T_he_model_max'] = nodata_val

            if calculate_resetting_depth is True and pybasin_params.calculate_thermochron_for_all_nodes is True:

                # modeled resetting depth
                ind_reset_min = he_age_nodes_min_array <= node_age
                ind_reset_max = he_age_nodes_max_array <= node_age
                if True in ind_reset_min:
                    model_results_series['partial_resetting_depth_he_model_min'] = \
                        z_nodes[-1][ind_reset_min].min()
                else:
                    model_results_series['partial_resetting_depth_he_model_min'] = nodata_val
                if True in ind_reset_max:
                    model_results_series['partial_resetting_depth_he_model_max'] = \
                        z_nodes[-1][ind_reset_max].min()
                else:
                    model_results_series['partial_resetting_depth_he_model_max'] = nodata_val

    ##################################
    # calculate model goodness of fit:
    ##################################
    # calculate model error temperature data
    ind = (T_data['well'] == well) & (T_data['depth'] < z_nodes[-1].max())
    T_data_well = T_data[ind].copy()

    if True in ind.values:

        T_gof, T_rmse, T_r2 = model_data_comparison.model_data_comparison_T(T_data_well, z_nodes,
                                                       T_nodes, active_nodes)

        T_model_data = (T_data_well['depth'].values,
                        T_data_well['temperature'].values,
                        T_data_well['temperature_unc_1sigma'].values,
                        T_data_well['data_type'],
                        T_gof, T_rmse)

    else:
        T_rmse = np.nan
        T_gof = np.nan
        T_r2 = np.nan
        T_model_data = None

    # calculate model error VR data
    vr_rmse = np.nan
    vr_gof = np.nan
    vr_r2 = np.nan
    if pybasin_params.simulate_VR is True:

        # calculate model error vitrinite reflectance data
        ind = ((vr_data_df['well'] == well)
               & (vr_data_df['depth'] < z_nodes[-1].max()))
        vr_data_well = vr_data_df[ind].copy()

        # interpolate vitrinite reflectance data
        if True in ind.values:

            vr_rmse, vr_gof, vr_r2, vr_data_well = model_data_comparison.model_data_comparison_VR(
                vr_data_well,
                z_nodes, vr_nodes,
                active_nodes,
                vr_unc_sigma=pybasin_params.vr_unc_sigma)

    # calculate model error AFT data
    aft_age_gof = np.nan
    aft_age_error = np.nan
    aft_age_r2 = np.nan
    if pybasin_params.simulate_AFT is True:

        # calculate model error fission track data
        ind = ((aft_samples['well'] == well) & (aft_samples['depth'] <= z_nodes[-1].max() + 1.0))
        aft_data_well = aft_samples[ind].copy()

        if True in ind.values:

            (aft_age_gof, aft_age_error, aft_age_r2,
             single_grain_aft_ages,
             single_grain_aft_ages_se_min,
             single_grain_aft_ages_se_plus,
             age_bins, age_pdfs, aft_data_well) = \
                model_data_comparison.model_data_comparison_AFT_age(aft_data_well, aft_ages,
                                              modeled_aft_age_samples_min,
                                              modeled_aft_age_samples_max,
                                              two_sided_gof=pybasin_params.two_sided_gof,
                                              gof_age_percentile=pybasin_params.gof_age_percentile)

    # simulate apatite (U-Th)/He data
    he_age_gof = np.nan
    he_age_error = np.nan
    he_age_r2 = np.nan
    he_ages_all_samples = None
    he_ages_all_samples_SE = None
    he_age_pdfs_all_samples = None
    he_samples_well = None

    if pybasin_params.simulate_He is True:

        # calculate model error fission track data
        ind = ((he_samples['location'] == well) & (he_samples['depth'] <= z_nodes[-1].max() + 1.0))
        he_samples_well = he_samples[ind].copy()

        he_age_bin = np.linspace(0, prov_start_nodes.max(), 1000)

        if True in ind.values:

            (he_age_gof, he_age_error, he_age_r2, he_ages_all_samples,
             he_ages_all_samples_SE,
             he_age_bin, he_age_pdfs_all_samples, he_samples_well) = \
                model_data_comparison.model_data_comparison_he(he_samples_well, he_data,
                                          he_age_bin,
                                          modeled_he_age_samples_min,
                                          modeled_he_age_samples_max,
                                          two_sided_gof=pybasin_params.two_sided_gof,
                                          gof_age_percentile=pybasin_params.gof_age_percentile)

        else:
            logger.warning(f'Warning: simulate_He is True but no He samples found within the depth range of well "{well}". He data will be skipped.')

    # calculate model error salinity data
    salinity_rmse = np.nan
    salinity_gof = np.nan
    salinity_r2 = np.nan
    if pybasin_params.simulate_salinity is True:

        ind = (salinity_data['well'] == well) & \
              (salinity_data['depth'] < z_nodes[-1].max())
        salinity_data_well = salinity_data[ind]

        if True in ind.values:
            salinity_gof, salinity_rmse, salinity_r2 = model_data_comparison.model_data_comparison_salinity(
                salinity_data_well, z_nodes, C_nodes, active_nodes)

    # calculate model error pressure data
    pressure_rmse = np.nan
    pressure_gof = np.nan
    pressure_r2 = np.nan
    Pressure_model_data = None
    if simulate_fluid_flow is True and pressure_data is not None:

        ind = ((pressure_data['well'] == well)
               & (pressure_data['depth_bottom'] < z_nodes[-1].max()))
        pressure_data_well = pressure_data[ind].copy()

        if True in ind.values:
            pressure_unc_sigma = getattr(pybasin_params, 'pressure_unc_sigma', 1.0)
            pressure_gof, pressure_rmse, pressure_r2, pressure_data_well = \
                model_data_comparison.model_data_comparison_pressure(
                    pressure_data_well, z_nodes, P_ex_nodes, active_nodes,
                    pressure_unc_sigma=pressure_unc_sigma)

            Pressure_model_data = [pressure_data_well['obs_depth'].values,
                                   pressure_data_well['FSIP_MPa'].values,
                                   pressure_data_well['simulated_pressure'].values,
                                   pressure_gof, pressure_rmse, pressure_r2,
                                   pressure_data_well]
        else:
            logger.info('no pressure data found for this location')

    # assemble output data
    if pybasin_params.simulate_AFT is True and location_has_AFT is True:

        AFT_data = [simulated_AFT_data,
                    aft_data_well['sample'].values,
                    aft_data_well['depth'].values,
                    aft_data_well['aft_age'].values,
                    aft_data_well['aft_age_stderr_min'].values,
                    aft_data_well['aft_age_stderr_plus'].values,
                    aft_data_well['length_mean'].values,
                    aft_data_well['length_std'].values,
                    modeled_aft_age_samples,
                    single_grain_aft_ages,
                    single_grain_aft_ages_se_min,
                    single_grain_aft_ages_se_plus,
                    age_bins,
                    age_pdfs,
                    aft_age_gof,
                    aft_age_error,
                    aft_sample_times,
                    aft_sample_temps,
                    time_array_bp,
                    z_aft_samples, T_samples,
                    aft_data_well]

    else:
        logger.info('no AFT data found for this location')
        AFT_data = None

    if pybasin_params.simulate_VR is True:
        VR_data = [vr_nodes,
                   vr_data_well['depth'].values,
                   vr_data_well['VR'].values,
                   vr_data_well['VR_min'].values,
                   vr_data_well['VR_max'].values,
                   vr_data_well['VR_unc_1sigma'].values,
                   vr_gof,
                   vr_rmse,
                   vr_data_well]
    else:
        VR_data = None

    if pybasin_params.simulate_salinity is True:
        C_data = [C_nodes, surface_salinity_array,
                  salinity_lwr_bnd,
                  salinity_data_well['depth'].values,
                  salinity_data_well['salinity'].values,
                  salinity_data_well['salinity_unc_1sigma'].values,
                  salinity_rmse,
                  q_solute_bottom, q_solute_top]
    else:
        C_data = None

    if (pybasin_params.simulate_He is True
            and location_has_He_data is True):
        He_model_data = [he_samples_well['depth'].values,
                          he_ages_all_samples,
                          he_ages_all_samples_SE,
                          he_age_bin,
                          he_age_pdfs_all_samples,
                          modeled_he_age_samples,
                          modeled_he_age_samples_min,
                          modeled_he_age_samples_max,
                          he_age_gof, he_age_error,
                          simulated_He_data,
                          he_samples_well]

    else:
        He_model_data = None

    model_run_data = [time_array_bp,
                      surface_temp_array, basal_hf_array,
                      z_nodes, active_nodes, T_nodes,
                      node_strat, node_age, porosity_nodes, rho_nodes,
                      P_ex_nodes if simulate_fluid_flow is True else None]

    return (model_run_data,
            T_model_data, T_gof, T_r2,
            C_data, salinity_gof, salinity_rmse, salinity_r2,
            vr_gof, vr_rmse, vr_r2, VR_data,
            aft_age_gof, aft_age_error, aft_age_r2, AFT_data,
            he_age_gof, he_age_error, he_age_r2,
            He_model_data,
            pressure_gof, pressure_rmse, pressure_r2, Pressure_model_data,
            model_results_series)


def update_model_params_and_run_model_new(model_scenario_number,
                                          pybasin_params,
                                          param_names, param_set,
                                          well_number, well,
                                          model_results_series,
                                          well_strat, well_strat_orig,
                                          strat_info_mod,
                                          surface_temp,
                                          surface_salinity_well,
                                          litho_props,
                                          T_data, vr_data_df,
                                          aft_samples, aft_ages,
                                          he_samples, he_data,
                                          salinity_data,
                                          csv_output_dir,
                                          output_dir,
                                          log_screen_output,
                                          pressure_data=None,
                                          record_data=True,
                                          save_burial_csv_files=True,
                                          show_progress=True):

    """
    update the model parameters class and run the model one time

    new version, uses beo style parameter space / sensitivity analysis runs

    returns a tuple with the model result variables
    """

    if log_screen_output is True:
        log_output_dir = os.path.join(output_dir, 'log')
        if os.path.exists(log_output_dir) is False:
            os.makedirs(log_output_dir)

        log_fn = 'log_well_%s_model_scen_%i_PID_%s.out' % (well, model_scenario_number, str(os.getpid()))
        log_path = os.path.join(log_output_dir, log_fn)
        logger.info('redirecting screen output to log file %s' % log_path)

        err_fn = 'error_log_well_%s_model_scen_%i_PID_%s.out' % (well, model_scenario_number, str(os.getpid()))
        err_path = os.path.join(log_output_dir, err_fn)
        logger.info('redirecting screen output for model errors to error log file %s' % err_path)

        sys.stdout = open(log_path, "w")
        sys.stderr = open(err_path, "w")
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s',
                            stream=sys.stdout, force=True)

    well_strat = well_strat_orig.copy()

    # update default parameters in pybasin_params class
    for scenario_param_name, scenario_parameter in zip(param_names, param_set):

        if scenario_parameter is not None:
            # find model parameter name to adjust
            model_param_name = scenario_param_name[:-2]

            if hasattr(pybasin_params, model_param_name) is False:
                msg = 'the parameter %s is not in the ModelParameters class in the pybasin_params.py file' \
                      % model_param_name
                msg += ', even though it should be updated for model sensitivity or parameter exploration '
                msg += 'according to the ParameterRanges class. Please check if the spelling of the parameter'

                logger.info(msg)

                raise AttributeError(msg)

            logger.info('updating parameter %s from %s to %s' % (model_param_name, str(getattr(pybasin_params, model_param_name)), str(scenario_parameter)))

            # update model parameter
            setattr(pybasin_params, model_param_name, scenario_parameter)

    # check exhumation timing:
    n_exh_phases = len(pybasin_params.exhumation_period_starts)
    pybasin_params.exhumation_period_starts = np.array(pybasin_params.exhumation_period_starts)
    pybasin_params.exhumation_period_ends = np.array(pybasin_params.exhumation_period_ends)

    # set up array for end of exhumation, if not specified directly
    if hasattr(pybasin_params, "exhumation_durations"):
        logger.info('using exhumation duration and using this to calculate end of exhumation phase')
        pybasin_params.exhuamtion_durations = np.array(pybasin_params.exhumation_durations)
        pybasin_params.exhumation_period_ends = pybasin_params.exhumation_period_starts \
                                                - pybasin_params.exhumation_durations
        logger.info(f"calculated end of exhumation period:  {pybasin_params.exhumation_period_ends}")

    ind_nok = pybasin_params.exhumation_period_ends < 0.0
    pybasin_params.exhumation_period_ends[ind_nok] = 0.0

    # calculate end of exhumation
    # for exhumation_phase_id in range(n_exh_phases):
        # exhumation_duration_temp = pybasin_params.exhumation_period_starts[exhumation_phase_id] - \
        #                           pybasin_params.exhumation_period_ends[exhumation_phase_id]

        # exhumation_duration_temp = pybasin_params.exhumation_durations[exhumation_phase_id]

        # check if exhumation duration exceeds starting age of exhumation
        # if (exhumation_duration_temp >= (pybasin_params.exhumation_period_starts[exhumation_phase_id]
        #                                  - pybasin_params.max_hf_timestep) / 1e6):
        #     exhumation_duration_temp = (pybasin_params.exhumation_period_starts[exhumation_phase_id]
        #                                 - (pybasin_params.max_hf_timestep * 3) / 1e6)

        # pybasin_params.exhumation_period_ends[exhumation_phase_id] = \
        #     (pybasin_params.exhumation_period_starts[exhumation_phase_id] - exhumation_duration_temp)

        # print('adjusted duration of exhumation = %0.2f My' \
        #       % exhumation_duration_temp)

    # check for new parameters if they are present in the parameter file
    # and if not, set a default value
    # this is for parameters that have been added since the last release to ensure backwards compatibility
    if hasattr(pybasin_params, 'vr_method') is False:
        pybasin_params.vr_method = 'easyRo'

    if hasattr(pybasin_params, "log_tT_paths") is False:
        pybasin_params.log_tT_paths = False

    if hasattr(pybasin_params, "provenance_start_temp") is False:
        logger.info('no provenance_start_temp specified in input file, using a value of 120 degrees C')
        pybasin_params.provenance_start_temp = 120.0

    if hasattr(pybasin_params, "simulate_flushing") is False:
        pybasin_params.simulate_flushing = False

    if hasattr(pybasin_params, "flushing_clay_fraction") is False:
        pybasin_params.flushing_clay_fraction = 0.4

    if hasattr(pybasin_params, "flushing_max_depth") is False:
        pybasin_params.flushing_max_depth = 200.0

    # get values of all input parameters in pybasin_params class
    attributes = inspect.getmembers(
        pybasin_params, lambda attribute: not (inspect.isroutine(attribute)))
    attribute_dict = [attribute for attribute in attributes
                      if not (attribute[0].startswith('__') and
                              attribute[0].endswith('__'))]

    # store input parameter values in dataframe:
    for a in attribute_dict:
        if a[0] in model_results_series.index:
            if type(a[1]) is list or type(a[1]) is np.ndarray:
                model_results_series[a[0]] = str(a[1])
                # store entries of lists or arrays in separate columns
                # TODO: fix this, this does not play nice with the model_results_df in the main loop yet
                # if len(a[1]) < 10:
                #    for i in range(len(a[1])):
                #        col_title = a[0] + '_' + str(i)
                #        model_results_series[col_title] = a[1][i]
            else:
                model_results_series[a[0]] = a[1]

    # run model:
    (model_run_data,
     T_model_data, T_gof, T_r2,
     C_data, salinity_gof, salinity_rmse, salinity_r2,
     vr_gof, vr_rmse, vr_r2, VR_model_data,
     aft_age_gof, aft_age_error, aft_age_r2, AFT_data,
     he_age_gof, he_age_error, he_age_r2,
     He_model_data,
     pressure_gof, pressure_rmse, pressure_r2, Pressure_model_data,
     model_results_series) = \
        run_model_and_compare_to_data(well_number, well, well_strat,
                                      strat_info_mod, pybasin_params,
                                      surface_temp, surface_salinity_well,
                                      litho_props,
                                      csv_output_dir,
                                      output_dir,
                                      model_scenario_number,
                                      model_results_series,
                                      T_data, vr_data_df,
                                      aft_samples, aft_ages,
                                      he_samples, he_data,
                                      salinity_data,
                                      pressure_data=pressure_data,
                                      vr_method=pybasin_params.vr_method,
                                      save_csv_files=save_burial_csv_files,
                                      show_progress=show_progress)

    if record_data is True:
        # store gof in model results dataframe
        model_results_series['well'] = well
        model_results_series['T_gof'] = T_gof
        model_results_series['T_r2'] = T_r2
        if pybasin_params.simulate_VR is True:
            model_results_series['vr_gof'] = vr_gof
            model_results_series['vr_rmse'] = vr_rmse
            model_results_series['vr_r2'] = vr_r2
        if pybasin_params.simulate_AFT is True:
            model_results_series['aft_age_gof'] = aft_age_gof
            model_results_series['aft_age_error'] = aft_age_error
            model_results_series['aft_age_r2'] = aft_age_r2
        if pybasin_params.simulate_He is True:
            model_results_series['he_gof'] = he_age_gof
            model_results_series['he_error'] = he_age_error
            model_results_series['he_r2'] = he_age_r2
        if pybasin_params.simulate_salinity is True:
            model_results_series['salinity_gof'] = salinity_gof
            model_results_series['salinity_rmse'] = salinity_rmse
            model_results_series['salinity_r2'] = salinity_r2
        if getattr(pybasin_params, 'simulate_fluid_flow', False) is True:
            model_results_series['pressure_gof'] = pressure_gof
            model_results_series['pressure_rmse'] = pressure_rmse
            model_results_series['pressure_r2'] = pressure_r2

        # calculate cumulative salt loss due to diffusion
        if pybasin_params.simulate_salinity is True:
            (C_nodes, surface_salinity_array, salinity_lwr_bnd,
             salinity_well_depth,
             salinity_well,
             salinity_well_sigma,
             salinity_rmse,
             q_solute_bottom, q_solute_top) = C_data

            (time_array_bp,
             surface_temp_array, basal_hf_array,
             z_nodes, active_nodes, T_nodes,
             node_strat, node_age, porosity_nodes, rho_nodes,
             P_ex_nodes) = model_run_data

            duration = -np.diff(time_array_bp)
            duration = np.append(duration, time_array_bp[-1])
            duration *= 365.25 * 24 * 60 * 60
            total_solute_flux_top = q_solute_top * duration
            total_solute_flux_bottom = q_solute_bottom * duration

            model_results_series['cumalative_solute_flux_top_kg_m-2'] = np.sum(total_solute_flux_top)
            model_results_series['cumalative_solute_flux_bottom_kg_m-2'] = np.sum(total_solute_flux_bottom)

    # restore well strat file to original
    well_strat = well_strat_orig.copy()

    # screen output GOF data
    logger.info('')
    logger.info('temperature GOF = %0.2f, R2 = %0.2f' % (T_gof, T_r2))
    if pybasin_params.simulate_VR is True:
        logger.info('vitrinite reflectance GOF = %0.2f, R2 = %0.2f' % (vr_gof, vr_r2))
    if pybasin_params.simulate_AFT is True:
        logger.info('AFT age GOF = %0.2f, R2 = %0.2f' % (aft_age_gof, aft_age_r2))
        logger.info('AFT age error = %0.2f' % aft_age_error)
    if pybasin_params.simulate_He is True:
        logger.info('He GOF = %0.2f, R2 = %0.2f' % (he_age_gof, he_age_r2))
        logger.info('He age error = %0.2f' % he_age_error)
    if pybasin_params.simulate_salinity is True:
        logger.info('salinity GOF = %0.2f, R2 = %0.2f' % (salinity_gof, salinity_r2))
        logger.info('salinity RMSE = %0.4f' % salinity_rmse)
    if getattr(pybasin_params, 'simulate_fluid_flow', False) is True:
        logger.info('pressure GOF = %0.2f, R2 = %0.2f' % (pressure_gof, pressure_r2))
        logger.info('pressure RMSE = %0.4f' % pressure_rmse)
    logger.info('')

    return (well_number, well, model_scenario_number,
            model_run_data,
            T_model_data, T_gof, T_r2,
            C_data,
            vr_gof, vr_r2, VR_model_data,
            aft_age_gof, aft_age_error, aft_age_r2, AFT_data,
            he_age_gof, he_age_error, he_age_r2,
            He_model_data,
            pressure_gof, pressure_r2, Pressure_model_data,
            model_results_series)


def setup_model_scenarios_new(ParameterRanges):

    """
    Generate list of model parameters for multiple model runs

    """

    pr = ParameterRanges

    # create list with param values for each model run
    scenario_param_names_raw = dir(pr)
    scenario_param_names = [m for m in scenario_param_names_raw
                            if '__' not in m and '_s' in m]

    scenario_parameter_list = [getattr(pr, p)
                               for p in scenario_param_names]

    # construct list with all parameter combinations
    if pr.parameter_combinations is True:
        scenario_parameter_combinations = \
            list(itertools.product(*scenario_parameter_list))
    else:
        nscens = np.sum(np.array([len(sp) for sp in scenario_parameter_list
                                  if sp is not None]))
        nparams = len(scenario_parameter_list)
        scenario_parameter_combinations = []

        if pr.initial_base_run is True:
            scenario_parameter_combinations.append([None] * nparams)

        for j, sl in enumerate(scenario_parameter_list):
            if sl[0] is not None:
                sc = [None] * nparams
                for sli in sl:
                    sci = list(sc)
                    sci[j] = sli
                    scenario_parameter_combinations.append(sci)

    return scenario_param_names, scenario_parameter_combinations


def setup_model_output_df(n_scenarios):

    """
    Set up a Pandas dataframe to store model results
    """

    model_scenario_numbers = np.arange(n_scenarios)

    cols = ['well',
            'exhumation_magnitude', 'exhumation_start', 'exhumation_duration',
            'exhumation_segment_factor', 'exhumation_duration_factor',
            'basal_heat_flow',
            'T_gof', 'vr_gof', 'aft_age_gof', 'aft_age_error',
            'he_gof', 'he_error',
            'mean_gof',
            'objective_function',
            'resetting_depth_model_min',
            'resetting_depth_model_max',
            'resetting_depth_data_min',
            'non-resetting_depth_data_max']

    # set up dataframe to store model results
    model_results_df = pd.DataFrame(index=model_scenario_numbers,
                                    columns=cols)
    model_results_df2 = pd.DataFrame(columns=cols, index=[1])

    return model_results_df, model_results_df2


def main():

    parser = argparse.ArgumentParser(description='PyBasin: Model burial, exhumation and thermal history '
                                                 'and low-temperature thermochronology')

    parser.add_argument('model_input_subfolder', metavar='directory', default=None, nargs='?',
                        help='directory containing input dataset for PyBasin')

    parser.add_argument('-w', dest='wells',
                        help='specify wells to include, separated by a comma for multiple wells')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='show detailed per-timestep and per-unit debug output on screen')

    args = parser.parse_args()

    if args.verbose is True:
        logging.getLogger().setLevel(logging.DEBUG)

    # check if script dir in python path
    scriptdir = os.path.realpath(sys.path[0])
    if scriptdir not in sys.path:
        sys.path.append(scriptdir)

    if args.model_input_subfolder is not None:
        model_input_subfolder = os.path.join(scriptdir, args.model_input_subfolder)
    else:
        # read default input folder
        fin = open(os.path.join(scriptdir, 'default_input_folder.txt'))
        d = fin.readline()
        fin.close()

        scenario_name = d.split()[-1]
        model_input_subfolder = os.path.join(scriptdir, d.rstrip())

    logger.info('running model input data from folder %s' % model_input_subfolder)

    if os.path.isdir(model_input_subfolder) is False:
        msg = 'could not find the input directory %s' % model_input_subfolder
        msg += '\ncheck that the directory name is spelled correctly and that you are running pybasin.py '
        msg += 'from the pybasin source code directory'
        raise FileNotFoundError(msg)

    mpath = os.path.join(model_input_subfolder, 'pybasin_params.py')

    if os.path.isfile(mpath) is False:
        msg = 'could not find a pybasin_params.py file in the input directory %s' % model_input_subfolder
        msg += '\nevery pybasin input directory needs a pybasin_params.py file, see the example datasets '
        msg += 'in input_data/example_dataset_1 to input_data/example_dataset_4 for examples'
        raise FileNotFoundError(msg)

    # add the input folder to sys.path so that 'pybasin_params' can be
    # located by name, both here and by worker processes that unpickle
    # ModelParameters instances during parallel model runs (multiprocessing
    # spawn propagates sys.path to child processes, which need to be able
    # to re-import the module the class was defined in)
    if model_input_subfolder not in sys.path:
        sys.path.append(model_input_subfolder)

    spec = importlib.util.spec_from_file_location('pybasin_params', mpath)
    param_module = importlib.util.module_from_spec(spec)
    sys.modules['pybasin_params'] = param_module
    spec.loader.exec_module(param_module)

    Parameters_original = param_module.ModelParameters
    # model_scenarios = param_module.model_scenarios
    ParameterRanges = param_module.ParameterRanges

    Parameters = Parameters_original()

    year = 365.25 * 24.0 * 60 * 60.

    input_dir = model_input_subfolder
    output_dir = Parameters.output_dir
    datafile_output_dir = Parameters.datafile_output_dir
    csv_output_dir = datafile_output_dir

    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir)

    # pck_output_dir = os.path.join(output_dir, 'model_run_data_files')
    # the datafile output dir is used for the pickled model run data
    # (only if save_model_run_data is True), but also for the burial
    # history csv files and the tT path log files, which are saved
    # regardless of save_model_run_data. Always create it so those
    # writes do not fail with a missing directory error.
    if os.path.exists(datafile_output_dir) is False:
        logger.info('creating directory %s to store model result datafiles' % datafile_output_dir)
        os.makedirs(datafile_output_dir)

    fig_output_dir = output_dir
    if os.path.exists(fig_output_dir) is False:
        logger.info('creating directory %s to store model-data comparison figures' % fig_output_dir)
        os.makedirs(fig_output_dir)

    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month, today.year)

    model_input_io.check_input_data_files(input_dir, Parameters)

    (well_strats, strat_info_mod, salinity_bnd_df,
     T_data_df, vr_data_df,
     aft_samples, aft_ages,
     he_samples, he_data,
     salinity_data, surface_temp, litho_props,
     pressure_data, porosity_data) \
        = model_input_io.read_model_input_data(input_dir, Parameters)

    # model_scenario_param_list, params_to_change = setup_model_scenarios(model_scenarios,
    #                                                  pybasin_params.correct_exhumation_duration)

    model_scenario_param_names, model_scenario_param_list = setup_model_scenarios_new(ParameterRanges)

    if len(model_scenario_param_names) == 0:
        model_scenario_param_names = [None]
        model_scenario_param_list = [[None]]
        logger.info('single model run, setting up parameter set with base case values')

    # check sys arguments to run a particular well
    if args.wells is not None:
        wells_str = args.wells
        wells = wells_str.split(",")
    else:
        wells = Parameters.wells

    # allow wells = "all" or wells = ["all"] (in pybasin_params.py, or via -w)
    # to run every well found in well_stratigraphy.csv
    if isinstance(wells, str):
        wells = [wells]

    if len(wells) == 1 and wells[0].lower() == 'all':
        wells = well_strats['well'].unique().tolist()
        logger.info(f"wells set to 'all', found {len(wells)} wells in well_stratigraphy.csv")

    logger.info(f"running the following wells:  {wells}")

    n_scenarios = len(wells) * len(model_scenario_param_list)

    # get attributes
    params = inspect.getmembers(
        Parameters, lambda attribute: not (inspect.isroutine(attribute)))
    param_names = [attribute[0] for attribute in params
                       if not (attribute[0].startswith('__') and
                               attribute[0].endswith('__'))]

    # model_results_df, model_results_df2 = setup_model_output_df(n_scenarios)
    columns = ['model_run', 'model_error', 'well', 'computational_time'] + param_names
    columns += ['well',
                'T_gof', 'vr_gof', 'aft_age_gof', 'aft_age_error',
                'he_gof', 'he_error', 'mean_gof',
                'resetting_depth_model_min',
                'resetting_depth_model_max',
                'resetting_depth_data_min',
                'non-resetting_depth_data_max']

    # set up dataframe to store model results
    index = np.arange(n_scenarios)
    model_results_df = pd.DataFrame(index=index, columns=columns)

    model_scenario_number = 0

    if ParameterRanges.parallel_model_runs is True:

        # check the number of processors available on this machine, and cap
        # max_number_of_processes to that number if it is set too high in
        # pybasin_params.py / ParameterRanges
        n_cpu_available = os.cpu_count()
        n_processes = ParameterRanges.max_number_of_processes

        if n_cpu_available is not None and n_processes > n_cpu_available:
            logger.warning('max_number_of_processes (%i) is higher than the number '
                          'of processors available on this machine (%i), '
                          'limiting the number of simultaneous processes to %i'
                          % (n_processes, n_cpu_available, n_cpu_available))
            n_processes = n_cpu_available

        pool = Pool(processes=n_processes)
        logger.info('initialized parallel model runs with max %i simultaneous processes' % n_processes)

        processes = []
        done_processing = []

    start_time = time.time()

    #######################
    # go through all wells
    #######################
    for well_number, well in enumerate(wells):

        logger.info('=== well %s (%i/%i) ===' % (well, well_number + 1, len(wells)))

        if np.any(well_strats['well'] == well) == False:
            raise IOError("could not find well %s in well_stratigraphy.csv in directory %s" % (well, input_dir))

        well_strat, well_strat_orig = model_input_io.select_well_strat(well, well_strats)

        # read well specific surface salinity bnd condition
        if (Parameters.simulate_salinity is True
                and Parameters.well_specific_surface_salinity_bnd is True):

            surface_salinity_well = model_input_io.select_well_salinity_bnd(well, salinity_bnd_df)
            # else:
            #     print 'could not find well %s in surface salintiy bnd condition file' % well
            #     cols = ['age_start', 'age_end', 'surface_salinity']
            #     surface_salinity_well = salinity_bnd_df[cols]
        else:
            surface_salinity_well = None

        # estimate max number of nodes
        well_strat['thickness'] = well_strat['depth_bottom'] - well_strat['depth_top']
        well_strat['n_nodes_est'] = np.ceil(well_strat['thickness'] / Parameters.max_thickness)
        # buffer of 20000 m for number of nodes to be safe
        n_nodes_est = int(np.sum(well_strat['n_nodes_est']) + 2e4 / Parameters.max_thickness)

        dfc = pd.DataFrame(columns=['well'], index=np.arange(n_nodes_est))
        dfc['well'] = well

        if Parameters.simulate_salinity is True:
            n_ts_est = int(strat_info_mod['age_bottom'].max() * 1e6 / Parameters.max_hf_timestep * 3)
            dfqs = pd.DataFrame(columns=['well'], index=np.arange(n_ts_est))
            dfqs['well'] = well

        save_counter = 0

        # determine interval for saving output to csv file
        save_counter_interval = Parameters.csv_save_interval

        # go through model scenarios specified in the model_scenarios.py file:
        for well_scenario_no, model_scenario_params \
                in enumerate(model_scenario_param_list):

            logger.info('--- model scenario %i / %i ---' % (model_scenario_number + 1, len(model_scenario_param_list)))

            # restore original parameter values
            Parameters = Parameters_original()

            # estimate total runtime and time left
            if model_scenario_number / 250 == model_scenario_number / 250.0 and model_scenario_number > 0:

                now = time.time()
                time_passed = (now - start_time)
                time_per_scenario = time_passed / model_scenario_number
                time_left = \
                    (n_scenarios - model_scenario_number) * time_per_scenario

                tekst = 'model scenario %i / %i\n' % (model_scenario_number,
                                                      n_scenarios)
                tekst += 'time passed = %s\n' \
                         % datetime.timedelta(seconds=time_passed)
                tekst += 'time left = %s\n' \
                         % datetime.timedelta(seconds=time_left)
                tekst += 'time per scenario = %s\n' \
                         % datetime.timedelta(seconds=time_per_scenario)

                logger.info(tekst)

                logger.info('writing estimated runtime to runtime.txt')

                fout = open('runtime_%s.txt' % well, 'w')
                fout.write(tekst)
                fout.close()

            # restore original well strat dataframe
            well_strat = well_strat_orig.copy()

            model_results_series = model_results_df.loc[model_scenario_number].copy()

            if ParameterRanges.parallel_model_runs is False:

                # run a single model scenario

                log_screen_output = False

                (well_number_check, well_check, model_scenario_number_check,
                 model_run_data,
                 T_model_data, T_gof, T_r2,
                 C_data,
                 vr_gof, vr_r2, VR_model_data,
                 aft_age_gof, aft_age_error, aft_age_r2, AFT_data,
                 he_age_gof, he_age_error, he_age_r2,
                 He_model_data,
                 pressure_gof, pressure_r2, Pressure_model_data,
                 model_results_series_updated) = update_model_params_and_run_model_new(
                    model_scenario_number,
                    Parameters,
                    model_scenario_param_names, model_scenario_params,
                    well_number, well,
                    model_results_series,
                    well_strat, well_strat_orig,
                    strat_info_mod,
                    surface_temp,
                    surface_salinity_well,
                    litho_props,
                    T_data_df, vr_data_df,
                    aft_samples, aft_ages,
                    he_samples, he_data,
                    salinity_data,
                    csv_output_dir,
                    output_dir,
                    log_screen_output,
                    pressure_data=pressure_data,
                    save_burial_csv_files=Parameters.save_model_run_data)

                well_number_store, well_store, model_scenario_number_store = well_number, well, model_scenario_number

                for col in model_results_series_updated.index:
                    if col not in model_results_df.columns:
                        model_results_df[col] = np.nan

                model_results_df.loc[model_scenario_number_store] = model_results_series_updated

                keep_on_processing = True

            else:

                # set up a new parallel model run
                log_screen_output = True

                p = pool.apply_async(update_model_params_and_run_model_new,
                                     (model_scenario_number,
                                      Parameters,
                                      model_scenario_param_names, model_scenario_params,
                                      well_number, well,
                                      model_results_series,
                                      well_strat, well_strat_orig,
                                      strat_info_mod,
                                      surface_temp,
                                      surface_salinity_well,
                                      litho_props,
                                      T_data_df, vr_data_df,
                                      aft_samples, aft_ages,
                                      he_samples, he_data,
                                      salinity_data,
                                      csv_output_dir,
                                      output_dir,
                                      log_screen_output),
                                     {'pressure_data': pressure_data,
                                      'save_burial_csv_files': Parameters.save_model_run_data,
                                      'show_progress': False})

                processes.append(p)

                done_processing.append(False)

                keep_on_processing = True

            while keep_on_processing is True:
                # if ParameterRanges.parallel_model_runs is False or process_parallel_runs_result is True:

                if ParameterRanges.parallel_model_runs is False:
                    model_result_ready = True

                else:
                    model_result_ready = False

                    # check if there is already a result to process:
                    for ip, p in enumerate(processes):
                        # print('checking if any runs are done')

                        # TODO: this is not so elegant, replace with conditional
                        #  loop instead and wrap output in a seperate function
                        if model_result_ready is False and p.ready() is True and done_processing[ip] is False:

                            logger.info('process %i is done' % ip)

                            p_result = p.get()

                            (well_number_store, well_store, model_scenario_number_store,
                             model_run_data,
                             T_model_data, T_gof, T_r2,
                             C_data,
                             vr_gof, vr_r2, VR_model_data,
                             aft_age_gof, aft_age_error, aft_age_r2, AFT_data,
                             he_age_gof, he_age_error, he_age_r2,
                             He_model_data,
                             pressure_gof, pressure_r2, Pressure_model_data,
                             model_results_series_updated) = p_result

                            for col in model_results_series_updated.index:
                                if col not in model_results_df.columns:
                                    model_results_df[col] = np.nan

                            model_results_df.loc[model_scenario_number_store] = model_results_series_updated

                            model_result_ready = True
                            done_processing[ip] = True

                if model_result_ready is False:
                    pass
                    
                else:
                    # process results of a single model run

                    today = datetime.datetime.now()
                    today_str = '%i-%i-%i' % (today.day, today.month, today.year)

                    # save salinity and T data
                    (time_array_bp,
                     surface_temp_array, basal_hf_array,
                     z_nodes, active_nodes, T_nodes,
                     node_strat, node_age,
                     porosity_nodes, rho_nodes, P_ex_nodes) = model_run_data

                    model_run_data_fig = pybasin_io.build_model_run_data(
                        time_array_bp, surface_temp_array, basal_hf_array,
                        z_nodes, active_nodes, T_nodes, node_strat, node_age,
                        T_model_data, C_data, VR_model_data, AFT_data, He_model_data,
                        porosity_nodes=porosity_nodes, rho_nodes=rho_nodes,
                        P_ex_nodes=P_ex_nodes,
                        pressure_data=Pressure_model_data, porosity_data=porosity_data)

                    # l = len(z_nodes[-1, active_nodes[-1]]) - 1
                    # dfc.loc[:l, 'depth_s%i' % model_scenario_number] = \
                    #    z_nodes[-1, active_nodes[-1]]
                    # dfc['depth_s%i' % model_scenario_number] = z_nodes[-1, active_nodes[-1]]

                    # dfc.loc[:l, 'T_s%i' % model_scenario_number] = \
                    #    T_nodes[-1, active_nodes[-1]]

                    if Parameters.save_model_run_data is True:

                        fn = os.path.join(output_dir,
                                      'model_data_%s_%s_ms%i.pck'
                                      % (well_store, today_str,
                                         model_scenario_number_store))
                        logger.info('saving all data for model run as %s' % fn)
                        fout = open(fn, 'wb')
                        pickle.dump(model_run_data_fig, fout)
                        fout.close()

                    ##############################
                    # save model results .csv file
                    ##############################
                    if save_counter == 0 or save_counter >= save_counter_interval:

                        if wells[0] == wells[-1]:
                            well_txt = wells[0]
                        else:
                            well_txt = '%s-%s' % (wells[0], wells[-1])
                        fn = os.path.join(output_dir, 'model_results_%s_%s_ms0-%i.csv'
                                          % (today_str, well_txt,
                                             n_scenarios))
                        logger.info('saving model results .csv file %s' % fn)
                        model_results_df.to_csv(fn, index_label='model_scenario_number')

                    ####################################
                    # detailed model output to csv files
                    ####################################
                    if Parameters.save_model_run_data is True:

                        # save modeled T-t paths, all nodes
                        _, n_nodes_store = T_nodes.shape

                        rs = Parameters.resample_timesteps
                        n_steps = len(time_array_bp[::rs])

                        cols = ['time_bp']

                        df_tt2 = pd.DataFrame(columns=cols,
                                              index=np.arange(n_steps))

                        df_tt2['time_bp'] = time_array_bp[::rs]

                        # new, store new cols in dicts first to avoid fragmented
                        # warning by pandas
                        z_cols = {}
                        T_cols = {}

                        for i in range(n_nodes_store):

                            # add depth
                            col_name = 'z_node_%i' % i
                            z_col = z_nodes[::rs, i]
                            z_col[active_nodes[::rs, i] == False] = np.nan
                            #df_tt2[col_name] = z_col
                            z_cols[col_name] = z_col

                            # add temperature
                            col_name = 'T_node_%i' % i
                            T_col = T_nodes[::rs, i]
                            T_col[active_nodes[::rs, i] == False] = np.nan
                            #df_tt2[col_name] = T_col
                            T_cols[col_name] = T_col
                            
                        # join the columns to the dataframe in one go to avoid fragmented warning by pandas
                        new_cols = pd.DataFrame({**z_cols, **T_cols}, index=df_tt2.index)
                        df_tt2 = pd.concat([df_tt2, new_cols], axis=1)

                        fn = os.path.join(csv_output_dir, 'time_depth_temp_%s_%s_ms%i.csv'
                                          % (well_store, today_str,
                                             model_scenario_number_store))

                        logger.info('saving time-temperature paths to %s' % fn)
                        df_tt2.to_csv(fn, index_label='timestep')

                        # salinity data:
                        if Parameters.simulate_salinity is True:
                            (C_nodes, surface_salinity_array, salinity_lwr_bnd,
                             salinity_well_depth,
                             salinity_well,
                             salinity_well_sigma,
                             salinity_rmse,
                             q_solute_bottom, q_solute_top) = C_data

                            l = len(z_nodes[-1, active_nodes[-1]]) - 1

                            dfc.loc[:l, 'salinity_s%i' % model_scenario_number_store] = \
                                C_nodes[-1, active_nodes[-1]]

                            # save solute flux data
                            columns = ['solute_flux_top_kg_m-2_s-1',
                                      'solute_flux_bottom_kg_m-2_s-1']

                            lqs = len(time_array_bp) - 1

                            dfqs.loc[:lqs, 'solute_flux_top_kg_m-2_s-1_s%i'
                                     % model_scenario_number_store] = q_solute_top
                            dfqs.loc[:lqs, 'solute_flux_bottom_kg_m-2_s-1_s%i'
                                     % model_scenario_number_store] = q_solute_bottom
                            dfqs.loc[:lqs, 'time_yr_s%i' % model_scenario_number_store] = \
                                time_array_bp
                            fn = os.path.join(csv_output_dir,
                                              'solute_flux_data_%s_%s_ms%i.csv'
                                              % (well_store, today_str, n_scenarios))
                            logger.info('saving solute flux data to %s' % fn)
                            dfqs.to_csv(fn, index=False)

                        ## VR and temperature
                        if Parameters.simulate_VR is True:

                            [vr_nodes,
                             vr_depth,
                             vr_obs,
                             vr_min,
                             vr_max,
                             vr_obs_sigma,
                             vr_GOF,
                             vr_rmse,
                             vr_data_well] = VR_model_data

                            if vr_nodes is not None:

                                nn = len(z_nodes[-1, active_nodes[-1]])
                                ind_nn = dfc.index[:nn]

                                dfc.loc[ind_nn, 'depth_s%i' % model_scenario_number_store] = \
                                    z_nodes[-1, active_nodes[-1]]
                                dfc.loc[ind_nn, 'T_s%i' % model_scenario_number_store] = T_nodes[-1, active_nodes[-1]]
                                dfc.loc[ind_nn, 'VR_s%i' % model_scenario_number_store] = vr_nodes[-1, active_nodes[-1]]

                                # save depth vs T and VR data
                                fn = os.path.join(csv_output_dir, 'modeled_depth_T_and_VR_%s_%s_ms%i.csv'
                                                  % (well_store, today_str, model_scenario_number_store))
                                logger.info('saving depth, temperature and VR data to %s' % fn)
                                dfc.to_csv(fn, index=False)

                                # save depth vs T and VR data
                                fn = os.path.join(csv_output_dir, 'model_data_comparison_VR_%s_%s_ms%i.csv'
                                                  % (well_store, today_str, model_scenario_number_store))
                                logger.info('saving depth, temperature and VR data to %s' % fn)
                                vr_data_well.to_csv(fn, index=False)

                        # AFT data:
                        if Parameters.simulate_AFT is True and AFT_data is not None:

                            [simulated_AFT_data,
                             aft_sample,
                             aft_age_depth,
                             aft_age,
                             aft_age_stderr_min,
                             aft_age_stderr_plus,
                             aft_length_mean,
                             aft_length_std,
                             aft_age_samples,
                             single_grain_aft_ages,
                             single_grain_aft_ages_se_min,
                             single_grain_aft_ages_se_plus,
                             aft_age_bins,
                             aft_age_pdfs,
                             aft_age_GOF,
                             aft_age_error,
                             aft_sample_times,
                             aft_sample_temps,
                             time_array_bp,
                             z_aft_samples, T_samples, aft_data_well] = AFT_data

                            # find max number of timesteps for AFT history
                            max_steps = np.max([len(tii) for ti in aft_sample_times for tii in ti])

                            n_models = len(aft_sample_times[0])

                            n_samples = len(aft_sample_times)

                            cols = []
                            for sample_i in range(n_samples):
                                for model_j in range(n_models):
                                    cols = cols + ['name_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['depth_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['aft_age_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['time_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['temp_sample_%i_model_%i' % (sample_i, model_j)]

                            df_tt_aft = pd.DataFrame(columns=cols, index=np.arange(max_steps))

                            for sample_i in range(n_samples):
                                df_tt_aft.loc[0, 'name_sample_%i_model_%i' % (sample_i, model_j)] = \
                                    aft_sample[sample_i]
                                df_tt_aft.loc[0, 'depth_sample_%i_model_%i' % (sample_i, model_j)] = \
                                    aft_age_depth[sample_i]
                                df_tt_aft.loc[0, 'aft_age_sample_%i_model_%i' % (sample_i, model_j)] = \
                                    aft_age[sample_i]

                                for model_j in range(n_models):
                                    steps = len(aft_sample_times[sample_i][model_j])

                                    df_tt_aft.loc[:(steps-1), 'time_sample_%i_model_%i' % (sample_i, model_j)] = \
                                        aft_sample_times[sample_i][model_j]
                                    df_tt_aft.loc[:(steps-1), 'temp_sample_%i_model_%i' % (sample_i, model_j)] = \
                                        aft_sample_temps[sample_i][model_j]

                            # save AFT time-temperature paths
                            fn = os.path.join(csv_output_dir,
                                              'aft_sample_time_temp_%s_%s_ms%i.csv'
                                              % (well_store, today_str,
                                                 model_scenario_number_store))
                            logger.info('saving time-temperature paths AFT samples to %s' % fn)
                            df_tt_aft.to_csv(fn, index=False)

                            ###################################################
                            # save csv file model-data comparison for AFT ages:
                            ###################################################
                            _, n_prov, n_kin = aft_age_samples.shape
                            for sample_i in range(n_samples):
                                for prov_model_i in range(n_prov):
                                    for kin_i in range(n_kin):
                                        aft_data_well.loc[sample_i,
                                                          'modeled_age_prov_%i_kinetic_param_%i'
                                                          % (prov_model_i, kin_i)] = \
                                            aft_age_samples[sample_i, prov_model_i, kin_i]

                            fn = os.path.join(csv_output_dir,
                                              'aft_model_vs_data_%s_%s_ms%i.csv'
                                              % (well_store, today_str,
                                                 model_scenario_number_store))
                            logger.info('saving modeled AFT data for samples to %s' % fn)
                            # df_aft.to_csv(fn)
                            aft_data_well.to_csv(fn)

                        #  U-Th/He data
                        if Parameters.simulate_He is True and He_model_data is not None:

                            [he_sample_depths,
                             he_ages_all_samples,
                             he_ages_all_samples_SE,
                             he_age_bin,
                             he_age_pdfs,
                             modeled_he_age_samples,
                             modeled_he_age_samples_min,
                             modeled_he_age_samples_max,
                             he_age_gof, he_age_error,
                             simulated_He_data,
                             he_data_samples] = He_model_data

                            n_samples_he = len(modeled_he_age_samples)

                            for sample_i, age_mod_sample in enumerate(modeled_he_age_samples):
                                for grain_i, age_mod_grains in enumerate(age_mod_sample):
                                    for prov_model_i, age_mod_prov in enumerate(age_mod_grains):
                                        he_data_samples.loc[sample_i, 'modeled_he_age_grain_%i_prov_%i'
                                                          % (grain_i, prov_model_i)] = age_mod_prov
                                fn = os.path.join(csv_output_dir,
                                                  'he_model_vs_data_%s_%s_ms%i.csv'
                                                  % (well_store, today_str,
                                                     model_scenario_number_store))
                                logger.info('saving modeled  U-Th/He data to %s' % fn)
                                he_data_samples.to_csv(fn)

                    #############################
                    # make a model vs data figure
                    #############################
                    if Parameters.make_model_data_fig is True:
                        fig = pybasin_figures.model_vs_data_figure(
                            model_run_data_fig,
                            contour_variable=Parameters.contour_variable,
                            show_strat_column=Parameters.show_strat_column,
                            show_thermochron_data=Parameters.show_thermochron_data,
                            show_gof_stats=getattr(Parameters, 'show_gof_stats', True))
                    #    vr_data['depth'], vr_data['VR'], vr_data['unc_range_sigma'])

                        # fn = os.path.join(fig_output_dir,
                        #                   'model_data_fig_%s_%s_ms%i.%s'
                        #                   % (well, today_str,
                        #                      model_scenario_number,
                        #                      pybasin_params.fig_adj))
                        # print 'saving model-data comparison figure %s' % fn
                        # fig.savefig(fn, dpi=200)
                        # pl.clf()

                        if type(Parameters.fig_adj) is list:
                            for fa in Parameters.fig_adj:
                                fn = os.path.join(fig_output_dir,
                                              'model_data_fig_%s_%s_ms%i.%s'
                                              % (well_store, today_str,
                                                 model_scenario_number_store,
                                                 fa))
                                logger.info('saving model-data comparison figure %s' % fn)
                                fig.savefig(fn, dpi=200)

                            pl.clf()
                        else:
                            fn = os.path.join(fig_output_dir,
                                              'model_data_fig_%s_%s_ms%i.%s'
                                              % (well_store, today_str,
                                                 model_scenario_number_store,
                                                 Parameters.fig_adj))
                            logger.info('saving model-data comparison figure %s' % fn)
                            fig.savefig(fn, dpi=200)
                            pl.clf()

                    save_counter += 1

                    if save_counter >= save_counter_interval:
                        save_counter = 0

                if ParameterRanges.parallel_model_runs is True:
                    # n_open_processes = len(processes) - np.sum(np.array([p.ready() for p in processes]))
                    n_open_processes = len(processes) - np.sum(done_processing)
                    last_model_scenario = model_scenario_number >= (n_scenarios - 1)
                    if n_open_processes >= n_processes or last_model_scenario is True:
                        keep_on_processing = True
                    else:
                        keep_on_processing = False

                    if n_open_processes == 0:
                        keep_on_processing = False
                else:
                    keep_on_processing = False

            model_scenario_number += 1
        logger.info('done with all model scenarios')

        # creating separate columns for exhumation rates and heat flow
        cols_to_separate = ['exhumed_thicknesses']
        for col in cols_to_separate:
            # new_df = pd.DataFrame(model_results_df[col].values.tolist(),
            #                      index=model_results_df.index)
            # _, n_new_cols = new_df.shape
            try:
                new_cols = [ast.literal_eval(','.join(c.split())) for c in model_results_df[col].values.tolist()]
                n_new_cols = len(new_cols[0])
                new_cols_a = np.array(new_cols)

                new_col_names = [col + '_' + str(i) for i in range(n_new_cols)]
                for new_col_name, new_col_value in zip(new_col_names, new_cols_a.T):
                    model_results_df[new_col_name] = new_col_value
            except:
                logger.info('failed to separate column %s in model output dataframe / csv file' % col)

        # saving model results to a .csv file
        if wells[0] == wells[-1]:
            well_txt = wells[0]
        else:
            well_txt = '%s-%s' % (wells[0], wells[-1])
        fn = os.path.join(output_dir, 'model_results_%s_%s_ms0-%i_final.csv'
                          % (today_str, well_txt,
                             n_scenarios))
        logger.info('saving model results .csv file %s' % fn)
        model_results_df.to_csv(fn, index_label='model_scenario_number')

    logger.info('done')


if __name__ == '__main__':
    main()
