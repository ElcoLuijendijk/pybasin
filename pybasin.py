"""
PyBasin, version 0.1

recoded in 2014-2015 to simplify the code and to use .csv files
+ pandas for model input

Elco Luijendijk, Goettingen University

<elco.luijendijk@geo.uni-goettingen.de>

"""

import sys
import os
import pdb
import datetime
import pickle
import itertools
import time
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as pl

import lib.pybasin_lib as pybasin_lib
import lib.pybasin_figures as pybasin_figures

# check if script dir in python path
scriptdir = os.path.realpath(sys.path[0])
if scriptdir not in sys.path:
    sys.path.append(scriptdir)

print ''

if len(sys.argv) > 1 and sys.argv[-1][-3:] != '.py':
    scenario_name = sys.argv[-1]
    model_input_subfolder = os.path.join(scriptdir, 'input_data',
                                         scenario_name)
    print 'running model input data from folder %s' % model_input_subfolder

else:
    # read default input folder
    fin = open(os.path.join(scriptdir, 'default_input_folder.txt'))
    d = fin.readline()
    fin.close()
    scenario_name = d.split()[0]
    model_input_subfolder = os.path.join(scriptdir, 'input_data',
                                         scenario_name)
    print 'running model input data from folder %s' % model_input_subfolder

# import model parameter and model functions scripts
sys.path.append(model_input_subfolder)
import pybasin_params
import model_scenarios

year = 365.25 * 24.0 * 60 * 60.

input_dir = pybasin_params.input_dir
output_dir = pybasin_params.output_dir

datafile_output_dir = pybasin_params.datafile_output_dir

if os.path.exists(output_dir) is False:
    os.mkdir(output_dir)

#pck_output_dir = os.path.join(output_dir, 'model_run_data_files')
if pybasin_params.save_model_run_data is True and os.path.exists(datafile_output_dir) is False:
    print 'creating directory %s to store model result datafiles' \
        % datafile_output_dir
    os.mkdir(datafile_output_dir)

csv_output_dir = os.path.join(output_dir, 'burial_history_csv_files')
if os.path.exists(csv_output_dir) is False:
    print 'creating directory %s to store burial history .csv datafiles' \
        % csv_output_dir
    os.mkdir(csv_output_dir)

fig_output_dir = os.path.join(output_dir, 'model_data_fig')
if os.path.exists(fig_output_dir) is False:
    print 'creating directory %s to store model-data comparison figures' \
        % fig_output_dir
    os.mkdir(fig_output_dir)


# read all input data
print 'reading input data'
# stratigraphy description
strat_info = pd.read_csv(os.path.join(input_dir, 'stratigraphy_info.csv'))
strat_info = strat_info.set_index('strat_unit')

# well stratigraphy
well_strats = pd.read_csv(os.path.join(input_dir, 'well_stratigraphy.csv'))

# surface temperature history
Ts = pd.read_csv(os.path.join(input_dir, 'surface_temperature.csv'))

# lithology properties
litho_props = pd.read_csv(os.path.join(input_dir, 'lithology_properties.csv'))
litho_props = litho_props.set_index('lithology')

# present-day temperature data
T_data = pd.read_csv(os.path.join(input_dir, 'temperature_data.csv'))

# vitrinite reflectance data
vr_data = pd.read_csv(os.path.join(input_dir, 'vitrinite_reflectance.csv'))

## fission track data
# fission track sample info
#aft_samples = pd.read_csv(os.path.join(input_dir, 'aft_samples.csv'))

# fission track data
aft_data_raw = pd.read_csv(os.path.join(input_dir, 'aft_data.csv'))

# filter out samples with low grain count
ind = ((aft_data_raw['n_grains'] > pybasin_params.min_grain_no)
        | (aft_data_raw['n_grains'].isnull()))
aft_data = aft_data_raw[ind]


########
# calculate porosity-depth and thermal parameters for each strat unit
# find lithology columns in stratigraphy dataframe
cols_temp = strat_info.columns[2:].tolist()
prov_cols = [col for col in cols_temp if 'provenance' in col]
litho_cols = [col for col in cols_temp if 'provenance' not in col]

litho_cols.sort()
litho_props = litho_props.sort_index()

# check if lithology data is given for each lithology
assert litho_props.index[:-1].tolist() == litho_cols

# create new copy of dataframe to store results
strat_info_mod = strat_info.copy()

# add new lithology properties columns to stratigraphy dataframe
for litho_prop_name in litho_props.columns:
    strat_info_mod[litho_prop_name] = 0

# go through all litho properties
for litho_prop_name in litho_props.columns:

    # go through all lithology types present in strat unit
    for col in litho_cols:
        # find fraction of lithology and multiply by lithology property
        s = strat_info[col].astype(float) * litho_props[litho_prop_name][col]

        # add fraction to column in strat dataframe
        strat_info_mod[litho_prop_name] = strat_info_mod[litho_prop_name] + s

# construct list with exhumation and heatflow estimates
model_scenario_param_list = \
    list(itertools.product(model_scenarios.exhumation_magnitudes,
                           model_scenarios.exhumation_starts_and_durations,
                           model_scenarios.basal_heat_flow_scenarios))

n_scenarios = len(model_scenarios.wells) * len(model_scenario_param_list)

model_scenario_numbers = np.arange(n_scenarios)

cols = ['well',
        'exhumation_magnitude', 'exhumation_start', 'exhumation_duration',
        'basal_heat_flow',
        'T_gof', 'vr_gof', 'aft_age_gof',
        'mean_gof',
        'resetting_depth_model_min',
        'resetting_depth_model_max',
        'resetting_depth_data_min',
        'non-resetting_depth_data_max']

# set up dataframe to store model results
model_results = pd.DataFrame(index=model_scenario_numbers,
                             columns=cols)

model_scenario_number = 0

start_time = time.time()

#######################
# go through all wells
#######################
for well_number, well in enumerate(model_scenarios.wells):

    print 'x' * 20
    print 'well %s, %i/%i' % (well, well_number + 1,
                              len(model_scenarios.wells))

    if np.any(well_strats['well'] == well) is False:
        #print 'warning, well %s not in well strat file'
        raise IOError('error, could not find well %s in well strat file')

    well_strat = well_strats[well_strats['well'] == well]
    well_strat = well_strat.set_index('strat_unit')

    # copy original well strat file
    well_strat_orig = well_strat.copy()

    # go through all model scenarios:
    for well_scenario_no, model_scenario_params \
            in enumerate(model_scenario_param_list):

        print '-' * 10
        print 'model scenario %i / %i' % (model_scenario_number,
                                          len(model_scenario_param_list))

        # estimate total runtime and time left
        if (model_scenario_number / 1000 == model_scenario_number / 1000.0
                and model_scenario_number > 0):

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

            print tekst

            print 'writing estimated runtime to runtime.txt'

            fout = open('runtime.txt', 'w')
            fout.write(tekst)
            fout.close()

        # update model scenario parameters:
        (exhumation_magnitude, exhumation_start_and_duration,
         basal_heat_flow) = model_scenario_params

        exhumation_start, exhumation_duration = exhumation_start_and_duration

        # record scenario params in dataframe:
        model_results.ix[model_scenario_number, 'exhumation_magnitude'] = \
            exhumation_magnitude
        model_results.ix[model_scenario_number, 'exhumation_start'] = \
            exhumation_start
        model_results.ix[model_scenario_number, 'exhumation_duration'] = \
            exhumation_duration
        model_results.ix[model_scenario_number, 'basal_heat_flow'] = \
            basal_heat_flow

        # update parameter file with new value of exhumation
        exhumation_phase_id = \
            pybasin_params.exhumation_phase_ids.index(
                model_scenarios.exhumation_scenarios_period)
        if exhumation_magnitude is not None:
            pybasin_params.exhumed_thicknesses[exhumation_phase_id] = \
                exhumation_magnitude
            print 'exhumation = %0.1f m' % exhumation_magnitude
        if exhumation_start is not None:
            pybasin_params.exhumation_period_starts[exhumation_phase_id] = \
                exhumation_start
            print 'start of exhumation = %0.2f Ma' % exhumation_start
        if exhumation_duration is not None:
            pybasin_params.exhumation_period_ends[exhumation_phase_id] = \
                (pybasin_params.exhumation_period_starts[exhumation_phase_id]
                 - exhumation_duration)
            print 'duration of exhumation = %0.2f My' % exhumation_duration

        # adjust entire heat flow history
        if basal_heat_flow is not None:
            if model_scenarios.basal_heat_flow_scenario_period == 'all':
                pybasin_params.heatflow_history[:] = basal_heat_flow

                print 'basal heat flow = %0.2e' % basal_heat_flow

        # restore original well strat dataframe
        well_strat = well_strat_orig.copy()

        # run burial history model
        (geohist_df, time_array, time_array_bp,
         surface_temp_array, basal_hf_array,
         z_nodes, T_nodes, active_nodes,
         n_nodes, n_cells,
         node_strat, node_age,
         prov_start_nodes, prov_end_nodes) = \
            pybasin_lib.run_burial_hist_model(well_number, well, well_strat,
                                              strat_info_mod, pybasin_params,
                                              Ts, litho_props,
                                              csv_output_dir,
                                              model_scenario_number)

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
            start = geohist_df.ix[exhumed_unit_block[-1], 'age_bottom']
            end = geohist_df.ix[exhumed_unit_block[0], 'age_top']

            if end < 0:
                pdb.set_trace()

            # find temperatures at start and end
            if start * 1e6 in time_array_bp and end * 1e6 in time_array_bp:
                start_ind = np.where(time_array_bp / 1e6 == start)[0][0]
                end_ind = np.where(time_array_bp / 1e6 == end)[0][0]

            else:
                print 'warning, could not find exact start and end of ' \
                      'exhumation in time array'
                print 'using closest time instead'
                start_ind = np.argmin(np.abs(time_array_bp / 1e6 - start))
                end_ind = np.argmin(np.abs(time_array_bp / 1e6 - end))

            # calculate cooling
            T_cooling = T_nodes[start_ind] - T_nodes[end_ind]

            # filter out eroded formations
            T_cooling_preserved = T_cooling[active_nodes[end_ind]]

            # store results
            model_results.ix[model_scenario_number,
                             'start_exhumation_phase_%i' % i] = start
            model_results.ix[model_scenario_number,
                             'end_exhumation_phase_%i' % i] = end
            model_results.ix[model_scenario_number,
                             'mean_cooling_exhumation_phase_%i' % i] = \
                T_cooling_preserved.mean()
            model_results.ix[model_scenario_number,
                             'min_cooling_exhumation_phase_%i' % i] = \
                T_cooling_preserved.min()
            model_results.ix[model_scenario_number,
                             'max_cooling_exhumation_phase_%i' % i] = \
                T_cooling_preserved.max()

        # record max temperature and depth
        model_results.ix[model_scenario_number,
                         'max_temperature'] = T_nodes.max()
        model_results.ix[model_scenario_number,
                         'max_present_temperature'] = \
            T_nodes[-1, active_nodes[-1]].max()
        model_results.ix[model_scenario_number,
                         'max_depth'] = z_nodes.max()

        #def model_data_comparison(time_array, time_array_bp, z_nodes, T_nodes,
        # active_nodes, n_nodes, pybasin_params):
        # model VR for all formations
        print 'calculating vitrinite reflectance for n=%i nodes' % n_nodes

        vr_nodes = pybasin_lib.calculate_vr(T_nodes,
                                            active_nodes,
                                            time_array,
                                            n_nodes)

        resample_t = pybasin_params.resample_AFT_timesteps
        nt_prov = pybasin_params.provenance_time_nt

        (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
         aft_ln_mean_nodes, aft_ln_std_nodes,
         aft_node_times_burial, aft_node_zs) =\
            pybasin_lib.simulate_aft(
                resample_t, nt_prov, n_nodes, time_array_bp,
                z_nodes, T_nodes, active_nodes,
                prov_start_nodes, prov_end_nodes,
                pybasin_params.annealing_kinetics_values,
                pybasin_params.annealing_kinetic_param,
                Ts)

        # store surface VR value
        model_results.ix[model_scenario_number, 'vr_surface'] = \
            vr_nodes[-1, active_nodes[-1]][0]

        ##################################
        # calculate model goodness of fit:
        ##################################
        # calculate model error temperature data
        ind = (T_data['well'] == well) & (T_data['depth'] < z_nodes[-1].max())
        T_data_well = T_data[ind]
        if True in ind.values:
            T_data_well['simulated_T'] = np.interp(T_data_well['depth'],
                                                   z_nodes[-1, active_nodes[-1]],
                                                   T_nodes[-1, active_nodes[-1]])

            # calculate model error temperature data
            #ind = [T_data_well['temperature_unc_1sigma'].isnull()]
            #T_data_well['temperature_unc_1sigma'][ind] = pybasin_params.vr_unc_sigma
            T_data_well['residual'] = (T_data_well['temperature']
                                       - T_data_well['simulated_T'])
            T_data_well['residual_norm'] = (T_data_well['residual']
                                            / T_data_well['temperature_unc_1sigma'])
            T_data_well['P_fit'] = \
                (1.0
                 - scipy.stats.norm.cdf(np.abs(T_data_well['residual_norm']))) * 2

            # assign P=0 for temperatures lower than uncorrected BHTs
            ind_bht_nofit = ((T_data_well['residual'] > 0)
                             & (T_data_well['data_type'] == 'BHT'))
            ind_bht_ok = ((T_data_well['residual'] <= 0)
                          & (T_data_well['data_type'] == 'BHT'))
            T_data_well['P_fit'][ind_bht_nofit] = 0
            T_data_well['residual'][ind_bht_nofit] = 15.0
            T_data_well['P_fit'][ind_bht_ok] = 1.00
            T_data_well['residual'][ind_bht_ok] = 0.0

            T_rmse = np.sqrt(np.mean(T_data_well['residual']**2))
            T_gof = np.mean(T_data_well['P_fit'])
        else:
            T_rmse = np.nan
            T_gof = np.nan

        # interpolate vitrinite reflectance data
        ind = ((vr_data['well'] == well)
               & (vr_data['depth'] < z_nodes[-1].max()))
        vr_data_well = vr_data[ind]
        if True in ind.values:
            vr_data_well['simulated_vr'] = \
                np.interp(vr_data_well['depth'],
                          z_nodes[-1, active_nodes[-1]],
                          vr_nodes[-1, active_nodes[-1]])

            # calculate model error vitrinite data
            #vr_data_well['VR_unc_1sigma'] = vr_data_well['VR_std']
            ind = vr_data_well['VR_unc_1sigma'].isnull()
            vr_data_well['VR_unc_1sigma'][ind] = pybasin_params.vr_unc_sigma
            vr_data_well['residual'] = (vr_data_well['VR']
                                        - vr_data_well['simulated_vr'])
            vr_data_well['residual_norm'] = (vr_data_well['residual']
                                             / vr_data_well['VR_unc_1sigma'])
            vr_data_well['P_fit'] = \
                (1.0
                 - scipy.stats.norm.cdf(np.abs(vr_data_well['residual_norm'])))*2

            vr_rmse = np.sqrt(np.mean(vr_data_well['residual']**2))
            vr_gof = np.mean(vr_data_well['P_fit'])
        else:
            vr_rmse = np.nan
            vr_gof = np.nan

        # calculate model error fission track data
        ind = ((aft_data['well'] == well)
               & (aft_data['depth'] <= z_nodes[-1].max() + 1.0))
        aft_data_well = aft_data[ind]
        if True in ind.values:

            aft_data_well['simulated_AFT_min'] = \
                np.interp(aft_data_well['depth'],
                          z_nodes[-1, active_nodes[-1]],
                          aft_age_nodes_min[active_nodes[-1]])

            aft_data_well['simulated_AFT_max'] = \
                np.interp(aft_data_well['depth'],
                          z_nodes[-1, active_nodes[-1]],
                          aft_age_nodes_max[active_nodes[-1]])

            # get pdf of ages
            age_bins, age_pdfs = \
                pybasin_lib.calculate_aft_ages_pdf(
                    aft_data_well['aft_age'],
                    aft_data_well['95ci_minus'] / 1.96,
                    aft_data_well['95ci_plus'] / 1.96)
            aft_data_well['age_peak_debug'] = \
                age_bins[np.argmax(age_pdfs, axis=1)]

            # go through samples and find out how much of age pdf is covered by
            #  min and max simulated age
            for i, sample_ix in enumerate(aft_data_well.index):

                # TODO: find more elegant solution for 0.0 simulated AFT age
                # and check if GOF for AFT ages of 0.0 Ma are correct
                if aft_data_well.ix[sample_ix, 'simulated_AFT_min'] == 0:
                    start_ind = 0
                else:
                    start_ind = np.where(
                        aft_data_well.ix[sample_ix, 'simulated_AFT_min']
                        >= age_bins)[0][-1]

                if aft_data_well.ix[sample_ix, 'simulated_AFT_max'] == 0.0:
                    end_ind = 0
                else:
                    # np.where(0.0 >= age_bins)[0]
                    end_ind = np.where(
                        aft_data_well.ix[sample_ix, 'simulated_AFT_max']
                        <= age_bins)[0][0]

                pdf_fit_sum = np.sum(age_pdfs[i][start_ind:end_ind])
                pdf_nofit_sum = np.sum(age_pdfs[i][:start_ind]) \
                    + np.sum(age_pdfs[i][end_ind:])
                aft_data_well.ix[sample_ix, 'GOF_aft_ages'] = pdf_fit_sum

            # make filters for central age and population data
            ca_ind = aft_data_well['data_type'] == 'central_age'
            pop_ind = aft_data_well['data_type'] == 'population_age'

            aft_age_gof = aft_data_well['GOF_aft_ages'][ca_ind].mean()
        else:
            aft_age_gof = np.nan

        # calculate mean gof temperature, vitrinite and aft age data
        gofs = np.array([T_gof, vr_gof, aft_age_gof], dtype=float)
        ind = np.isnan(gofs) == False
        gof_weights = np.array(pybasin_params.gof_weights)[ind]
        # make sure weights add up to 1
        gof_weights = gof_weights / gof_weights.sum()
        gof_mean = np.sum(gofs[ind] * gof_weights)

        # screen output GOF data
        print 'temperature GOF = %0.2f' % T_gof
        print 'vitrinite reflectance GOF = %0.2f' % vr_gof
        print 'AFT age GOF = %0.2f' % aft_age_gof
        print 'weighted mean GOF = %0.2f' % gof_mean

        # store gof in model results dataframe
        model_results.ix[model_scenario_number, 'well'] = well
        model_results.ix[model_scenario_number, 'T_gof'] = T_gof
        model_results.ix[model_scenario_number, 'vr_gof'] = vr_gof
        model_results.ix[model_scenario_number, 'aft_age_gof'] = aft_age_gof
        model_results.ix[model_scenario_number, 'mean_gof'] = gof_mean

        ############################
        # calculate resetting depth
        ############################

        # modeled resetting depth
        ind_reset_min = aft_age_nodes_min <= node_age
        ind_reset_max = aft_age_nodes_max <= node_age
        if True in ind_reset_min:
            model_results.ix[model_scenario_number,
                             'resetting_depth_model_min'] = \
                z_nodes[-1][ind_reset_min].min()
        if True in ind_reset_max:
            model_results.ix[model_scenario_number,
                             'resetting_depth_model_max'] = \
                z_nodes[-1][ind_reset_max].min()

        # aft data resetting depth:
        aft_data_well['strat_age_interp'] = \
            np.interp(aft_data_well['depth'], z_nodes[-1], node_age)
        aft_data_well['reset'] = \
            aft_data_well['aft_age'] < aft_data_well['strat_age_interp']

        if True in aft_data_well['reset'].values:
            model_results.ix[model_scenario_number,
                             'resetting_depth_data_min'] = \
                np.min(aft_data_well['depth'][aft_data_well['reset']])
        if False in aft_data_well['reset'].values:
            model_results.ix[model_scenario_number,
                             'non-resetting_depth_data_max'] = \
                np.max(aft_data_well['depth'][aft_data_well['reset']==False])

        # save model run data to .pck file
        model_run_data = [
            time_array_bp,
            surface_temp_array, basal_hf_array,
            z_nodes, active_nodes, T_nodes,
            node_strat, node_age,
            aft_node_times_burial, aft_node_zs,
            aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
            aft_ln_mean_nodes, aft_ln_std_nodes,
            vr_nodes,
            T_data_well['depth'], T_data_well['temperature'],
            T_data_well['temperature_unc_1sigma'],
            vr_data_well['depth'], vr_data_well['VR'],
            vr_data_well['VR_unc_1sigma'],
            aft_data_well['depth'], aft_data_well['aft_age'],
            aft_data_well['95ci_minus'], aft_data_well['95ci_plus'],
            aft_data_well['length_mean'], aft_data_well['length_std'],
            aft_data_well['data_type'].values,
            T_gof, vr_gof, aft_age_gof]

        today = datetime.datetime.now()
        today_str = '%i-%i-%i' % (today.day, today.month, today.year)

        if pybasin_params.save_model_run_data is True:

            fn = os.path.join(datafile_output_dir, 'model_data_%s_%s_ms%i.pck'
                              % (well, today_str, model_scenario_number))

            print 'saving model run results to %s' % fn

            fout = open(fn, 'w')
            pickle.dump(model_run_data, fout)
            fout.close()

        #############################
        # make a model vs data figure
        #############################
        if pybasin_params.make_model_data_fig is True:
            fig = pybasin_figures.model_vs_data_figure(
                time_array_bp,
                surface_temp_array, basal_hf_array,
                z_nodes, active_nodes, T_nodes,
                node_strat, node_age,
                aft_node_times_burial, aft_node_zs,
                aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
                aft_ln_mean_nodes, aft_ln_std_nodes,
                vr_nodes,
                T_data_well['depth'].values,
                T_data_well['temperature'].values,
                T_data_well['temperature_unc_1sigma'].values,
                vr_data_well['depth'].values,
                vr_data_well['VR'].values,
                vr_data_well['VR_unc_1sigma'].values,
                aft_data_well['depth'],
                aft_data_well['aft_age'],
                aft_data_well['95ci_minus'],
                aft_data_well['95ci_plus'],
                aft_data_well['length_mean'],
                aft_data_well['length_std'],
                aft_data_well['data_type'].values,
                T_gof, vr_gof, aft_age_gof,
                model_data_fig_bw=pybasin_params.model_data_fig_bw)
        #    vr_data['depth'], vr_data['VR'], vr_data['unc_range_sigma'])

            fn = os.path.join(fig_output_dir, 'model_data_fig_%s_%s_ms%i.%s'
                              % (well, today_str,
                                 model_scenario_number,
                                 pybasin_params.fig_adj))
            print 'saving model-data comparison figure %s' % fn
            fig.savefig(fn, dpi=200)
            pl.clf()

        # save model results .csv file
        fn = os.path.join(output_dir, 'model_results_all_wells_%s_ms0-%i.csv'
                          % (today_str,
                             n_scenarios))
        print 'saving model results .csv file %s' % fn
        model_results.to_csv(fn, index_label='model_scenario_number')

        model_scenario_number += 1

print 'done'
