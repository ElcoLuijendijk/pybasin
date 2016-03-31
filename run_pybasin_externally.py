"""
experimental minimal example for running PyBasin as an external module

"""


__author__ = 'Elco Luijendijk'

import sys
import pdb
import numpy as np
import pandas as pd
import PyBasin as pb

input_dir = 'input_data/AFT_benchmarks_Texas'
well = 'Frio'

sys.path.append(input_dir)
import pybasin_params

output_dir = 'model_output/salinity'
return_objective_function = True
record_data = True
model_scenario_number = 1
well_number = 1
csv_output_dir = None
n_scenarios = 1

(well_strats, strat_info_mod,
 salinity_bnd_df,
 T_data_df, vr_data_df,
 aft_samples, aft_ages,
 ahe_samples, ahe_data,
 salinity_data, surface_temp, litho_props) \
    = pb.read_model_input_data(input_dir, pybasin_params)

well_strat, well_strat_orig = pb.select_well_strat(well, well_strats)

if (pybasin_params.simulate_salinity is True
        and pybasin_params.well_specific_surface_salinity_bnd is True):

    surface_salinity_well = \
        pb.select_well_salinity_bnd(well, salinity_bnd_df)
else:
    surface_salinity_well = False

model_results_df, model_results_df2 = \
    pb.setup_model_output_df(n_scenarios)

calibration_target = ['AFT_age']

params_to_change = ['AFT_C0', 'AFT_C1',
                    'AFT_C2', 'AFT_C3',
                    'AFT_alpha']

model_scenario_params = [0.39528, 0.01073, -65.12969, -7.91715, 0.04672]

param_bounds_min = None
param_bounds_max = None

exhumation_scenarios_period = 'dummy_exhumation'
basal_heat_flow_scenario_period = 'all'

args = (model_scenario_params,
        params_to_change,
        exhumation_scenarios_period,
        basal_heat_flow_scenario_period,
        model_results_df,
        model_results_df2,
        well_number, well,
        well_strat.copy(), well_strat_orig,
        strat_info_mod,
        pybasin_params,
        surface_temp,
        surface_salinity_well,
        litho_props,
        csv_output_dir,
        output_dir,
        model_scenario_number,
        return_objective_function,
        calibration_target,
        record_data,
        param_bounds_min,
        param_bounds_max,
        T_data_df, vr_data_df,
        aft_samples, aft_ages,
        ahe_samples, ahe_data,
        salinity_data)

obj_function = pb.update_model_params_and_run_model(*args)

print 'done'




