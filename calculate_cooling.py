"""
read pybasin model results and calculate cooling at latest exhumation phase
from recorded max and present temperatures

also records the most likely basal heat flow if temperature data are available
and otherwise uses a best estimate heat flow

"""

__author__ = 'elco'


import os
import numpy as np
import pandas as pd


min_T_gof = 0.5
default_heat_flow = 0.065

model_result_fn = "/home/elco/python_scripts/pybasin/model_output/MB/5apr2016_results_two_stage_cooling/model_results_merged.csv"

df = pd.read_csv(model_result_fn)

wells = np.unique(df['well'])

# select cooling values and set NaN to 0 cooling:
df['cooling'] = df['mean_cooling_exhumation_phase_0']
df['cooling'][df['cooling'].isnull()] = 0.0

df['present_max_T'] = np.nan
df['estimated_present_day_heat_flow'] = np.nan
df['cooling_with_fixed_present_HF'] = np.nan
df['max_T_gof'] = np.nan

for well in wells:

    print 'well ', well

    well_ind = df['well'] == well

    # find best fit present-day temperature
    if True in df['T_gof'][well_ind].notnull().values and df['T_gof'][well_ind].max() > min_T_gof:
        best_T_gof_ind = df['T_gof'][well_ind] == df['T_gof'][well_ind].max()

        df['max_T_gof'] = df['T_gof'][well_ind].max()

        print 'using T data, with max gof of ', df['T_gof'][well_ind].max()

    else:
        best_T_gof_ind = df['basal_heat_flow'][well_ind] == default_heat_flow
        df['estimated_present_day_heat_flow'][well_ind] = default_heat_flow

        print 'using default heat flow ', default_heat_flow

    present_max_Ts = (df['max_temperature'][well_ind][best_T_gof_ind] -
                      df['cooling'][well_ind][best_T_gof_ind])
    present_max_T = present_max_Ts.mean()

    df['present_max_T'][well_ind] = present_max_T

    df['cooling_with_fixed_present_HF'][well_ind] = df['max_temperature'][well_ind] - present_max_T

#
df['cooling'] = df['cooling_with_fixed_present_HF']

# save results
mod_fn = model_result_fn.split('.csv')[0] + '_mod.csv'
print 'saving results to ', mod_fn
df.to_csv(mod_fn, index=False)
