"""
summarize pybasin model results
min, max and best estimates for cooling and exhumation for each well or surface sample

"""

__author__ = 'elco'


import os
import sys
import itertools
import itertools
import numpy as np
import pandas as pd
import matplotlib.cm
import matplotlib.pyplot as pl
import useful_functions

x_data = 'max_cooling'
y_data = 'aft_age_gof'
y_data2 = 'vr_gof'
z_data = None

# minimum value for gof for acceptable model fits
gof_cutoff = 0.7

# best estimate of heat flow if no temperature data are available
default_heat_flow = 0.065

# read model result data
model_result_fn = "model_output/MB/final_results_17feb2016/model_results_all_wells_17-2-2016_ms0-14880_mod.csv"
df = pd.read_csv(model_result_fn)

#
df['cooling'] = df['mean_cooling_exhumation_phase_0']
df['cooling'][df['cooling'].isnull()] = 0.0

#
wells = np.unique(df['well'])
#

# calculate gof aft and temperature
df['aft_age_and_T_gof'] = (df['aft_age_gof'] + df['T_gof']) / 2.0

cols = ['gof_vr_max', 'gof_aft_max', 'gof_T_max', 'gof_aft_T_max',
        'cooling_vr_best', 'cooling_vr_min', 'cooling_vr_max',
        'cooling_aft_best', 'cooling_aft_min', 'cooling_aft_max',
        'cooling_at_best', 'cooling_at_min', 'cooling_at_max',
        'cooling_at_est_best', 'cooling_at_est_min', 'cooling_at_est_max',
        'exhumation_vr_best', 'exhumation_vr_min', 'exhumation_vr_max',
        'exhumation_aft_best', 'exhumation_aft_min', 'exhumation_aft_max',
        'exhumation_at_best', 'exhumation_at_min', 'exhumation_at_max',
        'exhumation_at_est_best', 'exhumation_at_est_min', 'exhumation_at_est_max',
        'heat_flow_best', 'heat_flow_min', 'heat_flow_max']

dfs = pd.DataFrame(columns=cols, index=wells)

for well in wells:
    ind = df['well'] == well

    dfs.ix[well, 'gof_vr_max'] = np.max(df['vr_gof'][ind])
    dfs.ix[well, 'gof_aft_max'] = np.max(df['aft_age_gof'][ind])
    dfs.ix[well, 'gof_T_max'] = np.max(df['T_gof'][ind])
    dfs.ix[well, 'gof_aft_T_max'] = np.max(df['aft_age_and_T_gof'][ind])

    # min., max values cooling & exhumation acc. to VR data:
    if True in df['vr_gof'][ind].notnull().values:
        ind_vr_max = (df['well'] == well) & (df['vr_gof'] == dfs.ix[well, 'gof_vr_max'])
        dfs.ix[well, 'cooling_vr_best'] = df['cooling'][ind_vr_max].values[0]

        ind_vr = (df['well'] == well) & (df['vr_gof'] >= gof_cutoff)
        dfs.ix[well, 'cooling_vr_min'] = np.min(df['cooling'][ind_vr])
        dfs.ix[well, 'cooling_vr_max'] = np.max(df['cooling'][ind_vr])
        dfs.ix[well, 'exhumation_vr_min'] = np.min(df['exhumation_magnitude'][ind_vr])
        dfs.ix[well, 'exhumation_vr_max'] = np.max(df['exhumation_magnitude'][ind_vr])

    ind_aft_max = (df['well'] == well) & (df['aft_age_gof'] == dfs.ix[well, 'gof_aft_max'])
    dfs.ix[well, 'cooling_aft_best'] = df['cooling'][ind_aft_max].values[0]

    ind_aft = (df['well'] == well) & (df['aft_age_gof'] >= gof_cutoff)
    dfs.ix[well, 'cooling_aft_min'] = np.min(df['cooling'][ind_aft])
    dfs.ix[well, 'cooling_aft_max'] = np.max(df['cooling'][ind_aft])
    dfs.ix[well, 'exhumation_aft_min'] = np.min(df['exhumation_magnitude'][ind_aft])
    dfs.ix[well, 'exhumation_aft_max'] = np.max(df['exhumation_magnitude'][ind_aft])

    # best, min, max fit for both AFT and temperature data:
    if True in df['aft_age_and_T_gof'][ind].notnull().values:
        ind_at_max = (df['well'] == well) & (df['aft_age_and_T_gof'] == dfs.ix[well, 'gof_aft_T_max'])
        dfs.ix[well, 'cooling_at_best'] = df['cooling'][ind_at_max].values[0]

        ind_at = (df['well'] == well) & (df['aft_age_and_T_gof'] >= gof_cutoff)
        dfs.ix[well, 'cooling_at_min'] = np.min(df['cooling'][ind_at])
        dfs.ix[well, 'cooling_at_max'] = np.max(df['cooling'][ind_at])
        dfs.ix[well, 'exhumation_at_min'] = np.min(df['exhumation_magnitude'][ind_at])
        dfs.ix[well, 'exhumation_at_max'] = np.max(df['exhumation_magnitude'][ind_at])

        # copy to columns with combined T data or default heat flow:
        dfs.ix[well, 'cooling_at_est_best'] = dfs.ix[well, 'cooling_at_best']
        dfs.ix[well, 'cooling_at_est_min'] = dfs.ix[well, 'cooling_at_min']
        dfs.ix[well, 'cooling_at_est_max'] = dfs.ix[well, 'cooling_at_max']
        dfs.ix[well, 'exhumation_at_est_best'] = dfs.ix[well, 'exhumation_at_best']
        dfs.ix[well, 'exhumation_at_est_min'] = dfs.ix[well, 'exhumation_at_min']
        dfs.ix[well, 'exhumation_at_est_max'] = dfs.ix[well, 'exhumation_at_max']

        # record good values of heat flow for both AFT and T data
        dfs.ix[well, 'heat_flow_min'] = np.min(df['basal_heat_flow'][ind_at])
        dfs.ix[well, 'heat_flow_max'] = np.max(df['basal_heat_flow'][ind_at])

    else:

        # find combination of good aft fit and estimated default value of heat flow:
        ind_at_est = (df['well'] == well) & (df['aft_age_gof'] >= gof_cutoff) \
                     & (df['basal_heat_flow'] == default_heat_flow)

        dfs.ix[well, 'heat_flow_min'] = default_heat_flow
        dfs.ix[well, 'heat_flow_max'] = default_heat_flow

        if ind_at_est.sum() == 0:
            # no good fit for default heat flow: check which values of heat flow do generate good aft fit
            ind_at_est = (df['well'] == well) & (df['aft_age_gof'] >= gof_cutoff)

            dfs.ix[well, 'heat_flow_min'] = np.min(df['basal_heat_flow'][ind_at_est])
            dfs.ix[well, 'heat_flow_max'] = np.max(df['basal_heat_flow'][ind_at_est])

        aft_gof_best = np.max(df['aft_age_gof'][ind_at_est])

        ind_at_est_best = ind_at_est & (df['aft_age_gof'] == aft_gof_best)

        dfs.ix[well, 'cooling_at_est_best'] = np.max(df['cooling'][ind_at_est_best])
        dfs.ix[well, 'cooling_at_est_min'] = np.min(df['cooling'][ind_at_est])
        dfs.ix[well, 'cooling_at_est_max'] = np.max(df['cooling'][ind_at_est])
        dfs.ix[well, 'exhumation_at_est_best'] = np.max(df['exhumation_magnitude'][ind_at_est_best])
        dfs.ix[well, 'exhumation_at_est_min'] = np.min(df['exhumation_magnitude'][ind_at_est])
        dfs.ix[well, 'exhumation_at_est_max'] = np.max(df['exhumation_magnitude'][ind_at_est])

# save summary
summary_fn = model_result_fn.split('.csv')[0] + '_cooling_and_exhumation_summary.csv'
print 'saving results to %s' % summary_fn
dfs.to_csv(summary_fn, index_label='well')

print 'done'