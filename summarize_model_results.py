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
import pdb

x_data = 'max_cooling'
y_data = 'aft_age_gof'
y_data2 = 'vr_gof'
z_data = None

# minimum value for gof for acceptable model fits
gof_cutoff = 0.7

# best estimate of heat flow if no temperature data are available
default_heat_flow = 0.065

#
default_exhumation = 2000.0

# read model result data
model_result_fn = "model_output/MB/final_results_19feb2016/model_results_merged_mod.csv"
df = pd.read_csv(model_result_fn)

#
#df['cooling'] = df['mean_cooling_exhumation_phase_0']
#df['cooling'][df['cooling'].isnull()] = 0.0

#
wells = np.unique(df['well'].dropna())
#

# calculate gof aft and temperature
df['aft_age_and_T_gof'] = (df['aft_age_gof'] + df['T_gof']) / 2.0

cols = ['gof_T_best', 'gof_vr_best', 'gof_aft_best', 'gof_ahe_best',
        'cooling_vr_best', 'cooling_vr_min', 'cooling_vr_max',
        'cooling_aft_best', 'cooling_aft_min', 'cooling_aft_max',
        'cooling_ahe_best', 'cooling_ahe_min', 'cooling_ahe_max',
        'exhumation_vr_best', 'exhumation_vr_min', 'exhumation_vr_max',
        'exhumation_aft_best', 'exhumation_aft_min', 'exhumation_aft_max',
        'exhumation_ahe_best', 'exhumation_ahe_min', 'exhumation_ahe_max',
        'heat_flow_vr_best', 'heat_flow_vr_min', 'heat_flow_vr_max',
        'heat_flow_aft_best', 'heat_flow_aft_min', 'heat_flow_aft_max',
        'heat_flow_ahe_best', 'heat_flow_ahe_min', 'heat_flow_ahe_max']

dfs = pd.DataFrame(columns=cols, index=wells)

for well in wells:
    ind = df['well'] == well

    dfs.ix[well, 'gof_T_best'] = np.max(df['T_gof'][ind])

    # min., max values cooling & exhumation acc. to VR data:
    if True in df['vr_gof'][ind].notnull().values:
        dfs.ix[well, 'gof_vr_best'] = np.max(df['vr_gof'][ind])

        ind_vr_max = (df['well'] == well) & (df['vr_gof'] == dfs.ix[well, 'gof_vr_best'])
        dfs.ix[well, 'cooling_vr_best'] = df['cooling'][ind_vr_max].values[0]

        ind_vr = (df['well'] == well) & (df['vr_gof'] >= gof_cutoff)
        dfs.ix[well, 'cooling_vr_min'] = np.min(df['cooling'][ind_vr])
        dfs.ix[well, 'cooling_vr_max'] = np.max(df['cooling'][ind_vr])
        dfs.ix[well, 'exhumation_vr_min'] = np.min(df['exhumation_magnitude'][ind_vr])
        dfs.ix[well, 'exhumation_vr_max'] = np.max(df['exhumation_magnitude'][ind_vr])
        dfs.ix[well, 'heat_flow_vr_min'] = np.min(df['basal_heat_flow'][ind_vr])
        dfs.ix[well, 'heat_flow_vr_max'] = np.max(df['basal_heat_flow'][ind_vr])


    ind_aft = (df['well'] == well) & (df['aft_age_gof'] >= gof_cutoff)
    if True in df['aft_age_gof'][ind].notnull().values:

        dfs.ix[well, 'gof_aft_best'] = np.max(df['aft_age_gof'][ind])
                
        ind_aft_max = (df['well'] == well) & (df['aft_age_gof'] == dfs.ix[well, 'gof_aft_best'])
        dfs.ix[well, 'cooling_aft_best'] = df['cooling'][ind_aft_max].values[0]
        dfs.ix[well, 'exhumation_aft_best'] = df['exhumation_magnitude'][ind_aft_max].values[0]
        dfs.ix[well, 'heat_flow_aft_best'] = df['basal_heat_flow'][ind_aft_max].values[0]

        dfs.ix[well, 'cooling_aft_min'] = np.min(df['cooling'][ind_aft])
        dfs.ix[well, 'cooling_aft_max'] = np.max(df['cooling'][ind_aft])
        dfs.ix[well, 'exhumation_aft_min'] = np.min(df['exhumation_magnitude'][ind_aft])
        dfs.ix[well, 'exhumation_aft_max'] = np.max(df['exhumation_magnitude'][ind_aft])
        dfs.ix[well, 'heat_flow_aft_min'] = np.min(df['basal_heat_flow'][ind_aft])
        dfs.ix[well, 'heat_flow_aft_max'] = np.max(df['basal_heat_flow'][ind_aft])

        # find value of heat flow for best GOF assuming exhumation is known:
        ind_aft_exh = (df['well'] == well) & (df['exhumation_magnitude'] == default_exhumation)
        if True in ind_aft_exh.values:
            ind_aft_exh_max_rel = np.argmax(df['aft_age_gof'][ind_aft_exh])
            dfs.ix[well, 'gof_aft_exh%0.0f_best' % default_exhumation] = np.max(df['aft_age_gof'][ind_aft_exh])
            dfs.ix[well, 'heat_flow_aft_exh%0.0f_best' % default_exhumation] = df.ix[ind_aft_exh_max_rel, 'basal_heat_flow']

            ind_default_exhumation = ind_aft & (df['exhumation_magnitude'] == default_exhumation)
            dfs.ix[well, 'heat_flow_aft_exh%0.0f_min' % default_exhumation] = np.min(df['basal_heat_flow'][ind_default_exhumation])
            dfs.ix[well, 'heat_flow_aft_exh%0.0f_max' % default_exhumation] = np.max(df['basal_heat_flow'][ind_default_exhumation])

            print well, dfs.ix[well, 'heat_flow_aft_exh%0.0f_min' % default_exhumation]
        else:
            print well, 'all gof for exh %0.0f below cutoff %0.2f' % (default_exhumation, gof_cutoff)

    ind_ahe = (df['well'] == well) & (df['ahe_gof'] >= gof_cutoff)
    dfs.ix[well, 'gof_ahe_best'] = np.max(df['ahe_gof'][ind])

    if True in df['ahe_gof'][ind].notnull().values:

        ind_ahe_max = (df['well'] == well) & (df['ahe_gof'] == dfs.ix[well, 'gof_ahe_best'])
        dfs.ix[well, 'cooling_ahe_best'] = df['cooling'][ind_ahe_max].values[0]
        dfs.ix[well, 'exhumation_ahe_best'] = df['exhumation_magnitude'][ind_ahe_max].values[0]
        dfs.ix[well, 'heat_flow_ahe_best'] = df['basal_heat_flow'][ind_ahe_max].values[0]

        dfs.ix[well, 'cooling_ahe_min'] = np.min(df['cooling'][ind_ahe])
        dfs.ix[well, 'cooling_ahe_max'] = np.max(df['cooling'][ind_ahe])
        dfs.ix[well, 'exhumation_ahe_min'] = np.min(df['exhumation_magnitude'][ind_ahe])
        dfs.ix[well, 'exhumation_ahe_max'] = np.max(df['exhumation_magnitude'][ind_ahe])
        dfs.ix[well, 'heat_flow_ahe_min'] = np.min(df['basal_heat_flow'][ind_ahe])
        dfs.ix[well, 'heat_flow_ahe_max'] = np.max(df['basal_heat_flow'][ind_ahe])

        # find value of heat flow for best GOF assuming exhumation is known:
        ind_ahe_exh = (df['well'] == well) & (df['exhumation_magnitude'] == default_exhumation)
        if True in ind_ahe_exh.values:
            ind_ahe_exh_max_rel = np.argmax(df['ahe_gof'][ind_ahe_exh])
            dfs.ix[well, 'gof_ahe_exh%0.0f_best' % default_exhumation] = np.max(df['ahe_gof'][ind_ahe_exh])
            dfs.ix[well, 'heat_flow_ahe_exh%0.0f_best' % default_exhumation] = df.ix[ind_ahe_exh_max_rel, 'basal_heat_flow']
            ind_default_exhumation = ind_ahe & (df['exhumation_magnitude'] == default_exhumation)
            dfs.ix[well, 'heat_flow_ahe_exh%0.0f_min' % default_exhumation] = np.min(df['basal_heat_flow'][ind_default_exhumation ])
            dfs.ix[well, 'heat_flow_ahe_exh%0.0f_max' % default_exhumation] = np.max(df['basal_heat_flow'][ind_default_exhumation ])

# save summary
summary_fn = model_result_fn.split('.csv')[0] + '_cooling_and_exhumation_summary_cutoff%0.2f.csv' % gof_cutoff
print 'saving results to %s' % summary_fn
dfs.to_csv(summary_fn, index_label='well')

print 'done'