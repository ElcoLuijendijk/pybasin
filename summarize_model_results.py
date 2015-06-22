"""
create a series of 2D contour plots of 2 parameter values vs goodness of
fit statistic

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

gof_cutoff = 0.7

default_heat_flow = 0.065

# read model result data
model_result_fn = 'model_output/MB/' \
                  'model_results_all_wells_29-3-2015_ms0--9828_mod.csv'
df = pd.read_csv(model_result_fn)

wells = np.unique(df['well'])
#

# calculate gof aft and temperature
df['aft_age_and_T_gof'] = (df['aft_age_gof'] + df['T_gof']) / 2.0


cols = ['gof_vr_max', 'gof_aft_max',
        'cooling_vr_best', 'cooling_vr_min', 'cooling_vr_max',
        'cooling_aft_best', 'cooling_aft_min', 'cooling_aft_max',
        'cooling_at_best', 'cooling_at_min', 'cooling_at_max',
        'exhumation_vr_best', 'exhumation_vr_min', 'exhumation_vr_max',
        'exhumation_aft_best', 'exhumation_aft_min', 'exhumation_aft_max',
        'exhumation_at_best', 'exhumation_at_min', 'exhumation_at_max']

dfs = pd.DataFrame(columns=cols, index=wells)

for well in wells:
    ind = df['well'] == well

    dfs.ix[well, 'gof_vr_max'] = np.max(df['vr_gof'][ind])
    dfs.ix[well, 'gof_aft_max'] = np.max(df['aft_age_gof'][ind])
    dfs.ix[well, 'gof_aft_T_max'] = np.max(df['aft_age_and_T_gof'][ind])

    if True in df['vr_gof'][ind].notnull().values:
        ind_vr_max = (df['well'] == well) & (df['vr_gof'] == dfs.ix[well, 'gof_vr_max'])
        dfs.ix[well, 'cooling_vr_best'] = df['max_cooling'][ind_vr_max].values[0]

        ind_vr = (df['well'] == well) & (df['vr_gof'] >= gof_cutoff)
        dfs.ix[well, 'cooling_vr_min'] = np.min(df['max_cooling'][ind_vr])
        dfs.ix[well, 'cooling_vr_max'] = np.max(df['max_cooling'][ind_vr])
        dfs.ix[well, 'exhumation_vr_min'] = np.min(df['exhumation_magnitude'][ind_vr])
        dfs.ix[well, 'exhumation_vr_max'] = np.max(df['exhumation_magnitude'][ind_vr])

    ind_aft_max = (df['well'] == well) & (df['aft_age_gof'] == dfs.ix[well, 'gof_aft_max'])
    dfs.ix[well, 'cooling_aft_best'] = df['max_cooling'][ind_aft_max].values[0]

    ind_aft = (df['well'] == well) & (df['aft_age_gof'] >= gof_cutoff)
    dfs.ix[well, 'cooling_aft_min'] = np.min(df['max_cooling'][ind_aft])
    dfs.ix[well, 'cooling_aft_max'] = np.max(df['max_cooling'][ind_aft])
    dfs.ix[well, 'exhumation_aft_min'] = np.min(df['exhumation_magnitude'][ind_aft])
    dfs.ix[well, 'exhumation_aft_max'] = np.max(df['exhumation_magnitude'][ind_aft])

    if True in df['aft_age_and_T_gof'][ind].notnull().values:
        ind_at_max = (df['well'] == well) & (df['aft_age_and_T_gof'] == dfs.ix[well, 'gof_aft_T_max'])
        dfs.ix[well, 'cooling_at_best'] = df['max_cooling'][ind_at_max].values[0]

        ind_at = (df['well'] == well) & (df['aft_age_and_T_gof'] >= gof_cutoff)
        dfs.ix[well, 'cooling_at_min'] = np.min(df['max_cooling'][ind_at])
        dfs.ix[well, 'cooling_at_max'] = np.max(df['max_cooling'][ind_at])
        dfs.ix[well, 'exhumation_at_min'] = np.min(df['exhumation_magnitude'][ind_at])
        dfs.ix[well, 'exhumation_at_max'] = np.max(df['exhumation_magnitude'][ind_at])


# save summary
summary_fn = model_result_fn.split('.csv')[0] + '_summary.csv'
dfs.to_csv(summary_fn, index_label='well')

print 'done'