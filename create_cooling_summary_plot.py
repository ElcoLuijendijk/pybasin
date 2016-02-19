"""
create three scatter plots of cooling, exhumation start and duration per sample
location

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

x_datas = ['cooling', 'exhumation_start', 'exhumation_duration']
y_datas = ['vr_gof', 'aft_age_gof', 'ahe_gof']
markers = ['v', 'o', 's']
colors = ['brown', 'lightblue', 'gray']

degree_symbol = unichr(176)
s = 15
alpha = 0.5

gof_cutoff = 0.7

default_heat_flow = 0.065

fn_adj = 'cooling_vs_aft'

model_result_fn = 'model_output/MB/final_results_18feb2016_v2/' \
                  'model_results_merged_mod.csv'

df_all = pd.read_csv(model_result_fn)

wells = np.unique(df_all['well'].dropna())
nwells = len(wells)

ncols = 3
nrows = len(wells)

width = 8.0
height = width * float(nrows)/ncols

fig, panels_all = pl.subplots(nrows, ncols, figsize=(width, height), sharey=True, sharex=False)

panels = panels_all.ravel()

vmin = 0
vmax = 1.0

plot_count = 0

for well in wells:

    print well

    df = df_all[df_all['well'] == well]

    panels[plot_count].set_ylabel('Goodness of fit')
    panels[plot_count].set_ylim(0, 1)
    panels[plot_count].set_title(well)

    for col_no, x_data in enumerate(x_datas):
        for y_data, marker, color in zip(y_datas, markers, colors):
            panels[plot_count].scatter(df[x_data], df[y_data],
                                       marker=marker, facecolor=color,
                                       s=s, lw=0.5, alpha=alpha)
        panels[plot_count].set_xlim(0, df[x_data].max() * 1.1)

        plot_count += 1

for panel in panels:
    panel.xaxis.grid()
    panel.yaxis.set_ticks_position('left')
    panel.xaxis.set_ticks_position('bottom')
    panel.spines['right'].set_visible(False)
    panel.spines['top'].set_visible(False)
    panel.axhline(y=0.7, color='gray', ls='--', lw=1.0, zorder=0)

for col_no, x_data in enumerate(x_datas[::-1]):

    if 'cooling' in x_data:
        x_label = 'Cooling (%s C)' % degree_symbol
    elif 'start' in x_data:
        x_label = 'Start of cooling (Ma)'
    elif 'duration' in x_data:
        x_label = 'Duration of cooling (My)'
    else:
        x_label = x_data

    panels[-(col_no + 1)].set_xlabel(x_label)

fig_dir = os.path.join(os.path.dirname(model_result_fn), 'model_data_fig')
if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)

base_fn = os.path.basename(model_result_fn)
base_fn_new = base_fn.split('.csv')[0] + '_%s.png' % fn_adj
fn = os.path.join(fig_dir, base_fn_new)

fig.savefig(fn, dpi=200)

print 'done'
