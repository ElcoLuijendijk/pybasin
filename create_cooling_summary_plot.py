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

degree_symbol = unichr(176)
s = 15
alpha = 0.5
gof_cutoff = 0.7
default_heat_flow = 0.065


#x_datas = ['cooling', 'exhumation_start']
x_data = 'cooling'
x_label = 'Cooling (%s C)' % degree_symbol
xticks = [0, 25, 50, 75, 100, 125, 150, 175]
fn_adj = 'cooling_vs_thermochron'

#x_data = 'exhumation_start'
#x_label = 'Start of cooling (Ma)'
#xticks = [0, 2.5, 5.0, 7.5, 10.0, 12.5]
#fn_adj = 'exh_start_vs_thermochron'

x_data = 'exhumation_time_factor'
x_label = 'Cooling (%s C)' % degree_symbol
xticks = [0, 25, 50, 75, 100, 125, 150, 175]
fn_adj = 'cooling_vs_thermochron'


x_datas = ['cooling', 'exhumation_start',
           'exhumation_segment_factor', 'exhumation_duration_factor']
x_labels = ['Cooling (%s C)' % degree_symbol, 'Start of cooling (Ma)',
            'exhumation segment factor', 'exhumation duration factor']
xticks_all = [[0, 25, 50, 75, 100, 125, 150, 175],
              [0, 2.5, 5.0, 7.5, 10.0, 12.5],
              [0.0, 0.25, 0.5, 0.75, 1.0],
              [0.0, 0.25, 0.5, 0.75, 1.0]]
fn_adjs = ['cooling_vs_thermochron',
           'exh_start_vs_thermochron',
           'exh_time_factor_vs_thermochron',
           'exhumation rate factor']

#model_result_fn = 'model_output/MB/5apr2016_results_two_stage_cooling/' \
#                  'model_results_merged_mod.csv'
#model_result_fn = '/home/elco/model_files/pybasin/MB/' \
#                  'final_results_1mar2016_rdaam/' \
#                  'model_results_merged_mod.csv'

model_result_fn = '/home/elco/model_files/pybasin/MB/30aug2016_2stage_new/' \
                  'model_results_merged_mod.csv'

df_all = pd.read_csv(model_result_fn)


for x_data, x_label, xticks, fn_adj in zip(x_datas, x_labels,
                                           xticks_all, fn_adjs):

    #
    #x_label = 'Duration of cooling (My)'
    y_datas = ['aft_age_gof', 'ahe_gof']
    markers = ['o', 's']
    colors = ['lightblue', 'gray']

    #model_result_fn = 'model_output/MB/5apr2016_results_two_stage_cooling/' \
    #                  'model_results_merged_mod.csv'
    #model_result_fn = '/home/elco/model_files/pybasin/MB/' \
    #                  'final_results_1mar2016_rdaam/' \
    #                  'model_results_merged_mod.csv'

    #df_all = pd.read_csv(model_result_fn)

    wells = np.unique(df_all['well'].dropna())
    nwells = len(wells)

    ncols = 4
    nrows = int(np.ceil(len(wells) / float(ncols)))

    width = 8.0
    height = width * float(nrows)/ncols

    fig, panels_all = pl.subplots(nrows, ncols, figsize=(width, height),
                                  sharey=True, sharex=False)

    panels = panels_all.ravel()

    vmin = 0
    vmax = 1.0

    plot_count = 0

    for plot_count, well in enumerate(wells):

        print well

        df = df_all[df_all['well'] == well]

        leg_ydatas = []

        for y_data, marker, color in zip(y_datas, markers, colors):
            ly = panels[plot_count].scatter(df[x_data], df[y_data],
                                            marker=marker, facecolor=color,
                                            s=s, lw=0.5, alpha=alpha)
            leg_ydatas.append(ly)

        panels[plot_count].set_xlim(0, df[x_data].max() * 1.1)

        panels[plot_count].set_title(well)

    for i in range(nrows):
        panels[i * ncols].set_ylabel('Goodness of fit')

    for panel, well in zip(panels, wells):
        panel.xaxis.grid()
        panel.yaxis.set_ticks_position('left')
        panel.xaxis.set_ticks_position('bottom')
        panel.spines['right'].set_visible(False)
        panel.spines['top'].set_visible(False)
        leg_gf = panel.axhline(y=0.5, color='black', ls='--', lw=1.0, zorder=0)
        panels[plot_count].set_ylim(0, 1.02)

        xmin, xmax = panel.get_xlim()
        #x_int = 25.0
        #xticks = np.arange(0, xmax+x_int, x_int).astype(int)

        panel.set_xticks(xticks)
        panel.set_xlim(0, xticks[-1])

    if len(panels) > len(wells):
        for panel in panels[len(wells):]:
            panel.axis('off')

    for col_no in range(ncols):

        panels[-(col_no + 1)].set_xlabel(x_label)

    legs = leg_ydatas + [leg_gf]
    labels = ['model fit, AFT data', 'model fit, AHe data', 'good fit']
    #fig.legend(legs, labels, loc='lower right', fontsize='small')
    panels[-1].legend(legs, labels, bbox_to_anchor=(1.05, 1),
                              fontsize='small', frameon=False,
                              handlelength=3)

    fig_dir = os.path.join(os.path.dirname(model_result_fn), 'model_data_fig')
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)

    base_fn = os.path.basename(model_result_fn)
    base_fn_new = base_fn.split('.csv')[0] + '_%s.png' % fn_adj
    fn = os.path.join(fig_dir, base_fn_new)

    fig.tight_layout()

    print 'saving %s' % fn

    fig.savefig(fn, dpi=200)

print 'done'
