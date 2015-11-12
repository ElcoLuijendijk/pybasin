"""
create a series of scatter plots of cooling vs goodness of fit, with one subplot for each combination and start
and exhumation.

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

x_data = 'cooling'
y_data = 'aft_age_gof'
y_data2 = 'vr_gof'
z_data = None

s = 40

gof_cutoff = 0.7

default_heat_flow = 0.065

fn_adj = 'cooling_vs_aft'

model_result_fn = 'model_output/MB/final_results_10aug2015/' \
                  'model_results_all_wells_16-7-2015_ms0-15360_mod.csv'

df_all = pd.read_csv(model_result_fn)

# find all combinations of age and duration of exhumation event
ad_comb = zip(df_all['exhumation_start'], df_all['exhumation_duration'])

unique_ad_comb = list(set(ad_comb))
unique_ad_comb.sort()

ncombs = len(unique_ad_comb)

ncols = 3
nrows = int(np.ceil((len(unique_ad_comb)) / ncols))

# cooling figure
if z_data is not None and 'gof' in z_data:
    vmin = 0
    vmax = 1.0
else:
    vmin = None
    vmax = None

wells = np.unique(df_all['well'])

for well in wells:

    print well

    df = df_all[df_all['well'] == well]

    fig = useful_functions.setup_figure(width='2col',
                                        height=float(nrows)/ncols)

    panels = [fig.add_subplot(nrows, ncols, i)
              for i in xrange(1, ncombs + 1)]

    for panel, ad_comb_i in zip(panels, unique_ad_comb):

        print ad_comb_i

        ind = np.logical_and(df['exhumation_start'] == ad_comb_i[0],
                             df['exhumation_duration'] == ad_comb_i[1])

        if z_data is not None and 'gof' in z_data:
            ind_pass = np.logical_and(df[z_data] >= gof_cutoff, ind)
            ind_fail = np.logical_and(df[z_data] < gof_cutoff, ind)

            sc = panel.scatter(df[x_data][ind_pass], df[y_data][ind_pass],
                               c=df[z_data][ind_pass],
                               edgecolor='None',
                               s=60,
                               marker='o',
                               vmin=vmin, vmax=vmax,
                               cmap=matplotlib.cm.jet_r)

            sc = panel.scatter(df[x_data][ind_fail], df[y_data][ind_fail],
                               c=df[z_data][ind_fail],
                               #edgecolor='None',
                               s=30,
                               marker='x',
                               vmin=vmin, vmax=vmax,
                               cmap=matplotlib.cm.jet_r)

        else:
            if z_data is not None:
                c = df[z_data][ind]
                edgecolor = 'None'
            else:
                c = 'black'
                edgecolor='black'
            sc = panel.scatter(df[x_data][ind], df[y_data][ind], c=c,
                               edgecolor=edgecolor,
                               s=s,
                               vmin=vmin, vmax=vmax,
                               cmap=matplotlib.cm.jet_r,
                               zorder=10)

            if y_data2 is not None:
                panelr = panel.twinx()
                sc2 = panel.scatter(df[x_data][ind].values, df[y_data2][ind].values, c='lightgrey',
                                     edgecolor=edgecolor,
                                     s=s * 2/3.,
                                     marker='s',
                                     vmin=vmin, vmax=vmax,
                                     cmap=matplotlib.cm.jet_r,
                                     zorder=1)

                #offset = 60
                #new_fixed_axis = panel.get_grid_helper().new_fixed_axis
                panelr.yaxis.tick_left()
                panelr.yaxis.set_label_position('left')
                #panelr.set_offset_position(self, position)
                #panelr.axis["left2"] = new_fixed_axis(loc="left",
                #                                      axes=panelr,
                #                                      offset=(offset, 0))

                #panelr.axis["left2"].toggle(all=True)
                panelr.set_ylabel(y_data2.replace('_', ' '), labelpad=12)
                #panelr.spines["left"].label.set_color('brown')
                panelr.set_ylim(0, 1.0)
                panelr.yaxis.label.set_color('gray')

                #panelr.spines['left'].set_position(("axes", -0.4))

                # place second y axis and dataset below first
                #panel.set_zorder(panelr.get_zorder()+1) # put ax in front of ax2
                #panel.patch.set_visible(False) # hide the 'canvas'

        panel.grid()

        panel.set_xlabel(x_data.replace('_', ' '))
        panel.set_ylabel(y_data.replace('_', ' '))

        tekst = 'start exhumation = %0.0f Ma\nduration exhumation = %0.0f Ma' \
                % (ad_comb_i[0], ad_comb_i[1])
        panel.set_title(tekst, fontsize='x-small')

        if 'gof' in y_data:
            panel.set_ylim(0, 1.0)
            panel.axhline(y=gof_cutoff, color='black', lw=1.0)

            gof_ind = df[y_data][ind] > gof_cutoff
            x_min = np.min(df[x_data][ind][gof_ind])
            x_max = np.max(df[x_data][ind][gof_ind])

            panel.fill_between((x_min, x_max), (1.0, 1.0),
                               color='lightgrey', zorder=0)
            #print bla

    fig.tight_layout()

    if z_data is not None:
        cax = fig.add_axes([0.3, 0.09, 0.4, 0.015])
        cb = fig.colorbar(sc, cax=cax, orientation='horizontal')
        cb_label = z_data.replace('_', ' ')
        cb.set_label(cb_label)
        #cb.set_ticks(cax.get_xticks()[::2])

        fig.subplots_adjust(bottom=0.18)
    
    for panel in panels:
        panel.set_xlim(df[x_data].min(), df[x_data].max())

    fig_dir = os.path.join(os.path.dirname(model_result_fn), 'model_data_fig')
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)

    base_fn = os.path.basename(model_result_fn)
    base_fn_new = base_fn.split('.csv')[0] + '_%s_%s.png' % (well, fn_adj)
    fn = os.path.join(fig_dir, base_fn_new)

    fig.savefig(fn, dpi=300)

