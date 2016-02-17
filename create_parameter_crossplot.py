"""
create a series of 2D contour plots of 2 parameter values vs goodness of
fit statistic

"""

__author__ = 'elco'


import os
import sys
import itertools
import numpy as np
import pandas as pd
import matplotlib.cm
import matplotlib.pyplot as pl
import useful_functions


# check if script dir in python path
scriptdir = os.path.realpath(sys.path[0])
if scriptdir not in sys.path:
    sys.path.append(scriptdir)

print ''

if len(sys.argv) > 1 and sys.argv[-1][-3:] != '.py':
    scenario_name = sys.argv[-1]
    model_output_subfolder = os.path.join(scriptdir, 'model_output',
                                          scenario_name)
    print 'running model input data from folder %s' % model_output_subfolder

else:
    # read default input folder
    fin = open(os.path.join(scriptdir, 'default_input_folder.txt'))
    d = fin.readline()
    fin.close()
    scenario_name = d.split()[0]
    model_output_subfolder = os.path.join(scriptdir, 'model_output',
                                          scenario_name)
    print 'running model input data from folder %s' % model_output_subfolder


# read all model result summary files in folder
all_files = os.listdir(model_output_subfolder)

model_results_files = [fn for fn in all_files
                       if 'model_results_all_wells' in fn
                       and fn[-4:] == '.csv']

model_result_fns = [os.path.join(model_output_subfolder, fn)
                    for fn in model_results_files]

model_result_fns = ['model_output/MB/final_results_17feb2016/model_results_all_wells_17-2-2016_ms0-14880_mod.csv']

# go through all model files:
for model_result_fn in model_result_fns:

    print 'reading %s' % model_result_fn

    df = pd.read_csv(model_result_fn)

    all_cols = df.columns.tolist()

    gof_cols = [col for col in all_cols if 'gof' in col]

    param_cols = [col for col in all_cols
                  if 'gof' not in col
                  and 'model' not in col
                  and 'well' not in col]

    # set units
    units = []
    multipliers = np.ones(len(param_cols))

    for i, param_col in enumerate(param_cols):
        if 'heat_flow' in param_col:
            units.append('mW m$^{-2}$')
            multipliers[i] = 1000.0
        elif 'exhumation' in param_col:
            units.append('km')
            multipliers[i] = 1e-3
        elif 'age' in param_col:
            units.append('Ma')
        elif 'duration' in param_col:
            units.append('My')
        else:
            units.append(None)

    for param_col, multiplier in zip(param_cols, multipliers):
        df[param_col] *= multiplier

    # override param cols to restrict the number of figures
    param_cols = ['exhumation_magnitude', 'basal_heat_flow']

    # find number of wells
    wells = np.unique(df['well'])
    n_wells = len(wells)

    n_cols = 3
    n_rows = int(np.ceil(float(n_wells + 1) / n_cols))

    # go through all parameter combinations
    for i, param_combination, unit_combination in \
            zip(itertools.count(),
                itertools.combinations(param_cols, 2),
                itertools.combinations(units, 2)):

        for gof_col in gof_cols:

            ############################

            width = 8.0
            height=width * float(n_rows)/n_cols
            #fig = useful_functions.setup_figure(width='2col',
            #                                    height=float(nrows)/ncols)
            fig = pl.figure(figsize=(width, height))
            #fig = useful_functions.setup_figure(width='2col',
            #                                    height=float(n_rows)/n_cols)

            panels = [fig.add_subplot(n_rows, n_cols, i)
                      for i in xrange(1, n_wells+2)]

            for panel, well in zip(panels, wells):

                ind = df['well'] == well

                sc = panel.scatter(df[param_combination[0]][ind],
                                   df[param_combination[1]][ind],
                                   c=df[gof_col][ind],
                                   edgecolor='None',
                                   s=60,
                                   vmin=0.0, vmax=1.0,
                                   cmap=matplotlib.cm.jet_r)

                xlabel = r'%s' % param_combination[0].replace('_', ' ')
                ylabel = r'%s' % param_combination[1].replace('_', ' ')

                if unit_combination[0] is not None:
                    xlabel += ' (%s)' % unit_combination[0]
                if unit_combination[1] is not None:
                    ylabel += ' (%s)' % unit_combination[1]

                panel.set_xlabel(xlabel)
                panel.set_ylabel(ylabel)
                panel.grid()

                panel.set_xlim(df[param_combination[0]][ind].min(),
                               df[param_combination[0]][ind].max())
                panel.set_ylim(df[param_combination[1]][ind].min(),
                               df[param_combination[1]][ind].max())

                panel.set_title(well)

            fig.tight_layout()

            #fig.subplots_adjust(bottom=0.12)

            print 'drawing colorbar'

            panels[-1].set_xticks([])
            panels[-1].set_yticks([])
            panels[-1].axis('off')

            ax_dimension = panels[-1].get_position().get_points()
            xmin, ymin = ax_dimension[0]
            xmax, ymax = ax_dimension[1]
            ydist = 0.01
            xdist = xmax - xmin
            cax = fig.add_axes((xmin, (ymin + ymax)/2.0, xdist, ydist))
            #cb = pl.colorbar(cf, ax = axp, cax = axl, orientation = 'horizontal', format = '%0.1f')


            #cax = fig.add_axes([0.3, 0.06, 0.4, 0.01])
            cb = fig.colorbar(sc, cax=cax, orientation='horizontal')
            cb_label = gof_col.replace('_', ' ') + '\n(0=bad fit, 1=perfect fit)'
            cb.set_label(cb_label)
            cb.set_ticks(cax.get_xticks()[::2])

            fn_adj = '_' + '_vs_'.join(param_combination) + '_%s' % gof_col + '.png'
            fn = model_result_fn.split('.csv')[0] + fn_adj

            print 'saving %s' % fn

            fig.savefig(fn)

            pl.clf()

print 'done'