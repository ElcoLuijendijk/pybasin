import pickle
import argparse
import os
import sys
import pdb

import lib.pybasin_figures as pf

parser = argparse.ArgumentParser(description='Make figures of PyBasin model runs')

parser.add_argument('model_file_or_directory', metavar='directory or file', default=None, nargs='?',
                    help='input dataset or directory containing input datasets for PyBasin')

parser.add_argument('-f', '--figure_type', dest='fig_type', default='png',
                    help='figure type, for instance jpg, png, pdf or svg')

parser.print_help()

print '\n'

args = parser.parse_args()

# check if script dir in python path
scriptdir = os.path.realpath(sys.path[0])
if scriptdir not in sys.path:
    sys.path.append(scriptdir)

if args.model_file_or_directory is not None:
    model_file_or_directory = os.path.join(scriptdir, args.model_file_or_directory)

    if os.path.isdir(model_file_or_directory) == False:
        model_file = os.path.join(scriptdir, args.model_file_or_directory)
        model_dir = None
    else:
        model_file = None
        model_dir = model_file_or_directory

else:
    print 'enter a directory name:'
    model_dir = raw_input()

    model_files = os.listdir(model_dir)
    model_files = [f for f in model_files if f[-4:] == '.pck']

    for i, mf in enumerate(model_files):
        print i, '\t', mf

    print 'enter the number of the file/model run you want to make a figure of:'

    mfi = int(raw_input())

    model_file = os.path.join(model_dir, model_files[mfi])


fin = open(model_file)
model_run_data_fig = pickle.load(fin)
fin.close()

fig = pf.model_vs_data_figure(model_run_data_fig)

fn_out = model_file[:-4] + '_figure.%s' % args.fig_type

print 'saving figure as %s' % fn_out
fig.savefig(fn_out, dpi=300)

print 'done'





