"""
set up parallel model runs of pybasin

limited to running a single borehole / outcrop location per core for now

modified version, instead of running fixed batches this starts up a new run
everytime one model run is done, after initially starting up a first batch
"""

__author__ = 'elco'

import os
import subprocess
import pdb
import time
import numpy as np
import pandas as pd
from subprocess import Popen

log_folder = 'log'

if os.path.exists(log_folder) is False:
    os.makedirs(log_folder)
    print 'created folder %s to store model output' % log_folder

wait_time = 0.5

ncores = 30

# Rigi dataset: 27 outcrop samples, 2 boreholes, plus Entlebuch for reference:
#locations = ['Huenenberg',
#             'MRP025', 'MRP170', 'MRP172', 'MRP174',
#             'RH10', 'RH12', 'RH15', 'RH17', 'RH20', 'RH23', 'RH30a', 'RH30b',
#             'RH35', 'RH40', 'RH45', 'RH50', 'RH60c', 'RH65', 'RH70',
#             'RV05', 'RV10b', 'RV15', 'RV20', 'RV25',
#             'RV30a', 'RV30b', 'RV30c',
#             'Weggis', 'Entlebuch']

#locations = ['RH17', 'RV05', 'RV10b', 'RV30a', 'RV30b']
#locations = ['RH30a', 'RH30b']

df_loc = pd.read_csv('parallel_run_wells.csv')
locations = list(np.unique(df_loc['well']))

processes = []
fouts = []
threads_occupied = 0
done = 0

for j, inp_file in enumerate(locations):

    print '-' * 20
    print 'at location %s, %i of %i' % (inp_file, j+1, len(locations))
    print '%i processes running' % len(processes)
    print '-' * 20

    if threads_occupied < ncores:
        inp_file_root = inp_file.split('.')[0]
        inp_file_root_short = inp_file_root[3:]

        outfile = os.path.join(log_folder, inp_file_root + '_runlog.txt')

        fout = open(outfile, 'w')
        fouts.append(fout)

        args = ['python', 'pybasin.py', inp_file]

        print '=' * 10
        print 'starting %s' % ' '.join(args)
        print '=' * 10

        p = Popen(args, stdout=fout, stderr=fout, shell=False)
        processes.append(p)

        threads_occupied += 1
        print '=' * 10
        print 'number of threads running = %i' % threads_occupied
        print 'ok, waiting %0.1f sec' % wait_time
        print '=' * 10
        time.sleep(wait_time)

    if len(processes) >= ncores:

        print '=' * 10
        print 'number of processes = %i, waiting for closure of next process' \
              % len(processes)
        print '=' * 10
        processes[0].wait()
        done += 1
        threads_occupied -= 1

        processes.pop(0)
        fouts[0].close()
        fouts.pop(0)

        print '=' * 10
        print 'another process done, %i open' % (len(processes))
        print '=' * 10

while len(processes) > 0:
    print 'number of processes = %i, waiting for closure of next process' \
          % len(processes)
    processes[0].wait()
    processes.pop(0)
    done += 1
    threads_occupied -= 1

    print 'another process done, %i open' % (len(processes))

# close output files
if len(fouts) > 0:
    for fout in fouts:
        fout.close()

print 'done'