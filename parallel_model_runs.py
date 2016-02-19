"""
set up parallel model runs of pybasin

limited to running a single borehole / outcrop location per core for now

feb 2016, Elco Luijendijk
"""

__author__ = 'elco'

import os
import subprocess
import time
import numpy as np
from subprocess import Popen

log_folder = 'log'

if os.path.exists(log_folder) is False:
    os.makedirs(log_folder)
    print 'created folder %s to store model output' % log_folder

wait_time = 2.0

ncores = 38

locations = ['Huenenberg',
             'MRP025', 'MRP170', 'MRP172', 'MRP174',
             'RH10', 'RH15', 'RH25', 'RH30a', 'RH30b',
             'RH35', 'RH40', 'RH45', 'RH50', 'RH60c',
             'RH65', 'RH70', 'RV15', 'RV20', 'RV25',
             'RV30a', 'RV30b', 'RV30c', 'Weggis']

locations = ['Huenenberg',
             'MRP025', 'MRP170', 'MRP172', 'MRP174',
             'RH10', 'RH12', 'RH15', 'RH17', 'RH20', 'RH23', 'RH30a', 'RH30b',
             'RH35', 'RH40', 'RH45', 'RH50', 'RH60c', 'RH65', 'RH70',
             'RV15', 'RV20', 'RV25', 'RV30a', 'RV30b', 'RV30c', 'Weggis',
             'RV05', 'RV10b', 'RV15', 'RV25', 'RV30c']

# new locations, 19 feb 2016
locations = ['RH12', 'RH15', 'RH17', 'RH20', 'RH23', 'RV05', 'RV10b', 'RV15']

n_batches = int(np.ceil(len(locations) / float(ncores)))

for batch in range(n_batches):

    start = batch * ncores
    end = start + ncores

    if end > len(locations):
        end = -1
        locations_batch = locations[start:]
    else:
        locations_batch = locations[start:end]

    print 'running batch ', locations[start:end]

    processes = []
    fouts = []

    for inp_file in locations_batch:

        inp_file_root = inp_file.split('.')[0]
        inp_file_root_short = inp_file_root[3:]

        outfile = os.path.join(log_folder, inp_file_root + '_runlog.txt')

        fout = open(outfile, 'w')
        fouts.append(fout)

        args = ['python', 'pybasin.py', inp_file]

        p = Popen(args, stdout=fout, shell=False)
        processes.append(p)

        print 'ok, waiting %0.1f sec' % wait_time

        time.sleep(wait_time)

    print 'waiting for completion of processes'
    done = 0
    for i, p in enumerate(processes):
        p.wait()
        done += 1
        print 'another process done %i to go' % (ncores - done)

    # close output files
    for fout in fouts:
        fout.close()

print 'done'