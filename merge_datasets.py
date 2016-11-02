"""
combine AFT and kinetic parameter datasets

"""

__author__ = 'elco'

import pandas as pd
import numpy as np


###########################
print 'loading AFT lab data'

fn = 'benchmark_data/Ketcham2007_appendix_mod_v2.csv'
df = pd.read_csv(fn)


# load kinetic params
print 'loading kinetic parameters for lab data'

dfk = pd.read_csv('benchmark_data/Ketcham2007_kinetic_params.csv', index_col=0)

print 'combining kinetic parameters and annealing data table'
# populate table with kinetic params
for col in dfk.columns:
    df[col] = np.nan

for apatite in dfk.index:
    for col in dfk.columns:
        df[col][df['apatite'] == apatite] = dfk.ix[apatite, col]

fn_out = 'benchmark_data/Ketcham2007_appendix_mod_v2_with_kinetic_data.csv'
print 'saving results to %s' % fn_out
df.to_csv(fn_out, index=False)

print 'done'