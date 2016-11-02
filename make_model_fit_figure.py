"""


"""

__author__ = 'elco'

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import lib.aft_benchmark_lib as abl
import useful_functions


annealing_eq = 'FC'

micro_sign = unichr(181)

# values for AFT annealing parameters, Ketcham (2007)
alpha_old = 0.04672
C0_old = 0.39528
C1_old = 0.01073
C2_old = -65.12969
C3_old = -7.91715

# new calibrated values for annealing parameters, 1 aug 2016
#C0 = 1.75813733e-01
#C1 = 1.21848358e-02
#C2 = -7.06831374e+01
#C3 = -7.80206753e+00
#alpha = 7.78903711e-02

#C0 = 8.66867530e-02
#C1 = 1.15792194e-02
#C2 = -3.25983518e+01
#C3 = -7.09923628e+00
#alpha = 1.43954146e-01

#C0 = 1.94320233e-01
#C1 = 1.22140118e-02
#C2 = -7.05018354e+01
#C3 = -7.82251091e+00
#alpha = 7.26880362e-02

#C0, C1, C2, C3, alpha =
C0 = 2.74620936e-01
C1 = 1.15666643e-02
C2 = -6.94099931e+01
C3 = -7.85643253e+00
alpha = 5.99629016e-02

params_old = [C0_old, C1_old, C2_old, C3_old, alpha_old]
params_new = [C0, C1, C2, C3, alpha]

###################################################

fn = 'benchmark_data/Ketcham2007_appendix_mod_v2_with_kinetic_data.csv'
df = pd.read_csv(fn)

# convert dt to seconds
df['time'] = df['time'] * 60.0 * 60.0
# and T to Kelvin
df['temp'] = df['temp'] + 273.15

# get rid of samples without observed tracks
#ind =   (df['N']>0) & (df['Lm']!=-99999) &\
#        (df['Lc']>0) & (df['dt']>0) & (df['L0']>0)

df = df[df['sigma_Lm'].notnull()]

# remove Cf irradiated samples
df = df[df['Cf']=='n']

# remove samples that have no kinetic data or 0 duration annealing exps
#df = df[(df[dfk.columns[0]].notnull()) & (df['time'] > 0) & (df['sigma_Lm'] > 0)]
# -> keep 0 track lengths, this is also important info....
df = df[(df['rmr0'].notnull()) & (df['time'] > 0)]

# remove samples without initial track length info
# try to find where initial track length came from, not in 2007 xls appendix
# copied from paper???
df = df.loc[df['L0'].notnull()]

# remove samples without any c-axis projected length info
df = df.loc[df['Lc_fit'].notnull()]

observed_lengths = df['Lc_fit'].values
observed_lengths_SE = df['SE_Lc_fit']


##########################

modeled_ln_old, me_old, mae_old, mswd_old = \
    abl.calculate_objective_function_AFT_exp(params_old,
                                             annealing_eq,
                                             df['time'].values,
                                             df['temp'].values,
                                             df['L0'].values,
                                             df['rmr0'].values,
                                             df['kappa'].values,
                                             observed_lengths,
                                             observed_lengths_SE,
                                             False)

modeled_ln_new, me_old, mae_old, mswd_new = \
    abl.calculate_objective_function_AFT_exp(params_new,
                                             annealing_eq,
                                             df['time'].values,
                                             df['temp'].values,
                                             df['L0'].values,
                                             df['rmr0'].values,
                                             df['kappa'].values,
                                             observed_lengths,
                                             observed_lengths_SE,
                                             False)

#############################################
golden_ratio = (1.0 + np.sqrt(5))/2.0

fig, panels = pl.subplots(1, 2, figsize=(8, 4), sharey=True)

panels[0].scatter(observed_lengths, modeled_ln_old,
                  color='gray', edgecolor='black', zorder=2, alpha=0.5)

panels[1].scatter(observed_lengths, modeled_ln_new,
                  color='gray', edgecolor='black', zorder=2, alpha=0.5)

tekst_old = 'Ketcham et al. (2007) model\nN=%i\nmswd=%0.1f' \
            % (len(modeled_ln_old), mswd_old)
tekst_new = 'recalibrated model\nN=%i\nmswd=%0.1f' \
            % (len(modeled_ln_new), mswd_new)

panels[0].text(0.03, 0.80, tekst_old, ha='left',
               transform=panels[0].transAxes)
panels[1].text(0.03, 0.80, tekst_new, ha='left',
               transform=panels[1].transAxes)

for panel in panels:
    panel.plot([0, 20], [0, 20], ls='--', color='black', zorder=0)

for panel in panels:
    useful_functions.simpleaxis(panel)
    panel.yaxis.grid(True)
    panel.xaxis.grid(False)

    panel.set_xlim(-1, 19)
    panel.set_ylim(-1, 19)

    panel.set_xlabel('Observed track length (%sm)' % micro_sign)

panels[0].set_ylabel('Modeled track length (%sm)' % micro_sign)

useful_functions.make_subplot_titles(panels)

fig.tight_layout()

print 'saving figures in directory /fig'
fig.savefig('fig/observed_vs_modeled_ln.png', dpi=200)
fig.savefig('fig/observed_vs_modeled_ln.pdf')

print 'done'