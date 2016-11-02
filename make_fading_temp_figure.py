"""
make figures of fading temperature vs fading time for old and new annealing algorithms

"""

__author__ = 'elco'

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

import useful_functions
import lib.AFTannealingLib as AFTannealingLib


def calculate_fading_temperature(alpha, C0, C1, C2, C3, rmr0=None, kappa=None,
                                 rcf=0.527,
                                 fading_time_My=100.0):

    """
    calculate fading temperature = approx. of temperature at which all tracks anneal after a certain time
    first mentioned in Ketcham (1999), eq. not supplied, but assumed equal to annealing eq. with a r of 0.41

    quote: "In the present study, we approximate TF as the temperature required to anneal
    the reduced mean track length to a value of 0.41, the lowest value observed in any published
    annealing experiment to date in which more than 10 tracks were measured."

    rmr0 and kappa for fluorapatite = 0.83 and 0.21

    value of rcf given in FTindex source code: rcf = 0.55
    however, using density vs mean -caxis length conversion given by Ketcham
     rho = 9.205 rc^2 -9.157 rc + 2.269
    zero track density should occur for rc values of 0.527



    Parameters:
    -----------
    alpha : float
        empirical parameter fission track annealing equation
    C0 : float
        empirical parameter fission track annealing equation
    C1 : float
        empirical parameter fission track annealing equation
    C2 : float
        empirical parameter fission track annealing equation
    C3 : float
        empirical parameter fission track annealing equation
    rmr0 : float
        empirical kinetic parameter that governs the relative resistance to
        fission track annealing
    kappa : float
        empirical kinetic parameter that governs the relative resistance to
        fission track annealing
    rmf : float
        reduced track length at which fission track density is assumed to be 0
        (um)
    fading_time : float
        fading time (My)
    """

    # eq. given to convert normal ot c-axis projected length in fig. 8 of
    # Ketcham (1999)
    #a = -1.499
    #b = 4.150
    #c = -1.656
    #x = rmf - c
    #discr = np.sqrt(4 * a * x + b**2)
    #rcf = (discr - b) / (2 * a)

    # alternative equation, provided by Ketcham (2005) Ch. 22:
    #rcf = (rmf + 0.4017) / 1.396

    # alternative, rcf = 0.55, given in FTindex source code
    #rcf = 0.55

    # convert fading time from My to sec
    fading_time = fading_time_My * 1e6 * 365.25 * 24 * 60 * 60

    if rmr0 is not None:
        rcf_mod = AFTannealingLib.kinetic_modifier_reduced_lengths_inverse(
            rcf, rmr0, kappa)
    else:
        rcf_mod = rcf

    # fading value for g
    gf = (1.0 / rcf_mod - 1.0) ** alpha

    # fading temperature:
    Tfa = (np.log(fading_time) - C2) * C1 / (gf-C0) + C3
    Tf = 1.0 / np.exp(Tfa)

    g_check = C0 + C1 * ((np.log(fading_time) - C2) / (np.log(1.0 / Tf) - C3))
    f1 = g_check**(1./alpha)
    rc_check = 1.0 / (f1 + 1.0)

    if np.max(np.abs(rc_check - rcf_mod)) > 1e-3:
        msg = 'error, something wrong with eq., ' \
              'rcf_mod = %0.2e, rc_check = % 0.2e' % (rcf_mod, rc_check)
        raise ValueError(msg)

    return Tf


K = 273.15
degree_symbol = unichr(176)

# values for AFT annealing parameters FC model, Ketcham (2007)
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

C0 = 2.74620936e-01
C1 = 1.15666643e-02
C2 = -6.94099931e+01
C3 = -7.85643253e+00
alpha = 5.99629016e-02

# calculate fading temperature. ie, the temperature at which a fission track
# population is completely annealed after a given time of isothermal heating
fading_times = np.arange(0.1, 100.1, 0.1)

# kinetic parameters fluorapatite UNK, durango and B2
# see Ketcham (2007), Table 5
apatite_names = ['UNK', 'Durango', 'B2']
apatite_colors = ['blue', 'gray', 'black']
kinetic_params = [[0.83220, 0.23596], [0.79757, 0.22127], [0.0, 0.97933]]

fading_temps_old = np.zeros((len(kinetic_params), len(fading_times)))
fading_temps_new = np.zeros((len(kinetic_params), len(fading_times)))


for i, kinetic_param in enumerate(kinetic_params):

    fading_temps_old[i] = calculate_fading_temperature(alpha_old, C0_old, C1_old, C2_old, C3_old,
                                                       rmr0=kinetic_param[0],
                                                       kappa=kinetic_param[1],
                                                       fading_time_My=fading_times) - K

    fading_temps_new[i] = calculate_fading_temperature(alpha, C0, C1, C2, C3,
                                                       rmr0=kinetic_param[0],
                                                       kappa=kinetic_param[1],
                                                       fading_time_My=fading_times) - K

for i, apatite_name in enumerate(apatite_names):
    print 'mean, min, max difference in fading temp for apatite %s' % apatite_name
    print (np.mean(fading_temps_old[i] - fading_temps_new[i]),
           np.min(fading_temps_old[i] - fading_temps_new[i]),
           np.max(fading_temps_old[i] - fading_temps_new[i]))



###########################################################################
# figure 1: fading temperature vs time for the old and new annealing models
# for three different apatite
###########################################################################
golden_ratio = (1.0 + np.sqrt(5))/2.0

fig, panel = pl.subplots(1, 1, figsize=(6, 6/golden_ratio))

legs_new = []

for fading_temp_old, fading_temp_new, apatite_name, apatite_color \
        in zip(fading_temps_old, fading_temps_new, apatite_names, apatite_colors):

    leg_old, = panel.plot(fading_times, fading_temp_old,
                          color=apatite_color, ls='--', lw=0.5)
    leg_new, = panel.plot(fading_times, fading_temp_new,
                          color=apatite_color, ls='-', lw=0.5)
    legs_new.append(leg_new)

panel.set_ylim(0, panel.get_ylim()[-1] * 1.10)

useful_functions.simpleaxis(panel)

panel.yaxis.grid(True)
panel.xaxis.grid(False)

panel.set_xlabel('Fading time (My)')
panel.set_ylabel('Fading temperature (%sC)' % degree_symbol)

labels = ['Ketcham et al. (2007) annealing model',
          'recalibrated annealing model']
labels += ['apatite %s' % apatite_name for apatite_name in apatite_names]
legs = [leg_old, leg_new] + legs_new
panel.legend(legs, labels, loc='upper right', fontsize='x-small',
             frameon=False, handlelength=3)

fig.tight_layout()

print 'saving figures in directory /fig'
fig.savefig('fig/fading_temperatures.png', dpi=200)
fig.savefig('fig/fading_temperatures.pdf')


print 'done'
