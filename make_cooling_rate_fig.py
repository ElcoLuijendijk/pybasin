"""
make a figure of cooling rate vs AFT age and AHe age for the Ketcham(2007) and the new AFT annealing models
He diffusivity is calculated using the fission track density following the RDAAM model (Flowers et al 2009)

"""

import numpy as np
import matplotlib.pyplot as pl
import lib.AFTannealingLib as AFTannealingLib
import lib.helium_diffusion_models as he

import useful_functions

__author__ = 'elcopone'

# values for AFT annealing parameters FC model, Ketcham (2007)
alpha_old = 0.04672
C0_old = 0.39528
C1_old = 0.01073
C2_old = -65.12969
C3_old = -7.91715

# new calibrated values for annealing parameters, 17 aug 2016
C0 = 2.74620936e-01
C1 = 1.15666643e-02
C2 = -6.94099931e+01
C3 = -7.85643253e+00
alpha = 5.99629016e-02

# AHe parameters:
radius = 65e-6
# U concentration
# Durango: U=8.98 ppm, Th=161.3 ppm
# Wolf (1996) table 3:
U = 8.98e-6
Th = 161.3e-6

Myr = 1e6 * 365.25 * 24 * 60 * 60
Kelvin = 273.15

# option to use faster fortran algorithm for calculated AFT and Helium
# diffusivity with the RDAAM model
use_fortran_algorithm = False

# cooling rates in degr C/My
#cooling_rates = np.concatenate((np.arange(0.1, 1.0, 0.1), 
#                                np.arange(1.0, 11.0, 1.0)))

#cooling_rates = [0.1, 1.0, 10.0, 100.0]
c_int = 0.1
cooling_rates = 10**(np.arange(-1, 2 + c_int, c_int))
AFT_ages_old = np.zeros_like(cooling_rates)
AFT_ages_new = np.zeros_like(cooling_rates)

AHe_ages_old = np.zeros_like(cooling_rates)
AHe_ages_new = np.zeros_like(cooling_rates)

#
rmr0 = 0.83220
kappa = 0.23596

#
for i, cooling_rate in enumerate(cooling_rates):
    
    
    #
    nsteps = 1000
    time_150deg = 140.0/cooling_rate
    #time = np.arange(0, time_150deg + 0.1, 0.1)
    time = np.linspace(0, time_150deg, nsteps)

    # linear cooling path from 150 C to 10 C
    #temperature = 150.0 - time * cooling_rate
    temperature = np.linspace(150.0, 10, nsteps)# - time * cooling_rate

    track_length_pdf, aft_age_myr_old, l_mean, l_mean_std, \
        rm, rc_max, rho_age, dt = AFTannealingLib.simulate_AFT_annealing(
            time, temperature, rmr0,
            kinetic_parameter='rmr0', kappa=kappa,
            annealing_eq='FC',
            alpha=alpha_old, C0=C0_old, C1=C1_old, C2=C2_old, C3=C3_old,
            use_fortran_algorithm=use_fortran_algorithm)

    track_length_pdf, aft_age_myr_new, l_mean, l_mean_std, \
        rm, rc_max, rho_age, dt = AFTannealingLib.simulate_AFT_annealing(
            time, temperature, rmr0,
            kinetic_parameter='rmr0', kappa=kappa,
            annealing_eq='FC',
            alpha=alpha, C0=C0, C1=C1, C2=C2, C3=C3,
            use_fortran_algorithm=use_fortran_algorithm)
    
    
    time_sec = time * Myr
    temperature_K = temperature + Kelvin
    ahe_age_old = he.calculate_he_age_meesters_dunai_2002(time_sec, temperature_K, 
                                                       radius, U, Th,
                                                       alpha=alpha_old, 
                                                       C0=C0_old, C1=C1_old, 
                                                       C2=C2_old, C3=C3_old,
                                                       use_fortran_algorithm=
                                                       use_fortran_algorithm)

    ahe_age_new = he.calculate_he_age_meesters_dunai_2002(time_sec, temperature_K, 
                                                       radius, U, Th,
                                                       alpha=alpha, C0=C0, 
                                                       C1=C1, C2=C2, C3=C3,
                                                       use_fortran_algorithm=
                                                       use_fortran_algorithm)
                                         
    AFT_ages_old[i] = aft_age_myr_old
    AFT_ages_new[i] = aft_age_myr_new
    
    AHe_ages_old[i] = ahe_age_old[-1] / Myr
    AHe_ages_new[i] = ahe_age_new[-1] / Myr


ahe_ages_int = np.linspace(1, int(AFT_ages_old.max()), int(AFT_ages_old.max()))
cooling_rates_old = np.interp(ahe_ages_int, AFT_ages_old[::-1], cooling_rates)
cooling_rates_new = np.interp(ahe_ages_int, AFT_ages_new[::-1], cooling_rates)

################
degree_symbol = unichr(176)
golden_ratio = (1.0 + np.sqrt(5))/2.0

fig, panels = pl.subplots(1, 2, figsize=(8, 8/golden_ratio))

labels = ['AFT age (My)', 'AHe age (My)']
ages_new = [AFT_ages_new, AHe_ages_new]
ages_old = [AFT_ages_old, AHe_ages_old]

panelsr = [panel.twinx() for panel in panels]

rcolor = 'blue'

for panel, panelr, age_new, age_old, label in zip(panels, panelsr, 
                                                  ages_new, ages_old,
                                                  labels):

    #panelr = panel.twinx()
    relative_diff = (age_new / age_old - 1.0) * 100.0
    print 'diff = %0.1f, %0.1f' % (relative_diff.min(), relative_diff.max())
    leg_diff, = panelr.plot(cooling_rates, relative_diff, lw=1.5,
                           ls='-', color=rcolor)

    leg_old, = panel.plot(cooling_rates, age_old, 
                         ls='--', lw=1.0, color='gray')
    leg_new, = panel.plot(cooling_rates, age_new, 
                         ls='-', lw=1.0, color='black')

    panel.set_yscale('log')
    panel.set_xscale('log')
    
    panel.set_ylim(0.1, 1400)

    useful_functions.simpleaxis(panel)

    panel.yaxis.grid(True)
    panel.xaxis.grid(False)

    panelr.spines['right'].set_color(rcolor)
    panelr.yaxis.label.set_color(rcolor)
    panelr.tick_params(axis='y', colors=rcolor)
    
    panel.spines['top'].set_visible(False)
    panelr.spines['top'].set_visible(False)

    panel.set_xlabel(r'Cooling rate (%sC My$^{-1}$)' % degree_symbol)
    panel.set_ylabel(label)
    panelr.set_ylabel('Age difference new and old model (%)')

    panelr.set_ylim(0, 50)
    
    
useful_functions.make_subplot_titles(panels)

fig.subplots_adjust(bottom=0.23, wspace=0.5, left=0.1, right=0.90, top=0.97)

legs = [leg_old, leg_new, leg_diff]
labels = ['Ketcham(2007) annealing model', 'recalibrated annealing model',
          'relative difference']
fig.legend(legs, labels, loc='lower center', frameon=False, 
           ncol=2, fontsize='small')


#fig.tight_layout()

print 'saving figs'
fig.savefig('fig/cooling_rate_vs_thermochron_ages.png', dpi=300)
fig.savefig('fig/cooling_rate_vs_thermochron_ages.pdf')

print 'done'

