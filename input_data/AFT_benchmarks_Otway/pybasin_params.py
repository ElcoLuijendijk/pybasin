"""
edit this file to change the model parameters for PyBasin
"""


import numpy as np

#######################################################################
#  model parameters for the 1D burial & temperature history model
#######################################################################

print '-' * 10
print 'AFT benchmark input'
print '-' * 10

# location of input data .csv files
input_dir = 'input_data/AFT_benchmarks_Otway'
output_dir = 'model_output/AFT_benchmarks_Otway'
datafile_output_dir = '../../heavy_data/AFT_benchmarks'

# option to calculate apatite fission track data
simulate_AFT = True
simulate_AHe = False
simulate_VR = True
simulate_salinity = False

# option to calculate AHe ages for all nodes rather than just the samples
# in a particular borehole
# is automatically set to True if make_model_data_fig = True
calculate_thermochron_for_all_nodes = False

# option to save model run data (approx 10-20 MB per model run)
save_model_run_data = False

# option to generate 1 figure for each model run:
make_model_data_fig = True

# variable to show color contours for in burial history panel
# choose either 'temperature' or 'salinity'
# to show evolution of temperature or salinity over time
contour_variable = 'temperature'

# type of figure file to save (pdf, png or jpg):
fig_adj = ['pdf', 'png']


###################
# model_calibration
###################

# turn model calibration on or off:
calibrate_model_params = False

# calibration method, see scipy optimize documentation for list of available methods:
# http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.optimize.minimize.html
opt_method = 'Nelder-Mead'

# list the parameters that should be updated by either the automatic
# calibration function or the grid model space search
# chooose any combination of 'exhumation_magnitude', 'exhumation_start',
# 'exhumation_duration', 'basal_heat_flow',
# or fission track annealing params:
# 'AFT_C0', 'AFT_C1', 'AFT_C2', 'AFT_C3', 'AFT_alpha'
#params_to_change = ['exhumation_magnitude',
#                    'exhumation_start',
#                    'exhumation_duration',
#                    'basal_heat_flow']
params_to_change = ['AFT_C0', 'AFT_C1',
                    'AFT_C2', 'AFT_C3',
                    'AFT_alpha']

# initial values for model parameters
start_param_values = [0.39528, 0.01073, -65.12969, -7.91715, 0.04672]
#start_param_values = [2000.0, 10.0, 7.0]

# read initial params from file
load_initial_params = False
initial_params_file = 'initial_param_values.csv'

# min. and max bounds for parameters
# set to None for unconstrained calibration
#param_bounds_min = [0.0, 1.0, 0.5, 40e-3]
#param_bounds_max = [6000.0, 12.0, 11.0, 100e-3]
param_bounds_min = None
param_bounds_max = None

# list of variables to calibrate model to
# choose any combination of 'T', 'VR', 'AFT_age' or 'AHe'
# for temperature, vitrinite reflectance, apatite fission track age and
# apatite (U-Th)/He age, respectively
calibration_target = ['AFT_age']

#################
# goodness of fit
################
# weights for calculating overall goodness of fit from the gof statistic for
# temperature, vitrinite reflectance, apatite fission track age and
# apatite (U-Th)/He data
gof_weights = [1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0]

provenance_time_nt = 100

######################
# heat flow parameters
######################
# heatflow_history: age in My and basal heat flow in W/m^2
heatflow_ages = np.array([0.0, 90.0, 120.0])
heatflow_history = np.array([75.0, 90.0, 75.0]) * 1e-3

# optimize heat flow:
optimize_heatflow = False

# max size of heatflow timestep (in yrs)
max_hf_timestep = 10000.0

# resample timesteps for AFT calculation, number of timesteps that
resample_AFT_timesteps = 10

#############################
# exhumation scenarios
#############################
# name of exhumation phase
exhumation_phase_ids = ['late_miocene_exhumation']
# start (Ma)
exhumation_period_starts = np.array([10.0])
# end of exhumation phase (Ma)
exhumation_period_ends = np.array([0.0])
# exhumed thickness
exhumed_thicknesses = np.array([0.0])

# determine last deposited unit before unfconformity:
# list of 'normal' stratigrpahic thicknesses
# exhumation will start at the lowest missing unit
exhumed_strat_units = [['Gellibrand_Marl', 'Port_Campbell_Lst']]

# thicknesses
original_thicknesses = [[750.0, 400.0]]

# or set pre-exhumation thickness of particular unit, if known
# exhumation will then be calculated to match the present-day thickness
# of this unit in each well
# set to None if you do not use this feature
#pre_exhumation_units = [None, 'AT', None, None]
#pre_exhumation_thicknesses = np.array([0.0, 1200.0, 0.0, 0.0])

# support for two-stage exhumation history, enables fast and slow exhumation segments
# switch for two-stage exhumation
two_stage_exhumation = False
# fraction of total duration of exhumation phase that separates the first and second segment
exhumation_segment_factor = 0.5
# fraction of exhumation that takes place in the first of two segments
exhumation_duration_factor = 0.5

# parameter to automatically reduce exhumation duration if end of
# exhumation is < 0 Ma
correct_exhumation_duration = True

###########################################
# max thickness of strat units
# units that exceed this are subdivided
# to keep the modeled temperatures accurate
###########################################
max_thickness = 100.0

############################################
# Apatite fission track model params:
############################################

#
min_grain_no = 2

# use C-axis correction for apatite fission track lengths
use_caxis_correction = False

# parameters for annealing characteristics of apatite grains
# options for kinetic params: 
# 'Clwt' : Chloride wt fractions 
# 'Dpar' : Dpar / etch pit size 
annealing_kinetic_param = 'Clwt'
# end member values for kinetic parameters (if no value given in input data)
# min and max Cl wt % Flaxmans-1 sample (Green 1989): 0.179 % - 2.08 %
annealing_kinetics_values = np.array([0.00179, 0.0208])

# size of bins of (simulated) AFT length histogram, default = 0.25 um 
binsize = 0.25    

# annealing equation to use
# 'FA' for fanning Arrhenius equation by Laslett (1987)
# 'FC' for fanning curvelinear equation used by Ketcham (1999, 2007)
annealing_equation = 'FC'

# empirical coefficients AFT annealing equation
# default values from Ketcham et al. (2007) American Mineralogist
# fanning curvelinear model values in Table 5
#alpha = 0.04672
#C0 = 0.39528
#C1 = 0.01073
#C2 = -65.12969
#C3 = -7.91715

# calibrated values for Frio pseudo-borehole data by Corrigan (1993):
#alpha = 0.0570171869542
#C0 = 0.465675271105
#C1 = 0.00893230883784
#C2 = -75.4037631023
#C3 = -8.20203665287

# calibrated values Frio dataset (Corrigan 1993)
# + short term annealing experiments, Vrolijk ea (1992) data and
# Fish Canyon Tuff
# cal results 17 aug 2016:
C0 = 2.74620936e-01
C1 = 1.15666643e-02
C2 = -6.94099931e+01
C3 = -7.85643253e+00
alpha = 5.99629016e-02


##################
# (U-Th)/He params
##################

# apatite U-Th/He equations to use
# 'Farley2000' for helium diffusion parameters of Durango apatite
#   acc. to Farley(2000) JGR 105
# 'RDAAM' for he diffusion that depends on radiation damage acc. to
#   Flowers et al. (2009) GCA 73
ahe_method = 'RDAAM'

# decay constants
decay_constant_238U = 4.916e-18
decay_constant_232Th = 1.57e-18
decay_constant_235U = 3.12e-17


###################################################
# compaction
# (see input data for porosity vs depth parameters)
###################################################
# number of iterations for the calculation of decompaction
NcompactionIterations = 5

# max error when decompacting
max_decompaction_error = 0.01

#######
# VR
#######
# default sigma of uncertainty range for VR data,
# if not specified in input file
vr_unc_sigma = 0.05
