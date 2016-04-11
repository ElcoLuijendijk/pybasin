"""
edit this file to change the model parameters for PyBasin
"""


import numpy as np

#######################################################################
#  model parameters for the 1D burial & temperature history model
#######################################################################

print '-' * 10
print 'Salinity diffusion input data'
print '-' * 10

# location of input data .csv files
input_dir = 'input_data/salinity_grid'
output_dir = 'model_output/salinity_grid'
datafile_output_dir = '../../heavy_data/pybasin_salinity'

# option to calculate apatite fission track data
simulate_AFT = False
simulate_AHe = False
simulate_VR = False
simulate_salinity = True

# option to save model run data (approx 10-20 MB per model run)
save_model_run_data = False

# option to generate 1 figure for each model run:
make_model_data_fig = False

# black and white figure
model_data_fig_bw = False

# varaible to show color contours for in burial history panel
# choose either 'temperature' or 'salinity'
# to show evolution of temperature or salinity over time
contour_variable = 'salinity'

# type of figure file to save (pdf, png or jpg):
fig_adj = 'png'

#################
# goodness of fit
################
# weights for calculating overall goodness of fit from the gof statistic for
# temperature, vitrinite reflectance and apatite fission track age data
gof_weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]

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

##############################################
# sediment provenance parameters, used for AFT
##############################################
# calibrate provenance ages
calibrateProvenanceScenarios = False
# provenance age
# set to 0 for a provenance age that is equal to the stratigraphic
# age of the sample
prov_ages_start = [70.0, 70.0, 25.0, 20.0]
# time that the sample reaches the surface:
prov_ages_end = [68.0, 1.0, 23.1, 1.0]

provenance_time_nt = 100

######################
# heat flow parameters
######################
# heatflow_periods: first 1,2,3 or more letters of stratigraphic period to 
# set the basal heat flow for. use heatflow_periods = 'all' to set a default 
# value for all strat. periods
heatflow_ages = np.array([0, 100.0])
# heatflow_history: heat flow in W/m^2
heatflow_history = np.array([65.0, 65.0]) * 1e-3

# optimize heat flow:
optimize_heatflow = False

# max size of heatflow timestep (in yrs)
max_hf_timestep = 100000.0

# resample timesteps for AFT calculation, number of timesteps that
resample_AFT_timesteps = 10

#############################
# exhumation scenarios
#############################
# name of exhumation phase
exhumation_phase_ids = ['late_cretaceous_exhumation']
# start (Ma)
exhumation_period_starts = np.array([85.8])
# end of exhumation phase (Ma)
exhumation_period_ends = np.array([70.0])
# exhumed thickness
exhumed_thicknesses = np.array([1.0])

# determine last deposited unit before unfconformity:
# list of 'normal' stratigraphic thicknesses
# exhumation will start at the lowest missing unit
# list units from old to young
# nested list, separate list for each exhumation phase
# for example if unit A and B predate exhumation phase 1 and
# unit E,F and G predate exhumation phase 2
# exhumed_strat_units = [['A', 'B'], ['E', 'F', G']]

exhumed_strat_units = [['AT', 'SL', 'KN']]

# thicknesses
original_thicknesses = [[0.0]]


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
annealing_kinetic_param = 'Dpar'
# end member values for kinetic parameters (if no value given in input data)
#annealing_kinetics_values = np.array([1.5, 1.8])
annealing_kinetics_values = np.array([1.2, 2.2])

# size of bins of (simulated) AFT length histogram, default = 0.25 um 
binsize = 0.25    

# empirical coefficients AFT annealing equation
# default values from Ketcham et al. (2007) American Mineralogist
# fanning curvelinear model values in Table 5
alpha = 0.04672
C0 = 0.39528
C1 = 0.01073
C2 = -65.12969
C3 = -7.91715

##################
# (U-Th)/He params
##################
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
# sigma of uncertainty range for VR data, if not specified in input file
vr_unc_sigma = 0.05

###########################
# salinity diffusion params
###########################
tortuosity_factor = -1/3.

Dw = 20.3e-10

fixed_lower_bnd_salinity = 0.30

constant_diffusivity = False

salinity_seawater = 0.035
salinity_freshwater = 0.0001

# read surface salinity bnd condition for each well from a separate .csv file
well_specific_surface_salinity_bnd = True