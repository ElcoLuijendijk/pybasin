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
input_dir = 'input_data/salinity'
output_dir = 'model_output/salinity'
datafile_output_dir = '../../heavy_data/pybasin_salinity'

# option to calculate apatite fission track data
simulate_AFT = False
simulate_AHe = False
simulate_VR = False
simulate_salinity = True

# option to save model run data (approx 10-20 MB per model run)
save_model_run_data = False

# option to generate 1 figure for each model run:
make_model_data_fig = True

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
exhumed_thicknesses = np.array([500.0])

# determine last deposited unit before unconformity:
exhumed_strat_units = ['SLDNA']

# determine last deposited unit before unfconformity:
# list of 'normal' stratigrpahic thicknesses
# exhumation will start at the lowest missing unit
# Molasse strat units, following Kemp et al. (1999)
exhumed_strat_units = [['ATAL']]

# thicknesses
original_thicknesses = [[500.0]]

# or set pre-exhumation thickness of particular unit, if known
# exhumation will then be calculated to match the present-day thickness
# of this unit in each well
# set to None if you do not use this feature
#pre_exhumation_units = [None, 'AT', None, None]
#pre_exhumation_thicknesses = np.array([0.0, 1200.0, 0.0, 0.0])

###########################################
# max thickness of strat units
# units that exceed this are subdivided
# to keep the modeled temperatures accurate
###########################################
max_thickness = 100.0

#################################################################
# change thickness of strat units
# use this to simulate fault movement
#################################################################
# well ids of wells that are drilled across faults:
change_thickness_wells = ['HSW-01', 'ALM-01',]
# stratigraphic units to change thickness of:
change_thickness_units = ['ATWD', 'ATBR2']
# time steps (in strat units) for thickness changes
# if a sequence of time steps is given the changes are distributed evenly
change_thickness_timing = [['ATWD', 'ATBR2', 'ATBRU', 'SL'],
                          ['ATBR2', 'NL', 'NU']]
# change in thickness (m) for each time step:
# 0 means that the unit is at the present-day thickness
# positive numbers mean that the thickness of the unit increases over time
# negative that the thickness decreases
change_thickness_value = [[300, 600, 500, 0],
                          [210, 210, 0]]

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