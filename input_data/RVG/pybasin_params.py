"""
edit this file to change the model parameters for PyBasin
"""


import numpy as np

#######################################################################
#  model parameters for the 1D burial & temperature history model
#######################################################################

print '-' * 10
print 'RVG input data'
print '-' * 10

# location of input data .csv files
input_dir = 'input_data/RVG'
output_dir = 'model_output/RVG'

##############################################
# select wells
##############################################
#wells =[   'ALM-01', 'AND-06', 'BKZ-01', 'BRAK-01', 'HSW-01',
#           'KDK-01', 'LOZ-01', 'NDW-01', 'SMG-01', 'SPC-01',
#           'VEH-01', 'WWK-01', 'WWN-01']

wells = ['KDK-01']

###########################################################################
# set model scenarios: multple model runs with varying values of exhumation
# or basal heat flow
###########################################################################
# strat period for which to change the basal heat flow 
# must match a period of heatflow_periods further down this file
basal_heat_flow_scenario_period = ['all']
# basal heat flow scenarios to try, heat flow in mW m^-2
basal_heat_flow_scenarios = np.array([66.0])

# strat period for which to change exhumation 
# must exhumation_phase_id further down this file
exhumation_scenarios_period = ['late_cretaceous_unc']
# exhumation (m) for each model scenario
exhumation_scenarios = np.array([0.0, 250.0, 750.0, 1000.0, 1250.0])


##############################################
# sediment provenance parameters, used for AFT
##############################################
# calibrate provenance ages
calibrateProvenanceScenarios = False
# provenance age
# set to 0 for a provenance age that is equal to the stratigraphic age of the sample
prov_ages_start = [430.0, 400.0, 328.0, 305.0]
# time that the sample reaches the surface:
prov_ages_end = [418.0, 10.0, 312.0, 10.0]

######################
# heat flow parameters
######################
# heatflow_periods: first 1,2,3 or more letters of stratigraphic period to 
# set the basal heat flow for. use heatflow_periods = 'all' to set a default 
# value for all strat. periods
heatflow_ages = np.array([0, 260.0, 305, 312])
# heatflow_history: heat flow in W/m^2
heatflow_history = np.array([65.0, 65.0, 130.0, 130.0]) * 1e-3

# optimize heat flow:
optimize_heatflow = False

# max size of heatflow timestep (in yrs)
max_hf_timestep = 10000.0

# resample timesteps for AFT calculation
#
resample_AFT_timesteps = 10

#############################
# exhumation scenarios
#############################
# name of exhumation phase
exhumation_phase_ids = ['late_cretaceous_unc'] 
# start (Ma)
exhumation_period_starts = np.array([85.8])
# end of exhumation phase (Ma)
exhumation_period_ends = np.array([70.0])
# exhumed thickness
exhumed_thicknesses = np.array([500.0])

# determine last deposited unit before unconformity:
exhumed_strat_units = ['SLDNA']

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
max_thickness = 200.0


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
change_thickness_value = [[300, 600, 500, 0], [210, 210, 0]]

############################################
# Apatite fission track model params:
############################################
# use C-axis correction for apatite fission track lengths
use_caxis_correction = False

# parameters for annealing characteristics of apatite grains
# options for kinetic params: 
# 'Clwt' : Chloride wt fractions 
# 'Dpar' : Dpar / etch pit size 
annealing_kinetic_param = 'Clwt'
# end member values for kinetic parameters (if no value given in input data)
annealing_kinetics_values = np.array([0.0001, 0.02])

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
