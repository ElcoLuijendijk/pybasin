"""
edit this file to change the model parameters for PyBasin
"""


import numpy as np

#######################################################################
#  model parameters for the 1D burial & temperature history model
#######################################################################

print '-' * 10
print 'example input dataset from the Molasse basin'
print '-' * 10

# location of input data .csv files
input_dir = 'input_data/example_dataset'
output_dir = 'model_output/example_dataset'
datafile_output_dir = 'model_output/example_dataset'

# option to save model run data (approx 10-20 MB per model run)
save_model_run_data = False

# option to generate 1 figure for each model run:
make_model_data_fig = True

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
prov_ages_start = [70.0, 70.0, 30.0, 30.0]
# time that the sample reaches the surface:
prov_ages_end = [68.0, 1.0, 28.0, 1.0]

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
max_hf_timestep = 10000.0

# resample timesteps for AFT calculation
#
resample_AFT_timesteps = 10

#############################
# exhumation scenarios
#############################
# name of exhumation phase
exhumation_phase_ids = ['pre-molasse_exhumation',
                        'molasse_exhumation']
# start (Ma)
exhumation_period_starts = np.array([80, 12.0])
# end of exhumation phase (Ma)
exhumation_period_ends = np.array([40.0, 2.0])
# exhumed thickness
exhumed_thicknesses = np.array([3000.0, 1000.0])

# determine last deposited unit before unfconformity:
# list of 'normal' stratigraphic sequence and thicknesses
# exhumation will start at the lowest missing unit
# Molasse strat units, following Kemp et al. (1999)
exhumed_strat_units = [['Kimm'],
                       ['UMM', 'USM-I', 'USM-II',
                       'OMM', 'OSM', 'thrust_sheet_1',
                       'thrust_sheet_2']]

# original thicknesses of strat units
original_thicknesses = [[3000.0],
                        [600.0, 2175.0, 1125.0,
                        525.0, 4000.0, 4000.0, 4000.0]]

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
