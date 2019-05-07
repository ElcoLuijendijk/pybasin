"""
model parameters for the 1D burial & temperature history model

edit this file to change the model parameters for PyBasin
"""

import numpy as np


print(('-' * 10))
print('Example PyBasin dataset from the Roer Valley Graben (Luijendijk et al. 2011)')
print(('-' * 10))


class ModelParameters:

    # location of input data .csv files
    output_dir = 'model_output/example_dataset_1'
    datafile_output_dir = 'model_output/example_dataset_1'

    # names of wells or surface outcrops to include in a single set of model runs:
    wells = ['BKZ-01']

    # option to calculate apatite fission track data
    simulate_AFT = True
    simulate_AHe = False
    simulate_VR = True
    simulate_salinity = False

    # option to calculate AHe ages for all nodes rather than just the samples
    # in a particular borehole
    # is automatically set to True if make_model_data_fig = True
    # note that this increases the amount of computational time quite a bit
    calculate_thermochron_for_all_nodes = True

    # option to save model run data (approx 10-20 MB per model run)
    save_model_run_data = True

    # use stratigraphy input data from stratigraphic maps instead of text files
    # this is still an experimental feature, no guarantee that it actually works. Future updates will make this more
    # user friendly and bug-free (hopefully)
    use_strat_map_input = False

    # save results to a .csv file for each x number of model runs
    csv_save_interval = 1

    ################
    # figure options
    ################
    # option to generate 1 figure for each model run:
    make_model_data_fig = True

    # variable to show color contours for in burial history panel
    # choose either 'temperature' or 'salinity'
    # to show evolution of temperature or salinity over time
    contour_variable = 'temperature'

    # add a stratigraphic column to the figure
    show_strat_column = False

    # option to hide thermochron results
    show_thermochron_data = True

    # type of figure file to save (pdf, png or jpg):
    fig_adj = ['png']

    ###########################################
    # max thickness of strat units
    # units that exceed this are subdivided
    # to keep the modeled temperatures accurate
    ###########################################
    max_thickness = 100.0

    ###################################################
    # compaction
    # (see input data for porosity vs depth parameters)
    ###################################################
    # number of iterations for the calculation of decompaction
    NcompactionIterations = 5

    # max error when decompacting
    max_decompaction_error = 0.01

    #############################
    # exhumation scenarios
    #############################
    # start of exhumation (Ma)
    exhumation_period_starts = np.array([85.8])
    # end of exhumation phase (Ma)
    exhumation_period_ends = np.array([71.0])
    # exhumed thickness (m)
    exhumed_thicknesses = np.array([1000.0])

    # determine last deposited units before unconformity:
    # this should be one list for each exhumation phase, with stratigraphic unit codes ordered from old to young
    # the model will add units starting from the oldest to the youngest, untill the additional thickness needed for
    # erosion is filled
    exhumed_strat_units = [['ATBR3', 'ATBRU', 'ATBRO', 'SLDNA']]

    # initial (pre-erosion) thicknesses:
    original_thicknesses = [[50, 50, 50, 1500.0]]

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

    ######################
    # heat flow parameters
    ######################
    # heatflow_history: heat flow in W/m^2, age in Ma
    heatflow_ages = np.array([0, 260.0, 305, 312])
    heatflow_history = np.array([68.0, 68.0, 130.0, 130.0]) * 1e-3

    # max size of heatflow timestep (in yrs)
    max_hf_timestep = 10000.0

    #################
    # goodness of fit
    ################
    # weights for calculating overall goodness of fit from the gof statistic for
    # temperature, vitrinite reflectance, apatite fission track age and
    # apatite (U-Th)/He data
    gof_weights = [1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0]

    provenance_time_nt = 100

    ############################################
    # Apatite fission track model params:
    ############################################
    # resample timesteps for AFT calculation
    #
    resample_AFT_timesteps = 10

    # use C-axis correction for apatite fission track lengths
    use_caxis_correction = False

    # parameters for annealing characteristics of apatite grains
    # options for kinetic params:
    # 'Clwt' : Chloride wt fractions
    # 'Dpar' : Dpar / etch pit size
    annealing_kinetic_param = 'Clwt'
    # end member values for kinetic parameters (if no value given in input dataset)
    annealing_kinetics_values = np.array([0.0001, 0.02])

    # size of bins of (simulated) AFT length histogram, default = 0.25 um
    binsize = 0.25

    # annealing equation to use
    # 'FA' for fanning Arrhenius equation by Laslett (1987)
    # 'FC' for fanning curvelinear equation used by Ketcham (1999, 2007)
    annealing_equation = 'FC'

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
    # apatite U-Th/He equations to use
    # 'Farley2000' for helium diffusion parameters of Durango apatite
    #   acc. to Farley(2000) JGR 105
    # 'RDAAM' for he diffusion that depends on radiation damage acc. to
    #   Flowers et al. (2009) GCA 73
    ahe_method = 'Farley2000'

    # decay constants
    decay_constant_238U = 4.916e-18
    decay_constant_232Th = 1.57e-18
    decay_constant_235U = 3.12e-17

    #######
    # VR
    #######
    # default sigma of uncertainty range for VR data,
    # if not specified in input file
    vr_unc_sigma = 0.05


class ParameterRanges:

    """
    parameter ranges for sensitivity or uncertainty analysis

    PyBasin will look for any variable ending with _s below and then look for the
    corresponding variable in the class pybasin_params above

    each _s variable should be a list of values. PyBasin will replace the variable
    in model_parameters.py with each item in the list consecutively
    """

    year = 365.25 * 24 * 60 * 60.0

    # option whether to vary one model parameter at a time
    # (ie for an one at a time sensitivity analysis)
    # or to run all parameter combinations, using the parameter ranges specified below
    parameter_combinations = False

    # option to add a first base run with unchanged parameters to the list of model
    # runs
    initial_base_run = False

    #
    exhumed_thicknesses_s = [[500.0], [1000.0]]
