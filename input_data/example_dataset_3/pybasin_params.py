"""
model parameters for the 1D burial & temperature history model

edit this file to change the model parameters for PyBasin
"""

import numpy as np


print('-' * 10)
print('Example PyBasin dataset from the Roer Valley Graben (well AST-02), solute diffusion example')
print('-' * 10)


class ModelParameters:

    # location of input data .csv files
    output_dir = 'model_output/example_dataset_3'

    # names of wells or surface outcrops to include in a single set of model runs:
    wells = ['AST-02']

    # option to calculate apatite fission track data
    # (no aft_samples.csv/aft_data.csv, vitrinite_reflectance.csv or he_samples.csv/he_data.csv
    # are included with this dataset, so these are all turned off)
    simulate_AFT = False
    simulate_He = False
    simulate_VR = False
    simulate_salinity = True

    # option to calculate AHe ages for all nodes rather than just the samples
    # in a particular borehole
    # is automatically set to True if make_model_data_fig = True
    # note that this increases the amount of computational time quite a bit
    calculate_thermochron_for_all_nodes = True

    # option to save detailled model run data (approx 10-20 MB per model run)
    save_model_run_data = True

    # save time-temperature paths for each sample, so that they can be used with other codes such as HeFTy or QtQt
    log_tT_paths = True

    # location of detailled output data
    datafile_output_dir = 'model_output/example_dataset_3/thermal_history_datafiles'

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
    contour_variable = 'salinity'

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
    # NOTE: well AST-02 only reaches down to CKHM (Cretaceous chalk group), it
    # does not preserve the deeper Triassic/Carboniferous units (ATAL, ATBR,
    # DCCU, etc.) that the exhumation phases in the NDW-01 example dataset
    # refer to, and no exhumation history specific to AST-02 is included with
    # this dataset. PyBasin requires at least one exhumation phase to be
    # defined, so this uses a placeholder phase at the well documented
    # base-Cenozoic (end-Cretaceous) regional unconformity, with a nominal
    # exhumed thickness below the 5 m threshold that the model treats as
    # negligible (see min_exh_thickness in lib/pybasin_lib.py). Replace this
    # with a real, calibrated exhumation history for AST-02 if you have one;
    # this placeholder is only meant to make the solute diffusion example run.
    # start of exhumation (Ma)
    exhumation_period_starts = np.array([61.7])
    # end of exhumation phase (Ma)
    exhumation_period_ends = np.array([59.2])
    # exhumed thickness (m)
    exhumed_thicknesses = np.array([1.0])

    # determine last deposited units before unconformity:
    # this should be one list for each exhumation phase, with stratigraphic unit codes ordered from old to young
    # the model will add units starting from the oldest to the youngest, untill the additional thickness needed for
    # erosion is filled
    exhumed_strat_units = [['CKHM']]

    # maximum initial (pre-erosion) thicknesses:
    # make sure the last unit is thick enough so that all values of exhumation that you want to test can be accomodated
    original_thicknesses = [[10.0]]

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
    heatflow_history = np.array([65.0, 65.0, 100.0, 100.0]) * 1e-3

    # max size of heatflow timestep (in yrs)
    max_hf_timestep = 10000.0

    #################
    # goodness of fit
    ################
    # weights for calculating overall goodness of fit from the gof statistic for
    # temperature, vitrinite reflectance, apatite fission track age and
    # apatite (U-Th)/He data
    gof_weights = [1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0]

    # use two-sided GOF for age data: also penalizes modelled age range wider
    # than measured age range
    two_sided_gof = False

    # percentile range of measured age PDF used to define the measured age range
    # for the reverse GOF component (fraction of modelled ages within measured range)
    gof_age_percentile = [5, 95]

    #############################
    # Thermochronology parameters
    #############################
    # number of timesteps to discretize the provenance thermal history prior to deposition
    provenance_time_nt = 100

    # temperature at the start of the provenance history of the samples (degr. C)
    provenance_start_temp = 120.0

    # resample timesteps for AFT, AHe calculation and saving modeled temperature histories
    resample_timesteps = 10

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

    #####################
    # solute diffusion
    #####################
    # exponent of the porosity-tortuosity relation used to calculate the
    # effective diffusion coefficient: tortuosity = porosity ** tortuosity_factor
    tortuosity_factor = -1 / 3.

    # option to use a single, constant molecular diffusion coefficient (Dw
    # below) instead of a temperature and concentration dependent value
    constant_diffusivity = False

    # constant molecular diffusion coefficient (m^2 s^-1), only used if
    # constant_diffusivity is True
    Dw = 20.3e-10

    # fixed salinity at the lower boundary of the model domain (kg/kg)
    fixed_lower_bnd_salinity = 0.30

    # salinity of seawater (kg/kg), used for stratigraphic units or time
    # periods marked as marine
    salinity_seawater = 0.035

    # salinity of fresh/meteoric water (kg/kg), used for stratigraphic units
    # or time periods marked as terrestrial
    salinity_freshwater = 0.0001

    # option to read a well specific surface salinity boundary condition
    # history from surface_salinity.csv, instead of only relying on the
    # marine/non-marine history recorded in stratigraphy_info.csv
    well_specific_surface_salinity_bnd = True


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

    # run model scenarios parallel
    parallel_model_runs = False

    # max number of simultaneous model runs:
    max_number_of_processes = 20

    # example for running multiple models
    #exhumed_thicknesses_s = [[1000.0, 750.0], [1500.0, 750.0], [2000.0, 750.0], [2500.0, 750.0],
    #                         [3000.0, 750.0], [3500.0, 750.0], [4000.0, 750.0], [4500.0, 750.0],
    #                         [5000.0, 750.0]]
