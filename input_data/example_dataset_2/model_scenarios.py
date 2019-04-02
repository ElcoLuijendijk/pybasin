"""
set model scenarios: multiple model runs with varying values of exhumation
or basal heat flow
"""

import numpy as np

# list of wells to include in model runs
wells = ['BKZ-01']


# strat period for which to change the basal heat flow
# must match a period of heatflow_periods in the param file
basal_heat_flow_scenario_period = 'all'

# basal heat flow scenarios to try, heat flow in mW m^-2
# example, for testing basal heat flows of 50, 70 and 90 x 10^-3 W/m2:
basal_heat_flow_scenarios = [None]

# strat period for which to change exhumation
# must exhumation_phase_id further down this file
exhumation_scenarios_period = 'late_cretaceous_unc'

# exhumation (m) for each model scenario
# example for testing exhumation of 500, 1000 and 1500 m:
exhumation_magnitudes = [250.0, 1000.0]

# exhumation phase start (Ma) and duration (My)
exhumation_starts = [None]
exhumation_durations = [None]

exhumation_segment_factors = [None]
exhumation_duration_factors = [None]

# AFT annealing eq. parameters, see Ketcham et al. (2007) Am. Min.
AFT_C0 = [None]
AFT_C1 = [None]
AFT_C2 = [None]
AFT_C3 = [None]
AFT_C4 = [None]
AFT_alpha = [None]
