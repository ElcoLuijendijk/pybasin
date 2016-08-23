"""
set model scenarios: multiple model runs with varying values of exhumation
or basal heat flow
"""

import numpy as np

wells = ['Flaxmans-1']

# strat period for which to change the basal heat flow
# must match a period of heatflow_periods in the param file
# 'all' for entire history
# 'last' for last timestep
# or a list for a set number of timesteps, like this: [0, 2]
basal_heat_flow_scenario_period = [0, 2]

# basal heat flow scenarios to try, heat flow in W m^-2
# example, for testing basal heat flows of 50, 70 and 90 x 10^-3 W/m2:
#basal_heat_flow_scenarios = [75.0e-3]
#basal_heat_flow_scenarios = np.arange(62, 72, 1.0) * 1e-3
basal_heat_flow_scenarios = [66e-3]

# strat period for which to change exhumation
# must exhumation_phase_id further down this file
exhumation_scenarios_period = 'late_miocene_exhumation'

# magnitude of exhumation (m) for each model scenario
# example for testing exhumatwion of 500, 1000 and 1500 m:
#exhumation_magnitudes = [0.0, 500.0]
#exhumation_magnitudes = np.arange(400.0, 650.0, 50)
exhumation_magnitudes = [0.0]

# exhumation phase start (Ma) and duration (My)
#exhumation_starts = np.arange(1.0, 14.0, 2.0)
#exhumation_durations = np.arange(1.0, 13.0, 2.0)
#exhumation_durations = np.arange(0.1, 1.0, 0.1)
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
