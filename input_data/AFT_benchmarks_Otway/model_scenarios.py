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
basal_heat_flow_scenario_period = 'first'

# basal heat flow scenarios to try, heat flow in W m^-2
# example, for testing basal heat flows of 50, 70 and 90 x 10^-3 W/m2:
#basal_heat_flow_scenarios = [50e-3, 70e-3, 90e-3]
#basal_heat_flow_scenarios = np.arange(70, 80, 1.0) * 1e-3
basal_heat_flow_scenarios = [None]

# strat period for which to change exhumation
# must exhumation_phase_id further down this file
exhumation_scenarios_period = 'late_miocene_exhumation'

# magnitude of exhumation (m) for each model scenario
# example for testing exhumation of 500, 1000 and 1500 m:
# exhumation_magnitudes = [500, 1000, 1500]
exhumation_magnitudes = np.arange(0.0, 1100, 100)
#exhumation_magnitudes = [0.0]

# exhumation phase start (Ma) and duration (My)
#exhumation_starts = np.arange(1.0, 14.0, 2.0)
#exhumation_durations = np.arange(1.0, 13.0, 2.0)
#exhumation_durations = np.arange(0.1, 1.0, 0.1)
exhumation_starts = [10.0]
exhumation_durations = [10.0]

# AFT annealing eq. parameters, see Ketcham et al. (2007) Am. Min.
AFT_C0 = [None]
AFT_C1 = [None]
AFT_C2 = [None]
AFT_C3 = [None]
AFT_C4 = [None]
AFT_alpha = [None]
