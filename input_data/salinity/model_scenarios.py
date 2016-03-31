"""
set model scenarios: multiple model runs with varying values of exhumation
or basal heat flow
"""

import numpy as np

wells = ['AST-02']

# strat period for which to change the basal heat flow
# must match a period of heatflow_periods in the param file
basal_heat_flow_scenario_period = 'all'

# basal heat flow scenarios to try, heat flow in mW m^-2
# example, for testing basal heat flows of 50, 70 and 90 x 10^-3 W/m2:
#basal_heat_flow_scenarios = [50e-3, 70e-3, 90e-3]
basal_heat_flow_scenarios = [65.0 * 1e-3]

# strat period for which to change exhumation
# must exhumation_phase_id further down this file
exhumation_scenarios_period = 'late_cretaceous_exhumation'

# exhumation (m) for each model scenario
# example for testing exhumation of 500, 1000 and 1500 m:
# exhumation_magnitudes = [500, 1000, 1500]
#exhumation_magnitudes = np.arange(0., 5000, 500)
exhumation_magnitudes = [500.0]

# exhumation phase start and duration
# list of lists or list of tuples. first value of inner list or tuple is start of
# exhumation (Ma), second is duration of exhumation (My)
# example to test 2 scenarios, with exhumation starting at 10 Ma and 7.5 Ma
# and duration of 5 My and 2.5 My:
#exhumation_starts_and_durations = [[None, None]]

exhumation_starts = [None]
exhumation_durations = [None]

exhumation_segment_factors = [None]
exhumation_duration_factors = [None]

AFT_C0 = [None]
AFT_C1 = [None]
AFT_C2 = [None]
AFT_C3 = [None]
AFT_C4 = [None]
AFT_alpha = [None]
