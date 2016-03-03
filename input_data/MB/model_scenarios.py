"""
set model scenarios: multiple model runs with varying values of exhumation
or basal heat flow
"""

import numpy as np

# list of wells to include in model runs
#wells = ['Altishofen', 'Boswil', 'Chapelle', 'Entlebuch',
#         'Huenenberg', 'Herdern', 'Lindau', 'Linden',
#         'Rigi', 'Romanens', 'Savigny', 'Thun', 'Tuggen', 'Weggis']

#wells = ['MRP025', 'MRP067', 'MRP087', 'MRP117',
#         'MRP170', 'MRP172', 'MRP174', 'MRP179', 'MRP187', 'MRP198', 'MRT052',
#         'E20', 'E30', 'E35', 'E40', 'E45', 'E50', 'E55',
#         'E60', 'E65', 'E70b',
#         'RH10', 'RH15', 'RH25', 'RH30a', 'RH30b', 'RH35', 'RH40',
#         'RH45', 'RH50', 'RH60c', 'RH65', 'RH70', 'RV15', 'RV20', 'RV25',
#         'RV30a', 'RV30b', 'RV30c',
#         'B25', 'B30', 'B55', 'B60']

# Rigi wells:
wells = ['Huenenberg',
         'MRP025', 'MRP170', 'MRP172', 'MRP174',
         'RH10', 'RH12', 'RH15', 'RH17', 'RH20', 'RH23', 'RH30a', 'RH30b',
         'RH35', 'RH40', 'RH45', 'RH50', 'RH60c', 'RH65', 'RH70',
         'RV05', 'RV10b', 'RV15', 'RV20', 'RV25', 'RV30c',
         'Weggis', 'Entlebuch']

#wells = ['RH30b']
#wells = ['RH30a', 'RH30b']
wells = ['RV30c']
# strat period for which to change the basal heat flow
# must match a period of heatflow_periods in the param file
basal_heat_flow_scenario_period = 'all'

# basal heat flow scenarios to try, heat flow in mW m^-2
# example, for testing basal heat flows of 50, 70 and 90 x 10^-3 W/m2:
#basal_heat_flow_scenarios = [50e-3, 70e-3, 90e-3]
#basal_heat_flow_scenarios = np.arange(40, 105, 5.0) * 1e-3
basal_heat_flow_scenarios = [65.0 * 1e-3]

# strat period for which to change exhumation
# must exhumation_phase_id further down this file
exhumation_scenarios_period = 'molasse_exhumation'

# exhumation (m) for each model scenario
# example for testing exhumation of 500, 1000 and 1500 m:
# exhumation_magnitudes = [500, 1000, 1500]
#exhumation_magnitudes = np.arange(0., 5000, 500)
exhumation_magnitudes = [2000.0]

# exhumation phase start and duration
# list of lists or list of tuples. first value of inner list or tuple is start of
# exhumation (Ma), second is duration of exhumation (My)
# example to test 2 scenarios, with exhumation starting at 10 Ma and 7.5 Ma
# and duration of 5 My and 2.5 My:
# exhumation_starts_and_durations = [[10.0, 5.0], [7.5, 2.5]]

#exhumation_starts_and_durations = [[12.0, 10.0],
#                                   [12.0, 5.5],
#                                   [12.0, 1.0],
#                                   [8.0, 5.5],
#                                   [8.0, 1.0],
#                                   [4.0, 1.0]]
#exhumation_starts_and_durations = [[10.0, 7.0]]
#exhumation_starts = np.arange(1.0, 14.0, 2.0)
#exhumation_durations = np.arange(1.0, 13.0, 2.0)
#exhumation_durations = np.arange(0.1, 1.0, 0.1)
exhumation_starts = [9.0]
exhumation_durations = [7.0]