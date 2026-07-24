"""
Benchmark test for the Ketcham (1999, 2000) apatite fission track
annealing model against reference values validated against HeFTy.

The four time-temperature histories below are the exact NumPy array
definitions used in tests/test_aft.py::test_forward_model_aft of the
open source package gdtchron (https://github.com/dyvasey/gdtchron),
whose own docstring states that they come from Figure 7 of Ketcham
(2005) and that the comparison data comes from HeFTy v1.9.3. Using
method='Ketcham2000' in simulate_AFT_annealing selects the same
Ketcham et al. (1999) fanning curvilinear kinetics gdtchron implements
(as opposed to the default method='Ketcham2007'), so this is a same
model comparison, not a cross model plausibility check.

See tests/analytical_solution_AFT_annealing_Ketcham2005.ipynb for the
full comparison, including the Ketcham (2007) results, plots and
discussion.

Run with:
    pytest tests/test_aft_annealing_ketcham2000.py
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lib.AFTannealingLib as AFTannealingLib


def build_history_7a():
    initial_temps = np.linspace(190.1, 35, 20, endpoint=False)
    initial_tsteps = np.linspace(93, 89.1, 20, endpoint=False)
    second_temps = np.linspace(35, 20, 96 * 5, endpoint=True)
    second_tsteps = np.linspace(89.1, 0, 96 * 5, endpoint=True)
    timesteps_age = np.concatenate((initial_tsteps, second_tsteps))
    temperatures = np.concatenate((initial_temps, second_temps))
    return timesteps_age, temperatures


def build_history_7b():
    timesteps_age = np.linspace(93, 0, 500, endpoint=True)
    temperatures = np.linspace(190.1, 20, 500, endpoint=True)
    return timesteps_age, temperatures


def build_history_7c():
    first_temps = np.linspace(190.1, 19.8, 11 * 5, endpoint=False)
    first_tsteps = np.linspace(93, 79.5, 11 * 5, endpoint=False)
    second_temps = np.linspace(19.8, 92, 45 * 5, endpoint=False)
    second_tsteps = np.linspace(79.5, 37.9, 45 * 5, endpoint=False)
    third_temps = np.linspace(92, 20, 43 * 5, endpoint=True)
    third_tsteps = np.linspace(37.9, 0, 43 * 5, endpoint=True)
    timesteps_age = np.concatenate((first_tsteps, second_tsteps, third_tsteps))
    temperatures = np.concatenate((first_temps, second_temps, third_temps))
    return timesteps_age, temperatures


def build_history_7d():
    initial_temps = np.linspace(190.1, 121.7, 13 * 5, endpoint=False)
    initial_tsteps = np.linspace(93, 81, 13 * 5, endpoint=False)
    second_temps = np.linspace(121.7, 101.1, 75 * 5, endpoint=False)
    second_tsteps = np.linspace(81, 10.9, 75 * 5, endpoint=False)
    third_temps = np.linspace(101.1, 20, 13 * 5, endpoint=True)
    third_tsteps = np.linspace(10.9, 0, 13 * 5, endpoint=True)
    timesteps_age = np.concatenate((initial_tsteps, second_tsteps, third_tsteps))
    temperatures = np.concatenate((initial_temps, second_temps, third_temps))
    return timesteps_age, temperatures


# name, history builder, Dpar, reference age (Ma), reference mean length (um)
KETCHAM2005_CASES = [
    ('7a fast cooling', build_history_7a, 1.75, 89.4, 15.03),
    ('7b constant cooling', build_history_7b, 1.75, 39.8, 14.24),
    ('7c reheating', build_history_7c, 1.75, 63.3, 13.38),
    ('7d two cooling rates, Dpar 1.75', build_history_7d, 1.75, 13.0, 13.95),
    ('7d two cooling rates, Dpar 2.50', build_history_7d, 2.50, 51.1, 12.91),
]


@pytest.mark.parametrize(
    'name, builder, dpar, reference_age, reference_mean_length',
    KETCHAM2005_CASES,
    ids=[case[0] for case in KETCHAM2005_CASES])
def test_ketcham2000_matches_gdtchron_reference(
        name, builder, dpar, reference_age, reference_mean_length):
    """
    method='Ketcham2000' should reproduce the gdtchron / HeFTy v1.9.3
    reference age and mean length for each of the four Ketcham (2005)
    Figure 7 histories to within a few percent, since this is now a
    same model comparison. History 7d at Dpar 2.50 is the strictest
    case: the Ketcham (2007) model used elsewhere in this repository
    disagrees with the reference by 60 percent there, since the two
    annealing models genuinely differ, not just in fitted constants.
    """
    timesteps_age, temperatures = builder()
    timesteps = timesteps_age.max() - timesteps_age

    result = AFTannealingLib.simulate_AFT_annealing(
        timesteps, temperatures, dpar,
        kinetic_parameter='Dpar',
        apply_c_axis_correction=True,
        method='Ketcham2000',
        surpress_resampling=False)
    age, mean_length = result[1], result[2]

    assert age == pytest.approx(reference_age, rel=0.02)
    assert mean_length == pytest.approx(reference_mean_length, rel=0.02)
