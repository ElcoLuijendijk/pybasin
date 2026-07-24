"""
Benchmark test for the Ketcham (1999, 2000) apatite fission track
annealing model against reference ages from Ketcham (2005).

The four time-temperature histories below reproduce Figure 7 of
Ketcham, R. A. (2005), Forward and Inverse Modeling of Low-Temperature
Thermochronometry Data, Reviews in Mineralogy and Geochemistry 58,
275-314. Using method='Ketcham2000' in simulate_AFT_annealing selects
the Ketcham et al. (1999) fanning curvilinear kinetics that Figure 7
was generated with (as opposed to the default method='Ketcham2007'),
so this is a same model comparison, not a cross model plausibility
check.

Only ages are checked here: Figure 7 does not report numeric mean
track lengths, only track length density plots, so there is no
manuscript sourced reference value to compare a modeled mean length
against.

See tests/analytical_solution_AFT_annealing_Ketcham2005.ipynb for the
full comparison, plots and discussion.

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
    initial_temps = np.linspace(190.0, 35, 20, endpoint=False)
    initial_tsteps = np.linspace(93.0, 89, 20, endpoint=False)
    second_temps = np.linspace(35.0, 20, 96 * 5, endpoint=True)
    second_tsteps = np.linspace(89, 0, 96 * 5, endpoint=True)
    timesteps_age = np.concatenate((initial_tsteps, second_tsteps))
    temperatures = np.concatenate((initial_temps, second_temps))
    return timesteps_age, temperatures


def build_history_7b():
    timesteps_age = np.linspace(93, 0, 500, endpoint=True)
    temperatures = np.linspace(190.0, 20, 500, endpoint=True)
    return timesteps_age, temperatures


def build_history_7c():
    first_temps = np.linspace(190.0, 20.0, 11 * 5, endpoint=False)
    first_tsteps = np.linspace(93.0, 80.0, 11 * 5, endpoint=False)
    second_temps = np.linspace(20.0, 90.0, 45 * 5, endpoint=False)
    second_tsteps = np.linspace(80.0, 38.0, 45 * 5, endpoint=False)
    third_temps = np.linspace(92.0, 20.0, 43 * 5, endpoint=True)
    third_tsteps = np.linspace(38.0, 0, 43 * 5, endpoint=True)
    timesteps_age = np.concatenate((first_tsteps, second_tsteps, third_tsteps))
    temperatures = np.concatenate((first_temps, second_temps, third_temps))
    return timesteps_age, temperatures


def build_history_7d():
    initial_temps = np.linspace(190.0, 120.0, 13 * 5, endpoint=False)
    initial_tsteps = np.linspace(93.0, 81.0, 13 * 5, endpoint=False)
    second_temps = np.linspace(120.0, 100.0, 75 * 5, endpoint=False)
    second_tsteps = np.linspace(81.0, 11.0, 75 * 5, endpoint=False)
    third_temps = np.linspace(101.1, 20.0, 13 * 5, endpoint=True)
    third_tsteps = np.linspace(11.0, 0, 13 * 5, endpoint=True)
    timesteps_age = np.concatenate((initial_tsteps, second_tsteps, third_tsteps))
    temperatures = np.concatenate((initial_temps, second_temps, third_temps))
    return timesteps_age, temperatures


# name, history builder, Dpar, reference age (Ma)
KETCHAM2005_CASES = [
    ('7a fast cooling', build_history_7a, 1.75, 89.4),
    ('7b constant cooling', build_history_7b, 1.75, 39.9),
    ('7c reheating', build_history_7c, 1.75, 64.7),
    ('7d two cooling rates, Dpar 1.75', build_history_7d, 1.75, 13.5),
    ('7d two cooling rates, Dpar 2.50', build_history_7d, 2.50, 52.3),
]


@pytest.mark.parametrize(
    'name, builder, dpar, reference_age',
    KETCHAM2005_CASES,
    ids=[case[0] for case in KETCHAM2005_CASES])
def test_ketcham2000_matches_ketcham2005_reference(
        name, builder, dpar, reference_age):
    """
    method='Ketcham2000' should reproduce the age read from Ketcham
    (2005) Figure 7 for each of the four histories to within a few
    percent, since this is now a same model comparison. History 7d at
    Dpar 2.50 is the strictest case: the Ketcham (2007) model used
    elsewhere in this repository disagrees with the reference by about
    60 percent there, since the two annealing models genuinely differ,
    not just in fitted constants.

    The reference ages are read directly off the published figure
    rather than taken from an exact machine readable source, so the
    tolerance here (3 percent) is looser than it would be for a digit
    for digit reference; the actual differences observed range from
    0.3 to 1.8 percent.
    """
    timesteps_age, temperatures = builder()
    timesteps = timesteps_age.max() - timesteps_age

    result = AFTannealingLib.simulate_AFT_annealing(
        timesteps, temperatures, dpar,
        kinetic_parameter='Dpar',
        apply_c_axis_correction=True,
        method='Ketcham2000',
        surpress_resampling=False)
    age = result[1]

    assert age == pytest.approx(reference_age, rel=0.03)


def test_ketcham2000_dpar_to_rmr0():
    """
    Direct check of the Ketcham et al. (1999) Dpar to rmr0 conversion
    used by method='Ketcham2000':

        rmr0 = 1 - exp(0.647 * (Dpar - 1.75) - 1.834)

    and kappa = 1 - rmr0.
    """
    for dpar in [1.75, 2.50]:
        rmr0, kappa = AFTannealingLib.calculate_kinetic_parameters(
            'Dpar', dpar, method='Ketcham2000')

        expected_rmr0 = 1.0 - np.exp(0.647 * (dpar - 1.75) - 1.834)

        assert rmr0 == pytest.approx(expected_rmr0)
        assert kappa == pytest.approx(1.0 - expected_rmr0)
