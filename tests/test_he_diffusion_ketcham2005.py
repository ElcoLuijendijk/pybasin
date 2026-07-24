"""
Benchmark test for PyBasin's (U-Th)/He diffusion model against reference
values validated against HeFTy.

The three time-temperature histories below are the exact NumPy array
definitions used in tests/test_he.py::test_forward_model_he of the open
source package gdtchron (https://github.com/dyvasey/gdtchron), whose own
docstring states that they come from Figure 10 of Ketcham (2005) and
that the comparison data comes from HeFTy v1.9.3.

gdtchron's AHe diffusion kinetics (frequency factor 50 cm2/sec,
activation energy 138 kJ/mol, attributed to Reiners and Brandon, 2006)
match PyBasin's method='Farley2000' kinetics almost exactly. Since
gdtchron's diffusion-production solver (finite difference) and
PyBasin's calculate_he_age_meesters_dunai_2002 (eigenmode expansion,
Meesters and Dunai, 2002) are two independent numerical solutions of
the same underlying physics, rather than two different fitted models,
close agreement here is a correctness check, not just a plausibility
comparison.

PyBasin's alpha_ejection=True models the physical loss of He from
alpha ejection directly in the diffusion solution, the same concept as
gdtchron's uncorrected age (both are raw ages from an ejection-affected
solution, with no separate geometric correction factor applied
afterward), so that is the pairing checked here.

See tests/analytical_solution_He_diffusion_Ketcham2005.ipynb for the
full comparison, including the alpha_ejection=False versus corrected
age comparison, plots and discussion.

Run with:
    pytest tests/test_he_diffusion_ketcham2005.py
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lib.helium_diffusion_models as he

YEAR = 365.25 * 24 * 3600.0
MA = 1.0e6 * YEAR

TIMES = np.arange(60, -0.001, -0.1)


def build_history_1():
    early_temps = np.linspace(120, 20, 101)
    late_temps = np.linspace(20, 20, 500)
    return np.append(early_temps, late_temps)


def build_history_2():
    return np.linspace(120, 20, 601)


def build_history_3():
    early_temps = np.linspace(120, 65, 576)
    late_temps = np.linspace(65, 20, 26)
    return np.append(early_temps, late_temps[1:])


# name, temperature builder, reference uncorrected age (Ma)
KETCHAM2005_HE_CASES = [
    ('history 1, rapid cooling then held at 20 C', build_history_1, 46.5),
    ('history 2, constant cooling', build_history_2, 22.5),
    ('history 3, slow cooling then rapid final drop', build_history_3, 5.84),
]


@pytest.mark.parametrize(
    'name, builder, reference_age_uncorrected',
    KETCHAM2005_HE_CASES,
    ids=[case[0] for case in KETCHAM2005_HE_CASES])
def test_he_age_matches_gdtchron_reference(
        name, builder, reference_age_uncorrected):
    """
    method='Farley2000' with alpha_ejection=True should reproduce the
    gdtchron / HeFTy v1.9.3 uncorrected reference age for each of the
    three Ketcham (2005) Figure 10 histories to within a few percent,
    since this compares two independent numerical solutions of the
    same diffusion-production equation.
    """
    temps_c = builder()
    t_ma = TIMES.max() - TIMES
    t_sec = t_ma * MA
    T_kelvin = temps_c + 273.15

    age = he.calculate_he_age_meesters_dunai_2002(
        t_sec, T_kelvin, radius=100e-6, U=100e-6, Th=100e-6,
        method='Farley2000', alpha_ejection=True, n_eigenmodes=50)

    age_ma = age[-1] / MA

    assert age_ma == pytest.approx(reference_age_uncorrected, rel=0.05)
