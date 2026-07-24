"""
Benchmark test for PyBasin's (U-Th)/He diffusion model against reference
ages from Ketcham (2005).

The three time-temperature histories below reproduce Figure 10 of
Ketcham, R. A. (2005), Forward and Inverse Modeling of Low-Temperature
Thermochronometry Data, Reviews in Mineralogy and Geochemistry 58,
275-314. That figure's own caption states the ages are "reported using
the alpha correction factor of Farley et al. (1996)," ie. the raw,
ejection affected age divided by a purely geometric correction factor,
not the raw age itself.

PyBasin's alpha_ejection=True models the physical loss of He from
alpha ejection directly in the diffusion solution (Meesters and Dunai,
2002, Part II, eq. 24). This is checked here against the reference
after additionally dividing by farley1996_alpha_correction, PyBasin's
implementation of the Farley et al. (1996) geometric correction factor
(following eqs. 27-30 of Ketcham, 2005, for a sphere with uniformly
distributed U and Th).

Figure 10's caption describes the grain as "a grain diameter of 100
micrometre," but the body text of the same chapter, describing the
same figure, calls it "a 100 micrometre radius apatite." These cannot
both be correct. Checked directly: a 100 micrometre radius (the body
text) gives 0.8 to 5.9 percent agreement here; a 50 micrometre radius
(taking the caption's "100 micrometre diameter" literally) gives 5 to
47 percent disagreement. The body text's radius is used here.

See tests/analytical_solution_He_diffusion_Ketcham2005.ipynb for the
full comparison, plots and discussion.

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

RADIUS = 100e-6
U = 100e-6
TH = 100e-6


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


# name, temperature builder, reference alpha corrected age (Ma)
KETCHAM2005_HE_CASES = [
    ('history 1, rapid cooling then held at 20 C', build_history_1, 54.6),
    ('history 2, constant cooling', build_history_2, 26.5),
    ('history 3, slow cooling then rapid final drop', build_history_3, 7.12),
]


@pytest.mark.parametrize(
    'name, builder, reference_corrected_age',
    KETCHAM2005_HE_CASES,
    ids=[case[0] for case in KETCHAM2005_HE_CASES])
def test_he_age_matches_ketcham2005_reference(
        name, builder, reference_corrected_age):
    """
    method='Farley2000' with alpha_ejection=True, followed by dividing
    by the Farley et al. (1996) correction factor, should reproduce the
    alpha corrected age reported in Ketcham (2005) Figure 10 for each of
    the three histories to within a few percent.
    """
    temps_c = builder()
    t_ma = TIMES.max() - TIMES
    t_sec = t_ma * MA
    T_kelvin = temps_c + 273.15

    raw_age = he.calculate_he_age_meesters_dunai_2002(
        t_sec, T_kelvin, RADIUS, U=U, Th=TH,
        method='Farley2000', alpha_ejection=True, n_eigenmodes=50)
    raw_age_ma = raw_age[-1] / MA

    tau = he.farley1996_alpha_correction(RADIUS, U, TH)
    corrected_age_ma = raw_age_ma / tau

    assert corrected_age_ma == pytest.approx(reference_corrected_age, rel=0.07)
