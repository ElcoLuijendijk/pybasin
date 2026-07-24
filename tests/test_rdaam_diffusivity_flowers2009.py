"""
Validation test for PyBasin's RDAAM helium diffusivity model against
Flowers, R. M., Ketcham, R. A., Shuster, D. L. and Farley, K. A. (2009),
Apatite (U-Th)/He thermochronometry using a radiation damage
accumulation and annealing model, Geochimica et Cosmochimica Acta 73,
2347-2365 (see tests/literature/flowers_RDAAM.pdf).

lib.helium_diffusion_models.calculate_RDAAM_diffusivity implements the
paper's Eq. (1) diffusivity equation with the Table 1, parameter set 2
constants, fed by the radiation damage proxy e-rho_s accumulated via
Eqs. (5)-(8), which in turn depend on the Ketcham et al. (2007)
equivalent-time annealing recurrence (the paper's Eqs. (2)-(4)) already
covered independently in test_aft_annealing_ketcham2000.py.

The functions below (_reference_reduced_length,
_reference_normalized_density, _reference_rdaam_diffusivity) are a
second, independent transcription of Eqs. (1)-(8) written directly from
the paper text, not copied from lib/helium_diffusion_models.py or
lib/AFTannealingLib.py. Comparing PyBasin's production function against
this independent reference catches transcription errors (wrong
constant, sign, or exponent) that a test calling PyBasin's own helper
functions internally could not catch.

As of this writing, all three cases below FAIL, because writing this
test surfaced two real bugs in calculate_RDAAM_diffusivity itself (not
in the reference above, which was cross-checked step by step against
lib.AFTannealingLib's own functions before being trusted as an oracle):

1. lib/helium_diffusion_models.py, use_fortran_algorithm=True path
   (the default): the call to calculate_reduced_AFT_lengths.reduced_ln
   omits the annealing_eq_f90 argument that lib/AFTannealingLib.py's
   own calls to the same Fortran subroutine always pass. Every
   positional argument after rmr0/kappa therefore lands one slot to the
   left of where the compiled extension expects it (confirmed via the
   extension's own f2py-generated signature), so alpha is read as
   annealing_eq_f90 (truncated to an int, selecting the wrong annealing
   equation branch entirely), and C0-C3 are each read one slot early.
   Concretely, this returns rc=0 (fully annealed) for tracks that
   should be nearly pristine, e.g. rc~0.94 for the youngest tracks in
   the "monotonic cooling" case below.
2. Same file, use_fortran_algorithm=False path: calculate_reduced_
   track_lengths is called with the full temperature array (length
   nsteps+1) rather than temperature_midpoint (length nsteps), unlike
   the Fortran path a few lines above it, which correctly uses
   temperature_midpoint. This silently uses the wrong (shifted-by-one)
   temperature for each annealing step whenever temperature actually
   varies over the history, giving a few percent to ~10 percent error
   in the cases below (it is invisible for an isothermal history, which
   is why the simplest case here does not expose it, though it still
   fails from Fortran-path bug 1 above when use_fortran_algorithm
   defaults to True).

Run with:
    pytest tests/test_rdaam_diffusivity_flowers2009.py
"""

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lib.helium_diffusion_models as hdm

YEAR = 365.25 * 24 * 3600.0
MA = 1.0e6 * YEAR

# Flowers et al. (2009) Table 1, parameter set 2, the set adopted for
# the modeling in their Section 5.
LOG_OMEGA_P = -22.0
LOG_PHI_P = -13.0
E_TRAP = 34.0e3            # J/mol
LN_D0L_A2 = 9.733          # ln(D_oL/a^2), 1/s
E_L = 122.3e3              # J/mol
GAS_CONSTANT = 8.3144621

# Ketcham et al. (2007) fanning curvilinear annealing parameters, used
# by the paper's Eqs. (2)-(4)
C0, C1, C2, C3, ALPHA = 0.39528, 0.01073, -65.12969, -7.91715, 0.04672

# Eq. (8) constants
ETA_Q = 0.91
L_ETCH = 8.1e-4            # cm, "~8.1 um for an unannealed track"
LAMBDA_F = 8.46e-17 / YEAR  # 1/s, Jonckheere (2003) value adopted in the paper
LAMBDA_D = 1.551e-4 / (YEAR * 1e6)  # 1/s, total 238U decay constant

DECAY_U238 = 4.916e-18
DECAY_U235 = 3.12e-17
DECAY_TH232 = 1.57e-18
ATOMIC_MASS_U238 = 238.05078826
ATOMIC_MASS_U235 = 235.0439299
ATOMIC_MASS_TH232 = 232.0377
APATITE_DENSITY = 3190.0   # kg/m3
AVOGADRO = 6.022e23


def _reference_reduced_length(dts, temperatures_k):
    """Eq. (2): Ketcham et al. (2007) equivalent-time annealing
    recurrence, for c-axis projected reduced length r_c,B2."""
    nsteps = len(dts)
    rc = np.zeros(nsteps)
    for j in range(nsteps):
        g_prev = 0.0
        for i in range(j, nsteps):
            T = temperatures_k[i]
            if i == j:
                dteq = 0.0
            else:
                dteq = math.exp(
                    ((g_prev - C0) / C1 * (math.log(1.0 / T) - C3)) + C2)
            g_prev = C0 + C1 * (
                (math.log(dts[i] + dteq) - C2) / (math.log(1.0 / T) - C3))
        rc[j] = 1.0 / (g_prev ** (1.0 / ALPHA) + 1.0)
    return rc


def _reference_kinetic_modifier(rc, rmr0, kappa):
    """Eq. (3): r_c,lr from r_c,B2 and the rmr0/kappa kinetic parameters.

    rc < rmr0 (a track more annealed than the totally-annealed
    threshold) is physically a zero reduced length, not a fractional
    power of a negative number. lib.AFTannealingLib.kinetic_modifier_
    reduced_lengths relies on its @jit(nopython=True) compiled '**'
    happening to return 0.0 for this case; plain NumPy returns NaN for
    the same formula. The clip below reproduces the intended (and
    Numba's actual) behaviour explicitly rather than relying on that
    compiled-vs-interpreted difference.
    """
    with np.errstate(invalid='ignore'):
        result = ((rc - rmr0) / (1.0 - rmr0)) ** kappa
    return np.where(rc > rmr0, result, 0.0)


def _reference_normalized_density(r, r_crit=0.765, r_min=0.5275):
    """Eq. (4): reduced length to normalized track density rho_r.

    Flowers et al. (2009) print only the two-piece r_crit split shown
    below. lib.AFTannealingLib.calculate_normalized_density additionally
    floors any r <= r_min (0.5275) to zero, citing the cutoff value
    reported by Green (1988) via the broader Ketcham (2000) annealing
    model that Eq. (4) is drawn from but does not spell out. Without
    this floor, the quadratic branch evaluated at r=0 returns its
    constant term (d5 = 2.269) rather than zero, which is wrong for
    fully annealed tracks (e.g. those already clipped to rc=0 by
    _reference_kinetic_modifier above).
    """
    d1, d2, d3, d4, d5 = 1.600, -0.600, 9.205, -9.157, 2.269
    rho = np.where(r >= r_crit, d1 * r + d2, d3 * r ** 2 + d4 * r + d5)
    rho = np.where(r <= r_min, 0.0, rho)
    return np.clip(rho, 0.0, None)


def _reference_rdaam_diffusivity(temperature_k, time_sec, U_frac, Th_frac,
                                  radius_m, rmr0):
    """Independent reimplementation of Eqs. (1), (5)-(8), combined with
    the annealing chain above, following the paper directly rather than
    lib/helium_diffusion_models.py."""

    U238_frac = (137.88 / 138.88) * U_frac
    U235_frac = (1.0 / 138.88) * U_frac
    Th232_frac = Th_frac

    def frac_to_atoms_cm3(frac, atomic_mass):
        atoms_per_kg = (frac * 1000.0 / atomic_mass) * AVOGADRO
        return atoms_per_kg * APATITE_DENSITY / (100.0 ** 3)

    U238_cm3 = frac_to_atoms_cm3(U238_frac, ATOMIC_MASS_U238)
    U235_cm3 = frac_to_atoms_cm3(U235_frac, ATOMIC_MASS_U235)
    Th232_cm3 = frac_to_atoms_cm3(Th232_frac, ATOMIC_MASS_TH232)

    t1, t2 = time_sec[:-1], time_sec[1:]

    # Eq. (7): volume density of tracks, weighted 8/8, 7/8, 6/8 by
    # alpha-production per parent decay
    rho_v = (
        U238_cm3 * (np.exp(DECAY_U238 * t2) - np.exp(DECAY_U238 * t1))
        + 7.0 / 8.0 * U235_cm3 * (np.exp(DECAY_U235 * t2) - np.exp(DECAY_U235 * t1))
        + 6.0 / 8.0 * Th232_cm3 * (np.exp(DECAY_TH232 * t2) - np.exp(DECAY_TH232 * t1))
    )

    dts = t2 - t1
    T_mid = (temperature_k[1:] + temperature_k[:-1]) / 2.0
    kappa = 1.04 - rmr0

    rc_b2 = _reference_reduced_length(dts, T_mid)
    rc_lr = _reference_kinetic_modifier(rc_b2, rmr0, kappa)
    rho_r = _reference_normalized_density(rc_lr)

    # Eq. (8), then cumulative sum ("sum the density contributions from
    # all time steps prior to the one in which diffusivity is being
    # calculated")
    e_rho_s = (LAMBDA_F / LAMBDA_D) * rho_v * ETA_Q * L_ETCH * rho_r
    e_rho_s_sum = np.cumsum(e_rho_s)

    # Eq. (1), with Psi_rho * e_rho_s + Omega_rho * e_rho_s^3 in place
    # of Psi * [4He] (Section 3.3)
    C = 10 ** LOG_PHI_P * e_rho_s_sum + 10 ** LOG_OMEGA_P * e_rho_s_sum ** 3
    D_div_a2 = (
        np.exp(LN_D0L_A2) * np.exp(-E_L / (GAS_CONSTANT * T_mid))
        / (C * np.exp(E_TRAP / (GAS_CONSTANT * T_mid)) + 1.0))

    D_midpoint = D_div_a2 * radius_m ** 2

    # calculate_RDAAM_diffusivity returns D averaged from midpoint
    # resolution back onto the original, full time-array resolution;
    # replicate that same step so both sides return the same shape.
    D_full = np.zeros_like(temperature_k)
    D_full[0] = D_midpoint[0]
    D_full[-1] = D_midpoint[-1]
    D_full[1:-1] = (D_midpoint[1:] + D_midpoint[:-1]) / 2.0
    return D_full


# name, radius (m), rmr0, U (ppm), Th (ppm), time-temperature history
RDAAM_CASES = [
    (
        'negligible annealing, short isothermal history',
        60e-6, 0.83, 100.0, 50.0,
        np.array([0.0, 1000.0 * YEAR]),
        np.array([323.15, 323.15]),
    ),
    (
        'monotonic cooling, 100 Myr, typical apatite',
        60e-6, 0.83, 28.0, 12.0,
        np.linspace(0.0, 100.0 * MA, 101),
        np.linspace(120.0, 20.0, 101) + 273.15,
    ),
    (
        'monotonic cooling, small high-eU grain, resistant kinetics',
        40e-6, 0.75, 150.0, 20.0,
        np.linspace(0.0, 50.0 * MA, 51),
        np.linspace(140.0, 10.0, 51) + 273.15,
    ),
]


@pytest.mark.parametrize(
    'name, radius, rmr0, U_ppm, Th_ppm, time_sec, temperature_k',
    RDAAM_CASES, ids=[case[0] for case in RDAAM_CASES])
def test_rdaam_diffusivity_matches_flowers2009(
        name, radius, rmr0, U_ppm, Th_ppm, time_sec, temperature_k):
    """
    calculate_RDAAM_diffusivity should match an independent
    transcription of Flowers et al. (2009) Eqs. (1)-(8) with the
    Table 1 parameter set 2 constants, to within floating point
    summation-order differences.
    """
    U_frac = U_ppm * 1e-6
    Th_frac = Th_ppm * 1e-6

    U238 = (137.88 / 138.88) * U_frac
    U235 = (1.0 / 138.88) * U_frac
    Th232 = Th_frac

    D_pybasin = hdm.calculate_RDAAM_diffusivity(
        temperature_k, time_sec, U238, U235, Th232, radius,
        kinetic_parameter='rmr0', kinetic_value=rmr0)

    D_reference = _reference_rdaam_diffusivity(
        temperature_k, time_sec, U_frac, Th_frac, radius, rmr0)

    np.testing.assert_allclose(D_pybasin, D_reference, rtol=1e-6)
