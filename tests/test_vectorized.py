"""
Comparison tests for vectorized performance-optimized functions.

Each test calls the original (reference) function and its vectorized
equivalent with identical inputs and asserts that results are numerically
identical (np.allclose with default tolerances: rtol=1e-5, atol=1e-8).

These tests guard the equivalence of the refactor, not the correctness of
the underlying physics: if the original function is wrong, the vectorized
one is expected to be wrong in the same way and the test still passes.

Wall-clock benchmarks for these function pairs live in the standalone
script benchmarks/benchmark_vectorized.py; they are intentionally not part
of the test suite because timing is environment dependent and not
diagnostic of model state.

Run with:
    MPLBACKEND=Agg python -m pytest tests/test_vectorized.py -v
"""

import os
import sys
import numpy as np

# Ensure the project root is on the path so lib imports work
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lib.AFTannealingLib as AFTannealingLib
from lib.helium_diffusion_models import (
    He_diffusion_Meesters_and_Dunai_2002,
    He_diffusion_Meesters_and_Dunai_2002_vectorized,
)
import lib.easyRo as easyRo


# ---------------------------------------------------------------------------
# Synthetic thermal history helpers
# ---------------------------------------------------------------------------

def make_thermal_history(n=200, seed=42):
    """
    Build a realistic synthetic thermal history for AFT testing.

    Returns
    -------
    dts : ndarray, shape (n,)
        Timestep durations in seconds (converted from My).
    temperatures : ndarray, shape (n,)
        Temperature in Kelvin at each timestep midpoint.
    """
    rng = np.random.default_rng(seed)
    Myr = 1.0e6 * 365.0 * 24.0 * 60.0 * 60.0

    # Timestep durations: uniform ~1 My steps
    dts = np.ones(n) * Myr + rng.uniform(-0.1, 0.1, n) * Myr

    # Temperature: starts at 120 °C, cools to 20 °C, in Kelvin
    T_celsius = np.linspace(120.0, 20.0, n) + rng.uniform(-2.0, 2.0, n)
    temperatures = T_celsius + 273.15

    return dts, temperatures


def make_simulate_aft_inputs(n=80, seed=7):
    """
    Build inputs suitable for simulate_AFT_annealing / _vectorized.

    Returns
    -------
    timesteps : ndarray, shape (n+1,)  in My
    temperature_input : ndarray, shape (n+1,)  in °C
    kinetic_value : float
    """
    rng = np.random.default_rng(seed)

    # timesteps in My (must be monotonically increasing)
    timesteps = np.sort(rng.uniform(0, 100, n + 1))
    timesteps[0] = 0.0

    # temperature: 100 → 20 °C with small noise, resampling will handle steps >3.5°
    T = np.linspace(100.0, 20.0, n + 1) + rng.uniform(-1.0, 1.0, n + 1)

    kinetic_value = 0.0  # Cl wt fraction = 0 (fluorapatite)

    return timesteps, T, kinetic_value


def make_he_diffusion_inputs(n=100, seed=13):
    """
    Build inputs for He_diffusion_Meesters_and_Dunai_2002 / _vectorized.

    Returns
    -------
    t : ndarray, shape (n,) — time in seconds (monotonically increasing)
    D : ndarray, shape (n,) — diffusivity in m² s⁻¹
    radius : float — grain radius in m
    Ur0 : float — effective radiogenic production rate
    """
    rng = np.random.default_rng(seed)

    Myr = 1.0e6 * 365.0 * 24.0 * 60.0 * 60.0
    t = np.linspace(0.0, 100.0 * Myr, n)

    # Diffusivity increases with temperature (grain cooled from 120 → 20 °C)
    T_K = np.linspace(393.15, 293.15, n)
    Ea = 32.9 * 4184.0   # J/mol
    R = 8.3144621
    D0 = 50.0 / 1e4       # m²/s (Farley 2000)
    D = D0 * np.exp(-Ea / (R * T_K))

    radius = 60e-6         # 60 µm grain radius
    Ur0 = 1e-18            # representative production rate

    return t, D, radius, Ur0


# ---------------------------------------------------------------------------
# Test 1: calculate_reduced_track_lengths vs _vectorized
# ---------------------------------------------------------------------------

class TestCalculateReducedTrackLengths:

    def test_identical_results_medium_history(self):
        """Vectorized and original produce identical rc arrays (n=200)."""
        dts, temperatures = make_thermal_history(n=200)

        rc_orig = AFTannealingLib.calculate_reduced_track_lengths(
            dts, temperatures)
        rc_vec = AFTannealingLib.calculate_reduced_track_lengths_vectorized(
            dts, temperatures)

        assert rc_orig.shape == rc_vec.shape, (
            f"Shape mismatch: orig {rc_orig.shape} vs vec {rc_vec.shape}")
        assert np.allclose(rc_orig, rc_vec), (
            f"Max abs diff: {np.abs(rc_orig - rc_vec).max():.3e}\n"
            f"orig[:5]={rc_orig[:5]}\nvec[:5]={rc_vec[:5]}")

    def test_identical_results_short_history(self):
        """Identical results for a short (n=20) thermal history."""
        dts, temperatures = make_thermal_history(n=20, seed=99)

        rc_orig = AFTannealingLib.calculate_reduced_track_lengths(
            dts, temperatures)
        rc_vec = AFTannealingLib.calculate_reduced_track_lengths_vectorized(
            dts, temperatures)

        assert np.allclose(rc_orig, rc_vec), (
            f"Max abs diff: {np.abs(rc_orig - rc_vec).max():.3e}")

    def test_identical_results_constant_temperature(self):
        """Isothermal history: vectorized matches original."""
        Myr = 1.0e6 * 365.0 * 24.0 * 60.0 * 60.0
        n = 50
        dts = np.ones(n) * Myr
        temperatures = np.full(n, 80.0 + 273.15)  # constant 80 °C in K

        rc_orig = AFTannealingLib.calculate_reduced_track_lengths(
            dts, temperatures)
        rc_vec = AFTannealingLib.calculate_reduced_track_lengths_vectorized(
            dts, temperatures)

        assert np.allclose(rc_orig, rc_vec), (
            f"Max abs diff: {np.abs(rc_orig - rc_vec).max():.3e}")

    def test_non_default_annealing_parameters(self):
        """Custom C0/C1/C2/C3/alpha values: vectorized still matches original."""
        dts, temperatures = make_thermal_history(n=60, seed=3)

        kwargs = dict(alpha=0.05, C0=0.4, C1=0.012, C2=-60.0, C3=-8.0)
        rc_orig = AFTannealingLib.calculate_reduced_track_lengths(
            dts, temperatures, **kwargs)
        rc_vec = AFTannealingLib.calculate_reduced_track_lengths_vectorized(
            dts, temperatures, **kwargs)

        assert np.allclose(rc_orig, rc_vec), (
            f"Max abs diff: {np.abs(rc_orig - rc_vec).max():.3e}")


# ---------------------------------------------------------------------------
# Test 2: simulate_AFT_annealing vs _vectorized
# ---------------------------------------------------------------------------

class TestSimulateAFTAnnealing:

    def _compare(self, timesteps, temperature_input, kinetic_value, **kwargs):
        """Helper: run both functions, compare all 8 return values."""
        result_orig = AFTannealingLib.simulate_AFT_annealing(
            timesteps, temperature_input, kinetic_value,
            use_fortran_algorithm=False, **kwargs)
        result_vec = AFTannealingLib.simulate_AFT_annealing_vectorized(
            timesteps, temperature_input, kinetic_value,
            use_fortran_algorithm=False, **kwargs)

        names = ['track_length_pdf', 'aft_age_myr', 'l_mean', 'l_mean_std',
                 'rm', 'rc', 'rho_age', 'dt']
        for name, v_orig, v_vec in zip(names, result_orig, result_vec):
            if np.isscalar(v_orig) or (hasattr(v_orig, 'ndim') and v_orig.ndim == 0):
                # Both NaN → identical degenerate case, passes
                if np.isnan(float(v_orig)) and np.isnan(float(v_vec)):
                    continue
                assert np.isclose(v_orig, v_vec), (
                    f"{name}: orig={v_orig:.6g}, vec={v_vec:.6g}")
            else:
                a = np.asarray(v_orig)
                b = np.asarray(v_vec)
                # NaN positions must match; non-NaN values must be close
                nan_orig = np.isnan(a)
                nan_vec = np.isnan(b)
                assert np.array_equal(nan_orig, nan_vec), (
                    f"{name}: NaN pattern differs between original and vectorized")
                assert np.allclose(a[~nan_orig], b[~nan_vec]), (
                    f"{name}: max abs diff = {np.abs(a[~nan_orig] - b[~nan_vec]).max():.3e}")

    def test_fluorapatite_clwt(self):
        """Default Cl wt fraction = 0 (fluorapatite), Clwt kinetic parameter."""
        timesteps, T, kv = make_simulate_aft_inputs(n=80)
        self._compare(timesteps, T, kv, kinetic_parameter='Clwt')

    def test_chlorapatite_clwt(self):
        """Higher Cl content (0.3 wt%), different kinetics."""
        timesteps, T, _ = make_simulate_aft_inputs(n=60, seed=17)
        self._compare(timesteps, T, 0.3, kinetic_parameter='Clwt')

    def test_rmr0_kinetic_parameter(self):
        """Using rmr0 directly as kinetic parameter."""
        timesteps, T, _ = make_simulate_aft_inputs(n=50, seed=55)
        self._compare(timesteps, T, 0.5, kinetic_parameter='rmr0')

    def test_c_axis_correction(self):
        """With c-axis correction enabled."""
        timesteps, T, kv = make_simulate_aft_inputs(n=60, seed=22)
        self._compare(timesteps, T, kv,
                      kinetic_parameter='Clwt', apply_c_axis_correction=True)


# ---------------------------------------------------------------------------
# Test 3: He_diffusion_Meesters_and_Dunai_2002 vs _vectorized
# ---------------------------------------------------------------------------

class TestHeDiffusionVectorized:

    def _compare_he(self, t, D, radius, Ur0, **kwargs):
        """Helper: run both functions and compare returned Cav/Ur0 array."""
        t_c_orig = He_diffusion_Meesters_and_Dunai_2002(
            t, D, radius, Ur0, **kwargs)
        t_c_vec = He_diffusion_Meesters_and_Dunai_2002_vectorized(
            t, D, radius, Ur0, **kwargs)

        assert t_c_orig.shape == t_c_vec.shape, (
            f"Shape mismatch: orig {t_c_orig.shape} vs vec {t_c_vec.shape}")
        assert np.allclose(t_c_orig, t_c_vec), (
            f"Max abs diff: {np.abs(t_c_orig - t_c_vec).max():.3e}\n"
            f"Max rel diff: {np.abs((t_c_orig - t_c_vec) / (t_c_orig + 1e-30)).max():.3e}")

    def test_all_timesteps_sphere_constant_U(self):
        """all_timesteps=True, sphere shape, constant U function."""
        t, D, radius, Ur0 = make_he_diffusion_inputs(n=80)
        self._compare_he(t, D, radius, Ur0,
                         U_function='constant',
                         all_timesteps=True,
                         alpha_ejection=False,
                         n_eigenmodes=15)

    def test_all_timesteps_with_alpha_ejection(self):
        """all_timesteps=True with alpha ejection correction."""
        t, D, radius, Ur0 = make_he_diffusion_inputs(n=80)
        self._compare_he(t, D, radius, Ur0,
                         U_function='constant',
                         all_timesteps=True,
                         alpha_ejection=True,
                         stopping_distance=20e-6,
                         n_eigenmodes=15)

    def test_exponential_U_function(self):
        """U_function='exponential' with all_timesteps=True."""
        t, D, radius, Ur0 = make_he_diffusion_inputs(n=60, seed=99)
        self._compare_he(t, D, radius, Ur0,
                         U_function='exponential',
                         all_timesteps=True,
                         alpha_ejection=False,
                         n_eigenmodes=10)

    def test_fewer_eigenmodes(self):
        """n_eigenmodes=5 (fast): vectorized matches original."""
        t, D, radius, Ur0 = make_he_diffusion_inputs(n=50, seed=77)
        self._compare_he(t, D, radius, Ur0,
                         U_function='constant',
                         all_timesteps=True,
                         alpha_ejection=False,
                         n_eigenmodes=5)


# ---------------------------------------------------------------------------
# Synthetic VR thermal history helpers
# ---------------------------------------------------------------------------

def make_vr_thermal_history(n=500, seed=42):
    """
    Build a realistic synthetic thermal history for VR testing.

    Returns
    -------
    times : ndarray, shape (n,)
        Time in Ma, monotonically decreasing (oldest first → 0 = present).
    temperatures : ndarray, shape (n,)
        Temperature in °C at each timestep.
    """
    rng = np.random.default_rng(seed)
    # times: from ~n*0.5 Ma down to 0, with small jitter so spacing is uneven
    times_raw = np.linspace(n * 0.5, 0.0, n) + rng.uniform(-0.05, 0.05, n)
    times_raw = np.clip(times_raw, 0.0, None)
    times_raw[-1] = 0.0
    times = np.sort(times_raw)[::-1].copy()  # largest → smallest (oldest first)

    # Temperature: starts at 20 °C, heats to 120 °C, then cools back to 20 °C
    half = n // 2
    T = np.concatenate([
        np.linspace(20.0, 120.0, half) + rng.uniform(-1.0, 1.0, half),
        np.linspace(120.0, 20.0, n - half) + rng.uniform(-1.0, 1.0, n - half),
    ])
    return times, T


# ---------------------------------------------------------------------------
# Test 4: easyRo vs easyRo_vectorized
# ---------------------------------------------------------------------------

class TestEasyRoVectorized:
    """Correctness tests: easyRo_vectorized must produce bit-identical results to easyRo."""

    def _compare(self, times, temperatures, vr_method='easyRo'):
        Ro_orig = easyRo.easyRo(times, temperatures, vr_method=vr_method)
        Ro_vec = easyRo.easyRo_vectorized(times, temperatures, vr_method=vr_method)

        assert Ro_orig.shape == Ro_vec.shape, (
            f"Shape mismatch: orig {Ro_orig.shape} vs vec {Ro_vec.shape}")
        assert np.allclose(Ro_orig, Ro_vec), (
            f"vr_method={vr_method!r}: max abs diff = "
            f"{np.abs(Ro_orig - Ro_vec).max():.3e}\n"
            f"orig[:5]={Ro_orig[:5]}\nvec[:5]={Ro_vec[:5]}")

    def test_easyro_burial_and_exhumation(self):
        """Heat-up then cool-down history: vectorized matches original."""
        times, T = make_vr_thermal_history(n=200, seed=1)
        self._compare(times, T, vr_method='easyRo')

    def test_easyro_constant_temperature(self):
        """Isothermal edge case (no heating rate variation)."""
        n = 100
        times = np.linspace(50.0, 0.0, n)
        T = np.full(n, 80.0)
        self._compare(times, T, vr_method='easyRo')

    def test_basinro_method(self):
        """basinRo method: vectorized matches original."""
        times, T = make_vr_thermal_history(n=300, seed=42)
        self._compare(times, T, vr_method='basinRo')

    def test_debug_mode_intermediate_arrays(self):
        """In debug mode, all intermediate arrays are identical."""
        times, T = make_vr_thermal_history(n=100, seed=5)

        out_orig = easyRo.easyRo(times, T, vr_method='easyRo', debug=True)
        out_vec = easyRo.easyRo_vectorized(times, T, vr_method='easyRo', debug=True)

        names = ['Ro', 'sumReacted', 'cumulativeReacted', 'deltaI', 'I', 'EdivRT']
        for name, v_orig, v_vec in zip(names, out_orig, out_vec):
            assert np.allclose(v_orig, v_vec), (
                f"debug array '{name}': max abs diff = "
                f"{np.abs(np.asarray(v_orig) - np.asarray(v_vec)).max():.3e}")
