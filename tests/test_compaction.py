"""
Validation tests for the compaction-driven excess pore pressure solver.

Both tests build a synthetic sediment column that grows one layer at a
time, using the same node-activation logic as run_burial_hist_model in
lib/pybasin_lib.py, and drive it with the same loading-source helper
(cumulative_buoyant_stress) and pore pressure solver
(solve_1D_pore_pressure) that run_burial_hist_model uses. This avoids
needing a full example dataset while still exercising the actual
production loading and solver code.

test_compaction_envelope checks that the excess pressure and porosity
stay within their expected physical bounds for a moderately permeable
column. test_gibson_undrained_limit checks a quantitative limit: for a
column with zero permeability, excess pressure should match the buoyant
weight of the solid grains scaled by the undrained loading efficiency,
which approaches the classic Gibson (1958) result as the elastic
specific storage becomes small compared with the compaction storage.
Both of these use run_growing_column() above, in which each node is
loaded once, at its own creation step; this gives the correct answer in
the zero-permeability limit these two tests use, because that limit is
path-independent (only the final buoyant stress matters, not the
loading history), but it is not a valid driver for a finite-permeability
transient test, where the loading history does matter.
test_gibson_transient_solution instead uses run_growing_column_gradual(),
in which every existing node is buried by one grid spacing at every
step, giving every node the continuous, constant loading rate Gibson's
equation assumes, and compares the finite-permeability transient
solution against Gibson's (1958) exact solution (his Eq. 2, for a layer
growing at a constant rate on an impervious base), evaluated by
numerical quadrature. This is the transient, finite-permeability
regime that the other two tests do not exercise.
test_compaction_storage_irreversibility is a direct unit test of
compaction_storage(): it checks that the returned specific storage is
the elastic form beta_matrix + porosity * beta_water in both loading
regimes, and that the virgin/elastic flag flips once effective stress
drops below its historical maximum, which is what switches off the
inelastic compaction source term in run_burial_hist_model().
test_compaction_source_and_storage_units checks that the porosity
change rate, the compaction fluid source and the specific storage are
dimensionally consistent, ie. that every term of the excess pressure
equation reduces to inverse seconds.

Run with:
    pytest tests/test_compaction.py
"""

import os
import sys

import numpy as np
from scipy.integrate import quad

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import lib.pybasin_lib as pybasin_lib


YEAR = 365.25 * 24 * 60 * 60.0


def run_growing_column(n_layers, dz, dt, phi_fn, permeability, viscosity,
                        rho_grain, rho_f, compressibility_stress,
                        alpha_skeleton=1.0e-9, beta_water=4.4e-10):
    """
    simulate a sediment column that grows by one layer of thickness dz
    every timestep dt, and solve for excess pore pressure at each step.

    node identities are fixed integers 0 .. n_layers, with node
    n_layers always the first node deposited and node 0 the last.
    at a given step, active nodes are ordered from shallow (the node
    just deposited, depth 0) to deep, matching the node ordering
    expected by solve_1D_pore_pressure.

    the elastic specific storage is computed with compaction_storage()
    and the inelastic compaction is added as a porosity change source
    term, exactly as run_burial_hist_model() does, so this harness
    exercises the actual production storage, source and irreversibility
    logic. alpha_skeleton and beta_water default to the same small
    values as run_burial_hist_model(); leaving them at exactly 0 is not
    realistic (and not what production code does), since the elastic
    specific storage would then be 0 wherever the compaction source
    switches off (the elastic/unloading regime), making the pore
    pressure solve singular there.

    the historical maximum effective stress used for the
    virgin/elastic decision is ratcheted using the excess pressure
    just solved for, not the pre-solve estimate used for the decision
    itself, matching run_burial_hist_model(): a freshly deposited
    node's pre-solve estimate always assumes 0 excess pressure, which
    would otherwise overstate its true effective stress and lock in a
    ceiling it never actually reached.

    returns final depth, excess pressure and porosity, one value per
    node, in the fixed node order (node 0 first).
    """
    n_nodes = n_layers + 1
    P_ex = np.zeros(n_nodes)
    sigma_buoyant_prev = np.zeros(n_nodes)
    sigma_eff_max_prev = np.zeros(n_nodes)

    for i in range(1, n_layers + 1):
        active_ids = np.arange(n_layers, n_layers - i - 1, -1)
        depths = (n_layers - active_ids) * dz

        phi = phi_fn(depths)
        rho_bulk = rho_grain * (1.0 - phi) + rho_f * phi

        sigma_buoyant = pybasin_lib.cumulative_buoyant_stress(
            depths, rho_bulk, rho_f)

        loading_rate = (sigma_buoyant - sigma_buoyant_prev[active_ids]) / dt

        P_init = P_ex[active_ids]
        sigma_eff = np.clip(sigma_buoyant - P_init, 0.0, None)
        sigma_eff_max_prev_active = sigma_eff_max_prev[active_ids]
        c_sigma = np.full(i + 1, compressibility_stress)

        storage_elastic, is_virgin = pybasin_lib.compaction_storage(
            phi, sigma_eff, sigma_eff_max_prev_active,
            alpha_skeleton, beta_water)

        dphi_dsigma_eff = np.where(is_virgin, -c_sigma * phi, 0.0)

        Q = -dphi_dsigma_eff * loading_rate
        storage = storage_elastic - dphi_dsigma_eff

        K_hydraulic = np.full(i, permeability / viscosity)

        P_new, _ = pybasin_lib.solve_1D_pore_pressure(
            P_init, depths, dt, K_hydraulic, storage, Q,
            None, 0.0, 0.0, None)

        P_ex[active_ids] = P_new
        sigma_buoyant_prev[active_ids] = sigma_buoyant

        sigma_eff_post = np.clip(sigma_buoyant - P_new, 0.0, None)
        sigma_eff_max_prev[active_ids] = np.maximum(
            sigma_eff_post, sigma_eff_max_prev_active)

    depths_final = (n_layers - np.arange(n_nodes)) * dz

    return depths_final, P_ex, phi_fn(depths_final)


def test_compaction_envelope():

    n_layers = 80
    total_depth = 1000.0
    dz = total_depth / n_layers
    total_time = 10.0e6 * YEAR
    dt = total_time / n_layers
    phi_o = 0.5
    rho_grain = 2500.0
    rho_f = 1000.0
    compaction_coeff = 1.0e-3
    compressibility_stress = 5.0e-8
    permeability = 1.0e-20
    viscosity = 1.0e-3

    def phi_fn(z):
        return phi_o * np.exp(-compaction_coeff * z)

    depths, P_ex, phi = run_growing_column(
        n_layers, dz, dt, phi_fn, permeability, viscosity,
        rho_grain, rho_f, compressibility_stress)

    order = np.argsort(depths)
    depths_s = depths[order]
    P_ex_s = P_ex[order]
    phi_s = phi[order]

    assert not np.isnan(P_ex_s).any()

    assert abs(P_ex_s[0]) < 1.0

    assert np.all(P_ex_s >= -1.0)

    ceiling = (rho_grain - rho_f) * pybasin_lib.G_ACCEL * depths_s
    assert np.all(P_ex_s <= ceiling + 1.0)

    deep = depths_s > 0.5 * depths_s.max()
    trend = np.polyfit(depths_s[deep], P_ex_s[deep], 1)[0]
    assert trend > 0.0

    assert phi_s.max() <= phi_o + 1.0e-6
    assert phi_s[-1] < phi_o


def test_gibson_undrained_limit():

    n_layers = 60
    total_depth = 1000.0
    dz = total_depth / n_layers
    total_time = 10.0e6 * YEAR
    dt = total_time / n_layers
    phi_o = 0.5
    rho_grain = 2500.0
    rho_f = 1000.0
    compressibility_stress = 5.0e-8
    alpha_skeleton = 1.0e-9
    beta_water = 4.4e-10
    permeability = 0.0
    viscosity = 1.0e-3

    def phi_fn(z):
        return np.full_like(z, phi_o)

    depths, P_ex, _ = run_growing_column(
        n_layers, dz, dt, phi_fn, permeability, viscosity,
        rho_grain, rho_f, compressibility_stress,
        alpha_skeleton=alpha_skeleton, beta_water=beta_water)

    order = np.argsort(depths)
    depths_s = depths[order]
    P_ex_s = P_ex[order]

    # undrained loading efficiency: with a finite elastic specific
    # storage the excess pressure carries a fraction of the buoyant
    # load, not all of it. that fraction is the ratio of the compaction
    # (inelastic) storage to the total storage on the solver, ie.
    # c_sigma * phi / (c_sigma * phi + beta_matrix + phi * beta_water).
    # in the limit beta_matrix, beta_water << c_sigma this approaches 1
    # and recovers the classic Gibson (1958) result
    inelastic_storage = compressibility_stress * phi_o
    elastic_storage = alpha_skeleton + phi_o * beta_water
    loading_efficiency = inelastic_storage / (inelastic_storage
                                              + elastic_storage)

    buoyant_slope = (1.0 - phi_o) * (rho_grain - rho_f) * pybasin_lib.G_ACCEL
    expected_slope = loading_efficiency * buoyant_slope

    grown = depths_s > 0.05 * depths_s.max()
    slope, intercept = np.polyfit(depths_s[grown], P_ex_s[grown], 1)
    resid = P_ex_s[grown] - (slope * depths_s[grown] + intercept)
    r_squared = 1.0 - np.var(resid) / np.var(P_ex_s[grown])

    assert abs(slope / expected_slope - 1.0) < 0.01
    assert r_squared > 0.999


def run_growing_column_gradual(n_layers, dz, dt, phi_o, permeability,
                               viscosity, rho_grain, rho_f,
                               compressibility_stress):
    """
    Simulate a sediment column that grows by depositing one node at the
    surface every timestep dt, burying every previously deposited node by
    one grid spacing dz. Unlike run_growing_column() above, this gives
    every active node a constant, continuous loading rate as the column
    grows, matching the assumption behind Gibson's (1958) equation. It is
    needed for a finite-permeability transient benchmark, where the
    loading history matters and run_growing_column()'s single, one-off
    load applied at each node's own creation step is not a valid driver.

    The elastic skeleton and pore water compressibilities are fixed at
    zero, so that compaction_storage() returns a purely compaction-driven
    storage term, ie. the single linear compressibility Gibson's theory
    assumes. Leaving them at nonzero production-like values would instead
    give an undrained loading efficiency below 1 (see
    test_gibson_undrained_limit above), which departs from Gibson's
    linear solution.

    Returns the depth below the current surface and the excess pressure
    of every node, in the same (arbitrary) order.
    """
    alpha_skeleton = 0.0
    beta_water = 0.0

    depths = np.array([0.0])
    P_ex = np.array([0.0])
    sigma_buoyant_prev = np.array([0.0])
    sigma_eff_max_prev = np.array([0.0])

    for _ in range(n_layers):

        depths = np.concatenate(([0.0], depths + dz))
        P_ex = np.concatenate(([0.0], P_ex))
        sigma_buoyant_prev = np.concatenate(([0.0], sigma_buoyant_prev))
        sigma_eff_max_prev = np.concatenate(([0.0], sigma_eff_max_prev))

        phi = np.full_like(depths, phi_o)
        rho_bulk = rho_grain * (1.0 - phi) + rho_f * phi

        sigma_buoyant = pybasin_lib.cumulative_buoyant_stress(
            depths, rho_bulk, rho_f)
        loading_rate = (sigma_buoyant - sigma_buoyant_prev) / dt

        sigma_eff = np.clip(sigma_buoyant - P_ex, 0.0, None)

        storage_elastic, is_virgin = pybasin_lib.compaction_storage(
            phi, sigma_eff, sigma_eff_max_prev, alpha_skeleton, beta_water)

        c_sigma = np.full_like(phi, compressibility_stress)
        dphi_dsigma_eff = np.where(is_virgin, -c_sigma * phi, 0.0)

        Q = -dphi_dsigma_eff * loading_rate
        storage = storage_elastic - dphi_dsigma_eff

        K_hydraulic = np.full(len(depths) - 1, permeability / viscosity)

        P_ex, _ = pybasin_lib.solve_1D_pore_pressure(
            P_ex, depths, dt, K_hydraulic, storage, Q,
            None, 0.0, 0.0, None)

        sigma_buoyant_prev = sigma_buoyant
        sigma_eff_post = np.clip(sigma_buoyant - P_ex, 0.0, None)
        sigma_eff_max_prev = np.maximum(sigma_eff_post, sigma_eff_max_prev)

    return depths, P_ex


def gibson_excess_pressure(z, t, m, cv, gamma_prime, xi_limit_factor=40.0):
    """
    Gibson's (1958) exact excess pore pressure solution for a layer
    growing at a constant rate m on an impervious base, with a pervious
    top, at elevation z above the base and time t. Evaluated by numerical
    quadrature over the free variable xi, using the image-source identity

        exp(-z^2/4cvt) cosh(z xi/2cvt) exp(-xi^2/4cvt)
            = 0.5 * [exp(-(z - xi)^2/4cvt) + exp(-(z + xi)^2/4cvt)]

    which keeps the integrand bounded (the literal cosh term overflows
    for large xi) and folds the leading exp(-z^2/4cvt) prefactor into the
    integral.
    """
    xi_limit = (xi_limit_factor * np.sqrt(4.0 * cv * t)
               + 20.0 * cv / m)

    def integrand(xi):
        tanh_term = np.tanh(m * xi / (2.0 * cv))
        image_pair = 0.5 * (np.exp(-(z - xi) ** 2 / (4.0 * cv * t))
                            + np.exp(-(z + xi) ** 2 / (4.0 * cv * t)))
        return xi * tanh_term * image_pair

    integral, _ = quad(integrand, 0.0, xi_limit, limit=400)
    prefactor = gamma_prime * (np.pi * cv * t) ** -0.5

    return gamma_prime * m * t - prefactor * integral


def test_gibson_transient_solution():

    n_layers = 60
    total_depth = 1000.0
    dz = total_depth / n_layers
    total_time = 10.0e6 * YEAR
    dt = total_time / n_layers
    m = dz / dt
    phi_o = 0.4
    rho_grain = 2650.0
    rho_f = 1000.0
    compressibility_stress = 1.0e-3
    viscosity = 1.0e-3

    gamma_prime = (1.0 - phi_o) * (rho_grain - rho_f) * pybasin_lib.G_ACCEL
    specific_storage = compressibility_stress * phi_o

    # permeability chosen so the consolidation time factor cv * t / H^2
    # is of order 1, ie. an intermediate regime where neither drainage
    # nor storage dominates, unlike the zero-permeability limit used by
    # test_gibson_undrained_limit above
    time_factor = 1.0
    cv = time_factor * total_depth ** 2 / total_time
    permeability = cv * specific_storage * viscosity

    depths, P_ex = run_growing_column_gradual(
        n_layers, dz, dt, phi_o, permeability, viscosity,
        rho_grain, rho_f, compressibility_stress)

    # Gibson's z is elevation above the impervious base; depths is depth
    # below the pervious surface, so z = total_depth - depth
    z_gibson = total_depth - depths
    P_ex_gibson = np.array([
        gibson_excess_pressure(z, total_time, m, cv, gamma_prime)
        for z in z_gibson])

    max_error = np.max(np.abs(P_ex - P_ex_gibson))
    max_relative_error = max_error / np.max(np.abs(P_ex_gibson))

    assert max_relative_error < 0.03


def test_compaction_storage_irreversibility():

    porosity = np.array([0.3, 0.3])
    beta_matrix = 1.0e-9
    beta_water = 4.4e-10

    sigma_eff_max_prev = np.array([1.0e7, 1.0e7])

    # the specific storage is now the elastic form beta_matrix +
    # porosity * beta_water and does not depend on the loading regime;
    # the regime instead switches the inelastic compaction source term
    # on or off in run_burial_hist_model()
    expected_storage = beta_matrix + porosity * beta_water

    # virgin case: effective stress at or above its historical maximum
    sigma_eff_virgin = np.array([1.2e7, 1.0e7])
    storage_virgin, is_virgin = pybasin_lib.compaction_storage(
        porosity, sigma_eff_virgin, sigma_eff_max_prev,
        beta_matrix, beta_water)

    assert np.all(is_virgin)
    assert np.allclose(storage_virgin, expected_storage)

    # unloading case: effective stress below its historical maximum
    # (eg. excess pressure has built up), so the node leaves the virgin
    # regime and its compaction source switches off, but the elastic
    # storage is unchanged
    sigma_eff_unload = np.array([0.5e7, 1.0e7])
    storage_unload, is_virgin_unload = pybasin_lib.compaction_storage(
        porosity, sigma_eff_unload, sigma_eff_max_prev,
        beta_matrix, beta_water)

    assert not is_virgin_unload[0]
    assert is_virgin_unload[1]
    assert np.allclose(storage_unload, expected_storage)


class Dimension:
    """
    Minimal physical dimension tracker for checking that the compaction
    source and storage equations are dimensionally consistent.

    A dimension is stored as exponents over the SI base units kilogram,
    metre and second. Multiplication and division add and subtract
    exponents; addition and subtraction require identical dimensions and
    otherwise raise, which is what turns the equations below into genuine
    unit consistency assertions.
    """

    def __init__(self, kg=0, m=0, s=0):
        self.kg = kg
        self.m = m
        self.s = s

    def __eq__(self, other):
        return (self.kg, self.m, self.s) == (other.kg, other.m, other.s)

    def __repr__(self):
        return "kg^%g m^%g s^%g" % (self.kg, self.m, self.s)

    def __mul__(self, other):
        if not isinstance(other, Dimension):
            return self
        return Dimension(self.kg + other.kg, self.m + other.m,
                         self.s + other.s)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return Dimension(self.kg - other.kg, self.m - other.m,
                         self.s - other.s)

    def __neg__(self):
        return self

    def __add__(self, other):
        if self != other:
            raise ValueError("cannot add %r and %r" % (self, other))
        return self

    def __sub__(self, other):
        return self.__add__(other)


DIMENSIONLESS = Dimension()
PASCAL = Dimension(kg=1, m=-1, s=-2)
PER_PASCAL = Dimension(kg=-1, m=1, s=2)
PER_SECOND = Dimension(s=-1)
METRE = Dimension(m=1)
SECOND = Dimension(s=1)


def test_compaction_source_and_storage_units():

    # sanity check on the dimension helper itself
    assert PASCAL * PER_PASCAL == DIMENSIONLESS
    assert DIMENSIONLESS / PASCAL == PER_PASCAL

    # dimensions of the physical inputs
    c_sigma = PER_PASCAL
    phi = DIMENSIONLESS
    beta_matrix = PER_PASCAL
    beta_water = PER_PASCAL
    loading_rate = PASCAL / SECOND
    excess_pressure = PASCAL
    depth = METRE
    dt = SECOND
    # hydraulic conductivity is permeability / viscosity, m^2 (Pa s)^-1
    hydraulic_conductivity = Dimension(m=2) / (PASCAL * SECOND)

    # porosity law phi = phi0 * exp(-c_sigma * sigma_eff), so
    # dphi/dsigma_eff = -c_sigma * phi is a compressibility (Pa^-1)
    dphi_dsigma_eff = -c_sigma * phi
    assert dphi_dsigma_eff == PER_PASCAL

    # compaction source Q = -dphi/dsigma_eff * loading_rate expels pore
    # water at a volumetric rate per bulk volume (s^-1)
    Q_loading = -dphi_dsigma_eff * loading_rate
    assert Q_loading == PER_SECOND

    # elastic specific storage Ss = beta_matrix + phi * beta_water; the
    # addition only succeeds if both terms are Pa^-1
    storage_elastic = beta_matrix + phi * beta_water
    assert storage_elastic == PER_PASCAL

    # storage on the solver, elastic plus the implicit -dphi/dsigma_eff
    # term; the subtraction only succeeds if both terms are Pa^-1
    storage_active = storage_elastic - dphi_dsigma_eff
    assert storage_active == PER_PASCAL

    # every term of storage * dP/dt = d/dz(K dP/dz) + Q must be s^-1
    accumulation = storage_active * (excess_pressure / dt)
    diffusion = hydraulic_conductivity * (excess_pressure / depth) / depth
    assert accumulation == PER_SECOND
    assert diffusion == PER_SECOND
    assert Q_loading == PER_SECOND
