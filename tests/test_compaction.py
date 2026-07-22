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
column with (near) zero permeability, excess pressure should match the
analytic Gibson (1958) undrained solution, ie. the buoyant weight of
the solid grains. test_compaction_storage_irreversibility is a direct
unit test of compaction_storage(): the other two tests only exercise
its virgin (monotonic loading) branch, since neither ever unloads, so
this test checks the elastic-only (unloading) branch separately, where
storage drops to the elastic value because effective stress is below
its historical maximum.

Run with:
    pytest tests/test_compaction.py
"""

import os
import sys

import numpy as np

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

    storage is computed with compaction_storage(), the same function
    run_burial_hist_model() uses, so this harness exercises the actual
    production storage and irreversibility logic. alpha_skeleton and
    beta_water default to the same small values as
    run_burial_hist_model(); leaving them at exactly 0 is not
    realistic (and not what production code does), since the lagged
    effective stress estimate used here and in run_burial_hist_model()
    can dip into the elastic branch briefly even during otherwise
    monotonic loading, which would make storage exactly 0 and the
    pore pressure solve singular if alpha_skeleton and beta_water
    were both 0 too.

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

        storage, _ = pybasin_lib.compaction_storage(
            phi, sigma_eff, sigma_eff_max_prev_active, c_sigma,
            alpha_skeleton, beta_water)

        Q = storage * loading_rate

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
    permeability = 0.0
    viscosity = 1.0e-3

    def phi_fn(z):
        return np.full_like(z, phi_o)

    depths, P_ex, _ = run_growing_column(
        n_layers, dz, dt, phi_fn, permeability, viscosity,
        rho_grain, rho_f, compressibility_stress)

    order = np.argsort(depths)
    depths_s = depths[order]
    P_ex_s = P_ex[order]

    expected_slope = (1.0 - phi_o) * (rho_grain - rho_f) * pybasin_lib.G_ACCEL

    grown = depths_s > 0.05 * depths_s.max()
    slope, intercept = np.polyfit(depths_s[grown], P_ex_s[grown], 1)
    resid = P_ex_s[grown] - (slope * depths_s[grown] + intercept)
    r_squared = 1.0 - np.var(resid) / np.var(P_ex_s[grown])

    assert abs(slope / expected_slope - 1.0) < 0.01
    assert r_squared > 0.999


def test_compaction_storage_irreversibility():

    porosity = np.array([0.3, 0.3])
    compressibility_stress = np.array([5.0e-8, 5.0e-8])
    alpha_skeleton = 1.0e-9
    beta_water = 4.4e-10

    sigma_eff_max_prev = np.array([1.0e7, 1.0e7])

    # virgin case: effective stress at or above its historical
    # maximum, so the irreversible compaction coefficient contributes
    # to storage
    sigma_eff_virgin = np.array([1.2e7, 1.0e7])
    storage_virgin, is_virgin = pybasin_lib.compaction_storage(
        porosity, sigma_eff_virgin, sigma_eff_max_prev,
        compressibility_stress, alpha_skeleton, beta_water)

    assert np.all(is_virgin)
    expected_storage_virgin = porosity * (
        beta_water + alpha_skeleton + compressibility_stress)
    assert np.allclose(storage_virgin, expected_storage_virgin)

    # unloading case: effective stress below its historical maximum
    # (eg. excess pressure has built up), so storage is elastic only
    # (much smaller)
    sigma_eff_unload = np.array([0.5e7, 1.0e7])
    storage_unload, is_virgin_unload = pybasin_lib.compaction_storage(
        porosity, sigma_eff_unload, sigma_eff_max_prev,
        compressibility_stress, alpha_skeleton, beta_water)

    assert not is_virgin_unload[0]
    assert is_virgin_unload[1]
    expected_storage_unload_0 = porosity[0] * (beta_water + alpha_skeleton)
    assert np.isclose(storage_unload[0], expected_storage_unload_0)
    assert storage_unload[0] < storage_virgin[0]
