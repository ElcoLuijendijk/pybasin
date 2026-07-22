"""
Validation test for the 1D conductive heat flow solver.

test_conduction_step_change_vs_analytical drives the production solver
solve_1D_heat_flow (lib/pybasin_lib.py) with a step temperature change
imposed at the top boundary of an initially isothermal, homogeneous
column, and compares the resulting temperature profiles against the
analytical solution for a semi infinite medium (Carslaw and Jaeger,
1959):

    dT(x, t) = dT0 * erfc(x / (2 sqrt(kappa t)))

with kappa = K / (rho c) the thermal diffusivity and erfc the
complementary error function. The domain is made long enough that the
far boundary is never reached over the simulated time, so it behaves as
semi infinite.

This mirrors the analytical check in
tests/analytical_solution_heat_flow.ipynb, but turns the eyeball
comparison into an automated regression guard on the model's central
conductive solver, analogous to how test_gibson_undrained_limit guards
the compaction solver.

Run with:
    pytest tests/test_heat_flow.py
"""

import os
import sys

import numpy as np
from scipy.special import erfc

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import lib.pybasin_lib as pybasin_lib  # noqa: E402


YEAR = 365.25 * 24 * 3600.0


def run_conduction_step_change(kappa, dT0, z, dt, n_steps, output_steps):
    """
    march solve_1D_heat_flow forward for a top boundary step change.

    the column starts isothermal at 0, with a fixed temperature dT0
    imposed at the top node and a fixed far field temperature of 0 at
    the base. thermal conductivity, density and heat capacity are chosen
    so that K / (rho c) equals the requested diffusivity kappa.

    returns a dict mapping each requested step index to a copy of the
    temperature profile at that step.
    """
    n_nodes = len(z)

    K_value = 2.5
    rho_value = 2500.0
    c_value = K_value / (kappa * rho_value)

    rho = np.ones(n_nodes) * rho_value
    c = np.ones(n_nodes) * c_value
    Q = np.zeros(n_nodes)

    T = np.zeros(n_nodes)
    T[0] = dT0

    A = None
    profiles = {}

    for step in range(1, n_steps + 1):
        T, A = pybasin_lib.solve_1D_heat_flow(
            T, z, dt, K_value, rho, c, Q,
            upper_bnd_flux=None,
            lower_bnd_flux=None,
            fixed_upper_temperature=dT0,
            fixed_lower_temperature=0.0,
            A=A)

        if step in output_steps:
            profiles[step] = T.copy()

    return profiles


def test_conduction_step_change_vs_analytical():

    kappa = 1.5e-6
    dT0 = 100.0

    total_depth = 40000.0
    dz = 100.0
    z = np.arange(0.0, total_depth + dz, dz)

    dt = 2000.0 * YEAR
    total_time = 2.0e6 * YEAR
    n_steps = int(round(total_time / dt))

    output_times_my = [0.5, 1.0, 2.0]
    output_steps = [
        int(round(t_my * 1.0e6 * YEAR / dt)) for t_my in output_times_my
    ]

    profiles = run_conduction_step_change(
        kappa, dT0, z, dt, n_steps, output_steps)

    max_error = 0.0
    for step in output_steps:
        t_seconds = step * dt
        analytical = dT0 * erfc(z / (2.0 * np.sqrt(kappa * t_seconds)))
        numerical = profiles[step]

        assert not np.isnan(numerical).any()

        assert np.isclose(numerical[0], dT0)

        assert np.all(np.diff(numerical) <= 1.0e-9)

        max_error = max(max_error, np.max(np.abs(numerical - analytical)))

    assert max_error < 0.5
