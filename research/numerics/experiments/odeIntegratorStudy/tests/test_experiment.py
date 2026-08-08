"""Validation tests for the ODE integrator comparison experiment.

These tests must pass before any results are trusted (laboratory rule 3).
"""

import numpy as np
import pytest

from experiment import (
    H_ENERGY,
    METHODS,
    STEPSIZES,
    T_PERIOD,
    THEORETICAL_ORDER,
    X0,
    V0,
    analytic_solution,
    convergence_table,
    energy_study,
    integrate,
    total_energy,
)


def test_analytic_solution_satisfies_ic():
    assert analytic_solution(0.0) == pytest.approx(X0)
    # v(0) = -omega x0 sin(0) + v0 cos(0) = V0; check numerical derivative
    t = 1e-6
    d = (analytic_solution(t) - analytic_solution(-t)) / (2 * t)
    assert d == pytest.approx(V0, rel=1e-6)


def test_analytic_solution_energy_const():
    t = np.linspace(0.0, 3 * T_PERIOD, 1000)
    e = total_energy(analytic_solution(t), -X0 * np.sin(t))
    assert np.allclose(e, e[0], atol=1e-12)


@pytest.mark.parametrize("h", [0.05, 0.025])
def test_all_methods_close_to_analytic_at_coarse_step(h):
    states = integrate("rk4", h, T_PERIOD)
    x = states[:, 0]
    t = np.arange(len(x)) * h
    assert np.allclose(x, analytic_solution(t), atol=0.1)


@pytest.mark.parametrize("m", METHODS)
def test_measured_order_close_to_theoretical(m):
    errors = []
    for h in STEPSIZES:
        states = integrate(m, h, T_PERIOD)
        t = np.arange(len(states)) * h
        err = np.max(np.abs(states[:, 0] - analytic_solution(t)))
        errors.append(err)
    slope = (np.log(errors[0]) - np.log(errors[-1])) / (
        np.log(STEPSIZES[0]) - np.log(STEPSIZES[-1])
    )
    expected = THEORETICAL_ORDER[m]
    assert abs(slope - expected) < 0.5, (
        f"{m}: measured order {slope:.2f}, theoretical {expected}"
    )


def _max_energy_deviation(method: str, t_end: float) -> float:
    states = integrate(method, H_ENERGY, t_end)
    e = total_energy(states[:, 0], states[:, 1])
    return float(np.max(np.abs(e - total_energy(X0, V0))))


def test_symplectic_energy_error_is_bounded():
    """Verlet (symplectic) energy error must not grow with integration time."""
    short = _max_energy_deviation("velocity_verlet", 40 * np.pi)
    long = _max_energy_deviation("velocity_verlet", 400 * np.pi)
    # Bounded, non-secularly-growing error: ratio stays near 1 (oscillatory
    # in time), far below linear-in-horizon scaling.
    assert long <= 4.0 * short, f"bounded lost: short={short:.3e}, long={long:.3e}"


def test_nonsymplectic_energy_error_grows_secularly():
    """Euler (non-symplectic) energy error grows linearly with horizon."""
    short = _max_energy_deviation("euler", 40 * np.pi)
    long = _max_energy_deviation("euler", 400 * np.pi)
    assert long > 5.0 * short, f"secular growth: short={short:.3e}, long={long:.3e}"


def test_energy_drift_symplectic_vs_dissipative():
    study = energy_study()
    dev = study["energy_deviations"]
    # Symplectic methods must preserve energy far better than non-symplectic.
    assert dev["velocity_verlet"]["max_deviation"] < dev["euler"]["max_deviation"]
    assert dev["symplectic_euler"]["max_deviation"] < dev["euler"]["max_deviation"]


def test_determinism():
    a = integrate("rk4", 0.05, T_PERIOD)
    b = integrate("rk4", 0.05, T_PERIOD)
    assert np.array_equal(a, b)