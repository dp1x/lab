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

# Per-method coarse-step tolerance: at least ~10x the error each method actually
# produces at h = 0.05 (euler 1.7e-01, rk2 2.0e-03, rk4 2.5e-07,
# symplectic_euler 2.5e-02, velocity_verlet 5.0e-04). The test catches gross
# implementation errors without demanding precision typical of fine grids.
METHOD_ATOL = {
    "euler": 0.5,
    "rk2_midpoint": 0.02,
    "rk4": 1e-5,
    "symplectic_euler": 0.05,
    "velocity_verlet": 5e-3,
}


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


@pytest.mark.parametrize("m", METHODS)
@pytest.mark.parametrize("h", [0.05, 0.025])
def test_all_methods_close_to_analytic_at_coarse_step(m, h):
    states = integrate(m, h, T_PERIOD)
    x = states[:, 0]
    # Grid ends exactly on T_PERIOD (see grid_h in experiment.py).
    t = np.arange(len(x)) * (T_PERIOD / (len(x) - 1))
    err = float(np.max(np.abs(x - analytic_solution(t))))
    assert err < METHOD_ATOL[m], (
        f"{m} at h={h}: max error {err:.2e} exceeds tolerance {METHOD_ATOL[m]:.1e} "
        f"(expected O(h^{THEORETICAL_ORDER[m]}) at this scale)"
    )


@pytest.mark.parametrize("m", METHODS)
def test_measured_order_close_to_theoretical(m):
    errors = []
    for h in STEPSIZES:
        states = integrate(m, h, T_PERIOD)
        # Reference on the exact integrator grid (ends on T_PERIOD).
        t = np.arange(len(states)) * (T_PERIOD / (len(states) - 1))
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


def test_determinism_across_processes():
    """A fresh interpreter must produce bit-identical results.

    Catches global/process-level state (RNG, clocks, platform coupling) that a
    same-process double call cannot detect. Cheap: one subprocess spawn.
    """
    import json
    import subprocess
    import sys
    from pathlib import Path

    exp_dir = Path(__file__).resolve().parents[1]
    script = (
        "import json\n"
        "import sys\n"
        f"sys.path.insert(0, {str(exp_dir)!r})\n"
        "import experiment\n"
        "print(json.dumps(experiment.convergence_table(), sort_keys=True))\n"
    )
    out = subprocess.check_output([sys.executable, "-c", script], text=True).strip()
    expected = json.dumps(convergence_table(), sort_keys=True)
    assert out == expected, "results differ between interpreter processes"