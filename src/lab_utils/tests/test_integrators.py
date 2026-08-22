"""Tests for lab_utils.integrators (generic dimension-agnostic RK4).

Doctrine mirrors the experiment suites: theory values duplicated inline,
order-of-accuracy proven against analytic truth, and equivalence pinned
against the donor per-experiment RK4 loops this module retires clones of.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np

from lab_utils.integrators import rk4_propagate, rk4_step
from lab_utils.metrics import convergence_rate

_EXPERIMENTS_DIR = Path(__file__).resolve().parents[3] / "research" / "orbital-mechanics" / "experiments"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# L1 -- analytic truth: harmonic oscillator x'' = -x
# --------------------------------------------------------------------------- #
def test_rk4_step_sho_single_step_accuracy():
    """One step from x=(1,0): error must match the classic 1 - (h^4/24)*cos(h) bound scale."""
    h = 0.1

    def f(t, x):
        return np.array([x[1], -x[0]])

    x1 = rk4_step(f, 0.0, np.array([1.0, 0.0]), h)
    exact = np.array([math.cos(h), -math.sin(h)])
    # RK4 local error on SHO is O(h^5): h=0.1 -> ~1e-9 scale (measured 1.39e-9)
    assert abs(x1[0] - exact[0]) < 1e-8


def test_rk4_propagate_sho_order_four():
    """Global error must converge at order 4 (log-log slopes ~ 4)."""

    def f(t, x):
        return np.array([x[1], -x[0]])

    T = 2.0 * math.pi
    ns = np.array([50, 100, 200, 400], dtype=float)
    errs = []
    for n in ns:
        t = np.linspace(0.0, T, int(n) + 1)
        traj = rk4_propagate(f, t, np.array([1.0, 0.0]))
        errs.append(abs(traj[-1, 0] - 1.0))
    # convergence_rate takes matching (errors, stepsizes); SHO shows clean order 5
    orders = convergence_rate(np.array(errs), T / ns)
    assert np.all(orders > 4.5), f"RK4 oscillator convergence degraded: {orders}"


def test_rk4_propagate_grid_and_shapes():
    """Row 0 is x0; length matches grid; non-uniform grids accepted."""

    def f(t, x):
        return np.array([x[1], -x[0]])

    t = np.array([0.0, 0.01, 0.03, 0.1, 0.25])  # deliberately non-uniform
    traj = rk4_propagate(f, t, np.array([1.0, 2.0]))
    assert traj.shape == (5, 2)
    assert np.array_equal(traj[0], np.array([1.0, 2.0]))
    # closed-form check at the last point; final coarse step h=0.15 gives an
    # RK4 local error of order h^5/120 ~ 7.6e-7 (measured total 1.45e-6)
    exact = np.array([
        math.cos(0.25) + 2 * math.sin(0.25),
        -math.sin(0.25) + 2 * math.cos(0.25),
    ])
    assert np.linalg.norm(traj[-1] - exact) < 1e-5


def test_rk4_batched_state_shape():
    """Batched (n_states, dim) input stays supported through the step."""

    def f(t, x):
        return np.stack([x[..., 1], -x[..., 0]], axis=-1)

    xs = np.array([[1.0, 0.0], [0.0, 1.0]])
    out = rk4_step(f, 0.0, xs, 0.05)
    assert out.shape == (2, 2)
    exact = np.stack([[math.cos(0.05), -math.sin(0.05)], [math.sin(0.05), math.cos(0.05)]])
    # single-step local error at h=0.05 is O(h^5) ~ 2.6e-9 (measured)
    assert np.max(np.abs(out - exact)) < 1e-7


# --------------------------------------------------------------------------- #
# L2 -- equivalence pin vs donor experiment RK4 loops (Exp 009 J2=0 path)
# --------------------------------------------------------------------------- #
def _two_body_rhs(mu):
    def f(t, x):
        r = x[:3]
        a = -mu * r / np.linalg.norm(r) ** 3
        return np.concatenate([x[3:], a])

    return f


def test_equivalence_with_donor_kepler_propagation():
    """Generic rk4_propagate vs Exp 002's dedicated Kepler RK4 on the same grid."""
    exp002 = _load("exp002_for_integrators", _EXPERIMENTS_DIR / "keplerOrbitValidation" / "experiment.py")
    mu = 398600.4418  # SI-ish km^3/s^2 value passed explicitly to both sides
    a, e = 8000.0, 0.3
    r0, v0 = exp002.initial_state(a, e, mu)
    T = 2 * math.pi * math.sqrt(a**3 / mu)
    n = 720
    t = np.linspace(0.0, T, n + 1)

    def f(tt, x):
        return np.concatenate([x[2:], -mu * x[:2] / np.linalg.norm(x[:2]) ** 3])

    mine = rk4_propagate(f, t, np.concatenate([np.asarray(r0), np.asarray(v0)]))
    theirs = np.asarray(exp002.propagate_rk4(np.asarray(r0), np.asarray(v0), mu, t))
    assert theirs.shape == mine.shape
    denom = np.max(np.abs(mine))
    rel = np.max(np.abs(mine - theirs)) / denom
    assert rel < 1e-11, f"donor-vs-generic RK4 mismatch rel={rel:.3e}"


def test_equivalence_bit_level_donor_6d():
    """Same arithmetic ordering as donor 6-D loops => expect near-bit-exact agreement."""
    exp009 = _load("exp009_for_integrators", _EXPERIMENTS_DIR / "j2Precession" / "experiment.py")
    a, e, inc = 7000.0, 0.01, 0.5
    r0, v0, _ = exp009.seed_state(a, e, inc, 0.0, 0.0, 0.0, exp009.MU_EARTH_KM3S2)
    T = exp009.orbital_period(a)
    t = np.linspace(0.0, T, 257)
    ref = exp009.propagate_3d_rk4_j2(np.asarray(r0), np.asarray(v0), exp009.MU_EARTH_KM3S2, t, 0.0)  # J2=0
    mine = rk4_propagate(_two_body_rhs(exp009.MU_EARTH_KM3S2), t, np.concatenate([np.asarray(r0), np.asarray(v0)]))
    scale = np.max(np.abs(mine))
    rel = np.max(np.abs(mine - ref)) / scale
    assert rel < 1e-13, f"6D donor-vs-generic mismatch rel={rel:.3e}"
