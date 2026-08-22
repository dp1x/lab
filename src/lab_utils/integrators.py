"""Generic fixed-step Runge-Kutta 4 integrators (dimension-agnostic).

Canonical copies retire the per-experiment RK4 loop clones (Exp 002/006/008/009/010
each carry a frozen local version; this module is the shared form going forward).
The RHS signature ``f(t, x)`` supports non-autonomous systems, which is required for
inertial-frame validators whose primary positions advance through the RK4 stages --
freezing stage-time-dependent terms silently reduces the scheme to first order
(Exp 011 doctrine).

Covered by tests in ``src/lab_utils/tests/test_integrators.py``.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

__all__ = ["rk4_step", "rk4_propagate"]


def rk4_step(f: Callable[[float, np.ndarray], np.ndarray], t: float, x: np.ndarray, h: float) -> np.ndarray:
    """One classical RK4 step for dx/dt = f(t, x). Returns the new state.

    x may have any shape supported by f and the arithmetic below (flat state
    vector or (n, m) batch); no copying beyond one new output array.
    """
    k1 = f(t, x)
    k2 = f(t + 0.5 * h, x + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, x + 0.5 * h * k2)
    k4 = f(t + h, x + h * k3)
    return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def rk4_propagate(
    f: Callable[[float, np.ndarray], np.ndarray],
    t: np.ndarray,
    x0: np.ndarray,
) -> np.ndarray:
    """Fixed-step RK4 trajectory over a time grid. Returns (len(t), len(x0)).

    Step sizes come from consecutive grid points (non-uniform grids safe,
    matching the Exp 006/009/010 loop convention). t[0] is the epoch;
    row k is the state at t[k] (row 0 is x0 itself).
    """
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or t.shape[0] < 1:
        raise ValueError("t must be a 1-D array with at least one point")
    x0 = np.asarray(x0, dtype=float)
    traj = np.empty((t.shape[0], x0.shape[0]), dtype=float)
    traj[0] = x0
    x = x0.copy()
    for k in range(1, t.shape[0]):
        x = rk4_step(f, float(t[k - 1]), x, float(t[k] - t[k - 1]))
        traj[k] = x
    return traj
