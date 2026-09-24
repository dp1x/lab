"""Independent averaged-element propagator (materially different algorithm).

Integrates ONLY the node angle (1D ODE), not Cowell Cartesian state:
  dOmega/dt = Omega_dot_J2(a,i,lambda_J2) + S0(a,i,lambda_3b) + sum_k A_k*cos(theta_k(t)),
  theta_k = m_k*(Omega(t)-Omega3(t)) + j_k*u3(t) + phi_k,
with Omega3(t)=Omega30+Omega3_dot*t, u3(t)=u30+n3*t (Sun/Moon circular means).
This is independent of mission_experiment.py Cowell RK4 (Cartesian + full
third-body accel + IAU precession + node detection). Agreement on sign/order
is independent evidence per constitution §3.2; shared ancestry (same S0/A_k
formulas from theory_modulation) is documented, independence comes from
propagation algorithm + estimator implementation.
"""
from __future__ import annotations

import math

import numpy as np


def propagate_averaged(a_km: float, inc_rad: float, *, W_day: float, dt_day: float = 1.0,
                       Omega0: float = 0.0, oj2_deg_day: float, s0_deg_day: float,
                       harmonics: list, n3_deg_day: float = 0.0, o3dot_deg_day: float = 0.0,
                       u30: float = 0.0, o30: float = 0.0) -> dict:
    """Euler/RK4 on 1D Omega ODE with prescribed harmonics.

    harmonics: list of (m, j, A_deg_day, phi). Returns t_day, Omega_unwrapped_deg.
    Deterministic.
    """
    n = int(round(W_day / dt_day)) + 1
    t = np.linspace(0.0, W_day, n)
    om = np.zeros(n)
    om[0] = Omega0
    o3dot = math.radians(o3dot_deg_day)
    n3 = math.radians(n3_deg_day)
    oj2r = math.radians(oj2_deg_day)
    # integrate with RK4 on dOm/dt = oj2 + s0 + sum A cos(m(Om-O3)+j u3+phi)
    for k in range(n - 1):
        dt = t[k + 1] - t[k]

        def rhs(om_deg: float, tk: float) -> float:
            o3 = o30 + o3dot * tk * 180.0 / math.pi
            u3 = u30 + n3 * tk * 180.0 / math.pi
            tot = oj2_deg_day + s0_deg_day
            for (m, j, A, phi) in harmonics:
                tot += A * math.cos(m * math.radians(om_deg - o3) + j * math.radians(u3) + phi)
            return tot

        k1 = rhs(om[k], t[k])
        k2 = rhs(om[k] + 0.5 * dt * k1, t[k] + 0.5 * dt)
        k3 = rhs(om[k] + 0.5 * dt * k2, t[k] + 0.5 * dt)
        k4 = rhs(om[k] + dt * k3, t[k] + dt)
        om[k + 1] = om[k] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return {"t_day": t, "Omega_deg": om}


def ols_slope_deg_day(t_day: np.ndarray, y_deg: np.ndarray) -> float:
    A = np.column_stack([np.ones_like(t_day), t_day])
    coef, *_ = np.linalg.lstsq(A, y_deg, rcond=None)
    return float(coef[1])
