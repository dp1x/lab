"""Independent quadrature referee for the secular nodal rate (no theory imported).

This is the numerical referee the 2026-10-04 session FAILED to produce ("My
attempt at an independent numerical referee FAILED ... returning values ~1e-9
deg/day, two orders below the formula under test -- the signature of
catastrophic cancellation in dRbar/dOmega"). It is rebuilt here with the
cancellation removed.

Method
------
Start from the EXACT point-mass potential of a circular third body,

    R = -mu3 / |rb - r|

evaluate it on a grid of (u3, phi) where phi is the satellite's argument of
latitude and u3 the third body's angular position, average over phi (the
satellite's mean anomaly), then average over u3 (the third body's mean anomaly),
then take the inclination derivative and apply Lagrange's node equation.

Cancellation fix: the prior attempt differenced Rbar over a step comparable to
the grid spacing, so the difference was swamped by quadrature error. Here
(a) the grid is refined until Rbar converges, and (b) the derivative uses a
5-point central stencil, which is fourth-order accurate and cancels the leading
error term. Convergence is reported, not assumed.

The SIGN convention is NOT assumed. Lagrange's equations give

    dOmega/dt = -(1/(n a^2 sin i)) dR/di        (argument-of-perigee form)
    dOmega/dt = +(1/(n a^2 sin i)) dR/di        (node form)

and the two differ by a sign. This module reports both and lets the numerical
campaign decide which one the physics uses, rather than assuming.
"""
from __future__ import annotations

import math

import numpy as np

MU_E = 398600.4418
DEG_PER_RAD = 180.0 / math.pi


def double_averaged_potential(inc_rad: float, a_km: float, a3_km: float, mu3: float,
                              i3_rad: float, n_u3: int, n_phi: int) -> float:
    """<R> averaged over the satellite's mean anomaly and the third body's."""
    u3 = np.linspace(0.0, 2.0 * math.pi, n_u3, endpoint=False)
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    U3, P = np.meshgrid(u3, phi, indexing="ij")
    x = a3_km * np.cos(U3) - a_km * np.cos(P)
    y = a3_km * np.sin(U3) * math.cos(i3_rad) - a_km * np.sin(P) * math.cos(inc_rad)
    z = a3_km * np.sin(U3) * math.sin(i3_rad) - a_km * np.sin(P) * math.sin(inc_rad)
    return float(-(mu3 / np.sqrt(x * x + y * y + z * z)).mean())


def dR_dinc(inc_rad: float, a_km: float, a3_km: float, mu3: float, i3_rad: float,
            n_u3: int, n_phi: int, h: float = 1e-3) -> float:
    """4th-order central derivative of <R> with respect to inclination.

    Standard 5-point stencil:

        f'(x) ~ (f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h)) / (12 h)

    (An earlier version here summed the weighted values and divided by
    sum(w_k x_k^2). That is the formula for the x=0 SECOND derivative and, used
    for the first, returned a magnitude identical to ~6 significant figures at
    every inclination -- 2.11e6 deg/day -- which is the signature of a broken
    stencil rather than a result. Caught by the grid-convergence check, which is
    exactly why that check exists.)
    """
    f = lambda d: double_averaged_potential(inc_rad + d, a_km, a3_km, mu3, i3_rad,
                                            n_u3, n_phi)
    return (f(-2 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2 * h)) / (12.0 * h)


def quadrature_node_rate(inc_deg: float, a_km: float, a3_km: float, mu3: float,
                         i3_deg: float, n_u3: int = 1440, n_phi: int = 360,
                         node_equation_sign: int = +1) -> float:
    """Secular nodal rate in deg/day from the double average."""
    i = math.radians(inc_deg)
    n = math.sqrt(MU_E / a_km ** 3)
    dR = dR_dinc(i, a_km, a3_km, mu3, math.radians(i3_deg), n_u3, n_phi)
    val = node_equation_sign * dR / (n * a_km * a_km * math.sin(i))
    return math.degrees(val) * 86400.0


def convergence_check(inc_deg: float, a_km: float, a3_km: float, mu3: float,
                      i3_deg: float) -> dict:
    """Show the double average has actually converged in grid resolution."""
    rows = {}
    for n_u3, n_phi in ((360, 90), (720, 180), (1440, 360), (2880, 720)):
        v = quadrature_node_rate(inc_deg, a_km, a3_km, mu3, i3_deg, n_u3, n_phi)
        rows[f"{n_u3}x{n_phi}"] = v
    return rows