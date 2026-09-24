"""Closed-form J2-precession-modulated lunisolar theory (no fitted constants).

Implements DERIVATION.md: doubly-averaged secular S0 (Convention B),
singly-averaged long-period amplitudes A_k via quadrupole quadrature
(derived from model geometry, NOT fitted to RK4 R), beat frequencies
omega_k containing Omega_dot_J2, and exact finite-window OLS bias g().

All rates in deg/day unless noted. No empirical coefficient.
"""
from __future__ import annotations

import math

import numpy as np

from lab_utils import J2_EARTH, MU_EARTH_KM3S2, R_EARTH_KM

SOLAR_GM = 132712440018.0
LUNAR_GM = 4902.8001
AU_KM = 149597870.7
LUNAR_A_KM = 384400.0
OBLIQUITY_DEG = 23.4392911
LUNAR_I3_MEAN_DEG = 28.584  # secular mean to equator (2026 actual ~18.3; bounded as variation)
DEG_PER_DAY = 180.0 / math.pi  # rad/day -> deg/day via *86400? handled explicitly
N_SUN_DEG_DAY = 0.98560912
OMEGA3_DOT_DEG_DAY = -360.0 / 6798.4  # lunar nodal regression
N_MOON_DEG_DAY = 360.0 / 27.321661


def mean_motion_rad_s(a_km: float) -> float:
    return math.sqrt(MU_EARTH_KM3S2 / a_km ** 3)


def omega_dot_j2_deg_day(a_km: float, inc_rad: float) -> float:
    n_rad_s = mean_motion_rad_s(a_km)
    n_deg_day = n_rad_s * 86400.0 * 180.0 / math.pi
    return -1.5 * n_deg_day * J2_EARTH * (R_EARTH_KM / a_km) ** 2 * math.cos(inc_rad)


def secular_S0_deg_day(a_km: float, inc_rad: float, mu3: float, a3_km: float, i3_rad: float) -> float:
    n_rad_s = mean_motion_rad_s(a_km)
    n_deg_day = n_rad_s * 86400.0 * 180.0 / math.pi
    return (3.0 / 8.0) * n_deg_day * (mu3 / MU_EARTH_KM3S2) * (a_km / a3_km) ** 3 * math.sin(2.0 * (inc_rad - i3_rad)) / math.sin(inc_rad)


def secular_total_deg_day(a_km: float, inc_rad: float) -> float:
    s_sun = secular_S0_deg_day(a_km, inc_rad, SOLAR_GM, AU_KM, math.radians(OBLIQUITY_DEG))
    s_moon = secular_S0_deg_day(a_km, inc_rad, LUNAR_GM, LUNAR_A_KM, math.radians(LUNAR_I3_MEAN_DEG))
    return s_sun + s_moon


# --------------------------------------------------------------------------- #
# Quadrupole quadrature for long-period amplitudes (model-derived, no fit)
# --------------------------------------------------------------------------- #
def _r_hat(omega: float, inc: float, u: float) -> np.ndarray:
    cO, sO = math.cos(omega), math.sin(omega)
    ci, si = math.cos(inc), math.sin(inc)
    cu, su = math.cos(u), math.sin(u)
    # R_z(Omega) R_x(i) R_z(u) [1,0,0]
    x1, y1, z1 = cu, su, 0.0
    x2, y2, z2 = x1, ci * y1, si * y1
    return np.array([cO * x2 - sO * y2, sO * x2 + cO * y2, z2])


def quadrupole_R(mu3: float, a_km: float, a3_km: float, omega: float, inc: float, u: float,
                 omega3: float, i3: float, u3: float) -> float:
    rhat = _r_hat(omega, inc, u)
    r3hat = _r_hat(omega3, i3, u3)
    c = float(np.dot(rhat, r3hat))
    p2 = 0.5 * (3.0 * c * c - 1.0)
    return (mu3 / a3_km) * (a_km / a3_km) ** 2 * p2


def singly_averaged_rate_grid(a_km: float, inc: float, mu3: float, a3_km: float, i3: float,
                              n_omega: int = 72, n_u3: int = 72) -> tuple:
    """Return (S0_quad, harmonics) where harmonics is list of (m, j, A_deg_day, phi).

    Computes R_bar(Omega; u3) by quadrature over satellite u, differentiates
    w.r.t i by central difference, converts to node rate, then Fourier
    analyzes in (Omega-Omega3) and u3. Model-derived only.
    """
    n_rad_s = mean_motion_rad_s(a_km)
    n_deg_day = n_rad_s * 86400.0 * 180.0 / math.pi
    denom = n_rad_s * a_km * a_km * math.sin(inc)  # SI: R in km^2/s^2 -> rate rad/s
    di = 1e-6
    omegas = np.linspace(0, 2 * math.pi, n_omega, endpoint=False)
    u3s = np.linspace(0, 2 * math.pi, n_u3, endpoint=False)
    us = np.linspace(0, 2 * math.pi, 144, endpoint=False)
    # R_bar[io, ju3]
    rbar = np.zeros((n_omega, n_u3))
    for io, om in enumerate(omegas):
        for ju, u3 in enumerate(u3s):
            tot_p = 0.0
            tot_m = 0.0
            for u in us:
                tot_p += quadrupole_R(mu3, a_km, a3_km, om, inc + di, u, 0.0, i3, u3)
                tot_m += quadrupole_R(mu3, a_km, a3_km, om, inc - di, u, 0.0, i3, u3)
            rbar[io, ju] = (tot_p - tot_m) / (2 * di * len(us))
    rate_rad_s = rbar / denom
    rate_deg_day = rate_rad_s * 86400.0 * 180.0 / math.pi
    s0 = float(np.mean(rate_deg_day))
    # Fourier in Omega (m=0..2) and u3 (j=-2..2): least squares for amplitudes
    harmonics = []
    t_om, t_u3 = np.meshgrid(omegas, u3s, indexing="ij")
    y = (rate_deg_day - s0).ravel()
    for m in (1, 2):
        for j in (-2, -1, 1, 2):
            c = np.cos(m * t_om + j * t_u3).ravel()
            s = np.sin(m * t_om + j * t_u3).ravel()
            Ac = float(np.dot(y, c) / np.dot(c, c))
            As = float(np.dot(y, s) / np.dot(s, s))
            amp = math.hypot(Ac, As)
            phi = math.atan2(-As, Ac)
            harmonics.append((m, j, amp, phi))
    # also m-only (u3-averaged long-period in Omega-Omega3)
    for m in (1, 2):
        c = np.cos(m * t_om).ravel()
        s = np.sin(m * t_om).ravel()
        Ac = float(np.dot(y, c) / np.dot(c, c))
        As = float(np.dot(y, s) / np.dot(s, s))
        amp = math.hypot(Ac, As)
        phi = math.atan2(-As, Ac)
        harmonics.append((m, 0, amp, phi))
    return s0, harmonics


def leading_amplitudes_deg_day(a_km: float, inc_rad: float) -> dict:
    """Sum Sun+Moon leading long-period rate amplitudes (model-derived)."""
    _, hs = singly_averaged_rate_grid(a_km, inc_rad, SOLAR_GM, AU_KM, math.radians(OBLIQUITY_DEG),
                                      n_omega=48, n_u3=48)
    _, hm = singly_averaged_rate_grid(a_km, inc_rad, LUNAR_GM, LUNAR_A_KM, math.radians(LUNAR_I3_MEAN_DEG),
                                      n_omega=48, n_u3=48)
    out = {"sun": hs, "moon": hm}
    return out


# --------------------------------------------------------------------------- #
# Exact OLS bias of A*sin(omega t + phi)/omega angle term
# --------------------------------------------------------------------------- #
def ols_bias_deg_day(A_deg_day: float, omega_rad_day: float, W_day: float, phi: float) -> float:
    """Exact dense-sampling OLS slope bias of rate harmonic A*cos(omega t+phi).

    Angle: s*t + (A/omega)*sin(omega t+phi). Returns bias (deg/day).
    Handles omega->0 limit (resonant -> A*cos phi).
    """
    if abs(omega_rad_day) < 1e-14:
        return A_deg_day * math.cos(phi)
    x = omega_rad_day * W_day
    # Exact projection: slope of f(t)=(A/omega) sin(omega t+phi) on [0,W]
    # s_hat = Cov(t,f)/Var(t); Var= W^2/12; Cov = (1/W)int_0^W (t-W/2) f(t) dt.
    # Closed form:
    # int (t-W/2) sin(omega t+phi) dt = ...
    om = omega_rad_day
    W = W_day
    # indefinite: int t sin(om t+ph) = -t cos(om t+ph)/om + sin(om t+ph)/om^2
    def J(a: float, b: float) -> float:
        Fa = -b * math.cos(om * b + phi) / om + math.sin(om * b + phi) / (om * om)
        F0 = -a * math.cos(om * a + phi) / om + math.sin(om * a + phi) / (om * om)
        return Fa - F0
    # Cov = (A/om) * (1/W) * [ J(0,W) - (W/2)*( -cos(om W+ph)/om + cos(ph)/om ) ]
    int_sin = (-math.cos(om * W + phi) + math.cos(phi)) / om
    cov = (A_deg_day / om) * (J(0.0, W) - 0.5 * W * int_sin) / W
    var = W * W / 12.0
    return cov / var


def predicted_R_deg_day(a_km: float, inc_rad: float, W_day: float, phi_sun: float = 0.0,
                        phi_moon: float = 0.0, omega_dot_j2_override: float | None = None) -> dict:
    """Predict finite-window extra slope B = sum B_k for H-mod (no fit).

    Beats: solar m=1 with Omega_dot_J2 - n_sun; lunar m=1 with
    Omega_dot_J2 - Omega3_dot; plus monthly sidebands. Amplitudes from
    quadrature (leading m=1, j=0 terms). Phases from epoch (passed in).
    Returns dict with S0, per-term B, total B (= predicted R up to H-alias).
    """
    oj2 = omega_dot_j2_deg_day(a_km, inc_rad) if omega_dot_j2_override is None else omega_dot_j2_override
    s0 = secular_total_deg_day(a_km, inc_rad)
    amps = leading_amplitudes_deg_day(a_km, inc_rad)
    # leading m=1,j=0 amplitudes
    def lead(h, m=1, j=0):
        for (mm, jj, A, _) in h:
            if mm == m and jj == j:
                return A
        return 0.0
    A_sun = lead(amps["sun"], 1, 0)
    A_moon = lead(amps["moon"], 1, 0)
    om_sun = math.radians(oj2 - N_SUN_DEG_DAY)  # rad/day (oj2,n_sun in deg/day)
    om_moon = math.radians(oj2 - OMEGA3_DOT_DEG_DAY)
    B_sun = ols_bias_deg_day(A_sun, om_sun, W_day, phi_sun)
    B_moon = ols_bias_deg_day(A_moon, om_moon, W_day, phi_moon)
    return {"S0": s0, "A_sun": A_sun, "A_moon": A_moon, "om_sun_deg_day": oj2 - N_SUN_DEG_DAY,
            "om_moon_deg_day": oj2 - OMEGA3_DOT_DEG_DAY, "B_sun": B_sun, "B_moon": B_moon,
            "B_total": B_sun + B_moon, "Omega_dot_J2": oj2}


def octupole_bound_deg_day(a_km: float, inc_rad: float) -> float:
    # |(a/a3)| suppression vs quadrupole, summed Sun+Moon, O(1) geometry bound
    n_deg_day = mean_motion_rad_s(a_km) * 86400.0 * 180.0 / math.pi
    q_sun = (n_deg_day * (SOLAR_GM / MU_EARTH_KM3S2) * (a_km / AU_KM) ** 4) / abs(math.sin(inc_rad))
    q_moon = (n_deg_day * (LUNAR_GM / MU_EARTH_KM3S2) * (a_km / LUNAR_A_KM) ** 4) / abs(math.sin(inc_rad))
    return q_sun + q_moon
