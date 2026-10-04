"""Streaming RK4 propagator emitting the instantaneous nodal-rate series.

One propagator, one force model, used by every arm of the campaign. The
instantaneous rate is recorded at every integration step (not at node
crossings), so there is no node detection, no unwrapping and no window.

Force modes follow the predecessor `mission_mean_element_leakage` EXACTLY so
that a lunisolar rate measured here is comparable to the chain's committed
numbers:

    kepler_only  - two-body only
    j2_only      - two-body + J2
    sun_moon_j2  - two-body + J2 + point-mass Sun + point-mass Moon

The lunisolar term is isolated as (sun_moon_j2 rate) - (j2_only rate), both
estimated with the SAME estimator on the SAME grid, which removes the J2
common-mode exactly. That subtraction is why this mission needs only two modes
per inclination rather than the five the coupling mission needed.

Determinism: fixed step, fixed initial conditions, no RNG. The only
convergence knob is ``--dt``, used for the dt-convergence gate.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from lab_utils import J2_EARTH, MU_EARTH_KM3S2, R_EARTH_KM
from lab_utils.ephemeris import (LUNAR_GM, SOLAR_GM, interp_snapshot,
                                 load_snapshot, precession_j2000_to_mod)
from lab_utils.integrators import rk4_step

from estimator import instantaneous_node_rate, node_conditioning, node_rate_to_deg_day

RE2 = R_EARTH_KM ** 2
MODES = ("kepler_only", "j2_only", "sun_moon_j2")

#: T0 = 2026-01-01 00:00:00 TDB, the epoch every mission in the lunisolar chain
#: uses. Pinned so arcs are comparable across missions.
T0_S = 820476800.0


def snapshot_paths() -> tuple[Path, Path]:
    ref = Path(__file__).resolve().parent / "reference"
    return (ref / "horizons_sun_geocentric_vectors_2026_to_2064_icrf_tdb_daily.txt",
            ref / "horizons_moon_geocentric_vectors_2026_to_2064_icrf_tdb_daily.txt")


def load_snapshots() -> tuple[dict, dict]:
    sun_p, moon_p = snapshot_paths()
    return load_snapshot(sun_p), load_snapshot(moon_p)


def check_coverage(t0_s: float, t_end_s: float, sun: dict, moon: dict,
                   start_slack_days: float = 1.0) -> dict:
    """Fail loudly if the arc runs past the pinned ephemeris.

    ``interp_snapshot`` CLAMPS outside its range rather than raising, so a
    too-long arc silently freezes the Sun and Moon and yields a plausible,
    completely invalid rate. This is the guard against exactly that; the
    campaign calls it before propagating anything and the tests assert it.

    ``start_slack_days`` exists because the lab's ``T0_S = 820476800.0`` lies
    0.24 d BEFORE the first pinned snapshot row (JD 2461041.5). That offset is
    pre-existing in every mission in the chain and is harmless (it falls inside
    the first interpolation segment, so the Sun/Moon are extrapolated by at most
    5.8 h across a 1-day table), but it must not let the guard report a healthy
    arc as uncovered. The slack is reported explicitly rather than assumed.
    """
    t_lo = float(sun["t_s"][0])
    t_hi = float(sun["t_s"][-1])
    assert float(moon["t_s"][-1]) == t_hi, "Sun/Moon snapshot ends disagree"
    assert float(moon["t_s"][0]) == t_lo, "Sun/Moon snapshot starts disagree"
    end_ok = t_end_s <= t_hi
    start_ok = t0_s >= t_lo - start_slack_days * 86400.0
    return {
        "snapshot_start_s": t_lo, "snapshot_end_s": t_hi,
        "arc_start_s": t0_s, "arc_end_s": t_end_s,
        "covered": bool(start_ok and end_ok),
        "start_slack_days_used": (0.0 if t0_s >= t_lo else start_slack_days),
        "margin_days": (t_hi - t_end_s) / 86400.0,
        "arc_days": (t_end_s - t0_s) / 86400.0,
    }


def circular_ic(a_km: float, inc_rad: float, t0_s: float, phase: float = 0.0) -> np.ndarray:
    """Circular two-body initial state, RAAN = 0, phase = argument of latitude."""
    n = math.sqrt(MU_EARTH_KM3S2 / a_km ** 3)
    r_pf = np.array([a_km * math.cos(phase), a_km * math.sin(phase), 0.0])
    v_pf = np.array([-a_km * n * math.sin(phase), a_km * n * math.cos(phase), 0.0])
    R = np.array([[1.0, 0.0, 0.0],
                  [0.0, math.cos(inc_rad), -math.sin(inc_rad)],
                  [0.0, math.sin(inc_rad), math.cos(inc_rad)]])
    return np.concatenate([R @ r_pf, R @ v_pf])


def accel(r: np.ndarray, t: float, sun: dict, moon: dict, mode: str) -> np.ndarray:
    """Two-body + (mode-dependent) J2 + (mode-dependent) point-mass Sun/Moon.

    Written with scalar math on purpose. Profiling (2026-10-04) showed the
    original numpy-idiomatic version spent ~47 % of wall-clock in
    ``interp_snapshot`` + ``precession_j2000_to_mod`` and a further large share in
    ``np.linalg.norm`` on 3-element arrays. ``np.linalg.norm`` alone accounted for
    360 005 calls and 4.7 s of an 18.6 s 5-day benchmark.

    The pre-closed-form interpolation in :func:`_body_position` is
    ALGEBRAICALLY IDENTICAL to ``interp_snapshot`` and is pinned bit-for-bit
    against it by ``tests/test_propagator.py``. What it deliberately does NOT do
    is cache the precession rotation across steps: the 2026-10-02 stack audit
    found that "precess the snapshot once at load" is a 0.05 deg frame error
    because the rotation depends on ``t``, and that proposal was rejected. Here
    the rotation is applied at every evaluation, exactly as before.
    """
    rm2 = float(r[0]) ** 2 + float(r[1]) ** 2 + float(r[2]) ** 2
    rm = math.sqrt(rm2)
    mu = MU_EARTH_KM3S2
    a = (-mu * r[0] / (rm2 * rm), -mu * r[1] / (rm2 * rm), -mu * r[2] / (rm2 * rm))
    if mode in ("j2_only", "sun_moon_j2"):
        z2 = (r[2] * r[2]) / rm2
        c = -1.5 * J2_EARTH * mu * RE2 / (rm2 * rm2 * rm)
        g = 1.0 - 5.0 * z2
        a = (a[0] + c * r[0] * g, a[1] + c * r[1] * g, a[2] + c * r[2] * (3.0 - 5.0 * z2))
    if mode == "sun_moon_j2":
        ax = ay = az = 0.0
        for snap, gm in ((sun, SOLAR_GM), (moon, LUNAR_GM)):
            rx, ry, rz = _body_position(t, snap)
            dx, dy, dz = rx - r[0], ry - r[1], rz - r[2]
            dd2 = dx * dx + dy * dy + dz * dz
            d = dd2 * math.sqrt(dd2)
            b2 = rx * rx + ry * ry + rz * rz
            b = b2 * math.sqrt(b2)
            ax += gm * (dx / d - rx / b)
            ay += gm * (dy / d - ry / b)
            az += gm * (dz / d - rz / b)
        a = (a[0] + ax, a[1] + ay, a[2] + az)
    return np.array([a[0], a[1], a[2]])


def _body_position(t: float, snap: dict) -> tuple:
    """Precessed third-body position, algebraically identical to `interp_snapshot`.

    Linear interpolation between the two bracketing daily rows (clamped outside
    the table, preserving `interp_snapshot`'s documented clamp behaviour so the
    coverage guard remains the thing that catches an over-long arc), followed by
    the IAU-1976 precession rotation evaluated at this instant `t`.
    """
    ts = snap["t_s"]
    rows = snap["r_eci_km"]
    n = ts.size
    if t <= ts[0]:
        rx, ry, rz = rows[0]
    elif t >= ts[-1]:
        rx, ry, rz = rows[-1]
    else:
        i = int(np.searchsorted(ts, t))
        t_lo = ts[i - 1]
        frac = (t - t_lo) / (ts[i] - t_lo)
        lo = rows[i - 1]
        hi = rows[i]
        rx = lo[0] + frac * (hi[0] - lo[0])
        ry = lo[1] + frac * (hi[1] - lo[1])
        rz = lo[2] + frac * (hi[2] - lo[2])
    P = precession_j2000_to_mod(t)
    return (P[0, 0] * rx + P[0, 1] * ry + P[0, 2] * rz,
            P[1, 0] * rx + P[1, 1] * ry + P[1, 2] * rz,
            P[2, 0] * rx + P[2, 1] * ry + P[2, 2] * rz)


def propagate_rate_series(sun: dict, moon: dict, x0: np.ndarray, mode: str,
                          t0_s: float, t_end_s: float, dt_s: float,
                          every: int = 1) -> dict:
    """Propagate and record the instantaneous nodal rate every ``every`` steps.

    RK4 is applied to the full 6-vector (f returns [v, a(x,t)]), exactly as
    ``lab_utils.integrators.rk4_step`` expects. A hand-rolled variant that
    advances position with stage velocities but velocity with raw accelerations
    is WRONG and drove the semi-major axis negative within one sample.
    """
    x, t = np.asarray(x0, float).copy(), t0_s
    ts, rates = [], []
    n_steps = int(round((t_end_s - t0_s) / dt_s))
    for k in range(n_steps + 1):
        a = accel(x[:3], t, sun, moon, mode)
        if k % every == 0:
            ts.append(t)
            rates.append(instantaneous_node_rate(x, a))
        if k == n_steps:
            break
        x = rk4_step(lambda tt, xx: np.concatenate([xx[3:], accel(xx[:3], tt, sun, moon, mode)]),
                     t, x, dt_s)
        t += dt_s
    return {"t_s": np.array(ts),
            "rate_deg_day": np.asarray(node_rate_to_deg_day(np.array(rates)))}


def arc_mean_j2_rate(a_km: float, inc_rad: float) -> float:
    """Analytic first-order J2 secular nodal rate, deg/day (used as a gate)."""
    n = math.sqrt(MU_EARTH_KM3S2 / a_km ** 3) * 86400.0 * 180.0 / math.pi
    return -1.5 * n * J2_EARTH * (R_EARTH_KM / a_km) ** 2 * math.cos(inc_rad)