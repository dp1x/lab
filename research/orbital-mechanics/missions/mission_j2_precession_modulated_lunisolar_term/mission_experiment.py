"""Streaming RK4 + 6-mode isolation + lambda scaling + prescribed-precession controls.

Reuses force-model constants from mission_j2_lunisolar_coupling / closure
missions. snapshot interpolation is CORRECTED here (see REMEDIATION NOTE):
prior missions used rv = r_lo + frac*(r_lo - r_hi) (sign-reversed weights);
this mission uses rv = r_lo + frac*(r_hi - r_lo), pinned by regression test.

Modes: kepler_only, j2_only, sun_only, moon_only, sun_moon, sun_moon_j2.
Corrected isolation: use_j2 = mode in ("j2_only","sun_moon_j2").
Controls (ACTIVE rigid-rotation semantics):
  omega_p_deg_day: rigid z-rotation of the full Cartesian state (r AND v)
    applied after each inertial RK4 step (kinematic plane precession
    without J2 force when use_j2=False).
  freeze_plane: open-loop counter-rotation at the analytic J2 secular rate
    for J2 modes (common rotation; alters Sun/Moon geometry, does not hold
    the physical plane fixed -- gate (f) is evaluated as written).
Deterministic: fixed dt, fixed snapshots, fixed ICs.
NOTE 2026-09-14: a passive rotating-frame formulation (Coriolis/centrifugal
  + measure_in_rotating_frame, tried 2026-09-12) never went green
  (observable sign defect) and is reverted. All campaigns use the active
  semantics documented here.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

from lab_utils import J2_EARTH, MU_EARTH_KM3S2, R_EARTH_KM
from lab_utils.earth_frames import JD_J2000
from lab_utils.integrators import rk4_step

HERE = Path(__file__).resolve().parent
SUN_SNAPSHOT = HERE.parent / "mission_lunisolar_closure" / "reference" / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MOON_SNAPSHOT = HERE.parent / "mission_lunisolar_closure" / "reference" / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"

SOLAR_GM = 132712440018.0
LUNAR_GM = 4902.8001
DT_S = 60.0
DEG = math.pi / 180.0


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _rot3(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot2(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])


def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    T = t_s / (86400.0 * 36525.0)
    sec = math.radians(1.0 / 3600.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * sec
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * sec
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * sec
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)


def load_snapshot(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    soe = eoe = None
    for i, ln in enumerate(lines):
        if "$$SOE" in ln:
            soe = i + 1
        if "$$EOE" in ln:
            eoe = i
            break
    rows = []
    for ln in lines[soe:eoe]:
        s = ln.strip()
        if not s:
            continue
        p = [x.strip() for x in s.split(",")]
        rows.append((float(p[0]), float(p[2]), float(p[3]), float(p[4])))
    arr = np.array(rows)
    return {"t_s": (arr[:, 0] - JD_J2000) * 86400.0, "r_eci_km": arr[:, 1:4],
            "sha256": _sha256(path), "n_points": len(rows)}


def interp_snapshot(t_q: float, snap: dict, apply_precession: bool = True) -> np.ndarray:
    """Correct linear interpolation (regression-pinned)."""
    ts = snap["t_s"]
    r = snap["r_eci_km"]
    if t_q <= ts[0]:
        rv = r[0]
    elif t_q >= ts[-1]:
        rv = r[-1]
    else:
        idx = int(np.searchsorted(ts, t_q))
        t_lo, t_hi = ts[idx - 1], ts[idx]
        frac = (t_q - t_lo) / (t_hi - t_lo)
        rv = r[idx - 1] + frac * (r[idx] - r[idx - 1])
    if apply_precession:
        return precession_j2000_to_mod(t_q) @ rv
    return rv


def interp_snapshot_buggy(t_q: float, snap: dict, apply_precession: bool = True) -> np.ndarray:
    """Legacy sign-reversed weights, preserved ONLY for bug-impact quantification."""
    ts = snap["t_s"]
    r = snap["r_eci_km"]
    if t_q <= ts[0]:
        rv = r[0]
    elif t_q >= ts[-1]:
        rv = r[-1]
    else:
        idx = int(np.searchsorted(ts, t_q))
        t_lo, t_hi = ts[idx - 1], ts[idx]
        frac = (t_q - t_lo) / (t_hi - t_lo)
        rv = r[idx - 1] + frac * (r[idx - 1] - r[idx])
    if apply_precession:
        return precession_j2000_to_mod(t_q) @ rv
    return rv


def third_body_accel(r_eci: np.ndarray, t_s: float, sun: dict, moon: dict, *,
                     lambda_3body: float = 1.0, include_sun: bool = True,
                     include_moon: bool = True, use_buggy_interp: bool = False) -> np.ndarray:
    interp = interp_snapshot_buggy if use_buggy_interp else interp_snapshot
    a = np.zeros(3)
    if include_sun:
        rs = interp(t_s, sun, True)
        d = rs - r_eci
        a += SOLAR_GM * (d / np.linalg.norm(d) ** 3 - rs / np.linalg.norm(rs) ** 3)
    if include_moon:
        rm = interp(t_s, moon, True)
        d = rm - r_eci
        a += LUNAR_GM * (d / np.linalg.norm(d) ** 3 - rm / np.linalg.norm(rm) ** 3)
    return lambda_3body * a


def circular_ic(a_km: float, inc_rad: float, t0_s: float, phase: float = 0.0) -> np.ndarray:
    n = math.sqrt(MU_EARTH_KM3S2 / a_km ** 3)
    M = phase
    r_pf = np.array([a_km * math.cos(M), a_km * math.sin(M), 0.0])
    v_pf = np.array([-a_km * n * math.sin(M), a_km * n * math.cos(M), 0.0])
    R = _rot3(0.0) @ np.array([[1, 0, 0], [0, math.cos(inc_rad), -math.sin(inc_rad)],
                               [0, math.sin(inc_rad), math.cos(inc_rad)]])
    # Omega0 = 0
    return np.concatenate([R @ r_pf, R @ v_pf])


def _rot_z_vec(v: np.ndarray, ang: float) -> np.ndarray:
    c, s = math.cos(ang), math.sin(ang)
    x, y, z = v
    return np.array([c * x - s * y, s * x + c * y, z])


def propagate(sun: dict, moon: dict, x0: np.ndarray, *, mode: str, t0_s: float, t_end_s: float,
              dt_s: float = DT_S, lambda_j2: float = 1.0, lambda_3body: float = 1.0,
              omega_p_deg_day: float = 0.0, freeze_plane: bool = False,
              freeze_plane_rate_deg_day: float | None = None,
              use_buggy_interp: bool = False, subsample_every: int = 100,
              **legacy_kwargs) -> dict:
    if legacy_kwargs:
        raise TypeError(f"removed control kwargs (passive-frame experiment): {sorted(legacy_kwargs)}")
    use_sun = mode in ("sun_only", "sun_moon", "sun_moon_j2")
    use_moon = mode in ("moon_only", "sun_moon", "sun_moon_j2")
    use_j2 = mode in ("j2_only", "sun_moon_j2")
    re2 = R_EARTH_KM ** 2

    # ACTIVE rigid-rotation control (reactivated 2026-09-14; see module docstring).
    # Post-step rigid rotation of the full Cartesian state by dtheta about z:
    #   - omega_p_deg_day: rate applied to EVERY mode (kinematic precession).
    #   - freeze_plane: rate = analytic J2 secular rate applied to J2 modes only.
    # A common rigid rotation shifts every mode's raw node rate by the same
    # amount, so R = (full-j2) - ((sun-kep)+(moon-kep)) is invariant under it;
    # the prescribed sweep probes how R depends on the *commanded* rate, and
    # the frozen run must be read as an open-loop counter-rotation control
    # (it keeps the Sun/Moon field live), not a physical plane freeze.
    dtheta_prescribed = math.radians(omega_p_deg_day) * dt_s / 86400.0
    freeze_rate_deg_day = 0.0
    dtheta_freeze = 0.0
    if freeze_plane and use_j2:
        if freeze_plane_rate_deg_day is None:
            a0 = float(np.linalg.norm(x0[:3]))
            h0n = np.cross(x0[:3], x0[3:])
            h0n = h0n / float(np.linalg.norm(h0n))
            ci0 = max(-1.0, min(1.0, float(h0n[2])))
            n0 = math.sqrt(MU_EARTH_KM3S2 / a0 ** 3) * 86400.0 * 180.0 / math.pi
            freeze_rate_deg_day = -1.5 * n0 * J2_EARTH * (R_EARTH_KM / a0) ** 2 * ci0 * lambda_j2
        else:
            freeze_rate_deg_day = freeze_plane_rate_deg_day
        # Counter-rotation: analytic J2 rate is +0.99 deg/day at SSO, so the
        # applied rigid rotation is its negation (measured = raw - freeze_rate).
        dtheta_freeze = -math.radians(freeze_rate_deg_day) * dt_s / 86400.0

    def f_inertial(t, x_i):
        r, v = x_i[:3], x_i[3:]
        rm = float(np.linalg.norm(r))
        a = -MU_EARTH_KM3S2 * r / rm ** 3
        if use_j2 and lambda_j2 != 0.0:
            z2 = (r[2] * r[2]) / (rm * rm)
            c = -1.5 * lambda_j2 * J2_EARTH * MU_EARTH_KM3S2 * re2 / rm ** 5
            g = 1.0 - 5.0 * z2
            a = a + c * np.array([r[0] * g, r[1] * g, r[2] * (3.0 - 5.0 * z2)])
        if (use_sun or use_moon) and lambda_3body != 0.0:
            a = a + third_body_accel(r, t, sun, moon, lambda_3body=lambda_3body,
                                     include_sun=use_sun, include_moon=use_moon,
                                     use_buggy_interp=use_buggy_interp)
        return np.concatenate([v, a])

    def apply_control(x_s):
        if dtheta_prescribed != 0.0:
            R = _rot3(dtheta_prescribed)
            x_s[:3] = R @ x_s[:3]
            x_s[3:] = R @ x_s[3:]
        if dtheta_freeze != 0.0:
            R = _rot3(dtheta_freeze)
            x_s[:3] = R @ x_s[:3]
            x_s[3:] = R @ x_s[3:]
        return x_s

    def node_of(x_s):
        h = np.cross(x_s[:3], x_s[3:])
        return math.atan2(-h[0], h[1])

    t_cross, om_cross, t_node, om_node = [], [], [t0_s], []
    om0 = node_of(x0)
    om_node.append(om0)
    om_prev = om0
    x, t = x0.copy(), t0_s
    z_prev = x[2]
    n_steps = 0
    while t < t_end_s:
        x_new = apply_control(rk4_step(f_inertial, t, x, dt_s))
        z_curr = x_new[2]
        if z_prev <= 0 < z_curr:
            frac = -z_prev / (z_curr - z_prev)
            rc = x[:3] + frac * (x_new[:3] - x[:3])
            t_cross.append(t + frac * dt_s)
            om_cross.append(math.atan2(rc[1], rc[0]))
        if n_steps % subsample_every == 0:
            o = node_of(x_new)
            while o < om_prev - math.pi:
                o += 2 * math.pi
            while o > om_prev + math.pi:
                o -= 2 * math.pi
            t_node.append(t + dt_s)
            om_node.append(o)
            om_prev = o
        z_prev, x, t = z_curr, x_new, t + dt_s
        n_steps += 1
    om_arr = np.array(om_cross)
    om_unw = np.unwrap(om_arr) if len(om_arr) > 1 else om_arr
    return {"t_cross": np.array(t_cross), "om_cross": om_unw,
            "t_node": np.array(t_node), "omega_node": np.array(om_node),
            "n_steps": n_steps, "mode": mode, "lambda_j2": lambda_j2,
            "lambda_3body": lambda_3body, "omega_p_deg_day": omega_p_deg_day,
            "freeze_plane": freeze_plane,
            "freeze_plane_rate_deg_day": freeze_rate_deg_day,
            "use_buggy_interp": use_buggy_interp}


def ols_slope(t_s: np.ndarray, y_rad: np.ndarray) -> float:
    """Slope in rad/s."""
    A = np.column_stack([np.ones_like(t_s), t_s])
    coef, *_ = np.linalg.lstsq(A, y_rad, rcond=None)
    return float(coef[1])


def rate_deg_day(t_s: np.ndarray, y_rad: np.ndarray) -> float:
    return ols_slope(t_s, y_rad) * 86400.0 * 180.0 / math.pi


def code_hashes() -> dict:
    import pathlib
    out = {}
    for p in sorted((HERE).glob("*.py")):
        out[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    return out
