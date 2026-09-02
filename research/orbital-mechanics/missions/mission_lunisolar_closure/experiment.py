"""Experiment runner for mission_lunisolar_closure.

Memory-efficient: streams RK4 step-by-step, detects ascending-node
crossings and node-vector samples on the fly. Does NOT store the full
trajectory (which would require ~470 MB float64 for an 18.6-yr arc at
dt=60s and overwhelms memory bandwidth on commodity hardware).

For the harmonic regression (Estimator f) and node-vector OLS
(Estimator n), we only need samples at ascending-node crossings plus
a sparse trajectory sample every N steps for the harmonic basis fit.
The phase-locked 2-window estimator uses only the ascending-node
crossings within each window.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from lab_utils import (
    J2_EARTH,
    MU_EARTH_KM3S2,
    R_EARTH_KM,
    j2_rhs,
)
from lab_utils.earth_frames import JD_J2000
from lab_utils.integrators import rk4_step
from lab_utils.results import save_json_result

EXP_NAME = "mission_lunisolar_closure_001"
SOLAR_GM_KM3_S2 = 132712440018.0
LUNAR_GM_KM3_S2 = 4902.8001
LUNAR_DISTANCE_KM_MEAN = 384400.0
LUNAR_INCLINATION_DEG = 5.145
SOLAR_OBLIQUITY_DEG = 23.439
AU_KM = 149597870.7
DEG = math.pi / 180.0

H_SSO_KM = 600.0
I_SSO_DEG = 97.7876
I_90_DEG = 90.0
I_30_DEG = 30.0

DT_S = 60.0

HORIZONS_DAYS = {
    "W_1yr": 365.25,
    "W_2yr": 2.0 * 365.25,
    "W_5yr": 5.0 * 365.25,
    "W_9p3yr": 9.3067 * 365.25,
    "W_18p6yr": 18.6 * 365.25,
}

HERE = Path(__file__).resolve().parent
SUN_SNAPSHOT = HERE / "reference" / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MOON_SNAPSHOT = HERE / "reference" / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MANIFEST_PATH = HERE / "reference" / "MANIFEST.json"

LUNAR_NODAL_PERIOD_DAYS = 6798.4
HALF_NODAL_DAYS = LUNAR_NODAL_PERIOD_DAYS / 2.0


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# IAU-1976 precession
# --------------------------------------------------------------------------- #
def _rot3(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot2(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])


def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    T = t_s / (86400.0 * 36525.0)
    sec = math.radians(1.0 / 3600.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * sec
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * sec
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * sec
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)


# --------------------------------------------------------------------------- #
# Snapshot loading (same as before)
# --------------------------------------------------------------------------- #
def _load_snapshot(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"snapshot missing: {path}")
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    soe_idx = None
    eoe_idx = None
    for i, line in enumerate(lines):
        if "$$SOE" in line:
            soe_idx = i + 1
        if "$$EOE" in line:
            eoe_idx = i
            break
    rows = []
    for line in lines[soe_idx:eoe_idx]:
        s = line.strip()
        if not s:
            continue
        parts = [p.strip() for p in s.split(",")]
        jd_tt = float(parts[0])
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        rows.append((jd_tt, x, y, z))
    arr = np.array(rows)
    t_s = (arr[:, 0] - JD_J2000) * 86400.0
    r_eci = arr[:, 1:4]
    return {
        "t_s": t_s,
        "r_eci_km": r_eci,
        "sha256": _sha256(path),
        "n_points": len(rows),
        "jd_start": float(arr[0, 0]),
        "jd_end": float(arr[-1, 0]),
        "duration_days": float(arr[-1, 0] - arr[0, 0]),
    }


def _interp_snapshot_precessed(t_query_s: float, snap: dict,
                                apply_precession: bool = True) -> np.ndarray:
    t_s = snap["t_s"]
    r = snap["r_eci_km"]
    if t_query_s <= t_s[0]:
        rv = r[0]
    elif t_query_s >= t_s[-1]:
        rv = r[-1]
    else:
        idx = int(np.searchsorted(t_s, t_query_s))
        t_lo = t_s[idx - 1]
        t_hi = t_s[idx]
        frac = (t_query_s - t_lo) / (t_hi - t_lo)
        rv = r[idx - 1] + frac * (r[idx - 1] - r[idx])
    if apply_precession:
        P = precession_j2000_to_mod(t_query_s)
        return P @ rv
    return rv


def _third_body_accel(r_eci_km: np.ndarray, t_s: float, sun_snap: dict,
                      moon_snap: dict, *, include_sun: bool = True,
                      include_moon: bool = True,
                      apply_precession: bool = True) -> np.ndarray:
    a_total = np.zeros(3)
    if include_sun:
        r_sun = _interp_snapshot_precessed(t_s, sun_snap, apply_precession)
        r_sat_to_sun = r_sun - r_eci_km
        r3s = np.linalg.norm(r_sat_to_sun)
        r3 = np.linalg.norm(r_sun)
        a_total += SOLAR_GM_KM3_S2 * (r_sat_to_sun / r3s ** 3 - r_sun / r3 ** 3)
    if include_moon:
        r_moon = _interp_snapshot_precessed(t_s, moon_snap, apply_precession)
        r_sat_to_moon = r_moon - r_eci_km
        r3s = np.linalg.norm(r_sat_to_moon)
        r3 = np.linalg.norm(r_moon)
        a_total += LUNAR_GM_KM3_S2 * (r_sat_to_moon / r3s ** 3 - r_moon / r3 ** 3)
    return a_total


# --------------------------------------------------------------------------- #
# Streaming RK4 propagation (memory-efficient)
# --------------------------------------------------------------------------- #
def propagate_streaming(sun_snap, moon_snap, *, mode: str, t0_s: float,
                          t_end_s: float, dt_s: float = DT_S,
                          subsample_every: int = 100) -> dict:
    """DEPRECATED wrapper: use `propagate_streaming_with_x0`, which requires
    the caller to supply the initial state vector `x0`. This stub is kept
    to avoid importing failures in downstream code that referenced it.
    """
    raise RuntimeError(
        "propagate_streaming is a deprecated stub. Use "
        "propagate_streaming_with_x0(sun_snap, moon_snap, x0, mode=..., "
        "t0_s=..., t_end_s=..., dt_s=..., subsample_every=...) instead.")


# --------------------------------------------------------------------------- #
# Phase-locked 2-window estimator (operates on streaming outputs)
# --------------------------------------------------------------------------- #
def phase_locked_two_window(t_cross: np.ndarray, om_cross: np.ndarray, *,
                              window_days: float, separation_days: float,
                              t_start_s: float = 0.0,
                              mode: str = "as_measured") -> dict:
    """Phase-locked 2-window RAAN estimator.

    Splits the ascending-node-crossing RAAN history into two windows of
    length `window_days` separated by `separation_days`. Computes the
    RAAN drift in each window; the average of the two drifts is the
    "phase-locked" estimator that, in principle, cancels the contribution
    of a slow harmonic whose period equals `2 * separation_days` (the
    two-window baseline).

    With `separation_days = HALF_NODAL_DAYS = 3399.2 d`, the two-window
    baseline equals one full lunar nodal cycle (6798.4 d), so the
    estimator is asymptotically insensitive to the dominant slow harmonic
    of the lunar evection/variation family.

    Args:
      t_cross: ascending-node-crossing times [seconds since some epoch].
      om_cross: unwrapped RAAN at each crossing [radians].
      window_days: length of each window [days].
      separation_days: gap between window midpoints [days].
      t_start_s: epoch of t_cross = 0 [seconds since some epoch].
      mode: 'as_measured' (default) returns each window's drift in rad/s;
            'rate_only' returns the average rate only.

    Returns:
      Dict with window A and window B drift rates (rad/day), the average
      rate, and the number of node-crossings used in each window. Returns
      NaN if either window has fewer than 4 nodes.
    """
    n = len(t_cross)
    if n < 8:
        return {
            "window_a_drift_deg_day": float("nan"),
            "window_b_drift_deg_day": float("nan"),
            "avg_drift_deg_day": float("nan"),
            "window_a_n_nodes": 0,
            "window_b_n_nodes": 0,
        }

    t_day = (t_cross - t_start_s) / 86400.0

    win_s = window_days
    sep_s = separation_days

    # Window A: [0, win_s]
    # Window B: [win_s + sep_s, 2*win_s + sep_s]
    a_lo, a_hi = 0.0, win_s
    b_lo, b_hi = win_s + sep_s, 2.0 * win_s + sep_s

    a_mask = (t_day >= a_lo) & (t_day <= a_hi)
    b_mask = (t_day >= b_lo) & (t_day <= b_hi)

    out = {
        "window_a_n_nodes": int(np.sum(a_mask)),
        "window_b_n_nodes": int(np.sum(b_mask)),
        "window_a_drift_deg_day": float("nan"),
        "window_b_drift_deg_day": float("nan"),
        "avg_drift_deg_day": float("nan"),
    }
    if out["window_a_n_nodes"] >= 4:
        _, a_rate = ols_slope(t_day[a_mask], om_cross[a_mask])
        out["window_a_drift_deg_day"] = math.degrees(a_rate)
    if out["window_b_n_nodes"] >= 4:
        _, b_rate = ols_slope(t_day[b_mask], om_cross[b_mask])
        out["window_b_drift_deg_day"] = math.degrees(b_rate)
    if (not math.isnan(out["window_a_drift_deg_day"]) and
            not math.isnan(out["window_b_drift_deg_day"])):
        out["avg_drift_deg_day"] = 0.5 * (
            out["window_a_drift_deg_day"] + out["window_b_drift_deg_day"])
    return out


def propagate_streaming_with_x0(sun_snap, moon_snap, x0: np.ndarray, *,
                                  mode: str, t0_s: float, t_end_s: float,
                                  dt_s: float = DT_S,
                                  subsample_every: int = 100) -> dict:
    """Stream RK4 step-by-step, collect ascending-node crossings +
    subsampled node-vector samples. NO full-trajectory storage.

    For an 18.6-yr arc at dt=60s (~9.8M steps), subsample_every=100 gives
    ~98k node-vector samples (similar count to ascending nodes).
    """
    f_j2 = j2_rhs(MU_EARTH_KM3S2, J2_EARTH, R_EARTH_KM)
    use_sun = mode in ("sun_only", "sun_moon", "sun_moon_j2")
    use_moon = mode in ("moon_only", "sun_moon", "sun_moon_j2")
    use_j2 = mode != "kepler_only"

    def f(t, x):
        r = x[:3]
        v = x[3:]
        if mode == "kepler_only":
            a_kep = -MU_EARTH_KM3S2 * r / np.linalg.norm(r) ** 3
            return np.concatenate([v, a_kep])
        if use_j2:
            a = f_j2(t, x)[3:]
        else:
            a = -MU_EARTH_KM3S2 * r / np.linalg.norm(r) ** 3
        if use_sun or use_moon:
            a_3b = _third_body_accel(r, t, sun_snap, moon_snap,
                                       include_sun=use_sun, include_moon=use_moon)
            a = a + a_3b
        return np.concatenate([v, a])

    t_cross_list = []
    om_cross_list = []
    t_node_list = []
    omega_node_list = []
    z_prev = x0[2]
    h = np.cross(x0[:3], x0[3:])
    n_x_prev, n_y_prev = h[1], -h[0]
    omega_node_prev = math.atan2(-h[0], h[1])
    t_node_list.append(t0_s)
    omega_node_list.append(omega_node_prev)
    x = x0.copy()
    t = t0_s
    n_steps = 0
    t_start_wall = time.time()
    while t < t_end_s:
        # RK4 step
        x_new = rk4_step(f, t, x, dt_s)
        # Ascending-node detection
        z_curr = x_new[2]
        if z_prev <= 0 < z_curr:
            frac = -z_prev / (z_curr - z_prev)
            t_cross = t + frac * dt_s
            r_cross = x[:3] + frac * (x_new[:3] - x[:3])
            om_cross = math.atan2(r_cross[1], r_cross[0])
            t_cross_list.append(t_cross)
            om_cross_list.append(om_cross)
        # Node-vector at subsampled steps
        if n_steps % subsample_every == 0:
            h = np.cross(x_new[:3], x_new[3:])
            n_x, n_y = h[1], -h[0]
            omega_node = math.atan2(-h[0], h[1])
            # Unwrap: compare to previous
            while omega_node < omega_node_prev - math.pi:
                omega_node += 2 * math.pi
            while omega_node > omega_node_prev + math.pi:
                omega_node -= 2 * math.pi
            t_node_list.append(t + dt_s)
            omega_node_list.append(omega_node)
            omega_node_prev = omega_node
        z_prev = z_curr
        x = x_new
        t = t + dt_s
        n_steps += 1
    # Final sample
    h = np.cross(x[:3], x[3:])
    omega_node_list.append(math.atan2(-h[0], h[1]))
    t_node_list.append(t)

    # Unwrap crossings using np.unwrap
    om_arr = np.array(om_cross_list)
    if len(om_arr) > 1:
        om_unwrapped = np.unwrap(om_arr)
    else:
        om_unwrapped = om_arr
    return {
        "t_cross": np.array(t_cross_list),
        "om_cross": om_unwrapped,
        "t_node": np.array(t_node_list),
        "omega_node": np.array(omega_node_list),
        "n_steps": n_steps,
        "wall_clock_s": time.time() - t_start_wall,
    }


# --------------------------------------------------------------------------- #
# Estimators (same as before, but accept streaming outputs)
# --------------------------------------------------------------------------- #
def ols_slope(t_s: np.ndarray, y_rad: np.ndarray) -> tuple:
    A = np.column_stack([np.ones_like(t_s), t_s])
    result = np.linalg.lstsq(A, y_rad, rcond=None)
    return float(result[0][0]), float(result[0][1])


HARMONIC_BASIS_PERIODS_DAYS = (
    6798.4,
    365.2422, 182.6211, 121.7474, 91.3106, 73.0484,
    27.5546, 14.7653,
    9.3067 * 365.2422,
)


def harmonic_regression(t_day: np.ndarray, omega_rad: np.ndarray) -> dict:
    n = len(t_day)
    cols = [np.ones(n), t_day]
    for T_d in HARMONIC_BASIS_PERIODS_DAYS:
        omega_k = 2.0 * math.pi / T_d
        cols.append(np.cos(omega_k * t_day))
        cols.append(np.sin(omega_k * t_day))
    A = np.column_stack(cols)
    result = np.linalg.lstsq(A, omega_rad, rcond=None)
    coeffs = result[0]
    b_rad_per_day = float(coeffs[1])
    fit_pred = A @ coeffs
    rms = math.sqrt(np.mean((omega_rad - fit_pred) ** 2))
    harmonic_coeffs = {}
    for i, T_d in enumerate(HARMONIC_BASIS_PERIODS_DAYS):
        c_k = float(coeffs[2 + 2 * i])
        s_k = float(coeffs[2 + 2 * i + 1])
        amp = math.sqrt(c_k * c_k + s_k * s_k)
        harmonic_coeffs[T_d] = {
            "cos": c_k, "sin": s_k, "amp_rad": amp,
            "amp_deg": math.degrees(amp),
        }
    return {
        "b_rad_per_day": b_rad_per_day,
        "b_deg_per_day": math.degrees(b_rad_per_day),
        "rms_residual_deg": math.degrees(rms),
        "n_points": n,
        "harmonic_amplitudes_deg": harmonic_coeffs,
    }


# --------------------------------------------------------------------------- #
# Synthetic oracle, force-level identity, idealized bridge (unchanged)
# --------------------------------------------------------------------------- #
def synthetic_oracle_test() -> dict:
    a_true_deg = 1.0e-4
    harmonics = [
        (6798.4, 0.5, 1.5),
        (365.2422, 0.05, 0.0),
        (182.6211, 0.02, 0.5),
        (121.7474, 0.01, 1.0),
        (91.3106, 0.005, 1.5),
        (73.0484, 0.003, 2.0),
        (27.5546, 0.002, 0.3),
        (14.7653, 0.001, 1.1),
    ]
    W_days = 18.6 * 365.25
    n_samples = int(W_days * 5.6)
    t_day = np.linspace(0, W_days, n_samples)
    a_true_rad_per_day = math.radians(a_true_deg)
    omega_rad = a_true_rad_per_day * t_day
    for T_d, A_deg, phi in harmonics:
        A_rad = math.radians(A_deg)
        omega_k = 2 * math.pi / T_d
        omega_rad += A_rad * np.sin(omega_k * t_day + phi)
    fit = harmonic_regression(t_day, omega_rad)
    bias_deg_day = fit["b_deg_per_day"] - a_true_deg
    _, b_ols_rad = ols_slope(t_day, omega_rad)
    b_ols_deg_day = math.degrees(b_ols_rad)
    bias_ols_deg_day = b_ols_deg_day - a_true_deg
    return {
        "a_true_deg_day": a_true_deg,
        "W_days": W_days,
        "n_samples": n_samples,
        "estimator_f_harmonic_regression_deg_day": fit["b_deg_per_day"],
        "estimator_f_bias_deg_day": bias_deg_day,
        "estimator_f_rms_residual_deg": fit["rms_residual_deg"],
        "estimator_a_direct_ols_deg_day": b_ols_deg_day,
        "estimator_a_bias_deg_day": bias_ols_deg_day,
        "harmonic_amplitudes_deg_recovered": fit["harmonic_amplitudes_deg"],
        "verdict": "f wins" if abs(bias_deg_day) < abs(bias_ols_deg_day) else "a wins",
    }


def force_level_identity_check(h_km: float = H_SSO_KM, n_states: int = 50,
                                seed: int = 42) -> dict:
    rng = np.random.default_rng(seed)
    a = R_EARTH_KM + h_km
    r3_sun = np.array([1.0, 0.0, 0.0]) * AU_KM
    r3_moon = np.array([0.0, 1.0, 0.0]) * LUNAR_DISTANCE_KM_MEAN
    max_diff_sun = 0.0
    max_diff_moon = 0.0
    for _ in range(n_states):
        v = rng.standard_normal(3)
        v /= np.linalg.norm(v)
        r_sat = a * v
        r_sat_to_sun = r3_sun - r_sat
        r3s_sun = np.linalg.norm(r_sat_to_sun)
        r3_mag_sun = np.linalg.norm(r3_sun)
        a_sun_a = SOLAR_GM_KM3_S2 * (r_sat_to_sun / r3s_sun ** 3 - r3_sun / r3_mag_sun ** 3)
        a_sun_b = SOLAR_GM_KM3_S2 * (-(r_sat - r3_sun) / r3s_sun ** 3 - r3_sun / r3_mag_sun ** 3)
        r_sat_to_moon = r3_moon - r_sat
        r3s_moon = np.linalg.norm(r_sat_to_moon)
        r3_mag_moon = np.linalg.norm(r3_moon)
        a_moon_a = LUNAR_GM_KM3_S2 * (r_sat_to_moon / r3s_moon ** 3 - r3_moon / r3_mag_moon ** 3)
        a_moon_b = LUNAR_GM_KM3_S2 * (-(r_sat - r3_moon) / r3s_moon ** 3 - r3_moon / r3_mag_moon ** 3)
        max_diff_sun = max(max_diff_sun, float(np.max(np.abs(a_sun_a - a_sun_b))))
        max_diff_moon = max(max_diff_moon, float(np.max(np.abs(a_moon_a - a_moon_b))))
    return {
        "n_states": n_states,
        "max_diff_sun_km_s2": float(max_diff_sun),
        "max_diff_moon_km_s2": float(max_diff_moon),
        "passes_sun": max_diff_sun < 1e-15,
        "passes_moon": max_diff_moon < 1e-15,
    }


def idealized_circular_perturber_bridge(h_km: float, i_deg: float,
                                          *, n_states: int = 4) -> dict:
    a = R_EARTH_KM + h_km
    n = mean_motion(a)
    i_rad = math.radians(i_deg)
    r0 = np.array([a, 0.0, 0.0])
    v_circ = math.sqrt(MU_EARTH_KM3S2 / a)
    v0 = np.array([0.0, v_circ * math.cos(i_rad), v_circ * math.sin(i_rad)])
    rng = np.random.default_rng(42)
    i3_moon_rad = math.radians(SOLAR_OBLIQUITY_DEG + LUNAR_INCLINATION_DEG)
    omega_samples = np.zeros(256)
    for k in range(256):
        f_sat = 2 * math.pi * k / 256
        r_sat = r0 * math.cos(f_sat) + (v0 / v_circ) * a * math.sin(f_sat)
        M_moon = 2 * math.pi * rng.random()
        r_moon_ecl = np.array([
            LUNAR_DISTANCE_KM_MEAN * math.cos(M_moon),
            LUNAR_DISTANCE_KM_MEAN * math.sin(M_moon),
            0.0,
        ])
        c, s = math.cos(i3_moon_rad), math.sin(i3_moon_rad)
        R_ecl_to_eq = np.array([
            [1.0, 0.0, 0.0],
            [0.0, c, -s],
            [0.0, s, c],
        ])
        r_moon_eq = R_ecl_to_eq @ r_moon_ecl
        r_sat_to_moon = r_moon_eq - r_sat
        r3s_moon = np.linalg.norm(r_sat_to_moon)
        r3_mag_moon = np.linalg.norm(r_moon_eq)
        a_3b = LUNAR_GM_KM3_S2 * (r_sat_to_moon / r3s_moon ** 3 - r_moon_eq / r3_mag_moon ** 3)
        h = np.cross(r_sat, v0)
        h_mag = np.linalg.norm(h)
        h_hat = h / h_mag
        n_hat = np.array([-h_hat[1], h_hat[0], 0.0])
        n_hat /= np.linalg.norm(n_hat)
        a_nodal = np.dot(a_3b, n_hat)
        omega_samples[k] = a_nodal
    mean_a_nodal = float(np.mean(omega_samples))
    mean_a_nodal_deg_day = math.degrees(mean_a_nodal) * 86400.0
    cf_lunar_rad = (3.0 / 8.0) * n * (LUNAR_GM_KM3_S2 / MU_EARTH_KM3S2) * (
        a / LUNAR_DISTANCE_KM_MEAN) ** 3 * math.sin(2.0 * (i_rad - i3_moon_rad)) / math.sin(i_rad)
    cf_lunar_deg_day = math.degrees(cf_lunar_rad) * 86400.0
    return {
        "i_deg": i_deg,
        "idealized_orbit_averaged_nodal_deg_day": mean_a_nodal_deg_day,
        "cf_lunar_component_deg_day": cf_lunar_deg_day,
        "ratio": mean_a_nodal_deg_day / cf_lunar_deg_day if cf_lunar_deg_day != 0 else float("nan"),
    }


def corrected_secular_lunisolar_raan_rate_rad_s(
    h_km: float, i_deg: float,
    *,
    i3_sun_rad: float = math.radians(SOLAR_OBLIQUITY_DEG),
    i3_moon_rad: float = math.radians(SOLAR_OBLIQUITY_DEG + LUNAR_INCLINATION_DEG),
    mu_3_sun: float = SOLAR_GM_KM3_S2,
    mu_3_moon: float = LUNAR_GM_KM3_S2,
    a_3_sun: float = AU_KM,
    a_3_moon: float = LUNAR_DISTANCE_KM_MEAN,
) -> dict:
    a = R_EARTH_KM + h_km
    n = mean_motion(a)
    i_rad = math.radians(i_deg)
    solar = (3.0 / 8.0) * n * (mu_3_sun / MU_EARTH_KM3S2) * (
        a / a_3_sun) ** 3 * math.sin(2.0 * (i_rad - i3_sun_rad)) / math.sin(i_rad)
    lunar = (3.0 / 8.0) * n * (mu_3_moon / MU_EARTH_KM3S2) * (
        a / a_3_moon) ** 3 * math.sin(2.0 * (i_rad - i3_moon_rad)) / math.sin(i_rad)
    total = solar + lunar
    return {
        "h_km": h_km,
        "i_deg": i_deg,
        "solar_deg_day": math.degrees(solar) * 86400.0,
        "lunar_deg_day": math.degrees(lunar) * 86400.0,
        "total_deg_day": math.degrees(total) * 86400.0,
    }


def mean_motion(a_km: float) -> float:
    return math.sqrt(MU_EARTH_KM3S2 / a_km ** 3)


def code_hashes() -> dict:
    here = Path(__file__).resolve().parent
    # Walk up to find the repo root (the directory that contains
    # `src/lab_utils/`). This is robust to the mission's location.
    lab_root = here
    while lab_root != lab_root.parent:
        if (lab_root / "src" / "lab_utils").exists():
            break
        lab_root = lab_root.parent
    files = {
        "experiment.py": here / "experiment.py",
        "lab_utils/orbits.py": lab_root / "src" / "lab_utils" / "orbits.py",
        "lab_utils/earth_frames.py": lab_root / "src" / "lab_utils" / "earth_frames.py",
        "lab_utils/integrators.py": lab_root / "src" / "lab_utils" / "integrators.py",
        "lab_utils/results.py": lab_root / "src" / "lab_utils" / "results.py",
        "lab_utils/__init__.py": lab_root / "src" / "lab_utils" / "__init__.py",
        "018_experiment.py": lab_root / "research" / "orbital-mechanics" / "experiments" / "lunisolarReconciliation" / "experiment.py",
        "020_experiment.py": lab_root / "research" / "orbital-mechanics" / "experiments" / "lunisolarSecularLimit" / "experiment.py",
    }
    files["sun_reference_snapshot.txt"] = SUN_SNAPSHOT
    files["moon_reference_snapshot.txt"] = MOON_SNAPSHOT
    files["MANIFEST.json"] = MANIFEST_PATH
    return {name: _sha256(p) for name, p in files.items()}