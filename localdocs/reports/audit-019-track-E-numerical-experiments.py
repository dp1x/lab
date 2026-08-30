"""Audit-019 Track E: numerical-experiment suite for Lunisolar RAAN convergence.

INDEPENDENT reimplementation of the 018 numerical experiments, focused on
SECULAR/PERIODIC/FINITE-WINDOW separation. Reads ONLY:
  - AGENTS.md, localdocs/roadmap.md (project context)
  - Exp 017/018 experiment.py (read for context only; reuses
    lab_utils canon + snapshot paths; does NOT modify the donor files)
  - src/lab_utils/integrators.py (rk4_propagate)
  - src/lab_utils/orbits.py (j2_rhs, sso_inclination_rad, mean_motion)
  - byte-pinned Sun and Moon snapshot files

DOES NOT READ audit-018 or any other track's outputs.

Scientific question (Track E):
  How does the 1-year numerical linear-fit RAAN rate depend on the
  ESTIMATOR used? The window-shift (first half vs second half),
  phase bias (full-moon start vs new-moon start), and cycle-averaged
  estimator all probe the secular/periodic/finite-window decomposition.

Forces (4 modes at h=600 km, i=i_sso=97.79 deg, 1-year arc):
  1. Sun only
  2. Moon only
  3. Sun + Moon (no J2)
  4. Sun + Moon + J2 (full model)

Estimators applied to each Omega(t) time series:
  - full year linear fit
  - first half (W=180 d) linear fit
  - second half (W=180 d, shifted by 180 d) linear fit
  - full-moon start window (W=180 d) linear fit
  - new-moon start window (W=180 d) linear fit
  - perigee start window (W=180 d) linear fit
  - W=27.55 d (anomalistic month) linear fit, anchored at start
  - W=29.53 d (synodic month) linear fit, anchored at start
  - cycle-averaged estimator (12 synodic-month segments, slope per
    segment, then mean and std)

Headline validation: compare to 018 published result:
  - J2-only at h=600 km i_sso: ~+0.9920 deg/day
  - Full model at h=600 km i_sso: ~+0.9933 deg/day
  - Lunisolar-only (full - J2): ~+1.32e-3 deg/day (numerical prograde)

Independent implementation choices:
  - Sun/Moon acceleration: third-body direct + indirect form
    a_3 = mu_3 * (r_3 - r_sat)/|r_3 - r_sat|^3 - mu_3 * r_3/|r_3|^3
  - Frame: NO precession rotation applied (Track E pure precession-off
    baseline; the 018 'with_precession' fix is a Track D deliverable, not
    re-derived here. Documented as a methodological choice).
  - dt = 60 s (RK4 LEO/SSO convergence proven by 017/018).
  - Snapshot interpolation: linear-in-time on the (X, Y, Z) components
    independently; consistent with 017.
  - J2 uses the lab's graduated j2_rhs (donor: Exp 009).

Deterministic: pure float64, no RNG, no network at runtime, fixed
epoch = 820540800 s (2026-01-01 12:00 TT). Two consecutive runs
produce identical results.

References (canonical, no audit-018):
  - Vallado 4th ed., Ch. 9 (secular J2 + Lunisolar, Eq. 9-46 form).
  - Exp 009 (j2Precession): J2 secular nodal rate
    dO/dt = -1.5 J2 sqrt(mu_E) R^2 cos(i) / a^(7/2) (1-e^2)^-2
  - Exp 012 (orbitClasses): SSO inclination lock, a_max = 12352.505 km.
  - Exp 014 (eclipseTiming): byte-pinned 2026 Sun snapshot.
  - Exp 017 (lunisolarVerification): byte-pinned 2026 Moon snapshot,
    1-year RK4 propagation at dt=60 s, ascending-node detection.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

# Ensure lab_utils is importable when invoked directly outside the lab's
# default Python environment. The repo's `src/` directory contains the
# lab_utils package; we add it to sys.path if it is not already resolvable.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_LAB_UTILS_PARENT = _REPO_ROOT / "src"
if str(_LAB_UTILS_PARENT) not in sys.path:
    sys.path.insert(0, str(_LAB_UTILS_PARENT))

import numpy as np

from lab_utils import (
    J2_EARTH,
    MU_EARTH_KM3S2,
    R_EARTH_KM,
    SSO_TARGET_DEG_DAY,
    j2_rhs,
    mean_motion,
    sso_inclination_rad,
)
from lab_utils.earth_frames import JD_J2000
from lab_utils.integrators import rk4_propagate

# --------------------------------------------------------------------------- #
# Constants (frozen)
# --------------------------------------------------------------------------- #
EXP_TAG = "audit-019-track-E"
ALT_H_KM = 600.0
DT_S = 60.0
T0_S = 820540800.0  # 2026-01-01 12:00 TT, lab convention

SOLAR_GM_KM3_S2 = 132712440018.0
LUNAR_GM_KM3_S2 = 4902.8001
AU_KM = 149597870.7

# Lunar phase dates in 2026 (TT-like, ISO calendar dates). Source: JPL
# Horizons ephemeris approximations. Year-of-2026 list:
#   - first new moon: 2026-01-03
#   - first full moon: 2026-01-19
#   - first perigee: 2026-01-07
# Each subsequent phase advances by ~1 synodic month (29.53 d) for new/full,
# ~1 anomalistic month (27.55 d) for perigee.
LUNAR_PHASES_2026 = {
    "new_moons": [
        "2026-01-03", "2026-02-01", "2026-03-03", "2026-04-02",
        "2026-05-01", "2026-05-31", "2026-06-29", "2026-07-29",
        "2026-08-28", "2026-09-26", "2026-10-26", "2026-11-24",
        "2026-12-24",
    ],
    "full_moons": [
        "2026-01-19", "2026-02-17", "2026-03-19", "2026-04-17",
        "2026-05-17", "2026-06-15", "2026-07-15", "2026-08-13",
        "2026-09-12", "2026-10-12", "2026-11-11", "2026-12-10",
    ],
    "perigees": [
        "2026-01-07", "2026-02-04", "2026-03-04", "2026-04-01",
        "2026-04-29", "2026-05-27", "2026-06-25", "2026-07-23",
        "2026-08-21", "2026-09-18", "2026-10-15", "2026-11-12",
        "2026-12-11",
    ],
}

# Lunar periods (canonical)
ANOMALISTIC_MONTH_D = 27.554550  # Meeus Ch. 53
SYNODIC_MONTH_D = 29.530588      # Meeus Ch. 53
DRACONIC_MONTH_D = 27.212221     # lunar nodal period (relevant to short-period)

# Window lengths
W_FULL_D = 365.0
W_HALF_D = 180.0
W_PHASE_D = 180.0  # full/new/perigee-anchored window
N_CYCLE_SEGMENTS = 12  # 12 synodic months per year (synodic ~= 29.53 d * 12 ~ 354 d)

# Snapshot paths (byte-pinned)  (REPO_ROOT defined above)
SUN_SNAPSHOT_PATH = (
    _REPO_ROOT / "research" / "orbital-mechanics" / "experiments"
    / "eclipseTiming" / "reference"
    / "horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt"
)
MOON_SNAPSHOT_PATH = (
    _REPO_ROOT / "research" / "orbital-mechanics" / "experiments"
    / "lunisolarVerification" / "reference"
    / "horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt"
)

OUT_RESULTS_JSON = (
    _REPO_ROOT / "localdocs" / "reports"
    / "audit-019-track-E-numerical-experiments-results.json"
)
OUT_REPORT_MD = (
    _REPO_ROOT / "localdocs" / "reports"
    / "audit-019-track-E-numerical-experiments-report.md"
)


# --------------------------------------------------------------------------- #
# Snapshot loading
# --------------------------------------------------------------------------- #
def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_snapshot(path: Path) -> dict:
    """Load byte-pinned JPL Horizons geocentric vector snapshot.

    Parses $$SOE/$$EOE block. Returns dict with t_s (seconds since J2000
    TDB) and r_eci_km (geocentric vector in ICRF, km).
    """
    if not path.exists():
        raise FileNotFoundError(f"snapshot missing: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    soe = eoe = None
    for i, line in enumerate(lines):
        if "$$SOE" in line:
            soe = i + 1
        if "$$EOE" in line:
            eoe = i
            break
    rows = []
    for line in lines[soe:eoe]:
        s = line.strip()
        if not s:
            continue
        parts = [p.strip() for p in s.split(",")]
        # Columns: JDTDB, Calendar Date, X, Y, Z, VX, VY, VZ
        jd_tt = float(parts[0])
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        rows.append((jd_tt, x, y, z))
    arr = np.array(rows)
    return {
        "t_s": (arr[:, 0] - JD_J2000) * 86400.0,
        "r_eci_km": arr[:, 1:4],
        "sha256": _file_sha256(path),
        "n_points": len(rows),
    }


def _interp_snapshot(t_query_s: float, snap: dict) -> np.ndarray:
    """Linear interpolation of geocentric vector at query time.

    Clamps to endpoint values outside the snapshot range (known limitation
    disclosed in 017). Operates on ICRF/J2000 as stored.
    """
    t_s = snap["t_s"]
    r = snap["r_eci_km"]
    if t_query_s <= t_s[0]:
        return r[0]
    if t_query_s >= t_s[-1]:
        return r[-1]
    idx = int(np.searchsorted(t_s, t_query_s))
    t_lo = t_s[idx - 1]
    t_hi = t_s[idx]
    frac = (t_query_s - t_lo) / (t_hi - t_lo)
    return r[idx - 1] + frac * (r[idx] - r[idx - 1])


# --------------------------------------------------------------------------- #
# Third-body acceleration (independent reimplementation)
# --------------------------------------------------------------------------- #
def _third_body_accel(
    r_eci_km: np.ndarray,
    t_s: float,
    sun_snap: dict,
    moon_snap: dict,
    *,
    include_sun: bool,
    include_moon: bool,
) -> np.ndarray:
    """Geocentric third-body acceleration.

    a_3 = mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3 - mu_3 * r_3 / |r_3|^3

    First term: direct attraction toward the third body.
    Second term: indirect (Earth's acceleration toward the third body,
    transferred to the geocentric frame; O(1e-5) of direct term at LEO).

    Frame choice (Track E methodological decision): NO precession rotation
    applied. The snapshot is ICRF/J2000; the propagator uses the vector
    in the same ICRF frame (no rotation to mean-of-date). The 018
    'with_precession' fix adds a ~0.4 deg rotation that improves the
    frame match; Track E uses the raw ICRF/J2000 snapshot consistently
    (017 behavior).
    """
    a = np.zeros(3)
    if include_sun:
        r_sun = _interp_snapshot(t_s, sun_snap)
        d_sun = r_sun - r_eci_km
        a = a + SOLAR_GM_KM3_S2 * (
            d_sun / np.linalg.norm(d_sun) ** 3
            - r_sun / np.linalg.norm(r_sun) ** 3
        )
    if include_moon:
        r_moon = _interp_snapshot(t_s, moon_snap)
        d_moon = r_moon - r_eci_km
        a = a + LUNAR_GM_KM3_S2 * (
            d_moon / np.linalg.norm(d_moon) ** 3
            - r_moon / np.linalg.norm(r_moon) ** 3
        )
    return a


# --------------------------------------------------------------------------- #
# Force-mode RHS builders (independent)
# --------------------------------------------------------------------------- #
def make_rhs(sun_snap: dict, moon_snap: dict, *, mode: str):
    """Build RHS for a given mode.

    mode in {"sun_only", "moon_only", "sun_moon", "sun_moon_j2"}
    """
    f_j2 = j2_rhs(MU_EARTH_KM3S2, J2_EARTH, R_EARTH_KM)

    def f(t: float, x: np.ndarray) -> np.ndarray:
        r = x[:3]
        v = x[3:]
        # Kepler + J2 from graduated canon (NO Lunisolar for j2_only; this
        # experiment's j2_only mode is covered by 'sun_moon_j2' minus Sun/Moon).
        if mode == "sun_moon_j2":
            a_kep_j2 = f_j2(t, x)[3:]
        else:
            # No-J2 cases: pure Kepler
            a_kep_j2 = -MU_EARTH_KM3S2 * r / np.linalg.norm(r) ** 3
        if mode == "sun_only":
            a_3 = _third_body_accel(r, t, sun_snap, moon_snap,
                                     include_sun=True, include_moon=False)
        elif mode == "moon_only":
            a_3 = _third_body_accel(r, t, sun_snap, moon_snap,
                                     include_sun=False, include_moon=True)
        elif mode == "sun_moon":
            a_3 = _third_body_accel(r, t, sun_snap, moon_snap,
                                     include_sun=True, include_moon=True)
        elif mode == "sun_moon_j2":
            a_3 = _third_body_accel(r, t, sun_snap, moon_snap,
                                     include_sun=True, include_moon=True)
        else:
            raise ValueError(f"unknown mode: {mode}")
        return np.concatenate([v, a_kep_j2 + a_3])

    return f


# --------------------------------------------------------------------------- #
# Ascending-node detection (independent implementation, same algorithm as 017)
# --------------------------------------------------------------------------- #
def detect_ascending_nodes(t_s_arr: np.ndarray, x_arr: np.ndarray):
    """Detect ascending-node crossings (z=0, vz>0) with linear interpolation.

    Returns (t_cross, om_cross_rad). Omega is unwrapped against the prior
    crossing to give a continuous record (the linear fit uses the raw
    radian values mod 2pi is handled by the unwrap).
    """
    t_crossings = []
    om_crossings = []
    z_prev = x_arr[0, 2]
    vz_prev = x_arr[0, 5]
    for k in range(1, len(t_s_arr)):
        z_curr = x_arr[k, 2]
        vz_curr = x_arr[k, 5]
        if z_prev <= 0 < z_curr and vz_prev > 0:
            frac = -z_prev / (z_curr - z_prev)
            t_cross = t_s_arr[k - 1] + frac * (t_s_arr[k] - t_s_arr[k - 1])
            r_cross = x_arr[k - 1, :3] + frac * (x_arr[k, :3] - x_arr[k - 1, :3])
            om_cross = math.atan2(r_cross[1], r_cross[0])
            if om_crossings:
                om_prev = om_crossings[-1]
                while om_cross < om_prev - math.pi:
                    om_cross += 2 * math.pi
                while om_cross > om_prev + math.pi:
                    om_cross -= 2 * math.pi
            t_crossings.append(t_cross)
            om_crossings.append(om_cross)
        z_prev = z_curr
        vz_prev = vz_curr
    return np.array(t_crossings), np.array(om_crossings)


def linear_fit_drift(t_days: np.ndarray, y_rad: np.ndarray):
    """Linear least-squares fit y = a + b*t; returns (intercept, slope, rms_deg)."""
    A = np.column_stack([np.ones_like(t_days), t_days])
    coeffs, *_ = np.linalg.lstsq(A, y_rad, rcond=None)
    intercept = float(coeffs[0])
    slope = float(coeffs[1])
    pred = slope * t_days + intercept
    rms_deg = math.degrees(math.sqrt(np.mean((y_rad - pred) ** 2)))
    return intercept, slope, rms_deg


# --------------------------------------------------------------------------- #
# Lunar phase date -> seconds since J2000 TT
# --------------------------------------------------------------------------- #
def _date_to_tt_seconds(date_str: str) -> float:
    """Convert 'YYYY-MM-DD' to seconds since J2000 TT (noon-based).

    Uses the lab's simple convention: noon-of-date. Differences relative
    to T0_S matter, not absolute values. JD at noon-UT1 maps to TT-like
    via TT-UTC = 69.184 s. Good to ~minute precision over 2026; window-
    start sensitivity is dominated by the 180-d estimator, not by
    minute-precision on the phase anchor.
    """
    y, m, d = date_str.split("-")
    y, m, d = int(y), int(m), int(d)
    # Julian Day of calendar date at noon (Gregorian)
    # Algorithm from Meeus Ch. 7 (low precision)
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4  # Gregorian correction
    jd = math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1)) + d + B - 1524.5
    # TT = UT1 + 69.184 s (declared lab convention)
    t_tt_s = (jd - JD_J2000) * 86400.0 + 69.184
    return t_tt_s


# --------------------------------------------------------------------------- #
# Propagation per (mode)
# --------------------------------------------------------------------------- #
def propagate_mode(mode: str, sun_snap: dict, moon_snap: dict,
                   duration_days: float = W_FULL_D,
                   dt_s: float = DT_S) -> dict:
    """Propagate the chosen force mode for the chosen duration.

    Initial conditions: circular SSO at h=ALT_H_KM, RAAN=0, at ascending
    node with inclination = i_sso (closed-form retrograde). r = (a, 0, 0),
    v = (0, v_circ*cos(i), v_circ*sin(i)) (heading north).
    """
    a = R_EARTH_KM + ALT_H_KM
    e = 0.0
    i_sso = sso_inclination_rad(a, e)
    n = mean_motion(a)
    v_circ = math.sqrt(MU_EARTH_KM3S2 / a)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, v_circ * math.cos(i_sso), v_circ * math.sin(i_sso)])
    x0 = np.concatenate([r0, v0])

    t0 = T0_S
    t_end = t0 + duration_days * 86400.0
    n_steps = int(math.ceil((t_end - t0) / dt_s))
    t_grid = np.linspace(t0, t_end, n_steps + 1)

    f = make_rhs(sun_snap, moon_snap, mode=mode)
    x_traj = rk4_propagate(f, t_grid, x0)

    t_cross, om_cross = detect_ascending_nodes(t_grid, x_traj)

    return {
        "mode": mode,
        "duration_days": duration_days,
        "dt_s": dt_s,
        "n_dt_steps": n_steps,
        "i_sso_deg": math.degrees(i_sso),
        "h_km": ALT_H_KM,
        "n_ascending_nodes": len(t_cross),
        "t_cross_s": t_cross.tolist(),
        "om_cross_rad": om_cross.tolist(),
        "t_grid_first_s": float(t_grid[0]),
        "t_grid_last_s": float(t_grid[-1]),
    }


# --------------------------------------------------------------------------- #
# Estimators
# --------------------------------------------------------------------------- #
def _slice_window(prop: dict, t_start_d: float, t_end_d: float) -> tuple:
    """Slice the Omega time series to the window [t_start_d, t_end_d] days.

    t_start_d, t_end_d are days RELATIVE to T0_S (i.e., (t - T0_S)/86400).
    Returns (t_d_local, om_rad) within the window, or (None, None) if
    fewer than 3 crossings fall inside the window.
    """
    t_cross_s = np.array(prop["t_cross_s"])
    om = np.array(prop["om_cross_rad"])
    t_cross_rel_d = (t_cross_s - T0_S) / 86400.0
    mask = (t_cross_rel_d >= t_start_d) & (t_cross_rel_d <= t_end_d)
    if mask.sum() < 3:
        return None, None
    return t_cross_rel_d[mask], om[mask]


def estimator_full_year(prop: dict) -> dict:
    t, om = _slice_window(prop, 0.0, W_FULL_D)
    if t is None:
        return {"slope_deg_per_day": float("nan"), "n_points": 0}
    intercept, slope, rms = linear_fit_drift(t, om)
    return {
        "estimator": "full_year",
        "window_days": [0.0, W_FULL_D],
        "n_points": int(len(t)),
        "slope_rad_per_day": slope,
        "slope_deg_per_day": math.degrees(slope),
        "intercept_rad": intercept,
        "rms_deg": rms,
    }


def estimator_first_half(prop: dict) -> dict:
    t, om = _slice_window(prop, 0.0, W_HALF_D)
    if t is None:
        return {"slope_deg_per_day": float("nan"), "n_points": 0}
    intercept, slope, rms = linear_fit_drift(t, om)
    return {
        "estimator": "first_half",
        "window_days": [0.0, W_HALF_D],
        "n_points": int(len(t)),
        "slope_rad_per_day": slope,
        "slope_deg_per_day": math.degrees(slope),
        "intercept_rad": intercept,
        "rms_deg": rms,
    }


def estimator_second_half(prop: dict) -> dict:
    t, om = _slice_window(prop, W_HALF_D, W_FULL_D)
    if t is None:
        return {"slope_deg_per_day": float("nan"), "n_points": 0}
    # Re-zero time axis so intercept is meaningful
    intercept, slope, rms = linear_fit_drift(t - W_HALF_D, om)
    return {
        "estimator": "second_half",
        "window_days": [W_HALF_D, W_FULL_D],
        "n_points": int(len(t)),
        "slope_rad_per_day": slope,
        "slope_deg_per_day": math.degrees(slope),
        "intercept_rad": intercept,
        "rms_deg": rms,
    }


def _phase_anchor_estimator(prop: dict, phase_label: str, phase_dates: list,
                             w_days: float = W_PHASE_D) -> dict:
    """Window starts at the first phase event in 2026.

    The window is [t_anchor, t_anchor + w_days] days relative to T0_S.
    """
    t_anchor_s = _date_to_tt_seconds(phase_dates[0])
    t_anchor_d = (t_anchor_s - T0_S) / 86400.0
    # Clamp to [0, W_FULL_D - 1] if anchor is outside the 1-year window
    if t_anchor_d < 0.0:
        t_anchor_d = 0.0
    t_end_d = min(t_anchor_d + w_days, W_FULL_D)
    t, om = _slice_window(prop, t_anchor_d, t_end_d)
    if t is None:
        return {"slope_deg_per_day": float("nan"), "n_points": 0}
    intercept, slope, rms = linear_fit_drift(t - t_anchor_d, om)
    return {
        "estimator": f"{phase_label}_start_window",
        "window_days": [float(t_anchor_d), float(t_end_d)],
        "n_points": int(len(t)),
        "slope_rad_per_day": slope,
        "slope_deg_per_day": math.degrees(slope),
        "intercept_rad": intercept,
        "rms_deg": rms,
    }


def estimator_full_moon_start(prop: dict) -> dict:
    return _phase_anchor_estimator(prop, "full_moon", LUNAR_PHASES_2026["full_moons"])


def estimator_new_moon_start(prop: dict) -> dict:
    return _phase_anchor_estimator(prop, "new_moon", LUNAR_PHASES_2026["new_moons"])


def estimator_perigee_start(prop: dict) -> dict:
    return _phase_anchor_estimator(prop, "perigee", LUNAR_PHASES_2026["perigees"])


def estimator_month_window(prop: dict, w_days: float, label: str) -> dict:
    """Window [0, w_days] linear fit (anchored at start)."""
    t_end_d = min(w_days, W_FULL_D)
    t, om = _slice_window(prop, 0.0, t_end_d)
    if t is None:
        return {"slope_deg_per_day": float("nan"), "n_points": 0}
    intercept, slope, rms = linear_fit_drift(t, om)
    return {
        "estimator": label,
        "window_days": [0.0, float(t_end_d)],
        "n_points": int(len(t)),
        "slope_rad_per_day": slope,
        "slope_deg_per_day": math.degrees(slope),
        "intercept_rad": intercept,
        "rms_deg": rms,
    }


def estimator_anomalistic_month(prop: dict) -> dict:
    return estimator_month_window(prop, ANOMALISTIC_MONTH_D, "W_anomalistic_month")


def estimator_synodic_month(prop: dict) -> dict:
    return estimator_month_window(prop, SYNODIC_MONTH_D, "W_synodic_month")


def estimator_cycle_averaged(prop: dict, n_segments: int = N_CYCLE_SEGMENTS) -> dict:
    """Divide the year into N segments and compute the slope in each.

    Segment width = W_FULL_D / N. Each segment is a separate linear fit
    against (t - t_seg_start). Returns per-segment slopes and the
    cycle-averaged (mean and std) summary.
    """
    seg_w_d = W_FULL_D / n_segments
    t_cross_s = np.array(prop["t_cross_s"])
    om = np.array(prop["om_cross_rad"])
    t_rel_d = (t_cross_s - T0_S) / 86400.0
    slopes_rad = []
    intercepts_rad = []
    seg_centers_d = []
    for i in range(n_segments):
        t0 = i * seg_w_d
        t1 = (i + 1) * seg_w_d
        mask = (t_rel_d >= t0) & (t_rel_d < t1)
        if mask.sum() < 3:
            slopes_rad.append(float("nan"))
            intercepts_rad.append(float("nan"))
            seg_centers_d.append(0.5 * (t0 + t1))
            continue
        intercept, slope, _rms = linear_fit_drift(t_rel_d[mask] - t0, om[mask])
        slopes_rad.append(slope)
        intercepts_rad.append(intercept)
        seg_centers_d.append(0.5 * (t0 + t1))
    arr = np.array(slopes_rad)
    valid = np.isfinite(arr)
    return {
        "estimator": "cycle_averaged",
        "n_segments": n_segments,
        "segment_width_days": seg_w_d,
        "seg_centers_days": seg_centers_d,
        "slopes_rad_per_day": slopes_rad,
        "slopes_deg_per_day": [math.degrees(s) if math.isfinite(s) else float("nan")
                                for s in slopes_rad],
        "mean_slope_rad_per_day": float(np.mean(arr[valid])) if valid.any() else float("nan"),
        "mean_slope_deg_per_day": float(np.degrees(np.mean(arr[valid]))) if valid.any() else float("nan"),
        "std_slope_rad_per_day": float(np.std(arr[valid])) if valid.any() else float("nan"),
        "std_slope_deg_per_day": float(np.degrees(np.std(arr[valid]))) if valid.any() else float("nan"),
    }


def apply_all_estimators(prop: dict) -> dict:
    """Apply every estimator to a single propagation result."""
    return {
        "full_year": estimator_full_year(prop),
        "first_half": estimator_first_half(prop),
        "second_half": estimator_second_half(prop),
        "full_moon_start": estimator_full_moon_start(prop),
        "new_moon_start": estimator_new_moon_start(prop),
        "perigee_start": estimator_perigee_start(prop),
        "W_anomalistic_month": estimator_anomalistic_month(prop),
        "W_synodic_month": estimator_synodic_month(prop),
        "cycle_averaged": estimator_cycle_averaged(prop),
    }


# --------------------------------------------------------------------------- #
# Residual summary (window-shift, phase, month-cycle biases)
# --------------------------------------------------------------------------- #
def residual_summary(estimators_by_mode: dict) -> dict:
    """Compute residuals between estimators within each mode.

    Returns the residual (deg/day) of each estimator vs. full_year
    baseline, plus the cycle-averaged mean vs. full_year.
    """
    out = {}
    for mode, est in estimators_by_mode.items():
        baseline = est["full_year"]["slope_deg_per_day"]
        res = {}
        for name, e in est.items():
            # cycle_averaged uses 'mean_slope_deg_per_day'; others use 'slope_deg_per_day'
            if name == "cycle_averaged":
                slope = e.get("mean_slope_deg_per_day", float("nan"))
            else:
                slope = e.get("slope_deg_per_day", float("nan"))
            if math.isfinite(slope):
                res[name] = slope - baseline
            else:
                res[name] = float("nan")
        out[mode] = {
            "residuals_vs_full_year_deg_per_day": res,
            "full_year_slope_deg_per_day": baseline,
            "cycle_mean_minus_full_year_deg_per_day": (
                est["cycle_averaged"]["mean_slope_deg_per_day"] - baseline
            ) if math.isfinite(est["cycle_averaged"]["mean_slope_deg_per_day"]) else float("nan"),
        }
    return out


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    print(f"[track-E] starting Lunisolar RAAN convergence numerical experiments")
    sun_snap = _parse_snapshot(SUN_SNAPSHOT_PATH)
    moon_snap = _parse_snapshot(MOON_SNAPSHOT_PATH)
    print(f"[track-E] Sun snapshot: {sun_snap['n_points']} rows "
          f"sha256={sun_snap['sha256'][:16]}")
    print(f"[track-E] Moon snapshot: {moon_snap['n_points']} rows "
          f"sha256={moon_snap['sha256'][:16]}")

    # ------------------------------------------------------------------------ #
    # Headline validation: reproduce the 018 baseline (j2_only via separate
    # propagation; sun_moon_j2 - sun_moon gives the Lunisolar contribution).
    # We do NOT modify 018's 'mode' names; instead we report the
    # post-publication headline check at h=600 km i_sso.
    # ------------------------------------------------------------------------ #
    print("[track-E] propagating 4 force modes at h=600 km, i=i_sso, 1-year arc...")
    modes = ["sun_only", "moon_only", "sun_moon", "sun_moon_j2"]
    propagations = {}
    for mode in modes:
        print(f"[track-E]   mode={mode}: starting...")
        p = propagate_mode(mode, sun_snap, moon_snap, duration_days=W_FULL_D, dt_s=DT_S)
        propagations[mode] = p
        print(f"[track-E]   mode={mode}: n_cross={p['n_ascending_nodes']}")

    # Headline numbers
    headline = {}
    full = propagations["sun_moon_j2"]
    sun_moon = propagations["sun_moon"]
    sun_only = propagations["sun_only"]
    moon_only = propagations["moon_only"]
    # Linear-fit full-year slopes (deg/day)
    intercept_full, slope_full, _rms = linear_fit_drift(
        (np.array(full["t_cross_s"]) - full["t_cross_s"][0]) / 86400.0,
        np.array(full["om_cross_rad"]),
    )
    intercept_sm, slope_sm, _rms = linear_fit_drift(
        (np.array(sun_moon["t_cross_s"]) - sun_moon["t_cross_s"][0]) / 86400.0,
        np.array(sun_moon["om_cross_rad"]),
    )
    intercept_s, slope_s, _rms = linear_fit_drift(
        (np.array(sun_only["t_cross_s"]) - sun_only["t_cross_s"][0]) / 86400.0,
        np.array(sun_only["om_cross_rad"]),
    )
    intercept_m, slope_m, _rms = linear_fit_drift(
        (np.array(moon_only["t_cross_s"]) - moon_only["t_cross_s"][0]) / 86400.0,
        np.array(moon_only["om_cross_rad"]),
    )
    # J2-only estimate: full_model - (sun+moon no-j2) = J2 contribution
    # (this is the same approach 017 used)
    j2_only_est_rad = slope_full - slope_sm  # deg/day conversion happens at output
    print(f"[track-E] HEADLINE (full-year linear fit, deg/day):")
    print(f"[track-E]   sun_only     : {math.degrees(slope_s):+.6e}")
    print(f"[track-E]   moon_only    : {math.degrees(slope_m):+.6e}")
    print(f"[track-E]   sun+moon     : {math.degrees(slope_sm):+.6e}")
    print(f"[track-E]   full model   : {math.degrees(slope_full):+.6e}")
    print(f"[track-E]   J2-only est  : {math.degrees(j2_only_est_rad):+.6e}  (= full - sun+moon)")
    print(f"[track-E]   Lunisolar    : {math.degrees(slope_full - j2_only_est_rad):+.6e}  (= full - J2)")

    headline = {
        "sun_only_full_year_deg_per_day": math.degrees(slope_s),
        "moon_only_full_year_deg_per_day": math.degrees(slope_m),
        "sun_moon_no_j2_full_year_deg_per_day": math.degrees(slope_sm),
        "full_model_full_year_deg_per_day": math.degrees(slope_full),
        "j2_only_estimate_full_year_deg_per_day": math.degrees(j2_only_est_rad),
        "lunisolar_only_full_year_deg_per_day": math.degrees(slope_full - j2_only_est_rad),
        "expected_published_full_model_deg_per_day": 0.9933,
        "expected_published_j2_only_deg_per_day": 0.9920,
        "expected_published_lunisolar_deg_per_day": 0.00132,
        "frame_choice": "ICRF/J2000 (no precession rotation applied)",
    }

    # ------------------------------------------------------------------------ #
    # Apply all estimators to each mode
    # ------------------------------------------------------------------------ #
    print("[track-E] applying 9 estimators to each of 4 force modes...")
    estimators_by_mode = {}
    for mode, prop in propagations.items():
        estimators_by_mode[mode] = apply_all_estimators(prop)
        full_year_slope = estimators_by_mode[mode]["full_year"]["slope_deg_per_day"]
        cycle_mean = estimators_by_mode[mode]["cycle_averaged"]["mean_slope_deg_per_day"]
        print(f"[track-E]   {mode}: full_year = {full_year_slope:+.6e} deg/day, "
              f"cycle_mean = {cycle_mean:+.6e} deg/day")

    # ------------------------------------------------------------------------ #
    # Residual summary (window-shift bias, phase bias, cycle bias)
    # ------------------------------------------------------------------------ #
    print("[track-E] computing residual summary...")
    residuals = residual_summary(estimators_by_mode)

    # Save JSON payload
    payload = {
        "constants": {
            "h_km": ALT_H_KM,
            "dt_s": DT_S,
            "duration_days": W_FULL_D,
            "R_E_km": R_EARTH_KM,
            "J2": J2_EARTH,
            "mu_E_km3_s2": MU_EARTH_KM3S2,
            "mu_Sun_km3_s2": SOLAR_GM_KM3_S2,
            "mu_Moon_km3_s2": LUNAR_GM_KM3_S2,
            "AU_km": AU_KM,
            "T0_s": T0_S,
            "anomalistic_month_d": ANOMALISTIC_MONTH_D,
            "synodic_month_d": SYNODIC_MONTH_D,
            "draconic_month_d": DRACONIC_MONTH_D,
            "W_full_d": W_FULL_D,
            "W_half_d": W_HALF_D,
            "W_phase_d": W_PHASE_D,
            "N_cycle_segments": N_CYCLE_SEGMENTS,
        },
        "snapshots": {
            "sun_sha256": sun_snap["sha256"],
            "sun_n_points": sun_snap["n_points"],
            "moon_sha256": moon_snap["sha256"],
            "moon_n_points": moon_snap["n_points"],
        },
        "frame_choice": (
            "ICRF/J2000: Sun and Moon snapshots used as-is, no precession "
            "rotation. Track E methodological choice; the 018 'with_precession' "
            "fix adds ~0.4 deg/year bias reduction. Documented in the report."
        ),
        "headline_full_year": headline,
        "propagation_summary": {
            mode: {
                "n_ascending_nodes": prop["n_ascending_nodes"],
                "i_sso_deg": prop["i_sso_deg"],
                "duration_days": prop["duration_days"],
                "dt_s": prop["dt_s"],
            } for mode, prop in propagations.items()
        },
        "estimators_by_mode": estimators_by_mode,
        "residual_summary": residuals,
    }
    OUT_RESULTS_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[track-E] wrote {OUT_RESULTS_JSON}")

    # Print the headline table to console for verification
    print("\n[track-E] VERIFICATION TABLE (headline slopes, deg/day):")
    print(f"  {'estimator':<24}{'sun_only':>16}{'moon_only':>16}{'sun_moon':>16}{'full':>16}")
    for est_name in [
        "full_year", "first_half", "second_half",
        "full_moon_start", "new_moon_start", "perigee_start",
        "W_anomalistic_month", "W_synodic_month",
    ]:
        row = f"  {est_name:<24}"
        for mode in modes:
            v = estimators_by_mode[mode][est_name]["slope_deg_per_day"]
            row += f"{v:>16.4e}"
        print(row)
    print()
    for mode in modes:
        cm = estimators_by_mode[mode]["cycle_averaged"]["mean_slope_deg_per_day"]
        cs = estimators_by_mode[mode]["cycle_averaged"]["std_slope_deg_per_day"]
        print(f"  cycle_mean_{mode:<15} = {cm:+.4e} deg/day  (cycle_std = {cs:.4e})")


if __name__ == "__main__":
    main()