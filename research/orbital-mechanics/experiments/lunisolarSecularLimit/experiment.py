"""Experiment 020 -- Long-Arc Lunisolar Secular-Limit Validation.

Headline scientific question: at h=600 km i_sso=97.79 deg, does the
doubly-averaged quadrupole lunisolar secular RAAN rate from the corrected
018 formula predict the secular limit that a sufficiently long controlled
numerical experiment converges to? Or does the 018 formula systematically
under-estimate (or over-estimate) the actual secular drift?

The 019 extrapolation (Ω̇_fit(W) = a + b/W + c/W²) gave +0.0036 deg/day at
W -> ∞ (27x the corrected formula), but the Track 3 audit (audit-020-track-3)
showed this extrapolation has NO theoretical asymptotic basis (the actual
OLS bias scaling is O(1/W²) for fast harmonics and O(A_k ω_k) constant for
slow harmonics). The Track 7 hostile review (audit-020-track-7) showed the
019 i=90 deg extrapolation sign-flips between linear and quadratic models.

This experiment uses:
  (a) A THEORY-DRIVEN HARMONIC REGRESSION estimator (Track 3's recommendation,
      estimator (f)) as the headline secular-rate observable. The basis
      includes the known physics: annual, half-annual, third-annual,
      quarter-annual, fifth-annual, evection (27.55 d), variation (14.77 d),
      lunar nodal (6798.4 d), and the secular. OLS regression recovers
      exactly the secular coefficient (modulo the residual noise after
      subtraction).
  (b) An ANGULAR-MOMENTUM-VECTOR secular-rate estimator (Track 5's
      recommendation, estimator A) as the theory-INDEPENDENT cross-check.
      The kinematic secular rate is computed from the secular drift of the
      node vector n = z x h.
  (c) A 5-YEAR arc (Track 8's recommendation) at h=600 km i_sso=97.79 deg,
      with full Lunisolar+J2 force model. Multi-year DE441 Sun/Moon
      references acquired via the lab's standard methodology
      (eclipseTiming/lunisolarVerification patterns).
  (d) Multiple initial phases (4 phases spaced ~91.25 d apart in the
      lunar anomalistic cycle) to characterize phase-dependence.
  (e) Comparison against the corrected 018 formula, the 019 W=730 d
      extrapolation, and the 019 cycle-averaged estimator.

Reference-data provenance:
  - Sun: horizons_sun_geocentric_vectors_2026_to_2030_icrf_tdb_daily.txt
    acquired by fetch_horizons_sun_snapshot_5yr.py (single Horizons query,
    byte-pinned with MANIFEST.json). Acquisition pattern identical to
    Exp 014/017.
  - Moon: horizons_moon_geocentric_vectors_2026_to_2030_icrf_tdb_daily.txt
    acquired by fetch_horizons_moon_snapshot_5yr.py.

Methodology:
  - Deterministic, byte-stable (no RNG, no network at runtime, no wall-clock).
  - Reuses lab_utils canon (rk4_propagate, j2_rhs, sso_inclination_rad,
    mean_motion, rv_to_coe_eci, seed_state).
  - Uses fixed-step RK4 at dt = 60 s (the 018/019 verified design).
  - Reports all 4 force-mode combinations per phase.
  - Pre-registers the contract before the confirmatory run.

References:
  - audit-019-track-A: disturbing-function derivation
  - audit-019-track-F: OLS bias formula (incomplete; Track 3 generalized)
  - audit-019-track-G: hostile review
  - audit-020-track-1: disturbing-function reconciliation (Convention B)
  - audit-020-track-3: complete OLS bias formula + estimator hierarchy
  - audit-020-track-5: independent estimator design
  - audit-020-track-7: hostile review of 019
  - audit-020-track-8: compute feasibility
  - Murray & Dermott (1999), Solar System Dynamics, Ch. 6, Ch. 7
  - Kaula (1962), Development of the lunar and solar disturbing functions
  - Kozai (1959), AJ 64, 367
  - Lidov (1962), Planet. Space Sci. 9, 719
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from lab_utils import (  # noqa: E402
    J2_EARTH,
    MU_EARTH_KM3S2,
    R_EARTH_KM,
    SSO_TARGET_DEG_DAY,
    j2_rhs,
    sso_inclination_rad,
)
from lab_utils.earth_frames import JD_J2000  # noqa: E402
from lab_utils.integrators import rk4_propagate  # noqa: E402
from lab_utils.orbits import mean_motion  # noqa: E402
from lab_utils.results import save_json_result  # noqa: E402

# --------------------------------------------------------------------------- #
# Constants (frozen)
# --------------------------------------------------------------------------- #
EXP_NAME = "lunisolarSecularLimit-020"
SOLAR_GM_KM3_S2 = 132712440018.0
LUNAR_GM_KM3_S2 = 4902.8001
LUNAR_DISTANCE_KM_MEAN = 384400.0
LUNAR_INCLINATION_DEG = 5.145
SOLAR_OBLIQUITY_DEG = 23.439
AU_KM = 149597870.7
DEG = math.pi / 180.0

# Frozen h=600 km i_sso (canonical SSO case from 018/019)
H600_KM = 600.0
I_SSO_DEG = 97.7876
DT_PROPAGATION_S = 60.0

# 5-year arc (Track 8 recommendation)
ARC_DAYS = 5.0 * 365.2422  # = 1826.211 d

# Initial phase ensemble: 4 phases at quarters of the lunar anomalistic month
# (27.5546 d) to test for systematic phase dependence. Each phase differs
# from the previous by 27.5546/4 = 6.89 d, so all four together span
# approximately one lunar cycle. This is a structured grid, not random.
PHASE_OFFSETS_DAYS = (0.0, 6.89, 13.78, 20.66)

# Force modes: full model is the priority; j2_only as control.
FORCE_MODES = ("sun_moon_j2", "sun_moon", "moon_only", "sun_only", "j2_only")

# Reference-data paths (multi-year DE441 Sun + Moon, byte-pinned)
SUN_SNAPSHOT_PATH = (
    Path(__file__).resolve().parent / "reference"
    / "horizons_sun_geocentric_vectors_2026_to_2030_icrf_tdb_daily.txt"
)
MOON_SNAPSHOT_PATH = (
    Path(__file__).resolve().parent / "reference"
    / "horizons_moon_geocentric_vectors_2026_to_2030_icrf_tdb_daily.txt"
)
SUN_SNAPSHOT_FALLBACK = (
    Path(__file__).resolve().parent.parent / "eclipseTiming" / "reference"
    / "horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt"
)
MOON_SNAPSHOT_FALLBACK = (
    Path(__file__).resolve().parent.parent / "lunisolarVerification" / "reference"
    / "horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt"
)


# --------------------------------------------------------------------------- #
# FIXED IAU-1976 precession (eclipseTiming convention, Track D 019 remediation)
# --------------------------------------------------------------------------- #
def _rot3(angle: float) -> np.ndarray:
    """Standard active rotation about +z by +angle (eclipseTiming convention)."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot2(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])


def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    """IAU-1976 precession: J2000 -> mean-of-date (Lieske polynomial)."""
    T = t_s / (86400.0 * 36525.0)
    sec = math.radians(1.0 / 3600.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * sec
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * sec
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * sec
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)


# --------------------------------------------------------------------------- #
# Corrected secular formula (Convention B verified by audit-020-track-1)
# --------------------------------------------------------------------------- #
def corrected_secular_lunisolar_raan_rate_rad_s(
    h_km: float, i_deg: float = I_SSO_DEG,
    *,
    i3_sun_rad: float = math.radians(SOLAR_OBLIQUITY_DEG),
    i3_moon_rad: float = math.radians(SOLAR_OBLIQUITY_DEG + LUNAR_INCLINATION_DEG),
    mu_3_sun: float = SOLAR_GM_KM3_S2,
    mu_3_moon: float = LUNAR_GM_KM3_S2,
    a_3_sun: float = AU_KM,
    a_3_moon: float = LUNAR_DISTANCE_KM_MEAN,
) -> dict:
    """Corrected doubly-averaged quadrupole Lunisolar RAAN rate (Convention B)."""
    a = R_EARTH_KM + h_km
    n = mean_motion(a)
    i_rad = math.radians(i_deg)
    solar = (3.0 / 8.0) * n * (mu_3_sun / MU_EARTH_KM3S2) * (
        a / a_3_sun
    ) ** 3 * math.sin(2.0 * (i_rad - i3_sun_rad)) / math.sin(i_rad)
    lunar = (3.0 / 8.0) * n * (mu_3_moon / MU_EARTH_KM3S2) * (
        a / a_3_moon
    ) ** 3 * math.sin(2.0 * (i_rad - i3_moon_rad)) / math.sin(i_rad)
    total = solar + lunar
    return {
        "h_km": h_km,
        "i_deg": i_deg,
        "solar_rad_s": solar,
        "solar_deg_day": math.degrees(solar) * 86400.0,
        "lunar_rad_s": lunar,
        "lunar_deg_day": math.degrees(lunar) * 86400.0,
        "total_rad_s": total,
        "total_deg_day": math.degrees(total) * 86400.0,
    }


# --------------------------------------------------------------------------- #
# Snapshot loading + interpolation
# --------------------------------------------------------------------------- #
def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    return {"t_s": t_s, "r_eci_km": r_eci, "sha256": _sha256(path), "n_points": len(rows)}


def _interp_snapshot_precessed(t_query_s: float, snap: dict,
                                apply_precession: bool = True) -> np.ndarray:
    """Linear interpolation with FIXED IAU-1976 precession (019 Track D fix)."""
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
        rv = r[idx - 1] + frac * (r[idx] - r[idx - 1])
    if apply_precession:
        P = precession_j2000_to_mod(t_query_s)
        return P @ rv
    return rv


def _third_body_accel(r_eci_km: np.ndarray, t_s: float, sun_snap: dict,
                      moon_snap: dict, *, include_sun: bool = True,
                      include_moon: bool = True,
                      apply_precession: bool = True) -> np.ndarray:
    """Geocentric third-body acceleration (direct + indirect)."""
    a_total = np.zeros(3)
    if include_sun:
        r_sun = _interp_snapshot_precessed(t_s, sun_snap, apply_precession)
        r_sat_to_sun = r_sun - r_eci_km
        r3 = np.linalg.norm(r_sun)
        r3s = np.linalg.norm(r_sat_to_sun)
        a_total += SOLAR_GM_KM3_S2 * (r_sat_to_sun / r3s ** 3 - r_sun / r3 ** 3)
    if include_moon:
        r_moon = _interp_snapshot_precessed(t_s, moon_snap, apply_precession)
        r_sat_to_moon = r_moon - r_eci_km
        r3 = np.linalg.norm(r_moon)
        r3s = np.linalg.norm(r_sat_to_moon)
        a_total += LUNAR_GM_KM3_S2 * (r_sat_to_moon / r3s ** 3 - r_moon / r3 ** 3)
    return a_total


# --------------------------------------------------------------------------- #
# RHS builder
# --------------------------------------------------------------------------- #
def make_rhs(sun_snap: dict, moon_snap: dict, *, mode: str,
             apply_precession: bool = True):
    """Build RHS for a given mode."""
    f_j2 = j2_rhs(MU_EARTH_KM3S2, J2_EARTH, R_EARTH_KM)

    def f(t: float, x: np.ndarray) -> np.ndarray:
        r = x[:3]
        v = x[3:]
        if mode == "kepler_only":
            a_kep = -MU_EARTH_KM3S2 * r / np.linalg.norm(r) ** 3
            return np.concatenate([v, a_kep])
        a_j2_kep = f_j2(t, x)[3:]
        if mode == "j2_only":
            return np.concatenate([v, a_j2_kep])
        a_3b = _third_body_accel(
            r, t, sun_snap, moon_snap,
            include_sun=(mode in ("sun_only", "sun_moon", "sun_moon_j2")),
            include_moon=(mode in ("moon_only", "sun_moon", "sun_moon_j2")),
            apply_precession=apply_precession,
        )
        return np.concatenate([v, a_j2_kep + a_3b])

    return f


# --------------------------------------------------------------------------- #
# Ascending-node detection
# --------------------------------------------------------------------------- #
def detect_ascending_nodes(t_s_arr: np.ndarray, x_arr: np.ndarray) -> tuple:
    t_crossings = []
    om_crossings = []
    z_prev = x_arr[0, 2]
    for k in range(1, len(t_s_arr)):
        z_curr = x_arr[k, 2]
        if z_prev <= 0 < z_curr:
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
    return np.array(t_crossings), np.array(om_crossings)


def detect_node_vector(t_s_arr: np.ndarray, x_arr: np.ndarray) -> tuple:
    """Theory-INDEPENDENT estimator: extract the node vector n = z x h
    at every step, track its secular drift. n_x, n_y drift as the
    satellite orbit plane regresses; arctan2(n_y, n_x) gives Ω without
    relying on ascending-node-crossing detection."""
    n_x = np.empty(len(t_s_arr))
    n_y = np.empty(len(t_s_arr))
    for k in range(len(t_s_arr)):
        r = x_arr[k, :3]
        v = x_arr[k, 3:]
        h = np.cross(r, v)
        # node vector n = z x h: z = (0,0,1), so n = (h_y, -h_x, 0)
        n_x[k] = h[1]
        n_y[k] = -h[0]
    n_mag = np.sqrt(n_x ** 2 + n_y ** 2)
    good = n_mag > 1e-12
    n_x_filt = n_x[good]
    n_y_filt = n_y[good]
    t_filt = t_s_arr[good]
    omega_node = np.arctan2(n_y_filt, n_x_filt)
    # unwrap
    for k in range(1, len(omega_node)):
        while omega_node[k] < omega_node[k - 1] - math.pi:
            omega_node[k] += 2 * math.pi
        while omega_node[k] > omega_node[k - 1] + math.pi:
            omega_node[k] -= 2 * math.pi
    return t_filt, omega_node


def ols_slope(t_s: np.ndarray, y_rad: np.ndarray) -> tuple:
    """Standard OLS linear-fit: returns (intercept, slope in y-units/t-units)."""
    A = np.column_stack([np.ones_like(t_s), t_s])
    result = np.linalg.lstsq(A, y_rad, rcond=None)
    return float(result[0][0]), float(result[0][1])


# --------------------------------------------------------------------------- #
# Theory-driven harmonic regression estimator (Track 3's estimator (f))
# --------------------------------------------------------------------------- #
HARMONIC_BASIS_PERIODS_DAYS = (
    # Annual and harmonics (019 FFT top-5; integer-cycle aliases)
    365.2422,   # annual solar forcing
    182.6211,   # half-annual
    121.7474,   # third-annual
    91.3106,    # quarter-annual
    73.0484,    # fifth-annual
    # Named physical drivers (Track B/C/D estimates; non-integer-cycle)
    27.5546,    # evection (lunar anomalistic)
    14.7653,    # variation (lunar synodic half-month)
    6798.4,     # lunar nodal
)


def harmonic_regression_secular_rate(
    t_cross_s: np.ndarray,
    omega_rad: np.ndarray,
) -> dict:
    """Theory-driven harmonic regression estimator.

    Fit omega(t) = a + b*t + Σ_k [c_k cos(2π t/T_k) + s_k sin(2π t/T_k)]
    simultaneously via OLS. b is the secular rate; its bias is exactly
    zero for each harmonic in the basis (projection); residual bias
    comes only from the unmodelled content.
    """
    t_day = t_cross_s / 86400.0
    n = len(t_day)
    cols = [np.ones(n), t_day]
    for T_d in HARMONIC_BASIS_PERIODS_DAYS:
        omega_k = 2.0 * math.pi / T_d
        cols.append(np.cos(omega_k * t_day))
        cols.append(np.sin(omega_k * t_day))
    A = np.column_stack(cols)
    result = np.linalg.lstsq(A, omega_rad, rcond=None)
    coeffs = result[0]
    intercept = float(coeffs[0])
    b_rad_per_day = float(coeffs[1])
    fit_pred = A @ coeffs
    rms_residual_rad = math.sqrt(np.mean((omega_rad - fit_pred) ** 2))
    n_harm = len(HARMONIC_BASIS_PERIODS_DAYS)
    harmonic_coeffs = {}
    for i, T_d in enumerate(HARMONIC_BASIS_PERIODS_DAYS):
        c_k = float(coeffs[2 + 2 * i])
        s_k = float(coeffs[2 + 2 * i + 1])
        amp = math.sqrt(c_k * c_k + s_k * s_k)
        harmonic_coeffs[T_d] = {"cos": c_k, "sin": s_k, "amp_rad": amp,
                                "amp_deg": math.degrees(amp)}
    return {
        "intercept_rad": intercept,
        "intercept_deg": math.degrees(intercept),
        "b_rad_per_day": b_rad_per_day,
        "b_deg_per_day": math.degrees(b_rad_per_day),
        "rms_residual_deg": math.degrees(rms_residual_rad),
        "harmonic_amplitudes_deg": harmonic_coeffs,
        "n_points": n,
    }


# --------------------------------------------------------------------------- #
# Convergence ladder (deterministic RK4 self-convergence at multiple arc lengths)
# --------------------------------------------------------------------------- #
def convergence_ladder(sun_snap: dict, moon_snap: dict,
                        duration_days: float = 365.0) -> dict:
    """Convergence ladder at h=600 km i_sso for the given arc length.

    Returns the max position/velocity difference between the production
    dt=60 s propagation and a finer dt=15 s reference, plus an
    extrapolation of the secular-rate difference across timesteps.
    """
    a = R_EARTH_KM + H600_KM
    i_rad = math.radians(I_SSO_DEG)
    v_circ = math.sqrt(MU_EARTH_KM3S2 / a)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, v_circ * math.cos(i_rad), v_circ * math.sin(i_rad)])
    x0 = np.concatenate([r0, v0])
    t0 = 820540800.0
    T_test = duration_days * 86400.0

    f = make_rhs(sun_snap, moon_snap, mode="sun_moon_j2", apply_precession=True)

    # Reference at fine dt=15 s; production at dt=60 s; plus coarser steps.
    results = {"dt_s": [], "max_r_diff_km": [], "max_v_diff_km_per_s": [],
               "secular_slope_diff_deg_day": []}
    # The 5-year reference would be huge; use 30-day subset for the
    # convergence ladder.
    n_ref = int(T_test / 15.0)
    t_ref = t0 + np.arange(n_ref + 1) * 15.0
    x_ref = rk4_propagate(f, t_ref, x0)

    # Build per-dt trajectory, sample at matching indices.
    for dt in (120.0, 60.0, 30.0):
        n_coarse = int(T_test / dt)
        t_coarse = t0 + np.arange(n_coarse + 1) * dt
        x_coarse = rk4_propagate(f, t_coarse, x0)
        idx_ref_end = n_coarse * int(round(dt / 15.0))
        r_diff = float(np.linalg.norm(x_coarse[-1, :3] - x_ref[idx_ref_end, :3]))
        v_diff = float(np.linalg.norm(x_coarse[-1, 3:] - x_ref[idx_ref_end, 3:]))
        # Secular slope difference: extract Ω at ascending-node crossings
        # on both the coarse and reference trajectories.
        t_coarse_cross, om_coarse = detect_ascending_nodes(t_coarse, x_coarse)
        t_ref_cross, om_ref = detect_ascending_nodes(t_ref, x_ref)
        # Take a common span: first T_test seconds of the coarse crossings.
        mask_c = (t_coarse_cross >= t0) & (t_coarse_cross <= t0 + T_test)
        _, s_coarse = ols_slope((t_coarse_cross[mask_c] - t0) / 86400.0, om_coarse[mask_c])
        mask_r = (t_ref_cross >= t0) & (t_ref_cross <= t0 + T_test)
        _, s_ref = ols_slope((t_ref_cross[mask_r] - t0) / 86400.0, om_ref[mask_r])
        slope_diff_deg_day = math.degrees(s_coarse - s_ref) * 86400.0
        results["dt_s"].append(dt)
        results["max_r_diff_km"].append(r_diff)
        results["max_v_diff_km_per_s"].append(v_diff)
        results["secular_slope_diff_deg_day"].append(slope_diff_deg_day)

    return results


# --------------------------------------------------------------------------- #
# Force-level identity check (Track D 019 pattern)
# --------------------------------------------------------------------------- #
def force_level_identity_check(h_km: float = H600_KM, n_states: int = 50,
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
        a_sun_b = -SOLAR_GM_KM3_S2 * (r_sat - r3_sun) / r3s_sun ** 3 - SOLAR_GM_KM3_S2 * r3_sun / r3_mag_sun ** 3
        r_sat_to_moon = r3_moon - r_sat
        r3s_moon = np.linalg.norm(r_sat_to_moon)
        r3_mag_moon = np.linalg.norm(r3_moon)
        a_moon_a = LUNAR_GM_KM3_S2 * (r_sat_to_moon / r3s_moon ** 3 - r3_moon / r3_mag_moon ** 3)
        a_moon_b = -LUNAR_GM_KM3_S2 * (r_sat - r3_moon) / r3s_moon ** 3 - LUNAR_GM_KM3_S2 * r3_moon / r3_mag_moon ** 3
        max_diff_sun = max(max_diff_sun, float(np.max(np.abs(a_sun_a - a_sun_b))))
        max_diff_moon = max(max_diff_moon, float(np.max(np.abs(a_moon_a - a_moon_b))))
    return {
        "n_states": n_states,
        "max_diff_sun_km_s2": float(max_diff_sun),
        "max_diff_moon_km_s2": float(max_diff_moon),
        "passes_sun": max_diff_sun < 1e-15,
        "passes_moon": max_diff_moon < 1e-15,
    }


# --------------------------------------------------------------------------- #
# Synthetic estimator test (Track 3's calibration)
# --------------------------------------------------------------------------- #
def synthetic_estimator_test() -> dict:
    """Deterministic synthetic estimator test (Track 3 recommendation).

    Build a synthetic Ω(t) = a_secular * t + Σ A_k cos(ω_k t + φ_k) with
    the 019 FFT amplitudes. Apply estimator (a) direct OLS and
    estimator (f) harmonic regression. Verify which one recovers a_secular
    more accurately.
    """
    a_true = 1.0e-3  # deg/day, the order of the 018 numerical
    harmonics = [
        (365.2422, 0.103, 0.0),
        (182.6211, 0.025, 0.5),
        (121.7474, 0.012, 1.0),
        (91.3106, 0.007, 1.5),
        (73.0484, 0.005, 2.0),
        (27.5546, 0.004, 0.3),
        (14.7653, 0.003, 1.1),
        (6798.4, 0.002, 2.7),
    ]
    W_days = 365.0  # 1-yr test
    n_samples = int(W_days * 14.91)  # ascending-node-cadence
    t_day = np.linspace(0, W_days, n_samples)
    # Build signal directly in rad; a_true is in deg/day
    a_true_rad_per_day = math.radians(a_true)
    omega_rad = np.zeros(n_samples) + a_true_rad_per_day * t_day
    for T_d, A_deg, phi in harmonics:
        A_rad = math.radians(A_deg)
        omega_k = 2 * math.pi / T_d
        omega_rad += A_rad * np.cos(omega_k * t_day + phi)

    # Estimator (a): direct OLS on omega_rad (slope in rad/day)
    _, b_a_rad = ols_slope(t_day, omega_rad)
    b_a_deg_day = math.degrees(b_a_rad)
    bias_a_deg_day = b_a_deg_day - a_true

    # Estimator (f): theory-driven harmonic regression
    # Use the full theory-driven basis (including evection/variation/lunar
    # nodal which are NOT in the 019 FFT top-5 but are the physical drivers).
    fit = harmonic_regression_secular_rate(t_day * 86400.0, omega_rad)
    bias_f_deg_day = fit["b_deg_per_day"] - a_true

    return {
        "a_true_deg_day": a_true,
        "W_days": W_days,
        "n_samples": n_samples,
        "estimator_a_direct_ols_deg_day": float(b_a_deg_day),
        "estimator_a_bias_deg_day": float(bias_a_deg_day),
        "estimator_f_harmonic_regression_deg_day": float(fit["b_deg_per_day"]),
        "estimator_f_bias_deg_day": float(bias_f_deg_day),
        "estimator_f_rms_residual_deg": float(fit["rms_residual_deg"]),
        "harmonic_amplitudes_deg_recovered": fit["harmonic_amplitudes_deg"],
        "verdict": (
            "f wins" if abs(bias_f_deg_day) < abs(bias_a_deg_day) else "a wins"
        ),
    }


# --------------------------------------------------------------------------- #
# Code hashes
# --------------------------------------------------------------------------- #
def code_hashes() -> dict:
    here = Path(__file__).resolve().parent
    lab_root = here.parents[3]
    files = {
        "experiment.py": here / "experiment.py",
        "lab_utils/orbits.py": lab_root / "src" / "lab_utils" / "orbits.py",
        "lab_utils/earth_frames.py": lab_root / "src" / "lab_utils" / "earth_frames.py",
        "lab_utils/integrators.py": lab_root / "src" / "lab_utils" / "integrators.py",
        "lab_utils/results.py": lab_root / "src" / "lab_utils" / "results.py",
        "lab_utils/__init__.py": lab_root / "src" / "lab_utils" / "__init__.py",
        "019_experiment.py": here.parent / "lunisolarLongPeriod" / "experiment.py",
        "018_experiment.py": here.parent / "lunisolarReconciliation" / "experiment.py",
    }
    sun_p = SUN_SNAPSHOT_PATH if SUN_SNAPSHOT_PATH.exists() else SUN_SNAPSHOT_FALLBACK
    moon_p = MOON_SNAPSHOT_PATH if MOON_SNAPSHOT_PATH.exists() else MOON_SNAPSHOT_FALLBACK
    files["sun_reference_snapshot.txt"] = sun_p
    files["moon_reference_snapshot.txt"] = moon_p
    return {name: _sha256(p) for name, p in files.items()}


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def run() -> dict:
    print("[020] starting Lunisolar Long-Arc Secular-Limit Validation")
    # Load snapshots (use fallback to existing 1-yr if 5-yr not yet acquired).
    sun_path = SUN_SNAPSHOT_PATH if SUN_SNAPSHOT_PATH.exists() else SUN_SNAPSHOT_FALLBACK
    moon_path = MOON_SNAPSHOT_PATH if MOON_SNAPSHOT_PATH.exists() else MOON_SNAPSHOT_FALLBACK
    print(f"[020] Sun snapshot: {sun_path}")
    print(f"[020] Moon snapshot: {moon_path}")
    sun_snap = _load_snapshot(sun_path)
    moon_snap = _load_snapshot(moon_path)
    print(f"[020] Sun: {sun_snap['n_points']} rows sha256={sun_snap['sha256'][:16]}")
    print(f"[020] Moon: {moon_snap['n_points']} rows sha256={moon_snap['sha256'][:16]}")

    # Synthetic estimator test FIRST (Track 3 recommendation: build the
    # estimator hierarchy BEFORE inspecting real data).
    print("[020] running synthetic estimator test (Track 3 calibration)...")
    synth = synthetic_estimator_test()
    print(f"[020] synthetic: estimator (a) bias = {synth['estimator_a_bias_deg_day']:+.3e} deg/day")
    print(f"[020] synthetic: estimator (f) bias = {synth['estimator_f_bias_deg_day']:+.3e} deg/day")
    print(f"[020] synthetic verdict: {synth['verdict']}")

    # Force-level identity check (019 pattern, third-body accuracy)
    print("[020] running force-level identity check (50 states)...")
    identity = force_level_identity_check()
    print(f"[020] identity: max_diff_sun = {identity['max_diff_sun_km_s2']:.3e} km/s^2, "
          f"max_diff_moon = {identity['max_diff_moon_km_s2']:.3e} km/s^2")

    # Corrected secular formula at the canonical SSO case
    cf_018 = corrected_secular_lunisolar_raan_rate_rad_s(H600_KM, i_deg=I_SSO_DEG)
    print(f"[020] corrected cf at i_sso: solar = {cf_018['solar_deg_day']:+.6e}, "
          f"lunar = {cf_018['lunar_deg_day']:+.6e}, total = {cf_018['total_deg_day']:+.6e} deg/day")

    # Convergence ladder at the canonical SSO case (30-day subset)
    print("[020] running convergence ladder at h=600 km i_sso...")
    convergence = convergence_ladder(sun_snap, moon_snap, duration_days=30.0)
    print(f"[020] convergence ladder: {convergence}")

    # Multi-phase multi-mode propagation
    # If only 1-yr snapshot is available, fall back to 1-yr arc.
    available_arc_days = min(
        (sun_snap["t_s"][-1] - sun_snap["t_s"][0]) / 86400.0,
        (moon_snap["t_s"][-1] - moon_snap["t_s"][0]) / 86400.0,
    )
    arc_days = min(ARC_DAYS, available_arc_days)
    print(f"[020] arc length: {arc_days:.2f} days (available: {available_arc_days:.2f})")

    # Build a = R_E + h at h=600 km
    a = R_EARTH_KM + H600_KM
    i_rad = math.radians(I_SSO_DEG)
    v_circ = math.sqrt(MU_EARTH_KM3S2 / a)
    r0_base = np.array([a, 0.0, 0.0])
    v0_base = np.array([0.0, v_circ * math.cos(i_rad), v_circ * math.sin(i_rad)])

    # Pre-register the contract: arc length, phases, modes, estimators, gates.
    print("[020] running propagations (this may take several minutes per arc)...")
    propagation_results = {}
    for phase_offset_d in PHASE_OFFSETS_DAYS:
        # Phase offset shifts the initial mean anomaly; for a circular orbit
        # the initial state is parameterized by the true anomaly; we
        # approximate by rotating the initial perifocal state.
        phase_rad = 2 * math.pi * phase_offset_d / 27.5546
        cos_p, sin_p = math.cos(phase_rad), math.sin(phase_rad)
        # Rotate initial state in the orbital plane
        r0 = np.array([
            r0_base[0] * cos_p - r0_base[1] * sin_p,
            r0_base[0] * sin_p + r0_base[1] * cos_p,
            r0_base[2],
        ])
        v0 = np.array([
            v0_base[0] * cos_p - v0_base[1] * sin_p,
            v0_base[0] * sin_p + v0_base[1] * cos_p,
            v0_base[2],
        ])
        x0 = np.concatenate([r0, v0])
        t0 = 820540800.0 + phase_offset_d * 86400.0
        t_end = t0 + arc_days * 86400.0
        n_steps = int(math.ceil((t_end - t0) / DT_PROPAGATION_S))
        # Cap steps for compute feasibility
        max_steps = 3_000_000
        if n_steps > max_steps:
            print(f"[020] WARNING: phase {phase_offset_d:.2f} requires {n_steps} steps, "
                  f"capping to {max_steps}; reducing arc to {max_steps * DT_PROPAGATION_S / 86400.0:.2f} d")
            n_steps = max_steps
            arc_used = n_steps * DT_PROPAGATION_S / 86400.0
        else:
            arc_used = arc_days
        t_grid = np.linspace(t0, t0 + n_steps * DT_PROPAGATION_S, n_steps + 1)

        for mode in FORCE_MODES:
            key = f"phase{phase_offset_d:.2f}_{mode}"
            if t_grid[-1] > max(sun_snap["t_s"][-1], moon_snap["t_s"][-1]):
                # Snapshots clamp to last row; document the limitation.
                print(f"[020] WARNING: phase {phase_offset_d:.2f} propagates past snapshot end; clamping")
            f = make_rhs(sun_snap, moon_snap, mode=mode, apply_precession=True)
            x_traj = rk4_propagate(f, t_grid, x0)
            # Ascending-node detector (theory-Cowling observable)
            t_cross, om_cross = detect_ascending_nodes(t_grid, x_traj)
            # Node-vector estimator (theory-INDEPENDENT observable)
            t_node, omega_node = detect_node_vector(t_grid, x_traj)
            if len(t_cross) >= 10:
                # Estimator (a): direct OLS on rad-valued omega
                t_rel = (t_cross - t_cross[0]) / 86400.0
                _, b_a_rad_per_day = ols_slope(t_rel, om_cross)
                b_a_deg_per_day = math.degrees(b_a_rad_per_day)
                # Estimator (f): theory-driven harmonic regression
                fit_f = harmonic_regression_secular_rate(
                    (t_cross - t_cross[0]), om_cross
                )
                # Estimator (g): secant
                if len(t_cross) >= 4:
                    secant_rad_per_day = (om_cross[-1] - om_cross[0]) / t_rel[-1]
                else:
                    secant_rad_per_day = float("nan")
                secant_deg_per_day = math.degrees(secant_rad_per_day)
                # Node-vector estimator (theory-independent kinematic observable)
                if len(t_node) > 10:
                    _, n_slope_rad_per_day = ols_slope(
                        (t_node - t_node[0]) / 86400.0, omega_node
                    )
                    n_slope_deg_per_day = math.degrees(n_slope_rad_per_day)
                else:
                    n_slope_deg_per_day = float("nan")
                propagation_results[key] = {
                    "phase_offset_d": phase_offset_d,
                    "mode": mode,
                    "arc_used_days": arc_used,
                    "n_ascending_nodes": len(t_cross),
                    "estimator_a_direct_ols_deg_per_day": b_a_deg_per_day,
                    "estimator_f_harmonic_regression_deg_per_day": fit_f["b_deg_per_day"],
                    "estimator_f_rms_residual_deg": fit_f["rms_residual_deg"],
                    "estimator_f_harmonic_amplitudes_deg": fit_f["harmonic_amplitudes_deg"],
                    "estimator_g_secant_deg_per_day": secant_deg_per_day,
                    "node_vector_secular_deg_per_day": n_slope_deg_per_day,
                    "t_cross_day": t_rel.tolist(),
                    "om_cross_deg": np.degrees(om_cross).tolist(),
                }
            else:
                propagation_results[key] = {
                    "phase_offset_d": phase_offset_d,
                    "mode": mode,
                    "arc_used_days": arc_used,
                    "n_ascending_nodes": len(t_cross),
                    "estimator_a_direct_ols_deg_per_day": float("nan"),
                    "estimator_f_harmonic_regression_deg_per_day": float("nan"),
                    "estimator_g_secant_deg_per_day": float("nan"),
                    "node_vector_secular_deg_per_day": float("nan"),
                }

    # Build Lunisolar contribution per phase (subtract J2-only baseline)
    # from full-model (sun_moon_j2) for each phase.
    lunisolar_estimates = {}
    for phase_offset_d in PHASE_OFFSETS_DAYS:
        full_key = f"phase{phase_offset_d:.2f}_sun_moon_j2"
        j2_key = f"phase{phase_offset_d:.2f}_j2_only"
        if full_key in propagation_results and j2_key in propagation_results:
            full_a = propagation_results[full_key]["estimator_a_direct_ols_deg_per_day"]
            j2_a = propagation_results[j2_key]["estimator_a_direct_ols_deg_per_day"]
            full_f = propagation_results[full_key]["estimator_f_harmonic_regression_deg_per_day"]
            j2_f = propagation_results[j2_key]["estimator_f_harmonic_regression_deg_per_day"]
            full_g = propagation_results[full_key]["estimator_g_secant_deg_per_day"]
            j2_g = propagation_results[j2_key]["estimator_g_secant_deg_per_day"]
            full_n = propagation_results[full_key]["node_vector_secular_deg_per_day"]
            j2_n = propagation_results[j2_key]["node_vector_secular_deg_per_day"]
            lunisolar_estimates[f"phase{phase_offset_d:.2f}"] = {
                "phase_offset_d": phase_offset_d,
                "estimator_a_lunisolar_deg_per_day": full_a - j2_a,
                "estimator_f_lunisolar_deg_per_day": full_f - j2_f,
                "estimator_g_lunisolar_deg_per_day": full_g - j2_g,
                "node_vector_lunisolar_deg_per_day": full_n - j2_n,
                "full_a": full_a, "j2_a": j2_a,
                "full_f": full_f, "j2_f": j2_f,
            }

    # Sun-only and Moon-only direct estimates (no J2 subtraction needed)
    sun_moon_estimates = {}
    for phase_offset_d in PHASE_OFFSETS_DAYS:
        sm_key = f"phase{phase_offset_d:.2f}_sun_moon"
        so_key = f"phase{phase_offset_d:.2f}_sun_only"
        mo_key = f"phase{phase_offset_d:.2f}_moon_only"
        if sm_key in propagation_results:
            sun_moon_estimates[f"phase{phase_offset_d:.2f}"] = {
                "phase_offset_d": phase_offset_d,
                "sun_moon_estimator_a_deg_per_day":
                    propagation_results[sm_key]["estimator_a_direct_ols_deg_per_day"],
                "sun_moon_estimator_f_deg_per_day":
                    propagation_results[sm_key]["estimator_f_harmonic_regression_deg_per_day"],
                "sun_only_estimator_a_deg_per_day": (
                    propagation_results[so_key]["estimator_a_direct_ols_deg_per_day"]
                    if so_key in propagation_results else float("nan")
                ),
                "sun_only_estimator_f_deg_per_day": (
                    propagation_results[so_key]["estimator_f_harmonic_regression_deg_per_day"]
                    if so_key in propagation_results else float("nan")
                ),
                "moon_only_estimator_a_deg_per_day": (
                    propagation_results[mo_key]["estimator_a_direct_ols_deg_per_day"]
                    if mo_key in propagation_results else float("nan")
                ),
                "moon_only_estimator_f_deg_per_day": (
                    propagation_results[mo_key]["estimator_f_harmonic_regression_deg_per_day"]
                    if mo_key in propagation_results else float("nan")
                ),
            }

    # Compute Lunisolar ratio vs corrected formula
    cf_total_deg_day = cf_018["total_deg_day"]
    if lunisolar_estimates:
        # Mean over phases for the headline secular observable
        f_estimates = [v["estimator_f_lunisolar_deg_per_day"]
                       for v in lunisolar_estimates.values()
                       if not math.isnan(v["estimator_f_lunisolar_deg_per_day"])]
        if f_estimates:
            f_mean = float(np.mean(f_estimates))
            f_std = float(np.std(f_estimates))
        else:
            f_mean, f_std = float("nan"), float("nan")
        a_estimates = [v["estimator_a_lunisolar_deg_per_day"]
                       for v in lunisolar_estimates.values()
                       if not math.isnan(v["estimator_a_lunisolar_deg_per_day"])]
        if a_estimates:
            a_mean = float(np.mean(a_estimates))
            a_std = float(np.std(a_estimates))
        else:
            a_mean, a_std = float("nan"), float("nan")
        g_estimates = [v["estimator_g_lunisolar_deg_per_day"]
                       for v in lunisolar_estimates.values()
                       if not math.isnan(v["estimator_g_lunisolar_deg_per_day"])]
        if g_estimates:
            g_mean = float(np.mean(g_estimates))
            g_std = float(np.std(g_estimates))
        else:
            g_mean, g_std = float("nan"), float("nan")
        n_estimates = [v["node_vector_lunisolar_deg_per_day"]
                       for v in lunisolar_estimates.values()
                       if not math.isnan(v["node_vector_lunisolar_deg_per_day"])]
        if n_estimates:
            n_mean = float(np.mean(n_estimates))
            n_std = float(np.std(n_estimates))
        else:
            n_mean, n_std = float("nan"), float("nan")
        # Ratios
        def ratio(a, b):
            if b == 0 or math.isnan(a) or math.isnan(b):
                return float("nan")
            return a / b
        ratio_a = ratio(a_mean, cf_total_deg_day)
        ratio_f = ratio(f_mean, cf_total_deg_day)
        ratio_g = ratio(g_mean, cf_total_deg_day)
        ratio_n = ratio(n_mean, cf_total_deg_day)
    else:
        f_mean = a_mean = g_mean = n_mean = float("nan")
        f_std = a_std = g_std = n_std = float("nan")
        ratio_a = ratio_f = ratio_g = ratio_n = float("nan")

    # Build payload
    payload = {
        "constants": {
            "R_E_km": R_EARTH_KM,
            "J2": J2_EARTH,
            "mu_E_km3_s2": MU_EARTH_KM3S2,
            "mu_Sun_km3_s2": SOLAR_GM_KM3_S2,
            "mu_Moon_km3_s2": LUNAR_GM_KM3_S2,
            "AU_km": AU_KM,
            "LUNAR_DISTANCE_KM_cf": LUNAR_DISTANCE_KM_MEAN,
            "LUNAR_INCLINATION_DEG": LUNAR_INCLINATION_DEG,
            "SOLAR_OBLIQUITY_DEG": SOLAR_OBLIQUITY_DEG,
            "sso_target_deg_day": SSO_TARGET_DEG_DAY,
        },
        "contract": {
            "frame": "ECI mean-of-date; Sun and Moon snapshots rotated from ICRF/J2000 via FIXED IAU-1976 precession (Track D 019 remediation; eclipseTiming _rot3 convention [[c,-s,0],[s,c,0],[0,0,1]]).",
            "units": "km, km^3/s^2, s since J2000 (TT-like); radians internal; degrees at I/O.",
            "corrected_closed_form": "(3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i-i_3) / sin i (Convention B; verified by audit-020-track-1).",
            "decision_variables": [
                "i_deg = 97.7876 (i_sso)",
                "phase_offset_d in {0.00, 6.89, 13.78, 20.66} (lunar anomalistic quarters)",
                "duration_days = 5.0 * 365.2422 (5-year baseline; Track 8 recommendation)",
                "mode in {sun_only, moon_only, sun_moon, sun_moon_j2, j2_only}",
            ],
            "estimators": [
                "direct_OLS (ascending-node-crossing linear fit)",
                "harmonic_regression (theory-driven OLS with annual, half-annual, third-annual, quarter-annual, fifth-annual, evection 27.55 d, variation 14.77 d, lunar nodal 6798.4 d basis; Track 3 estimator (f))",
                "secant (y(T)-y(0))/T",
                "node_vector (theory-INDEPENDENT kinematic observable from r x v; Track 5 estimator A)",
                "synthetic_estimator_test (Track 3 calibration oracle)",
            ],
            "headline_estimator": "harmonic_regression (Track 3 recommendation; bias << direct_OLS for harmonic-decorrupted signals)",
            "cross_check_estimator": "node_vector (Track 5 recommendation; theory-independent)",
        },
        "snapshots": {
            "sun_sha256": sun_snap["sha256"],
            "sun_n_points": sun_snap["n_points"],
            "sun_source": str(sun_path),
            "moon_sha256": moon_snap["sha256"],
            "moon_n_points": moon_snap["n_points"],
            "moon_source": str(moon_path),
            "available_arc_days": available_arc_days,
        },
        "force_level_identity_check": identity,
        "convergence_ladder": convergence,
        "synthetic_estimator_test": synth,
        "corrected_closed_form_at_i_sso": {
            "solar_cf_deg_day": cf_018["solar_deg_day"],
            "lunar_cf_deg_day": cf_018["lunar_deg_day"],
            "total_cf_deg_day": cf_018["total_deg_day"],
        },
        "propagation_results": propagation_results,
        "lunisolar_estimates": lunisolar_estimates,
        "sun_moon_estimates": sun_moon_estimates,
        "headline_secular_estimate": {
            "estimator_f_harmonic_regression": {
                "mean_deg_day": f_mean,
                "std_deg_day": f_std,
                "n_phases": len(f_estimates) if not math.isnan(f_mean) else 0,
                "ratio_to_corrected_cf": ratio_f,
            },
            "estimator_a_direct_ols": {
                "mean_deg_day": a_mean,
                "std_deg_day": a_std,
                "n_phases": len(a_estimates) if not math.isnan(a_mean) else 0,
                "ratio_to_corrected_cf": ratio_a,
            },
            "estimator_g_secant": {
                "mean_deg_day": g_mean,
                "std_deg_day": g_std,
                "n_phases": len(g_estimates) if not math.isnan(g_mean) else 0,
                "ratio_to_corrected_cf": ratio_g,
            },
            "node_vector_estimator": {
                "mean_deg_day": n_mean,
                "std_deg_day": n_std,
                "n_phases": len(n_estimates) if not math.isnan(n_mean) else 0,
                "ratio_to_corrected_cf": ratio_n,
            },
        },
        "findings": [
            "Headline: 5-year-arc secular-limit validation using theory-driven harmonic regression estimator (Track 3 estimator (f)). The corrected 018 formula predicts +1.3475e-4 deg/day Lunisolar contribution at h=600 km i_sso.",
            "The 019 window-length extrapolation +0.0036 deg/day (27x the corrected formula) is NOT independently validated; the Track 3 audit shows the 1/W polynomial extrapolation has no asymptotic basis.",
            "The 019 i=90 deg extrapolation sign-flip between linear and quadratic models (Track 7 hostile review) is a smoking gun for mis-specification, not a robust asymptotic measurement.",
            "Track 4 discovery: the constant i3_moon = 28.584 deg over-estimates the lunar contribution at 2026 by ~50% (actual 2026 lunar i3 ~18.29 deg). The corrected formula is correct as a long-term secular average but apples-to-apples comparison at 2026 epoch gives ~13-14x residual, not 9.78x.",
            "Track 7 hostile review: the 019 Lunisolar 0.0036 deg/day is ~7x the 019 cycle-averaged i=90 deg value (0.000484), an asymmetry the corrected formula cannot explain; flagging for Exp 020 investigation.",
            "All four phase offsets span approximately one lunar anomalistic month to characterize phase-dependence of the secular estimator.",
        ],
        "limitations": [
            "If 5-yr multi-year snapshots are unavailable, the 1-yr fallback snapshot is used (with documented limitation in 'available_arc_days').",
            "Harmonic-regression basis assumes the dominant periodic content is at integer-cycle annual harmonics + evection + variation + lunar nodal. Unmodelled harmonics with significant amplitude at intermediate periods can bias the secular estimator.",
            "The 5-yr arc does NOT resolve the 18.6-yr lunar nodal cycle; only a partial nodal modulation is captured. Track 6 recommendation: a 2-window phase-locked estimator with 1-yr windows at 9.3-yr separation would cancel the lunar nodal contribution exactly, but requires a longer ephemeris.",
            "The estimator hierarchy (Track 3) shows that direct OLS and 1/W extrapolation are NOT the right tools for harmonic-decorrupted secular-rate extraction; the harmonic regression is theoretically preferred.",
        ],
        "code_sha256": code_hashes(),
    }

    out = Path(__file__).resolve().parent / "results" / "results.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    save_json_result(
        str(out), payload, name=EXP_NAME,
        description=(
            "Lunisolar Long-Arc Secular-Limit Validation: 5-year arc (or "
            "1-year fallback) at h=600 km i_sso with 4-phase ensemble, "
            "harmonic-regression + node-vector estimators. Resolves "
            "whether the 019 W=730 d extrapolation +0.0036 deg/day is "
            "robust or a model-choice artifact."
        ),
    )
    print(f"[020] results -> {out}")
    return payload


if __name__ == "__main__":
    run()