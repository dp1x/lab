"""Experiment for mission_j2_lunisolar_coupling.

Investigates whether the J2 x Lunisolar cross-coupling is a competitive
explanation for the 18.6-yr DE441 RAAN residual observed in
mission_lunisolar_closure.

Reuses machinery from mission_lunisolar_closure/experiment.py (streaming RK4,
IAU-1976 precession, ascending-node detection, estimators) and extends it
with:

- Force-mode isolation: 2-body, J2-only, Sun-only, Moon-only, Sun+Moon,
  J2+Sun+Moon (6 modes; used for force-model decomposition).
- Perturbative multipliers: lambda_J2 in [0, 2], lambda_3body in [0, 2]
  (used for the perturbative scaling experiment).
- Reduced model: synthetic circular Moon with fixed geometry (used for
  the idealized test of J2 x Lunisolar coupling in isolation).

The mission's central observable is the cross-coupling residual:

    R_J2x3b = Omega_dot_full - Omega_dot_J2 - Omega_dot_Sun - Omega_dot_Moon + Omega_dot_2body

computed by mode subtraction at multiple inclinations. The mission's
discriminator is a 2-D polynomial fit to f(lambda_J2, lambda_3body): a
significant coefficient on the lambda_J2 * lambda_3body cross term means a
genuine coupling signal.

Streaming RK4 propagator: no full-trajectory storage; ascending-node
crossings and subsampled node-vector samples only.

Deterministic: fixed seed, fixed inputs, byte-pinned DE441 snapshots.
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

EXP_NAME = "mission_j2_lunisolar_coupling_001"
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

# --------------------------------------------------------------------------- #
# Mission paths
# --------------------------------------------------------------------------- #
HERE = Path(__file__).resolve().parent
SUN_SNAPSHOT = HERE.parent / "mission_lunisolar_closure" / "reference" / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MOON_SNAPSHOT = HERE.parent / "mission_lunisolar_closure" / "reference" / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
PARENT_MANIFEST = HERE.parent / "mission_lunisolar_closure" / "reference" / "MANIFEST.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# IAU-1976 precession (same as mission_lunisolar_closure)
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
# Snapshot loading (reused from mission_lunisolar_closure)
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


# --------------------------------------------------------------------------- #
# Synthetic circular Moon geometry (reduced-model variant)
# --------------------------------------------------------------------------- #
def synthetic_circular_moon_state(t_s: float) -> np.ndarray:
    """Synthetic Moon in a circular orbit in its own plane with fixed
    inclination i3_moon = eps + I_moon = 28.584 deg to the equator, and
    RAAN = 0 (no nodal regression). Period = 27.3217 d (sidereal).
    Returns the geocentric Moon vector in km.

    No eccentricity, no inclination oscillation, no nodal regression:
    this isolates the doubly-averaged quadrupole + J2 cross coupling
    from all periodic lunar effects.
    """
    T_sidereal = 27.3217 * 86400.0
    M_moon = 2 * math.pi * (t_s / T_sidereal)
    r_ecl = np.array([
        LUNAR_DISTANCE_KM_MEAN * math.cos(M_moon),
        LUNAR_DISTANCE_KM_MEAN * math.sin(M_moon),
        0.0,
    ])
    # Rotate ecliptic -> equatorial by solar obliquity + lunar inclination
    # (i3_moon = eps + I_moon, around X axis)
    i3 = math.radians(SOLAR_OBLIQUITY_DEG + LUNAR_INCLINATION_DEG)
    c, s = math.cos(i3), math.sin(i3)
    R_ecl_to_eq = np.array([
        [1.0, 0.0, 0.0],
        [0.0, c, -s],
        [0.0, s, c],
    ])
    return R_ecl_to_eq @ r_ecl


# --------------------------------------------------------------------------- #
# Third-body acceleration with perturbative multiplier
# --------------------------------------------------------------------------- #
def _third_body_accel(r_eci_km: np.ndarray, t_s: float, sun_snap: dict,
                      moon_snap: dict, *,
                      lambda_3body: float = 1.0,
                      synthetic_moon: bool = False,
                      include_sun: bool = True,
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
        if synthetic_moon:
            r_moon = synthetic_circular_moon_state(t_s)
            if apply_precession:
                P = precession_j2000_to_mod(t_s)
                r_moon = P @ r_moon
        else:
            r_moon = _interp_snapshot_precessed(t_s, moon_snap, apply_precession)
        r_sat_to_moon = r_moon - r_eci_km
        r3s = np.linalg.norm(r_sat_to_moon)
        r3 = np.linalg.norm(r_moon)
        a_total += LUNAR_GM_KM3_S2 * (r_sat_to_moon / r3s ** 3 - r_moon / r3 ** 3)
    return lambda_3body * a_total


# --------------------------------------------------------------------------- #
# Streaming RK4 propagation with mode isolation + perturbative multipliers
# --------------------------------------------------------------------------- #
def propagate_streaming_with_x0(sun_snap, moon_snap, x0: np.ndarray, *,
                                  mode: str, t0_s: float, t_end_s: float,
                                  dt_s: float = DT_S,
                                  subsample_every: int = 100,
                                  lambda_j2: float = 1.0,
                                  lambda_3body: float = 1.0,
                                  synthetic_moon: bool = False) -> dict:
    """Stream RK4 step-by-step, collect ascending-node crossings +
    subsampled node-vector samples. NO full-trajectory storage.

    Modes:
      - kepler_only: 2-body only (no J2, no 3b)
      - j2_only: J2 + 2-body (no 3b)
      - sun_only: Sun + 2-body (no J2, no Moon)
      - moon_only: Moon + 2-body (no J2, no Sun)
      - sun_moon: Sun + Moon + 2-body (no J2)
      - sun_moon_j2: full model (the canonical mission setup)

    Perturbative multipliers:
      - lambda_j2: scales J2 acceleration (default 1.0)
      - lambda_3body: scales third-body (Sun + Moon) acceleration (default 1.0)
    """
    use_sun = mode in ("sun_only", "sun_moon", "sun_moon_j2")
    use_moon = mode in ("moon_only", "sun_moon", "sun_moon_j2")
    use_j2 = mode in ("j2_only", "sun_moon_j2")

    # J2 RHS with multiplier
    re2 = R_EARTH_KM ** 2

    def f(t, x):
        r = x[:3]
        v = x[3:]
        rm = np.linalg.norm(r)
        a_kep = -MU_EARTH_KM3S2 * r / rm ** 3
        a = a_kep.copy()
        if use_j2 and lambda_j2 != 0.0:
            z2r2 = (r[2] * r[2]) / (rm * rm)
            c = -1.5 * lambda_j2 * J2_EARTH * MU_EARTH_KM3S2 * re2 / rm ** 5
            g = 1.0 - 5.0 * z2r2
            a_j2 = c * np.array([r[0] * g, r[1] * g, r[2] * (3.0 - 5.0 * z2r2)])
            a = a + a_j2
        if (use_sun or use_moon) and lambda_3body != 0.0:
            a_3b = _third_body_accel(r, t, sun_snap, moon_snap,
                                       lambda_3body=lambda_3body,
                                       synthetic_moon=synthetic_moon,
                                       include_sun=use_sun,
                                       include_moon=use_moon)
            a = a + a_3b
        return np.concatenate([v, a])

    t_cross_list = []
    om_cross_list = []
    t_node_list = []
    omega_node_list = []
    z_prev = x0[2]
    h = np.cross(x0[:3], x0[3:])
    omega_node_prev = math.atan2(-h[0], h[1])
    t_node_list.append(t0_s)
    omega_node_list.append(omega_node_prev)
    x = x0.copy()
    t = t0_s
    n_steps = 0
    t_start_wall = time.time()
    while t < t_end_s:
        x_new = rk4_step(f, t, x, dt_s)
        z_curr = x_new[2]
        if z_prev <= 0 < z_curr:
            frac = -z_prev / (z_curr - z_prev)
            t_cross = t + frac * dt_s
            r_cross = x[:3] + frac * (x_new[:3] - x[:3])
            om_cross = math.atan2(r_cross[1], r_cross[0])
            t_cross_list.append(t_cross)
            om_cross_list.append(om_cross)
        if n_steps % subsample_every == 0:
            h = np.cross(x_new[:3], x_new[3:])
            omega_node = math.atan2(-h[0], h[1])
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
    h = np.cross(x[:3], x[3:])
    omega_node_list.append(math.atan2(-h[0], h[1]))
    t_node_list.append(t)

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
        "mode": mode,
        "lambda_j2": lambda_j2,
        "lambda_3body": lambda_3body,
        "synthetic_moon": synthetic_moon,
    }


# --------------------------------------------------------------------------- #
# Estimators (same as mission_lunisolar_closure)
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
# Forced-secular lunar nodal mode (analytical, alternative explanation)
# --------------------------------------------------------------------------- #
def forced_secular_lunar_nodal_node_rate_deg_day(h_km: float, i_deg: float) -> dict:
    """Analytical forced-secular lunar nodal mode contribution to RAAN rate.

    Derived from Kaula's expansion: the term proportional to cos(i) * cos(Omega - Omega_3)
    that survives short-period averaging but does NOT average to zero under
    the doubly-averaged theory unless the slow nodal regression is correctly
    accounted for.

    Following Murray & Dermott Sec. 7, the secular Hamiltonian contains
    H_3b = (3/4) * (n3/n)^2 * n * a^2 * (1 - e^2)^(-1/2) *
           [ (3 cos^2(i) - 1) * (...) + (3/2) sin^2(i) * cos(2(Omega - Omega_3)) ]

    The forced-secular part is the cos(2(Omega - Omega_3)) term, which has
    period pi/(Omega_dot - Omega_3_dot) ~ 18.6 yr at LEO (since Omega_3 is
    the lunar node, which regresses with period 18.6 yr, the relevant period
    is the synodic period between the satellite's mean Omega and the lunar
    node; at h=600 km i_sso, this is ~6798.4 d / 2 = 3399.2 d).

    Returns the secular rate contribution from this term in deg/day.
    """
    a = R_EARTH_KM + h_km
    n = math.sqrt(MU_EARTH_KM3S2 / a ** 3)
    i_rad = math.radians(i_deg)
    i3_rad = math.radians(SOLAR_OBLIQUITY_DEG + LUNAR_INCLINATION_DEG)

    # Use the Moon's mean motion at its orbital radius (circular approx)
    # n_moon = sqrt(GM_moon / a_moon^3) — but for a point-mass third body,
    # the relevant frequency is the third body's mean motion n3 = sqrt(mu3/a3^3)
    n3_moon = math.sqrt(LUNAR_GM_KM3_S2 / LUNAR_DISTANCE_KM_MEAN ** 3)

    # Forced-secular amplitude (rough; the exact coefficient depends on
    # the Legendre polynomial expansion and inclination functions)
    # Following Murray & Dermott Sec. 7.2:
    # The cos(2(Omega - Omega_3)) term in H_3b yields
    # dOmega/dt contribution ~ -(1/(n a^2 sin i)) * (3/2) sin^2(i3) * sin(i) *
    #   (n3/n)^2 * n * (...) * cos(2(Omega - Omega_3))
    # For the SECULAR average over 18.6 yr, cos(2(Omega - Omega_3)) averages
    # to zero IF the satellite's Omega precesses linearly; in practice the
    # 18.6-yr finite-window fit captures a residual proportional to the
    # window length and the satellite's nodal rate.
    coef = -(3.0 / 2.0) * (n3_moon / n) ** 2 * n
    sin_i = math.sin(i_rad)
    sin2_i3 = math.sin(i3_rad) ** 2

    # The forced-secular contribution to dOmega/dt is roughly:
    # coef * sin2_i3 * sin_i * <cos(2(Omega - Omega_3))>_W
    # where the bracket is the window-averaged value.
    # At W = 18.6 yr, the average is exactly 0 (one full nodal cycle).
    # At W = 1 yr, the average depends on the phase.

    # For this mission's analytical test, we use the standard
    # closed-form secular part WITHOUT the forced-secular term (the
    # 018 corrected formula):
    standard_secular_rad_s = (3.0 / 8.0) * n * (LUNAR_GM_KM3_S2 / MU_EARTH_KM3S2) * (
        a / LUNAR_DISTANCE_KM_MEAN) ** 3 * math.sin(2.0 * (i_rad - i3_rad)) / sin_i
    standard_secular_deg_day = math.degrees(standard_secular_rad_s) * 86400.0

    # Forced-secular upper-bound estimate: maximum amplitude the term
    # can contribute over the 18.6-yr cycle. The exact value depends on
    # phase; we report the BOUND as the diagnostic.
    # Maximum of cos(2(Omega - Omega_3)) = 1, so the bound is the amplitude.
    forced_secular_amplitude_rad_s = coef * sin2_i3 * sin_i
    forced_secular_amplitude_deg_day = math.degrees(forced_secular_amplitude_rad_s) * 86400.0

    return {
        "h_km": h_km,
        "i_deg": i_deg,
        "standard_secular_lunar_deg_day": standard_secular_deg_day,
        "forced_secular_amplitude_bound_deg_day": forced_secular_amplitude_deg_day,
        "ratio_bound_to_standard": (
            forced_secular_amplitude_deg_day / standard_secular_deg_day
            if abs(standard_secular_deg_day) > 1e-30 else float("nan")
        ),
        "note": "amplitude is the max value of the forced-secular term over the 18.6-yr nodal cycle; average over a full nodal cycle is zero; finite-window fits capture a phase-dependent residual",
    }


# --------------------------------------------------------------------------- #
# Octupole (a/a3)^4 third-body secular contribution (alternative explanation)
# --------------------------------------------------------------------------- #
def octupole_lunisolar_raan_rate_rad_s(h_km: float, i_deg: float) -> dict:
    """Octupole (Legendre P_3) term of the third-body disturbing function.

    The l=3 term in the Legendre expansion contributes (Kaula, Murray & Dermott):
        R_3 = (G m_3 a^3 / (8 r_3^4)) * [5 cos^3(i - i_3) - 3 cos(i - i_3)]
    Differentiating w.r.t. i gives the secular octupole nodal rate.

    The leading octupole contribution to dOmega/dt is (Convention B):
        dOmega/dt|_octupole = +(3/8) n (m_3/m_E) (a/a_3)^4 * ...
    The exact coefficient depends on the inclination function; for the
    i_sso geometry at h=600 km, the octupole is ~ (a/a_moon) ≈ 1.8e-2 of
    the quadrupole. At h=600 km, (a/a_moon) = 0.01815, so the octupole is
    ~ 5e-7 of the quadrupole (a fifth-power scaling).

    Returns the secular octupole contribution to RAAN rate in deg/day.
    """
    a = R_EARTH_KM + h_km
    n = math.sqrt(MU_EARTH_KM3S2 / a ** 3)
    i_rad = math.radians(i_deg)
    i3_rad = math.radians(SOLAR_OBLIQUITY_DEG + LUNAR_INCLINATION_DEG)
    delta = i_rad - i3_rad

    # Murray & Dermott Eq. 7.10 form, octupole correction
    # R_3 = (3/8) (mu_3 a^3 / r_3^4) (1/2) [5 cos^3(delta) - 3 cos(delta)]
    # dR_3/di = (3/8) (mu_3 a^3 / r_3^4) (1/2) [-3 sin(delta) (5 cos^2(delta) - 1)]
    # dOmega/dt = -(1/(n a^2 sin i)) dR_3/di
    # = +(3/16) (mu_3 a / (n r_3^4)) (sin(delta) (5 cos^2(delta) - 1) / sin i)
    # Note: the sign convention is Convention B with the minus sign in
    # dOmega/dt = -(1/(n a^2 sin i)) dR/di

    # Simplified: the octupole secular RAAN rate scales as (a/a3)^4 ~ 1e-7 of the
    # quadrupole. We compute it for diagnostic purposes; the expected magnitude
    # is too small to explain the 018/021 discrepancy.
    a3 = LUNAR_DISTANCE_KM_MEAN
    coef_rad_s = (3.0 / 16.0) * (LUNAR_GM_KM3_S2 / (n * a3 ** 4)) * a / math.sin(i_rad)
    octupole_rad_s = coef_rad_s * math.sin(delta) * (5 * math.cos(delta) ** 2 - 1)
    octupole_deg_day = math.degrees(octupole_rad_s) * 86400.0

    return {
        "h_km": h_km,
        "i_deg": i_deg,
        "octupole_lunar_deg_day": octupole_deg_day,
        "scaling_factor_vs_quadrupole": (a / a3) ** 1,  # (a/a3)^1 additional factor vs quadrupole (a/a3)^3
        "note": "octupole is (a/a3) ~ 1.8% of quadrupole at h=600 km; expected to be too small to explain 170x discrepancy",
    }


# --------------------------------------------------------------------------- #
# Code hashes for provenance
# --------------------------------------------------------------------------- #
def code_hashes() -> dict:
    here = Path(__file__).resolve().parent
    lab_root = here
    while lab_root != lab_root.parent:
        if (lab_root / "src" / "lab_utils").exists():
            break
        lab_root = lab_root.parent
    files = {
        "mission_experiment.py": here / "mission_experiment.py",
        "lab_utils/orbits.py": lab_root / "src" / "lab_utils" / "orbits.py",
        "lab_utils/earth_frames.py": lab_root / "src" / "lab_utils" / "earth_frames.py",
        "lab_utils/integrators.py": lab_root / "src" / "lab_utils" / "integrators.py",
        "lab_utils/results.py": lab_root / "src" / "lab_utils" / "results.py",
        "lab_utils/__init__.py": lab_root / "src" / "lab_utils" / "__init__.py",
        "parent_mission_experiment.py": here.parent / "mission_lunisolar_closure" / "experiment.py",
        "sun_reference_snapshot.txt": SUN_SNAPSHOT,
        "moon_reference_snapshot.txt": MOON_SNAPSHOT,
    }
    return {name: _sha256(p) for name, p in files.items() if p.exists()}
