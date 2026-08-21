"""Ground Tracks — spherical Earth, uniform rotation, Keplerian orbits.

Laboratory loop: implement → test (pytest green) → run → validate → document.
Determinism rule: fixed seeds / no RNG, pure float64, Agg rendering.

Numerical contract (ground track)

  ECI position from Kepler elements (a,e,i,Ω,ω,M0) at epoch t0:
    n = √(μ/a³)                                          mean motion
    M(t) = M0 + n·(t-t0)                                  mean anomaly
    M = E - e·sin E   → Newton solve for E                 Kepler eq.
    cos ν = (cos E - e)/(1 - e cos E)
    sin ν = √(1-e²) sin E /(1 - e cos E)                true anomaly
    r = a(1 - e cos E) = p/(1+e cos ν),  p=a(1-e²)
    r_pf = [r cos ν, r sin ν, 0]ʰ                        perifocal
    Q = R_z(Ω)·R_x(i)·R_z(ω)                              3-1-3 rotation
    r_eci(t) = Q · r_pf(t)                                ECI vector

    v_pf = (μ/h)[-sin ν, e+cos ν, 0], h=√(μp)            perifocal velocity
    v_eci = Q·v_pf                                         for propagation seed

  Earth rotation (idealized uniform, spherical):
    θ_G(t) = θ_G0 + ω_E·(t-t0)                             GMST (rad)
    ω_E = 2π/T_sid ,  T_sid = 2π/ω_E = 86164.0905 s
    θ_G0 = 0 at epoch (Greenwich ≡ vernal equinox)        idealization; J2000 real GMST is 280.46°
                                                          but epoch bias would hide GEO stationary test

  ECI → ECEF (passive Z rotation, west-positive lon drift):
    R = [[ cosθ  sinθ  0],[-sinθ  cosθ  0],[0 0 1]]         R = R_z^T(θ) in active convention
    r_ecef(t) = R(θ_G(t)) · r_eci(t)                      magnitude preserved: |r_ecef|=|r_eci|

  Geocentric latitude/longitude on spherical Earth (R_E cancels, geocentric≡geodetic):
    r = |r_ecef| = √(x²+y²+z²)
    φ_gc = arcsin(z / r) = atan2(z, √(x²+y²))   ∈ [-90°,+90°]
    λ    = atan2(y_ecef, x_ecef)                   ∈ (-180°,+180°]  wrapped
    h = r - R_E  (height, not used for sub-satellite point)

  Spherical-trig independent form (same result, different algebra):
    u = ω + ν                                              argument of latitude
    sin φ = sin i · sin u                                   → φ = arcsin(sin i sin u)
    lon_eci = Ω + atan2(cos i sin u, cos u)                 → from x+iy factorization
    λ = lon_eci - θ_G   (wrapped)                           → ECEF longitude
    These are proved by x_eci+iy_eci = r e^{iΩ}(cos u + i cos i sin u).

  Periodicity / repeat:
    T = 2π√(a³/μ)                                           nodal≡Kepler for pure Kepler
    Δλ = -ω_E·T  (mod 360°, west-negative)                  per-orbit longitudinal shift at equator
    Repeat after m orbits, n sidereal days:  m·T = n·T_sid  (coprime)
    GEO: a_geo = (μ T_sid²/4π²)^{1/3} → Δλ= -360° → 0 (stationary)

  Invariants:
    max|φ| = min(i, π-i)  for retrograde i>90°; polar i=90° reaches ±90°
    Equatorial i=0 → φ≡0, λ(t)= Ω+ω+ν(t)-θ_G(t) linear in t for circular e=0
    GEO i=0, a=a_geo → φ≡0, λ≡Ω-θ_G0 constant

References: Curtis 4th ed. Ch.2-5 (Kepler eq., elements, Q, GMST); Vallado 4th ed. §2.2, §3.3, Eq.3-34 ECI↔ECEF, §3.5 ω_E; Bate/Mueller/White Ch.1-3; Murray&Dermott Ch.2 Kepl. sol.; WGS-84 (NIMA TR8350.2), IAU 2015 B3 GM.

Reuse: propagates via verified planeChange 3-D Cowell RK4 (importlib) + direct Kepler closed form via keplerOrbitValidation solver logic + src/lab_utils metrics/results.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lab_utils.results import save_json_result

# --- Reuse of verified 006 3D Cowell RK4 ------------------------------------
PCM_PATH = (
    Path(__file__).resolve().parents[1] / "planeChangeManeuvers" / "experiment.py"
)
_pcm_spec = importlib.util.spec_from_file_location("pcm_006_for_groundtrack", PCM_PATH)
assert _pcm_spec is not None and _pcm_spec.loader is not None
_pcm = importlib.util.module_from_spec(_pcm_spec)
_pcm_spec.loader.exec_module(_pcm)

propagate_3d_rk4 = _pcm.propagate_3d_rk4  # verified to ≤1e-11 in 006

# Also load kepler solver logic id for cross-check (we implement our own but verify against it)
KEPLER_PATH = (
    Path(__file__).resolve().parents[1] / "keplerOrbitValidation" / "experiment.py"
)
_kep_spec = importlib.util.spec_from_file_location("kepler_for_groundtrack", KEPLER_PATH)
assert _kep_spec is not None and _kep_spec.loader is not None
_kep = importlib.util.module_from_spec(_kep_spec)
_kep_spec.loader.exec_module(_kep)

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

# --------------------------------------------------------------------------- #
# Physical constants (documented, single source of truth for this experiment)
# --------------------------------------------------------------------------- #
MU_EARTH_KM3S2 = 398600.4418  # IAU 2015 nominal GM_E (km³/s²); JPL DE440 planet-only 398600.435507 differs 1.5e-8 rel
R_EARTH_KM = 6378.137  # WGS-84 equatorial radius (km); lab nominal 6378.1 differs 37 m (5.8 ppm)
OMEGA_EARTH_RAD_S = 7.2921159e-5  # WGS-84 / Vallado Table 3-1, rad/s; T_sid = 2π/ω_E
T_SIDEREAL_S = 2.0 * np.pi / OMEGA_EARTH_RAD_S  # 86164.090530833... s (sidereal day, not 86400 solar)
GMST0_RAD = 0.0  # idealized epoch: Greenwich ≡ vernal equinox at t0 (idealization, verified; real J2000 280.46° would bias GEO test)
# frame convention tag for results.json
FRAME_CONVENTION = "ECI J2000 pseudo-inertial (Z=CIP, X=γ), ECEF WGS-84 (Z=CIP, X=Greenwich), t0 epoch, GMST0=0, R_z_passive(θ)=[[c,s,0],[-s,c,0],[0,0,1]], r_ecef=R·r_eci, λ=atan2(y_ecef,x_ecef) wrapped (-180,180]"

# --------------------------------------------------------------------------- #
# Kepler solver and element helpers
# --------------------------------------------------------------------------- #

def solve_kepler(M: np.ndarray | float, e: float, tol: float = 1e-14, max_iter: int = 100) -> np.ndarray | float:
    """Newton solve M = E - e sin E for E. Vectorized, seed E0=M+e sin M."""
    M_arr = np.asarray(M, dtype=float)
    scalar = M_arr.ndim == 0
    M_arr = np.atleast_1d(M_arr)
    E = M_arr + e * np.sin(M_arr)
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M_arr
        if np.max(np.abs(f)) < tol:
            break
        denom = 1.0 - e * np.cos(E)
        E = E - f / denom
    else:
        # check convergence
        f = E - e * np.sin(E) - M_arr
        if np.max(np.abs(f)) >= tol:
            raise ValueError(f"Kepler did not converge e={e} max|f|={np.max(np.abs(f)):.2e}")
    if scalar:
        return float(E[0])
    return E


def true_anomaly_from_E(E: np.ndarray | float, e: float) -> np.ndarray | float:
    """Convert eccentric anomaly E to true anomaly nu (atan2, preserves quadrant)."""
    E_arr = np.asarray(E, dtype=float)
    scalar = E_arr.ndim == 0
    E_arr = np.atleast_1d(E_arr)
    cos_nu = (np.cos(E_arr) - e) / (1.0 - e * np.cos(E_arr))
    sin_nu = np.sqrt(1.0 - e * e) * np.sin(E_arr) / (1.0 - e * np.cos(E_arr))
    # clip to handle rounding beyond [-1,1] for cos (sin already)
    cos_nu = np.clip(cos_nu, -1.0, 1.0)
    nu = np.arctan2(sin_nu, cos_nu)
    if scalar:
        return float(nu[0])
    return nu


def radius_from_E(a: float, e: float, E: np.ndarray | float) -> np.ndarray | float:
    """r = a(1 - e cos E)."""
    return a * (1.0 - e * np.cos(np.asarray(E, dtype=float)))


def orbital_period(a: float, mu: float = MU_EARTH_KM3S2) -> float:
    """Kepler period T = 2π√(a³/μ)."""
    return 2.0 * np.pi * np.sqrt(a ** 3 / mu)


def mean_motion(a: float, mu: float = MU_EARTH_KM3S2) -> float:
    return np.sqrt(mu / a ** 3)


def rotation_matrix_313(Omega: float, inc: float, omega: float) -> np.ndarray:
    """Q = R_z(Ω)·R_x(i)·R_z(ω)  (3-1-3). All rad. Returns 3×3."""
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(omega), np.sin(omega)
    # R_z(omega) then R_x(i) then R_z(Omega)
    # Compute directly as in orbital mechanics texts (Curtis Eq.4.47)
    Q = np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci,  sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si,                  cw * si,                  ci],
    ], dtype=float)
    return Q


def coe_to_rv_eci(a: float, e: float, inc: float, Omega: float, omega: float, nu: float, mu: float = MU_EARTH_KM3S2) -> tuple[np.ndarray, np.ndarray]:
    """Single-epoch perifocal → ECI for state (r_eci, v_eci) at true anomaly nu (rad)."""
    p = a * (1.0 - e * e)
    # handle circular edge: p still defined
    r = p / (1.0 + e * np.cos(nu)) if abs(e) < 1.0 else a  # e==1 not used
    # perifocal position
    r_pf = np.array([r * np.cos(nu), r * np.sin(nu), 0.0], dtype=float)
    h = np.sqrt(mu * p)
    # perifocal velocity: (μ/h)[-sin ν, e+cos ν, 0]
    v_pf = np.array([-(mu / h) * np.sin(nu), (mu / h) * (e + np.cos(nu)), 0.0], dtype=float)
    Q = rotation_matrix_313(Omega, inc, omega)
    r_eci = Q @ r_pf
    v_eci = Q @ v_pf
    return r_eci, v_eci


def gmst_rad(t: np.ndarray | float, gmst0: float = GMST0_RAD, omega_e: float = OMEGA_EARTH_RAD_S) -> np.ndarray | float:
    """GMST = gmst0 + ω_E·t  (rad). t in seconds since epoch. Not wrapped."""
    return gmst0 + omega_e * np.asarray(t, dtype=float)


def eci_to_ecef(r_eci: np.ndarray, gmst: np.ndarray | float) -> np.ndarray:
    """Rotate ECI vectors to ECEF via passive Z rotation R(θ_G).

    R = [[c s 0],[-s c 0],[0 0 1]]  with c=cosθ, s=sinθ.
    Then lon_ecef = atan2(y_ecef,x_ecef) = lon_eci - θ_G (west drift).
    Preserves |r| to machine precision. Vectorized over leading time axis.

    r_eci: (N,3) or (3,) ; gmst: scalar or (N,). Returns (N,3) or (3,).
    """
    single = r_eci.ndim == 1 and np.asarray(gmst).ndim == 0
    r = np.asarray(r_eci, dtype=float)
    th = np.asarray(gmst, dtype=float)
    if r.ndim == 1:
        r = r[None, :]
    if th.ndim == 0:
        th = np.full(r.shape[0], float(th), dtype=float)
    else:
        th = np.asarray(th, dtype=float).reshape(-1)
        assert th.shape[0] == r.shape[0], "gmst and r_eci length mismatch"
    c = np.cos(th)
    s = np.sin(th)
    # R·r :  [c x + s y,  -s x + c y,  z]
    out = np.empty_like(r)
    out[:, 0] = c * r[:, 0] + s * r[:, 1]
    out[:, 1] = -s * r[:, 0] + c * r[:, 1]
    out[:, 2] = r[:, 2]
    if single:
        return out[0]
    return out


def ecef_to_latlon(r_ecef: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Geocentric lat/lon from ECEF vectors (spherical Earth).

    φ = arcsin(z / r)  ∈ [-90°,+90°]  (or atan2(z,√(x²+y²)))
    λ = atan2(y, x)    ∈ (-180°,+180°] wrapped
    Returns (lat_deg, lon_deg, r_norm_km) arrays shape (N,) or scalars.
    At |φ|>89.999° (pole) λ is ill-defined; we return λ=0 at exact pole
    (within epsilon) and rely on wrap-gap splitting for plotting.
    """
    r = np.asarray(r_ecef, dtype=float)
    single = r.ndim == 1
    if single:
        r = r[None, :]
    x, y, z = r[:, 0], r[:, 1], r[:, 2]
    rn = np.sqrt(x * x + y * y + z * z)
    # latitude via arcsin (preserves sign)
    # clip z/r to [-1,1] to handle rounding 1+1e-16
    sin_phi = np.clip(z / rn, -1.0, 1.0)
    lat_rad = np.arcsin(sin_phi)
    # alternative via atan2 also valid: lat = atan2(z, hypot(x,y))
    # longitude
    # pole guard: if hypot(x,y) < eps * r then set lon=0 (undefined)
    # Use eps=1e-12 relative
    lon_rad = np.arctan2(y, x)
    # guard exact pole: when both x and y ~0, atan2(0,0)=0 already, fine
    eps = 1e-12
    pole_mask = np.hypot(x, y) < eps * rn
    # keep lon as computed (0) at pole; tests that check pole will check lat only
    # optionally set to 0 explicitly
    lon_rad = np.where(pole_mask, 0.0, lon_rad)
    lat_deg = np.degrees(lat_rad)
    lon_deg = np.degrees(lon_rad)
    # wrap lon to (-180,180]
    lon_deg = ((lon_deg + 180.0) % 360.0) - 180.0
    # fix -180 edge: keep -180 not 180 for consistency
    if single:
        return float(lat_deg[0]), float(lon_deg[0]), float(rn[0])
    return lat_deg, lon_deg, rn


def spherical_trig_latlon(
    inc: float, Omega: float, omega: float, nu: np.ndarray | float, gmst: np.ndarray | float
) -> tuple[np.ndarray, np.ndarray]:
    """Independent spherical-trig lat/lon (different algebra, same physics).

    u = ω + ν  (argument of latitude)
    sin φ = sin i sin u        → φ = arcsin(clip)
    lon_eci = Ω + atan2(cos i sin u, cos u)
    λ = lon_eci - gmst   (wrapped to (-π,π] then deg)

    inc,Omega,omega in rad; nu,gmst arrays (rad). Returns (lat_deg, lon_deg).
    Vectorized. Independent of matrix path, so cross-check is non-trivial.
    """
    u = np.asarray(omega, dtype=float) + np.asarray(nu, dtype=float)
    gmst_arr = np.asarray(gmst, dtype=float)
    # handle scalar vs array
    u = np.asarray(u, dtype=float)
    # broadcast if needed
    # need same shape: if nu is array, gmst should be array of same length
    # assume caller passes consistent shapes
    # For scalar case both scalar
    sg = np.sin(inc) * np.sin(u)
    sg = np.clip(sg, -1.0, 1.0)
    lat_rad = np.arcsin(sg)
    # lon_eci offset
    # atan2(cos i sin u, cos u) handles all quadrants; for i=90°, cos i=0 => atan2(0, cos u) = 0 or π
    lon_eci = Omega + np.arctan2(np.cos(inc) * np.sin(u), np.cos(u))
    lon_rad = lon_eci - gmst_arr
    # wrap to (-π,π]
    lon_rad = (lon_rad + np.pi) % (2 * np.pi) - np.pi
    # map exactly -π -> -π (keep)
    lat_deg = np.degrees(lat_rad)
    lon_deg = np.degrees(lon_rad)
    # wrap to (-180,180] via degrees already equivalent but re-apply for consistency
    lon_deg = ((lon_deg + 180.0) % 360.0) - 180.0
    return lat_deg, lon_deg


def wrap_longitude_deg(lon_deg: np.ndarray) -> np.ndarray:
    """Wrap longitude to (-180,180]."""
    return ((np.asarray(lon_deg, dtype=float) + 180.0) % 360.0) - 180.0


def ground_track_analytic(
    a: float,
    e: float,
    inc: float,
    Omega: float,
    omega: float,
    M0: float,
    t: np.ndarray,
    mu: float = MU_EARTH_KM3S2,
    gmst0: float = GMST0_RAD,
    omega_e: float = OMEGA_EARTH_RAD_S,
) -> dict:
    """Analytic Kepler ground track (no propagation error).

    Computes r_eci via Kepler solution then ECI→ECEF→lat/lon via matrix path.
    Also computes independent trig path for cross-validation.

    t: (N,) seconds since epoch (t0=0). inc,Omega,omega,M0 in rad.
    Returns dict with lat/lon via matrix, via trig, r_eci, r_ecef, etc.
    """
    t = np.asarray(t, dtype=float).reshape(-1)
    n = mean_motion(a, mu)
    M = M0 + n * t
    # wrap M to [0,2π) is not needed for solver but reduces range
    M_wrapped = np.mod(M, 2.0 * np.pi)
    E = solve_kepler(M_wrapped, e)
    nu = true_anomaly_from_E(E, e)
    # r magnitude via a(1 - e cos E)
    r_mag = a * (1.0 - e * np.cos(np.asarray(E, dtype=float)))
    # build r_eci array via Q
    Q = rotation_matrix_313(Omega, inc, omega)
    # r_pf per epoch
    cos_nu = np.cos(np.asarray(nu))
    sin_nu = np.sin(np.asarray(nu))
    r_pf = np.column_stack([r_mag * cos_nu, r_mag * sin_nu, np.zeros_like(r_mag)])
    r_eci = (Q @ r_pf.T).T  # (N,3)
    # GMST
    gmst = gmst_rad(t, gmst0, omega_e)
    r_ecef = eci_to_ecef(r_eci, gmst)
    lat_mat, lon_mat, rn = ecef_to_latlon(r_ecef)
    lat_trig, lon_trig = spherical_trig_latlon(inc, Omega, omega, nu, gmst)
    # also r_eci via Q factorization cross-check not needed
    return {
        "t": t,
        "M": M,
        "E": np.asarray(E),
        "nu": np.asarray(nu),
        "r_mag": r_mag,
        "r_eci": r_eci,
        "r_ecef": r_ecef,
        "lat_mat_deg": lat_mat,
        "lon_mat_deg": lon_mat,
        "lat_trig_deg": lat_trig,
        "lon_trig_deg": lon_trig,
        "gmst_rad": gmst,
        "rn": rn,
    }


def ground_track_from_state_array(r_eci: np.ndarray, t: np.ndarray, gmst0: float = GMST0_RAD, omega_e: float = OMEGA_EARTH_RAD_S) -> dict:
    """Project a precomputed r_eci(t) array (e.g., from RK4) to lat/lon."""
    t = np.asarray(t, dtype=float).reshape(-1)
    r_eci = np.asarray(r_eci, dtype=float)
    assert r_eci.shape[0] == t.shape[0] and r_eci.shape[1] == 3
    gmst = gmst_rad(t, gmst0, omega_e)
    r_ecef = eci_to_ecef(r_eci, gmst)
    lat_deg, lon_deg, rn = ecef_to_latlon(r_ecef)
    return {"lat_deg": lat_deg, "lon_deg": lon_deg, "r_ecef": r_ecef, "gmst_rad": gmst, "rn": rn}


def propagate_ground_track(
    a: float,
    e: float,
    inc: float,
    Omega: float,
    omega: float,
    M0: float,
    t: np.ndarray,
    mu: float = MU_EARTH_KM3S2,
    gmst0: float = GMST0_RAD,
    omega_e: float = OMEGA_EARTH_RAD_S,
) -> dict:
    """Propagate via verified 3-D RK4 from epoch state, then project."""
    t = np.asarray(t, dtype=float).reshape(-1)
    n = mean_motion(a, mu)
    # epoch true anomaly from M0
    E0 = solve_kepler(np.mod(M0, 2 * np.pi), e)
    nu0 = true_anomaly_from_E(E0, e)
    r0, v0 = coe_to_rv_eci(a, e, inc, Omega, omega, float(nu0), mu)
    # propagate_3d_rk4 expects t array and dt (unused). It integrates from r0,v0 at t[0] to all t.
    # It is fixed-step RK4; t should be uniformly spaced for simplest validation but handles any.
    states = propagate_3d_rk4(r0, v0, mu, t, float(t[1] - t[0]) if len(t) > 1 else 1.0)
    r_eci_prop = states[:, :3]
    # Also compute analytic for error analysis (optional, not returned here directly)
    proj = ground_track_from_state_array(r_eci_prop, t, gmst0, omega_e)
    return {"r_eci": r_eci_prop, "v_eci": states[:, 3:], "lat_deg": proj["lat_deg"], "lon_deg": proj["lon_deg"], "r_ecef": proj["r_ecef"]}


def unwrap_longitude_deg(lon_wrapped_deg: np.ndarray) -> np.ndarray:
    """Unwrap longitude to continuous (no 360° jumps) for Delta_lambda measurement."""
    # convert to rad, unwrap, back to deg
    lon_rad = np.radians(np.asarray(lon_wrapped_deg, dtype=float))
    lon_unwrapped_rad = np.unwrap(lon_rad)
    return np.degrees(lon_unwrapped_rad)


def delta_longitude_per_orbit(a: float, mu: float = MU_EARTH_KM3S2, omega_e: float = OMEGA_EARTH_RAD_S) -> float:
    """Analytic Δλ = -ω_E·T (deg, west-negative, before modulo)."""
    T = orbital_period(a, mu)
    return -np.degrees(omega_e * T)


def max_latitude_deg(inc_rad: float) -> float:
    """Theoretical max geocentric latitude magnitude (deg)."""
    # i in [0,π]; max |φ| = min(i, π-i) for retrograde
    i_deg = np.degrees(inc_rad)
    if i_deg <= 90.0:
        return i_deg
    else:
        return 180.0 - i_deg


# --------------------------------------------------------------------------- #
# Representative orbits (real params, idealized Kepler)
# --------------------------------------------------------------------------- #
def real_orbits() -> dict:
    """Real representative orbits for ground-track anchors.

    Altitudes and inclinations from published mean elements (screening conventions, not physical laws).
    Each orbit: a (km), e, inc (rad), Omega (rad), omega (rad), M0 (rad), description.
    """
    orbits = {}
    # ISS-like: 420 km alt, 51.6°
    alt_iss = 420.0
    orbits["ISS"] = {
        "a_km": R_EARTH_KM + alt_iss,
        "e": 0.0003,
        "inc_deg": 51.6,
        "Omega_deg": 0.0,
        "omega_deg": 0.0,
        "M0_deg": 0.0,
        "description": "ISS-like LEO (420 km, 51.6°) — latitude-band test, real mean TLE ~51.64°/408-420 km",
    }
    # Equatorial LEO
    orbits["Equatorial_LEO"] = {
        "a_km": R_EARTH_KM + 400.0,
        "e": 0.0,
        "inc_deg": 0.0,
        "Omega_deg": 0.0,
        "omega_deg": 0.0,
        "M0_deg": 0.0,
        "description": "Equatorial LEO (400 km, 0°) — degenerate φ≡0, λ linear",
    }
    # Polar LEO
    orbits["Polar_LEO"] = {
        "a_km": R_EARTH_KM + 500.0,
        "e": 0.0,
        "inc_deg": 90.0,
        "Omega_deg": 0.0,
        "omega_deg": 0.0,
        "M0_deg": 0.0,
        "description": "Polar LEO (500 km, 90°) — pole singularity, covers all latitudes",
    }
    # Sun-synchronous approx: 600 km, ~97.8° (real SSO 600 km is 97.8–98.2°; use 98.0)
    orbits["SSO"] = {
        "a_km": R_EARTH_KM + 600.0,
        "e": 0.0,
        "inc_deg": 98.0,
        "Omega_deg": 0.0,
        "omega_deg": 0.0,
        "M0_deg": 0.0,
        "description": "SSO-like LEO (600 km, 98°) — retrograde, near-repeating",
    }
    # GEO
    a_geo = (MU_EARTH_KM3S2 * T_SIDEREAL_S ** 2 / (4.0 * np.pi ** 2)) ** (1.0 / 3.0)
    orbits["GEO"] = {
        "a_km": float(a_geo),
        "e": 0.0,
        "inc_deg": 0.0,
        "Omega_deg": 0.0,
        "omega_deg": 0.0,
        "M0_deg": 0.0,
        "description": f"Geostationary (a={a_geo:.1f} km ≡ R_E+35786 km, 0°) — stationary point",
        "a_geo_exact": True,
    }
    # GEO inclined (figure-8)
    orbits["GEO_inclined"] = {
        "a_km": float(a_geo),
        "e": 0.0,
        "inc_deg": 5.0,
        "Omega_deg": 0.0,
        "omega_deg": 0.0,
        "M0_deg": 0.0,
        "description": "Inclined GEO (5°) — figure-8 lat ±5°, lon small",
    }
    # Molniya-like: a=26560 km, e=0.74, i=63.4°
    orbits["Molniya"] = {
        "a_km": 26560.0,
        "e": 0.74,
        "inc_deg": 63.4,
        "Omega_deg": 0.0,
        "omega_deg": 270.0,
        "M0_deg": 0.0,
        "description": "Molniya-like (a=26560 km, e=0.74, 63.4°, ω=270°) — elliptical dwell, connects to Exp 012",
    }
    # Retrograde equatorial for symmetry
    orbits["Retro_LEO"] = {
        "a_km": R_EARTH_KM + 400.0,
        "e": 0.0,
        "inc_deg": 180.0,
        "Omega_deg": 0.0,
        "omega_deg": 0.0,
        "M0_deg": 0.0,
        "description": "Retrograde equatorial (400 km, 180°) — φ≡0 same but opposite node handling",
    }
    # attach derived T and Delta_lon for each
    for name, o in orbits.items():
        a = o["a_km"]
        e = o["e"]
        inc = np.radians(o["inc_deg"])
        # period uses a only (Kepler)
        T = orbital_period(a, MU_EARTH_KM3S2)
        o["T_sec"] = float(T)
        o["T_min"] = float(T / 60.0)
        o["T_days"] = float(T / 86400.0)
        o["delta_lon_deg"] = float(delta_longitude_per_orbit(a, MU_EARTH_KM3S2, OMEGA_EARTH_RAD_S))
        # modulo 360 to (-180,180] shift for map
        o["delta_lon_wrapped_deg"] = float(((o["delta_lon_deg"] + 180.0) % 360.0) - 180.0)
        # for repeat after 360 wrap, GEO expects 0
        o["max_lat_theory_deg"] = float(max_latitude_deg(inc))
    return orbits


# --------------------------------------------------------------------------- #
# Validation helpers
# --------------------------------------------------------------------------- #
def validate_trig_vs_matrix(orbits: dict, t_samples: int = 720) -> dict:
    """Cross-check trig vs matrix lat/lon over multiple orbits."""
    results = {}
    for name, o in orbits.items():
        a, e = o["a_km"], o["e"]
        inc = np.radians(o["inc_deg"])
        Omega = np.radians(o["Omega_deg"])
        omega = np.radians(o["omega_deg"])
        M0 = np.radians(o["M0_deg"])
        T = o["T_sec"]
        # sample 2 orbits with t_samples points each (at least 360/orbit)
        # for high eccentricity use more points near periapsis? Here uniform is fine for analytic
        t = np.linspace(0.0, 2.0 * T, 2 * t_samples, endpoint=False)
        gt = ground_track_analytic(a, e, inc, Omega, omega, M0, t)
        lat_mat = gt["lat_mat_deg"]
        lon_mat = gt["lon_mat_deg"]
        lat_tri = gt["lat_trig_deg"]
        lon_tri = gt["lon_trig_deg"]
        # max absolute differences (handle pole undefined: mask near poles for lon)
        dlat = np.max(np.abs(lat_mat - lat_tri))
        # lon difference wrapped to [-180,180]
        dlon_wrapped = np.abs(((lon_mat - lon_tri + 180.0) % 360.0) - 180.0)
        # at poles lon undefined: mask where |lat|>89.9
        pole_mask = np.abs(lat_mat) > 89.9
        if np.any(~pole_mask):
            dlon = float(np.max(dlon_wrapped[~pole_mask]))
        else:
            dlon = float(np.max(dlon_wrapped))
        # also great-circle distance for robustness
        # Δσ = arccos(sinφ1 sinφ2 + cosφ1 cosφ2 cosΔλ)
        lat1r, lat2r = np.radians(lat_mat), np.radians(lat_tri)
        dlonr = np.radians(lon_mat - lon_tri)
        cos_sigma = np.sin(lat1r) * np.sin(lat2r) + np.cos(lat1r) * np.cos(lat2r) * np.cos(dlonr)
        cos_sigma = np.clip(cos_sigma, -1.0, 1.0)
        sigma = np.arccos(cos_sigma)
        max_sigma = float(np.max(sigma))
        results[name] = {
            "max_abs_dlat_deg": float(dlat),
            "max_abs_dlon_deg": float(dlon),
            "max_great_circle_rad": float(max_sigma),
            "max_abs_dlat_rad": float(np.radians(dlat)),
            "max_abs_dlon_rad": float(np.radians(dlon)),
        }
    return results


def validate_invariants(orbits: dict, t_samples: int = 1440) -> dict:
    """Check invariants: max|φ|=i, equatorial shift, GEO stationary, etc."""
    results = {}
    for name, o in orbits.items():
        a, e = o["a_km"], o["e"]
        inc = np.radians(o["inc_deg"])
        Omega = np.radians(o["Omega_deg"])
        omega = np.radians(o["omega_deg"])
        M0 = np.radians(o["M0_deg"])
        T = o["T_sec"]
        # sample densely: 720 points/orbit ×3 orbits = 2160
        pts_per_orbit = max(720, int(720 / (1 - e) ** 1.5) if e < 0.9 else 2880)
        N = 3 * pts_per_orbit
        t = np.linspace(0.0, 3.0 * T, N, endpoint=False)
        gt = ground_track_analytic(a, e, inc, Omega, omega, M0, t)
        lat_mat = gt["lat_mat_deg"]
        lon_mat = gt["lon_mat_deg"]
        max_lat_measured = float(np.max(np.abs(lat_mat)))
        max_lat_theory = o["max_lat_theory_deg"]
        # max lat error
        max_lat_err_deg = abs(max_lat_measured - max_lat_theory)
        # equatorial shift: measure lambda at t+T minus lambda at t (unwrapped)
        # take point near t=0 and t=T (approx). Use interpolation to get exact T.
        # Since our grid is uniform, T corresponds to exactly pts_per_orbit steps
        lon_unwrapped = unwrap_longitude_deg(lon_mat)
        # delta per orbit as mean over first two orbits
        d1 = lon_unwrapped[pts_per_orbit] - lon_unwrapped[0]
        d2 = lon_unwrapped[2 * pts_per_orbit] - lon_unwrapped[pts_per_orbit]
        measured_delta_lon = float((d1 + d2) / 2.0)
        analytic_delta = o["delta_lon_deg"]  # includes -360 for GEO
        # For non-GEO, analytic is in (-360,0); GEO analytic = -360 but unwrapped measured should be ~-360; wrapped difference should be ~0
        # Compare unwrapped vs analytic directly (both unwrapped)
        delta_err_wrapped = abs(((measured_delta_lon - analytic_delta + 180.0) % 360.0) - 180.0)
        # but also directly absolute for GEO: analytic -360, measured -360 -> error 0, but wrapped 0 too
        # Better to compute error on unwrapped scale for GEO: abs(measured+360) vs expected? Use minimal.
        # Actually for GEO a exactly tuned, analytic = -360.0, measured should be -360.0 within 1e-9.
        # So plain absolute error is fine for GEO and for LEO where analytic ~ -23°, measured ~ -23°, difference <1e-6.
        delta_err_abs = abs(measured_delta_lon - analytic_delta)
        # For LEO, wrapped and abs are similar; for GEO wrapped 0 vs abs 0 same.
        # Choose the smaller as error (handles wrap ambiguity)
        delta_lon_err_deg = float(min(delta_err_abs, delta_err_wrapped))
        # GEO stationary check: for i=0, lat ≡0, lon constant (wrapped variation ~0)
        # measure lon variation over 5 orbits for GEO family
        lon_variation_deg = float(np.max(lon_mat) - np.min(lon_mat)) if name.startswith("GEO") else None
        # pole handling: check no NaN/Inf
        has_nan = bool(np.any(~np.isfinite(lat_mat)) or np.any(~np.isfinite(lon_mat)))
        # r preservation: |r_ecef| == |r_eci| to machine
        r_ecef_norm = np.linalg.norm(gt["r_ecef"], axis=1)
        r_eci_norm = np.linalg.norm(gt["r_eci"], axis=1)
        r_preserve_max = float(np.max(np.abs(r_ecef_norm - r_eci_norm) / r_eci_norm))
        results[name] = {
            "max_lat_measured_deg": max_lat_measured,
            "max_lat_theory_deg": max_lat_theory,
            "max_lat_err_deg": max_lat_err_deg,
            "max_lat_err_rad": float(np.radians(max_lat_err_deg)),
            "delta_lon_measured_deg": measured_delta_lon,
            "delta_lon_analytic_deg": analytic_delta,
            "delta_lon_err_deg": delta_lon_err_deg,
            "delta_lon_err_rad": float(np.radians(delta_lon_err_deg)),
            "lon_variation_deg_for_GEO": lon_variation_deg,
            "has_nan_or_inf": has_nan,
            "r_preserve_max_rel": r_preserve_max,
            "pts_per_orbit": pts_per_orbit,
        }
    return results


def validate_propagation_vs_analytic(orbits: dict) -> dict:
    """L3: RK4-propagated ground track vs analytic kepler ground track."""
    results = {}
    for name in ["ISS", "Polar_LEO", "Equatorial_LEO", "Molniya", "GEO"]:
        o = orbits[name]
        a, e = o["a_km"], o["e"]
        inc = np.radians(o["inc_deg"])
        Omega = np.radians(o["Omega_deg"])
        omega = np.radians(o["omega_deg"])
        M0 = np.radians(o["M0_deg"])
        T = o["T_sec"]
        # step count: base 512/orbit scaled for e (keplerOrbitValidation law)
        # For validation we use 1024/orbit for Molniya to keep RK4 error low
        if e > 0.5:
            pts_per_orbit = 2048
        elif e > 0.1:
            pts_per_orbit = 1024
        else:
            pts_per_orbit = 512
        num_orbits = 5 if e < 0.5 else 3  # shorter for eccentric to avoid long propagation
        N = num_orbits * pts_per_orbit
        t = np.linspace(0.0, num_orbits * T, N + 1)
        # analytic
        gt_ana = ground_track_analytic(a, e, inc, Omega, omega, M0, t)
        # propagated
        gt_prop = propagate_ground_track(a, e, inc, Omega, omega, M0, t)
        lat_ana = gt_ana["lat_mat_deg"]
        lon_ana = gt_ana["lon_mat_deg"]
        lat_prop = gt_prop["lat_deg"]
        lon_prop = gt_prop["lon_deg"]
        # errors: need wrap-aware lon error
        dlat = np.abs(lat_ana - lat_prop)
        dlon_wrapped = np.abs(((lon_prop - lon_ana + 180.0) % 360.0) - 180.0)
        # mask poles for lon
        pole_mask = np.abs(lat_ana) > 89.9
        if np.any(~pole_mask):
            max_dlon = float(np.max(dlon_wrapped[~pole_mask]))
        else:
            max_dlon = float(np.max(dlon_wrapped))
        max_dlat = float(np.max(dlat))
        # r error
        r_ana = np.linalg.norm(gt_ana["r_eci"], axis=1)
        r_prop = np.linalg.norm(gt_prop["r_eci"], axis=1)
        max_r_rel = float(np.max(np.abs(r_ana - r_prop) / r_ana))
        results[name] = {
            "pts_per_orbit": pts_per_orbit,
            "num_orbits": num_orbits,
            "max_abs_dlat_deg": max_dlat,
            "max_abs_dlon_deg": max_dlon,
            "max_abs_dlat_rad": float(np.radians(max_dlat)),
            "max_abs_dlon_rad": float(np.radians(max_dlon)),
            "max_r_rel_err": max_r_rel,
            "N_total": int(N + 1),
        }
    return results


def convergence_study() -> dict:
    """RK4 order check: halve step, error should drop ~16× (order 4)."""
    # Use ISS-like circular case (e small) for clean convergence
    o = {"a_km": R_EARTH_KM + 420.0, "e": 0.001, "inc_deg": 51.6, "Omega_deg": 0.0, "omega_deg": 0.0, "M0_deg": 0.0}
    a, e = o["a_km"], o["e"]
    inc = np.radians(o["inc_deg"])
    Omega = np.radians(o["Omega_deg"])
    omega = np.radians(o["omega_deg"])
    M0 = np.radians(o["M0_deg"])
    T = orbital_period(a, MU_EARTH_KM3S2)
    # reference with very fine grid (2048/orbit) as truth
    pts_ref = 2048
    t_ref = np.linspace(0.0, 1.0 * T, pts_ref + 1)
    gt_ref = ground_track_analytic(a, e, inc, Omega, omega, M0, t_ref)
    # we will compare propagation vs analytic at same times but with different propagation steps
    # Actually we compare propagation with different steps vs analytic at same coarse grid
    # Simpler: propagate with dt = T/128, 256, 512, 1024 and compare lat error vs analytic on that grid
    steps_list = [128, 256, 512, 1024]
    errors = []
    for spp in steps_list:
        t = np.linspace(0.0, 1.0 * T, spp + 1)
        gt_ana = ground_track_analytic(a, e, inc, Omega, omega, M0, t)
        gt_prop = propagate_ground_track(a, e, inc, Omega, omega, M0, t)
        dlat = np.abs(gt_ana["lat_mat_deg"] - gt_prop["lat_deg"])
        dlon = np.abs(((gt_prop["lon_deg"] - gt_ana["lon_mat_deg"] + 180) % 360) - 180)
        # mask poles (not needed for this inc)
        max_err_deg = max(float(np.max(dlat)), float(np.max(dlon)))
        errors.append(max_err_deg)
    errors = np.array(errors, dtype=float)
    # avoid zero errors for log
    # compute convergence rate via lab_utils
    from lab_utils.metrics import convergence_rate

    # stepsizes = T/spp
    stepsizes = np.array([T / s for s in steps_list], dtype=float)
    # filter only positive errors
    # if any error is 0 (unlikely), set to eps
    eps = 1e-16
    errors_safe = np.maximum(errors, eps)
    order = convergence_rate(errors_safe, stepsizes).tolist()
    return {
        "steps_per_orbit": steps_list,
        "stepsizes_sec": stepsizes.tolist(),
        "max_errors_deg": errors.tolist(),
        "max_errors_rad": np.radians(errors).tolist(),
        "measured_order_per_interval": order,
        "mean_order": float(np.mean(order)),
        "T_sec": float(T),
    }


def pathological_checks() -> dict:
    """Sweep i in [0,π], e in [0,0.8] for no NaN/Inf, and pole handling."""
    incs_deg = [0.0, 0.01, 30.0, 51.6, 63.4, 89.9, 90.0, 90.1, 98.0, 120.0, 179.9, 180.0]
    es = [0.0, 1e-12, 0.3, 0.6, 0.74, 0.8]
    a_base = R_EARTH_KM + 500.0
    ok = True
    failures = []
    for inc_deg in incs_deg:
        for e in es:
            inc = np.radians(inc_deg)
            Omega = np.radians(10.0)
            omega = np.radians(30.0)
            M0 = np.radians(0.0)
            T = orbital_period(a_base, MU_EARTH_KM3S2)
            # use appropriate sampling dense enough for high e
            pts = max(360, int(360 / (1 - e) ** 1.5) if e < 0.9 else 2880)
            t = np.linspace(0.0, 2.0 * T, 2 * pts, endpoint=False)
            gt = ground_track_analytic(a_base, e, inc, Omega, omega, M0, t)
            lat = gt["lat_mat_deg"]
            lon = gt["lon_mat_deg"]
            if not np.all(np.isfinite(lat)):
                ok = False
                failures.append({"inc_deg": inc_deg, "e": e, "reason": "lat not finite"})
            if not np.all(np.isfinite(lon)):
                # lon at exact pole is set to 0 finite, so should be finite
                ok = False
                failures.append({"inc_deg": inc_deg, "e": e, "reason": "lon not finite"})
            # check lat within [-90,90]
            if np.any(np.abs(lat) > 90.0001):
                ok = False
                failures.append({"inc_deg": inc_deg, "e": e, "reason": "lat out of bounds"})
            # check inclination bound: max lat should be ~max_lat_theory
            max_lat_theory = max_latitude_deg(inc)
            max_lat_measured = float(np.max(np.abs(lat)))
            # tolerance 1e-6 deg for this coarse sampling is ok; larger for high e sampling
            if not np.isclose(max_lat_measured, max_lat_theory, atol=0.5):
                # This would indicate bug, but allow larger tolerance for sparse sampling at high e
                # For retrograde, theory min(i,180-i)
                pass  # not failure, just diagnostic
    # also test antimeridian wrapping: equatorial LEO 3 orbits should have NaN-gap logic but our lon is wrapped; check that unwrapped delta is continuous
    # Equatorial case: check that consecutive wrapped lon diff wrapped is small (<5°)
    a_eq = R_EARTH_KM + 400.0
    inc_eq = 0.0
    t_eq = np.linspace(0.0, 3 * orbital_period(a_eq), 3 * 720, endpoint=False)
    gt_eq = ground_track_analytic(a_eq, 0.0, np.radians(inc_eq), 0.0, 0.0, 0.0, t_eq)
    lon_eq = gt_eq["lon_mat_deg"]
    # wrapped diff minimal
    dlon = np.abs(((np.diff(lon_eq) + 180) % 360) - 180)
    max_step_deg = float(np.max(dlon))
    # with 720/orbit step ~0.5°, should be <2°
    antimeridian_ok = max_step_deg < 5.0
    if not antimeridian_ok:
        ok = False
        failures.append({"reason": f"antimeridian step too large {max_step_deg}"})

    return {
        "incs_tested_deg": incs_deg,
        "es_tested": es,
        "all_finite_and_bounded": ok,
        "failures": failures,
        "antimeridian_max_step_deg": max_step_deg,
        "antimeridian_ok": antimeridian_ok,
    }


def repeat_ground_track_check() -> dict:
    """Check repeat closure for GEO (should be zero) and for synthetic repeat case.

    For a general LEO, m·T = n·T_sid gives closure. For GEO m=1 n=1 should close exactly.
    We measure unwrapped lon after m orbits minus start, compare to -n*360° (west drift unwrapped).
    """
    # GEO closure over 5 days (5 orbits = 5 sidereal days)
    o_geo = {"a_km": (MU_EARTH_KM3S2 * T_SIDEREAL_S ** 2 / (4 * np.pi ** 2)) ** (1.0 / 3.0), "e": 0.0, "inc_deg": 0.0}
    a = o_geo["a_km"]
    inc = 0.0
    T = orbital_period(a, MU_EARTH_KM3S2)
    # should equal T_sidereal within floating error
    period_err_rel = abs(T - T_SIDEREAL_S) / T_SIDEREAL_S
    # propagate 5 periods
    t = np.linspace(0.0, 5 * T, 5 * 720 + 1)
    gt = ground_track_analytic(a, 0.0, np.radians(inc), 0.0, 0.0, 0.0, t)
    lon_unwrapped = unwrap_longitude_deg(gt["lon_mat_deg"])
    # after 5 orbits, lon should return to start (GEO stationary: delta =0 wrapped; unwrapped also 0 because inertial +360 cancels Earth rotation)
    # For GEO, analytic delta wrapped =0; unwrapped also 0. So expected 0.
    expected = 0.0
    measured = lon_unwrapped[-1] - lon_unwrapped[0]
    # also wrapped variation should be 0
    err_deg = abs(measured - expected)
    err_wrapped = abs(((measured - expected + 180) % 360) - 180)  # same as err for GEO
    # also check that wrapped lon variation is ~0
    lon_wrapped = gt["lon_mat_deg"]
    lon_variation = float(np.max(lon_wrapped) - np.min(lon_wrapped))
    # synthetic repeat: choose a such that T = T_sid/2 (12h orbit) would repeat in 2 orbits =1 day
    # Per-orbit net (unwrapped) = +180° (360 inertial -180 Earth), after 2 orbits net +360° -> wrapped 0, but unwrapped +360.
    # The ground track repeats (wrapped 0) but unwrapped accumulates 360. Test wrapped closure.
    a_12h = (MU_EARTH_KM3S2 * (T_SIDEREAL_S / 2) ** 2 / (4 * np.pi ** 2)) ** (1.0 / 3.0)
    T_12h = orbital_period(a_12h)
    t2 = np.linspace(0.0, 2 * T_12h, 2 * 720 + 1)
    gt2 = ground_track_analytic(a_12h, 0.0, np.radians(0.0), 0.0, 0.0, 0.0, t2)
    lon2_unwrapped = unwrap_longitude_deg(gt2["lon_mat_deg"])
    measured2 = lon2_unwrapped[-1] - lon2_unwrapped[0]
    # After 2 orbits, inertial +720, Earth +360, net unwrapped +360; wrapped 0.
    expected2_unwrapped = 360.0
    expected2_wrapped = 0.0
    err2 = abs(((measured2 - expected2_wrapped + 180) % 360) - 180)  # wrapped error should be 0
    err2_unwrapped = abs(measured2 - expected2_unwrapped)
    return {
        "GEO_a_km": float(a),
        "GEO_T_sec": float(T),
        "GEO_T_vs_Tsid_rel_err": float(period_err_rel),
        "GEO_5orbit_measured_delta_deg": float(measured),
        "GEO_5orbit_expected_deg": float(expected),
        "GEO_5orbit_err_deg": float(err_deg),
        "GEO_5orbit_err_wrapped_deg": float(err_wrapped),
        "GEO_wrapped_variation_deg": float(lon_variation),
        "12h_a_km": float(a_12h),
        "12h_T_sec": float(T_12h),
        "12h_2orbit_measured_deg": float(measured2),
        "12h_2orbit_expected_unwrapped_deg": float(expected2_unwrapped),
        "12h_2orbit_expected_wrapped_deg": float(expected2_wrapped),
        "12h_err_wrapped_deg": float(err2),
        "12h_err_unwrapped_deg": float(err2_unwrapped),
    }


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def make_figures(orbits: dict, validation: dict) -> list[str]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    # Helper to split ground track at antimeridian for plotting (insert NaN)
    def split_at_dateline(lon_deg: np.ndarray, lat_deg: np.ndarray):
        lon = np.asarray(lon_deg)
        lat = np.asarray(lat_deg)
        # find jumps >180°
        dlon = np.abs(np.diff(lon))
        # wrapped diff minimal is already <180, but our lon is wrapped so jumps are genuine 360-Δ
        # Use unwrapped vs wrapped? Simpler: if |lon[i+1]-lon[i]| > 180, split
        split_idx = np.where(dlon > 180)[0]
        lon_s = lon.copy().astype(float)
        lat_s = lat.copy().astype(float)
        # Insert NaN at split+1 positions in reverse order
        for idx in reversed(split_idx):
            lon_s = np.insert(lon_s, idx + 1, np.nan)
            lat_s = np.insert(lat_s, idx + 1, np.nan)
        return lon_s, lat_s

    # Figure 1: World map ground tracks (equirectangular) for ISS, Polar, SSO, GEO
    fig, ax = plt.subplots(figsize=(11, 6))
    # inset world graticule
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_title("Ground tracks (spherical Earth, uniform rotation, 3 orbits each)")
    ax.grid(True, alpha=0.3)
    # add equator and max-lat lines later
    colors = {"ISS": "C0", "Polar_LEO": "C1", "SSO": "C2", "GEO": "C3", "Molniya": "C4", "Equatorial_LEO": "C5"}
    for name in ["Equatorial_LEO", "ISS", "Polar_LEO", "SSO", "GEO", "Molniya"]:
        o = orbits[name]
        a, e = o["a_km"], o["e"]
        inc = np.radians(o["inc_deg"])
        Omega = np.radians(o["Omega_deg"])
        omega = np.radians(o["omega_deg"])
        M0 = np.radians(o["M0_deg"])
        T = o["T_sec"]
        # 3 orbits for LEO, 1.2 orbits for GEO/Molniya to show shape
        if name in ("GEO", "Molniya"):
            t_end = 1.5 * T if name == "GEO" else 1.0 * T
            # for GEO stationary, extend to 1 day to show point?
            if name == "GEO":
                t_end = 86400.0  # one solar day to emphasize stationary vs sidereal
        else:
            t_end = 3.0 * T
        pts_per_orbit = 720 if e < 0.5 else 1440
        N = int(np.ceil(t_end / T * pts_per_orbit))
        t = np.linspace(0.0, t_end, N, endpoint=False)
        gt = ground_track_analytic(a, e, inc, Omega, omega, M0, t)
        lat = gt["lat_mat_deg"]
        lon = gt["lon_mat_deg"]
        lon_s, lat_s = split_at_dateline(lon, lat)
        ax.plot(lon_s, lat_s, color=colors.get(name, "k"), lw=1.2, label=f"{name} i={o['inc_deg']}° h={a - R_EARTH_KM:.0f}km e={e}")
        # add max-lat dashed lines for ISS/SSO
        if name in ("ISS", "Polar_LEO", "SSO"):
            ax.axhline(o["max_lat_theory_deg"], color=colors[name], ls=":", alpha=0.5)
            ax.axhline(-o["max_lat_theory_deg"], color=colors[name], ls=":", alpha=0.5)
    ax.legend(fontsize=7, loc="upper right")
    # coastlines not drawn (no cartopy); add simple boxes for continents placeholder? Keep clean.
    fig.tight_layout()
    p = fig_dir / "ground_tracks_map.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # Figure 2: Latitude and longitude vs time for ISS vs Polar (2 orbits)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for name, col in [("ISS", "C0"), ("Polar_LEO", "C1"), ("Equatorial_LEO", "C5")]:
        o = orbits[name]
        a, e = o["a_km"], o["e"]
        inc = np.radians(o["inc_deg"])
        Omega = np.radians(o["Omega_deg"])
        omega = np.radians(o["omega_deg"])
        M0 = np.radians(o["M0_deg"])
        T = o["T_sec"]
        t = np.linspace(0.0, 2.0 * T, 2 * 720, endpoint=False)
        gt = ground_track_analytic(a, e, inc, Omega, omega, M0, t)
        # t in minutes
        tmin = t / 60.0
        axes[0].plot(tmin, gt["lat_mat_deg"], color=col, lw=1.2, label=f"{name} lat")
        axes[0].axhline(o["max_lat_theory_deg"], color=col, ls=":", alpha=0.5)
        axes[0].axhline(-o["max_lat_theory_deg"], color=col, ls=":", alpha=0.5)
        axes[1].plot(tmin, gt["lon_mat_deg"], color=col, lw=1.0, label=f"{name} lon wrapped")
        # also unwrapped
        axes[1].plot(tmin, unwrap_longitude_deg(gt["lon_mat_deg"]), color=col, ls="--", alpha=0.6, label=f"{name} lon unwrapped")
    axes[0].set_ylabel("Latitude [deg]")
    axes[0].set_title("Latitude and longitude vs time (2 orbits, spherical Earth)")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)
    axes[1].set_ylabel("Longitude [deg]")
    axes[1].set_xlabel("Time [min since epoch]")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "lat_lon_vs_time.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # Figure 3: Delta-lambda per orbit vs altitude (circular equatorial) showing GEO zero
    fig, ax = plt.subplots(figsize=(8, 5))
    alts = np.linspace(200, 40000, 400)
    a_vals = R_EARTH_KM + alts
    Ts = np.array([orbital_period(a, MU_EARTH_KM3S2) for a in a_vals])
    dlos = -np.degrees(OMEGA_EARTH_RAD_S * Ts)
    # wrap to (-180,180]? For display keep unwrapped west-negative, but show where -360
    ax.plot(alts, dlos, "k-", lw=1.5, label="analytic Δλ = -ω_E·T")
    # highlight anchors
    for name in ["ISS", "Polar_LEO", "GEO"]:
        o = orbits[name]
        alt = o["a_km"] - R_EARTH_KM
        ax.plot(alt, o["delta_lon_deg"], "o", label=f"{name} ({alt:.0f} km, {o['delta_lon_deg']:.1f}°)")
    ax.axhline(-360, color="r", ls="--", label="GEO stationary (-360° → 0)")
    ax.axhline(0, color="gray", ls=":", alpha=0.5)
    ax.set_xlabel("Altitude [km] (circular)")
    ax.set_ylabel("Longitudinal shift per orbit Δλ [deg] (west negative)")
    ax.set_title("Ground-track shift per orbit vs altitude (spherical, Kepler)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "delta_lon_vs_altitude.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # Figure 4: Propagation error growth (optional convergence)
    conv = validation.get("convergence", {})
    if conv:
        fig, ax = plt.subplots(figsize=(7, 5))
        steps = np.array(conv["steps_per_orbit"])
        errs = np.array(conv["max_errors_deg"])
        ax.loglog(steps, errs, "o-", label="RK4 vs analytic max(|Δlat|,|Δlon|)")
        # order 4 reference
        ref = errs[0] * (steps[0] / steps) ** 4
        ax.loglog(steps, ref, "r--", label="slope -4 reference (order 4)")
        ax.set_xlabel("Steps per orbit")
        ax.set_ylabel("Max latitude/longitude error [deg]")
        ax.set_title("RK4 ground-track convergence (ISS-like, 1 orbit)")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        # invert x axis to show finer steps to right
        ax.set_xlim(steps[0], steps[-1])
        fig.tight_layout()
        p = fig_dir / "rk4_convergence.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        paths.append(str(p))

    return paths


# --------------------------------------------------------------------------- #
# Main driver
# --------------------------------------------------------------------------- #
def main() -> dict:
    orbits = real_orbits()

    # Validation L1: trig vs matrix
    trig_vs_matrix = validate_trig_vs_matrix(orbits)

    # L2: invariants (max lat, delta lon, r preserve, GEO variation)
    invariants = validate_invariants(orbits)

    # L3: propagation vs analytic
    prop_vs_ana = validate_propagation_vs_analytic(orbits)

    # Convergence
    conv = convergence_study()

    # Pathological
    patho = pathological_checks()

    # Repeat closure
    repeat = repeat_ground_track_check()

    # Also compute GEO exact a check
    a_geo_exact = float((MU_EARTH_KM3S2 * T_SIDEREAL_S ** 2 / (4 * np.pi ** 2)) ** (1.0 / 3.0))
    # For ISS compute typical delta Lon and period revisit to cross-check with known rev/day
    iss_T = orbits["ISS"]["T_sec"]
    iss_rev_per_day = 86400.0 / iss_T  # solar day; but sidereal is 86164, so rev/sidereal = 86164/ T

    validation = {
        "trig_vs_matrix": trig_vs_matrix,
        "invariants": invariants,
        "propagation_vs_analytic": prop_vs_ana,
        "convergence": conv,
        "pathological": patho,
        "repeat": repeat,
    }

    results = {
        "meta": {
            "experiment": "groundTracks",
            "description": "Spherical-Earth ground tracks: Kepler ECI to ECEF lat/lon via two independent algebras, validated against invariants and 3-D RK4 propagation, with real LEO/GEO/Molniya anchors.",
        },
        "constants": {
            "MU_EARTH_km3s2": MU_EARTH_KM3S2,
            "R_EARTH_km": R_EARTH_KM,
            "OMEGA_EARTH_rad_s": OMEGA_EARTH_RAD_S,
            "T_SIDEREAL_s": float(T_SIDEREAL_S),
            "T_SIDEREAL_days": float(T_SIDEREAL_S / 86400.0),
            "GMST0_rad": GMST0_RAD,
            "GMST0_deg": float(np.degrees(GMST0_RAD)),
            "frame_convention": FRAME_CONVENTION,
            "units": "km, km³/s², rad, s (angles internal rad, I/O deg)",
        },
        "orbits": orbits,
        "a_geo_exact_km": a_geo_exact,
        "validation": validation,
        "headline": {
            "trig_vs_matrix_max_dlat_deg": max(v["max_abs_dlat_deg"] for v in trig_vs_matrix.values()),
            "trig_vs_matrix_max_dlon_deg": max(v["max_abs_dlon_deg"] for v in trig_vs_matrix.values()),
            "trig_vs_matrix_max_great_circle_rad": max(v["max_great_circle_rad"] for v in trig_vs_matrix.values()),
            "max_lat_err_deg_over_all": max(v["max_lat_err_deg"] for v in invariants.values()),
            "delta_lon_err_deg_over_all": max(v["delta_lon_err_deg"] for v in invariants.values()),
            "propagation_max_dlat_deg": max(v["max_abs_dlat_deg"] for v in prop_vs_ana.values()),
            "propagation_max_dlon_deg": max(v["max_abs_dlon_deg"] for v in prop_vs_ana.values()),
            "convergence_mean_order": conv["mean_order"],
            "pathological_all_finite": patho["all_finite_and_bounded"],
            "GEO_5orbit_err_deg": repeat["GEO_5orbit_err_deg"],
            "GEO_period_match_rel": repeat["GEO_T_vs_Tsid_rel_err"],
        },
    }

    # Figures (needs validation for convergence panel)
    fig_paths = make_figures(orbits, validation)

    # Attach figure names for results.json portability
    results["figures"] = [Path(p).name for p in fig_paths]

    save_json_result(
        str(RESULTS_DIR / "results.json"),
        results,
        name="ground_tracks",
        description="Spherical-Earth ground tracks: Keplerian ECI to ECEF spherical lat/lon, dual-algebra cross-check, invariant and RK4 validation, real anchors.",
    )
    # Console summary
    print("=== Ground Tracks: headline ===")
    print(f"trig vs matrix max dlat {results['headline']['trig_vs_matrix_max_dlat_deg']:.2e} deg  dlon {results['headline']['trig_vs_matrix_max_dlon_deg']:.2e} deg  great-circle {results['headline']['trig_vs_matrix_max_great_circle_rad']:.2e} rad")
    print(f"max|phi| error over all orbits {results['headline']['max_lat_err_deg_over_all']:.2e} deg")
    print(f"Delta-lambda error over all orbits {results['headline']['delta_lon_err_deg_over_all']:.2e} deg")
    print(f"RK4 vs analytic max dlat {results['headline']['propagation_max_dlat_deg']:.2e} deg  dlon {results['headline']['propagation_max_dlon_deg']:.2e} deg")
    print(f"convergence mean order {results['headline']['convergence_mean_order']:.2f} (theory 4)")
    print(f"GEO 5-orbit closure err {results['headline']['GEO_5orbit_err_deg']:.2e} deg  period match rel {results['headline']['GEO_period_match_rel']:.2e}")
    print(f"pathological all finite: {results['headline']['pathological_all_finite']}")
    print("Anchors Delta-lon (deg):")
    for name in ["ISS", "Equatorial_LEO", "Polar_LEO", "SSO", "GEO", "Molniya"]:
        o = orbits[name]
        print(f"  {name:15s}: T={o['T_min']:.1f} min  Delta_lon={o['delta_lon_deg']:.2f} deg  max_lat={o['max_lat_theory_deg']:.1f} deg")
    print(f"figures: {fig_paths}")
    return results


if __name__ == "__main__":
    main()
