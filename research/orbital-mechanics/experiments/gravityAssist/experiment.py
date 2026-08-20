"""Planetary gravity assist: the two-body hyperbolic flyby (patched conic).

An unpowered planetary flyby conserves the planet-frame hyperbolic-excess
SPEED (|v_inf,out| = |v_inf,in|) and rotates only its DIRECTION by the turn
angle delta.  Assembled back into the heliocentric frame,

    v_in  = V_p + v_inf,in
    v_out = V_p + v_inf,out
    dv_helio = v_out - v_in                      (vector change, NOT a burn)
    d_epsilon = 0.5 (v_out^2 - v_in^2) = V_p . dv_helio

the spacecraft exchanges heliocentric orbital energy with the moving planet.

Numerical contract (x = r_p v_inf^2 / mu_p; cancellation-safe forms):

    e = 1 + x
    delta = 2 atan2(1, sqrt(x(x+2)))             ==  2 arcsin(1/e)
    pi - delta = 2 atan2(sqrt(x(x+2)), 1)        (near-parabolic form)
    b = r_p sqrt(1 + 2/x)                        (impact parameter; the
      simple form 2 atan(mu/(b v_inf^2)) uses b, NEVER r_p)
    v_inf,out = R_nhat(delta) v_inf,in           (Rodrigues, B-plane normal)

Orientation landscape (alpha = angle(v_inf,in, V_p), phi = flyby-plane phase,
h = delta/2):

    d_eps(alpha, phi) = 2 V_p v_inf sin(h) [ -sin(h) cos(alpha)
                                            + cos(h) sin(alpha) cos(phi) ]
    global max at alpha* = pi/2 + h, phi* = 0, value 2 V_p v_inf sin(h)
      (bend vector EXACTLY parallel to +V_p by the Cauchy-Schwarz equality)
    global min at alpha = pi/2 - h, phi = pi (antiparallel bend)

Canonical anchors are reconstructed from NASA/JPL published hyperbolic
encounter elements (a, e) with planet-only GMs, and the heliocentric energy
change is checked against the published pre/post-encounter semimajor axes.

The independent numerical validation (L3) integrates the raw two-body EOM with
the verified Experiment 006 3D Cowell RK4 from an exact inbound conic state,
then RECOVERS the hyperbola from the final state (energy, angular momentum,
eccentricity vector) and computes the asymptotic directions analytically.
Finite-radius velocity directions are never compared directly: at 5-10 SOI
they misstate the turn angle by 1e-3..1e-1 rad, while element recovery is
patch-radius-insensitive.

References: Curtis Ch. 8 (2021); Bate/Mueller/White Ch. 1 (1971); JPL Pub
82-43 (B-plane); NASA/JPL Voyager encounter elements; JPL SSD astrodynamic
tables (DE440 solar GM; planet-only Jupiter GM from JUP365); NASA Trajectory
Browser (patched-two-body model and its documented limitations).
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

from lab_utils.results import save_json_result

# --- Reuse of the verified 006 machinery (3D Cowell RK4) --------------------
#
# The laboratory rule is "reuse verified code, never rebuild scaffolding".
# Experiment 006 owns a 3D fixed-step RK4 two-body integrator validated to
# <=1e-11 there; it is loaded here by explicit path (same importlib pattern
# the other experiments use).

_pcm_path = (
    Path(__file__).resolve().parents[1] / "planeChangeManeuvers" / "experiment.py"
)
_pcm_spec = importlib.util.spec_from_file_location("plane_change_exp006", _pcm_path)
assert _pcm_spec is not None and _pcm_spec.loader is not None
pcm = importlib.util.module_from_spec(_pcm_spec)
_pcm_spec.loader.exec_module(pcm)

propagate_3d_rk4 = pcm.propagate_3d_rk4

# --- Physical constants ------------------------------------------------------
#
# Solar GM from the DE440-based JPL astrodynamic-parameter table.  Planet GMs
# are PLANET-ONLY values (Jupiter 126686531.9 from JUP365; Saturn
# 37931206.23 from the satellite-ephemeris solution): the major moons must
# not be silently folded into the central body of a planet-centered flyby.
# The DE440 *system* GMs (Jupiter 126712764.1) are deliberately NOT used.
# Radii are JPL equatorial values; a_p are JPL mean heliocentric distances.

MU_SUN_KM3S2 = 132712440041.279  # DE440

PLANETS: dict[str, dict[str, float]] = {
    "Earth": {
        "mu": 398600.435507,
        "R_eq": 6378.1366,
        "a_p": 1.495978707e8,
    },
    "Venus": {
        "mu": 324858.592,
        "R_eq": 6051.8,  # Venus is nearly spherical; mean = equatorial
        "a_p": 1.0821e8,
    },
    "Mars": {
        "mu": 42828.375816,
        "R_eq": 3396.19,
        "a_p": 2.2794e8,
    },
    "Jupiter": {
        "mu": 126686531.9,  # planet-only (JUP365), NOT the DE440 system GM
        "R_eq": 71492.0,
        "a_p": 7.7848e8,
    },
    "Saturn": {
        "mu": 37931206.23,  # planet-only, NOT the DE440 system GM
        "R_eq": 60268.0,
        "a_p": 1.4337e9,
    },
}

K0 = np.array([0.0, 0.0, 1.0])  # reference pole: J2000 ecliptic north

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

# Deterministic sweep grids
V_INF_KMS = np.logspace(np.log10(0.5), np.log10(30.0), 32)
R_P_FACTOR = np.logspace(np.log10(1.02), np.log10(50.0), 64)
ORIENTATION_ALPHA_DEG = np.arange(0.0, 181.0, 1.0)
ORIENTATION_PHI_DEG = np.arange(0.0, 360.0, 1.0)


# --------------------------------------------------------------------------- #
# Closed-form core (the numerical contract)
# --------------------------------------------------------------------------- #
def x_parameter(r_p_km: float, v_inf_kms: float, mu_p_km3s2: float) -> float:
    """x = r_p v_inf^2 / mu_p (keep x directly; never rebuild e-1 from e)."""
    return r_p_km * v_inf_kms**2 / mu_p_km3s2


def eccentricity(x: float) -> float:
    """e = 1 + x."""
    return 1.0 + x


def turn_angle_rad(x: float) -> float:
    """Cancellation-safe turn angle delta = 2 atan2(1, sqrt(x(x+2)))."""
    return 2.0 * np.arctan2(1.0, np.sqrt(x * (x + 2.0)))


def turn_angle_naive_rad(x: float) -> float:
    """Reference form delta = 2 arcsin(1/(1+x)) (loses precision for x << 1)."""
    return 2.0 * np.arcsin(1.0 / eccentricity(x))


def turn_angle_near_pi_rad(x: float) -> float:
    """Near-parabolic complement: pi - delta = 2 atan2(sqrt(x(x+2)), 1)."""
    return 2.0 * np.arctan2(np.sqrt(x * (x + 2.0)), 1.0)


def impact_parameter_km(r_p_km: float, x: float) -> float:
    """b = r_p sqrt(1 + 2/x)."""
    return r_p_km * np.sqrt(1.0 + 2.0 / x)


def periapsis_speed_kms(r_p_km: float, v_inf_kms: float, mu_p_km3s2: float) -> float:
    """v_p = sqrt(v_inf^2 + 2 mu / r_p) (vis-viva at periapsis)."""
    return np.sqrt(v_inf_kms**2 + 2.0 * mu_p_km3s2 / r_p_km)


def planet_speed_kms(planet: dict[str, float]) -> float:
    """Circular heliocentric speed V_p = sqrt(mu_sun / a_p) (idealization)."""
    return np.sqrt(MU_SUN_KM3S2 / planet["a_p"])


def soi_radius_km(planet: dict[str, float]) -> float:
    """Laplace sphere of influence r_SOI = a_p (mu_p / mu_sun)^(2/5)."""
    return planet["a_p"] * (planet["mu"] / MU_SUN_KM3S2) ** 0.4


# --------------------------------------------------------------------------- #
# B-plane geometry and the closed-form flyby
# --------------------------------------------------------------------------- #
def bplane_basis(s_hat: np.ndarray, k0: np.ndarray = K0) -> tuple[np.ndarray, np.ndarray]:
    """Right-handed B-plane basis (t_hat, r_hat) perpendicular to s_hat.

    t_hat = (k0 x s)/|k0 x s|, r_hat = s x t_hat.  beta = 0 means +t_hat.
    """
    t = np.cross(k0, s_hat)
    nrm = np.linalg.norm(t)
    if nrm < 1e-12:  # s_hat parallel to the pole: pick any orthogonal axis
        fallback = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(fallback, s_hat)) > 0.9:
            fallback = np.array([0.0, 1.0, 0.0])
        t = np.cross(fallback, s_hat)
        nrm = np.linalg.norm(t)
    t_hat = t / nrm
    r_hat = np.cross(s_hat, t_hat)
    return t_hat, r_hat


def beta_for_bend(q_hat_desired: np.ndarray, s_hat: np.ndarray) -> float:
    """B-plane phase beta realizing a desired unit bend direction q (q . s = 0).

    q = -(cos(beta) t + sin(beta) r)  =>  cos(beta) = -q.t, sin(beta) = -q.r.
    """
    t_hat, r_hat = bplane_basis(s_hat)
    return float(np.arctan2(-np.dot(q_hat_desired, r_hat), -np.dot(q_hat_desired, t_hat)))


def flyby_v_inf_out(
    v_inf_in_kms: np.ndarray, delta_rad: float, beta_rad: float, b_km: float
) -> np.ndarray:
    """Outbound excess-velocity vector via the B-plane construction.

    B = b(cos(beta) t + sin(beta) r); bend direction q = -B/b (gravity pulls
    toward the planet); n = s x q; v_out = R_n(delta) v_in (Rodrigues).
    """
    v_inf = float(np.linalg.norm(v_inf_in_kms))
    s_hat = v_inf_in_kms / v_inf
    t_hat, r_hat = bplane_basis(s_hat)
    B_vec = b_km * (np.cos(beta_rad) * t_hat + np.sin(beta_rad) * r_hat)
    q_hat = -B_vec / b_km
    n_hat = np.cross(s_hat, q_hat)
    return (
        v_inf_in_kms * np.cos(delta_rad)
        + np.cross(n_hat, v_inf_in_kms) * np.sin(delta_rad)
    )


def assemble_heliocentric(
    planet: dict[str, float], v_inf_in_kms: np.ndarray, v_inf_out_kms: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """v_in = V_p + v_inf,in ; v_out = V_p + v_inf,out (V_p along +z here)."""
    V_p = planet_speed_kms(planet)
    V_vec = np.array([0.0, 0.0, V_p])
    return V_vec + v_inf_in_kms, V_vec + v_inf_out_kms


# --------------------------------------------------------------------------- #
# Orientation landscape (closed form + grid cross-check)
# --------------------------------------------------------------------------- #
def orientation_factor(alpha_rad: np.ndarray, phi_rad: np.ndarray, h_rad: float):
    """F(alpha, phi) = -sin(h) cos(alpha) + cos(h) sin(alpha) cos(phi).

    d_eps = A * F with A = 2 V_p v_inf sin(h).  Vectorized over meshes.
    """
    return -np.sin(h_rad) * np.cos(alpha_rad) + np.cos(h_rad) * np.sin(
        alpha_rad
    ) * np.cos(phi_rad)


def analytic_optimum(h_rad: float) -> dict[str, float]:
    """Exact global extrema of the orientation landscape.

    Max at alpha* = pi/2 + h (phi = 0), min at alpha = pi/2 - h (phi = pi);
    the extremal bend vector is exactly parallel/antiparallel to V_p.
    """
    return {
        "alpha_max_rad": np.pi / 2.0 + h_rad,
        "phi_max_rad": 0.0,
        "alpha_min_rad": np.pi / 2.0 - h_rad,
        "phi_min_rad": np.pi,
        "F_max": 1.0,
        "F_min": -1.0,
    }


def orientation_grid_max(h_rad: float) -> tuple[float, float, float]:
    """Grid cross-check: max of F on the 1 deg x 1 deg orientation grid.

    Returns (F_grid_max, alpha_at_max_rad, phi_at_max_rad).
    """
    am = np.radians(ORIENTATION_ALPHA_DEG)
    pm = np.radians(ORIENTATION_PHI_DEG)
    amg, pmg = np.meshgrid(am, pm, indexing="ij")
    F = orientation_factor(amg, pmg, h_rad)
    idx = np.unravel_index(np.argmax(F), F.shape)
    return float(F[idx]), float(am[idx[0]]), float(pm[idx[1]])


# --------------------------------------------------------------------------- #
# Canonical anchors (published Voyager encounter elements)
# --------------------------------------------------------------------------- #
# NASA/JPL published hyperbolic encounter elements (NASA Science Voyager
# planetary-elements table, provenance attributed to JPL trajectory
# engineering) and heliocentric pre/post semimajor axes (km).
ANCHORS = [
    {
        "name": "Voyager 1 - Jupiter 1979",
        "planet": "Jupiter",
        "a_enc": -1092356.0,
        "e_enc": 1.318976,
        "a_helio_in": 745761000.0,
        "a_helio_out": -593237000.0,
        "published_ca_km": 348890.0,  # NASA 1979 mission account
        "published_ca_rounded_km": 350000.0,  # JPL DESCANSO history
    },
    {
        "name": "Voyager 2 - Jupiter 1979",
        "planet": "Jupiter",
        "a_enc": -2184140.0,
        "e_enc": 1.330279,
        "a_helio_in": 544470000.0,
        "a_helio_out": -2220315000.0,
        "published_ca_km": 730000.0,  # JPL historical (rounded)
    },
    {
        "name": "Voyager 1 - Saturn 1980",
        "planet": "Saturn",
        "a_enc": -166152.0,
        "e_enc": 2.107561,
        "a_helio_in": -593237000.0,
        "a_helio_out": -480926000.0,
        "published_ca_km": None,
    },
]


def reconstruct_anchor(anchor: dict) -> dict:
    """Derive the flyby quantities from the published (a, e) elements.

    v_inf = sqrt(mu_p/|a|); r_p = |a|(e-1); delta = 2 arcsin(1/e);
    |dv_helio| = 2 v_inf sin(delta/2);
    d_eps = mu_sun/2 (1/a_in - 1/a_out);
    scalar heliocentric speed change at the planet's distance (circular
    approximation for the encounter radius r_enc ~ a_p).
    """
    p = PLANETS[anchor["planet"]]
    mu_p = p["mu"]
    a_abs = abs(anchor["a_enc"])
    e = anchor["e_enc"]
    v_inf = float(np.sqrt(mu_p / a_abs))
    r_p = a_abs * (e - 1.0)
    delta = 2.0 * float(np.arcsin(1.0 / e))
    dv_vec_mag = 2.0 * v_inf * np.sin(delta / 2.0)
    d_eps = 0.5 * MU_SUN_KM3S2 * (
        1.0 / anchor["a_helio_in"] - 1.0 / anchor["a_helio_out"]
    )
    r_enc = p["a_p"]  # circular idealization of the encounter distance
    v_in = float(np.sqrt(MU_SUN_KM3S2 * (2.0 / r_enc - 1.0 / anchor["a_helio_in"])))
    v_out = float(np.sqrt(MU_SUN_KM3S2 * (2.0 / r_enc - 1.0 / anchor["a_helio_out"])))
    out = {
        "v_inf_kms": v_inf,
        "r_p_km": float(r_p),
        "delta_deg": float(np.degrees(delta)),
        "dv_vector_kms": float(dv_vec_mag),
        "d_epsilon_km2s2": float(d_eps),
        "v_helio_in_kms": v_in,
        "v_helio_out_kms": v_out,
        "dv_scalar_kms": v_out - v_in,
        "x": x_parameter(r_p, v_inf, mu_p),
        "soi_km": soi_radius_km(p),
        "r_p_over_R_eq": float(r_p / p["R_eq"]),
    }
    if anchor.get("published_ca_km") is not None:
        out["ca_spread_pct"] = 100.0 * (anchor["published_ca_km"] - r_p) / r_p
    return out


# --------------------------------------------------------------------------- #
# L3: independent hyperbolic propagation + element recovery
# --------------------------------------------------------------------------- #
def hyperbolic_state_at_radius(
    r_p_km: float,
    v_inf_kms: float,
    mu_p_km3s2: float,
    R_km: float,
    inbound: bool,
    e_hat: np.ndarray,
    h_hat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Exact conic state (r, v) at radius R on the inbound/outbound leg.

    Built in the flyby plane (e_hat = periapsis direction, h_hat = normal).
    This is SETUP for the numerical test, not the thing under test.
    """
    e = eccentricity(x_parameter(r_p_km, v_inf_kms, mu_p_km3s2))
    p = (mu_p_km3s2 / v_inf_kms**2) * (e**2 - 1.0)  # = h^2/mu
    cos_nu = np.clip((p / R_km - 1.0) / e, -1.0, 1.0)
    nu = np.arccos(cos_nu)
    if inbound:
        nu = -nu
    q_hat = np.cross(h_hat, e_hat)
    r_vec = (p / (1.0 + e * np.cos(nu))) * (
        np.cos(nu) * e_hat + np.sin(nu) * q_hat
    )
    v_vec = np.sqrt(mu_p_km3s2 / p) * (
        -np.sin(nu) * e_hat + (e + np.cos(nu)) * q_hat
    )
    return r_vec, v_vec


def recover_hyperbola(
    r_vec: np.ndarray, v_vec: np.ndarray, mu_p_km3s2: float
) -> dict[str, float]:
    """Recover the hyperbola from a state: eps, v_inf, h, e, asymptotes, delta.

    The asymptotic directions follow from the recovered orbital elements, so
    the result carries NO finite-patch truncation (only integrator error).
    """
    r = float(np.linalg.norm(r_vec))
    v = float(np.linalg.norm(v_vec))
    eps = 0.5 * v**2 - mu_p_km3s2 / r
    v_inf_rec = float(np.sqrt(2.0 * eps))
    h_vec = np.cross(r_vec, v_vec)
    h = float(np.linalg.norm(h_vec))
    e_vec = np.cross(v_vec, h_vec) / mu_p_km3s2 - r_vec / r
    e_rec = float(np.linalg.norm(e_vec))
    e_hat = e_vec / e_rec
    h_hat = h_vec / h
    q_hat = np.cross(h_hat, e_hat)
    s = np.sqrt(e_rec**2 - 1.0)
    v_inf_in_dir = (e_hat + s * q_hat) / e_rec
    v_inf_out_dir = (-e_hat + s * q_hat) / e_rec
    delta_rec = float(np.arccos(np.clip(np.dot(v_inf_in_dir, v_inf_out_dir), -1, 1)))
    r_p_rec = h**2 / (mu_p_km3s2 * (1.0 + e_rec))
    return {
        "eps": eps,
        "v_inf_kms": v_inf_rec,
        "e": e_rec,
        "r_p_km": float(r_p_rec),
        "delta_rad": delta_rec,
        "h": h,
        "v_inf_in_dir": v_inf_in_dir,
        "v_inf_out_dir": v_inf_out_dir,
    }


def propagate_flyby_l3(
    r_p_km: float,
    v_inf_kms: float,
    mu_p_km3s2: float,
    R0_factor: float = 100.0,
    theta: float = 0.005,
) -> dict[str, float]:
    """Independent L3 validation of one flyby.

    Starts from the exact conic state at R0 = R0_factor * r_p on the inbound
    leg, integrates the raw two-body EOM through periapsis with the verified
    006 3D Cowell RK4 using the deterministic angular step rule
    dt = theta * r / v, and recovers the hyperbola from the outbound state.
    """
    e = eccentricity(x_parameter(r_p_km, v_inf_kms, mu_p_km3s2))
    # flyby plane: periapsis along +x, normal along +z (any deterministic frame)
    e_hat = np.array([1.0, 0.0, 0.0])
    h_hat = np.array([0.0, 0.0, 1.0])
    R0 = R0_factor * r_p_km
    r0, v0 = hyperbolic_state_at_radius(
        r_p_km, v_inf_kms, mu_p_km3s2, R0, inbound=True, e_hat=e_hat, h_hat=h_hat
    )
    r = r0.copy()
    v = v0.copy()
    n_steps = 0
    max_steps = 200_000
    while True:
        rr = float(np.linalg.norm(r))
        vv = float(np.linalg.norm(v))
        outbound = bool(np.dot(r, v) > 0.0)
        if outbound and rr >= R0:
            break
        if n_steps >= max_steps:
            raise RuntimeError("L3 propagation did not exit the flyby")
        dt = theta * rr / vv
        seg = propagate_3d_rk4(r, v, mu_p_km3s2, np.array([0.0, dt]), dt)
        r, v = seg[-1, :3], seg[-1, 3:]
        n_steps += 1
    rec = recover_hyperbola(r, v, mu_p_km3s2)
    rec["n_steps"] = n_steps
    rec["R0_factor"] = R0_factor
    rec["eps_planet_conservation_rel"] = abs(
        (0.5 * float(np.linalg.norm(v)) ** 2 - mu_p_km3s2 / float(np.linalg.norm(r)))
        - rec["eps"]
    ) / abs(rec["eps"])
    return rec


# --------------------------------------------------------------------------- #
# Pathological regimes
# --------------------------------------------------------------------------- #
def pathological_checks() -> dict:
    """Stable-vs-naive delta, monotonicity, and mpmath verification at extremes."""
    xs = np.logspace(-18, 12, 61)
    d_stable = np.array([turn_angle_rad(x) for x in xs])
    d_naive = np.array([turn_angle_naive_rad(x) for x in xs])
    d_near_pi = np.array([turn_angle_near_pi_rad(x) for x in xs])
    finite = bool(np.all(np.isfinite(d_stable)))
    monotone = bool(np.all(np.diff(d_stable) < 0.0))
    rel_err = float(np.max(np.abs(d_stable - d_naive)[xs > 1e-8] / d_naive[xs > 1e-8]))
    near_pi_err = float(np.max(np.abs((np.pi - d_stable) - d_near_pi)[xs < 1e-2]))
    mp.mp.dps = 50
    worst_mp = 0.0
    for x in (1e-18, 1e-12, 1e-6, 1.0, 1e6, 1e12):
        xm = mp.mpf(x)
        d_mp = 2 * mp.atan2(mp.mpf(1), mp.sqrt(xm * (xm + 2)))
        d_fl = turn_angle_rad(x)
        worst_mp = max(worst_mp, float(abs(mp.mpf(d_fl) - d_mp) / d_mp))
    return {
        "x_grid_log10min": float(np.log10(xs[0])),
        "x_grid_log10max": float(np.log10(xs[-1])),
        "all_finite": finite,
        "delta_monotone_decreasing_in_x": monotone,
        "stable_vs_naive_max_rel_err_above_1e-8": rel_err,
        "near_pi_form_max_abs_err_below_x_1e-2": near_pi_err,
        "stable_vs_mpmath_max_rel_err": worst_mp,
        "naive_form_fails_below_x": float(xs[np.argmax(np.abs(d_stable - d_naive) > 1e-12 * d_naive)])
        if np.any(np.abs(d_stable - d_naive) > 1e-12 * d_naive)
        else None,
    }


# --------------------------------------------------------------------------- #
# Experiment driver
# --------------------------------------------------------------------------- #
def run_sweep() -> dict:
    """Analytic ceiling |d_eps|_max over planets x v_inf x r_p (10,240 cases)."""
    sweep = {}
    for name, p in PLANETS.items():
        V_p = planet_speed_kms(p)
        vv, rr_factor = np.meshgrid(V_INF_KMS, R_P_FACTOR, indexing="ij")
        rr = rr_factor * p["R_eq"]  # factor -> km for THIS planet
        x = x_parameter(rr, vv, p["mu"])
        h = 0.5 * turn_angle_rad(x)
        d_eps_max = 2.0 * V_p * vv * np.sin(h)
        dv_vec = 2.0 * vv * np.sin(h)
        sweep[name] = {
            "V_p_kms": float(V_p),
            "soi_km": float(soi_radius_km(p)),
            "d_eps_max_mean_km2s2": float(np.mean(d_eps_max)),
            "d_eps_max_grid_km2s2": d_eps_max,
            "dv_vector_grid_kms": dv_vec,
            "x_grid": x,
        }
    return sweep


def representative_cases() -> list[dict]:
    """The L3 validation matrix: moderate, strong, weak, and near-parabolic."""
    cases = []
    for name, rpf, vinf in [
        ("Earth", 1.10, 5.0),
        ("Earth", 5.0, 2.0),
        ("Jupiter", 1.02, 10.0),
        ("Jupiter", 1.05, 0.7),  # near-parabolic (x ~ 2.8e-4)
        ("Saturn", 1.10, 15.0),
        ("Mars", 1.30, 3.0),
    ]:
        p = PLANETS[name]
        r_p = rpf * p["R_eq"]
        x = x_parameter(r_p, vinf, p["mu"])
        cases.append(
            {
                "planet": name,
                "r_p_km": r_p,
                "r_p_over_R_eq": rpf,
                "v_inf_kms": vinf,
                "x": x,
                "delta_analytic_rad": turn_angle_rad(x),
            }
        )
    return cases


def main() -> dict:
    results: dict = {"meta": {"experiment": "gravityAssist"}}

    # 1. constants and SOI table
    results["constants"] = {
        name: {
            "mu_planet_only_km3s2": p["mu"],
            "R_eq_km": p["R_eq"],
            "a_p_km": p["a_p"],
            "V_p_circular_kms": float(planet_speed_kms(p)),
            "soi_laplace_km": float(soi_radius_km(p)),
        }
        for name, p in PLANETS.items()
    }

    # 2. orientation landscape: analytic optimum vs 1 deg grid (representative)
    ori = {}
    for name, rpf, vinf in [("Jupiter", 1.02, 10.0), ("Earth", 1.10, 3.0)]:
        p = PLANETS[name]
        r_p = rpf * p["R_eq"]
        x = x_parameter(r_p, vinf, p["mu"])
        h = 0.5 * turn_angle_rad(x)
        V_p = planet_speed_kms(p)
        A = 2.0 * V_p * vinf * np.sin(h)
        opt = analytic_optimum(h)
        f_grid, a_grid, p_grid = orientation_grid_max(h)
        # closed-form check: at the analytic optimum (phi = 0, i.e. the bend
        # lies in the (V_p, s) plane) the bend vector is parallel to V_p.
        # Construct the outbound vector through the B-plane machinery by
        # solving for the beta that realizes the desired bend direction.
        alpha_star = opt["alpha_max_rad"]
        k_hat = np.array([0.0, 0.0, 1.0])
        s_hat = np.array([np.sin(alpha_star), 0.0, np.cos(alpha_star)])
        q_desired = (k_hat - np.cos(alpha_star) * s_hat) / np.sin(alpha_star)
        beta_star = beta_for_bend(q_desired, s_hat)
        v_in = vinf * s_hat
        v_out = flyby_v_inf_out(
            v_in, 2.0 * h, beta_star, impact_parameter_km(r_p, x)
        )
        dv = v_out - v_in
        dv_parallel_err = float(
            np.linalg.norm(dv - (np.dot(dv, k_hat)) * k_hat)
        ) / float(np.linalg.norm(dv))
        vi, vo = assemble_heliocentric(p, v_in, v_out)
        d_eps_direct = 0.5 * (float(np.dot(vo, vo)) - float(np.dot(vi, vi)))
        ori[f"{name}_rp{rpf}_v{vinf}"] = {
            "delta_deg": float(np.degrees(2.0 * h)),
            "alpha_star_deg": float(np.degrees(alpha_star)),
            "A_km2s2": float(A),
            "F_grid_max": f_grid,
            "F_grid_max_alpha_deg": float(np.degrees(a_grid)),
            "F_grid_max_phi_deg": float(np.degrees(p_grid)),
            "grid_loss_pct": 100.0 * (1.0 - f_grid),
            "dv_parallel_err_at_optimum": dv_parallel_err,
            "d_eps_direct_at_optimum": d_eps_direct,
            "d_eps_direct_vs_A_rel_err": abs(d_eps_direct - A) / A,
        }
    results["orientation_landscape"] = ori

    # 3. anchors
    results["anchors"] = {a["name"]: reconstruct_anchor(a) for a in ANCHORS}

    # 4. L3 independent propagation
    l3 = {}
    for case in representative_cases():
        p = PLANETS[case["planet"]]
        rec = propagate_flyby_l3(case["r_p_km"], case["v_inf_kms"], p["mu"])
        d_an = case["delta_analytic_rad"]
        l3[f"{case['planet']}_rp{case['r_p_over_R_eq']}_v{case['v_inf_kms']}"] = {
            "delta_analytic_deg": float(np.degrees(d_an)),
            "delta_recovered_deg": float(np.degrees(rec["delta_rad"])),
            "delta_abs_err_rad": abs(rec["delta_rad"] - d_an),
            "delta_rel_err": abs(rec["delta_rad"] - d_an) / d_an,
            "r_p_rel_err": abs(rec["r_p_km"] - case["r_p_km"]) / case["r_p_km"],
            "v_inf_rel_err": abs(rec["v_inf_kms"] - case["v_inf_kms"]) / case["v_inf_kms"],
            "n_steps": rec["n_steps"],
            "R0_factor": rec["R0_factor"],
        }
    # patch-radius insensitivity on one case
    p = PLANETS["Jupiter"]
    conv = {}
    for r0f in (50.0, 100.0, 200.0):
        rec = propagate_flyby_l3(
            1.02 * p["R_eq"], 10.0, p["mu"], R0_factor=r0f
        )
        conv[str(r0f)] = {
            "delta_recovered_deg": float(np.degrees(rec["delta_rad"])),
            "n_steps": rec["n_steps"],
        }
    d_ref = conv["100.0"]["delta_recovered_deg"]
    conv["max_spread_deg"] = max(abs(c["delta_recovered_deg"] - d_ref) for c in conv.values() if isinstance(c, dict))
    # integrator order verification on the hardest case (near-parabolic):
    # halving theta must cut the delta error by ~16 (RK4 order 4).
    p = PLANETS["Jupiter"]
    r_p_np = 1.05 * p["R_eq"]
    d_an = turn_angle_rad(x_parameter(r_p_np, 0.7, p["mu"]))
    errs = {}
    for th in (0.02, 0.01, 0.005):
        rec_t = propagate_flyby_l3(r_p_np, 0.7, p["mu"], theta=th)
        errs[str(th)] = {
            "delta_abs_err_rad": abs(rec_t["delta_rad"] - d_an),
            "n_steps": rec_t["n_steps"],
        }
    conv["order_verification_near_parabolic"] = {
        "errors": errs,
        "ratio_0.02_to_0.01": errs["0.02"]["delta_abs_err_rad"]
        / errs["0.01"]["delta_abs_err_rad"],
        "ratio_0.01_to_0.005": errs["0.01"]["delta_abs_err_rad"]
        / errs["0.005"]["delta_abs_err_rad"],
        "expected_ratio_order4": 16.0,
    }
    l3["patch_radius_convergence_jupiter"] = conv
    results["L3_propagation"] = l3

    # 5. sweep (analytic ceiling; compact summaries only)
    sweep = run_sweep()
    sweep_summary = {}
    for name, s in sweep.items():
        sweep_summary[name] = {
            "V_p_kms": s["V_p_kms"],
            "soi_km": s["soi_km"],
            "d_eps_max_at_rmin_over_v": [
                float(s["d_eps_max_grid_km2s2"][j, 0])
                for j in range(len(V_INF_KMS))
            ],
            "v_inf_grid_kms": [float(v) for v in V_INF_KMS],
            "d_eps_max_best_km2s2": float(np.max(s["d_eps_max_grid_km2s2"])),
            "d_eps_max_best_at": {
                "v_inf_kms": float(V_INF_KMS[
                    np.unravel_index(np.argmax(s["d_eps_max_grid_km2s2"]),
                                     s["d_eps_max_grid_km2s2"].shape)[0]]),
                "r_p_over_R_eq": float(R_P_FACTOR[
                    np.unravel_index(np.argmax(s["d_eps_max_grid_km2s2"]),
                                     s["d_eps_max_grid_km2s2"].shape)[1]]),
            },
        }
    results["sweep"] = sweep_summary

    # 6. pathological regimes
    results["pathological"] = pathological_checks()

    # figures
    make_figures(sweep, ori)

    save_json_result(
        str(RESULTS_DIR / "results.json"),
        results,
        name="gravity_assist_flyby",
        description=(
            "Patched-conic planetary gravity assist: orientation landscape, "
            "analytic energy-exchange ceiling, Voyager anchors, L3 Cowell "
            "validation with element recovery, sweep, pathological regimes."
        ),
    )
    return results


def make_figures(sweep: dict, ori: dict) -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    # Figure 1: analytic ceiling vs v_inf at r_p = 1.02 R_eq (axis 0 = v_inf)
    fig, ax = plt.subplots(figsize=(7, 5))
    for name, s in sweep.items():
        ax.loglog(V_INF_KMS, s["d_eps_max_grid_km2s2"][:, 0], label=name)
    ax.set_xlabel("v_inf [km/s]")
    ax.set_ylabel(r"$|\Delta\varepsilon|_{\max}$ [km$^2$/s$^2$]")
    ax.set_title(r"Max heliocentric energy exchange, $r_p = 1.02\,R_{eq}$")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "d_eps_max_vs_vinf.png", dpi=150)
    plt.close(fig)

    # Figure 2: turn-angle map (Jupiter)
    x = sweep["Jupiter"]["x_grid"]
    delta = np.degrees(turn_angle_rad(x))
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.pcolormesh(
        V_INF_KMS, R_P_FACTOR, delta.T, shading="auto", cmap="viridis"
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("v_inf [km/s]")
    ax.set_ylabel(r"$r_p / R_{eq}$")
    ax.set_title("Jupiter turn angle [deg]")
    fig.colorbar(im, ax=ax, label=r"$\delta$ [deg]")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "turn_angle_map_jupiter.png", dpi=150)
    plt.close(fig)

    # Figure 3: orientation landscape (Jupiter, r_p = 1.02 R_eq, v_inf = 10)
    key = "Jupiter_rp1.02_v10.0"
    h = 0.5 * np.radians(ori[key]["delta_deg"])
    am = np.radians(ORIENTATION_ALPHA_DEG)
    pm = np.radians(ORIENTATION_PHI_DEG)
    amg, pmg = np.meshgrid(am, pm, indexing="ij")
    F = orientation_factor(amg, pmg, h)
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.pcolormesh(
        ORIENTATION_ALPHA_DEG, ORIENTATION_PHI_DEG, F.T, shading="auto", cmap="RdBu_r"
    )
    a_star = np.degrees(np.pi / 2 + h)
    ax.plot([a_star], [0.0], "w*", markersize=14, label="analytic max")
    ax.plot([180.0 - a_star], [180.0], "k*", markersize=14, label="analytic min")
    ax.set_xlabel(r"$\alpha$ [deg]  (angle of $v_{\infty,in}$ from $V_p$)")
    ax.set_ylabel(r"$\phi$ [deg]  (flyby-plane phase)")
    ax.set_title(f"Jupiter flyby: d_eps landscape (delta = {ori[key]['delta_deg']:.1f} deg)")
    fig.colorbar(im, ax=ax, label=r"$F = \Delta\varepsilon / A$")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "orientation_landscape_jupiter.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    res = main()
    print("Anchors:")
    for k, v in res["anchors"].items():
        print(
            f"  {k}: v_inf={v['v_inf_kms']:.4f} r_p={v['r_p_km']:.0f} km "
            f"delta={v['delta_deg']:.3f} deg |dv|={v['dv_vector_kms']:.3f} "
            f"d_eps={v['d_epsilon_km2s2']:.2f} km2/s2 dv_scalar={v['dv_scalar_kms']:+.3f}"
        )
    print("Orientation (Jupiter low):")
    o = res["orientation_landscape"]["Jupiter_rp1.02_v10.0"]
    print(
        f"  delta={o['delta_deg']:.2f} deg alpha*={o['alpha_star_deg']:.2f} deg "
        f"A={o['A_km2s2']:.2f} grid_loss={o['grid_loss_pct']:.2e}% "
        f"dv_parallel_err={o['dv_parallel_err_at_optimum']:.2e} "
        f"d_eps_vs_A={o['d_eps_direct_vs_A_rel_err']:.2e}"
    )
    print("L3 (max delta rel err):",
          max(v["delta_rel_err"] for v in res["L3_propagation"].values()
              if isinstance(v, dict) and "delta_rel_err" in v))
    print("Pathological:", res["pathological"])
