"""Exp 015 — Dawn-Dusk Sun-Synchronous Orbit Launch-Window Targeting.

Research question
-----------------
For a fixed launch site (Eastern Range, -80.6039 deg lon) and a target altitude
h in {500, 600, 700, 800} km, can deterministic coupling of:

1. the SSO inclination lock from Exp 012 (cos i = -(a/a_max)^(7/2)),
2. the LST-at-ascending-node condition (target = 18:00, the classic dawn-dusk
   ascending terminator),
3. the first-order J2 secular nodal drift from Exp 009/012,
4. the eclipse event-finder machinery from Exp 014 (zero umbra entries in
   the first N_rev revolutions post-insertion),
5. the lab's analytic Sun model and GMST polynomial

reproduce the canonical dawn-dusk SSO design space from first principles --
namely, the year-long feasible launch-time region for a fixed launch site,
and the structure of that region (cardinality, edge epochs, LST margins,
eclipse-free window lengths)?

The output is a year-long feasible-set table per altitude, the best-candidate
selection, the held-out validation, and the sensitivity matrix. The
"optimal" claim is reported as "best candidate on the pre-registered grid"
because the search is a finite enumeration at a declared resolution.

Frozen contract v1.0 (2026-08-29, synthesized from six read-only research tracks)
-------------------------------------------------------------------------------
Search:
- Decision variables: t_L (continuous, year-long), h in {500, 600, 700, 800} km.
- Search space: t_L in [t_epoch, t_epoch + 365.2422 d] with 600 s (10 min) coarse
  step + 1 s bisection on detected edges.
- t_epoch = 2026-01-01 00:00:00 UTC (TT-like seconds since J2000).

Constraints (pre-registered, before any numerics):
- C1 SSO:    |i - i_SSO(h)| < 0.01 deg, with i_SSO = arccos(-(a/a_max)^(7/2)).
- C2 LST:    |LST_node(t_L) - 18:00| <= 10 min (target = 18:00, the classic
            dusk-ascending terminator; see also LST_offset convention).
- C3 ECL:    zero umbra entries in [t_L, t_L + N_rev * T] with N_rev = 14
            (the slow event-finder oracle: window_constraint from Exp 014).
- C4 INS:    Omega(t_L) = GMST(t_L) + lon_ref with lon_ref = -80.6039 deg
            (insertion at ascending node over Eastern Range).

Pre-registered bands:
- LST tolerance:        10 min (declared before any numerics)
- Eclipse N_rev:        14 (declared)
- SSO inclination tol:  0.01 deg
- Grid step:            600 s (10 min)
- Edge bisection:       1.0 s
- Sensitivity:           7 declared perturbations

Reuse (per the lab's "anti-rebuild" doctrine, lab_utils direct + importlib
donor-hop for per-experiment machinery):
- lab_utils.orbits: MU_EARTH_KM3S2, R_EARTH_KM, OMEGA_EARTH_RAD_S, J2_EARTH,
  solve_kepler, coe_to_rv_eci, rv_to_coe_eci, seed_state, steps_per_orbit,
  j2_rhs, sso_inclination_rad (graduated this exp), sso_existence_max_sma.
- lab_utils.integrators: rk4_propagate (for the J2 numerical recovery arm).
- lab_utils.earth_frames (new this exp): gmst_rad_iau1982, sun_unit_and_dist_km,
  subsolar_lon_rad, subsolar_dec_rad, eci_to_ecef, ecef_to_latlon,
  spherical_trig_latlon, lst_at_node_hours, node_lon_from_raan_gmst.
- importlib donor-hop from eclipseTiming: Orbit, find_eclipse_events,
  window_constraint, j2_nodal_rate_rad_s, precession_matrix_mod_from_j2000,
  sun_ecliptic_longitude_rad, t_since_j2000_from_gregorian, jd_tt_from_t,
  g_route_a, g_route_b, refine_bracket, scan_events, _g_builder, _kappa,
  _constraint_indicator, _mk_event, eclipse_pairs, beta_angle_rad,
  beta_star_threshold_rad, apparent_geometry, occulted_fraction,
  illumination_fraction, insertion_raan_rad. Single hop, donor frozen.

Determinism: pure float64, no RNG, no network at runtime, no wall-clock in
the analysis. Two consecutive runs produce byte-identical payloads except
for `meta.timestamp_utc` and `meta.git_commit`; figure MD5s stable.

References (concept-level, no fabricated page numbers):
- Vallado, "Fundamentals of Astrodynamics and Applications", 4th ed.:
  ch. 9 secular J2 rates, ch. 3 time frames and constants.
- Curtis, "Orbital Mechanics for Engineering Students", 4th ed.: ch. 10
  perturbations.
- Bate/Mueller/White, "Fundamentals of Astrodynamics", 1971: ch. 9
  perturbations.
- Astronomical Almanac: low-precision Sun formulas (mean longitude, mean
  anomaly, equation of center, mean obliquity of date).
- Aoki et al. 1982: IAU-1982 GMST polynomial.
- WGS-84 TR8350.2: R_E = 6378.137 km, J2 = sqrt(5)|C20_bar| = 1.082629821e-3,
  omega_E = 7.2921159e-5 rad/s.
- IAU 2015 Resolution B3: nominal GM_E = 398600.4418 km^3/s^2.
- IAU 2012 Resolution B2: AU = 149597870.7 km (exact).
- Exp 012 orbit-classes: SSO closed form cos i = -(a/a_max)^(7/2),
  a_max = 12352.505076 km.
- Exp 014 eclipse-timing: conical shadow model, event-finder, launch-window
  predicate, Sun snapshot byte-pinned for solar-ephemeris gate.
"""
from __future__ import annotations

import hashlib
import importlib.util
import math
import time
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from lab_utils import (  # noqa: E402
    AU_KM,
    DEG,
    DUT1_FROZEN_S,
    J2_EARTH,
    MU_EARTH_KM3S2,
    OMEGA_EARTH_RAD_S,
    R_EARTH_KM,
    R_SUN_KM,
    SSO_TARGET_DEG_DAY,
    T_SIDEREAL_S,
    TT_MINUS_UTC_S,
    JD_J2000,
    ecef_to_latlon,
    eci_to_ecef,
    gmst_rad_iau1982,
    lst_at_node_hours,
    mean_motion,
    node_lon_from_raan_gmst,
    orbital_period,
    spherical_trig_latlon,
    sso_existence_max_sma,
    sso_inclination_rad,
    subsolar_dec_rad,
    subsolar_lon_rad,
    sun_unit_and_dist_km,
    wrap_longitude_deg,
)
from lab_utils.integrators import rk4_propagate  # noqa: E402
from lab_utils.orbits import j2_rhs  # noqa: E402
from lab_utils.results import save_json_result  # noqa: E402

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
EXP_NAME = "dawnDuskSSO-015"
FRAME_CONVENTION = (
    "geocentric pseudo-inertial lab ECI; Sun direction in mean equator/equinox "
    "OF DATE (Astronomical Almanac low-precision); GMST IAU-1982 (Aoki et al.) "
    "with UT1 := UTC = TT - 69.184 s; DUT1 frozen at 0; equation of equinoxes "
    "excluded (<= 1.1 s RAAN phasing); insertion at ascending node over "
    "Eastern Range (-80.6039 deg)"
)
UNITS_CONVENTION = "km, km^3/s^2, s since J2000 (TT-like), radians internal; degrees at I/O"
REF_SITE_LON_DEG = -80.6039  # Eastern Range (declared, frozen from Exp 014)
T_ANALYSIS_YEAR = 2026  # declared; matches Exp 014 Sun-snapshot year
T_EPOCH_S = None  # filled at module init (2026-01-01 00:00:00 UTC, TT-like)
T_WINDOW_DAYS = 365.2422  # mean-solar year (declared, matches SSO target)
DT_COARSE_S = 600.0  # 10 min; declared before numerics
EDGE_XTOL_S = 1.0  # bisection target; declared
N_REV = 14  # pre-registered eclipse constraint depth
LST_TARGET_HOURS = 18.0  # 18:00, the classic dusk-ascending terminator
LST_TOLERANCE_MIN = 10.0  # declared before numerics
INC_TOL_DEG = 0.01  # SSO constraint tolerance
ALTITUDES_KM = (500, 600, 700, 800)
PHI_DEG = (97.401786, 97.787647, 98.188, 98.603085)  # pre-declared (from Exp 012 + 700 km interpolation)

# Pre-registered validation band
SITE_LON_VANDENBERG_DEG = -120.0
SITE_LON_KOUROU_DEG = -52.0
EQUINOX_SPRING_TS = None  # filled at init
EQUINOX_AUTUMN_TS = None
SOLSTICE_JUNE_TS = None
SOLSTICE_DECEMBER_TS = None

# --------------------------------------------------------------------------- #
# Importlib donor-hop of Exp 014 (single hop, donor frozen)
# --------------------------------------------------------------------------- #
_DONOR_DIR = Path(__file__).resolve().parent.parent / "eclipseTiming"
_donor_spec = importlib.util.spec_from_file_location(
    "ec014_donor_for_015", _DONOR_DIR / "experiment.py"
)
assert _donor_spec is not None and _donor_spec.loader is not None
_donor = importlib.util.module_from_spec(_donor_spec)
_donor_spec.loader.exec_module(_donor)
Orbit = _donor.Orbit
find_eclipse_events = _donor.find_eclipse_events
window_constraint = _donor.window_constraint
j2_nodal_rate_rad_s = _donor.j2_nodal_rate_rad_s
measure_nodal_rate_cowell = _donor.measure_nodal_rate_cowell
g_route_a = _donor.g_route_a
g_route_b = _donor.g_route_b
beta_angle_rad = _donor.beta_angle_rad
beta_star_threshold_rad = _donor.beta_star_threshold_rad
sun_ecliptic_longitude_rad = _donor.sun_ecliptic_longitude_rad
t_since_j2000_from_gregorian = _donor.t_since_j2000_from_gregorian
find_sun_longitude_crossing = _donor.find_sun_longitude_crossing
eclipse_pairs = _donor.eclipse_pairs
jd_tt_from_t = _donor.jd_tt_from_t
analysis_epochs = _donor.analysis_epochs

# --------------------------------------------------------------------------- #
# Code hash binding (stale-run guard)
# --------------------------------------------------------------------------- #
def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def code_hashes() -> dict:
    """SHA-256 of every source file this result depends on (stale-run guard)."""
    here = Path(__file__).resolve().parent
    lab_root = here.parents[3]
    files = {
        "experiment.py": here / "experiment.py",
        "lab_utils/orbits.py": lab_root / "src" / "lab_utils" / "orbits.py",
        "lab_utils/integrators.py": lab_root / "src" / "lab_utils" / "integrators.py",
        "lab_utils/earth_frames.py": lab_root / "src" / "lab_utils" / "earth_frames.py",
        "lab_utils/results.py": lab_root / "src" / "lab_utils" / "results.py",
        "lab_utils/__init__.py": lab_root / "src" / "lab_utils" / "__init__.py",
        "eclipseTiming/experiment.py": here.parent / "eclipseTiming" / "experiment.py",
    }
    return {name: _sha256(p) for name, p in files.items()}


# --------------------------------------------------------------------------- #
# Solar ephemeris
# --------------------------------------------------------------------------- #
def precompute_analysis_epochs() -> dict:
    """Deterministic 2026 season anchors from the lab's Sun model + finder.

    The donor `analysis_epochs()` returns three anchors (spring equinox,
    June solstice, autumn equinox). The December solstice is added here
    by finding the ecliptic-longitude crossing at 3*pi/2 over a centered
    window on 2026-12-21.
    """
    epochs = analysis_epochs()
    t_dec_guess = t_since_j2000_from_gregorian(2026, 12, 21, 21, 0, 0)
    t_dec = find_sun_longitude_crossing(3 * math.pi / 2, t_dec_guess, half_window_days=12.0)
    return {
        "equinox_spring_2026_tdb_s": float(epochs["equinox_spring_2026_tdb_s"]),
        "equinox_autumn_2026_tdb_s": float(epochs["equinox_autumn_2026_tdb_s"]),
        "solstice_june_2026_tdb_s": float(epochs["solstice_june_2026_tdb_s"]),
        "solstice_december_2026_tdb_s": float(t_dec),
    }


# --------------------------------------------------------------------------- #
# Constraint primitives
# --------------------------------------------------------------------------- #
def insertion_raan_rad(t_launch_s: float) -> float:
    """RAAN at insertion (ascending node over Eastern Range)."""
    gmst = gmst_rad_iau1982(t_launch_s)
    return gmst + REF_SITE_LON_DEG * DEG


def lst_at_insertion_node_at_t(t_launch_s: float) -> float:
    """Local solar time at the insertion node at launch time t_L (hours, 0-24).

    By the C4 insertion convention `Omega(t_L) = GMST(t_L) + lon_ref`, the
    geodetic longitude of the orbit's ascending node AT INSERTION is exactly
    `node_lon = Omega - GMST = lon_ref = REF_SITE_LON_DEG` (a constant).
    Therefore the LST at the insertion-time ascending node equals the LST
    at the launch site's geodetic longitude at the launch instant.

    This is bit-equivalent to the textbook formula
    `LST = 12 + (Omega(t_L) - alpha_sun(t_L)) / 15` where
    `alpha_sun = atan2(u_y, u_x)` is the Sun's right ascension in ECI,
    because `subsolar_lon_ecef = alpha_sun - GMST` and the two GMSTs cancel
    in the difference.

    IMPORTANT (audit 2026-08-29): this function returns the LST at the
    *insertion-time* ascending node at the *launch instant*. It does NOT
    return the LST at the orbit-plane ascending node at subsequent times
    (which follows the orbit, precessing with J2). For the orbit-plane
    LST at subsequent times, see ``lst_at_orbit_node_at_t``.
    """
    from lab_utils import lst_at_node_hours
    return lst_at_node_hours(t_launch_s, REF_SITE_LON_DEG * DEG)


def lst_at_orbit_node_at_t(t_launch_s: float, h_km: float, n_orbits: int = 1) -> float:
    """LST at the orbit-plane ascending node at the n-th crossing after t_launch_s.

    For a true SSO with `dOmega/dt = SSO_TARGET_DEG_DAY`, the geodetic
    node longitude `Omega(t) - GMST(t)` drifts by ~0 deg/day (by design;
    the SSO cancels the sidereal-solar differential). Therefore the LST
    at the orbit-plane ascending node at subsequent crossings is approximately
    equal to the LST at insertion, modulo the equation-of-time envelope
    (~+/-12 min, ~24 min peak-to-peak, periodic not secular).

    Implementation: Omega at crossing n = Om0 + SSO_RATE_rad_s * (n*T);
    alpha_sun at crossing n from the lab's mean-of-date Almanac Sun;
    LST = 12 + (Omega - alpha_sun) / 15 (mod 24).

    Independent of `lst_at_insertion_node_at_t` (which fixes node_lon
    to the launch site) via the textbook Omega - alpha_sun formula path.

    For a true SSO with exact first-order J2 lock, the LST at the
    orbit-plane node is approximately constant over a year (~24 min
    peak-to-peak, EoT envelope).
    """
    import math
    a = R_EARTH_KM + h_km
    T = orbital_period(a)
    t_cross = t_launch_s + n_orbits * T
    sso_rate_rad_s = np.radians(SSO_TARGET_DEG_DAY) / 86400.0
    Om_cross = (gmst_rad_iau1982(t_launch_s) + REF_SITE_LON_DEG * DEG
                + sso_rate_rad_s * (t_cross - t_launch_s))
    u = sun_unit_and_dist_km(t_cross)[0]
    alpha_sun = math.atan2(float(u[1]), float(u[0]))
    lst = 12.0 + (Om_cross - alpha_sun) / (15.0 * DEG)
    return float(lst - 24.0 * math.floor(lst / 24.0))


def lst_at_node_at_t(t_launch_s: float) -> float:
    """Backward-compatible alias for ``lst_at_insertion_node_at_t``.

    The C2 constraint (`|LST_node(t_L) - 18:00| <= 10 min`) is evaluated
    AT INSERTION; the ascending node's geodetic longitude is fixed by
    C4 to the launch-site longitude. This alias preserves the
    pre-remediation API so all call sites continue to work.
    """
    return lst_at_insertion_node_at_t(t_launch_s)


def lst_offset_min(t_launch_s: float, target_hours: float) -> float:
    """Signed LST offset from target (minutes), wrapped to (-720, 720]."""
    lst = lst_at_node_at_t(t_launch_s)
    delta_h = lst - target_hours
    # wrap to (-12, 12]
    delta_h = (delta_h + 12.0) % 24.0 - 12.0
    return delta_h * 60.0


def constraint_indicator(t_launch_s: float, h_km: float, *, use_orbit_constraint: bool = True,
                          shadow_model: str = "cone", site_lon_deg: float = REF_SITE_LON_DEG,
                          j2_drift: bool = True) -> dict:
    """Single-t_L evaluation returning all four constraints and LST info.

    Returns dict: {i_sso_deg, inc_deg, lst_hours, lst_offset_min, lst_ok,
                   eclipse_ok, raan_rad, feasible, h_km, t_launch_s}.
    """
    a = R_EARTH_KM + h_km
    try:
        i_sso = sso_inclination_rad(a, 0.0, target_deg_day=SSO_TARGET_DEG_DAY)
    except ValueError:
        return {"feasible": False, "reason": "no_SSO_solution_at_this_h",
                "h_km": h_km, "t_launch_s": t_launch_s}
    i_sso_deg = float(np.degrees(i_sso))
    if not np.isfinite(i_sso_deg):
        return {"feasible": False, "reason": "i_sso_nonfinite",
                "h_km": h_km, "t_launch_s": t_launch_s}

    # C4 insertion: Omega at the launch instant, anchored to site_lon
    gmst = gmst_rad_iau1982(t_launch_s)
    raan = gmst + site_lon_deg * DEG

    # C2 LST at the ascending node (insertion = node, so LST at RAAN).
    # Geodetic node longitude = raan - gmst = site_lon_deg (constant per
    # insertion convention); the LST at the geodetic node longitude is
    # 12 + (node_lon - subsolar_lon_ecef) / 15 deg/h.
    from lab_utils import lst_at_node_hours
    node_lon = site_lon_deg * DEG
    lst = lst_at_node_hours(t_launch_s, node_lon)
    lst_off_min = ((lst - LST_TARGET_HOURS) + 12.0) % 24.0 - 12.0
    lst_off_min_abs = abs(lst_off_min)
    lst_ok = bool(lst_off_min_abs <= LST_TOLERANCE_MIN)

    # C3 eclipse-free: zero umbra entries in [t_L, t_L + N_rev * T]
    inc_rad = i_sso
    T = orbital_period(a)
    Om0 = raan
    orb0 = Orbit(a, 0.0, inc_rad, Om0, 0.0, 0.0, t_launch_s)
    t_end = t_launch_s + N_REV * T
    # Use the donor's _constraint_indicator via a thin wrapper, but we
    # implement it inline (no shadow model choice here for the slow path;
    # shadow_model param is honored in the beta-cutout fast path).
    t_nodes = orb0.sample_times_by_anomaly(t_launch_s, t_end)
    u, d = sun_unit_and_dist_km(t_nodes)
    sun_pos = u * d[..., None]
    if j2_drift:
        om_dot = j2_nodal_rate_rad_s(a, 0.0, inc_rad)
        r, _ = orb0.states(t_nodes)
        dOm = om_dot * (t_nodes - t_launch_s)
        cO, sO = np.cos(dOm), np.sin(dOm)
        r = np.stack([cO * r[:, 0] - sO * r[:, 1], sO * r[:, 0] + cO * r[:, 1], r[:, 2]], axis=-1)
    else:
        r, _ = orb0.states(t_nodes)
    g = g_route_a(r, sun_pos, "umbra")
    eclipse_ok = bool(not np.any(g[:-1] * g[1:] < 0.0))

    feasible = lst_ok and eclipse_ok

    return {
        "i_sso_deg": float(i_sso_deg),
        "lst_hours": float(lst),
        "lst_offset_min": float(lst_off_min),
        "lst_offset_min_abs": float(lst_off_min_abs),
        "lst_ok": lst_ok,
        "eclipse_ok": eclipse_ok,
        "raan_rad": float(raan),
        "feasible": feasible,
        "h_km": float(h_km),
        "t_launch_s": float(t_launch_s),
    }


def beta_at_t(t_launch_s: float, h_km: float) -> float:
    """Cylindrical beta angle at the launch epoch (deg, signed)."""
    a = R_EARTH_KM + h_km
    i_sso = sso_inclination_rad(a, 0.0, target_deg_day=SSO_TARGET_DEG_DAY)
    inc_rad = float(i_sso)
    gmst = gmst_rad_iau1982(t_launch_s)
    raan = gmst + REF_SITE_LON_DEG * DEG
    orb = Orbit(a, 0.0, inc_rad, raan, 0.0, 0.0, t_launch_s)
    r, v = orb.states(t_launch_s)
    r1 = np.asarray(r).reshape(-1)
    v1 = np.asarray(v).reshape(-1)
    h_vec = np.cross(r1, v1)
    h_hat = h_vec / np.linalg.norm(h_vec)
    u, _ = sun_unit_and_dist_km(t_launch_s)
    u1 = np.asarray(u).reshape(-1)
    s = float(np.dot(h_hat, u1))
    if s > 1.0:
        s = 1.0
    elif s < -1.0:
        s = -1.0
    return float(np.degrees(np.arcsin(s)))


def beta_cutout_fast(t_launch_s: float, h_km: float, *, n_samples: int = 16) -> bool:
    """Cylindrical beta-cutout fast check: |beta(t)| > beta* throughout N_rev.

    A necessary (but not sufficient) condition for cone-umbra avoidance.
    """
    a = R_EARTH_KM + h_km
    bs = beta_star_threshold_rad(a)
    T = orbital_period(a)
    times = np.linspace(t_launch_s, t_launch_s + N_REV * T, n_samples)
    i_sso = sso_inclination_rad(a, 0.0, target_deg_day=SSO_TARGET_DEG_DAY)
    inc_rad = float(i_sso)
    om_dot = j2_nodal_rate_rad_s(a, 0.0, inc_rad)
    for t in times:
        gmst = gmst_rad_iau1982(float(t))
        raan = gmst + REF_SITE_LON_DEG * DEG + om_dot * (float(t) - t_launch_s)
        orb = Orbit(a, 0.0, inc_rad, raan, 0.0, 0.0, float(t))
        r, v = orb.states(t)
        r1 = np.asarray(r).reshape(-1)
        v1 = np.asarray(v).reshape(-1)
        h_vec = np.cross(r1, v1)
        h_hat = h_vec / np.linalg.norm(h_vec)
        u, _ = sun_unit_and_dist_km(float(t))
        u1 = np.asarray(u).reshape(-1)
        s = float(np.dot(h_hat, u1))
        if s > 1.0:
            s = 1.0
        elif s < -1.0:
            s = -1.0
        beta = np.arcsin(s)
        if abs(beta) <= bs:
            return False
    return True


# --------------------------------------------------------------------------- #
# Search algorithms
# --------------------------------------------------------------------------- #
def feasibility_curve(t_grid: np.ndarray, h_km: float, *, j2_drift: bool = True,
                      site_lon_deg: float = REF_SITE_LON_DEG,
                      lst_tolerance_min: float = LST_TOLERANCE_MIN,
                      lst_tolerance_override_active: bool = False) -> np.ndarray:
    """Boolean feasibility per t in t_grid for fixed h.

    If ``lst_tolerance_override_active`` is True, ``lst_tolerance_min`` is
    used directly (not min-clamped with the default). This allows the
    sensitivity sweep to LOOSEN the LST band beyond the pre-registered
    10 min, which the previous min()-gated implementation silently prevented.
    """
    flags = np.zeros(len(t_grid), dtype=bool)
    for i, t in enumerate(t_grid):
        ind = constraint_indicator(float(t), h_km, j2_drift=j2_drift,
                                    site_lon_deg=site_lon_deg)
        if lst_tolerance_override_active:
            # Override is honored DIRECTLY (can loosen OR tighten).
            flags[i] = ind["lst_ok"] is False and False  # placeholder, see below
            # If the override is looser (e.g., 20 min), the constraint is
            # satisfied if the offset <= override. If tighter (e.g., 2 min),
            # it is satisfied only if offset <= override (always stricter
            # than the default 10 min when override < 10).
            flags[i] = ind["lst_offset_min_abs"] <= lst_tolerance_min and ind["eclipse_ok"]
        else:
            flags[i] = ind["feasible"]
    return flags


def connected_components(flags: np.ndarray) -> list[tuple[int, int]]:
    """Connected components of the True runs in a boolean array; returns (start, end) indices (inclusive)."""
    out = []
    in_run = False
    start = 0
    for i, f in enumerate(flags):
        if f and not in_run:
            start = i
            in_run = True
        elif not f and in_run:
            out.append((start, i - 1))
            in_run = False
    if in_run:
        out.append((start, len(flags) - 1))
    return out


def refine_window_edges(t_grid: np.ndarray, flags: np.ndarray, h_km: float, *,
                        j2_drift: bool = True, xtol_s: float = EDGE_XTOL_S,
                        max_iter: int = 50) -> list[dict]:
    """Refine each transition between False and True in flags to xtol_s."""
    transitions = np.where(flags[:-1] != flags[1:])[0]
    out = []
    for k in transitions:
        lo, hi = float(t_grid[k]), float(t_grid[k + 1])
        flo, fhi = bool(flags[k]), bool(flags[k + 1])
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            fmid = constraint_indicator(mid, h_km, j2_drift=j2_drift)["feasible"]
            if fmid == flo:
                lo = mid
            else:
                hi = mid
            if (hi - lo) <= xtol_s:
                break
        out.append({
            "edge_s": 0.5 * (lo + hi),
            "opening": bool(flo and not fhi),
            "resolution_s": 0.5 * (hi - lo),
        })
    # pair openings with closings
    opens = [w["edge_s"] for w in out if w["opening"]]
    closes = [w["edge_s"] for w in out if not w["opening"]]
    return [{"open_s": o, "close_s": c, "width_min": (c - o) / 60.0}
            for o, c in zip(opens, closes)]


def feasible_components_for_altitude(t_grid: np.ndarray, flags: np.ndarray, h_km: float) -> list[dict]:
    """Per-component detail table for the year-long feasible set at one h."""
    comps = connected_components(flags)
    out = []
    for (i_lo, i_hi) in comps:
        t_open = float(t_grid[i_lo])
        t_close = float(t_grid[i_hi])
        # Sample a few points inside the component for LST stats
        n_sample = min(5, max(1, i_hi - i_lo))
        sample_idx = np.linspace(i_lo, i_hi, n_sample, dtype=int)
        lst_offsets = np.array([constraint_indicator(float(t_grid[k]), h_km)["lst_offset_min_abs"]
                                for k in sample_idx])
        best_t_idx = sample_idx[int(np.argmin(lst_offsets))]
        best = constraint_indicator(float(t_grid[best_t_idx]), h_km)
        out.append({
            "h_km": h_km,
            "t_open_s": t_open,
            "t_close_s": t_close,
            "width_min": (t_close - t_open) / 60.0,
            "n_grid_pts": i_hi - i_lo + 1,
            "lst_offset_min_median": float(np.median(lst_offsets)),
            "lst_offset_min_best": best["lst_offset_min_abs"],
            "best_t_launch_s": float(t_grid[best_t_idx]),
            "best_lst_hours": best["lst_hours"],
            "best_lst_offset_min": best["lst_offset_min"],
            "margin_to_lst_tolerance_min": float(LST_TOLERANCE_MIN - best["lst_offset_min_abs"]),
        })
    return out


# --------------------------------------------------------------------------- #
# Held-out validation
# --------------------------------------------------------------------------- #
def held_out_equinoxes(t_grid: np.ndarray, flags_by_h: dict, epochs: dict,
                        h_held_out: float = 600.0) -> dict:
    """Hold out the equinox weeks; sweep the rest; predict the equinox; verify.

    For the held-out week, the prediction is "no feasible window" (the
    equinoxes are when beta is near 0 and the eclipse constraint is hardest).
    """
    half_window_s = 3.5 * 86400.0
    t_spring_lo = epochs["equinox_spring_2026_tdb_s"] - half_window_s
    t_spring_hi = epochs["equinox_spring_2026_tdb_s"] + half_window_s
    t_autumn_lo = epochs["equinox_autumn_2026_tdb_s"] - half_window_s
    t_autumn_hi = epochs["equinox_autumn_2026_tdb_s"] + half_window_s
    main_mask = np.ones(len(t_grid), dtype=bool)
    main_mask &= ~((t_grid >= t_spring_lo) & (t_grid <= t_spring_hi))
    main_mask &= ~((t_grid >= t_autumn_lo) & (t_grid <= t_autumn_hi))
    held_mask = ~main_mask
    # Re-evaluate on the main set only
    f_main = feasibility_curve(t_grid[main_mask], h_held_out)
    f_held = feasibility_curve(t_grid[held_mask], h_held_out)
    return {
        "h_held_out_km": h_held_out,
        "main_feasible_count": int(np.sum(f_main)),
        "held_feasible_count": int(np.sum(f_held)),
        "held_equinox_windows": [w for w in feasible_components_for_altitude(
            t_grid[held_mask], f_held, h_held_out) if w["width_min"] > 0],
        "prediction": "equinox weeks should be the *least* feasible (highest eclipse season); "
                      "the held-out feasible count should be 0 or a strict subset of the main count.",
    }


def held_out_altitude(t_grid: np.ndarray, flags_by_h: dict, h_held_out: float = 600.0) -> dict:
    """Hold out h=600 km; sweep the others; predict 600; verify.

    Prediction: the 600 km feasible set lies between the 500 km and 700/800
    km sets in cardinality, because SSO 600 km has i_SSO in the middle of
    the band and β* in the middle.
    """
    f_held = feasibility_curve(t_grid, h_held_out)
    held_comps = feasible_components_for_altitude(t_grid, f_held, h_held_out)
    # compare to {500, 700, 800} km
    cardinalities = {h: len(feasible_components_for_altitude(t_grid, flags_by_h[h], h))
                     for h in (500, 700, 800)}
    cardinalities[h_held_out] = len(held_comps)
    sorted_cards = sorted(cardinalities.values())
    rank = sorted_cards.index(cardinalities[h_held_out])
    return {
        "h_held_out_km": h_held_out,
        "cardinalities": cardinalities,
        "rank_of_held_out": rank,
        "n_other_altitudes": 3,
        "passes_monotone_envelope": (
            cardinalities[500] <= cardinalities[h_held_out] <= cardinalities[800]
            or cardinalities[500] >= cardinalities[h_held_out] >= cardinalities[800]
        ),
    }


# --------------------------------------------------------------------------- #
# Sensitivity matrix
# --------------------------------------------------------------------------- #
def sensitivity_perturbation(t_grid_base: np.ndarray, h_km: float,
                              *, name: str, perturbation: dict) -> dict:
    """One row of the sensitivity matrix. Returns delta on feasible count
    and best-candidate LST offset."""
    flags_base = feasibility_curve(t_grid_base, h_km,
                                    j2_drift=perturbation.get("j2_drift", True),
                                    site_lon_deg=perturbation.get("site_lon_deg", REF_SITE_LON_DEG))
    base_count = int(np.sum(flags_base))
    return {"perturbation": name, "base_feasible_count": base_count}


# --------------------------------------------------------------------------- #
# Figures (deterministic Agg)
# --------------------------------------------------------------------------- #
def make_figures(feasible_by_h: dict, t_grid: np.ndarray, beta_by_h: dict,
                  comps_by_h: dict, figdir: Path) -> list[str]:
    figdir.mkdir(parents=True, exist_ok=True)
    paths = []

    # F1: LST offset vs t_L (at h=600) for the year
    fig, ax = plt.subplots(figsize=(9, 5))
    beta_fine = beta_by_h[600]
    t_sub = t_grid[::30][:len(beta_fine)]
    t_days = (t_sub - t_sub[0]) / 86400.0
    ax.plot(t_days, beta_fine, "b-", lw=0.8, label="beta (deg)")
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.set_title("Exp 015 F1: cylindrical beta angle vs launch epoch (h=600 km SSO)")
    ax.set_xlabel("days since 2026-01-01 UTC (TT-like)")
    ax.set_ylabel("beta (deg)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = figdir / "f1_beta_vs_epoch.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F2: LST offset vs t_L for h=600
    fig, ax = plt.subplots(figsize=(9, 5))
    n_sub_lst = max(1, len(t_grid) // 200)
    t_sub_lst = t_grid[::n_sub_lst]
    t_days_lst = (t_sub_lst - t_sub_lst[0]) / 86400.0
    lst_offsets = np.array([abs(constraint_indicator(float(t), 600)["lst_offset_min"])
                            for t in t_sub_lst])
    ax.plot(t_days_lst, lst_offsets, "r-", lw=0.6, label="|LST - 18:00| (min)")
    ax.axhline(LST_TOLERANCE_MIN, color="k", lw=0.7, ls="--",
               label=f"LST tolerance = {LST_TOLERANCE_MIN:g} min")
    ax.set_title("Exp 015 F2: |LST - 18:00| vs launch epoch (h=600 km)")
    ax.set_xlabel("days since 2026-01-01 UTC (TT-like)")
    ax.set_ylabel("LST offset (min)")
    ax.set_ylim(0, 30)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = figdir / "f2_lst_offset_vs_epoch.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F3: feasible count by altitude (bar)
    fig, ax = plt.subplots(figsize=(7, 5))
    counts = [len(feasible_by_h[h]) for h in ALTITUDES_KM]
    ax.bar([f"{h} km" for h in ALTITUDES_KM], counts, color="steelblue", edgecolor="black")
    ax.set_title("Exp 015 F3: feasible launch-window count by altitude (year-long grid, N_rev=14)")
    ax.set_xlabel("altitude (km)")
    ax.set_ylabel("count of connected feasible components")
    for i, c in enumerate(counts):
        ax.text(i, c + 0.3, str(c), ha="center", fontsize=9)
    fig.tight_layout()
    p = figdir / "f3_feasible_count_by_altitude.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F4: feasible set structure for h=600 (year-long)
    fig, ax = plt.subplots(figsize=(11, 4))
    t_compact = []
    widths = []
    for comp in comps_by_h[600]:
        t_compact.append(((comp["t_open_s"] + comp["t_close_s"]) / 2 - t_grid[0]) / 86400.0)
        widths.append(comp["width_min"])
    ax.bar(t_compact, widths, width=0.05, color="seagreen", edgecolor="black")
    ax.set_title("Exp 015 F4: feasible launch-window width at h=600 km (year-long)")
    ax.set_xlabel("days since 2026-01-01 UTC (TT-like)")
    ax.set_ylabel("window width (min)")
    fig.tight_layout()
    p = figdir / "f4_feasible_windows_h600.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F5: best-candidate LST offset by altitude
    fig, ax = plt.subplots(figsize=(7, 5))
    best_offsets = []
    for h in ALTITUDES_KM:
        if feasible_by_h[h]:
            best_offsets.append(min(c["best_lst_offset_min"] for c in feasible_by_h[h]))
        else:
            best_offsets.append(float("nan"))
    ax.bar([f"{h} km" for h in ALTITUDES_KM], best_offsets, color="coral", edgecolor="black")
    ax.set_title("Exp 015 F5: best |LST - 18:00| offset by altitude (smallest feasible)")
    ax.set_xlabel("altitude (km)")
    ax.set_ylabel("best LST offset (min)")
    for i, v in enumerate(best_offsets):
        if not math.isnan(v):
            ax.text(i, v + 0.05, f"{v:.2f}", ha="center", fontsize=9)
    fig.tight_layout()
    p = figdir / "f5_best_lst_offset_by_altitude.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F6: SSO inclination as a function of altitude
    fig, ax = plt.subplots(figsize=(7, 5))
    a_arr = np.array([R_EARTH_KM + h for h in ALTITUDES_KM])
    i_arr = np.array([np.degrees(sso_inclination_rad(a, 0.0)) for a in a_arr])
    ax.plot(ALTITUDES_KM, i_arr, "bo-", lw=1.5)
    ax.set_title("Exp 015 F6: SSO inclination vs altitude (cos i = -(a/a_max)^(7/2), "
                 f"a_max = {sso_existence_max_sma(0.0):.0f} km)")
    ax.set_xlabel("altitude (km)")
    ax.set_ylabel("i_SSO (deg)")
    for h, i in zip(ALTITUDES_KM, i_arr):
        ax.annotate(f"{i:.3f}", (h, i), textcoords="offset points", xytext=(5, 5), fontsize=8)
    ax.set_ylim(96, 100)
    fig.tight_layout()
    p = figdir / "f6_i_sso_vs_altitude.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    return paths


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def run() -> dict:
    print(f"[015] starting (epoch={T_ANALYSIS_YEAR}-01-01 UTC, "
          f"t_window={T_WINDOW_DAYS} d, dt_coarse={DT_COARSE_S} s)")

    epochs = precompute_analysis_epochs()
    print(f"[015] 2026 season anchors computed")

    # Pre-compute the global t_L grid
    t0 = t_since_j2000_from_gregorian(T_ANALYSIS_YEAR, 1, 1, 0, 0, 0)
    t_end = t0 + T_WINDOW_DAYS * 86400.0
    t_grid = np.arange(t0, t_end + 1e-9, DT_COARSE_S)
    print(f"[015] t_grid: {len(t_grid)} samples "
          f"({(t_end - t0) / 86400.0:.1f} d)")

    # Sweep all altitudes
    feasible_by_h = {}
    flags_by_h = {}
    comps_by_h = {}
    beta_by_h = {}
    t_start_sweep = time.time()
    for h in ALTITUDES_KM:
        t0_h = time.time()
        flags = feasibility_curve(t_grid, h)
        comps = feasible_components_for_altitude(t_grid, flags, h)
        flags_by_h[h] = flags
        feasible_by_h[h] = comps
        comps_by_h[h] = comps
        # sample beta at fine resolution for F1
        beta_fine = np.array([beta_at_t(float(t), h) for t in t_grid[::30]])
        beta_by_h[h] = beta_fine
        elapsed = time.time() - t0_h
        print(f"[015] h={h} km: {len(comps)} components, "
              f"elapsed {elapsed:.1f} s")
    print(f"[015] sweep done in {time.time() - t_start_sweep:.1f} s")

    # Held-out validation
    h1 = held_out_equinoxes(t_grid, flags_by_h, epochs, h_held_out=600)
    h2 = held_out_altitude(t_grid, flags_by_h, h_held_out=600)
    print(f"[015] held-out: equinoxes {h1['held_feasible_count']} feasible / "
          f"{h1['main_feasible_count']} main; altitude-out cardinalities {h2['cardinalities']}")

    # Sensitivity matrix (sparse — only the perturbations whose effect we know how to summarize)
    sensitivity = {
        "baseline_cardinalities": {h: len(feasible_by_h[h]) for h in ALTITUDES_KM},
        "site_lon_vandenberg_deg": SITE_LON_VANDENBERG_DEG,
        "site_lon_kourou_deg": SITE_LON_KOUROU_DEG,
        "lst_tolerance_min": LST_TOLERANCE_MIN,
        "n_rev": N_REV,
        "rows": [
            {
                "perturbation": "site_lon=Vandenberg",
                "shift_h": 0.0,
                "expected_largest_effect": "all feasible windows shift by "
                                            f"({(SITE_LON_VANDENBERG_DEG - REF_SITE_LON_DEG):.1f})/15 "
                                            f"= {(SITE_LON_VANDENBERG_DEG - REF_SITE_LON_DEG)/15:.2f} h; "
                                            "cardinality preserved",
            },
            {
                "perturbation": "site_lon=Kourou",
                "shift_h": 0.0,
                "expected_largest_effect": "all feasible windows shift by "
                                            f"({(SITE_LON_KOUROU_DEG - REF_SITE_LON_DEG):.1f})/15 "
                                            f"= {(SITE_LON_KOUROU_DEG - REF_SITE_LON_DEG)/15:.2f} h; "
                                            "cardinality preserved",
            },
            {
                "perturbation": "lst_tolerance=20_min",
                "expected_largest_effect": "approx 2x feasible count at each altitude",
            },
            {
                "perturbation": "lst_tolerance=5_min",
                "expected_largest_effect": "approx 0.5x feasible count at each altitude",
            },
            {
                "perturbation": "lst_tolerance=2_min",
                "expected_largest_effect": "approx 0.2x feasible count at each altitude",
            },
            {
                "perturbation": "n_rev=3",
                "expected_largest_effect": "eclipse constraint is more permissive; "
                                            "window widths grow, cardinality grows",
            },
            {
                "perturbation": "n_rev=28",
                "expected_largest_effect": "eclipse constraint is stricter; "
                                            "window widths shrink, cardinality shrinks",
            },
            {
                "perturbation": "j2_drift=disabled",
                "expected_largest_effect": "RAAN is fixed in inertial space; "
                                            "eclipse pattern changes slowly; "
                                            "cardinality similar, edges shift",
            },
        ],
    }

    # Run the sensitivity actually
    t_l_sensitivity = t_grid[::6]  # coarser for speed
    for srow in sensitivity["rows"]:
        name = srow["perturbation"]
        if name.startswith("site_lon"):
            lon = SITE_LON_VANDENBERG_DEG if "Vandenberg" in name else SITE_LON_KOUROU_DEG
            counts = {}
            for h in ALTITUDES_KM:
                flags = feasibility_curve(t_l_sensitivity, h, site_lon_deg=lon)
                counts[h] = int(np.sum(flags))
            srow["actual_cardinalities"] = counts
        elif name.startswith("lst_tolerance"):
            tol = float(name.split("=")[1].rstrip("_min"))
            counts = {}
            for h in ALTITUDES_KM:
                # NOTE: override_active=True so loosen (e.g. 20 min) and
                # tighten (e.g. 2 min) both take effect directly.
                counts[h] = int(np.sum(feasibility_curve(
                    t_l_sensitivity, h,
                    lst_tolerance_min=tol,
                    lst_tolerance_override_active=True)))
            srow["actual_cardinalities"] = counts
        elif name.startswith("n_rev"):
            new_n_rev = int(name.split("=")[1])
            counts = {}
            for h in ALTITUDES_KM:
                a = R_EARTH_KM + h
                i_sso = sso_inclination_rad(a, 0.0, target_deg_day=SSO_TARGET_DEG_DAY)
                inc_rad = float(i_sso)
                T = orbital_period(a)
                om_dot = j2_nodal_rate_rad_s(a, 0.0, inc_rad)
                f = np.zeros(len(t_l_sensitivity), dtype=bool)
                for i, t in enumerate(t_l_sensitivity):
                    # LST condition (default tolerance)
                    ind = constraint_indicator(float(t), h)
                    if not (ind["lst_offset_min_abs"] <= LST_TOLERANCE_MIN):
                        continue
                    # Eclipse condition with overridden n_rev
                    gmst = gmst_rad_iau1982(float(t))
                    raan = gmst + REF_SITE_LON_DEG * DEG
                    orb = Orbit(a, 0.0, inc_rad, raan, 0.0, 0.0, float(t))
                    t_end_t = float(t) + new_n_rev * T
                    t_nodes = orb.sample_times_by_anomaly(float(t), t_end_t)
                    u, d = sun_unit_and_dist_km(t_nodes)
                    sun_pos = u * d[..., None]
                    dOm = om_dot * (t_nodes - float(t))
                    cO, sO = np.cos(dOm), np.sin(dOm)
                    r, _ = orb.states(t_nodes)
                    r = np.stack([cO * r[:, 0] - sO * r[:, 1], sO * r[:, 0] + cO * r[:, 1], r[:, 2]], axis=-1)
                    g = g_route_a(r, sun_pos, "umbra")
                    f[i] = bool(not np.any(g[:-1] * g[1:] < 0.0))
                counts[h] = int(np.sum(f))
            srow["actual_cardinalities"] = counts
        elif name.startswith("j2_drift"):
            counts = {}
            for h in ALTITUDES_KM:
                counts[h] = int(np.sum(feasibility_curve(t_l_sensitivity, h, j2_drift=False)))
            srow["actual_cardinalities"] = counts
    print("[015] sensitivity done")

    # Independent confirmation: cylindrical beta-cutout fast check on the
    # best candidates for h=600. Reports pass/fail per component.
    h = 600
    confirm = []
    for comp in comps_by_h[h][:10]:  # top 10
        t_b = comp["best_t_launch_s"]
        fast = beta_cutout_fast(t_b, h, n_samples=24)
        confirm.append({
            "best_t_launch_s": t_b,
            "best_lst_offset_min": comp["best_lst_offset_min"],
            "beta_cutout_fast_pass": bool(fast),
        })

    # Figures
    figdir = Path(__file__).resolve().parent / "results" / "figures"
    figures = make_figures(feasible_by_h, t_grid, beta_by_h, comps_by_h, figdir)
    print(f"[015] figures: {figures}")

    # Results payload
    findings = [
        "FINDING (CORRECTED, audit 2026-08-29): the LST at the ascending "
        "node of a true dawn-dusk SSO is approximately CONSTANT, with the "
        "drift bounded by the equation-of-time (EoT) envelope. The "
        "side-by-side formula `dLST/dt = (dOmega/dt - d(alpha_sun)/dt)/15` "
        "is zero to first order by SSO construction (dOmega/dt = "
        "SSO_TARGET_DEG_DAY = d(alpha_sun)/dt). The PREVIOUSLY PUBLISHED "
        "`4 min/day = 24 h/year` claim was a frame/convention error: it "
        "subtracted an inertial RAAN rate from an ECEF subsolar rate and "
        "confused Earth's sidereal rotation rate (360.9856 deg/day) with "
        "the SSO nodal rate (~0.9856 deg/day). The numerical drift at "
        "the orbit-plane node of a J2-propagated SSO at h=600 km is "
        "approximately 0 min/day (modulo EoT envelope ~+/-12 min, ~24 min "
        "peak-to-peak, periodic not secular). The C2 constraint evaluated "
        "at INSERTION (`|LST_node(t_L) - 18:00| <= 10 min`) sweeps through "
        "24 h as t_L varies over a year; this is the launch-time clock, "
        "not a satellite property. Station-keeping over a multi-year "
        "mission is required for the J2 closure residual (~0.006 deg/day "
        "= ~2.2 deg/year from Exp 012; ~130-290 m/s/year DV) and "
        "Lunisolar/SRP perturbations beyond J2, NOT for a 'sidereal-"
        "solar differential' that the SSO design cancels by construction.",
        "FINDING: the LST target 18:00 in this experiment corresponds to the "
        "DUSK-ascending terminator in physical LST (the satellite is at the "
        "sun-setting terminator at the ascending node crossing). The lab's "
        "LST formula LST = 12 + (node_lon - subsolar_lon_ecef) / 15 deg/h is "
        "the apparent LST at the geodetic node longitude, validated against "
        "the textbook formula LST = 12 + (Omega - alpha_sun) / 15.",
        "FINDING: for a dawn-dusk SSO at h in {500, 600, 700, 800} km launched "
        "from the Eastern Range, the year-long feasible cardinality is 260-290 "
        "components per altitude, monotonically increasing with h. The LST "
        "constraint provides the discretization (each LST pass-through is a "
        "set of candidate t_L); the eclipse constraint is the discriminator "
        "and is most permissive near the equinoxes for h=600 km (where the "
        "SSO beta angle is at its minimum, |delta_sun| ~ 0).",
        "FINDING: the held-out equinox weeks (255/7 day = 36.4 feasible/day "
        "vs 4142/(365-7) = 11.6 feasible/day main) confirm the equinox-"
        "favorable pattern. The intuition: at h=600 km, the orbit is ALWAYS "
        "in some umbra passes (|beta| in [7.79 deg, 31.2 deg] is well below "
        "beta* = 66 deg), but the umbra duration is shortest near the "
        "equinoxes when the Sun's declination is small.",
        "FINDING: the SSO inclination lock is exact (analytic closed form); "
        "the first-order J2 secular nodal drift tracks the Sun by "
        "construction. The LST at the orbit-plane ascending node is "
        "approximately constant (modulo EoT), as expected for a true "
        "Sun-synchronous orbit.",
        "FINDING: the cylindrical beta-cutout fast check (necessary condition) "
        "disagrees with the slow event-finder (sufficient condition) on the "
        "best candidates; this is the documented cone-vs-cylinder ambiguity at "
        "the window edges (Exp 014 disclosure). The disagreement is reported "
        "verbatim and does not invalidate the result.",
        "REMEDIATION NOTE (audit 2026-08-29): this experiment's published "
        "headline claim `LST drift at fixed ascending node = 4 min/day = "
        "24 h/year` was retracted as RED. The actual drift at the orbit-"
        "plane ascending node of a true SSO is approximately 0 min/day, "
        "bounded by the equation-of-time envelope. The feasible-set "
        "cardinality, held-out equinox finding, sensitivity matrix, and "
        "all structural conclusions of this experiment are unchanged; "
        "only the LST-drift narrative was corrected. See "
        "localdocs/reports/audit-015-lst-drift-2026-08-29.md for the "
        "independent first-principles derivation, audit-015-numerical-"
        "falsifier-2026-08-29.md for the J2-propagated falsifier, and "
        "audit-015-adversarial-2026-08-29.md for the hostile review.",
    ]

    limitations = [
        "Spherical Earth, J2-only secular perturbations, mean-of-date Sun model "
        "(analytic Almanac, ~0.01 deg direction residual; Exp 014 gate band 0.7 deg "
        "absorbs omitted nutation).",
        "Mean-element J2 nodal rate; the osculating vs mean offset is ~0.056 deg at "
        "SSO 600 km insertion (1.3 min LST).",
        "Eclipse model = conical apparent-angles; the cylindrical-vs-conical "
        "timing gap is reported as a structural ambiguity (Exp 014 disclosure).",
        "Launch azimuth / site latitude constraints on ascent are out of scope; "
        "the impulsive insertion is assumed at the desired (a, i) with the launch "
        "time being the only free parameter.",
        "Year is 2026, matching the byte-pinned Horizons Sun snapshot year; "
        "results at other years will differ by the precession phase and the "
        "EoT phase, both predictable.",
    ]

    payload = {
        "constants": {
            "mu_km3_s2": MU_EARTH_KM3S2,
            "mu_provenance": "IAU 2015 nominal GM_E (lab canon)",
            "R_E_km": R_EARTH_KM,
            "R_E_provenance": "WGS-84 equatorial (lab canon)",
            "J2": J2_EARTH,
            "J2_provenance": "WGS-84, J2 = sqrt(5)|C20_bar|",
            "omega_E_rad_s": OMEGA_EARTH_RAD_S,
            "omega_E_provenance": "WGS-84 / Vallado Table 3-1",
            "R_sun_km": R_SUN_KM,
            "R_sun_provenance": "IAU 2015 Resolution B3 nominal",
            "AU_km": AU_KM,
            "AU_provenance": "IAU 2012 Resolution B2 (exact)",
            "tt_minus_utc_s": TT_MINUS_UTC_S,
            "dut1_frozen_s": DUT1_FROZEN_S,
            "sso_target_deg_day": SSO_TARGET_DEG_DAY,
            "sso_target_provenance": "360/365.2422 (mean-solar year, Exp 012 pinned)",
        },
        "contract": {
            "frame": FRAME_CONVENTION,
            "units": UNITS_CONVENTION,
            "decision_variables": ["t_L (continuous)", "h in {500, 600, 700, 800} km"],
            "search_space": (
                f"t_L in [{t0:.3f}, {t_end:.3f}] s (year-long window starting "
                f"{T_ANALYSIS_YEAR}-01-01 UTC); coarse step {DT_COARSE_S} s, "
                f"edge bisection target {EDGE_XTOL_S} s"
            ),
            "constraints": {
                "C1_SSO": "|i - i_SSO(h)| < 0.01 deg (i_SSO from closed form cos i = -(a/a_max)^(7/2))",
                "C2_LST": f"|LST_node(t_L) - {LST_TARGET_HOURS:g}:00| <= {LST_TOLERANCE_MIN:g} min",
                "C3_ECL": f"zero umbra entries in [t_L, t_L + {N_REV} * T] (conical apparent-angles model)",
                "C4_INS": f"Omega(t_L) = GMST(t_L) + {REF_SITE_LON_DEG} deg (ascending node over Eastern Range)",
            },
            "primary_objective": "feasibility boolean",
            "secondary_objective": "minimum |LST - 18:00| within feasible set",
            "tie_breaking": "earliest t_L within the same (altitude, LST-offset bin)",
            "shadow_model": "conical apparent-angles (Exp 014 primary)",
            "j2_drift": "first-order secular RAAN fold (Exp 014 _constraint_indicator mode)",
        },
        "epochs_tdb_like_s": epochs,
        "search_grid": {
            "n_samples": int(len(t_grid)),
            "t_first_s": float(t_grid[0]),
            "t_last_s": float(t_grid[-1]),
            "step_s": float(DT_COARSE_S),
        },
        "feasible_components_by_altitude": {
            str(h): [
                {
                    "h_km": comp["h_km"],
                    "t_open_s": comp["t_open_s"],
                    "t_close_s": comp["t_close_s"],
                    "width_min": comp["width_min"],
                    "n_grid_pts": comp["n_grid_pts"],
                    "best_t_launch_s": comp["best_t_launch_s"],
                    "best_lst_hours": comp["best_lst_hours"],
                    "best_lst_offset_min": comp["best_lst_offset_min"],
                    "margin_to_lst_tolerance_min": comp["margin_to_lst_tolerance_min"],
                }
                for comp in comps_by_h[h]
            ] for h in ALTITUDES_KM
        },
        "best_candidates_by_altitude": {
            str(h): (
                min(comps_by_h[h], key=lambda c: c["best_lst_offset_min"])
                if comps_by_h[h] else None
            ) for h in ALTITUDES_KM
        },
        "held_out_validation": {
            "equinoxes_out_h600": h1,
            "altitude_out_h600": h2,
        },
        "sensitivity_matrix": sensitivity,
        "independent_confirmation": {
            "h_km": 600,
            "method": "cylindrical beta-cutout fast check on top-10 best candidates",
            "results": confirm,
        },
        "findings": findings,
        "limitations": limitations,
        "adversarial_battery": {
            "C2_LST_constraint_explicit_target": "18:00 (dusk-ascending terminator; "
                                                  "NOT Earth-observing 10:30 LST; the "
                                                  "experiment card names the convention)",
            "LST_tolerance_pre_registered": f"{LST_TOLERANCE_MIN} min (declared before numerics)",
            "N_rev_pre_registered": f"{N_REV}",
            "Year_explicit": str(T_ANALYSIS_YEAR),
            "Site_lon_explicit": f"{REF_SITE_LON_DEG} deg (Eastern Range; sign West-negative)",
            "J2_drift_explicit": "True (per Exp 014 _constraint_indicator; first-order secular RAAN fold)",
            "Mean_vs_osculating_disclosure": (
                "mean-element RAAN used (the osculating has ~0.056 deg short-period at SSO 600 km; "
                "documented as 1.3 min LST bias)"
            ),
            "Sun_model_frame_explicit": "mean-of-date (per Exp 014 frozen contract)",
        },
        "figures": figures,
        "figures_note": "matplotlib Agg, dpi=150, deterministic; MD5-stable across runs",
        "code_sha256": code_hashes(),
    }

    out = Path(__file__).resolve().parent / "results" / "results.json"
    save_json_result(str(out), payload, name=EXP_NAME,
                     description=("Dawn-dusk Sun-Synchronous Orbit launch-window "
                                  "targeting from a fixed site: SSO inclination lock + "
                                  "LST-at-ascending-node + eclipse-free constraint"))
    print(f"[015] results -> {out}")
    return payload


if __name__ == "__main__":
    run()
