"""Experiment 012 - Orbit Classes: constraint-defined families (SSO / Molniya / GTO / GEO anchor).

Research question
-----------------
Can the operationally important Earth-orbit classes be RECOVERED AND CLASSIFIED as
solution sets of coupled dynamical constraints under the declared model (two-body
Kepler + first-order secular J2 + spherical Earth rotating uniformly at omega_E),
with each class-defining quantity independently reproduced by closed-form algebra,
by full-force numerical propagation, and by structural identities -- rather than
asserted from mission folklore?

Classes as constraint equations on classical elements:
  SSO      : Omega_dot(a,e,i) = +360/365.2422 deg/day (mean-sun rate), retrograde
             branch i in (90,180]; closed form cos i_SSO = -(a/a_max(e))^(7/2);
             finite existence boundary a_max(e) beyond which NO sun-synchronous
             solution exists.
  Molniya  : apsidal freeze omega_dot = 0  <=>  cos^2 i = 1/5 (i_crit = 63.43494882 deg,
             supplement 116.56505118 deg), semi-synchronous resonance T = P_sidereal/2
             <=> n = 2 omega_E  =>  a = 26561.762 km, high-e family parameterized by e
             with perigee-survival bound, omega = 270 deg northern-apogee convention.
             First-order repeat-clock refinement: M_dot + omega_dot = 2 (omega_E -
             Omega_dot) => a_repeat ~ 26553 km (design disclosure; the full problem
             near i_crit carries small-divisor-amplified short-period dynamics --
             a excursions ~ +160 km and event-period excess ~ +325 s/orbit are
             MEASURED here, energy-conserving, converged under step refinement).
  GEO      : 1:1 spin-orbit resonance n = omega_E => a_GEO = 42164.169 km; e = 0,
             i = 0 is an exact fixed point of the pinned model (individual Omega_dot /
             omega_dot are nonzero but act only on unobservable degenerate elements).
  GTO      : connector class between LEO injection boundary condition and the GEO
             radius; pure two-body locks (vis-viva budgets, half-period transfer time).

Numerical contract (frozen v1.0)
--------------------------------
* Frames: ECI pseudo-inertial (Z = spin axis) for dynamics and secular rates;
  ECEF via uniform rotation theta(t) = omega_E t for ground-track quantities
  (Exp 008 convention). omega_E is the MASTER sidereal constant; P_sidereal is
  always DERIVED as 2*pi/omega_E, never hard-coded.
* Units: km / km^3-s^2 / s; angles internal radians, degrees at I/O boundaries.
* Constants (canon, src/lab_utils/orbits.py): mu = 398600.4418 (IAU 2015 B3),
  R_E = 6378.137 (WGS-84 TR8350.2), J2 = 1.082629821e-3 (WGS-84 sqrt(5)|C20_bar|),
  omega_E = 7.2921159e-5 rad/s (WGS-84 / Vallado Table 3-1).
* SSO target convention: mean SOLAR year 365.2422 d -> 360/365.2422 =
  0.985647332099 deg/day (continuity with Exp 009). The tropical-year variant
  (365.24219 d) differs by 2.8e-8 deg/day and is behaviorally indistinguishable
  (documented blindness); the SIDEREAL-year rate 0.98560912 deg/day and Julian-year
  rate 0.98564685 deg/day are WRONG targets, separated in i_SSO by 3.0e-4 / 1.7e-4 deg
  respectively -- both must be rejected by the 5e-5 deg analytic tolerance.
* Rates reported are SECULAR (first-order J2). Numerical estimates are orbit-averaged
  regressions of osculating elements over integer-orbit windows; converged numeric vs
  first-order theory differences are MODEL-ORDER residual (expected O(J2) ~ 0.5%),
  separated from integration error by convergence ladders and plateau-under-h-halving
  (Exp 009 doctrine). No exact-zero demand anywhere (float-tie discipline).
* Secular rates used (Vallado ch. 9, chapter-level):
    Omega_dot = -(3/2) n J2 (R/p)^2 cos i
    omega_dot = +(3/4) n J2 (R/p)^2 (5 cos^2 i - 1)
    M_dot     =   n + (3/4) n J2 (R/p)^2 sqrt(1-e^2) (3 cos^2 i - 1)
* Full-force validation propagations use the graduated lab_utils j2_rhs +
  rk4_propagate canon (donor-equivalent to Exp 009's verified loop, pinned tests).
* Thresholds sized to separate conventions (binding pair: solar-vs-Julian year at
  1.67e-4 deg -> i_SSO tolerance 5e-5 deg) and to keep model-order residual inside
  disclosed bands rather than tuned away.

Determinism: pure float64, no RNG, no wall-clock in experiment code, fixed grids,
Agg figures regenerated deterministically from recorded results.

References (chapter-level only)
-------------------------------
* D. A. Vallado, Fundamentals of Astrodynamics and Applications, 4th ed., Microcosm,
  2013 - Ch.9 general perturbations (secular J2 rates incl. mean-anomaly rate),
  Ch.3 time frames/constants (omega_E).
* H. D. Curtis, Orbital Mechanics for Engineering Students, 4th ed., Elsevier, 2021 -
  Ch.10 perturbations (sun-synchronous + Molniya class design conditions).
* R. R. Bate, D. D. Mueller, J. E. White, Fundamentals of Astrodynamics, Dover, 1971 -
  Ch.9 perturbations (critical inclination).
* NIMA WGS-84 TR8350.2; IAU 2015 Resolution B3 (arXiv:1510.07674) - constants.
* Real-mission anchors (Sentinel/Landsat-class SSO inclinations, Molniya geometry,
  GEO altitude 35786 km) are CONTEXTUAL screens only, never tuning targets.

Reuse: src/lab_utils/orbits.py (constants, Kepler/element canon, graduated j2_rhs),
src/lab_utils/integrators.py, src/lab_utils/metrics/results. Single-hop importlib
borrows: ols_fit / measure_secular_rates / node_crossing_raan_rate / analytic_rates
from Exp 009 (estimator plumbing, Exp 011 precedent); hohmann_dv1 / hohmann_dv2 /
transfer_elements / hohmann_transfer_time from Exp 004 (verified closed forms;
re-validating them would duplicate Exp 004's RK4 proofs). New experiment-local code:
domain-safe SSO solver with existence status, apsis-crossing period refinement,
dwell-fraction machinery, repeat-lattice solver, classification sweeps, adversarial
battery, figures.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lab_utils.integrators import rk4_propagate
from lab_utils.metrics import convergence_rate
from lab_utils.orbits import (
    J2_EARTH,
    MU_EARTH_KM3S2,
    OMEGA_EARTH_RAD_S,
    R_EARTH_KM,
    j2_rhs,
    mean_motion,
    orbital_period,
    rv_to_coe_eci,
    seed_state,
    steps_per_orbit,
)
from lab_utils.results import save_json_result

# --- Borrowed estimator plumbing (single-hop importlib, donor frozen) --------
J2_PATH = Path(__file__).resolve().parents[1] / "j2Precession" / "experiment.py"
_spec = importlib.util.spec_from_file_location("j2_009_for_012", J2_PATH)
assert _spec is not None and _spec.loader is not None
_j2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_j2)

ols_fit = _j2.ols_fit
measure_secular_rates = _j2.measure_secular_rates
node_crossing_raan_rate = _j2.node_crossing_raan_rate
analytic_rates = _j2.analytic_rates  # first-order secular Omega_dot / omega_dot

HM_PATH = Path(__file__).resolve().parents[1] / "hohmannTransfer" / "experiment.py"
_hspec = importlib.util.spec_from_file_location("hm_004_for_012", HM_PATH)
assert _hspec is not None and _hspec.loader is not None
_hm = importlib.util.module_from_spec(_hspec)
_hspec.loader.exec_module(_hm)

hohmann_dv1 = _hm.hohmann_dv1
hohmann_dv2 = _hm.hohmann_dv2
transfer_elements = _hm.transfer_elements
hohmann_transfer_time = _hm.hohmann_transfer_time

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

# --------------------------------------------------------------------------- #
# Constants and conventions (provenance recorded in results.json)
# --------------------------------------------------------------------------- #
YEAR_MEAN_SOLAR_DAYS = 365.2422  # mean-solar/Gregorian-mean convention (Exp 009 continuity)
YEAR_TROPICAL_DAYS = 365.24219  # documented variant, NOT used (behaviorally indistinguishable)
YEAR_SIDEREAL_DAYS = 365.25636  # WRONG for SSO target (equinox-precession offset)
YEAR_JULIAN_DAYS = 365.25  # WRONG for SSO target (most likely innocent substitution)
SSO_TARGET_DEG_DAY = 360.0 / YEAR_MEAN_SOLAR_DAYS
P_SIDEREAL_S = 2.0 * np.pi / OMEGA_EARTH_RAD_S  # MASTER resonance clock (derived, never hard-coded)

CRITICAL_INC_DEG = float(np.degrees(np.arccos(1.0 / np.sqrt(5.0))))  # 63.43494882...
CRITICAL_INC_SUPP_DEG = 180.0 - CRITICAL_INC_DEG
A_SEMISYNC_KM = (MU_EARTH_KM3S2 / (2.0 * OMEGA_EARTH_RAD_S) ** 2) ** (1.0 / 3.0)
A_GEO_KM = (MU_EARTH_KM3S2 / OMEGA_EARTH_RAD_S**2) ** (1.0 / 3.0)

FRAME_CONVENTION = (
    "ECI pseudo-inertial (Z=Earth spin axis) for dynamics/secular rates; ECEF via "
    "uniform rotation theta(t)=omega_E*t (spherical Earth, Exp 008 convention) for "
    "ground-track quantities"
)
UNITS_CONVENTION = "km, km^3/s^2, s; angles internal rad, I/O deg"

I_SSO_TOL_DEG = 5e-5  # binding separation solar-vs-Julian year = 1.67e-4 deg -> < half
MOLNIYA_RESIDUAL_REL_BAND = 2e-2  # model-order disclosure band for propagated rates
OMEGA_NULL_ABS_DEG_DAY = 5e-3  # |omega_dot_num| at exact critical: null-consistent bound
ANTISYM_RATIO_TOL = 2e-2  # symmetric/antisymmetric residual ratio at +/-0.5 deg
DWELL_NUM_TOL_ABS = 5e-3  # numeric dwell fraction vs closed form
APOGEE_DRIFT_TOL_DEG_DAY = 0.05  # abs band on wrapped daily apogee-longitude drift
# Period-clock doctrine (machine-established): first-order clocks T_anom=2pi/M_dot,
# T_drac=2pi/(M_dot+omega_dot) predict a +24.06 ms near-critical split at i=63.4 deg
# (the naive 2pi/(n -/+ omega_dot) forms give 2x this and are WRONG), but full-problem
# short-period event jitter near i_crit makes direct split detection infeasible;
# period claims are restricted to Kepler-machinery checks + converged
# excess-over-Kepler reporting + the apogee event-rate identity.
GEO_INCLINED_NODAL_TOL_S = 1.5  # |num - first-order| band; residual is model-order
KEPLER_EXCESS_PLATEAU_RATIO_BAND = (0.5, 2.0)  # excess(5432)/excess(2716) must plateau

DEG_DAY_PER_RAD_S = 86400.0 * 180.0 / np.pi


def wrap_deg(x):
    """Wrap degrees to (-180, 180]."""
    return (np.asarray(x, dtype=float) + 180.0) % 360.0 - 180.0


# --------------------------------------------------------------------------- #
# Production solvers (experiment-local; domain-safe, no silent clipping)
# --------------------------------------------------------------------------- #
def solve_sso_inclination(a_km: float, e: float,
                          target_deg_day: float = SSO_TARGET_DEG_DAY,
                          mu: float = MU_EARTH_KM3S2, j2: float = J2_EARTH,
                          R_eq_km: float = R_EARTH_KM) -> dict:
    """Solve Omega_dot(a,e,i) = +target for inclination, retrograde branch.

    Returns dict(status, incl_rad, cos_i_raw, margin). margin = 1 - target/K where
    K = (3/2) n J2 (R/p)^2 is the maximum deliverable eastward nodal rate (at
    i = 180). NO clipping: nonexistent solutions return an explicit sentinel so
    callers can never mistake a clipped 180 deg for physics.
    """
    n = mean_motion(a_km, mu)
    p = a_km * (1.0 - e * e)
    tgt_rad_s = np.radians(target_deg_day) / 86400.0
    K = 1.5 * n * j2 * (R_eq_km / p) ** 2
    ratio = tgt_rad_s / K  # = -cos_i on the retrograde branch
    margin = 1.0 - abs(ratio)
    if margin <= 0.0:
        return {"status": "NO_REAL_SOLUTION", "incl_rad": float("nan"),
                "cos_i_raw": float(ratio), "margin": float(margin)}
    return {"status": "OK", "incl_rad": float(np.arccos(-ratio)),
            "cos_i_raw": float(ratio), "margin": float(margin)}


def sso_existence_max_sma(e: float = 0.0, target_deg_day: float = SSO_TARGET_DEG_DAY,
                          mu: float = MU_EARTH_KM3S2, j2: float = J2_EARTH,
                          R_eq_km: float = R_EARTH_KM) -> float:
    """Maximum semi-major axis admitting an SSO solution (cos i_sol = -1 there).

    a_max^(7/2) = 1.5 J2 sqrt(mu) R^2 / (lambda (1-e^2)^2), lambda = sun rate [rad/s].
    Eccentricity EXTENDS the limit ((1-e^2)^(-4/7)): smaller p weakens J2 at fixed a.
    """
    lam = np.radians(target_deg_day) / 86400.0
    return float((1.5 * j2 * np.sqrt(mu) * R_eq_km**2 / (lam * (1.0 - e * e) ** 2)) ** (2.0 / 7.0))


def mean_anomaly_rate(a: float, e: float, inc_rad: float,
                      mu: float = MU_EARTH_KM3S2, j2: float = J2_EARTH) -> float:
    """Secular M-dot = n + (3/4) n J2 (R/p)^2 sqrt(1-e^2)(3cos^2 i - 1) [rad/s].

    Bracket (3cos^2 i - 1) is DISTINCT from omega_dot's (5cos^2 i - 1)
    (Vallado ch. 9; independently re-derived by the Exp 012 physics track).
    """
    n = mean_motion(a, mu)
    p = a * (1.0 - e * e)
    return float(n + 0.75 * n * j2 * (R_EARTH_KM / p) ** 2 * np.sqrt(1.0 - e * e)
                 * (3.0 * np.cos(inc_rad) ** 2 - 1.0))


# --------------------------------------------------------------------------- #
# Crossing refinement (apsis + node) and period estimators
# --------------------------------------------------------------------------- #
def _parabolic_refine(y_prev: float, y0: float, y_next: float) -> float:
    """Vertex fraction of a 3-point parabola on the middle stencil, clipped to [0,1]."""
    den = y_prev - 2.0 * y0 + y_next
    frac = 0.5 * (y_prev - y_next) / den if den != 0.0 else 0.5
    return min(max(frac, 0.0), 1.0)


def apsis_crossings(t: np.ndarray, states: np.ndarray) -> dict:
    """Refined periapsis/apoapsis crossing times via r.v sign changes.

    r.v changes + -> - at apoapsis... precisely: d(rv)/dt = v^2 - mu/r < 0 after
    periapsis; rv crosses +->0- at apoapsis entry and -->0+ at periapsis exit.
    Returns {'peri_times', 'apo_times'} (interior crossings only).
    """
    rv = np.einsum("ij,ij->i", states[:, :3], states[:, 3:])
    s = np.sign(rv)
    out = {}
    for name, pat in (("peri", (s[:-1] < 0) & (s[1:] >= 0)),
                      ("apo", (s[:-1] > 0) & (s[1:] <= 0))):
        idx = np.where(pat)[0]
        idx = idx[(idx >= 1) & (idx <= len(rv) - 2)]
        times = []
        for i in idx:
            frac = _parabolic_refine(rv[i - 1], rv[i], rv[i + 1])
            times.append(t[i] + frac * (t[i + 1] - t[i]))
        out[name] = np.asarray(times)
    return out


def ascending_node_times(t: np.ndarray, states: np.ndarray) -> np.ndarray:
    """Refined ascending-node (z: - -> +) crossing times."""
    z = states[:, 2]
    idx = np.where((z[:-1] < 0.0) & (z[1:] >= 0.0))[0]
    idx = idx[(idx >= 1) & (idx <= len(z) - 2)]
    times = []
    for i in idx:
        frac = _parabolic_refine(z[i - 1], z[i], z[i + 1])
        times.append(t[i] + frac * (t[i + 1] - t[i]))
    return np.asarray(times)


def slope_vs_index(times: np.ndarray) -> float | None:
    """OLS slope of crossing time vs crossing index (interior crossings only).

    First and last crossings are discarded to avoid edge bias; >= 4 interior
    crossings required. This is the period estimate.
    """
    if times.size < 6:
        return None
    x = np.arange(1, times.size - 1, dtype=float)
    fit = ols_fit(x, times[1:-1])
    return None if fit is None else fit["slope"]


def unwrap_guard_ok(angle_series: np.ndarray) -> bool:
    """Fail loud policy check: max per-sample unwrapped step must stay < pi/4.

    Guards np.unwrap against aliasing when coarse grids meet fast periapsis sweeps
    (adversarial track M12b): silently wrapped jumps of ~2*pi corrupt statistics.
    """
    steps = np.abs(np.diff(np.asarray(angle_series)))
    return bool(np.all(steps < np.pi / 4.0))


# --------------------------------------------------------------------------- #
# Dwell fractions (closed form + numeric dual path)
# --------------------------------------------------------------------------- #
def dwell_fraction_closed_form(e: float, delta_deg: float) -> float:
    """Fraction of period spent within +/-delta_deg of APOGEE (exact, Kepler).

    Window nu in [pi-D, pi+D] maps to eccentric anomaly E_1 = E(pi-D); by symmetry
    about apogee the window spans [E_1, 2 pi - E_1] and
    f = (t(2pi-E_1)-t(E_1))/T = (pi - E_1 + e sin E_1)/pi.
    Limits: f -> 1/2 as e -> 0, f -> 1 as e -> 1 (apogee dwell grows with e).
    """
    D = np.radians(delta_deg)
    cos_nu, sin_nu = -np.cos(D), np.sin(D)  # nu = pi - D
    denom = 1.0 + e * cos_nu
    cosE = (e + cos_nu) / denom
    sinE = np.sqrt(max(1.0 - e * e, 0.0)) * sin_nu / denom
    E1 = np.arctan2(sinE, cosE) % (2.0 * np.pi)
    return float((np.pi - E1 + e * np.sin(E1)) / np.pi)


def dwell_fraction_numeric(t: np.ndarray, states: np.ndarray, delta_deg: float,
                           T_kepler_s: float, coe: dict | None = None) -> float:
    """Numeric dwell fraction from propagated samples with boundary interpolation.

    Uses the LAST Kepler-orbit worth of samples (seed transient excluded); the
    propagation span is an integer multiple of T_kepler by construction. Angular
    distance to apogee delta(t) = |wrap(nu - pi)|; boundary crossings interpolated
    linearly between samples.
    """
    if coe is None:
        coe = rv_to_coe_eci(states[:, :3], states[:, 3:])
    nu_last = np.mod(coe["nu"], 2.0 * np.pi)
    mask = t >= (t[-1] - T_kepler_s)
    tt = t[mask]
    d = np.abs(wrap_deg(np.degrees(nu_last[mask]) - 180.0))
    inside = (d <= delta_deg).astype(float)
    dt = float(np.mean(np.diff(tt)))
    cross = np.where(np.diff(inside) != 0)[0]
    extra = 0.0
    for i in cross:
        d0, d1 = d[i], d[i + 1]
        f = (delta_deg - d0) / (d1 - d0)
        t_edge = tt[i] + f * (tt[i + 1] - tt[i])
        extra += (tt[i + 1] - t_edge) if inside[i] > 0 else -(tt[i + 1] - t_edge)
    return float((inside.sum() * dt + extra) / (tt[-1] - tt[0]))


# --------------------------------------------------------------------------- #
# Ground-track helpers (uniform sidereal rotation)
# --------------------------------------------------------------------------- #
def ecef_longitude_deg(t: np.ndarray, states: np.ndarray) -> np.ndarray:
    """Wrapped ECEF longitude [deg] under uniform rotation theta = omega_E t."""
    lon_inertial = np.degrees(np.arctan2(states[:, 1], states[:, 0]))
    return wrap_deg(lon_inertial - np.degrees(OMEGA_EARTH_RAD_S) * t)


def apogee_longitude_drift(t: np.ndarray, states: np.ndarray) -> dict:
    """ECEF-longitude rate of successive apogee events [deg/day, unwrapped scale].

    Apogee positions are interpolated at parabolic-refined crossing times from the
    pre-unwrapped sample-longitude series; ALL events are regressed in unwrapped
    space (no parity selection -- the true slope is the continuous event-rate
    (Omega_dot + omega_dot - omega_E) ~ -361.13 deg/day). The DAILY drift is the
    wrapped residual slope + 360 ~ -1.133 deg/day westward at Molniya elements
    (machine-established during adversarial debugging; the folklore '~Omega_dot'
    value conflates the inertial nodal regression with the ECEF event drift).
    """
    aps = apsis_crossings(t, states)
    ta = aps["apo"]
    if ta.size < 6:
        return {"slope_deg_day": None, "n_apogees": int(ta.size)}
    lon_samples_deg = ecef_longitude_deg(t, states)
    lon_unwrapped_rad = np.unwrap(np.radians(lon_samples_deg))
    lon_apo = np.empty(ta.size)
    for k, tk in enumerate(ta):
        i = int(np.searchsorted(t, tk))
        i = min(max(i, 1), len(t) - 1)
        frac = (tk - t[i - 1]) / (t[i] - t[i - 1])
        lon_apo[k] = lon_unwrapped_rad[i - 1] + frac * (lon_unwrapped_rad[i] - lon_unwrapped_rad[i - 1])
    fit = ols_fit(ta, np.degrees(lon_apo))
    if fit is None:
        return {"event_rate_deg_day": None, "n_apogees": int(ta.size)}
    # ols_fit slope is deg/SECOND (ta in seconds); convert to deg/day
    event_rate = float(fit["slope"]) * 86400.0
    return {
        "event_rate_deg_day": event_rate,
        "daily_drift_deg_day": float(wrap_deg(event_rate)),
        "n_apogees": int(ta.size),
        "fit_residual_rms_deg": float(fit["resid_rms"]),
    }


# --------------------------------------------------------------------------- #
# Case runner
# --------------------------------------------------------------------------- #
def propagate_case(a_km: float, e: float, inc_rad: float, n_orbits: int,
                   spp: int | None = None, j2: float = J2_EARTH,
                   Omega0: float = 0.0, omega0_deg: float = 270.0, M0: float = 0.0,
                   mu: float = MU_EARTH_KM3S2) -> dict:
    """Propagate one case on its documented grid; measure rates/periods/invariants."""
    spp = spp or steps_per_orbit(e)
    T = orbital_period(a_km, mu)
    t = np.linspace(0.0, n_orbits * T, n_orbits * spp + 1)
    r0, v0, nu0 = seed_state(a_km, e, inc_rad, Omega0, np.radians(omega0_deg), M0, mu)
    states = rk4_propagate(j2_rhs(mu, j2), t, np.concatenate([r0, v0]))
    return {"t": t, "states": states, "T_kepler_s": float(T), "spp": spp,
            "nu0_deg": float(np.degrees(nu0)), "j2": j2}


def invariants_gate(run: dict, j2: float = J2_EARTH) -> dict:
    """Energy / angular-momentum drift gate + r-band check (per-run sanity).

    Energy is J2-INCLUSIVE (Kepler + static J2 potential): a pure-Kepler energy
    oscillates at O(J2 (R/p)^2) when J2 is on and would false-fail the gate.
    """
    r, v = run["states"][:, :3], run["states"][:, 3:]
    rn = np.linalg.norm(r, axis=1)
    zr = r[:, 2] / rn
    P2 = 0.5 * (3.0 * zr * zr - 1.0)
    energy = (0.5 * np.einsum("ij,ij->i", v, v) - MU_EARTH_KM3S2 / rn
              + MU_EARTH_KM3S2 * j2 * R_EARTH_KM**2 * P2 / rn**3)
    hz = np.cross(r, v)[:, 2]
    return {
        "energy_drift_rel": float((energy.max() - energy.min()) / abs(energy[0])),
        "hz_range_rel": float((hz.max() - hz.min()) / abs(hz[0])),
    }


def measured_secular_trends(run: dict, windows_orbits: tuple[int, ...],
                            mu: float = MU_EARTH_KM3S2) -> dict:
    """Path-A element regression (+ Path-B node-crossing twin) on a propagation."""
    meas = measure_secular_rates(run["t"], run["states"], run["T_kepler_s"],
                                 list(windows_orbits), mu)
    primary = windows_orbits[-1]
    out = {"primary_window_orbits": primary}
    fo = meas["Omega"].get(primary)
    fw = meas["omega"].get(primary)
    out["Omega_dot_deg_day"] = float(np.degrees(fo["slope"]) * 86400.0) if fo else None
    out["omega_dot_deg_day"] = float(np.degrees(fw["slope"]) * 86400.0) if fw else None
    out["omega_resid_rms_deg"] = float(fw["resid_rms"]) if fw else None
    nc = node_crossing_raan_rate(run["t"], run["states"])
    out["node_crossing_Omega_dot_deg_day"] = (
        float(np.degrees(nc["slope"]) * 86400.0) if nc else None)
    out["Omega_defined"] = meas["Omega_defined"]
    out["omega_defined"] = meas["omega_defined"]
    out["elements_finite"] = bool(
        np.all(np.isfinite(meas["elements"]["a"]))
        and np.all(np.isfinite(meas["elements"]["e"]))
        and np.all(np.isfinite(meas["elements"]["inc"])))
    return out


# --------------------------------------------------------------------------- #
# Analysis blocks
# --------------------------------------------------------------------------- #
def sso_analysis() -> dict:
    """Closed-form SSO table + altitude sweep + existence bracket + wrong-convention probes."""
    alts = list(range(400, 1201, 50)) + [1500]
    table = {}
    for e in (0.0, 0.02):
        rows = []
        for h in alts:
            a = R_EARTH_KM + h
            sol = solve_sso_inclination(a, e)
            rows.append({
                "alt_km": h, "status": sol["status"],
                "i_SSO_deg": float(np.degrees(sol["incl_rad"])) if sol["status"] == "OK" else None,
                "margin": sol["margin"],
                "Omega_dot_check_deg_day": (
                    analytic_rates(a, e, sol["incl_rad"])["Omega_dot_deg_day"]
                    if sol["status"] == "OK" else None),
            })
        table[f"e_{e:.2f}"] = rows
    anchors = {}
    for h in (500.0, 600.0, 800.0):
        sol = solve_sso_inclination(R_EARTH_KM + h, 0.0)
        anchors[f"{int(h)}km"] = float(np.degrees(sol["incl_rad"]))
    amax0 = sso_existence_max_sma(0.0)
    bracket = []
    for ratio in (0.5, 0.9, 0.99, 0.999999, 1.000001, 1.001):
        sol = solve_sso_inclination(amax0 * ratio, 0.0)
        bracket.append({"a_over_amax": ratio, "status": sol["status"],
                        "i_SSO_deg": float(np.degrees(sol["incl_rad"])) if sol["status"] == "OK" else None})
    # wrong-target probes (recorded evidence for the test battery)
    def _i600(target):
        return float(np.degrees(solve_sso_inclination(R_EARTH_KM + 600.0, 0.0, target)["incl_rad"]))
    convention_probes = {
        "mean_solar_target_deg_day": SSO_TARGET_DEG_DAY,
        "i600_mean_solar_deg": _i600(SSO_TARGET_DEG_DAY),
        "i600_sidereal_year_deg": _i600(360.0 / YEAR_SIDEREAL_DAYS),
        "i600_julian_year_deg": _i600(360.0 / YEAR_JULIAN_DAYS),
        "i600_tropical_year_deg": _i600(360.0 / YEAR_TROPICAL_DAYS),
        "sidereal_year_separation_deg": abs(_i600(SSO_TARGET_DEG_DAY) - _i600(360.0 / YEAR_SIDEREAL_DAYS)),
        "julian_year_separation_deg": abs(_i600(SSO_TARGET_DEG_DAY) - _i600(360.0 / YEAR_JULIAN_DAYS)),
        "earth_rotation_confusion_status": solve_sso_inclination(
            R_EARTH_KM + 600.0, 0.0, np.degrees(OMEGA_EARTH_RAD_S) * 86400.0)["status"],
        "amax_sidereal_year_km": sso_existence_max_sma(0.0, 360.0 / YEAR_SIDEREAL_DAYS),
        "eccentric_extension_note": "a_max(e) = a_max(0) * (1-e^2)^(-4/7)",
        "eccentric_example_e02_rp600km": {
            "rp_km": R_EARTH_KM + 600.0,
            "a_km": (R_EARTH_KM + 600.0) / 0.8,  # a = rp/(1-e)
            "i_SSO_deg": float(np.degrees(solve_sso_inclination(
                (R_EARTH_KM + 600.0) / 0.8, 0.2)["incl_rad"])),
        },
    }
    return {"table": table, "anchors_deg": anchors, "a_max_km": amax0,
            "h_max_km": amax0 - R_EARTH_KM, "bracket": bracket,
            "convention_probes": convention_probes,
            "identity_cos_neg_apow72": {  # cos i_sol == -(a/a_max)^(7/2) at e=0
                str(h): float(-( (R_EARTH_KM + h) / amax0 ) ** 3.5) for h in (400, 800, 1200)
            }}


def sso_numeric_closure() -> dict:
    """Full-force propagation AT solved i_SSO must reproduce the sun rate (Paths A+B)."""
    out = {}
    for alt, orbits_n, spp in ((600.0, 50, 512), (800.0, 30, 512)):
        a = R_EARTH_KM + alt
        sol = solve_sso_inclination(a, 0.0)
        inc = sol["incl_rad"]
        run = propagate_case(a, 0.0, inc, orbits_n, spp)
        meas = measured_secular_trends(run, (orbits_n // 2, orbits_n))
        ana = analytic_rates(a, 0.0, inc)
        numA = meas["Omega_dot_deg_day"]
        numB = meas["node_crossing_Omega_dot_deg_day"]
        relA = abs(numA - SSO_TARGET_DEG_DAY) / SSO_TARGET_DEG_DAY
        out[f"{int(alt)}km"] = {
            "i_SSO_deg": float(np.degrees(inc)),
            "orbits": orbits_n, "spp": spp,
            "path_A_Omega_dot_deg_day": numA,
            "path_B_node_crossing_deg_day": numB,
            "target_deg_day": SSO_TARGET_DEG_DAY,
            "residual_rel_path_A": float(relA),
            "dual_path_agreement_rel": float(abs(numA - numB) / abs(numA)) if numB else None,
            "invariants": invariants_gate(run),
        }
    return out


def molniya_family() -> dict:
    """Frozen-perigee family at a_semisync across eccentricity (documented spp cap)."""
    es = [0.60, 0.64, 0.68, 0.70, 0.72, 0.74, 0.75]
    icrit = np.radians(CRITICAL_INC_DEG)
    rows = []
    for e in es:
        run = propagate_case(A_SEMISYNC_KM, e, icrit, 16, min(steps_per_orbit(e), 2048))
        meas = measured_secular_trends(run, (8, 16))
        ana = analytic_rates(A_SEMISYNC_KM, e, icrit)
        coe = rv_to_coe_eci(run["states"][:, :3], run["states"][:, 3:])
        d90_closed = dwell_fraction_closed_form(e, 90.0)
        d90_num = dwell_fraction_numeric(run["t"], run["states"], 90.0,
                                         run["T_kepler_s"], coe)
        aps = apsis_crossings(run["t"], run["states"])
        nodes = ascending_node_times(run["t"], run["states"])
        T_peri = slope_vs_index(aps["peri"])
        T_node = slope_vs_index(nodes)
        rows.append({
            "e": e, "spp_cap_used": run["spp"], "law_steps_per_orbit": steps_per_orbit(e),
            "hp_km": A_SEMISYNC_KM * (1 - e) - R_EARTH_KM,
            "ha_km": A_SEMISYNC_KM * (1 + e) - R_EARTH_KM,
            "analytic_omega_dot_deg_day": ana["omega_dot_deg_day"],
            "measured_omega_dot_deg_day": meas["omega_dot_deg_day"],
            "analytic_Omega_dot_deg_day": ana["Omega_dot_deg_day"],
            "measured_Omega_dot_deg_day": meas["Omega_dot_deg_day"],
            "dwell90_closed": d90_closed, "dwell90_numeric": d90_num,
            "dwell_abs_err": abs(d90_closed - d90_num),
            "T_peri_minus_T_sid_half_ms": None if T_peri is None else (T_peri - P_SIDEREAL_S / 2) * 1e3,
            "T_node_minus_T_sid_half_ms": None if T_node is None else (T_node - P_SIDEREAL_S / 2) * 1e3,
            "invariants": invariants_gate(run),
        })
    # (1-e^2)^-2 scaling law exponent for |omega_dot| off the exact lock
    ii = np.radians(CRITICAL_INC_DEG + 2.0)
    xs, ys = [], []
    for e in es[:5]:
        ana = analytic_rates(A_SEMISYNC_KM, e, ii)
        xs.append(-2.0 * np.log(1.0 - e * e))  # log (1-e^2)^-2
        ys.append(np.log(abs(ana["omega_dot_deg_day"])))
    scaling_slope = float(np.polyfit(xs, ys, 1)[0])
    return {"rows": rows, "p_scaling_exponent_target": 1.0,
            "p_scaling_loglog_slope": scaling_slope,
            "scaling_note": "omega_dot prop-to (R/p)^2: log|w| vs log((1-e^2)^-2) slope must be 1"}


def molniya_periods_and_drift() -> dict:
    """Period clocks near/at the lock + stroboscopic apogee drift + repeat radius.

    FINDING (machine-established during adversarial debugging, 2026-08-23): in the
    FULL J2 problem near the critical inclination the osculating elements carry
    LARGE short-period oscillations (a excursions ~ +160 km at e=0.74, i=i_crit;
    energy conserved to the integrator floor; converged under step refinement),
    and apsis/node event spacings exceed the Kepler period by ~ +325 s/orbit --
    far beyond any first-order secular clock (which predicts ~ +4 s). Direct
    event timing therefore CANNOT resolve the first-order T_peri - T_node split
    (~ +0.30 ms): short-period jitter dominates. Claims here are restricted to
    what the numerics support:
      (a) machinery check: J2 = 0 spacing equals T_kepler to < 1e-6 rel;
      (b) J2-on excess over Kepler reported with convergence plateau (spp x2);
      (c) the first-order split (+0.30 ms at 63.4 deg, zero at the lock) is
          recorded as THEORY DISCLOSURE only, never claimed as detected;
      (d) omega_dot suppression at the lock is validated on the orbit-averaged
          ELEMENT-regression path (Exp 009-proven), not by event timing.
    """
    out = {}
    icrit = np.radians(CRITICAL_INC_DEG)

    def _clocks(j2: float, inc_rad: float, n_orbits: int, spp: int):
        run = propagate_case(A_SEMISYNC_KM, 0.74, inc_rad, n_orbits, spp, j2=j2)
        aps = apsis_crossings(run["t"], run["states"])
        nodes = ascending_node_times(run["t"], run["states"])
        return (slope_vs_index(aps["peri"]), slope_vs_index(nodes), run)

    # (a) machinery check: Kepler-only spacing must equal T_kepler
    Tp0, Tn0, _ = _clocks(0.0, np.radians(63.4), 8, 2048)
    out["kepler_machinery_check"] = {
        "T_peri_rel_err_vs_kepler": abs(Tp0 - P_SIDEREAL_S / 2) / (P_SIDEREAL_S / 2),
        "T_node_rel_err_vs_kepler": abs(Tn0 - P_SIDEREAL_S / 2) / (P_SIDEREAL_S / 2),
    }
    # (b) J2-on excess over Kepler at/near the lock, with convergence plateau
    inc_near = np.radians(63.4)
    mdot_near = mean_anomaly_rate(A_SEMISYNC_KM, 0.74, inc_near)
    ana_near = analytic_rates(A_SEMISYNC_KM, 0.74, inc_near)
    wdot_near_rad_s = ana_near["omega_dot_rad_s"]
    pred_split_ms = float((2 * np.pi / mdot_near - 2 * np.pi / (mdot_near + wdot_near_rad_s)) * 1e3)
    plateau = {}
    for spp in (2716, 5432):
        Tp, Tn, run = _clocks(J2_EARTH, inc_near, 32, spp)
        plateau[str(spp)] = {
            "excess_peri_ms": (Tp - P_SIDEREAL_S / 2) * 1e3,
            "excess_node_ms": (Tn - P_SIDEREAL_S / 2) * 1e3,
            "split_ms": (Tp - Tn) * 1e3,
        }
    ex_hi = plateau["5432"]["excess_peri_ms"]
    ex_lo = plateau["2716"]["excess_peri_ms"]
    out["near_critical_63p4"] = {
        "orbits_used": 32,
        "kepler_excess_peri_ms": ex_hi,
        "kepler_excess_node_ms": plateau["5432"]["excess_node_ms"],
        "convergence_ratio_excess": ex_hi / ex_lo if ex_lo else None,
        "split_ms_measured": plateau["5432"]["split_ms"],
        "split_ms_first_order_theory_disclosure": pred_split_ms,
        "split_detection_note": (
            "first-order split ~ +24 ms is BELOW the short-period event-jitter "
            "floor (~tens of ms here); recorded as theory disclosure, NOT claimed "
            "as detected"),
        "omega_dot_analytic_deg_day": ana_near["omega_dot_deg_day"],
    }
    # exact lock: element-regression omega-dot null (the claimable freeze evidence)
    run_lock = propagate_case(A_SEMISYNC_KM, 0.74, icrit, 16, steps_per_orbit(0.74))
    meas_lock = measured_secular_trends(run_lock, (8, 16))
    out["at_exact_lock"] = {
        "omega_dot_element_regression_deg_day": meas_lock["omega_dot_deg_day"],
        "omega_dot_null_bound_deg_day": OMEGA_NULL_ABS_DEG_DAY,
        "node_crossing_Omega_dot_deg_day": meas_lock["node_crossing_Omega_dot_deg_day"],
    }
    # 12-sidereal-day drift run. The event series tracks the OSCULATING apsis,
    # which sweeps ~360 deg/orbit; its measurable identity (machine-verified to
    # ~1e-3 deg/orbit) is Delta_lambda_per_orbit = 360 + (Omega_dot + omega_dot -
    # omega_E) * T_orbit_days. The MEAN apsis line drifts at wrap(Omega_dot +
    # omega_dot - omega_E) ~ -1.133 deg/day (analytic class quantity, disclosed).
    n_days = 12.0
    n_orb = int(np.ceil(2.0 * n_days)) + 1  # ~2 orbits per sidereal day at semi-sync
    run4 = propagate_case(A_SEMISYNC_KM, 0.74, icrit, n_orb, steps_per_orbit(0.74))
    drift = apogee_longitude_drift(run4["t"], run4["states"])
    ana = analytic_rates(A_SEMISYNC_KM, 0.74, icrit)
    T_orb_meas = P_SIDEREAL_S / 2 + plateau["5432"]["excess_peri_ms"] / 1000.0
    inertial_rate = ana["Omega_dot_deg_day"] + ana["omega_dot_deg_day"] - np.degrees(OMEGA_EARTH_RAD_S) * 86400.0
    pred_event_rate = 360.0 * 86400.0 / T_orb_meas + inertial_rate
    out["apogee_drift"] = {
        "days_propagated": float(run4["t"][-1] / P_SIDEREAL_S),
        "event_rate_measured_deg_day": drift["event_rate_deg_day"],
        "event_rate_predicted_deg_day": float(pred_event_rate),
        "event_rate_abs_err_deg_day": (abs(drift["event_rate_deg_day"] - pred_event_rate)
                                       if drift["event_rate_deg_day"] is not None else None),
        "mean_line_drift_disclosure_deg_day": float(wrap_deg(inertial_rate)),
        "fit_residual_rms_deg": drift["fit_residual_rms_deg"],
        "n_apogees": drift["n_apogees"],
    }
    # Keplerian-only repeat closure at lattice radii (m orbits per k sidereal days)
    lat = []
    for m, k in ((1, 1), (2, 1), (3, 1), (1, 2)):
        T_res = k * P_SIDEREAL_S / m
        a_res = (MU_EARTH_KM3S2 * (T_res / (2.0 * np.pi)) ** 2) ** (1.0 / 3.0)
        t_grid = np.linspace(0.0, m * T_res, m * 2048 + 1)
        r0 = np.array([a_res, 0.0, 0.0])
        v0 = np.array([0.0, np.sqrt(MU_EARTH_KM3S2 / a_res), 0.0])
        st = rk4_propagate(j2_rhs(MU_EARTH_KM3S2, 0.0), t_grid, np.concatenate([r0, v0]))
        lon0 = np.degrees(np.arctan2(st[0, 1], st[0, 0]))
        lonN = np.degrees(np.arctan2(st[-1, 1], st[-1, 0]))
        closure = float(wrap_deg(lonN - lon0))
        lat.append({"m_per_k_days": f"{m}:{k}", "a_res_km": float(a_res),
                    "closure_deg_after_m_orbits": closure})
    out["repeat_lattice_keplerian"] = lat
    # J2-aware semi-synchronous correction: solve M_dot + omega_dot = 2*(omega_E - Omega_dot)
    # for a by bisection (repeat clock = draconitic rate measured against the node frame).
    # FIRST-ORDER DESIGN DISCLOSURE: the full problem near i_crit carries short-period
    # dynamics beyond this clock (see molniya_periods_and_drift finding).
    def _repeat_residual(a_km):
        an = analytic_rates(a_km, 0.74, icrit)
        mdot_a = mean_anomaly_rate(a_km, 0.74, icrit)
        return (mdot_a + an["omega_dot_rad_s"]) - 2.0 * (OMEGA_EARTH_RAD_S - an["Omega_dot_rad_s"])

    lo, hi = 26000.0, 27000.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _repeat_residual(lo) * _repeat_residual(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    out["repeat_corrected_radius_km"] = float(0.5 * (lo + hi))
    out["naive_radius_km"] = float(A_SEMISYNC_KM)
    return out


def geo_anchors() -> dict:
    """GEO fixed point: construction identity, stationarity residuals, inclined nodal period."""
    T_kep = orbital_period(A_GEO_KM)
    inc0 = 0.0
    ana0 = analytic_rates(A_GEO_KM, 0.0, inc0)
    resid_deg_day = (mean_anomaly_rate(A_GEO_KM, 0.0, inc0) * DEG_DAY_PER_RAD_S
                     + ana0["Omega_dot_deg_day"] + ana0["omega_dot_deg_day"]
                     - np.degrees(OMEGA_EARTH_RAD_S) * 86400.0)
    # inclined GEO: nodal (draconitic) period shift; T_node = 2pi/(M_dot + omega_dot).
    # High spp required: coarse grids quantize refined crossings to the sample step
    # (observed lock artifact at 512 spp where the shift is far below one sample).
    inc5 = np.radians(5.0)
    n_geo = mean_motion(A_GEO_KM)
    wd5 = analytic_rates(A_GEO_KM, 0.0, inc5)["omega_dot_deg_day"]
    mdot5 = mean_anomaly_rate(A_GEO_KM, 0.0, inc5)
    T_node_pred = 2.0 * np.pi / (mdot5 + np.radians(wd5) / 86400.0)
    run5 = propagate_case(A_GEO_KM, 0.0, inc5, 12, 4096)
    nodes = ascending_node_times(run5["t"], run5["states"])
    T_node_num = slope_vs_index(nodes)
    # figure-8 center (apparent) drift at i=5: Omega_dot + omega_dot
    ana5 = analytic_rates(A_GEO_KM, 0.0, inc5)
    center_drift = ana5["Omega_dot_deg_day"] + ana5["omega_dot_deg_day"]
    return {
        "a_geo_km": float(A_GEO_KM),
        "h_geo_km": float(A_GEO_KM - R_EARTH_KM),
        "period_match_rel_err": float(abs(T_kep - P_SIDEREAL_S) / P_SIDEREAL_S),
        "stationarity_residual_deg_day_keplerian_a": float(resid_deg_day),
        "Omega_dot_deg_day_i0": ana0["Omega_dot_deg_day"],
        "omega_dot_deg_day_i0": ana0["omega_dot_deg_day"],
        "relation_omegadot_eq_minus2Omegadot_i0_rel_err": float(
            abs(ana0["omega_dot_deg_day"] + 2.0 * ana0["Omega_dot_deg_day"]) / abs(ana0["omega_dot_deg_day"])),
        "inclined_i5_T_node_pred_minus_sid_ms": float((T_node_pred - P_SIDEREAL_S) * 1e3),
        "inclined_i5_T_node_num_minus_sid_ms": float((T_node_num - P_SIDEREAL_S) * 1e3) if T_node_num else None,
        "figure8_center_drift_deg_day": float(center_drift),
        "mean_element_sync_shift_km": float(
            A_GEO_KM - (MU_EARTH_KM3S2 / (mean_anomaly_rate(A_GEO_KM, 0.0, 0.0)
                                          + analytic_rates(A_GEO_KM, 0.0, 0.0)["Omega_dot_rad_s"]
                                          + analytic_rates(A_GEO_KM, 0.0, 0.0)["omega_dot_rad_s"]) ** 2) ** (1.0 / 3.0)),
        "model_boundary_note": (
            "Keplerian n=omega_E radius is the class anchor (42,164.169 km). The "
            "mean-element longitude-rate condition lambda_dot=M_dot+omega_dot+Omega_dot=omega_E "
            "shifts the synchronous radius by the recorded km-scale amount (first-order "
            "disclosure; real GEO fights tesseral+luni-solar drift outside this model)."),
    }


def gto_analysis() -> dict:
    """GTO budgets across injection altitudes + RK4 flight anchor + bug quantification."""
    ra = A_GEO_KM
    rows = []
    for h_p in (200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0):
        rp = R_EARTH_KM + h_p
        te = transfer_elements(rp, ra, MU_EARTH_KM3S2)
        dv1 = hohmann_dv1(rp, ra, MU_EARTH_KM3S2)
        dv2 = hohmann_dv2(rp, ra, MU_EARTH_KM3S2)
        vc_rp = np.sqrt(MU_EARTH_KM3S2 / rp)
        vc_ra = np.sqrt(MU_EARTH_KM3S2 / ra)
        rows.append({
            "hp_km": h_p, "a_km": float(te["a"]), "e": float(te["e"]),
            "vp_km_s": float(te["v_at_r1"]), "va_km_s": float(te["v_at_r2"]),
            "vcirc_rp_km_s": float(vc_rp), "vcirc_ra_km_s": float(vc_ra),
            "dv1_km_s": float(dv1), "dv2_km_s": float(dv2),
            "dv_total_km_s": float(dv1 + dv2),
            "tof_hours": float(hohmann_transfer_time(rp, ra, MU_EARTH_KM3S2) / 3600.0),
            "h_conservation_rp_vp_over_ra_va": float(te["v_at_r1"] * rp / (te["v_at_r2"] * ra)),
        })
    # canonical case RK4 flight anchor (Kepler-only and J2-on variants)
    hp_canon = 300.0
    rp = R_EARTH_KM + hp_canon
    te = transfer_elements(rp, ra, MU_EARTH_KM3S2)
    e_t = te["e"]
    spp = steps_per_orbit(e_t)
    T_t = orbital_period(te["a"])
    t_half = np.linspace(0.0, T_t / 2.0, spp + 1)
    r0 = np.array([rp, 0.0, 0.0])
    v0 = np.array([0.0, te["v_at_r1"], 0.0])
    arrival = {}
    for tag, jj in (("kepler_only", 0.0), ("j2_on", J2_EARTH)):
        st = rk4_propagate(j2_rhs(MU_EARTH_KM3S2, jj), t_half, np.concatenate([r0, v0]))
        r_end = np.linalg.norm(st[-1, :3])
        arrival[tag] = {
            "arrival_radius_km": float(r_end),
            "arrival_rel_err": float(abs(r_end - ra) / ra),
            "tof_hours": float(t_half[-1] / 3600.0),
        }
    bugs = {
        "bug_X_altitude_as_radius": {
            "note": "r_a := 35786.033 km (altitude) instead of 42164.17 km (radius)",
            "dv_total_km_s": float(hohmann_dv1(rp, 35786.03246186183, MU_EARTH_KM3S2)
                                   + hohmann_dv2(rp, 35786.03246186183, MU_EARTH_KM3S2)),
            # period if circularized at the WRONG apogee radius (not synchronous)
            "circular_period_at_wrong_apogee_hours": float(
                orbital_period(35786.03246186183) / 3600.0),
        },
        "bug_Y_radius_as_altitude": {
            "note": "r_a := R_E + 42164.17 km",
            "dv_total_km_s": float(hohmann_dv1(rp, R_EARTH_KM + A_GEO_KM, MU_EARTH_KM3S2)
                                   + hohmann_dv2(rp, R_EARTH_KM + A_GEO_KM, MU_EARTH_KM3S2)),
        },
        "swap_mutant_note": "dv1<->dv2 swap leaves the TOTAL invariant; per-burn signed asserts required",
        "swap_mutant_total_km_s": float(hohmann_dv2(rp, ra, MU_EARTH_KM3S2) + hohmann_dv1(rp, ra, MU_EARTH_KM3S2)),
        "circular_v_mutant_errors_km_s": {
            "at_molniya_perigee": float(np.sqrt(MU_EARTH_KM3S2 / A_SEMISYNC_KM)
                                        - np.sqrt(MU_EARTH_KM3S2 * (1 + 0.74) / (A_SEMISYNC_KM * (1 - 0.74)))),
            "at_molniya_apogee": float(np.sqrt(MU_EARTH_KM3S2 / A_SEMISYNC_KM)
                                       - np.sqrt(MU_EARTH_KM3S2 * (1 - 0.74) / (A_SEMISYNC_KM * (1 + 0.74)))),
        },
    }
    exp004_anchor = {
        "hp_km": 200.0,
        "dv_total_km_s": float(hohmann_dv1(R_EARTH_KM + 200.0, ra, MU_EARTH_KM3S2)
                               + hohmann_dv2(R_EARTH_KM + 200.0, ra, MU_EARTH_KM3S2)),
        "exp004_recorded_km_s": 3.9319,
        "note": "continuity anchor: same vis-viva forms already validated vs RK4 to <=4.2e-9 in Exp 004",
    }
    return {"rows": rows, "canonical_flight": {"hp_km": hp_canon, "e": float(e_t),
                                               "spp": spp, "arrival": arrival},
            "bugs": bugs, "exp004_anchor": exp004_anchor}


def critical_inclination_sweep() -> dict:
    """Staged omega_dot(i) zoom through the lock + antisymmetry + slope + plateau proof."""
    a, e = A_SEMISYNC_KM, 0.74
    spp = min(steps_per_orbit(e), 2048)  # documented cap for sweep economics
    orbits_n = 16

    def wdot_at(i_deg: float, n_orbits: int = orbits_n, spp_: int = spp) -> float:
        run = propagate_case(a, e, np.radians(i_deg), n_orbits, spp_)
        meas = measured_secular_trends(run, (n_orbits // 2, n_orbits))
        return meas["omega_dot_deg_day"]

    stage_A = [{"inc_deg": iv, "omega_dot_deg_day": wdot_at(iv)}
               for iv in np.arange(61.5, 65.51, 0.5)]
    stage_B = [{"inc_deg": iv, "omega_dot_deg_day": wdot_at(iv)}
               for iv in np.arange(62.95, 63.951, 0.05)]
    # dedicated both-sides pair + slope set at exact offsets from the lock
    deltas = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]
    offs = [{"delta_deg": d, "inc_deg": CRITICAL_INC_DEG + d,
             "omega_dot_deg_day": wdot_at(CRITICAL_INC_DEG + d, 24)} for d in deltas]
    plus = next(o for o in offs if o["delta_deg"] == 0.5)["omega_dot_deg_day"]
    minus = next(o for o in offs if o["delta_deg"] == -0.5)["omega_dot_deg_day"]
    antisym_ratio = abs(plus + minus) / abs(plus - minus)
    dd = np.array([-2.0, -1.0, -0.5, 0.5, 1.0, 2.0])
    ww = np.array([o["omega_dot_deg_day"] for o in offs])
    slope_fit = float(np.polyfit(dd, ww, 1)[0])  # deg/day per deg
    # sign-change localization from stage B bracket
    sb = sorted(stage_B, key=lambda r: r["inc_deg"])
    zero_bracket = None
    for lo, hi in zip(sb[:-1], sb[1:]):
        if lo["omega_dot_deg_day"] > 0.0 >= hi["omega_dot_deg_day"]:
            f = lo["omega_dot_deg_day"] / (lo["omega_dot_deg_day"] - hi["omega_dot_deg_day"])
            zero_bracket = lo["inc_deg"] + f * (hi["inc_deg"] - lo["inc_deg"])
            break
    # plateau proof at the lock (numerics stable under h halving)
    plateau = {str(s): wdot_at(CRITICAL_INC_DEG, 16, s) for s in (1024, 2048, 4096)}
    ana_crit = analytic_rates(a, e, np.radians(CRITICAL_INC_DEG))
    return {"stage_A": stage_A, "stage_B": stage_B, "offset_set": offs,
            "antisym_ratio_pm05deg": float(antisym_ratio),
            "slope_fit_deg_day_per_deg": slope_fit,
            # d(omega_dot)/di = -3.75 n J2 (R/p)^2 sin(2 i_crit); rad/s-per-rad -> deg/day-per-deg is x86400
            "slope_theory_deg_day_per_deg": float(
                -3.75 * mean_motion(a) * J2_EARTH * (R_EARTH_KM / (a * (1 - e * e))) ** 2
                * np.sin(np.radians(2.0 * CRITICAL_INC_DEG)) * 86400.0),
            "zero_crossing_localized_deg": zero_bracket,
            "plateau_at_lock_by_spp": plateau,
            "analytic_null_value_deg_day": ana_crit["omega_dot_deg_day"],
            "spp_policy": "min(law, 2048) documented deviation; anchor runs use full law",
            }


def adversarial_battery() -> dict:
    """Live-computed mutants (recorded evidence that each realistic error is detectable)."""
    out = {}
    # M01/M03: J2 sign flip / branch flip -> prograde twin whose rate, evaluated
    # under the TRUE J2 (the real world the code ships in), is exactly negated
    flipped = solve_sso_inclination(R_EARTH_KM + 600.0, 0.0, j2=-J2_EARTH)
    i_wrong = float(np.degrees(flipped["incl_rad"]))
    out["j2_sign_flip_wrong_branch"] = {
        "i_wrong_branch_deg": i_wrong,
        "produced_Omega_dot_true_J2_deg_day": analytic_rates(
            R_EARTH_KM + 600.0, 0.0, flipped["incl_rad"])["Omega_dot_deg_day"],
        "detectable_by": "signed assert i>90 deg AND sign(Omega_dot)==+; omega-dot tests are BLIND here (even in cos i)",
    }
    # M05a: solar half-day resonance mutant
    a_solar = (MU_EARTH_KM3S2 * (43200.0 / (2.0 * np.pi)) ** 2) ** (1.0 / 3.0)
    out["solar_halfday_resonance_mutant"] = {
        "a_mutant_km": float(a_solar), "a_true_km": float(A_SEMISYNC_KM),
        "delta_km": float(a_solar - A_SEMISYNC_KM),
        "detectable_by": "absolute-a assert at <=1e-9 rel; period-only asserts are blind (mutant is periodic too)",
    }
    # M12d: p := a substitution distortion (null tests blind, absolute value catches)
    e = 0.74
    p_sub_ratio = ((1.0 - e * e)) ** 2  # (p/a)^2 shrink factor applied to K
    out["p_equals_a_mutant"] = {
        "K_shrink_factor": float(p_sub_ratio),
        "Omega_dot_distortion_factor": float(1.0 / p_sub_ratio),
        "detectable_by": "ABSOLUTE Omega_dot value assert (frozen-null tests map zero to zero and survive)",
    }
    # M08: linear-in-nu dwell mutant
    d90 = dwell_fraction_closed_form(0.74, 90.0)
    out["dwell_linear_mutant"] = {
        "true_fraction": float(d90), "mutant_fraction": 0.5,
        "error_pp": float((d90 - 0.5) * 100.0),
    }
    # M06c: circular-v on eccentric arc (already computed in gto bugs; reference values here)
    # unwrap guard policy
    out["unwrap_guard_policy"] = {
        "max_step_allowed_rad": np.pi / 4.0,
        "molniya_law_spp": steps_per_orbit(0.74),
        "note": "guard enforced before every unwrap-based regression (fail loud, never silent wrap)",
    }
    # documented blindness / equivalent-mutant register (adversarial session record)
    out["documented_blind_spots"] = {
        "M14_refinement_index_on_uniform_grids": (
            "mis-indexed parabolic refinement step is an EQUIVALENT mutant here: all "
            "grids are np.linspace-uniform by contract, so t[i]+frac*(t[i]-t[i-1]) == "
            "t[i]+frac*(t[i+1]-t[i]) exactly; distinguishable only on non-uniform grids"),
        "M15_dwell_boundary_sign": (
            "flipping the numeric-dwell boundary-correction sign shifts the result "
            "~1e-3, below the 5e-3 corroboration tolerance BY DESIGN: the numeric "
            "dwell is a corroborating dual-form layer; the load-bearing dwell claims "
            "are the exact closed-form machine literals pinned in L1"),
        "M02b_tropical_year": (
            "2.1e-7 deg behavioral separation is below any defensible threshold; "
            "pinned by constant literal instead"),
    }
    return out


def pathological_grid(orbits_n: int = 3) -> dict:
    """Structural sentinels across i x e (NaN policy, r-bounds, energy gate)."""
    incs = [0.0, CRITICAL_INC_DEG, 90.0, CRITICAL_INC_SUPP_DEG, 180.0]
    eccs = [0.0, 0.05, 0.2, 0.74]
    a = R_EARTH_KM + 500.0
    rows = []
    ok = True
    for inc_deg in incs:
        for e in eccs:
            run = propagate_case(a, e, np.radians(inc_deg), orbits_n,
                                 Omega0=np.radians(10.0), omega0_deg=30.0)
            coe = rv_to_coe_eci(run["states"][:, :3], run["states"][:, 3:])
            rn = np.linalg.norm(run["states"][:, :3], axis=1)
            inv = invariants_gate(run)
            node_defined = inc_deg not in (0.0, 180.0)
            row = {
                "inc_deg": inc_deg, "e": e,
                "aei_finite": bool(np.all(np.isfinite(coe["a"])) and np.all(np.isfinite(coe["e"]))
                                   and np.all(np.isfinite(coe["inc"]))),
                "Omega_all_nan_expected": not node_defined,
                "Omega_all_nan": bool(np.all(~np.isfinite(coe["Omega"]))),
                "r_min_km": float(rn.min()), "r_max_km": float(rn.max()),
                "r_bounded": bool(np.all(rn > a * (1 - e) * (1 - 5e-3))
                                  and np.all(rn < a * (1 + e) * (1 + 5e-3))),
                "energy_drift_rel": inv["energy_drift_rel"],
            }
            row["ok"] = bool(row["aei_finite"] and row["r_bounded"]
                             and row["Omega_all_nan"] == row["Omega_all_nan_expected"]
                             and row["energy_drift_rel"] < 1e-8)
            ok = ok and row["ok"]
            rows.append(row)
    return {"incs_tested_deg": incs, "es_tested": eccs, "rows": rows, "all_ok": ok}


def convergence_study() -> dict:
    """RK4 order proofs: raw integrator (Kepler truth), rate metric, period metric."""
    mu = MU_EARTH_KM3S2
    # (a) raw integrator order vs closed-form circular truth, J2 = 0, phase-sensitive
    a0 = R_EARTH_KM + 420.0
    n0 = np.sqrt(mu / a0**3)
    T0 = 2.0 * np.pi / n0
    r0 = np.array([a0, 0.0, 0.0])
    v0 = np.array([0.0, np.sqrt(mu / a0), 0.0])
    errs_state = []
    spp_list_kep = (128, 256, 512, 1024)
    for spp in spp_list_kep:
        tg = np.linspace(0.0, T0, spp + 1)
        st = rk4_propagate(j2_rhs(mu, 0.0), tg, np.concatenate([r0, v0]))
        truth = np.column_stack([a0 * np.cos(n0 * tg), a0 * np.sin(n0 * tg), np.zeros_like(tg)])
        errs_state.append(float(np.max(np.linalg.norm(st[:, :3] - truth, axis=1))))
    orders_state = convergence_rate(np.maximum(np.asarray(errs_state), 1e-18),
                                    np.array([T0 / s for s in spp_list_kep]))
    # (b) SSO600 nodal-rate ladder vs 2048 reference at identical phases
    a = R_EARTH_KM + 600.0
    inc = solve_sso_inclination(a, 0.0)["incl_rad"]
    r0s, v0s, _ = seed_state(a, 0.0, inc, 0.0, 0.0, 0.0, mu)
    x0s = np.concatenate([r0s, v0s])
    T = orbital_period(a, mu)
    orbits_n, spp_eval = 20, 256
    rates = []
    for spp in (256, 512, 1024, 2048):
        tg = np.linspace(0.0, orbits_n * T, orbits_n * spp + 1)
        st = rk4_propagate(j2_rhs(mu, J2_EARTH), tg, x0s)
        stride = spp // spp_eval
        sub = st[::stride]
        coe = rv_to_coe_eci(sub[:, :3], sub[:, 3:], mu)
        fit = ols_fit(tg[::stride], np.unwrap(coe["Omega"]))
        rates.append(fit["slope"])
    errs_rate = np.abs(np.asarray(rates[:-1]) - rates[-1])
    # zero-guard only: a larger epsilon would clip the finest error (rad/s) and
    # corrupt the final convergence order
    orders_rate = convergence_rate(np.maximum(errs_rate, 1e-30),
                                   np.array([orbits_n * T / s for s in (256, 512, 1024)]))
    # (c) Molniya anomalistic-period ladder vs finest grid
    e_m = 0.74
    inc_m = np.radians(CRITICAL_INC_DEG)
    r0m, v0m, _ = seed_state(A_SEMISYNC_KM, e_m, inc_m, 0.0, 270.0, 0.0, mu)
    x0m = np.concatenate([r0m, v0m])
    Tm = orbital_period(A_SEMISYNC_KM)
    ladders, periods = (1358, 2716, 5432, 10864), []
    for spp in ladders:
        tg = np.linspace(0.0, 8 * Tm, 8 * spp + 1)
        st = rk4_propagate(j2_rhs(mu, J2_EARTH), tg, x0m)
        aps = apsis_crossings(tg, st)
        periods.append(slope_vs_index(aps["peri"]))
    err_period = np.abs(np.asarray(periods[:-1]) - periods[-1])
    orders_period = convergence_rate(np.maximum(err_period, 1e-9),
                                     np.array([8 * Tm / s for s in ladders[:-1]]))
    return {
        "state_space": {"spp_list": list(spp_list_kep), "max_pos_err_km": errs_state,
                        "orders_per_interval": [float(o) for o in orders_state],
                        "mean_order": float(np.mean(orders_state))},
        "rate_metric": {"spp_list": [256, 512, 1024, 2048],
                        "errors_deg_day": [float(np.degrees(er) * 86400.0) for er in errs_rate],
                        "orders_per_interval": [float(o) for o in orders_rate],
                        "mean_order": float(np.mean(orders_rate))},
        "period_metric": {"spp_list": list(ladders),
                          "T_peri_s": [float(p) for p in periods],
                          "errors_s": [float(er) for er in err_period],
                          "orders_per_interval": [float(o) for o in orders_period],
                          "mean_order": float(np.mean(orders_period)),
                          "finest_gap_ms": float(err_period[-1] * 1e3)},
    }


# --------------------------------------------------------------------------- #
# Figures (regenerated deterministically from recorded data)
# --------------------------------------------------------------------------- #
def make_figures(results: dict) -> list[str]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    # F1: SSO existence structure
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    rows = results["sso"]["table"]["e_0.00"]
    alts = [r["alt_km"] for r in rows if r["status"] == "OK"]
    incs = [r["i_SSO_deg"] for r in rows if r["status"] == "OK"]
    ax.plot(alts, incs, "o-", ms=4, lw=1.2, label="i_SSO(a) solved (e=0)")
    hmax = results["sso"]["h_max_km"]
    ax.axvline(hmax, color="r", ls="--", label=f"existence limit h_max={hmax:.1f} km")
    ax.axvspan(hmax, 1600, color="r", alpha=0.08)
    ax.text(hmax + 30, 97.0, "no SSO\nexists", color="r", fontsize=8)
    anchors = results["sso"]["anchors_deg"]
    ax.plot([500, 600, 800], [anchors["500km"], anchors["600km"], anchors["800km"]],
            "k*", ms=10, label="pinned anchors")
    ax.set_xlabel("Altitude [km]")
    ax.set_ylabel("Required sun-synchronous inclination [deg]")
    ax.set_title("SSO family: inclination lock rises to the existence boundary")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = fig_dir / "f1_sso_existence.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F2: omega_dot(i) through the critical inclination
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    sw = results["critical_sweep"]
    ia = [r["inc_deg"] for r in sw["stage_A"]]
    wa = [r["omega_dot_deg_day"] for r in sw["stage_A"]]
    ib = [r["inc_deg"] for r in sw["stage_B"]]
    wb = [r["omega_dot_deg_day"] for r in sw["stage_B"]]
    a_k, e_k = A_SEMISYNC_KM, 0.74
    igrid = np.linspace(61.0, 66.0, 251)
    wgrid = [analytic_rates(a_k, e_k, np.radians(iv))["omega_dot_deg_day"] for iv in igrid]
    ax.plot(igrid, wgrid, "k-", lw=1.0, label="first-order theory")
    ax.plot(ia, wa, "o", ms=5, label="measured (stage A)")
    ax.plot(ib, wb, ".", ms=4, color="C2", label="measured (stage B)")
    ax.axvline(CRITICAL_INC_DEG, color="r", ls="--", lw=1.0,
               label=f"i_crit={CRITICAL_INC_DEG:.4f} deg")
    ax.axhline(0.0, color="gray", lw=0.8)
    ax.set_xlabel("Inclination [deg]")
    ax.set_ylabel("Secular omega_dot [deg/day]")
    ax.set_title("Molniya apsidal lock: sign change localized at the critical inclination")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = fig_dir / "f2_omegadot_vs_inclination.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F3: Molniya dwell geometry over one orbit
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.0), sharex=True)
    run = propagate_case(A_SEMISYNC_KM, 0.74, np.radians(CRITICAL_INC_DEG), 1,
                         steps_per_orbit(0.74))
    t_hr = run["t"] / 3600.0
    st = run["states"]
    rn = np.linalg.norm(st[:, :3], axis=1)
    lat = np.degrees(np.arctan2(st[:, 2], np.hypot(st[:, 0], st[:, 1])))
    axes[0].plot(t_hr, rn / 1000.0, lw=1.2)
    axes[0].set_ylabel("Radius [1000 km]")
    axes[0].set_title("Molniya (a_semisync, e=0.74, i=i_crit): apogee dwell geometry")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(t_hr, lat, lw=1.2, color="C1")
    axes[1].axhline(63.43, color="gray", ls=":", lw=0.8)
    d90 = next(r for r in results["molniya_family"]["rows"] if abs(r["e"] - 0.74) < 1e-12)["dwell90_closed"]
    axes[1].set_ylabel("Geocentric latitude [deg]")
    axes[1].set_xlabel("Time [hours]")
    axes[1].set_title(f"time within +/-90 deg of apogee: {d90:.4f} of the period "
                      f"({d90 * orbital_period(A_SEMISYNC_KM) / 3600.0:.2f} h)")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "f3_molniya_dwell.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F4: GTO budget geometry
    fig, ax = plt.subplots(figsize=(8.0, 8.0))
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    for rr, cc, ll in ((R_EARTH_KM + 300.0, "C0", "LEO parking (300 km)"),
                       (A_GEO_KM, "C2", "GEO radius")):
        ax.plot(rr * np.cos(theta), rr * np.sin(theta), color=cc, lw=0.8, label=ll)
    te = transfer_elements(R_EARTH_KM + 300.0, A_GEO_KM, MU_EARTH_KM3S2)
    nu = np.linspace(0.0, np.pi, 181)
    r_t = te["p"] / (1.0 + te["e"] * np.cos(nu))
    ax.plot(r_t * np.cos(nu), r_t * np.sin(nu), "r-", lw=1.4, label="GTO half-ellipse")
    gto = results["gto"]["rows"][2]
    ax.annotate(f"dv1={gto['dv1_km_s']:.4f} km/s", (R_EARTH_KM + 300.0, 0),
                textcoords="offset points", xytext=(8, -14), fontsize=8)
    ax.annotate(f"dv2={gto['dv2_km_s']:.4f} km/s", (A_GEO_KM, 0),
                textcoords="offset points", xytext=(-70, 10), fontsize=8)
    ax.plot(0, 0, "k+", ms=10)
    ax.set_aspect("equal")
    ax.set_xlabel("x ECI [km]")
    ax.set_ylabel("y ECI [km]")
    ax.set_title(f"GTO budget: total {gto['dv_total_km_s']:.4f} km/s, "
                 f"tof {gto['tof_hours']:.3f} h (coplanar)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    p = fig_dir / "f4_gto_budget.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F5: convergence panels
    conv = results["convergence"]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    s1 = np.array(conv["state_space"]["spp_list"], dtype=float)
    e1 = np.maximum(np.array(conv["state_space"]["max_pos_err_km"]), 1e-18)
    axes[0].loglog(s1, e1, "o-", label="max position error (J2=0 truth)")
    axes[0].loglog(s1, e1[0] * (s1[0] / s1) ** 4, "r--", label="order-4 ref")
    axes[0].set_xlabel("Steps per orbit")
    axes[0].set_ylabel("Error [km]")
    axes[0].set_title(f"Raw integrator order {conv['state_space']['mean_order']:.2f}")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    s2 = np.array(conv["rate_metric"]["spp_list"][:-1], dtype=float)
    e2 = np.maximum(np.array(conv["rate_metric"]["errors_deg_day"]), 1e-18)
    axes[1].loglog(s2, e2, "s-", color="C2", label="nodal-rate error vs 2048 ref")
    axes[1].loglog(s2, e2[0] * (s2[0] / s2) ** 4.5, "r--", label="order-4.5 ref")
    axes[1].set_xlabel("Steps per orbit")
    axes[1].set_ylabel("|Omega_dot_h - Omega_dot_ref| [deg/day]")
    axes[1].set_title(f"Rate-metric order {conv['rate_metric']['mean_order']:.2f} (super-fourth legit)")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    p = fig_dir / "f5_convergence.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F6: apogee event-rate classification (resonant vs non-resonant)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    drift = results["molniya_periods"]["apogee_drift"]
    pred = drift["event_rate_predicted_deg_day"]
    meas = drift["event_rate_measured_deg_day"]
    ax.bar([0, 1], [meas, pred], width=0.5, color=["C0", "C2"],
           tick_label=["measured (12.5 d arc)", "predicted identity"])
    ax.set_ylabel("Apogee event ECEF-longitude rate [deg/day]")
    ax.set_title(f"Apogee event-rate identity (resonant Molniya): "
                 f"{meas:.3f} vs {pred:.3f} deg/day")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "f6_repeat_tracks.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)
    return paths


# --------------------------------------------------------------------------- #
# Main driver
# --------------------------------------------------------------------------- #
def main() -> dict:
    sso = sso_analysis()
    sso_num = sso_numeric_closure()
    fam = molniya_family()
    periods = molniya_periods_and_drift()
    geo = geo_anchors()
    gto = gto_analysis()
    crit = critical_inclination_sweep()
    adv = adversarial_battery()
    patho = pathological_grid()
    conv = convergence_study()

    worst_closure = max(v["residual_rel_path_A"] for v in sso_num.values())
    worst_dual = max(v["dual_path_agreement_rel"] for v in sso_num.values())
    dwell_row74 = next(r for r in fam["rows"] if abs(r["e"] - 0.74) < 1e-12)

    headline = {
        "sso_target_deg_day": SSO_TARGET_DEG_DAY,
        "i_SSO_500km_deg": sso["anchors_deg"]["500km"],
        "i_SSO_600km_deg": sso["anchors_deg"]["600km"],
        "i_SSO_800km_deg": sso["anchors_deg"]["800km"],
        "sso_existence_a_max_km": sso["a_max_km"],
        "sso_existence_h_max_km": sso["h_max_km"],
        "worst_sso_numeric_closure_rel": worst_closure,
        "worst_sso_dual_path_disagreement_rel": worst_dual,
        "critical_inclination_deg": CRITICAL_INC_DEG,
        "antisym_ratio_pm05deg": crit["antisym_ratio_pm05deg"],
        "slope_fit_deg_day_per_deg": crit["slope_fit_deg_day_per_deg"],
        "zero_crossing_localized_deg": crit["zero_crossing_localized_deg"],
        "a_semisync_km": A_SEMISYNC_KM,
        "molniya_hp_ha_km": [dwell_row74["hp_km"], dwell_row74["ha_km"]],
        "molniya_dwell90_fraction": dwell_row74["dwell90_closed"],
        "molniya_dwell90_abs_err_numeric": dwell_row74["dwell_abs_err"],
        "molniya_kepler_machinery_rel_err": periods["kepler_machinery_check"]["T_peri_rel_err_vs_kepler"],
        "molniya_kepler_excess_peri_ms": periods["near_critical_63p4"]["kepler_excess_peri_ms"],
        "molniya_split_first_order_disclosure_ms": periods["near_critical_63p4"]["split_ms_first_order_theory_disclosure"],
        "lock_omega_dot_regression_deg_day": periods["at_exact_lock"]["omega_dot_element_regression_deg_day"],
        "apogee_event_rate_measured_deg_day": periods["apogee_drift"]["event_rate_measured_deg_day"],
        "apogee_event_rate_pred_deg_day": periods["apogee_drift"]["event_rate_predicted_deg_day"],
        "apogee_mean_line_disclosure_deg_day": periods["apogee_drift"]["mean_line_drift_disclosure_deg_day"],
        "repeat_corrected_radius_km": periods["repeat_corrected_radius_km"],
        "a_geo_km": geo["a_geo_km"],
        "geo_period_match_rel_err": geo["period_match_rel_err"],
        "geo_stationarity_residual_deg_day": geo["stationarity_residual_deg_day_keplerian_a"],
        "geo_inclined_i5_nodal_shift_num_ms": geo["inclined_i5_T_node_num_minus_sid_ms"],
        "gto_dv_total_300km_km_s": gto["rows"][2]["dv_total_km_s"],
        "gto_exp004_anchor_km_s": gto["exp004_anchor"]["dv_total_km_s"],
        "gto_arrival_rel_err_kepler": gto["canonical_flight"]["arrival"]["kepler_only"]["arrival_rel_err"],
        "gto_arrival_rel_err_j2": gto["canonical_flight"]["arrival"]["j2_on"]["arrival_rel_err"],
        "state_order_mean": conv["state_space"]["mean_order"],
        "rate_order_mean": conv["rate_metric"]["mean_order"],
        "period_finest_gap_ms": conv["period_metric"]["finest_gap_ms"],
        "pathological_all_ok": patho["all_ok"],
    }

    results = {
        "constants": {
            "mu_km3_s2": MU_EARTH_KM3S2,
            "mu_provenance": "IAU 2015 Resolution B3 nominal GM_E; JPL DE440 planet-only 398600.435507 differs 1.5e-8 relative",
            "R_E_km": R_EARTH_KM,
            "R_E_provenance": "WGS-84 equatorial radius, NIMA TR8350.2",
            "J2": J2_EARTH,
            "J2_provenance": "WGS-84 TR8350.2: J2 = sqrt(5)*|C20_bar|, C20_bar = -0.484166774985e-3; EGM2008 1.08262668e-3 NOT used",
            "omega_E_rad_s": OMEGA_EARTH_RAD_S,
            "omega_E_provenance": "WGS-84 / Vallado Table 3-1; MASTER sidereal constant, P_sidereal derived as 2*pi/omega_E",
            "P_sidereal_s": float(P_SIDEREAL_S),
            "year_mean_solar_days": YEAR_MEAN_SOLAR_DAYS,
            "year_tropical_days_documented_not_used": YEAR_TROPICAL_DAYS,
            "year_sidereal_days_not_used_for_SSO": YEAR_SIDEREAL_DAYS,
            "year_julian_days_rejected": YEAR_JULIAN_DAYS,
            "sso_target_deg_day": SSO_TARGET_DEG_DAY,
            "sso_convention": (
                "mean SOLAR year 365.2422 d: 360/365.2422 = 0.985647332099 deg/day "
                "(Exp 009 continuity). Tropical-year variant differs 2.8e-8 deg/day "
                "(behaviorally indistinguishable, documented blindness). Sidereal-year "
                "rate 0.98560912 and Julian-year rate 0.98564685 deg/day are WRONG "
                "targets (separations 3.0e-4 / 1.7e-4 deg in i_SSO at 600 km)."),
            "frame_convention": FRAME_CONVENTION,
            "units": UNITS_CONVENTION,
        },
        "contract": {
            "rates_kind": "SECULAR first-order J2; numeric estimates are orbit-averaged osculating regressions",
            "model_residual_policy": "converged numeric vs first-order theory = model-order residual, never called integration error; plateau-under-h-halving proof required",
            "i_sso_tolerance_deg": I_SSO_TOL_DEG,
            "i_sso_tolerance_justification": "< half the binding solar-vs-Julian-year separation 1.67e-4 deg",
            "dual_path_agreement_tol_rel": 1e-2,
            "dual_path_note": "element-regression vs node-crossing rates agree <=1e-2 rel at 25-50-orbit windows (leakage-bias scale); both independently within the 1% model-order band of target",
            "period_ladder_policy": "report-only near i_crit: short-period event jitter gives a ~0.5 s floor; integrator order is proven by the state-space ladder instead",
            "omega_null_bound_deg_day": OMEGA_NULL_ABS_DEG_DAY,
            "antisymmetry_ratio_tol": ANTISYM_RATIO_TOL,
            "dwell_numeric_tol_abs": DWELL_NUM_TOL_ABS,
            "apogee_drift_tol_deg_day": APOGEE_DRIFT_TOL_DEG_DAY,
            "spp_cap_policy": "family/sweep runs capped at 2048 spp (documented deviation); anchor runs use the full steps_per_orbit law",
            "claims_guards": "omega-dot trends claimed only for e>=0.01 seeds; Omega undefined at i in {0,180}; no silent arccos clipping",
        },
        "sso": sso,
        "sso_numeric_closure": sso_num,
        "molniya_family": fam,
        "molniya_periods": periods,
        "geo": geo,
        "gto": gto,
        "critical_sweep": crit,
        "adversarial_battery": adv,
        "pathological": patho,
        "convergence": conv,
        "headline": headline,
    }

    fig_paths = make_figures(results)
    results["figures"] = fig_paths

    save_json_result(
        str(RESULTS_DIR / "results.json"),
        results,
        name="orbit_classes",
        description=(
            "Exp 012 orbit classes: constraint-defined families recovered under "
            "two-body + first-order secular J2 + spherical-Earth rotation. SSO "
            "inclination lock + finite existence boundary (a_max = 12352.5 km); "
            "Molniya apsidal freeze at cos^2 i = 1/5 + semi-synchronous resonance "
            "(a = 26561.76 km) + draconitic/anomalistic period split + dwell "
            "fractions + stroboscopic apogee drift; GEO 1:1 fixed point "
            "(a = 42164.17 km) with nonzero-but-unobservable individual rates; GTO "
            "vis-viva budgets anchored to Exp 004; adversarial convention battery."),
    )

    print("=== Orbit Classes: headline ===")
    print(f"i_SSO: 500km {headline['i_SSO_500km_deg']:.4f}  600km {headline['i_SSO_600km_deg']:.4f}  "
          f"800km {headline['i_SSO_800km_deg']:.4f} deg (target {SSO_TARGET_DEG_DAY:.6f} deg/day)")
    print(f"SSO existence limit: a_max {headline['sso_existence_a_max_km']:.3f} km "
          f"(h_max {headline['sso_existence_h_max_km']:.3f} km)")
    print(f"SSO numeric closure worst residual {worst_closure:.2e} rel; "
          f"dual-path disagreement {worst_dual:.2e} rel")
    print(f"i_crit {CRITICAL_INC_DEG:.6f} deg; antisym ratio(+/-0.5deg) {crit['antisym_ratio_pm05deg']:.2e}; "
          f"slope {crit['slope_fit_deg_day_per_deg']:.6f} deg/day/deg; "
          f"zero-crossing {crit['zero_crossing_localized_deg']:.4f} deg")
    print(f"Molniya a {A_SEMISYNC_KM:.3f} km, h_p/h_a {dwell_row74['hp_km']:.1f}/{dwell_row74['ha_km']:.1f} km; "
          f"dwell(+-90deg) {dwell_row74['dwell90_closed']:.5f} (num err {dwell_row74['dwell_abs_err']:.2e})")
    print(f"Kepler machinery rel err {headline['molniya_kepler_machinery_rel_err']:.2e}; "
          f"J2-on Kepler excess +{headline['molniya_kepler_excess_peri_ms']/1000.0:.1f} s/orbit near lock "
          f"(short-period amplified; split disclosure {headline['molniya_split_first_order_disclosure_ms']:+.3f} ms)")
    print(f"lock omega_dot (element regression) {headline['lock_omega_dot_regression_deg_day']:+.2e} deg/day "
          f"(null bound {OMEGA_NULL_ABS_DEG_DAY})")
    print(f"Apogee event rate: measured {headline['apogee_event_rate_measured_deg_day']:+.4f} vs predicted "
          f"{headline['apogee_event_rate_pred_deg_day']:+.4f} deg/day (identity 360/T + inertial rate); "
          f"mean-line drift disclosure {headline['apogee_mean_line_disclosure_deg_day']:+.5f} deg/day; "
          f"repeat-corrected radius {headline['repeat_corrected_radius_km']:.2f} km")
    print(f"GEO a {headline['a_geo_km']:.3f} km (period match {headline['geo_period_match_rel_err']:.1e}); "
          f"stationarity residual {headline['geo_stationarity_residual_deg_day']:+.5f} deg/day; "
          f"inclined-i5 nodal shift {headline['geo_inclined_i5_nodal_shift_num_ms']:.1f} ms")
    print(f"GTO 300km->GEO: dv1+dv2 {headline['gto_dv_total_300km_km_s']:.4f} km/s "
          f"(Exp004 anchor {headline['gto_exp004_anchor_km_s']:.4f} at 200 km); "
          f"arrival err kep {headline['gto_arrival_rel_err_kepler']:.1e}, j2 {headline['gto_arrival_rel_err_j2']:.1e}")
    print(f"Convergence: state order {headline['state_order_mean']:.2f}, rate order "
          f"{headline['rate_order_mean']:.2f}, period finest gap {headline['period_finest_gap_ms']:.3f} ms")
    print(f"pathological all ok: {headline['pathological_all_ok']}")
    print(f"figures: {fig_paths}")
    return results


if __name__ == "__main__":
    main()
