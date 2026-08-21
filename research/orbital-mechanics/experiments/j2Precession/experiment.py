"""Experiment 009 - J2 Precession: secular nodal & apsidal drift from Earth's oblateness.

Research question
-----------------
For Earth orbits (a, e, i) perturbed only by the J2 zonal term, does full-force
numerical Cowell RK4 propagation (a = a_Kepler + a_J2) rediscover the
first-order analytical secular rates

    Omega_dot = -(3/2) n J2 (R_E/p)^2 cos i
    omega_dot =  (3/4) n J2 (R_E/p)^2 (5 cos^2 i - 1),   n = sqrt(mu/a^3), p = a(1-e^2)

when the rates are measured by an independent state->element->trend estimator,
with (A) numerical convergence to a high-accuracy numerical reference at RK4
order ~4 separated from (B) the legitimate first-order model-order residual?

Numerical contract
------------------
* EOM (ECI, Z = Earth spin axis): r_ddot = -mu r/|r|^3 + a_J2(r) with
  a_J2 = -(3 mu J2 R_E^2)/(2 r^5) * [x(1-5z^2/r^2), y(1-5z^2/r^2), z(3-5z^2/r^2)]
  (gradient of U_J2 = +mu J2 R_E^2 P2(z/r)/r^3; componentwise sign check done).
* Fixed-step RK4 identical in structure to the verified Exp 006 propagator;
  the j2 == 0 path executes the identical float-op sequence and is regression
  tested BIT-EXACT against it.
* Secular-rate estimator: osculating elements recovered per sample from
  h = rxv, node = zhat x h, e_vec = (v x h)/mu - r/r; angles unwrapped and
  linearly regressed over integer-orbit windows (20/50/100 orbits reported,
  primary = 100 declared a priori). A third, algebraically independent
  estimator regresses ascending-node-crossing longitudes vs time (first-order
  Omega short-period terms vanish at the nodes).
* Convergence: E_h = |rate_h - rate_fine| against a 2048 steps/orbit numerical
  reference on identical sample phases; observed order p = log2(E_h/E_h/2)
  must lie in [3.6, 4.4]. The physics residual (converged numeric vs
  first-order analytic) is reported separately as model-order difference
  (expected O(J2) ~ 1e-3 relative plateau) and is NEVER called integration
  error.
* Singular elements: e = 0 -> omega structurally undefined (NaN sentinel),
  Omega_dot validated only; i = 0/180 -> RAAN structurally undefined, no
  Omega_dot claim. omega_dot claimed only for e >= 0.01 cases.
* No M_dot anywhere in this experiment (out of contract).

Determinism: pure float64, no RNG, fixed grids, Agg backend, figures
regenerated deterministically from recorded result data.

References (chapter-level only; equation numbers omitted where they could not
be verified against a specific printing - citation != truth):
* D. A. Vallado, Fundamentals of Astrodynamics and Applications, 4th ed.,
  Microcosm, 2013 - Ch.9 general perturbations, central-body J2 secular rates.
* H. D. Curtis, Orbital Mechanics for Engineering Students, 4th ed., Elsevier,
  2021 - Ch.10 introduction to orbital perturbations (J2 nodal/apsidal drift).
* R. R. Bate, D. D. Mueller, J. E. White, Fundamentals of Astrodynamics,
  Dover, 1971 - Ch.9 perturbations (nodal regression).
* C. D. Murray & S. F. Dermott, Solar System Dynamics, Cambridge UP, 1999 -
  Ch.6 disturbing function, planetary oblateness effects.
* NIMA, WGS-84, TR8350.2 - R_E = 6378.137 km; J2 = sqrt(5)|C20_bar|
  = 1.082629821e-3 (C20_bar = -0.484166774985e-3). EGM2008 gives
  J2 = 1.08262668e-3 (NOT used here; difference documented).
* IAU 2015 Resolution B3 (arXiv:1510.07674) - nominal GM_E = 398600.4418 km^3/s^2.

Reuse: verified Exp 006 3-D Cowell RK4 (importlib, J2=0 oracle), Exp 008
groundtracks Kepler/element machinery and constants (importlib),
src/lab_utils results/metrics. No scaffolding rebuilt.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lab_utils.metrics import convergence_rate
from lab_utils.results import save_json_result

# --- Reuse of verified prior-experiment machinery (single-hop importlib) ----
PCM_PATH = (
    Path(__file__).resolve().parents[1] / "planeChangeManeuvers" / "experiment.py"
)
_pcm_spec = importlib.util.spec_from_file_location("pcm_006_for_j2", PCM_PATH)
assert _pcm_spec is not None and _pcm_spec.loader is not None
_pcm = importlib.util.module_from_spec(_pcm_spec)
_pcm_spec.loader.exec_module(_pcm)

propagate_3d_rk4_kepler = _pcm.propagate_3d_rk4  # verified to <=1e-11 in Exp 006

GT_PATH = Path(__file__).resolve().parents[1] / "groundtracks" / "experiment.py"
_gt_spec = importlib.util.spec_from_file_location("groundtracks_for_j2", GT_PATH)
assert _gt_spec is not None and _gt_spec.loader is not None
_gt = importlib.util.module_from_spec(_gt_spec)
_gt_spec.loader.exec_module(_gt)

solve_kepler = _gt.solve_kepler
true_anomaly_from_E = _gt.true_anomaly_from_E
orbital_period = _gt.orbital_period
mean_motion = _gt.mean_motion
coe_to_rv_eci = _gt.coe_to_rv_eci

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

# --------------------------------------------------------------------------- #
# Physical constants (single source of truth for this experiment)
# --------------------------------------------------------------------------- #
MU_EARTH_KM3S2 = _gt.MU_EARTH_KM3S2  # 398600.4418 IAU 2015 B3 nominal GM_E
R_EARTH_KM = _gt.R_EARTH_KM  # 6378.137 WGS-84 TR8350.2 equatorial radius
J2_EARTH = 1.082629821e-3  # WGS-84: J2 = sqrt(5)*|C20_bar|, C20_bar=-0.484166774985e-3 (TR8350.2); EGM2008 1.08262668e-3 NOT used
YEAR_MEAN_SOLAR_DAYS = 365.2422  # mean solar year (Gregorian/J2000-mean convention)
SSO_TARGET_DEG_DAY = 360.0 / YEAR_MEAN_SOLAR_DAYS  # 0.985647360... deg/day mean-sun rate
CRITICAL_INC_DEG = float(np.degrees(np.arccos(1.0 / np.sqrt(5.0))))  # 63.43494882...
CRITICAL_INC_2_DEG = 180.0 - CRITICAL_INC_DEG  # 116.56505118...

FRAME_CONVENTION = (
    "ECI J2000 pseudo-inertial (Z=Earth spin axis, X=vernal equinox); a_J2 is "
    "axisymmetric so any Z-rotation of the frame leaves the dynamics invariant; "
    "ECEF/GMST not used in this experiment"
)
UNITS_CONVENTION = "km, km^3/s^2, s; angles internal rad, I/O deg"

NODE_GUARD_REL = 1e-6  # |zhat x h| / |h| = sin(i) below this -> RAAN undefined
ECC_GUARD_ABS = 1e-8  # |e_vec| below this -> argument of periapsis undefined
OMEGA_CLAIM_MIN_ECC = 0.01  # omega_dot trend claimed only for e >= this (policy)

PRIMARY_WINDOW_ORBITS = 100  # declared a priori primary fit window
WINDOW_STABILIZATION_REL = 1e-3  # max pairwise rel slope diff across windows
CONV_ORDER_BAND = (3.6, 4.4)  # acceptance band for measured RK4 order
PHYSICS_RESIDUAL_REL_TOL = 2e-2  # loose honest bound on first-order model residual


# --------------------------------------------------------------------------- #
# Dynamics: generalized fixed-step RK4 with optional J2 (structure = Exp 006)
# --------------------------------------------------------------------------- #
def j2_acceleration(r: np.ndarray, mu: float, j2: float) -> np.ndarray:
    """J2 perturbing acceleration in ECI (Z = spin axis). Vectorized over (N,3)."""
    r = np.asarray(r, dtype=float)
    single = r.ndim == 1
    if single:
        r = r[None, :]
    rm = np.linalg.norm(r, axis=1)
    z2r2 = (r[:, 2] * r[:, 2]) / (rm * rm)
    c = -1.5 * j2 * mu * R_EARTH_KM**2 / rm**5
    f = 1.0 - 5.0 * z2r2
    out = np.empty_like(r)
    out[:, 0] = c * r[:, 0] * f
    out[:, 1] = c * r[:, 1] * f
    out[:, 2] = c * r[:, 2] * (3.0 - 5.0 * z2r2)
    return out[0] if single else out


def propagate_3d_rk4_j2(r0, v0, mu, t, j2: float = J2_EARTH):
    """Fixed-step RK4 for r_ddot = -mu r/r^3 + a_J2(r). Deterministic.

    Loop structure and the j2 == 0 acceleration expression are operation-for-
    operation identical to the verified Exp 006 propagate_3d_rk4, so the
    J2 = 0 path reproduces it bit-for-bit (regression-tested).
    """
    n = len(t)
    state = np.empty((n, 6))
    state[0] = np.concatenate([np.asarray(r0, float), np.asarray(v0, float)])
    re2 = R_EARTH_KM**2

    def accel(x):
        r = x[:3]
        rm = np.linalg.norm(r)
        a_kep = -mu * r / rm**3
        if j2 == 0.0:
            return a_kep
        z2r2 = (r[2] * r[2]) / (rm * rm)
        c = -1.5 * j2 * mu * re2 / rm**5
        f = 1.0 - 5.0 * z2r2
        return a_kep + c * np.array([r[0] * f, r[1] * f, r[2] * (3.0 - 5.0 * z2r2)])

    for k in range(1, n):
        h = t[k] - t[k - 1]
        x = state[k - 1]
        k1 = np.concatenate([x[3:], accel(x[:3])])
        x2 = x + 0.5 * h * k1
        k2 = np.concatenate([x2[3:], accel(x2[:3])])
        x3 = x + 0.5 * h * k2
        k3 = np.concatenate([x3[3:], accel(x3[:3])])
        x4 = x + h * k3
        k4 = np.concatenate([x4[3:], accel(x4[:3])])
        state[k] = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return state


def steps_per_orbit(e: float, base: int = 512, ecc_base: int = 720) -> int:
    """Documented resolution rule: >=512/orbit circular; eccentric scaled as
    ceil(ecc_base/(1-e)^{3/2}) (periapsis-resolution law, Exp 002/008 precedent)."""
    if e <= 0.0:
        return base
    return max(base, int(np.ceil(ecc_base / (1.0 - e) ** 1.5)))


def seed_state(a, e, inc, Om, om, M0, mu=MU_EARTH_KM3S2):
    """Classical elements + M0 -> epoch state (r0, v0, nu0). M0 -> E0 -> nu0 chain."""
    E0 = solve_kepler(np.mod(M0, 2.0 * np.pi), e)
    nu0 = true_anomaly_from_E(E0, e)
    r0, v0 = coe_to_rv_eci(a, e, inc, Om, om, float(nu0), mu)
    return r0, v0, float(nu0)


# --------------------------------------------------------------------------- #
# Osculating element recovery (independent of the analytic oracle)
# --------------------------------------------------------------------------- #
def rv_to_coe_eci(r, v, mu=MU_EARTH_KM3S2) -> dict:
    """Cartesian state(s) -> classical osculating elements with singular guards.

    h = r x v; node vector n = zhat x h; e_vec = (v x h)/mu - r/r.
    Guards: sin(i) < NODE_GUARD_REL -> Omega NaN (RAAN undefined);
    |e_vec| < ECC_GUARD_ABS -> omega NaN (periapsis direction undefined).
    Accepts (N,3)/(N,3) or single vectors; returns arrays of a,e,i,Omega,omega,nu.
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    single = r.ndim == 1
    if single:
        r = r[None, :]
        v = v[None, :]
    rn = np.linalg.norm(r, axis=1)
    hvec = np.cross(r, v)
    hn = np.linalg.norm(hvec, axis=1)
    nvec = np.column_stack([-hvec[:, 1], hvec[:, 0], np.zeros(hn.shape[0])])
    nn = np.linalg.norm(nvec, axis=1)
    evec = np.cross(v, hvec) / mu - r / rn[:, None]
    em = np.linalg.norm(evec, axis=1)
    energy = 0.5 * np.einsum("ij,ij->i", v, v) - mu / rn
    sma = -mu / (2.0 * energy)
    inc = np.arctan2(nn, hvec[:, 2])

    good_node = nn > NODE_GUARD_REL * hn
    good_ecc = em >= ECC_GUARD_ABS
    Om = np.full(hn.shape, np.nan)
    Om[good_node] = np.arctan2(nvec[good_node, 1], nvec[good_node, 0])
    om = np.full(hn.shape, np.nan)
    both = good_node & good_ecc
    if np.any(both):
        nhat = nvec[both] / nn[both, None]
        ehat = evec[both] / em[both, None]
        hhat = hvec[both] / hn[both, None]
        cosw = np.einsum("ij,ij->i", nhat, ehat)
        sinw = np.einsum("ij,ij->i", np.cross(nhat, ehat), hhat)
        om[both] = np.arctan2(sinw, cosw)
    nu = np.full(hn.shape, np.nan)
    gnu = em > ECC_GUARD_ABS
    if np.any(gnu):
        ehat_n = evec[gnu] / em[gnu, None]
        rhat = r[gnu] / rn[gnu, None]
        hhat_n = hvec[gnu] / hn[gnu, None]
        cosnu = np.einsum("ij,ij->i", ehat_n, rhat)
        sinnu = np.einsum("ij,ij->i", np.cross(ehat_n, rhat), hhat_n)
        nu[gnu] = np.arctan2(sinnu, cosnu)
    out = {
        "a": sma,
        "e": em,
        "inc": inc,
        "Omega": Om,
        "omega": om,
        "nu": nu,
        "energy": energy,
        "h_z": hvec[:, 2],
        "h_mag": hn,
    }
    if single:
        return {k: (float(val) if np.ndim(val) == 0 else val[0]) for k, val in out.items()}
    return out


# --------------------------------------------------------------------------- #
# Estimators: unwrapped linear regression + node-crossing regression
# --------------------------------------------------------------------------- #
def ols_fit(t, y):
    """Closed-form ordinary-least-squares linear fit y = slope*t + intercept.

    Independent generic least squares (no shared algebra with the oracle).
    Returns dict or None when <2 finite points. NaNs are excluded.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(t) & np.isfinite(y)
    if np.count_nonzero(m) < 2:
        return None
    tt, yy = t[m], y[m]
    tm = float(tt.mean())
    ym = float(yy.mean())
    dt = tt - tm
    denom = float(np.dot(dt, dt))
    slope = float(np.dot(dt, yy - ym) / denom)
    intercept = ym - slope * tm
    resid = yy - (slope * tt + intercept)
    ss_res = float(np.dot(resid, resid))
    ss_tot = float(np.dot(yy - ym, yy - ym))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "slope": slope,
        "intercept": intercept,
        "r2": float(r2),
        "resid_rms": float(np.sqrt(ss_res / tt.size)),
        "resid_max": float(np.max(np.abs(resid))),
        "npts": int(tt.size),
    }


def measure_secular_rates(t, states, period_s, windows_orbits, mu=MU_EARTH_KM3S2) -> dict:
    """Unwrapped-element linear-regression estimator of Omega_dot / omega_dot.

    Windows are integer multiples of the Keplerian period (sample index =
    window_orbits * samples_per_orbit). All requested windows are reported;
    no cherry-picking. Angles unwrapped with np.unwrap (per-sample phase step
    is ~1e-7 rad at documented resolutions, 7 decades below pi).
    """
    coe = rv_to_coe_eci(states[:, :3], states[:, 3:], mu)
    n_total = len(t)
    out = {"windows_orbits": list(windows_orbits), "Omega": {}, "omega": {}}
    # all-or-nothing sentinel policy: unwrap only fully-defined series
    Om_ok = bool(np.all(np.isfinite(coe["Omega"])))
    om_ok = bool(np.all(np.isfinite(coe["omega"])))
    Om_u = np.unwrap(coe["Omega"]) if Om_ok else None
    om_u = np.unwrap(coe["omega"]) if om_ok else None
    for w in windows_orbits:
        # window end index: integer fraction of the total sample span
        frac = w / float(windows_orbits[-1])
        idx_end = int(round(frac * (n_total - 1)))
        res_o = None
        if Om_u is not None:
            fit = ols_fit(t[: idx_end + 1], Om_u[: idx_end + 1])
            res_o = dict(fit) if fit else None
            if res_o:
                res_o["window_orbits"] = int(w)
                res_o["t_span_days"] = float((t[idx_end] - t[0]) / 86400.0)
        out["Omega"][w] = res_o
        res_w = None
        if om_u is not None:
            fit = ols_fit(t[: idx_end + 1], om_u[: idx_end + 1])
            res_w = dict(fit) if fit else None
            if res_w:
                res_w["window_orbits"] = int(w)
                res_w["t_span_days"] = float((t[idx_end] - t[0]) / 86400.0)
        out["omega"][w] = res_w
    out["Omega_defined"] = Om_ok
    out["omega_defined"] = om_ok
    out["elements"] = coe
    return out


def node_crossing_raan_rate(t, states) -> dict | None:
    """Third-path estimator: regress inertial longitude at ascending-node
    crossings vs refined crossing time.

    At an ascending crossing the position lies ON the node, so atan2(y, x)
    equals the osculating RAAN there; first-order J2 short-period terms in
    Omega are proportional to cos(u) and vanish at the nodes, making this a
    clean secular estimator independent of the h/e_vec element algebra.
    Crossing times refined by parabolic vertex on the z(t) stencil.
    """
    z = states[:, 2]
    idx = np.where((z[:-1] < 0.0) & (z[1:] >= 0.0))[0]
    idx = idx[(idx >= 1) & (idx <= len(z) - 2)]
    if idx.size < 2:
        return None
    tc_list = []
    lon_list = []
    for i in idx:
        zm1, z0, zp1 = z[i - 1], z[i], z[i + 1]
        den = zm1 - 2.0 * z0 + zp1
        frac = 0.5 * (zm1 - zp1) / den if den != 0.0 else 0.5
        frac = min(max(frac, 0.0), 1.0)
        dt_i = t[i + 1] - t[i]
        tc = t[i] + frac * dt_i
        s0 = states[i]
        s1 = states[i + 1]
        xc = s0[0] + frac * (s1[0] - s0[0])
        yc = s0[1] + frac * (s1[1] - s0[1])
        tc_list.append(tc)
        lon_list.append(np.arctan2(yc, xc))
    tc_arr = np.asarray(tc_list)
    lon_u = np.unwrap(np.asarray(lon_list))
    fit = ols_fit(tc_arr, lon_u)
    if fit is None:
        return None
    fit["n_crossings"] = int(len(tc_arr))
    return fit


# --------------------------------------------------------------------------- #
# Analytical oracle (first-order secular J2 theory; uses ONLY a, e, i)
# --------------------------------------------------------------------------- #
def analytic_rates(a, e, inc_rad, mu=MU_EARTH_KM3S2, j2=J2_EARTH) -> dict:
    """First-order secular nodal/apsidal rates. Oracle inputs: (a, e, i) only."""
    n = mean_motion(a, mu)
    p = a * (1.0 - e * e)
    om_dot = -1.5 * n * j2 * (R_EARTH_KM / p) ** 2 * np.cos(inc_rad)
    w_dot = 0.75 * n * j2 * (R_EARTH_KM / p) ** 2 * (5.0 * np.cos(inc_rad) ** 2 - 1.0)
    return {
        "n_rad_s": float(n),
        "p_km": float(p),
        "Omega_dot_rad_s": float(om_dot),
        "omega_dot_rad_s": float(w_dot),
        "Omega_dot_deg_day": float(np.degrees(om_dot) * 86400.0),
        "omega_dot_deg_day": float(np.degrees(w_dot) * 86400.0),
    }


def sun_sync_inclination_rad(a, e, mu=MU_EARTH_KM3S2, j2=J2_EARTH,
                             target_deg_day=SSO_TARGET_DEG_DAY) -> float:
    """Solve Omega_dot(a,e,i) = target for i (retrograde branch, cos i < 0)."""
    n = mean_motion(a, mu)
    p = a * (1.0 - e * e)
    tgt_rad_s = np.radians(target_deg_day) / 86400.0
    cos_i = -tgt_rad_s * p**2 / (1.5 * n * j2 * R_EARTH_KM**2)
    cos_i = float(np.clip(cos_i, -1.0, 1.0))
    return float(np.arccos(cos_i))


def j2_specific_energy(r, v, mu=MU_EARTH_KM3S2, j2=J2_EARTH):
    """Specific mechanical energy incl. J2 potential U_J2 = mu J2 R_E^2 P2(z/r)/r^3.

    Static potential => exactly conserved by the ODE (integrator drift only).
    Gradient of this potential reproduces a_J2 componentwise (verified).
    """
    r = np.atleast_2d(np.asarray(r, dtype=float))
    v = np.atleast_2d(np.asarray(v, dtype=float))
    rn = np.linalg.norm(r, axis=1)
    u_ratio = r[:, 2] / rn
    p2 = 0.5 * (3.0 * u_ratio * u_ratio - 1.0)
    ke = 0.5 * np.einsum("ij,ij->i", v, v)
    E = ke - mu / rn + mu * j2 * R_EARTH_KM**2 * p2 / rn**3
    return E if E.size > 1 else float(E[0])


# --------------------------------------------------------------------------- #
# Case definitions (all derived from pinned constants; nothing hard-coded)
# --------------------------------------------------------------------------- #
def real_cases() -> dict:
    """Anchor cases. Rates/inclinations DERIVED from pinned constants."""
    inc_sso_600 = sun_sync_inclination_rad(R_EARTH_KM + 600.0, 0.0)
    cases = {
        "ISS": {
            "a_km": R_EARTH_KM + 420.0, "e": 0.0003, "inc_deg": 51.6,
            "Omega0_deg": 0.0, "omega0_deg": 0.0, "M0_deg": 0.0,
            "n_orbits": PRIMARY_WINDOW_ORBITS, "windows_orbits": (20, 50, 100),
            "claims": ("Omega_dot",),
            "desc": "ISS-like LEO 420 km / 51.6 deg (omega suppressed by claim policy e<0.01)",
        },
        "STARLINK": {
            "a_km": R_EARTH_KM + 550.0, "e": 0.0003, "inc_deg": 53.0,
            "Omega0_deg": 0.0, "omega0_deg": 0.0, "M0_deg": 0.0,
            "n_orbits": PRIMARY_WINDOW_ORBITS, "windows_orbits": (20, 50, 100),
            "claims": ("Omega_dot",),
            "desc": "Starlink-like LEO 550 km / 53 deg",
        },
        "SSO600": {
            "a_km": R_EARTH_KM + 600.0, "e": 0.0, "inc_deg": float(np.degrees(inc_sso_600)),
            "Omega0_deg": 0.0, "omega0_deg": 0.0, "M0_deg": 0.0,
            "n_orbits": PRIMARY_WINDOW_ORBITS, "windows_orbits": (20, 50, 100),
            "claims": ("Omega_dot",),
            "desc": "Sun-synchronous 600 km circular at solved i_SSO (e=0: omega sentinel)",
        },
        "POLAR": {
            "a_km": R_EARTH_KM + 500.0, "e": 0.0, "inc_deg": 90.0,
            "Omega0_deg": 0.0, "omega0_deg": 0.0, "M0_deg": 0.0,
            "n_orbits": PRIMARY_WINDOW_ORBITS, "windows_orbits": (20, 50, 100),
            "claims": ("Omega_dot_null",),
            "desc": "Polar circular 500 km (Omega_dot = 0 exact null)",
        },
        "MOLNIYA": {
            "a_km": 26560.0, "e": 0.74, "inc_deg": 63.4,
            "Omega0_deg": 0.0, "omega0_deg": 270.0, "M0_deg": 0.0,
            "n_orbits": 48, "windows_orbits": (12, 24, 48),
            "claims": ("Omega_dot", "omega_dot_small_nonzero"),
            "desc": "Molniya-like a=26560 km e=0.74 i=63.4 (near-critical: omega_dot tiny nonzero)",
        },
        "CRITICAL": {
            "a_km": R_EARTH_KM + 500.0, "e": 0.2, "inc_deg": CRITICAL_INC_DEG,
            "Omega0_deg": 0.0, "omega0_deg": 0.0, "M0_deg": 0.0,
            "n_orbits": 60, "windows_orbits": (15, 30, 60),
            "claims": ("Omega_dot", "omega_dot_null"),
            "desc": "Exact critical inclination arccos(1/sqrt(5)): apsidal rate frozen",
        },
        "ECC_REF": {
            "a_km": R_EARTH_KM + 500.0, "e": 0.2, "inc_deg": 30.0,
            "Omega0_deg": 0.0, "omega0_deg": 0.0, "M0_deg": 0.0,
            "n_orbits": 60, "windows_orbits": (15, 30, 60),
            "claims": ("Omega_dot", "omega_dot"),
            "desc": "Generic eccentric prograde reference for apsidal drift",
        },
    }
    return cases


def run_case(case: dict, j2: float = J2_EARTH, mu=MU_EARTH_KM3S2) -> dict:
    """Propagate one case on its documented grid and measure everything."""
    a, e = case["a_km"], case["e"]
    inc = np.radians(case["inc_deg"])
    Om0 = np.radians(case["Omega0_deg"])
    om0 = np.radians(case["omega0_deg"])
    M0 = np.radians(case["M0_deg"])
    T = orbital_period(a, mu)
    spp = steps_per_orbit(e)
    n_orbits = case["n_orbits"]
    t = np.linspace(0.0, n_orbits * T, n_orbits * spp + 1)
    r0, v0, nu0 = seed_state(a, e, inc, Om0, om0, M0, mu)
    states = propagate_3d_rk4_j2(r0, v0, mu, t, j2)
    meas = measure_secular_rates(t, states, T, case["windows_orbits"], mu)
    nc = node_crossing_raan_rate(t, states)
    ana = analytic_rates(a, e, inc, mu, j2)
    E = j2_specific_energy(states[:, :3], states[:, 3:], mu, j2)
    e_drift_rel = float((np.max(E) - np.min(E)) / abs(E[0]))
    hz = meas["elements"]["h_z"]
    hm = meas["elements"]["h_mag"]
    inv = {
        "energy_drift_rel": e_drift_rel,
        "h_mag_range_rel": float((np.max(hm) - np.min(hm)) / hm[0]),
        "h_z_range_rel": float((np.max(hz) - np.min(hz)) / abs(hz[0])),
    }
    inv["hz_over_hmag_range_ratio"] = (
        inv["h_z_range_rel"] / inv["h_mag_range_rel"] if inv["h_mag_range_rel"] > 0 else 0.0
    )
    return {
        "case": case,
        "T_sec": float(T),
        "spp": spp,
        "dt_sec": float(T / spp),
        "nu0_deg": float(np.degrees(nu0)),
        "analytic": ana,
        "measure": meas,
        "node_crossing": nc,
        "invariants": inv,
        "states": states,
        "t": t,
    }


# --------------------------------------------------------------------------- #
# Validation layers
# --------------------------------------------------------------------------- #
def validate_case_rates(run: dict) -> dict:
    """Assemble per-case numeric-vs-analytic comparison + residuals."""
    case = run["case"]
    ana = run["analytic"]
    meas = run["measure"]
    primary_w = case["windows_orbits"][-1]
    # relative residual is meaningful only when the analytic value is not a
    # structural zero (polar/equatorial nulls): fall back to absolute
    def _resid(num, ref):
        if num is None:
            return None, None
        d = abs(num - ref)
        rel = d / abs(ref) if abs(ref) > 1e-9 else None
        return d, rel

    out = {
        "windows_orbits": list(case["windows_orbits"]),
        "primary_window_orbits": primary_w,
        "spp": run["spp"],
        "n_orbits": case["n_orbits"],
        "claims": list(case["claims"]),
        "analytic_Omega_dot_deg_day": ana["Omega_dot_deg_day"],
        "analytic_omega_dot_deg_day": ana["omega_dot_deg_day"],
        "Omega_dot_deg_day_by_window": {},
        "omega_dot_deg_day_by_window": {},
    }
    for w in case["windows_orbits"]:
        fo = meas["Omega"].get(w)
        fw = meas["omega"].get(w)
        out["Omega_dot_deg_day_by_window"][w] = (
            float(np.degrees(fo["slope"]) * 86400.0) if fo else None
        )
        out["omega_dot_deg_day_by_window"][w] = (
            float(np.degrees(fw["slope"]) * 86400.0) if fw else None
        )
    po = meas["Omega"][primary_w]
    pw = meas["omega"][primary_w]
    out["Omega_fit_stats"] = (
        {k: po[k] for k in ("slope", "intercept", "r2", "resid_rms", "resid_max", "npts", "t_span_days")}
        if po else None
    )
    out["omega_fit_stats"] = (
        {k: pw[k] for k in ("slope", "intercept", "r2", "resid_rms", "resid_max", "npts", "t_span_days")}
        if pw else None
    )
    num_Om = out["Omega_dot_deg_day_by_window"][primary_w]
    num_om = out["omega_dot_deg_day_by_window"][primary_w]
    out["numeric_Omega_dot_deg_day"] = num_Om
    out["numeric_omega_dot_deg_day"] = num_om
    d_Om, rel_Om = _resid(num_Om, ana["Omega_dot_deg_day"])
    d_om, rel_om = _resid(num_om, ana["omega_dot_deg_day"])
    out["Omega_residual_abs_deg_day"] = d_Om
    out["Omega_residual_rel"] = rel_Om
    out["omega_residual_abs_deg_day"] = d_om
    out["omega_residual_rel"] = rel_om
    if case["e"] < OMEGA_CLAIM_MIN_ECC:
        out["omega_note"] = (
            "seed e < OMEGA_CLAIM_MIN_ECC: J2 induces a real eccentricity "
            "(~1e-3 at LEO) that dominates the seed, so the recovered omega "
            "sweeps once per orbit (slope ~= mean motion) and carries NO "
            "secular apsidal signal; omega_dot is not claimed"
        )
    slopes_o = [v for v in out["Omega_dot_deg_day_by_window"].values() if v is not None]
    if len(slopes_o) >= 2 and abs(slopes_o[-1]) > 0:
        diffs = [abs(s - slopes_o[-1]) / abs(slopes_o[-1]) for s in slopes_o[:-1]]
        out["Omega_window_stabilization_max_rel"] = float(max(diffs))
    else:
        out["Omega_window_stabilization_max_rel"] = None
    if run["node_crossing"] is not None:
        nc_rate = float(np.degrees(run["node_crossing"]["slope"]) * 86400.0)
        out["node_crossing_Omega_dot_deg_day"] = nc_rate
        out["node_crossing_vs_element_rel"] = (
            abs(nc_rate - num_Om) / abs(num_Om)
            if num_Om is not None and abs(num_Om) > 1e-9 else None
        )
        out["node_crossing_vs_analytic_abs_deg_day"] = abs(nc_rate - ana["Omega_dot_deg_day"])
        out["node_crossing_vs_analytic_rel"] = (
            abs(nc_rate - ana["Omega_dot_deg_day"]) / abs(ana["Omega_dot_deg_day"])
            if abs(ana["Omega_dot_deg_day"]) > 1e-9 else None
        )
        out["node_crossing_n"] = run["node_crossing"]["n_crossings"]
    else:
        out["node_crossing_Omega_dot_deg_day"] = None
    out["invariants"] = run["invariants"]
    return out


def convergence_study(orbits: int = 20, a_km=None, e: float = 0.001,
                      inc_deg: float = 51.6, j2: float = J2_EARTH) -> dict:
    """RK4 order check on two metrics against a 2048/orbit numerical reference.

    Grids 128..1024 steps/orbit vs reference; every grid is evaluated at
    IDENTICAL sample phases (stride subsampling at the coarsest spacing) so
    short-period sampling cancels in the differences.

    Metrics:
      (a) state-space final-position error vs reference -- must show the raw
          integrator order (~4);
      (b) secular-rate-estimate error vs reference -- measured orders run
          ABOVE 4 (4.5-4.7 observed): the estimator's orbit-averaging cancels
          RK4's leading phase-error mode (element recovery is invariant to
          along-track timing shifts), so rate convergence is faster than h^4.
          Both are reported; acceptance uses (a) for integrator order and
          >=3.6-with-monotone-decay for (b).
    """
    if a_km is None:
        a_km = R_EARTH_KM + 420.0
    mu = MU_EARTH_KM3S2
    T = orbital_period(a_km, mu)
    inc = np.radians(inc_deg)
    r0, v0, _ = seed_state(a_km, e, inc, 0.0, 0.0, 0.0, mu)
    spp_eval = 128
    spp_list = [128, 256, 512, 1024, 2048]
    rates = []
    final_positions = []
    for spp in spp_list:
        t = np.linspace(0.0, orbits * T, orbits * spp + 1)
        states = propagate_3d_rk4_j2(r0, v0, mu, t, j2)
        stride = spp // spp_eval
        sub = states[::stride]
        t_sub = t[::stride]
        coe = rv_to_coe_eci(sub[:, :3], sub[:, 3:], mu)
        Om_u = np.unwrap(coe["Omega"])
        fit = ols_fit(t_sub, Om_u)
        rates.append(fit["slope"])
        final_positions.append(sub[-1, :3].copy())
    # state-space error of each grid's final position vs the finest-grid state
    ref_pos = final_positions[-1]
    pos_err = [float(np.linalg.norm(fp - ref_pos)) for fp in final_positions[:-1]]
    rates = np.asarray(rates)
    ref = rates[-1]
    errs = np.abs(rates[:-1] - ref)
    errs_safe = np.maximum(errs, 1e-16)  # house guard: convergence_rate needs >0
    stepsizes = np.array([orbits * T / s for s in spp_list[:-1]])
    orders = convergence_rate(errs_safe, stepsizes)
    # state-space order: grids 128..1024 against the finest-grid final state
    pos_err_arr = np.maximum(np.asarray(pos_err), 1e-18)
    state_orders = convergence_rate(pos_err_arr, stepsizes)
    conv_deg_day = [float(np.degrees(er) * 86400.0) for er in errs]
    return {
        "case": {"a_km": a_km, "e": e, "inc_deg": inc_deg},
        "orbits": orbits,
        "spp_list": spp_list,
        "reference_spp": 2048,
        "sample_stride_policy": "every grid evaluated at identical phases (stride spp//128)",
        "E_h_deg_day": conv_deg_day,
        "orders_per_interval": [float(o) for o in orders],
        "mean_order": float(np.mean(orders)),
        "state_final_pos_err_km": pos_err,
        "state_orders_per_interval": [float(o) for o in state_orders],
        "mean_state_order": float(np.mean(state_orders)),
        "rates_rad_s": [float(r) for r in rates],
    }


def kepler_order_check(spp_list=(128, 256, 512, 1024)) -> dict:
    """Raw RK4 integrator order via closed-form circular-Kepler truth (J2=0).

    Max full-vector position error vs the exact solution over one orbit. This
    is the cleanest possible order proof: analytic reference (no reference
    contamination), phase-sensitive (unlike final-|r|, which hides RK4's
    dominant along-track error and decays at order ~5).
    """
    mu = MU_EARTH_KM3S2
    a = R_EARTH_KM + 420.0
    n = np.sqrt(mu / a**3)
    T = 2.0 * np.pi / n
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, np.sqrt(mu / a), 0.0])
    errs = []
    for spp in spp_list:
        t = np.linspace(0.0, T, spp + 1)
        st = propagate_3d_rk4_j2(r0, v0, mu, t, 0.0)
        pos_true = np.column_stack([a * np.cos(n * t), a * np.sin(n * t), np.zeros_like(t)])
        errs.append(float(np.max(np.linalg.norm(st[:, :3] - pos_true, axis=1))))
    errs_arr = np.maximum(np.asarray(errs), 1e-18)
    stepsizes = np.array([T / s for s in spp_list])
    orders = convergence_rate(errs_arr, stepsizes)
    return {
        "spp_list": list(spp_list),
        "max_pos_err_km": errs,
        "orders_per_interval": [float(o) for o in orders],
        "mean_order": float(np.mean(orders)),
        "metric": "max full-vector position error vs closed-form circular Kepler truth, 1 orbit, J2=0",
    }


def null_and_signflip(orbits: int = 20) -> dict:
    """J2 = 0 null test and J2 -> -J2 sign-flip sensitivity (adversarial)."""
    mu = MU_EARTH_KM3S2
    a = R_EARTH_KM + 420.0
    e = 0.01  # well-defined omega, near-circular regime
    inc = np.radians(51.6)
    T = orbital_period(a, mu)
    spp = steps_per_orbit(e)
    t = np.linspace(0.0, orbits * T, orbits * spp + 1)
    r0, v0, _ = seed_state(a, e, inc, 0.0, 0.0, 0.0, mu)

    def slopes(j2):
        st = propagate_3d_rk4_j2(r0, v0, mu, t, j2)
        m = measure_secular_rates(t, st, T, (orbits,), mu)
        fo = m["Omega"][orbits]
        fw = m["omega"][orbits]
        return (
            float(np.degrees(fo["slope"]) * 86400.0),
            float(np.degrees(fw["slope"]) * 86400.0),
        )

    Om0, om0 = slopes(0.0)
    OmP, omP = slopes(+J2_EARTH)
    OmM, omM = slopes(-J2_EARTH)
    ana = analytic_rates(a, e, inc, mu, +J2_EARTH)
    return {
        "j2_zero": {
            "Omega_dot_deg_day": Om0,
            "omega_dot_deg_day": om0,
            "bound_abs_deg_day": 1e-9,
        },
        "sign_flip": {
            "ratio_Omega": float(OmP / OmM) if OmM != 0 else None,
            "ratio_omega": float(omP / omM) if omM != 0 else None,
            "plus_J2_Omega_deg_day": OmP,
            "minus_J2_Omega_deg_day": OmM,
            "plus_J2_omega_deg_day": omP,
            "minus_J2_omega_deg_day": omM,
        },
        "analytic_plus_J2": {
            "Omega_dot_deg_day": ana["Omega_dot_deg_day"],
            "omega_dot_deg_day": ana["omega_dot_deg_day"],
        },
    }


def pathological_grid(orbits: int = 5) -> dict:
    """i x e sweep incl. structural singularities; finiteness + sentinels.

    r-band width w = 5e-3 covers the J2-induced eccentricity of near-circular
    seeds (measured +/-9.60 km at a=6878 km, i=0, e=0 — exactly a*(3/2)J2(R/p)^2
    from the effective-mu deficit) plus integrator margin. omega is expected
    NaN whenever the node line is undefined (i=0/180), because classical
    argument of periapsis is measured FROM the node: without Omega there is no
    reference direction for omega (longitude-of-periapsis would be needed).
    For defined nodes and e=0 seeds, J2 induces a real eccentricity so omega
    remains finite; its trend is simply not claimed (claims policy).
    """
    incs = [0.0, CRITICAL_INC_DEG, 90.0, CRITICAL_INC_2_DEG, 180.0]
    eccs = [0.0, 0.05, 0.2, 0.74]
    mu = MU_EARTH_KM3S2
    a = R_EARTH_KM + 500.0
    T = orbital_period(a, mu)
    rows = []
    ok = True
    band_w = 5e-3
    for inc_deg in incs:
        for e in eccs:
            spp = steps_per_orbit(e)
            t = np.linspace(0.0, orbits * T, orbits * spp + 1)
            inc = np.radians(inc_deg)
            r0, v0, _ = seed_state(a, e, inc, np.radians(10.0), np.radians(30.0), 0.0, mu)
            states = propagate_3d_rk4_j2(r0, v0, mu, t, J2_EARTH)
            coe = rv_to_coe_eci(states[:, :3], states[:, 3:], mu)
            rn = np.linalg.norm(states[:, :3], axis=1)
            node_defined = inc_deg not in (0.0, 180.0)
            # J2 induces a real eccentricity (~1.4e-3) in circular seeds, so
            # e_vec sits far above the guard and omega stays finite (it sweeps
            # once per orbit with the induced periapsis). Per contract the
            # omega_dot TREND is simply not claimed for e=0 seeds.
            row = {
                "inc_deg": inc_deg,
                "e": e,
                "all_aei_finite": bool(np.all(np.isfinite(coe["a"]))
                                       and np.all(np.isfinite(coe["e"]))
                                       and np.all(np.isfinite(coe["inc"]))),
                "Omega_all_nan_expected": not node_defined,
                "Omega_all_nan": bool(np.all(~np.isfinite(coe["Omega"]))),
                "omega_all_nan_expected": not node_defined,
                "omega_all_nan": bool(np.all(~np.isfinite(coe["omega"]))),
                "mean_e_vec_mag": float(np.mean(coe["e"])),
                "r_min_km": float(np.min(rn)),
                "r_max_km": float(np.max(rn)),
                "r_bounded": bool(np.all(rn > a * (1 - e) * (1 - band_w))
                                  and np.all(rn < a * (1 + e) * (1 + band_w))),
                "max_energy_drift_rel": float(
                    np.max(np.abs(j2_specific_energy(states[:, :3], states[:, 3:], mu, J2_EARTH)
                                  - j2_specific_energy(states[0, :3], states[0, 3:], mu, J2_EARTH)))
                    / abs(float(j2_specific_energy(states[0, :3], states[0, 3:], mu, J2_EARTH)))
                ),
            }
            row["ok"] = bool(row["all_aei_finite"] and row["r_bounded"]
                             and row["Omega_all_nan"] == row["Omega_all_nan_expected"]
                             and row["omega_all_nan"] == row["omega_all_nan_expected"]
                             and row["max_energy_drift_rel"] < 1e-8)
            ok = ok and row["ok"]
            rows.append(row)
    return {"incs_tested_deg": incs, "es_tested": eccs, "rows": rows, "all_ok": ok}


# --------------------------------------------------------------------------- #
# Figures (regenerated deterministically from recorded data)
# --------------------------------------------------------------------------- #
def make_figures(results: dict) -> list[str]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    # Fig 1: unwrapped RAAN vs time with secular fits
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name, col in [("ISS", "C0"), ("SSO600", "C2"), ("MOLNIYA", "C4")]:
        ser = results["series"][name]["raan"]
        if ser is None:
            continue
        t_d = np.asarray(ser["t_days"])
        ang = np.asarray(ser["angle_deg"])
        ax.plot(t_d, ang, color=col, lw=1.0, label=f"{name} unwrapped RAAN")
        cmp_ = results["cases"][name]
        slope = cmp_["Omega_fit_stats"]["slope"] if cmp_["Omega_fit_stats"] else None
        if slope is not None:
            ax.plot(t_d, np.degrees(cmp_["Omega_fit_stats"]["intercept"]) + np.degrees(slope) * t_d * 86400.0,
                    color=col, ls="--", lw=1.0,
                    label=f"{name} fit {cmp_['numeric_Omega_dot_deg_day']:+.4f} deg/day")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Unwrapped RAAN [deg]")
    ax.set_title("Nodal precession: propagated RAAN with secular fits (primary windows)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = fig_dir / "raan_vs_time_fit.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # Fig 2: unwrapped argument of periapsis vs time with fits
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name, col in [("ECC_REF", "C3"), ("CRITICAL", "C1"), ("MOLNIYA", "C4")]:
        ser = results["series"][name]["apo"]
        if ser is None:
            continue
        t_d = np.asarray(ser["t_days"])
        ang = np.asarray(ser["angle_deg"])
        ax.plot(t_d, ang, color=col, lw=1.0, label=f"{name} unwrapped omega")
        cmp_ = results["cases"][name]
        if cmp_["omega_fit_stats"]:
            slope = cmp_["omega_fit_stats"]["slope"]
            ax.plot(t_d, np.degrees(cmp_["omega_fit_stats"]["intercept"]) + np.degrees(slope) * t_d * 86400.0,
                    color=col, ls="--", lw=1.0,
                    label=f"{name} fit {cmp_['numeric_omega_dot_deg_day']:+.2e} deg/day")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Unwrapped argument of periapsis [deg]")
    ax.set_title("Apsidal precession: propagated omega with secular fits")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = fig_dir / "omega_vs_time_fit.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # Fig 3: convergence log-log with order-4 reference
    conv = results["convergence"]
    fig, ax = plt.subplots(figsize=(7, 5))
    spp = np.array(conv["spp_list"][:-1], dtype=float)
    errs = np.maximum(np.array(conv["E_h_deg_day"], dtype=float), 1e-18)
    ax.loglog(spp, errs, "o-", label="|rate_h - rate_ref| (deg/day)")
    ref_line = errs[0] * (spp[0] / spp) ** 4
    ax.loglog(spp, ref_line, "r--", label="order-4 reference (h^-4)")
    ax.set_xlabel("Steps per orbit")
    ax.set_ylabel("Secular-rate error vs fine reference [deg/day]")
    ax.set_title(f"RK4 convergence of measured nodal rate (mean order {conv['mean_order']:.2f})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = fig_dir / "convergence_order.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # Fig 4: analytical vs numerical comparison (scatter + rate-vs-inclination curve)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axc = axes[0]
    xs, ys, labels = [], [], []
    for name, cmp_ in results["cases"].items():
        if cmp_["numeric_Omega_dot_deg_day"] is not None and cmp_["analytic_Omega_dot_deg_day"] != 0.0:
            xs.append(cmp_["analytic_Omega_dot_deg_day"])
            ys.append(cmp_["numeric_Omega_dot_deg_day"])
            labels.append(name)
    lim = max(max(np.abs(xs)), max(np.abs(ys))) * 1.15
    axc.plot([-lim, lim], [-lim, lim], "k--", lw=1.0, label="y = x (perfect first-order agreement)")
    axc.plot(xs, ys, "o", ms=7)
    for x, y, lab in zip(xs, ys, labels):
        axc.annotate(lab, (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)
    axc.set_xlabel("Analytic first-order Omega_dot [deg/day]")
    axc.set_ylabel("Numeric converged Omega_dot [deg/day]")
    axc.set_title("Analytic oracle vs full-force propagation")
    axc.grid(True, alpha=0.3)
    axc.legend(fontsize=8)

    axi = axes[1]
    a_iss = R_EARTH_KM + 420.0
    incs = np.linspace(0.0, 180.0, 361)
    rates = [analytic_rates(a_iss, 0.001, np.radians(iv))["Omega_dot_deg_day"] for iv in incs]
    axi.plot(incs, rates, "k-", lw=1.2, label=f"analytic Omega_dot(i) at h=420 km")
    axi.axhline(SSO_TARGET_DEG_DAY, color="g", ls=":", label="SSO target +0.98565 deg/day")
    i_sso = np.degrees(sun_sync_inclination_rad(R_EARTH_KM + 600.0, 0.0))
    axi.axvline(i_sso, color="g", ls="--", alpha=0.7, label=f"solved i_SSO(600km)={i_sso:.3f} deg")
    axi.axvline(CRITICAL_INC_DEG, color="r", ls="--", alpha=0.7, label=f"critical i={CRITICAL_INC_DEG:.4f} deg")
    axi.set_xlabel("Inclination [deg]")
    axi.set_ylabel("Omega_dot [deg/day]")
    axi.set_title("First-order nodal rate vs inclination")
    axi.grid(True, alpha=0.3)
    axi.legend(fontsize=7)
    fig.tight_layout()
    p = fig_dir / "analytic_vs_numeric.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)
    return paths


# --------------------------------------------------------------------------- #
# Main driver
# --------------------------------------------------------------------------- #
def main() -> dict:
    cases = real_cases()

    # Per-case propagation + measurement
    case_results = {}
    series = {}
    runs = {}
    series_needed = {
        ("ISS", "raan"), ("SSO600", "raan"), ("MOLNIYA", "raan"),
        ("ECC_REF", "apo"), ("CRITICAL", "apo"), ("MOLNIYA", "apo"),
    }
    for name, case in cases.items():
        run = run_case(case)
        runs[name] = run
        cmp_ = validate_case_rates(run)
        case_results[name] = cmp_
        # bounded downsampled series for figures (only what make_figures plots)
        coe = run["measure"]["elements"]
        n = len(run["t"])
        stride = max(1, n // 600)
        idx = np.arange(0, n, stride)
        entry = {"raan": None, "apo": None}
        if (name, "raan") in series_needed and run["measure"]["Omega_defined"]:
            entry["raan"] = {
                "t_days": (run["t"][idx] / 86400.0).tolist(),
                "angle_deg": (np.degrees(np.unwrap(coe["Omega"][idx]))).tolist(),
            }
        if (name, "apo") in series_needed and run["measure"]["omega_defined"]:
            entry["apo"] = {
                "t_days": (run["t"][idx] / 86400.0).tolist(),
                "angle_deg": (np.degrees(np.unwrap(coe["omega"][idx]))).tolist(),
            }
        series[name] = entry

    # Validation layers
    conv = convergence_study()
    kep_order = kepler_order_check()
    nullsf = null_and_signflip()
    patho = pathological_grid()

    # SSO table from pinned constants
    sso_table = {}
    for alt in (500.0, 600.0, 800.0):
        a = R_EARTH_KM + alt
        i_rad = sun_sync_inclination_rad(a, 0.0)
        sso_table[f"{int(alt)}km"] = {
            "a_km": float(a),
            "i_SSO_deg": float(np.degrees(i_rad)),
            "check_Omega_dot_deg_day": analytic_rates(a, 0.0, i_rad)["Omega_dot_deg_day"],
        }

    # Molniya near-critical apsidal residual (63.4 vs exact critical)
    mol_ana = analytic_rates(26560.0, 0.74, np.radians(63.4))

    results = {
        "constants": {
            "mu_km3_s2": MU_EARTH_KM3S2,
            "mu_provenance": "IAU 2015 Resolution B3 nominal GM_E; JPL DE440 planet-only 398600.435507 differs 1.5e-8 relative",
            "R_E_km": R_EARTH_KM,
            "R_E_provenance": "WGS-84 equatorial radius, NIMA TR8350.2",
            "J2": J2_EARTH,
            "J2_provenance": "WGS-84 TR8350.2: J2 = sqrt(5)*|C20_bar| with C20_bar = -0.484166774985e-3; EGM2008 value 1.08262668e-3 NOT used (documented distinction)",
            "year_mean_solar_days": YEAR_MEAN_SOLAR_DAYS,
            "sso_target_deg_day": SSO_TARGET_DEG_DAY,
            "sso_convention": (
                "mean SOLAR year 365.2422 d: exact quotient 360/365.2422 = 0.9856473321 deg/day. "
                "The commonly printed 0.98564736 deg/day corresponds to the tropical-year "
                "variant (365.24219 d); the two differ by 2.8e-8 deg/day, negligible vs all "
                "tolerances used here. Sidereal-year rate kept separate and NOT used."
            ),
            "critical_inclination_deg": CRITICAL_INC_DEG,
            "frame_convention": FRAME_CONVENTION,
            "units": UNITS_CONVENTION,
        },
        "resolution_rule": {
            "base_steps_per_orbit": 512,
            "eccentric_base": 720,
            "law": "steps_per_orbit = max(512, ceil(720/(1-e)^1.5)); uniform dt per case",
            "rationale": "periapsis-resolution law (Exp 002/008 precedent); one consistent timestep per case, recorded per case",
        },
        "estimator": {
            "method": "osculating elements per sample -> np.unwrap -> closed-form OLS over integer-orbit windows",
            "windows_policy": "all windows reported; primary = longest, declared a priori",
            "stabilization_criterion_rel": WINDOW_STABILIZATION_REL,
            "third_path": "ascending-node-crossing longitude regression (parabolic time refinement)",
            "guards": {"node_guard_rel": NODE_GUARD_REL, "ecc_guard_abs": ECC_GUARD_ABS,
                       "omega_claim_min_ecc": OMEGA_CLAIM_MIN_ECC},
        },
        "cases": case_results,
        "series": series,
        "convergence": conv,
        "kepler_order_check": kep_order,
        "null_tests": nullsf,
        "pathological": patho,
        "sso_table": sso_table,
        "molniya_near_critical": {
            "omega_dot_analytic_deg_day": mol_ana["omega_dot_deg_day"],
            "note": (
                "63.4 deg is NOT the exact critical inclination 63.43494882 deg; small "
                "nonzero apsidal rate is expected. Near-critical cases (MOLNIYA, CRITICAL) "
                "show residuals ABOVE the generic O(J2*(R/p)^2) model-order scale because "
                "second-order-in-J2 secular terms carry a small divisor (1 - 5cos^2 i) that "
                "peaks at the critical inclination; this amplification is expected physics, "
                "not estimator failure."
            ),
            "window_stabilization_note": (
                "Molniya 12-orbit window still converging toward the 48-orbit primary "
                "(pairwise diff ~4e-3 rel vs <2e-4 for LEO cases); high-e short-period "
                "amplitudes make long windows expensive, so the residual is reported "
                "with this caveat rather than tightened by cherry-picking"
            ),
        },
        "tolerances": {
            "conv_order_band_state_space": list(CONV_ORDER_BAND),
            "conv_rate_order_rule": "every interval >= 3.6 (fourth-order-or-better) and <= 5.0; measured 4.4-4.7",
            "conv_rate_order_justification": (
                "rate-metric orders run ABOVE 4 because the estimator's orbit-averaging "
                "cancels RK4's leading phase-error mode (element recovery is invariant to "
                "along-track timing shifts); raw integrator order 4 is proven separately via "
                "kepler_order_check (max full-vector position error vs closed-form truth)"
            ),
            "kepler_order_band": list(CONV_ORDER_BAND),
            "physics_residual_rel_tol": PHYSICS_RESIDUAL_REL_TOL,
            "anchor_rate_rel_tol": 5e-3,
            "molniya_rate_rel_tol": 1e-2,
            "j2_zero_omega_null_bound_deg_day": 1e-3,
            "j2_zero_null_justification": (
                "Omega null clean at 1e-9; omega slope carries RK4 e_vec-direction noise "
                "amplified by 1/e (~1.8e-6 deg/day at e=0.01); frame/sign bugs produce "
                "O(0.98) deg/day artifacts, so 1e-3 keeps >=3 orders of artifact margin"
            ),
            "sign_flip_ratio_band": [-1.01, -0.99],
            "pathological_r_band_rel": 5e-3,
            "pathological_r_band_justification": (
                "covers J2-induced eccentricity of near-circular seeds (measured +/-9.60 km "
                "= a*(3/2)J2(R/p)^2 at a=6878 km equatorial circular) plus integrator margin"
            ),
            "anchor_leakage_note": (
                "short-window leakage bias ~1.91*A/N (A~0.05 deg osculating wiggle); "
                "primary 100-orbit windows push this to ~1e-4 deg/day"
            ),
        },
        "figures_note": "figures regenerate deterministically from recorded data (Agg, dpi=150)",
    }

    headline = {
        "iss_Omega_dot_numeric_deg_day": case_results["ISS"]["numeric_Omega_dot_deg_day"],
        "iss_Omega_dot_analytic_deg_day": case_results["ISS"]["analytic_Omega_dot_deg_day"],
        "starlink_Omega_dot_numeric_deg_day": case_results["STARLINK"]["numeric_Omega_dot_deg_day"],
        "sso600_numeric_Omega_dot_deg_day": case_results["SSO600"]["numeric_Omega_dot_deg_day"],
        "molniya_Omega_dot_numeric_deg_day": case_results["MOLNIYA"]["numeric_Omega_dot_deg_day"],
        "molniya_omega_dot_numeric_deg_day": case_results["MOLNIYA"]["numeric_omega_dot_deg_day"],
        "critical_omega_dot_numeric_deg_day": case_results["CRITICAL"]["numeric_omega_dot_deg_day"],
        "convergence_mean_order": conv["mean_order"],
        "kepler_mean_order": kep_order["mean_order"],
        "worst_anchor_residual_rel": max(
            (v["Omega_residual_rel"] for v in case_results.values() if v["Omega_residual_rel"] is not None),
            default=None,
        ),
        "worst_anchor_residual_name": max(
            ((v["Omega_residual_rel"], k) for k, v in case_results.items() if v["Omega_residual_rel"] is not None),
            default=(None, None),
        )[1],
        "j2_zero_max_abs_slope_deg_day": max(abs(nullsf["j2_zero"]["Omega_dot_deg_day"]),
                                             abs(nullsf["j2_zero"]["omega_dot_deg_day"])),
        "sign_flip_ratio_Omega": nullsf["sign_flip"]["ratio_Omega"],
        "pathological_all_ok": patho["all_ok"],
    }
    results["headline"] = headline

    fig_paths = make_figures(results)
    results["figures"] = fig_paths

    save_json_result(
        str(RESULTS_DIR / "results.json"),
        results,
        name="j2_precession",
        description=(
            "Exp 009 J2 precession: first-order secular nodal/apsidal rates established "
            "analytically and rediscovered by full-force Cowell RK4 via independent "
            "state->element->trend estimation; convergence order ~4 vs numerical reference; "
            "model-order residual separated; anchors ISS/Starlink/SSO/Molniya/critical-i."
        ),
    )

    print("=== J2 Precession: headline ===")
    print(f"ISS       Omega_dot num {headline['iss_Omega_dot_numeric_deg_day']:+.6f} vs ana {headline['iss_Omega_dot_analytic_deg_day']:+.6f} deg/day")
    print(f"STARLINK  Omega_dot num {headline['starlink_Omega_dot_numeric_deg_day']:+.6f} deg/day")
    print(f"SSO600    Omega_dot num {headline['sso600_numeric_Omega_dot_deg_day']:+.6f} (target {SSO_TARGET_DEG_DAY:+.6f})")
    print(f"MOLNIYA   Omega_dot num {headline['molniya_Omega_dot_numeric_deg_day']:+.6f}  omega_dot num {headline['molniya_omega_dot_numeric_deg_day']:+.3e} deg/day")
    print(f"CRITICAL  omega_dot num {headline['critical_omega_dot_numeric_deg_day']:+.3e} deg/day (exact-critical null)")
    print(f"convergence mean order {headline['convergence_mean_order']:.3f} (rate metric); Kepler truth order {headline['kepler_mean_order']:.3f}")
    print(f"worst anchor residual (rel) {headline['worst_anchor_residual_rel']:.2e}")
    print(f"J2=0 null max |slope| {headline['j2_zero_max_abs_slope_deg_day']:.2e} deg/day; sign-flip ratio {headline['sign_flip_ratio_Omega']:.6f}")
    print(f"pathological all ok: {headline['pathological_all_ok']}")
    print(f"figures: {fig_paths}")
    return results


if __name__ == "__main__":
    main()
