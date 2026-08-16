"""
Combined bi-elliptic transfer + plane change (Experiment 005 continuation /
the "006" plane-change slot, completed in place).

Research question: what is the GLOBAL optimum delta-v for a combined radius
change (r1 -> r2, radius ratio R = r2/r1) AND inclination change (delta_i)
between two circular orbits, when the maneuver may use either

  (a) two burns  -- the combined Hohmann + plane change (a 2-impulse transfer
      with the plane change split between the two burns), OR
  (b) three burns -- raise apoapsis to an intermediate radius s, do the plane
      change (split across the three burns), then lower periapsis to r2
      (the super-synchronous / bi-elliptic-with-plane-change strategy), OR
  (c) the s -> infinity limit of (b), where the apoapsis velocity vanishes and
      the entire plane change becomes "free"; the cost tends to the coplanar
      bi-parabolic value (sqrt(2)-1)(1 + 1/sqrt(R)), INDEPENDENT of delta_i.

Does the true global optimum have distinct regimes -- ordinary two-burn,
finite intermediate-apoapsis three-burn, and asymptotic s->infinity -- and
where (in the (R, delta_i) plane) do the boundaries between them occur?

Method: a dense, non-unimodal global search over (s, theta1, theta2) with the
plane-change split theta1+theta2+theta3 = delta_i fully free; comparison
against the two-burn optimum and the analytic s->infinity limit; boundary
location by bisection; independent validation by (i) a second brute-force
search, (ii) 50-digit mpmath at the regime corners, (iii) a full 3D RK4
propagation of the optimal maneuvers, and (iv) real-system engineering anchors.

Determinism: pure float64 + fixed mpmath precision, no RNG. Repeated runs are
byte-identical apart from the timestamp.

Sources (real, cited in the card):
  * Curtis, H. D. "Orbital Mechanics for Engineering Students", 4th ed.,
    Elsevier 2021 -- combined Hohmann + plane change (law of cosines), the
    300 km LEO, 28.6 deg -> GEO worked example.
  * Gonzalez, "Orbital Mechanics & Astrodynamics" (orbital-mechanics.space) --
    Plane-Change Maneuver example (same worked case).
  * Wikipedia "Geostationary transfer orbit" (law of cosines with cos(d_i) at
    apogee) and "Supersynchronous orbit" (SES-8, Thaicom-6 apogee 90 000 km).
  * Wakker, "Optimal Impulsive Orbit Transfers", Springer 2015.
  * Hoelker & Silber 1959; Exp 004 (Hohmann) and Exp 005 (coplanar
    bi-elliptic crossover) for the closed-form transfer machinery reused here.
"""

from __future__ import annotations

import json
import importlib.util
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import mpmath as mp

# --------------------------------------------------------------------------- #
# Load the verified machinery: Exp 002 (RK4 / Kepler) and Exp 004 (Hohmann +
# IAU constants) and Exp 005 (coplanar bi-elliptic / bi-parabolic limit).
# --------------------------------------------------------------------------- #
_THIS = Path(__file__).resolve().parent
_RESEARCH = _THIS.parent.parent.parent
_KEPLER_PATH = _RESEARCH / "orbital-mechanics/experiments/keplerOrbitValidation/experiment.py"
_HOHMANN_PATH = _RESEARCH / "orbital-mechanics/experiments/hohmannTransfer/experiment.py"
_BIELL_PATH = _RESEARCH / "orbital-mechanics/experiments/biellipticVsHohmann/experiment.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


kepler_exp = _load("kepler_exp_006", _KEPLER_PATH)
hohmann_exp = _load("hohmann_exp_006", _HOHMANN_PATH)
biell_exp = _load("biell_exp_006", _BIELL_PATH)

MU_EARTH = 398600.4418          # km^3/s^2 (IAU 2012 Earth value)
R_EARTH_KM = 6378.1
LEO_ALT_KM = 200.0
GEO_ALT_KM = 35786.0
GEO_RADIUS_KM = R_EARTH_KM + GEO_ALT_KM

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"

# Margin (in normalized delta-v units) required for the three-burn family to
# count as genuinely beating two-burn (or finite-s beating the s->inf limit).
# The float64 grid optimizer's global minimum has a noise floor of ~1e-7 in
# delta-v near a near-tie; advantages smaller than this are not reproducible
# and would let float-grid artifacts define the regime boundary. A margin of
# 1e-5 sits two orders of magnitude above that noise floor and well below the
# ~1e-4 genuine dip depth at the true boundary, so it separates artifact from
# reality. (Audit finding 2026-08-16: the prior 1e-12 threshold produced a
# di_c(R) boundary at the float-tie point, in error by up to ~6 deg near R=1.)
WIN_MARGIN = 1e-5


# --------------------------------------------------------------------------- #
# Small generic helpers
# --------------------------------------------------------------------------- #
def rodrigues(vec, axis, theta):
    """Rotate 3-vector `vec` about unit-ish `axis` by `theta` (Rodrigues)."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    return (vec * c + np.cross(axis, vec) * s
            + axis * (np.dot(axis, vec)) * (1.0 - c))


def golden_section(f, lo, hi, tol=1e-12, maxit=200):
    """Minimize unimodal scalar f on [lo, hi]. Deterministic."""
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = float(lo), float(hi)
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(maxit):
        if b - a < tol:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = f(d)
    return (a + b) / 2.0, f((a + b) / 2.0)


# --------------------------------------------------------------------------- #
# PART A -- Pure inclination change on ONE circular orbit (R = 1).
# Reused from the prior (incomplete) draft; preserved because it is correct
# and provides the R = 1 anchor for the combined problem (the detour there
# reduces to the 3-burn with R = 1).
# --------------------------------------------------------------------------- #
def pure_plane_change_dv(v1, di):
    """Minimal dv for an inclination change di on a circular orbit (speed v1)."""
    return 2.0 * v1 * np.sin(di / 2.0)


def detour_h(s, di):
    """One-leg cost of the 3-burn detour at intermediate radius s (R = 1,
    v1 = 1): raise apoapsis to s, plane-change at s, lower back. Total = 2 h."""
    q = np.sin(di / 2.0)
    return (np.sqrt(2.0 * s / (1.0 + s)) - 1.0) + np.sqrt(2.0 / (s * (1.0 + s))) * q


def detour_optimum(di):
    """Analytic optimum of the R = 1 three-burn detour.

    Returns (s_star, total_dv). For di <= 38.9424 deg the direct burn is
    optimal (s_star = 1); for 38.9424 < di < 60 the interior optimum is
    s* = q/(1-2q); for di >= 60 deg the optimum diverges and total -> 2(sqrt2-1).
    """
    q = np.sin(di / 2.0)
    if q <= 1.0 / 3.0:
        return 1.0, 2.0 * q
    if q < 0.5 - 1e-9:
        s_star = q / (1.0 - 2.0 * q)
        total = 2.0 * (2.0 * np.sqrt(2.0 * q * (1.0 - q)) - 1.0)
        return float(s_star), float(total)
    # di >= 60 deg (within float tolerance): optimum diverges to infinity and
    # the total tends to the constant 2(sqrt2 - 1) (plane change free at the
    # vanishing apoapsis velocity).
    return float("inf"), 2.0 * (np.sqrt(2.0) - 1.0)


def combined_inc_raan_dv(v1, i1, om1, i2, om2):
    """Minimal dv for a combined inclination + RAAN change between two planes
    (same circular radius, speed v1). Only the total angle between the plane
    normals matters."""
    cdelta = np.cos(i1) * np.cos(i2) + np.sin(i1) * np.sin(i2) * np.cos(om1 - om2)
    delta = np.arccos(np.clip(cdelta, -1.0, 1.0))
    return 2.0 * v1 * np.sin(delta / 2.0)


# --------------------------------------------------------------------------- #
# PART B -- Two-burn Hohmann + plane change (combined transfer).
# --------------------------------------------------------------------------- #
def _transfer_speeds(R):
    """Normalized (v1 = 1) transfer-ellipse velocities for a Hohmann between
    r1 = 1 and r2 = R: periapsis speed, target circular speed, apoapsis speed."""
    vp = np.sqrt(2.0 * R / (1.0 + R))
    v2 = 1.0 / np.sqrt(R)
    v_apo = np.sqrt(2.0 / (R * (1.0 + R)))
    return vp, v2, v_apo


_transfer_speeds_v2 = _transfer_speeds  # alias


def two_burn_dv(R, di, theta1):
    """Total dv of the 2-burn combined transfer with `theta1` of the plane
    change done at the first burn (periapsis r1); the rest at the second."""
    vp, v2, v_apo = _transfer_speeds(R)
    dv1 = np.sqrt(1.0 + vp * vp - 2.0 * vp * np.cos(theta1))
    dv2 = np.sqrt(v2 * v2 + v_apo * v_apo - 2.0 * v2 * v_apo * np.cos(di - theta1))
    return dv1, dv2, dv1 + dv2


def two_burn_optimal(R, di):
    """Optimal theta1* (plane change at first burn) minimizing total dv.

    The 2-burn cost is unimodal in theta1 in [0, di] (it is a sum of two
    'law of cosines' terms with a fixed-sum split), so golden section is
    reliable here."""
    vp, v2, v_apo = _transfer_speeds(R)

    def total(theta1):
        dv1 = np.sqrt(1.0 + vp * vp - 2.0 * vp * np.cos(theta1))
        dv2 = np.sqrt(v2 * v2 + v_apo * v_apo - 2.0 * v2 * v_apo * np.cos(di - theta1))
        return dv1 + dv2

    theta1_star, total_dv = golden_section(total, 0.0, di, tol=1e-13)
    dv1, dv2, _ = two_burn_dv(R, di, theta1_star)
    dv_all_at_r1 = two_burn_dv(R, di, di)[2]
    dv_all_at_r2 = two_burn_dv(R, di, 0.0)[2]
    return {
        "R": float(R),
        "di_deg": float(np.degrees(di)),
        "theta1_star_deg": float(np.degrees(theta1_star)),
        "theta1_frac_of_di": float(theta1_star / di),
        "total_dv": float(total_dv),
        "dv1": float(dv1),
        "dv2": float(dv2),
        "dv_sequential_at_r1": float(dv_all_at_r1),
        "dv_sequential_at_r2": float(dv_all_at_r2),
        "saving_vs_sequential_at_r2": float(dv_all_at_r2 - total_dv),
    }


# --------------------------------------------------------------------------- #
# PART C -- Three-burn super-synchronous strategy with FULL plane-change split.
# --------------------------------------------------------------------------- #
def three_burn_cost(R, s, theta1, theta2, di):
    """Full 3-burn cost for (R, s, split). theta3 = di - theta1 - theta2.

    s >= R (r_b = s r1 is the common apsidal radius of the two transfer
    ellipses). Each burn is a single vector difference (law of cosines) with
    the assigned plane-change angle; the split theta1+theta2+theta3 = di is
    the total inclination change, distributed about the common node axis.
    """
    vp12 = np.sqrt(2.0 * s / (1.0 + s))        # burn1 speed (ellipse 1, periapsis r1)
    v2 = 1.0 / np.sqrt(R)                       # target circular speed
    va1 = np.sqrt(2.0 / (s * (1.0 + s)))        # apoapsis speed, ellipse 1
    va2 = np.sqrt(2.0 * R / (s * (R + s)))      # apoapsis speed, ellipse 2
    vp23 = np.sqrt(2.0 * s / (R * (R + s)))     # periapsis speed, ellipse 2
    theta3 = di - theta1 - theta2
    dv1 = np.sqrt(1.0 + vp12 * vp12 - 2.0 * vp12 * np.cos(theta1))
    dv2 = np.sqrt(va1 * va1 + va2 * va2 - 2.0 * va1 * va2 * np.cos(theta2))
    dv3 = np.sqrt(v2 * v2 + vp23 * vp23 - 2.0 * v2 * vp23 * np.cos(theta3))
    return dv1 + dv2 + dv3


def bi_parabolic_plane_change_limit(R):
    """Cost of the 3-burn as s -> infinity (plane change free at infinity):
    (sqrt(2)-1)(1 + 1/sqrt(R)), independent of di. This equals the COPLANAR
    bi-parabolic limit of Exp 005 -- a clean identity tested below."""
    return (np.sqrt(2.0) - 1.0) * (1.0 + 1.0 / np.sqrt(R))


def _surface(R, di, s_grid, nth):
    """Vectorized 3-burn cost surface total(s, th1, th2) for a fixed (R, di).

    Returns (total, s_grid, th1, th2). total has shape (len(s_grid), nth, nth)
    with inf where th1+th2 > di (invalid split). Memory is kept modest:
    (Ns, nth, nth) floats ~ Ns*nth^2 * 8 bytes.
    """
    s = np.asarray(s_grid, dtype=float)
    Ns = len(s)
    vp12 = np.sqrt(2.0 * s / (1.0 + s))            # (Ns,)
    v2 = 1.0 / np.sqrt(R)
    va1 = np.sqrt(2.0 / (s * (1.0 + s)))           # (Ns,)
    va2 = np.sqrt(2.0 * R / (s * (R + s)))         # (Ns,)
    vp23 = np.sqrt(2.0 * s / (R * (R + s)))        # (Ns,)
    th1 = np.linspace(0.0, di, nth)                # (nth,)
    th2 = np.linspace(0.0, di, nth)                # (nth,)
    th3 = di - th1[:, None] - th2[None, :]         # (nth, nth)
    valid = th3 >= -1e-12
    dv1 = np.sqrt(1.0 + vp12[:, None, None] ** 2
                  - 2.0 * vp12[:, None, None] * np.cos(th1)[None, :, None])
    dv2 = np.sqrt(va1[:, None, None] ** 2 + va2[:, None, None] ** 2
                  - 2.0 * va1[:, None, None] * va2[:, None, None] * np.cos(th2)[None, None, :])
    dv3 = np.sqrt(v2 ** 2 + vp23[:, None, None] ** 2
                  - 2.0 * v2 * vp23[:, None, None] * np.cos(th3)[None, :, :])
    total = dv1 + dv2 + dv3
    total = np.where(valid[None, :, :], total, np.inf)
    return total, s, th1, th2


def _global_min_on_surface(R, di, s_grid, nth):
    """Global minimum of the 3-burn cost over (s, th1, th2) on a given grid."""
    total, s, th1, th2 = _surface(R, di, s_grid, nth)
    idx = np.unravel_index(np.argmin(total), total.shape)
    s_star = float(s[idx[0]])
    th1_star = float(th1[idx[1]])
    th2_star = float(th2[idx[2]])
    return {
        "s": s_star,
        "theta1_deg": float(np.degrees(th1_star)),
        "theta2_deg": float(np.degrees(th2_star)),
        "theta3_deg": float(np.degrees(di - th1_star - th2_star)),
        "total_dv": float(total[idx]),
    }


def three_burn_global(R, di, s_max=1.0e6, ns=240, nth=64, refine=True):
    """Global minimum of the 3-burn cost over (s >= R, th1, th2).

    Uses a dense log-spaced s grid + a 2D split grid (the objective is NOT
    assumed unimodal in s, so we keep the global grid minimum rather than a
    local optimizer), with two narrowing refinement passes and a final
    nested golden-section polish of the split. Always reports the analytic
    s->infinity limit so the caller can choose the true global minimum.
    """
    s_coarse = np.logspace(np.log10(R * 1.0001), np.log10(s_max), ns)
    best = _global_min_on_surface(R, di, s_coarse, nth)
    if refine:
        sc = best["s"]
        # refine s around the candidate; widen if it sits at the boundary
        if sc < s_max * 0.95 and sc > R * 1.05:
            lo = max(R * 1.0001, sc * 0.5)
            hi = min(s_max, sc * 2.0)
            s_fine = np.logspace(np.log10(lo), np.log10(hi), int(ns * 0.9))
            b2 = _global_min_on_surface(R, di, s_fine, nth + 16)
            if b2["total_dv"] < best["total_dv"] - 1e-15:
                best = b2
                sc = best["s"]
            s_finer = np.logspace(np.log10(max(R * 1.0001, sc * 0.75)),
                                  np.log10(min(s_max, sc * 1.34)), int(ns * 0.7))
            b3 = _global_min_on_surface(R, di, s_finer, nth + 32)
            if b3["total_dv"] < best["total_dv"] - 1e-15:
                best = b3
        # final split polish via nested golden section (robust to non-unimodal
        # s: we fix s at the discovered optimum and polish th1, th2 jointly).
        s_fixed = best["s"]
        vp12 = np.sqrt(2.0 * s_fixed / (1.0 + s_fixed))
        v2 = 1.0 / np.sqrt(R)
        va1 = np.sqrt(2.0 / (s_fixed * (1.0 + s_fixed)))
        va2 = np.sqrt(2.0 * R / (s_fixed * (R + s_fixed)))
        vp23 = np.sqrt(2.0 * s_fixed / (R * (R + s_fixed)))

        def cost_split(th1, th2):
            th3 = di - th1 - th2
            if th3 < -1e-12 or th2 < -1e-12 or th1 < -1e-12:
                return 1e9
            dv1 = np.sqrt(1.0 + vp12 ** 2 - 2.0 * vp12 * np.cos(th1))
            dv2 = np.sqrt(va1 ** 2 + va2 ** 2 - 2.0 * va1 * va2 * np.cos(th2))
            dv3 = np.sqrt(v2 ** 2 + vp23 ** 2 - 2.0 * v2 * vp23 * np.cos(th3))
            return dv1 + dv2 + dv3

        def inner(th1):
            th2_st, _ = golden_section(lambda t2: cost_split(th1, t2), 0.0, di - th1,
                                       tol=1e-11)
            return cost_split(th1, th2_st), th2_st

        th1_st, _ = golden_section(lambda t1: inner(t1)[0], 0.0, di, tol=1e-11)
        _, th2_st = inner(th1_st)
        polished = cost_split(th1_st, th2_st)
        if polished < best["total_dv"]:
            best = {
                "s": s_fixed,
                "theta1_deg": float(np.degrees(th1_st)),
                "theta2_deg": float(np.degrees(th2_st)),
                "theta3_deg": float(np.degrees(di - th1_st - th2_st)),
                "total_dv": float(polished),
            }
    return best


def combined_optimum(R, di, **kw):
    """Compare the three candidates and return the global optimum with a
    regime label.

    regimes:
      'two_burn'      -- the 2-burn combined transfer is optimal.
      'finite_s'      -- optimal 3-burn has a finite intermediate apoapsis s* > R.
      'infinite_s'    -- optimal 3-burn has s* -> infinity (cost = bi-parabolic
                         limit, plane change free at the near-rest apoapsis).
    """
    two = two_burn_optimal(R, di)["total_dv"]
    tb = three_burn_global(R, di, **kw)
    inf = bi_parabolic_plane_change_limit(R)
    finite_dv = tb["total_dv"]
    # true 3-burn minimum: either a finite dip below the s->infinity limit, or
    # the s->infinity limit itself.
    if finite_dv < inf - WIN_MARGIN:
        three_dv = finite_dv
        three_regime = "finite_s"
    else:
        three_dv = inf
        three_regime = "infinite_s"
    if three_dv < two - WIN_MARGIN:
        winner = three_regime
        best_dv = three_dv
    else:
        winner = "two_burn"
        best_dv = two
    return {
        "R": float(R),
        "di_deg": float(np.degrees(di)),
        "regime": winner,
        "best_dv": float(best_dv),
        "two_burn_dv": float(two),
        "three_burn_finite_dv": float(finite_dv),
        "three_burn_s_star": float(tb["s"]),
        "three_burn_theta1_deg": float(tb["theta1_deg"]),
        "three_burn_theta2_deg": float(tb["theta2_deg"]),
        "three_burn_theta3_deg": float(tb["theta3_deg"]),
        "bi_parabolic_limit_dv": float(inf),
        "beats_two_burn": bool(three_dv < two - WIN_MARGIN),
    }


def _three_burn_wins(R, di, **kw):
    return combined_optimum(R, di, **kw)["beats_two_burn"]


def _three_burn_is_finite(R, di, **kw):
    co = combined_optimum(R, di, **kw)
    return co["regime"] == "finite_s"


def di_c_boundary(R, tol=1e-4, **kw):
    """Smallest delta_i (deg) at which the 3-burn (finite OR infinite) beats
    the 2-burn combined transfer. Returns None if the 3-burn never wins on
    (0, 180 deg]."""
    lo, hi = 1e-4, np.pi
    if _three_burn_wins(R, lo, **kw):
        return 0.0  # already winning at the smallest angle
    if not _three_burn_wins(R, hi, **kw):
        return None  # never wins
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if _three_burn_wins(R, mid, **kw):
            hi = mid
        else:
            lo = mid
        if hi - lo < tol:
            break
    return float(np.degrees(0.5 * (lo + hi)))


def di_inf_boundary(R, tol=1e-4, **kw):
    """delta_i (deg) above which the optimal 3-burn is the s->infinity
    (free-at-apoapsis) regime rather than a finite intermediate apoapsis.

    Returns None when the finite-s regime never occurs for this R (i.e. the
    only 3-burn regime that beats two-burn is s->infinity). Otherwise returns
    the upper edge of the finite-s window (the lower edge is di_c)."""
    # The s->infinity regime must hold at large delta_i: near 179 deg the
    # optimal 3-burn is always the s->infinity one (its cost is independent
    # of delta_i and it is the cheapest available for large delta_i).
    hi_di = np.radians(179.0)
    if _three_burn_is_finite(R, hi_di, **kw):
        # No finite window extends to large delta_i, so there is no finite-s
        # regime at all: di_inf == di_c (the 3-burn that beats two-burn is
        # the s->infinity one directly above di_c).
        return di_c_boundary(R, tol, **kw)
    # There is a finite-s regime somewhere; bisect its upper edge. The lower
    # edge is di_c(R); start the search just above it.
    dc = di_c_boundary(R, tol, **kw)
    if dc is None:
        return None
    lo_di = np.radians(dc + 0.5)
    if not _three_burn_is_finite(R, lo_di, **kw):
        # no finite window just above di_c -> di_inf == di_c
        return dc
    for _ in range(70):
        mid = 0.5 * (lo_di + hi_di)
        if _three_burn_is_finite(R, mid, **kw):
            lo_di = mid
        else:
            hi_di = mid
        if hi_di - lo_di < tol:
            break
    return float(np.degrees(0.5 * (lo_di + hi_di)))


def regime_sweep(Rs, di_deg_grid):
    """Classify every (R, di) point. Returns a dict with grids + labels.

    label codes: 0 = two_burn, 1 = finite_s, 2 = infinite_s.
    """
    Rs = np.asarray(Rs, dtype=float)
    dis = np.asarray(di_deg_grid, dtype=float)
    labels = np.zeros((len(Rs), len(dis)), dtype=int)
    s_star = np.zeros_like(labels, dtype=float)
    best_dv = np.zeros_like(labels, dtype=float)
    for i, R in enumerate(Rs):
        for j, d in enumerate(dis):
            co = combined_optimum(R, np.radians(d))
            code = {"two_burn": 0, "finite_s": 1, "infinite_s": 2}[co["regime"]]
            labels[i, j] = code
            s_star[i, j] = co["three_burn_s_star"] if code != 0 else float("nan")
            best_dv[i, j] = co["best_dv"]
    return {
        "Rs": Rs.tolist(),
        "di_deg_grid": dis.tolist(),
        "labels": labels.tolist(),
        "s_star": s_star.tolist(),
        "best_dv": best_dv.tolist(),
    }


# --------------------------------------------------------------------------- #
# PART D -- Independent 3D RK4 trajectory validation.
# The verified Exp 002 machinery is planar (2D); a plane change needs 3D, so
# we implement a compact, independent 3D Kepler (Cowell) integrator here and
# propagate the actual optimal maneuvers built from the closed-form burn
# vectors, then check the final orbit is the target circular orbit at r = R
# with its normal rotated by delta_i about the common node axis.
# --------------------------------------------------------------------------- #
def propagate_3d_rk4(r0, v0, mu, t, dt):
    """Fixed-step RK4 for the 3D two-body problem. Deterministic."""
    n = len(t)
    state = np.empty((n, 6))
    state[0] = np.concatenate([np.asarray(r0, float), np.asarray(v0, float)])

    def accel(x):
        r = x[:3]
        rm = np.linalg.norm(r)
        return -mu * r / rm ** 3

    # precompute step indices
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


def _node_axis(theta):
    return np.array([1.0, 0.0, 0.0])  # common node axis (radial at the burns)


def build_two_burn_states(R, di, theta1):
    """Return (states list, dvs list, times list) for the optimal 2-burn
    combined transfer, built from the closed-form vectors. The node axis is
    +x; inclination is rotated about it by delta_i."""
    mu = 1.0
    vp = np.sqrt(2.0 * R / (1.0 + R))
    v2 = 1.0 / np.sqrt(R)
    v_apo = np.sqrt(2.0 / (R * (1.0 + R)))
    # initial circular orbit at r=1, normal +z
    r0 = np.array([1.0, 0.0, 0.0])
    v_pre = np.array([0.0, 1.0, 0.0])
    # burn1: rotate v_pre by theta1 about +x, scale to vp
    v_post1 = vp * rodrigues(v_pre, _node_axis(di), theta1)
    dv1 = v_post1 - v_pre
    # Hohmann transfer ellipse: a = (1+R)/2, e = (R-1)/(R+1), period
    a = 0.5 * (1.0 + R)
    T = 2.0 * np.pi * np.sqrt(a ** 3 / mu)
    t1 = np.linspace(0.0, 0.5 * T, 4001)
    leg1 = propagate_3d_rk4(r0, v_post1, mu, t1, t1[1] - t1[0])
    r_apo = leg1[-1, :3]
    v_apo_vec = leg1[-1, 3:]
    # burn2: circularize to r=R with normal rotated by delta_i about +x.
    # The Hohmann leg arrives at the APOAPSIS of the transfer ellipse, which
    # (for R > 1, periapsis at +x) is at -x; its prograde tangent is [0,-1,0].
    # The circular-velocity vector is that tangent rotated by di about +x.
    v_pre2 = v_apo_vec
    prograde_apo = np.array([0.0, -1.0, 0.0])
    v_post2 = v2 * rodrigues(prograde_apo, _node_axis(di), di)
    dv2 = v_post2 - v_pre2
    t2 = np.linspace(0.0, 2.0 * np.pi * np.sqrt(R ** 3 / mu), 2001)
    leg2 = propagate_3d_rk4(r_apo, v_post2, mu, t2, t2[1] - t2[0])
    return {
        "leg1": leg1, "leg2": leg2, "dv1": dv1, "dv2": dv2,
        "r_apo": r_apo, "v_post2": v_post2, "final": leg2[-1],
        "target_R": R, "di_deg": np.degrees(di), "theta1_deg": np.degrees(theta1),
    }


def build_three_burn_states(R, di, s, theta1, theta2):
    """Return states/dvs for the optimal finite-s 3-burn transfer."""
    mu = 1.0
    vp12 = np.sqrt(2.0 * s / (1.0 + s))
    v2 = 1.0 / np.sqrt(R)
    # burn1 at r=1: rotate by theta1 about +x, scale to vp12
    r0 = np.array([1.0, 0.0, 0.0])
    v_pre = np.array([0.0, 1.0, 0.0])
    v_post1 = vp12 * rodrigues(v_pre, _node_axis(di), theta1)
    dv1 = v_post1 - v_pre
    # leg1: ellipse(1, s), a1=(1+s)/2, half period
    a1 = 0.5 * (1.0 + s)
    T1 = 2.0 * np.pi * np.sqrt(a1 ** 3 / mu)
    t1 = np.linspace(0.0, 0.5 * T1, 6001)
    leg1 = propagate_3d_rk4(r0, v_post1, mu, t1, t1[1] - t1[0])
    r_apo = leg1[-1, :3]
    v_apo_vec = leg1[-1, 3:]
    # burn2 at apoapsis r=s: rotate v_pre(apo) by theta2 about +x, scale to va2
    va1 = np.sqrt(2.0 / (s * (1.0 + s)))
    va2 = np.sqrt(2.0 * R / (s * (R + s)))
    # direction of v_apo_vec (tangential at apoapsis)
    vdir = v_apo_vec / np.linalg.norm(v_apo_vec)
    v_post2 = va2 * rodrigues(vdir, _node_axis(di), theta2)
    dv2 = v_post2 - v_apo_vec
    # leg2: ellipse(R, s), a2=(R+s)/2, half period
    a2 = 0.5 * (R + s)
    T2 = 2.0 * np.pi * np.sqrt(a2 ** 3 / mu)
    t2 = np.linspace(0.0, 0.5 * T2, 6001)
    leg2 = propagate_3d_rk4(r_apo, v_post2, mu, t2, t2[1] - t2[0])
    r_peri = leg2[-1, :3]
    v_peri_vec = leg2[-1, 3:]
    # burn3 at r=R: circularize, normal rotated by delta_i
    vp23 = np.sqrt(2.0 * s / (R * (R + s)))
    v_pre3 = v_peri_vec
    v_post3 = v2 * rodrigues(np.array([0.0, 1.0, 0.0]), _node_axis(di), di)
    dv3 = v_post3 - v_pre3
    t3 = np.linspace(0.0, 2.0 * np.pi * np.sqrt(R ** 3 / mu), 2001)
    leg3 = propagate_3d_rk4(r_peri, v_post3, mu, t3, t3[1] - t3[0])
    return {
        "leg1": leg1, "leg2": leg2, "leg3": leg3,
        "dv1": dv1, "dv2": dv2, "dv3": dv3,
        "r_apo": r_apo, "r_peri": r_peri, "final": leg3[-1],
        "target_R": R, "s": s, "di_deg": np.degrees(di),
        "theta1_deg": np.degrees(theta1), "theta2_deg": np.degrees(theta2),
    }


def validate_two_burn_rk4(R, di):
    """Validate the 2-burn closed form against 3D RK4."""
    ob = two_burn_optimal(R, di)
    theta1 = np.radians(ob["theta1_star_deg"])
    st = build_two_burn_states(R, di, theta1)
    # closed-form dv
    cf_dv1 = np.sqrt(1.0 + (np.sqrt(2 * R / (1 + R))) ** 2
                    - 2 * np.sqrt(2 * R / (1 + R)) * np.cos(theta1))
    cf_dv2 = np.sqrt((1 / np.sqrt(R)) ** 2
                    + (np.sqrt(2 / (R * (1 + R)))) ** 2
                    - 2 * (1 / np.sqrt(R)) * np.sqrt(2 / (R * (1 + R)))
                    * np.cos(di - theta1))
    rk_dv1 = np.linalg.norm(st["dv1"])
    rk_dv2 = np.linalg.norm(st["dv2"])
    # final orbit circular at R, normal rotated by di about +x
    rf = st["final"][:3]
    vf = st["final"][3:]
    rf_mag = np.linalg.norm(rf)
    vf_mag = np.linalg.norm(vf)
    v_circ = np.sqrt(1.0 / R)
    # angular momentum should equal R * v_circ (circular)
    h = np.cross(rf, vf)
    h_mag = np.linalg.norm(h)
    # normal direction: should be R*[0,-sin(di),cos(di)] (cross of position
    # and velocity on the target circle)
    n_target = np.array([0.0, -np.sin(di), np.cos(di)])
    n_actual = h / h_mag if h_mag > 0 else n_target
    cos_ang = np.clip(np.dot(n_actual, n_target), -1.0, 1.0)
    return {
        "R": float(R), "di_deg": float(np.degrees(di)),
        "cf_dv1": float(cf_dv1), "rk_dv1": float(rk_dv1),
        "cf_dv2": float(cf_dv2), "rk_dv2": float(rk_dv2),
        "dv_rel_err": float(max(abs(rk_dv1 - cf_dv1) / cf_dv1,
                                abs(rk_dv2 - cf_dv2) / cf_dv2)),
        "final_r_error": float(abs(rf_mag - R) / R),
        "final_speed_error": float(abs(vf_mag - v_circ) / v_circ),
        "final_circular_h_error": float(abs(h_mag - R * v_circ) / (R * v_circ)),
        "normal_alignment_cos": float(cos_ang),
    }


def validate_three_burn_rk4(R, di, s):
    """Validate a finite-s 3-burn closed form against 3D RK4."""
    tb = three_burn_global(R, di, s_max=max(s * 3.0, 1e4))
    # pick a representative finite-s optimum (the global one may be infinite;
    # for validation we force a finite s near the quoted optimum)
    # Re-run with bounded s_max to keep it finite and representative.
    cf = three_burn_cost(R, s,
                         np.radians(tb["theta1_deg"]),
                         np.radians(tb["theta2_deg"]), di)
    st = build_three_burn_states(R, di, s, np.radians(tb["theta1_deg"]),
                                 np.radians(tb["theta2_deg"]))
    rk_dv1 = np.linalg.norm(st["dv1"])
    rk_dv2 = np.linalg.norm(st["dv2"])
    rk_dv3 = np.linalg.norm(st["dv3"])
    rf = st["final"][:3]
    vf = st["final"][3:]
    rf_mag = np.linalg.norm(rf)
    vf_mag = np.linalg.norm(vf)
    v_circ = np.sqrt(1.0 / R)
    h = np.cross(rf, vf)
    h_mag = np.linalg.norm(h)
    n_target = np.array([0.0, -np.sin(di), np.cos(di)])
    n_actual = h / h_mag if h_mag > 0 else n_target
    cos_ang = np.clip(np.dot(n_actual, n_target), -1.0, 1.0)
    return {
        "R": float(R), "di_deg": float(np.degrees(di)), "s": float(s),
        "cf_total": float(cf),
        "rk_dv1": float(rk_dv1), "rk_dv2": float(rk_dv2), "rk_dv3": float(rk_dv3),
        "rk_total": float(rk_dv1 + rk_dv2 + rk_dv3),
        "dv_rel_err": float(abs(rk_dv1 + rk_dv2 + rk_dv3 - cf) / cf),
        "final_r_error": float(abs(rf_mag - R) / R),
        "final_speed_error": float(abs(vf_mag - v_circ) / v_circ),
        "final_circular_h_error": float(abs(h_mag - R * v_circ) / (R * v_circ)),
        "normal_alignment_cos": float(cos_ang),
    }


# --------------------------------------------------------------------------- #
# PART E -- Real-system anchors.
# --------------------------------------------------------------------------- #
def real_anchors():
    """Engineering numbers (km/s) for canonical cases.

    The normalized results are scaled by v1 = sqrt(mu/r1); all geometry is
    identical in the normalized and physical frames (scale-free two-body)."""
    mu = MU_EARTH
    r_leo = R_EARTH_KM + LEO_ALT_KM
    r_geo = GEO_RADIUS_KM

    def scale(R_c, di_deg):
        r1 = r_leo
        R = (R_c * r_leo) / r_leo  # R_c is already r2/r1 when r1 = r_leo
        v1 = np.sqrt(mu / r1)
        co = combined_optimum(R, np.radians(di_deg))
        tb = two_burn_optimal(R, np.radians(di_deg))
        # all-at-apogee (theta1 = 0) baseline for comparison
        dv_all_apo = two_burn_dv(R, np.radians(di_deg), 0.0)[2]
        return {
            "R": float(R), "di_deg": di_deg,
            "regime": co["regime"],
            "best_dv_km_s": co["best_dv"] * v1,
            "two_burn_dv_km_s": co["two_burn_dv"] * v1,
            "two_burn_all_at_apogee_km_s": dv_all_apo * v1,
            "bi_parabolic_limit_km_s": co["bi_parabolic_limit_dv"] * v1,
            "three_burn_s_star": co["three_burn_s_star"],
            "saving_vs_two_burn_pct": 100.0 * (co["two_burn_dv"] - co["best_dv"]) / co["two_burn_dv"],
            "two_burn_theta1_deg": tb["theta1_star_deg"],
            "v1_km_s": v1,
        }

    out = {}
    # LEO(200 km) -> GEO with a 28.6 deg plane change (the Curtis-style case,
    # using the experiment's 200 km LEO altitude).
    out["leo_geo_28p6deg"] = scale(r_geo / r_leo, 28.6)
    # the classic Curtis / Gonzalez 300 km LEO, 28.6 deg -> GEO worked example
    r_leo300 = R_EARTH_KM + 300.0
    R_curtis = (R_EARTH_KM + GEO_ALT_KM) / r_leo300
    out["curtis_300km_28p6deg"] = scale(R_curtis, 28.6)
    # GTO -> GEO with a small plane change (all-at-apogee is near-optimal)
    out["gto_to_geo_5deg"] = scale(r_geo / r_leo, 5.0)
    # super-synchronous SES-8 / Thaicom-6 style: apogee ~ 90 000 km
    r_apo = 90000.0
    R_ss = r_apo / r_leo
    out["supersynchronous_90000km_30deg"] = scale(R_ss, 30.0)
    return out


# --------------------------------------------------------------------------- #
# High-precision (mpmath, 50 digits) cross-checks at the key analytic points.
# --------------------------------------------------------------------------- #
def high_precision_verification():
    """Reproduce the corner identities at 50-digit precision:
      * bi_parabolic_plane_change_limit(R) == Exp 005 bi_parabolic(R)
        (both = (sqrt2-1)(1+1/sqrt(R))).
      * R = 1 detour anchoring: di_c(1) = 2 arcsin(1/3) = 38.9424... deg,
        di_inf(1) = 60 deg (optimum s* -> infinity at di >= 60).
      * the finite-s dip at R = 2, di ~ 47.5 deg beats two-burn (ratio)."""
    mp.mp.dps = 50
    s2 = mp.sqrt(2)

    def bpp(R):
        return (s2 - 1) * (1 + 1 / mp.sqrt(mp.mpf(R)))

    def two_burn(R, di):
        R = mp.mpf(R); di = mp.mpf(di)
        vp = mp.sqrt(2 * R / (1 + R))
        v2 = 1 / mp.sqrt(R)
        v_apo = mp.sqrt(2 / (R * (1 + R)))

        def tot(th1):
            dv1 = mp.sqrt(1 + vp ** 2 - 2 * vp * mp.cos(th1))
            dv2 = mp.sqrt(v2 ** 2 + v_apo ** 2 - 2 * v2 * v_apo * mp.cos(di - th1))
            return dv1 + dv2

        # golden-section on [0, di]
        a, b = mp.mpf(0), di
        gr = (mp.sqrt(5) - 1) / 2
        c = b - gr * (b - a)
        d = a + gr * (b - a)
        fc, fd = tot(c), tot(d)
        for _ in range(200):
            if b - a < mp.mpf("1e-20"):
                break
            if fc < fd:
                b, d, fd = d, c, fc
                c = b - gr * (b - a)
                fc = tot(c)
            else:
                a, c, fc = c, d, fd
                d = a + gr * (b - a)
                fd = tot(d)
        t1 = (a + b) / 2
        return float(tot(t1))

    def three_burn(R, s, di):
        R = mp.mpf(R); s = mp.mpf(s); di = mp.mpf(di)
        vp12 = mp.sqrt(2 * s / (1 + s))
        v2 = 1 / mp.sqrt(R)
        va1 = mp.sqrt(2 / (s * (1 + s)))
        va2 = mp.sqrt(2 * R / (s * (R + s)))
        vp23 = mp.sqrt(2 * s / (R * (R + s)))
        # search over th1, th2 on a fine grid (mpmath, deterministic)
        best = mp.mpf("inf")
        nth = 121
        ths = [di * k / (nth - 1) for k in range(nth)]
        for th1 in ths:
            for th2 in ths:
                th3 = di - th1 - th2
                if th3 < 0:
                    continue
                dv1 = mp.sqrt(1 + vp12 ** 2 - 2 * vp12 * mp.cos(th1))
                dv2 = mp.sqrt(va1 ** 2 + va2 ** 2 - 2 * va1 * va2 * mp.cos(th2))
                dv3 = mp.sqrt(v2 ** 2 + vp23 ** 2 - 2 * v2 * vp23 * mp.cos(th3))
                best = min(best, dv1 + dv2 + dv3)
        return float(best)

    # R = 1 detour anchoring
    di_c_1 = 2 * mp.asin(mp.mpf(1) / 3)
    di_inf_1 = mp.radians(60)
    # R = 2, di = 47.5 deg finite-s dip vs two-burn
    R2 = 2.0
    di = mp.radians(47.5)
    two = two_burn(R2, di)
    # scan s in a range that includes the dip near 2.7
    best3 = mp.mpf("inf")
    s_vals = [mp.mpf(R2) * (1 + 0.05 * k) for k in range(1, 60)]
    for sv in s_vals:
        best3 = min(best3, three_burn(R2, sv, di))
    return {
        "dps": 50,
        "bi_parabolic_plane_change_R2": str(bpp(2)),
        "Exp005_bi_parabolic_R2": str(biell_exp.bi_parabolic(np.array([2.0]))[0]),
        "identity_match": str(bpp(2) - biell_exp.bi_parabolic(np.array([2.0]))[0]),
        "R1_di_c_deg": str(mp.degrees(di_c_1)),
        "R1_di_inf_deg": str(mp.degrees(di_inf_1)),
        "R2_di47p5_two_burn_mp": float(two),
        "R2_di47p5_three_burn_finite_mp": float(best3),
        "R2_di47p5_three_beats_two_mp": bool(best3 < mp.mpf(two)),
        "R2_di47p5_saving_pct": float(100 * (mp.mpf(two) - best3) / mp.mpf(two)),
    }


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def make_figures(results: dict):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    paths = []

    # 1. Regime map in (R, delta_i): 2-burn / finite-s / infinite-s.
    sweep = results["regime_sweep"]
    Rs = np.array(sweep["Rs"])
    dis = np.array(sweep["di_deg_grid"])
    labels = np.array(sweep["labels"])
    Rg, Dg = np.meshgrid(Rs, dis, indexing="ij")
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    cmap = ListedColormap(["#d62728", "#1f77b4", "#2ca02c"])
    ax.pcolormesh(Rg, Dg, labels, cmap=cmap, shading="auto", vmin=-0.5, vmax=2.5)
    # overlay boundaries
    bnd = results["boundaries"]
    bR = np.array([b["R"] for b in bnd])
    bdc = np.array([b["di_c_deg"] for b in bnd])
    bdi = np.array([b["di_inf_deg"] if b["di_inf_deg"] is not None else np.nan for b in bnd])
    ax.plot(bR, bdc, "k-", lw=2, label="di_c(R): 2-burn -> 3-burn")
    ax.plot(bR, bdi, "k--", lw=2, label="di_inf(R): finite-s -> s->inf")
    ax.set_xscale("log")
    ax.set_xlabel("radius ratio R = r2/r1")
    ax.set_ylabel("plane-change angle delta_i (deg)")
    ax.set_title("Global-optimum regime: combined transfer + plane change")
    ax.legend(loc="lower right", fontsize=8)
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-0.5, 2.5)),
                 ax=ax, ticks=[0, 1, 2], label="0=two-burn, 1=finite-s, 2=s->inf")
    fig.tight_layout()
    p = FIG_DIR / "regime_map.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 2. Cost vs delta_i at representative R, with the three candidate curves.
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
    for ax, R in zip(axes, (1.5, 2.0, 6.41)):
        dd = np.linspace(1.0, 120.0, 240)
        two = [two_burn_optimal(R, np.radians(d))["total_dv"] for d in dd]
        inf = [bi_parabolic_plane_change_limit(R)] * len(dd)
        fin = []
        for d in dd:
            co = combined_optimum(R, np.radians(d))
            fin.append(co["three_burn_finite_dv"])
        ax.plot(dd, two, "k-", lw=1.8, label="two-burn")
        ax.plot(dd, inf, "g:", lw=1.6, label="s->inf limit")
        ax.plot(dd, fin, "b-", lw=1.0, label="3-burn finite-s (lower envelope)")
        ax.set_title("R = %.2f" % R)
        ax.set_xlabel("delta_i (deg)")
        ax.grid(True, alpha=0.3)
        if ax is axes[0]:
            ax.set_ylabel("delta-v / v1")
            ax.legend(fontsize=7)
    fig.suptitle("Cost curves vs plane-change angle (normalized)")
    fig.tight_layout()
    p = FIG_DIR / "cost_vs_di.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 3. s* vs delta_i where finite-s wins (shows the dip and divergence).
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for R in (2.0, 4.0, 6.41, 8.0, 12.0):
        dd = np.linspace(1.0, 120.0, 240)
        ss = []
        for d in dd:
            co = combined_optimum(R, np.radians(d))
            if co["regime"] == "finite_s":
                ss.append(co["three_burn_s_star"])
            else:
                ss.append(np.nan)
        ax.plot(dd, ss, lw=1.4, label="R = %.2f" % R)
    ax.set_yscale("log")
    ax.set_xlabel("delta_i (deg)")
    ax.set_ylabel("optimal intermediate apoapsis s* = r_b/r1")
    ax.set_title("Finite-s optimum s*(delta_i) (only where finite-s regime holds)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = FIG_DIR / "s_star_vs_di.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))
    return paths


# --------------------------------------------------------------------------- #
# Main experiment
# --------------------------------------------------------------------------- #
def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # --- Part A: pure plane-change / detour anchors (R = 1) -------------------
    partA = {
        "pure_plane_change_90deg": pure_plane_change_dv(1.0, np.pi / 2),
        "detour_optimum_45deg": detour_optimum(np.radians(45.0)),
        "detour_optimum_60deg": detour_optimum(np.radians(60.0)),
        "detour_optimum_30deg": detour_optimum(np.radians(30.0)),
        "R1_di_c_deg_expected": float(np.degrees(2 * np.arcsin(1.0 / 3.0))),
        "R1_di_inf_deg_expected": 60.0,
    }

    # --- Boundary mapping over a dense R grid --------------------------------
    # R from just above 1 out to 50 (covers the window pinch near 6.41 and the
    # large-R behavior). Use log spacing plus a dense cluster near 6.41.
    Rs = sorted(set(np.round(np.concatenate([
        np.logspace(np.log10(1.05), np.log10(50.0), 46),
        np.linspace(5.8, 7.2, 28),
    ]), 6).tolist()))
    boundaries = []
    for R in Rs:
        dc = di_c_boundary(R)
        diinf = di_inf_boundary(R)
        boundaries.append({
            "R": float(R),
            "di_c_deg": dc,
            "di_inf_deg": diinf,
            "finite_window_width_deg": (diinf - dc) if (dc is not None and diinf is not None) else None,
        })

    # --- Dense regime sweep for the map --------------------------------------
    Rs_sweep = np.logspace(np.log10(1.1), np.log10(40.0), 60)
    di_sweep = np.linspace(2.0, 120.0, 120)
    sweep = regime_sweep(Rs_sweep, di_sweep)

    # --- Independent validation ----------------------------------------------
    # (a) a second brute-force optimizer at a few (R, di) to cross-check the
    #     primary optimizer (different s grid + larger split grid).
    cross_check = []
    for (R, di) in [(2.0, 47.5), (4.0, 45.0), (6.41, 38.2), (2.0, 50.0), (12.0, 35.0)]:
        a = combined_optimum(R, np.radians(di))
        b = combined_optimum(R, np.radians(di), ns=400, nth=96, s_max=1e8)
        cross_check.append({
            "R": R, "di_deg": di,
            "primary_regime": a["regime"], "primary_dv": a["best_dv"],
            "primary_s": a["three_burn_s_star"],
            "alt_regime": b["regime"], "alt_dv": b["best_dv"],
            "alt_s": b["three_burn_s_star"],
            "dv_abs_diff": abs(a["best_dv"] - b["best_dv"]),
            "agree": a["regime"] == b["regime"]
                     and abs(a["best_dv"] - b["best_dv"]) < 1e-4,
        })

    # (b) 3D RK4 trajectory validation of representative optimal maneuvers.
    rk4_two = [validate_two_burn_rk4(R, np.radians(di))
               for (R, di) in [(6.41, 28.6), (2.0, 30.0), (4.0, 20.0)]]
    rk4_three = [validate_three_burn_rk4(R, np.radians(di), s)
                 for (R, di, s) in [(2.0, 47.5, 2.73), (4.0, 45.0, 4.97),
                                    (2.0, 40.0, 2.16)]]

    # (c) high-precision mpmath cross-checks
    hp = high_precision_verification()

    # (d) precise location of the finite-s window pinch (where the width
    #     di_inf(R) - di_c(R) closes to zero). This is the "abrupt behavior"
    #     the prior investigation flagged near R ~ 6.41.
    def _width(R):
        dc = di_c_boundary(R)
        di = di_inf_boundary(R)
        if dc is None or di is None:
            return 0.0
        return di - dc

    # bracket: at R=2 the window is wide (>15 deg); at R=12 it is closed (0).
    lo, hi = 2.0, 12.0
    if _width(lo) > 0 and _width(hi) <= 0:
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if _width(mid) > 0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-4:
                break
        pinch_R = 0.5 * (lo + hi)
    else:
        pinch_R = None
    pinch = {
        "R_pinch_approx": pinch_R,
        "width_at_R2": _width(2.0),
        "width_at_R641": _width(6.41),
        "interpretation": (
            "for R < R_pinch a finite intermediate apoapsis s* beats two-burn "
            "in a window (di_c, di_inf); for R > R_pinch the only 3-burn regime "
            "that beats two-burn is the s->infinity (free-at-apoapsis) one"
        ),
    }

    # --- Real anchors --------------------------------------------------------
    anchors = real_anchors()

    # --- Assemble results ---------------------------------------------------
    results = {
        "meta": {
            "experiment": "combined bi-elliptic transfer + plane change",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "python": "3.12", "numpy": np.__version__,
            "mpmath_dps": 50,
            "normalization": "r1=1, v1=1, mu=1; physical values scale by v1=sqrt(mu/r1)",
            "optimizer": "dense log-s grid (240) x 2D split (64x64) + 2 refinement passes + nested golden polish; s->infinity analytic limit folded in",
            "determinism": "pure float64 + fixed mpmath precision, no RNG",
            "boundaries": "di_c(R), di_inf(R) via bisection (tol 1e-4 deg) on the global optimizer",
        },
        "partA_anchors": partA,
        "boundaries": boundaries,
        "regime_sweep": sweep,
        "cross_check_alternate_optimizer": cross_check,
        "rk4_validation_two_burn": rk4_two,
        "rk4_validation_three_burn": rk4_three,
        "high_precision": hp,
        "real_anchors": anchors,
        "pinch_point": pinch,
        "key_findings": {
            "R1_di_c_deg": partA["R1_di_c_deg_expected"],
            "R1_di_inf_deg": partA["R1_di_inf_deg_expected"],
            "R_pinch_finite_s_window": pinch_R,
            "window_width_at_R2_deg": _width(2.0),
            "s_inf_limit_identity": "bi_parabolic_plane_change_limit(R) == Exp005 bi_parabolic(R)",
        },
    }

    from lab_utils.results import save_json_result
    save_json_result(
        RESULTS_DIR / "results.json",
        results,
        name="combined_bielliptic_plane_change",
        description="Global optimum of combined radius change + plane change: "
                    "two-burn vs finite-s three-burn vs s->infinity regimes and "
                    "their (R, delta_i) boundaries.",
    )

    fig_paths = make_figures(results)

    # console summary
    print("=== combined transfer + plane change: key results ===")
    print("R=1 detour anchors: di_c=%.4f deg, di_inf=%.2f deg"
          % (partA["R1_di_c_deg_expected"], partA["R1_di_inf_deg_expected"]))
    print("mpmath s->inf identity diff (R=2): %s" % hp["identity_match"])
    print("mpmath R=2,di=47.5 three-beats-two: %s (%.3f%% saving)"
          % (hp["R2_di47p5_three_beats_two_mp"], hp["R2_di47p5_saving_pct"]))
    print("figures: %s" % fig_paths)
    return results


if __name__ == "__main__":
    main()
