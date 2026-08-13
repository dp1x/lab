"""Bi-elliptic transfer vs Hohmann: the three-impulse crossover structure.

Research question: between two coplanar circular orbits of radii r1 < r2 about
a central body mu, when is the three-impulse bi-elliptic transfer (r1 -> r_b ->
r2, r_b > r2) cheaper than the Hohmann transfer, as a function of the radius
ratio R = r2/r1 and the intermediate apoapsis ratio s = r_b/r1?

The three burns, from the vis-viva equation v^2 = mu(2/r - 1/a) (ellipse 1:
periapsis r1, apoapsis r_b, a1 = (r1+r_b)/2; ellipse 2: periapsis r2,
apoapsis r_b, a2 = (r2+r_b)/2), in units of v1 = sqrt(mu/r1):

  dv1/v1 = sqrt(2s/(1+s)) - 1                                (at r1)
  dv2/v1 = (1/sqrt(s)) * [ sqrt(2R/(R+s)) - sqrt(2/(1+s)) ]  (at r_b)
  dv3/v1 = (1/sqrt(R)) * [ sqrt(2s/(R+s)) - 1 ]              (at r2)

Key analytic results derived and verified in this experiment:

  1. Degeneracy: as s -> R the bi-elliptic reduces to the Hohmann transfer
     (dv3 -> 0, dv1, dv2 -> the Hohmann burns): the bi-elliptic family is the
     continuous generalization of Hohmann with one extra impulse.
  2. Bi-parabolic limit (s -> infinity, the two ellipses -> parabolas):
     dv/v1 -> (sqrt(2)-1)(1 + 1/sqrt(R)); it beats Hohmann iff R > R_bp with
     R_bp ~ 11.9388 (the classic "11.94" crossover of Hoelker & Silber 1959).
  3. Corner identity: d/ds [dv_biell/v1] at s = R equals d/dR [dv_H/v1]
     (proved by algebra in the card, verified numerically and at 50-digit
     precision here). The Hohmann cost curve has its maximum at R* ~ 15.5817
     (Experiment 004), so: for R > R* the bi-elliptic is cheaper than Hohmann
     for EVERY finite s > R ("always cheaper" regime, Hoelker & Silber), for
     R < R_bp Hohmann is cheaper for every s, and in (R_bp, R*) there is a
     unique crossover s_c(R): bi-elliptic wins iff s > s_c(R), with
     s_c -> infinity as R -> R_bp+ and s_c -> R as R -> R*-.
  4. The crossover curve reproduces the classical table (Escobal 1968;
     Gobetz & Doll 1969): s_c(12) ~ 815.8, s_c(13) ~ 48.9, s_c(14) ~ 26.1,
     s_c(15) ~ 18.2, s_c(15.58) ~ 15.58. (One widely copied table entry,
     "12 -> 15.81", is shown by independent computation and 50-digit
     arithmetic to be a dropped-digit transcription error for 815.81.)
  5. The low-family bi-elliptic (intermediate radius BELOW both orbits,
     s < min(1, R)) is never cheaper than Hohmann for any R > 1; inward
     transfers are governed by the same thresholds via time-reversal
     (bi-elliptic wins inward iff outer/inner > R_bp with s large enough).
  6. The fuel saving over Hohmann peaks near R ~ 50 (dv/v1 ~ 0.041, ~4 % of
     the inner circular speed) and decays as (2-sqrt(2))/sqrt(R) at large R.
  7. Flight time is strictly worse: t = pi*sqrt((r1+r_b)^3/(8 mu)) +
     pi*sqrt((r2+r_b)^3/(8 mu)) grows ~ s^{3/2} against the Hohmann half
     period; the bi-parabolic limit has infinite transfer time.

The full three-burn trajectory is propagated with the verified RK4 machinery
of Experiment 002 (loaded through Experiment 004, same explicit-path
importlib pattern) to validate every closed-form quantity. The Hohmann
reference machinery and IAU nominal constants are reused from Experiment 004.

References:
  - R. F. Hoelker, R. Silber, "The bi-elliptical transfer between co-planar
    circular orbits", Planetary and Space Science 7, 164-175 (1961) (the
    primary source; R > 15.58 "always" claim, 11.94 crossover).
  - P. R. Escobal, "Methods of Astrodynamics", Wiley 1968 (crossover table).
  - F. W. Gobetz, J. R. Doll, "A Survey of Impulsive Trajectories",
    AIAA Journal 7(5), 801-834 (1969) (the table source cited by Wikipedia).
  - R. R. Bate, D. D. Mueller, J. E. White, "Fundamentals of Astrodynamics",
    Dover 1971, Ch. 6.
  - H. D. Curtis, "Orbital Mechanics for Engineering Students", 4th ed.,
    Elsevier 2021, Ch. 6 (worked bi-elliptic example numbers).
  - D. A. Vallado, "Fundamentals of Astrodynamics and Applications", 4th ed.,
    Microcosm 2013, Sec. 6.5.
  - Wikipedia "Bi-elliptic transfer" (rev. 1233203053, cited Curtis 2005 and
    Vallado 2001): formula set, 11.94/15.58 structure and the worked R = 14,
    r_b = 40 r0 example (3061.04 / 608.825 / 447.662 m/s, total 4117.53 m/s)
    - reproduced exactly here.
  - "Orbital Mechanics & Astrodynamics" open textbook (orbital-mechanics.space,
    bi-elliptic Hohmann transfer section): 11.94/15.58 statements.
  - poliastro docs, "Comparing Hohmann and bielliptic transfers": the
    15.58 boundary coincides with the Hohmann cost maximum (R[idx_max]).
  - IAU 2015 Resolution B3 nominal constants (Mamajek et al. arXiv:1510.07674)
    and IAU 2012 Resolution B2 au - via Experiment 004.

Structure of the results:

  - theory_constants: the two boundary radius ratios (R_bp ~ 11.9388,
    R* ~ 15.5817) and the corner identity verification.
  - crossover_table: s_c(R) on (R_bp, R*), plus the classical table points.
  - region_verification: adversarial grid checks of the three-regime claims.
  - shape_diagnostics: the f(s) shape (monotone increasing / single hump /
    monotone decreasing) across R.
  - saving_curve: dv_H - dv_bi-parabolic, its peak and the 1/sqrt(R)
    asymptote; time penalty t_biell/t_H.
  - low_family / inward: the never-wins low family and the inward
    time-reversal equivalence.
  - rk4_validation: full three-burn propagation of six (R, s) cases.
  - high_precision: mpmath 50-digit cross-checks of the roots and the
    corner identity.
  - real_anchors: LEO -> GEO, LEO -> 14 x LEO (Wikipedia example), LEO ->
    50 x LEO (max-saving ratio), GEO -> lunar, GEO -> 15.58 x GEO, Earth ->
    Mars, with delta-v, saving and flight-time budgets.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from lab_utils.results import save_json_result

# --- Reuse of the verified 004 machinery ----------------------------------
#
# Experiment 004 (hohmannTransfer) owns the verified Hohmann closed forms,
# the digit-safe variants, the (R-1)/2 and sqrt(2)-1 asymptote checks and the
# IAU nominal constants; it in turn loads Experiment 002's verified RK4
# propagator / Kepler solver / elements machinery (kov). Loading 004 by its
# explicit path (the same importlib pattern the test suites use) gives this
# experiment a single source of truth for the Hohmann reference costs.

_ht_path = Path(__file__).resolve().parents[1] / "hohmannTransfer" / "experiment.py"
_ht_spec = importlib.util.spec_from_file_location("hohmann_transfer_exp004", _ht_path)
assert _ht_spec is not None and _ht_spec.loader is not None
ht = importlib.util.module_from_spec(_ht_spec)
_ht_spec.loader.exec_module(ht)

hohmann_dv_total = ht.hohmann_dv_total
hohmann_transfer_time = ht.hohmann_transfer_time
hohmann_split = ht.hohmann_split
dv_over_v1 = ht.dv_over_v1
peak_of_cost_curve = ht.peak_of_cost_curve
kov = ht.kov
propagate_rk4 = kov.propagate_rk4
kepler_solution = kov.kepler_solution
orbital_elements = kov.orbital_elements

MU = 1.0  # canonical gravitational parameter for the core studies
MU_EARTH_KM3S2 = ht.MU_EARTH_KM3S2
R_EARTH_KM = ht.R_EARTH_KM
LEO_ALT_KM = ht.LEO_ALT_KM
GEO_ALT_KM = ht.GEO_ALT_KM
AU_KM = ht.AU_KM
MU_SUN_KM3S2 = ht.MU_SUN_KM3S2
MARS_A_AU = ht.MARS_A_AU
G0 = ht.G0

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

S2_MINUS_1 = np.sqrt(2.0) - 1.0  # the escape burn, in units of the local v_circ

# --- Closed-form bi-elliptic transfer -------------------------------------


def ellipse_speed(r: float, a: float, mu: float = MU) -> float:
    """Vis-viva speed at radius r on an orbit of semi-major axis a."""
    return np.sqrt(mu * (2.0 / r - 1.0 / a))


def bielliptic_burns(
    r1: float, r2: float, r_b: float, mu: float = MU
) -> tuple[float, float, float]:
    """The three impulse magnitudes of a bi-elliptic transfer between
    circular orbits r1 and r2 via the common apsidal radius r_b.

    r_b may lie above both orbits (the high family: both transfer ellipses
    have apoapsis r_b) or below both (the low family: both have periapsis
    r_b); the formula is direction- and family-agnostic. Burn 1 and burn 3
    are the speed mismatch between the circular orbit and the transfer
    ellipse at the departure/arrival radius; burn 2 is the tangential speed
    difference between the two transfer ellipses at r_b.
    """
    a1 = 0.5 * (r1 + r_b)
    a2 = 0.5 * (r2 + r_b)
    v1c = np.sqrt(mu / r1)
    v2c = np.sqrt(mu / r2)
    dv1 = abs(ellipse_speed(r1, a1, mu) - v1c)
    dv2 = abs(ellipse_speed(r_b, a2, mu) - ellipse_speed(r_b, a1, mu))
    dv3 = abs(ellipse_speed(r2, a2, mu) - v2c)
    return dv1, dv2, dv3


def bielliptic_dv_total(r1: float, r2: float, r_b: float, mu: float = MU) -> float:
    return sum(bielliptic_burns(r1, r2, r_b, mu))


def bielliptic_transfer_time(r1: float, r2: float, r_b: float, mu: float = MU) -> float:
    """Two half-periods: periapsis -> apoapsis on ellipse 1, apoapsis ->
    periapsis on ellipse 2."""
    a1 = 0.5 * (r1 + r_b)
    a2 = 0.5 * (r2 + r_b)
    return np.pi * (np.sqrt(a1**3 / mu) + np.sqrt(a2**3 / mu))


def f_high(R: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Normalized high-family bi-elliptic cost in units of v1 = sqrt(mu/r1),
    r1 = 1, r2 = R > 1, r_b = s > R:

        dv1 = sqrt(2s/(1+s)) - 1
        dv2 = (1/sqrt(s)) [ sqrt(2R/(R+s)) - sqrt(2/(1+s)) ]
        dv3 = (1/sqrt(R)) [ sqrt(2s/(R+s)) - 1 ]
    """
    return (
        (np.sqrt(2.0 * s / (1.0 + s)) - 1.0)
        + (1.0 / np.sqrt(s))
        * (np.sqrt(2.0 * R / (R + s)) - np.sqrt(2.0 / (1.0 + s)))
        + (1.0 / np.sqrt(R)) * (np.sqrt(2.0 * s / (R + s)) - 1.0)
    )


def f_low(R: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Normalized low-family bi-elliptic cost (intermediate radius below both
    orbits, s < min(1, R)), r1 = 1, r2 = R, r_b = s:

        dv1 = 1 - sqrt(2s/(1+s))                     (brake at apoapsis r1)
        dv2 = (1/sqrt(s)) [ sqrt(2R/(R+s)) - sqrt(2/(1+s)) ]   (boost at r_b)
        dv3 = (1/sqrt(R)) [ 1 - sqrt(2s/(R+s)) ]     (raise at apoapsis r2)

    The deep-space term is a BOOST: the periapsis speed of ellipse 2
    (periapsis s, apoapsis R) exceeds that of ellipse 1 (periapsis s,
    apoapsis 1) because a2 > a1; at the degenerate corner s -> min(1, R) the
    sum reduces exactly to the Hohmann cost dv_H(R) (ellipse 2 becomes the
    Hohmann ellipse and ellipse 1 the initial circle).
    """
    return (
        (1.0 - np.sqrt(2.0 * s / (1.0 + s)))
        + (1.0 / np.sqrt(s))
        * (np.sqrt(2.0 * R / (R + s)) - np.sqrt(2.0 / (1.0 + s)))
        + (1.0 / np.sqrt(R)) * (1.0 - np.sqrt(2.0 * s / (R + s)))
    )


def bi_parabolic(R: np.ndarray) -> np.ndarray:
    """Bi-parabolic limit of the bi-elliptic cost (s -> infinity): the first
    burn tends to the escape burn and the third to (sqrt(2)-1)/sqrt(R); the
    deep-space burn vanishes."""
    return S2_MINUS_1 * (1.0 + 1.0 / np.sqrt(R))


# --- Digit-safe forms near the Hohmann corner (s -> R) ---------------------


def burn1_stable(s: np.ndarray) -> np.ndarray:
    """sqrt(2s/(1+s)) - 1 = (s-1) / ((1+s)(1 + sqrt(2s/(1+s)))) - free of
    cancellation at s -> 1."""
    return (s - 1.0) / ((1.0 + s) * (1.0 + np.sqrt(2.0 * s / (1.0 + s))))


def burn3_stable(R: np.ndarray, s: np.ndarray) -> np.ndarray:
    """(1/sqrt(R)) (sqrt(2s/(R+s)) - 1) = (1/sqrt(R)) (s-R) /
    ((R+s)(1 + sqrt(2s/(R+s)))) - free of cancellation at s -> R."""
    return (
        (1.0 / np.sqrt(R))
        * (s - R)
        / ((R + s) * (1.0 + np.sqrt(2.0 * s / (R + s))))
    )


def burn2_stable(R: np.ndarray, s: np.ndarray) -> np.ndarray:
    """The deep-space burn with the difference of the two square roots
    rearranged to (2 s (R-1)) / ((R+s)(1+s)(sqrt(2R/(R+s)) + sqrt(2/(1+s))))
    - free of cancellation at s -> R."""
    return (
        (1.0 / np.sqrt(s))
        * (2.0 * s * (R - 1.0))
        / (
            (R + s)
            * (1.0 + s)
            * (np.sqrt(2.0 * R / (R + s)) + np.sqrt(2.0 / (1.0 + s)))
        )
    )


def f_high_stable(R: np.ndarray, s: np.ndarray) -> np.ndarray:
    """f_high evaluated without cancellation anywhere in s >= R."""
    return burn1_stable(s) + burn2_stable(R, s) + burn3_stable(R, s)


# --- Crossover machinery ---------------------------------------------------


def d_dv_ho_dR(R: np.ndarray) -> np.ndarray:
    """Analytic derivative of the Hohmann cost dv_H/v1 with respect to R:

        d/dR = (1/2) sqrt(2R/(1+R)) / (R(1+R))
               - (1/2) R^{-3/2} (1 - sqrt(2/(1+R)))
               + (1/2) R^{-1/2} sqrt(2/(1+R)) / (1+R)

    Its zero is the Hohmann cost maximum (R* ~ 15.58, Experiment 004).
    """
    A = 0.5 * np.sqrt(2.0 * R / (1.0 + R)) / (R * (1.0 + R))
    B = 0.5 * R ** (-1.5) * (1.0 - np.sqrt(2.0 / (1.0 + R)))
    C = 0.5 * R ** (-0.5) * np.sqrt(2.0 / (1.0 + R)) / (1.0 + R)
    return A - B + C


def _bisect_root(f, lo: float, hi: float, iters: int = 200) -> float:
    flo, fhi = f(lo), f(hi)
    assert flo * fhi < 0.0, (lo, hi, flo, fhi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm > 0.0:
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def R_bp_crossover() -> float:
    """Root of bi_parabolic(R) - dv_H(R) = 0 in (1, 100): the radius ratio
    where the bi-parabolic limit ties the Hohmann transfer (~11.9388)."""
    return _bisect_root(
        lambda R: float(bi_parabolic(R) - dv_over_v1(R)), 1.01, 100.0
    )


def R_always_cheaper() -> float:
    """Root of d/dR dv_H = 0 in (10, 30): the Hohmann cost maximum, which is
    also the boundary R* above which the bi-elliptic is cheaper than Hohmann
    for every intermediate apoapsis (corner identity, Sec. 3 of the card)."""
    return _bisect_root(d_dv_ho_dR, 10.0, 30.0)


def crossover_s(R: float, s_hi: float = 1e12) -> float:
    """The intermediate-apoapsis ratio s_c(R) at which the high-family
    bi-elliptic cost equals the Hohmann cost, for R in (R_bp, R*).

    g = f_high - dv_H is positive at s = R (rising corner, R < R*) and
    negative at s = infinity (R > R_bp); bisection in log s.
    """
    lo = np.log(R * 1.0000001)
    hi = np.log(s_hi)
    glo = f_high_stable(np.array([R]), np.exp(np.array([lo]))) - dv_over_v1(
        np.array([R])
    )
    glo = float(glo[0])
    ghi = float(
        (f_high(np.array([R]), np.array([s_hi])) - dv_over_v1(np.array([R])))[0]
    )
    assert glo * ghi < 0.0, (R, glo, ghi)
    for _ in range(250):
        mid = 0.5 * (lo + hi)
        gm = float(
            (
                f_high_stable(np.array([R]), np.exp(np.array([mid])))
                - dv_over_v1(np.array([R]))
            )[0]
        )
        if glo * gm > 0.0:
            lo = mid
        else:
            hi = mid
    return float(np.exp(0.5 * (lo + hi)))


def crossover_curve(
    R_bp: float, R_star: float, n: int = 61, s_hi: float = 1e12
) -> dict:
    """s_c(R) on (R_bp, R*): grid uniform in R, plus the two endpoints whose
    limits are known analytically (infinity and R*)."""
    R_grid = np.linspace(R_bp + 1e-4, R_star - 1e-4, n)
    s_c = np.array([crossover_s(float(R), s_hi) for R in R_grid])
    # classical table points (Escobal 1968 / Gobetz & Doll 1969)
    table_R = np.array([12.0, 13.0, 14.0, 15.0, 15.58])
    table_s = np.array([crossover_s(float(R)) for R in table_R])
    return {
        "R_grid": R_grid.tolist(),
        "s_c_grid": s_c.tolist(),
        "endpoint_R_bp": {"R": R_bp, "s_c": float("inf")},
        "endpoint_R_star": {"R": R_star, "s_c": R_star},
        "classical_table_R": table_R.tolist(),
        "classical_table_s_c": table_s.tolist(),
        "note_wikipedia_R12_entry": (
            "independent computation gives s_c(12) ~ 815.82; the widely "
            "copied table entry '12 -> 15.81' (Wikipedia, from Gobetz & Doll "
            "1969) is a dropped-digit transcription error for 815.81, "
            "consistent with this experiment's float64 and 50-digit results"
        ),
    }


# --- Region verification (adversarial checks of the textbook claims) -------


def region_verification(
    R_bp: float, R_star: float
) -> dict:
    """Grid attacks on the three claims:

      (i)   R < R_bp:  g = f_high - dv_H > 0 for every s in (R, 1e6 R)
      (ii)  R > R_star: g < 0 for every s in (R, 1e6 R)
      (iii) R_bp < R < R_star: exactly one crossing of g in log s

    The corner s = R itself (g = 0) is excluded; the nearest grid line sits
    at (1 + 1e-6) R.
    """
    Rs_lo = np.logspace(np.log10(1.001), np.log10(R_bp * 0.9999), 90)
    Rs_hi = np.logspace(np.log10(R_star * 1.0001), 3.0, 90)
    Rs_mid = np.linspace(R_bp * 1.0001, R_star * 0.9999, 60)

    def scan(Rs: np.ndarray, sign_expected: float) -> dict:
        worst = float("inf")
        worst_at = None
        crossings = []
        for R in Rs:
            s = np.logspace(np.log10(R * 1.000001), np.log10(R * 1e6), 400)
            g = f_high_stable(R * np.ones_like(s), s) - dv_over_v1(R * np.ones_like(s))
            if sign_expected > 0.0:
                assert np.all(g > -1e-12), (R, g.min())
                if g.min() < worst:
                    worst, worst_at = float(g.min()), float(R)
            else:
                assert np.all(g < 1e-12), (R, g.max())
                if -g.min() < worst:
                    worst, worst_at = float(-g.min()), float(R)
            ncross = int(np.sum(np.sign(g[:-1]) * np.sign(g[1:]) < 0))
            crossings.append(ncross)
        return {
            "n_R": len(Rs),
            "n_s_per_R": 400,
            "worst_margin": worst,          # closest approach to g = 0
            "worst_margin_at_R": worst_at,
            "crossings_counted": crossings,
        }

    lo_scan = scan(Rs_lo, +1.0)
    hi_scan = scan(Rs_hi, -1.0)
    mid_crossings = []
    for R in Rs_mid:
        s = np.logspace(np.log10(R * 1.000001), np.log10(R * 1e7), 800)
        g = f_high_stable(R * np.ones_like(s), s) - dv_over_v1(R * np.ones_like(s))
        ncross = int(np.sum(np.sign(g[:-1]) * np.sign(g[1:]) < 0))
        mid_crossings.append(ncross)
    return {
        "R_lo_region": lo_scan,
        "R_hi_region": hi_scan,
        "mid_crossing_counts": {
            "n_R": len(Rs_mid),
            "max_crossings_seen": max(mid_crossings),
            "min_crossings_seen": min(mid_crossings),
            "all_exactly_one": all(c == 1 for c in mid_crossings),
        },
    }


# --- Shape diagnostics of f(s) ---------------------------------------------


def shape_diagnostics() -> dict:
    """Classify the shape of f(s) as a function of s for fixed R:

      monotone increasing  (R below the hump onset R_m)
      single interior hump (R_m < R < R_star)
      monotone decreasing  (R > R_star)

    plus, in the hump regime, confirmation that f decreases monotonically
    after the hump (this is what makes the crossover unique).
    """
    Rs = np.geomspace(1.2, 200.0, 80)
    shapes = []
    hump_onset = None
    for R in Rs:
        s = np.logspace(np.log10(R * 1.000001), np.log10(R * 1e8), 4000)
        f = f_high(R * np.ones_like(s), s)
        df = np.diff(f)
        rising = np.all(df > 0.0)
        falling = np.all(df < 0.0)
        nmin = int(np.sum((df[:-1] < 0) & (df[1:] > 0)))
        nmax = int(np.sum((df[:-1] > 0) & (df[1:] < 0)))
        if rising:
            kind = "monotone_increasing"
        elif falling:
            kind = "monotone_decreasing"
        else:
            kind = "single_hump" if (nmax == 1 and nmin == 0) else "other"
        if hump_onset is None and kind == "single_hump":
            hump_onset = float(R)
        # monotone decreasing after the hump?
        if kind == "single_hump":
            i_hump = int(np.argmax(f))
            assert np.all(np.diff(f[i_hump:]) < 0.0), R
        shapes.append({"R": float(R), "shape": kind, "n_max": nmax, "n_min": nmin})
    return {
        "shapes": shapes,
        "hump_onset_R_approx": hump_onset,
        "note": (
            "the hump exists for intermediate R but never dips below dv_H "
            "for R < R_bp; the unique crossing for R_bp < R < R_star follows "
            "from the hump-plus-monotone-decrease structure"
        ),
    }


# --- Saving curve and time penalty -----------------------------------------


def saving_and_time() -> dict:
    """dv_H - dv_bi-parabolic (the asymptotic fuel saving) and the transfer
    time ratio t_biell/t_H."""
    R_save = np.logspace(np.log10(1.0001), 4.0, 200001)
    saving = dv_over_v1(R_save) - bi_parabolic(R_save)
    i = int(np.argmax(saving))
    R_peak = float(R_save[i])
    peak = float(saving[i])
    # large-R asymptote: dv_H - dv_bp ~ (2 - sqrt(2))/sqrt(R)
    R_asym = 1e4
    asym_pred = (2.0 - np.sqrt(2.0)) / np.sqrt(R_asym)
    asym_meas = float(
        (dv_over_v1(np.array([R_asym])) - bi_parabolic(np.array([R_asym])))[0]
    )
    # time ratio for representative (R, s)
    R_t = np.array([2.0, 6.41, 14.0, 20.0])
    s_t = np.array([3.0, 30.0, 40.0, 100.0])
    t_ho = np.pi * np.sqrt((1.0 + R_t) ** 3 / 8.0)
    t_be = np.pi * (np.sqrt((1.0 + s_t) ** 3 / 8.0) + np.sqrt((R_t + s_t) ** 3 / 8.0))
    return {
        "peak_saving_over_v1": peak,
        "peak_saving_at_R": R_peak,
        "large_R_asymptote": "(2 - sqrt(2))/sqrt(R)",
        "asymptote_pred_at_1e4": asym_pred,
        "asymptote_meas_at_1e4": asym_meas,
        "asymptote_ratio": asym_meas / asym_pred,
        "time_ratio_cases": [
            {"R": float(R_t[k]), "s": float(s_t[k]),
             "t_biell_over_t_hohmann": float(t_be[k] / t_ho[k])}
            for k in range(len(R_t))
        ],
    }


# --- Low family and inward case --------------------------------------------


def low_family_and_inward(R_bp: float, R_star: float) -> dict:
    """(a) Low-family bi-elliptic (r_b below both orbits): verify it never
    beats Hohmann for R > 1.  (b) Inward transfers: the cost of the inward
    high-family bi-elliptic equals the outward high-family cost of the same
    two ellipses (time reversal), so the same thresholds apply with the roles
    of the inner/outer radius fixed: bi-elliptic wins inward iff
    r_outer/r_inner > R_bp (large s)."""
    # (a) low family grid: R in (1.001, 100), s in (1e-6 R, min(1, R))
    worst_low = 0.0
    R_grid = np.geomspace(1.001, 100.0, 40)
    for R in R_grid:
        s = np.logspace(np.log10(R * 1e-6), np.log10(min(1.0, R) * 0.9999), 200)
        g = f_low(R * np.ones_like(s), s) - dv_over_v1(R * np.ones_like(s))
        assert np.all(g > -1e-12), R
        worst_low = max(worst_low, float(g.min()))
    # (b) inward: r1 = 20, r2 = 1 (R_in = 0.05, r_outer/r_inner = 20 > 11.94)
    #     high family r_b = s r1, s > 1; compare with inward Hohmann.
    r1, r2 = 20.0, 1.0
    dv_ho_in = hohmann_dv_total(r1, r2, MU)
    cases = []
    for s in (1.1, 2.0, 5.0, 20.0, 100.0):
        r_b = s * r1
        dv_be = bielliptic_dv_total(r1, r2, r_b, MU)
        cases.append(
            {"s": s, "dv_bielliptic": dv_be, "dv_hohmann": dv_ho_in,
             "wins": bool(dv_be < dv_ho_in)}
        )
    # time-reversal magnitude identity at three radii sets: reversing the
    # transfer swaps the first and third burns (the same two ellipses, the
    # burns applied in the opposite order), so fwd[k] == bwd[2-k].
    rev_ok = True
    for (a, b, c) in ((1.0, 14.0, 40.0), (1.0, 20.0, 100.0), (14.0, 1.0, 0.025)):
        fwd = bielliptic_burns(a, b, c, MU)
        bwd = bielliptic_burns(b, a, c, MU)
        rev_ok = rev_ok and all(
            abs(fwd[k] - bwd[2 - k]) < 1e-12 for k in range(3)
        )
    return {
        "low_family_never_wins": {
            "grid_R": R_grid.tolist(),
            "worst_margin_g_min": worst_low,
        },
        "inward_high_family": {
            "r1_outer": r1,
            "r2_inner": r2,
            "R_in": r2 / r1,
            "dv_hohmann_inward": dv_ho_in,
            "cases": cases,
            "note": (
                "inward bi-elliptic wins when r_outer/r_inner > R_bp and s is "
                "large enough - the same thresholds as outward, by "
                "time-reversal symmetry"
            ),
        },
        "time_reversal_burn_identity": bool(rev_ok),
        "crossover_thresholds": {
            "R_bp": R_bp,
            "R_star": R_star,
            "inward_equivalent_outer_over_inner_R_bp": R_bp,
        },
    }


# --- Corner identity verification ------------------------------------------


def corner_identity_check(Rs: np.ndarray) -> dict:
    """Verify the corner identity numerically:

        d/ds f_high(R, s) |_{s=R}  ==  d/dR dv_H(R)

    (finite differences; analytic proof in the card). This identity is what
    makes the 'always cheaper' boundary coincide with the Hohmann maximum.
    R* itself is excluded: d/dR dv_H = 0 there, so the relative form is
    ill-conditioned (the 50-digit identity below covers R* instead).
    """
    out = []
    for R in Rs:
        h = 1e-6 * R
        dds = (f_high_stable(np.array([R]), np.array([R + h]))
               - f_high_stable(np.array([R]), np.array([R - h]))) / (2.0 * h)
        ddr = d_dv_ho_dR(np.array([R]))
        out.append({"R": float(R), "d_f_ds_at_corner": float(dds[0]),
                    "d_dvH_dR": float(ddr[0]),
                    "rel_diff": abs(dds[0] - ddr[0]) / abs(ddr[0])})
    return {"checks": out, "max_rel_diff": max(o["rel_diff"] for o in out)}


# --- High-precision verification (mpmath, 50 digits) -----------------------


def high_precision_verification() -> dict:
    """Recompute the two boundary ratios and the crossover curve table with
    mpmath at 50-digit precision; cross-check the float64 values."""
    import mpmath as mp

    mp.mp.dps = 50
    s2 = mp.sqrt(2)

    def dvH(R: mp.mpf) -> mp.mpf:
        R = mp.mpf(R)
        return (
            mp.sqrt(2 * R / (1 + R)) - 1
            + (1 / mp.sqrt(R)) * (1 - mp.sqrt(2 / (1 + R)))
        )

    def dvbp(R: mp.mpf) -> mp.mpf:
        return (s2 - 1) * (1 + 1 / mp.sqrt(R))

    def dH(R: mp.mpf) -> mp.mpf:
        R = mp.mpf(R)
        A = mp.mpf("0.5") * mp.sqrt(2 * R / (1 + R)) / (R * (1 + R))
        B = mp.mpf("0.5") * R ** mp.mpf("-1.5") * (
            1 - mp.sqrt(2 / (1 + R))
        )
        C = mp.mpf("0.5") * R ** mp.mpf("-0.5") * mp.sqrt(2 / (1 + R)) / (1 + R)
        return A - B + C

    def root(f, lo, hi):
        lo, hi = mp.mpf(lo), mp.mpf(hi)
        flo, fhi = f(lo), f(hi)
        assert flo * fhi < 0
        for _ in range(400):
            mid = (lo + hi) / 2
            if flo * f(mid) > 0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    R_bp_mp = root(lambda x: dvbp(x) - dvH(x), 11.9, 12.0)
    R_star_mp = root(dH, 15.5, 15.7)

    def fbe(R: mp.mpf, s: mp.mpf) -> mp.mpf:
        R, s = mp.mpf(R), mp.mpf(s)
        return (
            mp.sqrt(2 * s / (1 + s)) - 1
            + (1 / mp.sqrt(s)) * (
                mp.sqrt(2 * R / (R + s)) - mp.sqrt(2 / (1 + s))
            )
            + (1 / mp.sqrt(R)) * (mp.sqrt(2 * s / (R + s)) - 1)
        )

    def sc(R: mp.mpf) -> mp.mpf:
        lo, hi = mp.log(R * (1 + mp.mpf("1e-14"))), mp.log(mp.mpf(10) ** 40)
        glo = fbe(R, mp.e**lo) - dvH(R)
        ghi = fbe(R, mp.e**hi) - dvH(R)
        assert glo * ghi < 0
        for _ in range(500):
            mid = (lo + hi) / 2
            if glo * (fbe(R, mp.e**mid) - dvH(R)) > 0:
                lo = mid
            else:
                hi = mid
        return mp.e ** ((lo + hi) / 2)

    tbl = {12.0: sc(mp.mpf(12)), 13.0: sc(mp.mpf(13)),
           14.0: sc(mp.mpf(14)), 15.0: sc(mp.mpf(15))}

    # corner identity at 50 digits (h = 1e-20 balances the finite-difference
    # roundoff ~1e-50/h against truncation O(h^2); measured rel diff ~1e-29)
    Rtest = mp.mpf("6.409676")
    h = mp.mpf("1e-20")
    dds = (fbe(Rtest, Rtest + h) - fbe(Rtest, Rtest - h)) / (2 * h)
    ddr = dH(Rtest)

    return {
        "dps": 50,
        "R_bp": {"mpmath": str(R_bp_mp), "float64": R_bp_crossover()},
        "R_star": {"mpmath": str(R_star_mp), "float64": R_always_cheaper()},
        "s_c_table": {
            "R": [12.0, 13.0, 14.0, 15.0],
            "mpmath": [str(v) for v in [tbl[12.0], tbl[13.0], tbl[14.0], tbl[15.0]]],
            "float64": [crossover_s(12.0), crossover_s(13.0),
                        crossover_s(14.0), crossover_s(15.0)],
            "wikipedia_12_entry_suspected_error": "15.81 (should be 815.81)",
        },
        "corner_identity_50_digits": {
            "d_f_ds_at_corner": str(dds),
            "d_dvH_dR": str(ddr),
            "rel_diff": str(abs(dds - ddr) / abs(ddr)),
        },
    }


# --- RK4 trajectory validation of the full three-burn transfer -------------


def transfer_case_steps(e: float, base_per_orbit: int = 512) -> int:
    """Steps for a half-orbit leg, from the periapsis-resolution law of
    Experiment 002: steps/orbit ~ base * (1-e)^(-3/2)."""
    return max(64, int(np.ceil(0.5 * base_per_orbit / (1.0 - e) ** 1.5)))


def validate_bielliptic_rk4(
    r1: float, r2: float, r_b: float, mu: float = MU
) -> dict:
    """Propagate the complete three-burn bi-elliptic transfer:

      burn 1 at r1 -> half-orbit of ellipse 1 (periapsis r1 -> apoapsis r_b)
      burn 2 at r_b -> half-orbit of ellipse 2 (apoapsis r_b -> periapsis r2)
      burn 3 at r2 -> one full circular orbit at r2

    Every closed-form quantity (burns, leg times, arrival states) is checked
    against the verified RK4 propagator of Experiment 002; the analytic
    reference (kepler_solution of the same ellipse) is phase-shifted by a
    half-period for legs that start at apoapsis (the Experiment 004 lesson).
    """
    a1 = 0.5 * (r1 + r_b)
    a2 = 0.5 * (r2 + r_b)
    e1 = (r_b - r1) / (r_b + r1)
    e2 = (r_b - r2) / (r_b + r2)
    dv1, dv2, dv3 = bielliptic_burns(r1, r2, r_b, mu)
    t1 = np.pi * np.sqrt(a1**3 / mu)
    t2 = np.pi * np.sqrt(a2**3 / mu)
    v1c = np.sqrt(mu / r1)
    v2c = np.sqrt(mu / r2)

    def leg(t: np.ndarray, r0, v0):
        return propagate_rk4(r0, v0, mu, t)

    # --- Leg 1: from circular r1, burn 1 to ellipse 1 (periapsis r1).
    v_ell1_r1 = ellipse_speed(r1, a1, mu)
    assert abs(abs(v_ell1_r1 - v1c) - dv1) < 1e-9 * v1c
    r0 = np.array([r1, 0.0])
    v0 = np.array([0.0, v_ell1_r1])
    n1 = transfer_case_steps(e1)
    t_leg1 = np.linspace(0.0, t1, n1 + 1)
    st1 = leg(t_leg1, r0, v0)
    r_ap = np.hypot(st1[-1, 0], st1[-1, 1])
    v_ap = st1[-1, 2:]
    v_ap_mag = float(np.hypot(*v_ap))
    ana1 = kepler_solution(a1, e1, mu, t_leg1)  # starts at periapsis, as flown
    # arrival must be the apoapsis at exactly t1
    r_along1 = np.hypot(st1[:, 0], st1[:, 1])
    apside_ok1 = int(np.argmax(r_along1)) == n1

    # --- Burn 2 at r_b: boost from ellipse-1 apoapsis speed to ellipse-2
    #     apoapsis speed (both tangential at the apside).
    v_ell1_rb = ellipse_speed(r_b, a1, mu)
    v_ell2_rb = ellipse_speed(r_b, a2, mu)
    assert abs(abs(v_ell2_rb - v_ell1_rb) - dv2) < 1e-9 * v_ell2_rb
    dirn = v_ap / v_ap_mag  # tangential unit vector at apoapsis
    v_after2 = dirn * v_ell2_rb
    dv2_measured = float(np.hypot(*(v_after2 - v_ap)))
    assert abs(dv2_measured - dv2) / dv2 < 1e-6

    # --- Leg 2: from apoapsis r_b on ellipse 2 down to periapsis r2.
    r0b = np.array([-r_b, 0.0])
    v0b = np.array([0.0, -v_ell2_rb])  # prograde; r_b on the -x axis
    n2 = transfer_case_steps(e2)
    t_leg2 = np.linspace(0.0, t2, n2 + 1)
    st2 = leg(t_leg2, r0b, v0b)
    r_per = np.hypot(st2[-1, 0], st2[-1, 1])
    v_per = st2[-1, 2:]
    v_per_mag = float(np.hypot(*v_per))
    # analytic reference: same ellipse, but this leg starts at apoapsis, i.e.
    # the same ellipse one half-period later (Experiment 004 lesson).
    ana2 = kepler_solution(a2, e2, mu, t_leg2 + t2)
    r_along2 = np.hypot(st2[:, 0], st2[:, 1])
    apside_ok2 = int(np.argmin(r_along2)) == n2

    # --- Burn 3 at r2: circularize. The leg-2 arrival is the periapsis on
    #     the +x axis with prograde velocity (0, +v_per), so the target
    #     circular state there is (0, +v2c) - same direction, magnitude
    #     difference dv3 (Experiment 004's apsidal-burn pattern).
    v_ell2_r2 = ellipse_speed(r2, a2, mu)
    assert abs(abs(v_ell2_r2 - v2c) - dv3) < 1e-9 * v_ell2_r2
    v_after3 = np.array([0.0, v2c])
    dv3_measured = float(np.hypot(*(v_after3 - v_per)))
    assert abs(dv3_measured - dv3) / dv3 < 1e-6

    # --- One full circular orbit at r2.
    t_circ = np.linspace(0.0, 2.0 * np.pi * np.sqrt(r2**3 / mu), 1025)
    circ = propagate_rk4(np.array([r2, 0.0]), v_after3, mu, t_circ)
    r_circ = np.hypot(circ[:, 0], circ[:, 1])

    # energy / angular momentum drift on both legs (RK4, sampled)
    dE, dH = 0.0, 0.0
    for st in (st1, st2):
        for i in range(0, len(st), 4):
            els = orbital_elements(st[i, :2], st[i, 2:], mu)
            a_leg = a1 if st is st1 else a2
            h_leg = np.sqrt(mu * a_leg * (1.0 - (e1 if st is st1 else e2) ** 2))
            dE = max(dE, abs(els["energy"] - (-mu / (2.0 * a_leg))) / abs(
                -mu / (2.0 * a_leg)))
            dH = max(dH, abs(els["angular_momentum"] - h_leg) / abs(h_leg))

    return {
        "case": {"r1": r1, "r2": r2, "r_b": r_b, "R": r2 / r1, "s": r_b / r1,
                 "mu": mu, "e1": e1, "e2": e2, "t1": t1, "t2": t2,
                 "steps_leg1": n1, "steps_leg2": n2},
        "burn1": {"dv1": dv1, "v_circ_r1": v1c, "v_ellipse_r1": v_ell1_r1},
        "arrival_apoapsis_rk4": {
            "r_final": float(r_ap), "r_b": r_b,
            "rel_r_error": abs(r_ap - r_b) / r_b,
            "v_final_mag": v_ap_mag,
            "v_ellipse1_at_rb": v_ell1_rb,
            "rel_v_error": abs(v_ap_mag - v_ell1_rb) / v_ell1_rb,
            "apsis_at_final": bool(apside_ok1),
        },
        "arrival_apoapsis_analytic": {
            "rel_r_error": abs(float(np.hypot(ana1[-1, 0], ana1[-1, 1])) - r_b) / r_b,
            "rel_v_error": abs(float(np.hypot(ana1[-1, 2], ana1[-1, 3])) - v_ell1_rb)
            / v_ell1_rb,
        },
        "burn2": {"dv2": dv2, "v_ellipse1_at_rb": v_ell1_rb,
                  "v_ellipse2_at_rb": v_ell2_rb,
                  "dv2_measured": dv2_measured,
                  "rel_dv2_error": abs(dv2_measured - dv2) / dv2},
        "arrival_periapsis_rk4": {
            "r_final": float(r_per), "r2": r2,
            "rel_r_error": abs(r_per - r2) / r2,
            "v_final_mag": v_per_mag,
            "v_ellipse2_at_r2": v_ell2_r2,
            "rel_v_error": abs(v_per_mag - v_ell2_r2) / v_ell2_r2,
            "apsis_at_final": bool(apside_ok2),
        },
        "arrival_periapsis_analytic": {
            "rel_r_error": abs(float(np.hypot(ana2[-1, 0], ana2[-1, 1])) - r2) / r2,
            "rel_v_error": abs(float(np.hypot(ana2[-1, 2], ana2[-1, 3])) - v_ell2_r2)
            / v_ell2_r2,
        },
        "burn3": {"dv3": dv3, "v_ellipse2_at_r2": v_ell2_r2, "v_circ_r2": v2c,
                  "dv3_measured": dv3_measured,
                  "rel_dv3_error": abs(dv3_measured - dv3) / dv3},
        "post_burn_circular_orbit": {
            "radius_max_rel_variation": float(
                np.max(np.abs(r_circ - r2)) / r2),
            "speed_rel_error_vs_sqrt_mu_r2": float(
                max(abs(np.hypot(circ[i, 2], circ[i, 3]) - v2c) / v2c
                    for i in range(0, len(circ), 8))),
        },
        "max_rel_drift": {"energy": dE, "angular_momentum": dH},
    }


# --- Real-system anchors ----------------------------------------------------


def real_anchors() -> dict:
    """Engineering numbers with IAU 2015 B3 nominal constants (Earth) and
    IAU Sun + JPL mean orbits (heliocentric), reusing Experiment 004's
    constants and Hohmann machinery."""
    r_leo = R_EARTH_KM + LEO_ALT_KM
    r_geo = R_EARTH_KM + GEO_ALT_KM
    v_leo = np.sqrt(MU_EARTH_KM3S2 / r_leo)
    mu_e = MU_EARTH_KM3S2

    def case(name, r1, r2, r_b, mu):
        R = r2 / r1
        s = r_b / r1
        v1 = np.sqrt(mu / r1)
        dv_ho = hohmann_dv_total(r1, r2, mu)
        dv1, dv2, dv3 = bielliptic_burns(r1, r2, r_b, mu)
        dv_be = dv1 + dv2 + dv3
        dv_bp = S2_MINUS_1 * (1.0 + 1.0 / np.sqrt(R)) * v1
        t_ho = hohmann_transfer_time(r1, r2, mu)
        t_be = bielliptic_transfer_time(r1, r2, r_b, mu)
        return {
            "name": name, "r1_km": r1, "r2_km": r2, "r_b_km": r_b,
            "R": R, "s": s, "v1_km_s": v1,
            "dv_hohmann_km_s": dv_ho,
            "dv_bielliptic_km_s": dv_be,
            "dv1_km_s": dv1, "dv2_km_s": dv2, "dv3_km_s": dv3,
            "dv_bi_parabolic_km_s": dv_bp,
            "saving_vs_hohmann_km_s": dv_ho - dv_be,
            "saving_percent_of_hohmann": 100.0 * (dv_ho - dv_be) / dv_ho,
            "t_hohmann_days": t_ho / 86400.0,
            "t_bielliptic_days": t_be / 86400.0,
            "t_ratio": t_be / t_ho,
        }

    out = {}
    out["leo_geo"] = case(
        "LEO(200 km) -> GEO", r_leo, r_geo, 100.0 * r_leo, mu_e)
    out["wiki_14x_r0_6700"] = case(
        "Wikipedia worked example (r0 = 6700 km, R = 14, s = 40)",
        6700.0, 14.0 * 6700.0, 40.0 * 6700.0, mu_e)
    out["leo_50x_best_saving"] = case(
        "LEO -> 50 x LEO (max-saving ratio)", r_leo, 50.0 * r_leo,
        1e6 * r_leo, mu_e)  # near bi-parabolic: use finite s = 1e6
    out["geo_lunar"] = case(
        "GEO -> lunar distance", r_geo, 3.844e5, 3.0 * 3.844e5, mu_e)
    out["geo_15_58x"] = case(
        "GEO -> 15.58 x GEO (boundary)", r_geo, 15.58 * r_geo,
        30.0 * r_geo, mu_e)
    out["earth_mars"] = case(
        "Earth -> Mars (heliocentric, circular mean orbits)",
        AU_KM, MARS_A_AU * AU_KM, 3.0 * MARS_A_AU * AU_KM, MU_SUN_KM3S2)
    return out


# --- Figures ----------------------------------------------------------------


def make_figures(
    R_bp: float, R_star: float, curve: dict, rk4: list[dict]
) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (RESULTS_DIR / "figures").mkdir(parents=True, exist_ok=True)
    paths = []

    # 1. Cost curves vs R (the classical figure): Hohmann + bi-elliptic at
    #    fixed intermediate apoapsis ratios; crossover markers.
    R = np.logspace(np.log10(1.0001), np.log10(75.0), 600)
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    ax.plot(R, dv_over_v1(R), "k-", lw=2.2, label="Hohmann")
    for alpha in (16.0, 20.0, 40.0, 60.0, 100.0, 1000.0):
        ax.plot(R, f_high(R, alpha * np.ones_like(R)), lw=1.0,
                label=f"bi-elliptic, r_b/r1 = {alpha:g}")
    ax.plot(R, bi_parabolic(R), ":", lw=1.6, color="0.45",
            label="bi-parabolic limit (r_b -> inf)")
    ax.axvline(R_bp, color="C3", ls="--", lw=1.2)
    ax.text(R_bp, 0.415, f"R_bp = {R_bp:.4f}", rotation=90, fontsize=8,
            color="C3", va="bottom")
    ax.axvline(R_star, color="C1", ls="--", lw=1.2)
    ax.text(R_star, 0.415, f"R* = {R_star:.4f}", rotation=90, fontsize=8,
            color="C1", va="bottom")
    ax.set_xlabel("radius ratio R = r2/r1")
    ax.set_ylabel("total delta-v / v1")
    ax.set_title("Hohmann vs bi-elliptic transfer cost")
    ax.legend(fontsize=7, ncol=2, loc="lower right")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "cost_curves.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 2. Crossover map in the (R, s) plane: g < 0 (bi-elliptic wins) region.
    Rs = np.logspace(np.log10(1.001), np.log10(50.0), 260)
    ss = np.logspace(0.0, np.log10(50.0 * 1e5), 300)
    Rg, Sg = np.meshgrid(Rs, ss, indexing="ij")
    g = f_high(Rg, Sg) - dv_over_v1(Rg)
    # mask the invalid region s < R (the bi-elliptic needs s > R for the
    # high family; s = R is the Hohmann corner)
    valid = Sg > Rg * (1.0 + 1e-4)
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    c = ax.contourf(Rg, Sg, g, levels=np.linspace(-0.05, 0.05, 41),
                    cmap="RdBu_r", extend="both")
    ax.contour(Rg, Sg, g, levels=[0.0], colors="k", linewidths=1.5)
    ax.plot(Rs, Rs, "k--", lw=0.8, label="s = R (Hohmann corner)")
    ax.plot(curve["classical_table_R"], curve["classical_table_s_c"], "o",
            color="C2", ms=5, label="classical table points")
    ax.plot(curve["R_grid"], curve["s_c_grid"], "C2-", lw=1.8,
            label="crossover s_c(R)")
    ax.axvline(R_bp, color="C3", ls="--", lw=1.2)
    ax.axvline(R_star, color="C1", ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("radius ratio R = r2/r1")
    ax.set_ylabel("intermediate apoapsis ratio s = r_b/r1")
    ax.set_title("Bi-elliptic minus Hohmann delta-v; blue = bi-elliptic wins")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    fig.colorbar(c, ax=ax, label="(dv_biell - dv_Hohmann) / v1")
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "crossover_map.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 3. Shape of f(s) at representative R.
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 8.0), sharey=True)
    for axx, R in zip(axes.ravel(), (2.0, 10.0, 13.0, 14.0, 15.0, 20.0)):
        s = np.logspace(np.log10(R * 1.0001), np.log10(R * 1e7), 2000)
        f = f_high(R * np.ones_like(s), s)
        dvh = float(dv_over_v1(np.array([R]))[0])
        axx.semilogx(s, f, "C0-", lw=1.5)
        axx.axhline(dvh, color="C3", ls="--", lw=1.0,
                    label=f"dv_Hohmann = {dvh:.5f}")
        axx.axhline(float(bi_parabolic(np.array([R]))[0]), color="C1", ls=":",
                    lw=1.0, label="bi-parabolic limit")
        axx.set_title(f"R = {R:g}")
        axx.grid(True, which="both", alpha=0.3)
        axx.legend(fontsize=7)
    fig.suptitle("Bi-elliptic cost vs intermediate apoapsis s (v1 units)")
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "shape_per_R.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 4. Fuel saving vs R (Hohmann minus bi-parabolic) with peak + asymptote.
    R_sav = np.logspace(np.log10(1.0001), 4.0, 400)
    sav = dv_over_v1(R_sav) - bi_parabolic(R_sav)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.semilogx(R_sav, sav, "C0-", lw=1.6)
    ax.semilogx(R_sav, (2.0 - np.sqrt(2.0)) / np.sqrt(R_sav), "C1--", lw=1.1,
                label="asymptote (2 - sqrt(2))/sqrt(R)")
    ax.set_xlabel("radius ratio R = r2/r1")
    ax.set_ylabel("(dv_Hohmann - dv_bi-parabolic) / v1")
    ax.set_title("Asymptotic fuel saving of the bi-parabolic transfer")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "saving_curve.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 5. Time penalty: t_biell / t_H vs s at representative R.
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    for R in (2.0, 6.41, 14.0, 20.0):
        s = np.logspace(np.log10(R * 1.01), np.log10(R * 1e5), 400)
        t_ho = np.pi * np.sqrt((1.0 + R) ** 3 / 8.0)
        t_be = np.pi * (np.sqrt((1.0 + s) ** 3 / 8.0)
                        + np.sqrt((R + s) ** 3 / 8.0))
        ax.loglog(s, t_be / t_ho, lw=1.4, label=f"R = {R:g}")
    ax.set_xlabel("intermediate apoapsis ratio s = r_b/r1")
    ax.set_ylabel("transfer time ratio t_biell / t_Hohmann")
    ax.set_title("Bi-elliptic flight-time penalty")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "time_penalty.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 6. RK4 trajectory geometry of the (R = 14, s = 40) case.
    r1, r2, r_b = 1.0, 14.0, 40.0
    a1, a2 = 0.5 * (r1 + r_b), 0.5 * (r2 + r_b)
    e1, e2 = (r_b - r1) / (r_b + r1), (r_b - r2) / (r_b + r2)
    mu = 1.0
    t1 = np.pi * np.sqrt(a1**3 / mu)
    t2 = np.pi * np.sqrt(a2**3 / mu)
    # ellipses via kepler_solution over full periods (periapsis at t = 0);
    # ellipse 2 phase-shifted by half a period so both trace from their
    # common apoapsis r_b.
    th = np.linspace(0.0, 2.0 * np.pi, 800)
    fig, ax = plt.subplots(figsize=(7.0, 7.0))
    ax.plot(np.cos(th), np.sin(th), "C0-", lw=1.0, label="r1 = 1")
    ax.plot(14.0 * np.cos(th), 14.0 * np.sin(th), "C2-", lw=1.0,
            label="r2 = 14")
    e1_full = kepler_solution(a1, e1, mu, np.linspace(0.0, 2 * t1, 1200))
    ax.plot(e1_full[:, 0], e1_full[:, 1], "C1-", lw=1.8,
            label="transfer ellipse 1 (r1 -> r_b)")
    e2_full = kepler_solution(a2, e2, mu, np.linspace(0.0, 2 * t2, 1200) + t2)
    ax.plot(e2_full[:, 0], e2_full[:, 1], "C4-", lw=1.8,
            label="transfer ellipse 2 (r_b -> r2)")
    for (x, y, lbl, col) in (
        (1.0, 0.0, "burn 1", "r"),
        (-40.0, 0.0, "burn 2", "g"),
        (14.0, 0.0, "burn 3", "b"),
    ):
        ax.plot(x, y, "o", ms=7, color=col)
        ax.annotate(lbl, (x, y), textcoords="offset points", xytext=(8, 8),
                    fontsize=9, color=col)
    ax.plot(0, 0, "k*", ms=12)
    ax.set_aspect("equal")
    ax.set_xlabel("x [r1]")
    ax.set_ylabel("y [r1]")
    ax.set_title("Bi-elliptic transfer geometry: R = 14, r_b = 40 r1")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "trajectory_geometry.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    return paths


# --- Main -------------------------------------------------------------------


def main() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    R_bp = R_bp_crossover()
    R_star = R_always_cheaper()
    # independent cross-check: the R_star root must equal Exp 004's peak
    R_star_004 = peak_of_cost_curve()["R_star"]
    assert abs(R_star - R_star_004) / R_star_004 < 1e-4

    curve = crossover_curve(R_bp, R_star)
    regions = region_verification(R_bp, R_star)
    shapes = shape_diagnostics()
    sav_time = saving_and_time()
    low_in = low_family_and_inward(R_bp, R_star)
    corner = corner_identity_check(np.array([1.5, 2.0, 6.41, 11.94, 15.0]))
    hp = high_precision_verification()

    rk4_cases = [
        validate_bielliptic_rk4(1.0, 2.0, 5.0),
        validate_bielliptic_rk4(1.0, 6.41, 100.0),
        validate_bielliptic_rk4(1.0, 12.5, 100.0),
        validate_bielliptic_rk4(1.0, 14.0, 40.0),
        validate_bielliptic_rk4(1.0, 20.0, 20.5),
        validate_bielliptic_rk4(1.0, 20.0, 200.0),
    ]

    anchors = real_anchors()
    figures = make_figures(R_bp, R_star, curve, rk4_cases)

    print("=== Bi-elliptic vs Hohmann: theory constants ===")
    print(f"R_bp (bi-parabolic crossover)  = {R_bp:.10f}   (literature 11.94)")
    print(f"R*   (always-cheaper boundary)  = {R_star:.10f}  "
          f"(literature 15.58; Hohmann peak of Exp 004 = {R_star_004:.10f})")
    print(f"corner identity max rel diff    = {max(o['rel_diff'] for o in corner['checks']):.2e}")
    print(f"50-digit identity rel diff      = {hp['corner_identity_50_digits']['rel_diff']}")
    print("=== Classical crossover table (s_c, float64 vs mpmath 50 d) ===")
    for k, R in enumerate([12.0, 13.0, 14.0, 15.0]):
        print(f"  R = {R:5.2f}: s_c = {curve['classical_table_s_c'][k]:.6f} "
              f" (mpmath {hp['s_c_table']['mpmath'][k]})")
    print("=== Region verification (adversarial) ===")
    lo = regions["R_lo_region"]
    hi = regions["R_hi_region"]
    print(f"R < R_bp:   g > 0 on all {lo['n_R']}x{lo['n_s_per_R']} grid "
          f"(worst margin {lo['worst_margin']:.2e})")
    print(f"R > R_star: g < 0 on all {hi['n_R']}x{hi['n_s_per_R']} grid "
          f"(worst margin {hi['worst_margin']:.2e})")
    print(f"R_bp < R < R_star: crossings per R all == 1: "
          f"{regions['mid_crossing_counts']['all_exactly_one']}")
    print("=== Shape diagnostics ===")
    onset = shapes["hump_onset_R_approx"]
    print(f"hump onset R ~ {onset:.2f}")
    print("=== Saving and time ===")
    print(f"max saving dv_H - dv_bp = {sav_time['peak_saving_over_v1']:.6f} v1 "
          f"at R = {sav_time['peak_saving_at_R']:.1f}")
    for tc in sav_time["time_ratio_cases"]:
        print(f"  R={tc['R']:5.2f} s={tc['s']:5.1f}: t_biell/t_H = {tc['t_biell_over_t_hohmann']:8.1f}")
    print("=== RK4 validation (3-burn trajectory) ===")
    for case in rk4_cases:
        c = case["case"]
        a1r, a2r = case["arrival_apoapsis_rk4"], case["arrival_periapsis_rk4"]
        print(
            f"R={c['R']:.4g} s={c['s']:.4g} e1={c['e1']:.4f} e2={c['e2']:.4f}: "
            f"apoapsis rel r-err {a1r['rel_r_error']:.2e} "
            f"| periapsis rel r-err {a2r['rel_r_error']:.2e} "
            f"| burn2 err {case['burn2']['rel_dv2_error']:.2e} "
            f"| burn3 err {case['burn3']['rel_dv3_error']:.2e}"
        )
    print("=== Real anchors ===")
    for key in ("leo_geo", "wiki_14x_r0_6700", "leo_50x_best_saving",
                "geo_lunar", "geo_15_58x", "earth_mars"):
        a = anchors[key]
        print(
            f"{key:20s}: dv_H = {a['dv_hohmann_km_s']:7.4f}  "
            f"dv_BE(s={a['s']:.2g}) = {a['dv_bielliptic_km_s']:7.4f}  "
            f"saving = {a['saving_percent_of_hohmann']:+.2f}%  "
            f"t_ratio = {a['t_ratio']:7.1f}"
        )

    result = {
        "theory_constants": {
            "R_bp_bi_parabolic_crossover": R_bp,
            "R_star_always_cheaper_boundary": R_star,
            "R_star_exp004_peak_cross_check": R_star_004,
            "identity_d_ds_f_at_s_eq_R_equals_d_dR_dvH": corner,
        },
        "crossover_curve": curve,
        "region_verification": regions,
        "shape_diagnostics": shapes,
        "saving_and_time": sav_time,
        "low_family_and_inward": low_in,
        "high_precision_mpmath": hp,
        "rk4_transfer_validation": rk4_cases,
        "real_anchors": anchors,
        "constants": {
            "mu_earth_km3_s2": MU_EARTH_KM3S2,
            "r_earth_km": R_EARTH_KM,
            "leo_alt_km": LEO_ALT_KM,
            "geo_alt_km": GEO_ALT_KM,
            "au_km": AU_KM,
            "mu_sun_km3_s2": MU_SUN_KM3S2,
            "mars_a_au": MARS_A_AU,
        },
        "figures": [Path(p).name for p in figures],
    }
    path = save_json_result(
        RESULTS_DIR / "results.json",
        result,
        name="bielliptic_vs_hohmann",
        description=(
            "Bi-elliptic vs Hohmann transfer crossover: closed-form three-burn "
            "costs, the bi-parabolic limit, the 11.94 and 15.58 boundaries "
            "(the latter equal to the Hohmann cost maximum via the corner "
            "identity d/ds dv_biell|_{s=R} = d/dR dv_H), the crossover curve "
            "s_c(R), adversarial region verification, shape diagnostics, "
            "RK4 validation of the full three-burn trajectory, mpmath 50-digit "
            "cross-checks, and real-system anchors."
        ),
    )
    print(f"\nSaved results -> {path}")
    return result


if __name__ == "__main__":
    main()
