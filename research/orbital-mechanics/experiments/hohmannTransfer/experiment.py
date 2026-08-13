"""Hohmann transfer: the least-fuel two-impulse transfer between coplanar
circular orbits, and its delta-v budget.

Given circular orbits of radii r1 < r2 about a central body (mu), the Hohmann
transfer connects them with two tangential impulses:

  dv1 = v1 * (sqrt(2R/(1+R)) - 1)              at departure (r1), R = r2/r1
  dv2 = v2 * (1 - sqrt(2/(1+R)))               at arrival (r2), v = sqrt(mu/r)
  t_tr = pi * sqrt((r1 + r2)^3 / (8 mu))       half-period of the transfer
                                                ellipse, a_t = (r1+r2)/2,
                                                e_t = (r2-r1)/(r2+r1)

This experiment verifies the closed forms from first principles (vis-viva,
energy, angular momentum), sweeps the full radius-ratio space including the
pathological regimes R -> 1 and R -> infinity, numerically validates the
transfer by propagating it with the verified RK4 machinery of Experiment 002
(imported by explicit path - no scaffolding rebuilt), scans the two-impulse
transfer family for optimality of the Hohmann ellipse, and reproduces the
canonical real-world budgets (LEO->GEO ~3.93 km/s; Earth->Mars ~259 days,
v_infinity ~2.94/2.65 km/s; trans-Mars injection ~3.6 km/s) from IAU 2015
Resolution B3 nominal constants and JPL mean orbits.

Structure of the results:

  - the closed-form dvs and transfer elements, cross-checked via independent
    algebra (vis-viva, h-conservation r1 v_p = r2 v_a, energy -mu/(2 a_t)).
  - R-sweep: dv_total/v1 vs R with the small-R asymptote (R-1)/2, the
    R -> infinity asymptote sqrt(2)-1 (the escape burn!), the interior
    maximum near R* = 15.58, and outward/inward symmetry.
  - RK4 propagation of the full transfer (burn, coast, burn) for cases of
    growing eccentricity; arrival radius/speed/timing and the post-burn
    circular orbit are verified against the closed forms.
  - optimality scan: over the two-impulse family whose transfer ellipse has
    periapsis r_p <= r1 and apoapsis r_a >= r2 (burns at r1 and r2, no
    tangency assumption), the total dv is minimized at the Hohmann ellipse.
  - real-system anchors: LEO(200 km) -> GEO, Earth -> Mars, Earth -> Venus,
    trans-Mars injection from LEO, and Tsiolkovsky propellant fractions.

References: W. Hohmann, "Die Erreichbarkeit der Himmelskoerper",
Oldenbourg, 1925; R. R. Bate, D. D. Mueller, J. E. White, "Fundamentals of
Astrodynamics", Dover, 1971, Ch. 6; H. D. Curtis, "Orbital Mechanics for
Engineering Students", 4th ed., Elsevier, 2021, Ch. 6; D. A. Vallado,
"Fundamentals of Astrodynamics and Applications", 4th ed., Microcosm, 2013,
Ch. 6; IAU 2015 Resolution B3 nominal constants (Mamajek et al.,
arXiv:1510.07674); JPL SSD "Approximate Positions of the Planets"
(Standish & Williams 1992) for the Mars and Venus mean semi-major axes.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from lab_utils.results import save_json_result

# --- Reuse of the verified 002 machinery -----------------------------------
#
# The laboratory rule is "reuse existing code, never rebuild scaffolding".
# Experiment 002 (keplerOrbitValidation) owns the verified RK4 propagator,
# Kepler solver and elements extraction; it is loaded here from its explicit
# path (the same importlib pattern the test suites use), so this experiment
# stays self-contained in the sense of one runnable entry point while the
# propagation code has a single source of truth. Experiments are closed, so
# this dependency is stable.

_kov_path = (
    Path(__file__).resolve().parents[1] / "keplerOrbitValidation" / "experiment.py"
)
_kov_spec = importlib.util.spec_from_file_location(
    "kepler_orbit_validation_exp002", _kov_path
)
assert _kov_spec is not None and _kov_spec.loader is not None
kov = importlib.util.module_from_spec(_kov_spec)
_kov_spec.loader.exec_module(kov)

propagate_rk4 = kov.propagate_rk4
kepler_solution = kov.kepler_solution
orbital_elements = kov.orbital_elements

# --- Physical constants (verified against primary sources, 2026-08-13) -----

# IAU 2015 Resolution B3 nominal conversion constants (ia u.org resolution
# text; Mamajek et al., arXiv:1510.07674): exact by definition.
MU_EARTH_KM3S2 = 3.986004e5  # 1 (GM)^N_E = 3.986004e14 m^3/s^2
R_EARTH_KM = 6.3781e3  # 1 R^N_eE = 6.3781e6 m (nominal equatorial radius)
# Canonical geostationary altitude over the equatorial radius (public
# standard, e.g., Curtis Ch. 6; again used here as the textbook value).
GEO_ALT_KM = 3.5786e4
LEO_ALT_KM = 2.0e2
AU_KM = kov.AU_KM  # 149597870.7 km, IAU 2012 Resolution B2 (exact)
MU_SUN_KM3S2 = kov.MU_SUN_KM3_S2  # (GM)^N_Sun IAU 2015 B3 nominal
# Mean heliocentric semi-major axes, JPL SSD "Approximate Positions of the
# Planets" (Standish & Williams 1992), mean ecliptic J2000, EM barycenter
# (same source family Experiment 002 used for Earth's eccentricity).
MARS_A_AU = 1.523679
VENUS_A_AU = 0.723332
G0 = 9.80665  # standard gravity, m/s^2 (SI definition)

MU = 1.0  # canonical gravitational parameter for the core studies

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

# --- Closed-form Hohmann transfer -------------------------------------------


def hohmann_dv1(r1: float, r2: float, mu: float = MU) -> float:
    """First impulse of the transfer, executed at radius r1.

    For r2 > r1 (outward): a prograde raise from the circular speed v1 to the
    transfer-ellipse periapsis speed. For r2 < r1 (inward): a retrograde drop
    from v1 to the transfer-ellipse apoapsis speed. The textbook formulas are
    written for r2 > r1; the inward branch follows from the same transfer
    ellipse (periapsis r_min, apoapsis r_max) evaluated at r1.
    """
    R = r2 / r1
    if R >= 1.0:
        return np.sqrt(mu / r1) * (np.sqrt(2.0 * R / (1.0 + R)) - 1.0)
    rho = r1 / r2
    return np.sqrt(mu / r1) * (1.0 - np.sqrt(2.0 / (1.0 + rho)))


def hohmann_dv2(r1: float, r2: float, mu: float = MU) -> float:
    """Second impulse of the transfer, executed at radius r2.

    Outward: a prograde raise from the transfer-ellipse apoapsis speed to v2.
    Inward: a prograde raise from v2 to the transfer-ellipse periapsis speed.
    """
    R = r2 / r1
    if R >= 1.0:
        return np.sqrt(mu / r2) * (1.0 - np.sqrt(2.0 / (1.0 + R)))
    rho = r1 / r2
    return np.sqrt(mu / r2) * (np.sqrt(2.0 * rho / (1.0 + rho)) - 1.0)


def hohmann_dv_total(r1: float, r2: float, mu: float = MU) -> float:
    return hohmann_dv1(r1, r2, mu) + hohmann_dv2(r1, r2, mu)


def hohmann_transfer_time(r1: float, r2: float, mu: float = MU) -> float:
    """Half the period of the transfer ellipse (periapsis to apoapsis)."""
    return np.pi * np.sqrt((r1 + r2) ** 3 / (8.0 * mu))


def transfer_elements(r1: float, r2: float, mu: float = MU) -> dict:
    """Transfer ellipse geometry from its two apsidal radii r1 and r2.

    The transfer ellipse has periapsis min(r1, r2) and apoapsis max(r1, r2);
    the state speeds at the two radii follow from vis-viva. Direction-agnostic
    (the ellipse is identical for both transfer directions).
    """
    a = 0.5 * (r1 + r2)
    e = abs(r2 - r1) / (r2 + r1)
    p = 2.0 * r1 * r2 / (r1 + r2)  # semi-latus rectum = a(1-e^2)
    h = np.sqrt(mu * p)  # specific angular momentum
    v_at_r1 = np.sqrt(mu * (2.0 / r1 - 1.0 / a))  # vis-viva at r1
    v_at_r2 = np.sqrt(mu * (2.0 / r2 - 1.0 / a))  # vis-viva at r2
    eps = -mu / (2.0 * a)  # specific energy of the transfer orbit
    return {
        "a": a,
        "e": e,
        "p": p,
        "h": h,
        "v_at_r1": v_at_r1,
        "v_at_r2": v_at_r2,
        "energy": eps,
        "v_periapsis": np.sqrt(mu * (2.0 / min(r1, r2) - 1.0 / a)),
        "v_apoapsis": np.sqrt(mu * (2.0 / max(r1, r2) - 1.0 / a)),
    }


def hohmann_split(r1: float, r2: float, mu: float = MU) -> dict:
    """Full closed-form budget for the transfer (r1 < r2 assumed)."""
    R = r2 / r1
    v1 = np.sqrt(mu / r1)
    v2 = np.sqrt(mu / r2)
    dv1 = hohmann_dv1(r1, r2, mu)
    dv2 = hohmann_dv2(r1, r2, mu)
    els = transfer_elements(r1, r2, mu)
    return {
        "r1": r1,
        "r2": r2,
        "R": R,
        "v1": v1,
        "v2": v2,
        "dv1": dv1,
        "dv2": dv2,
        "dv_total": dv1 + dv2,
        "dv_total_over_v1": (dv1 + dv2) / v1,
        "transfer_time": hohmann_transfer_time(r1, r2, mu),
        "transfer": els,
        "v_at_r1": els["v_at_r1"],
        "v_at_r2": els["v_at_r2"],
        "dv1_magnitude_consistency": abs(els["v_at_r1"] - v1),
        "dv2_magnitude_consistency": abs(els["v_at_r2"] - v2),
    }


def hohmann_dv_stable(r1: float, r2: float, mu: float = MU) -> tuple[float, float]:
    """Digit-safe forms of dv1, dv2 for R -> 1 (no cancellation).

    sqrt(2R/(1+R)) - 1 = (R-1)/((1+R)(1 + sqrt(2R/(1+R))))
      1 - sqrt(2/(1+R)) = (R-1)/((1+R)(1 + sqrt(2/(1+R))))
    """
    R = r2 / r1
    v1 = np.sqrt(mu / r1)
    v2 = np.sqrt(mu / r2)
    d1 = v1 * (R - 1.0) / ((1.0 + R) * (1.0 + np.sqrt(2.0 * R / (1.0 + R))))
    d2 = v2 * (R - 1.0) / ((1.0 + R) * (1.0 + np.sqrt(2.0 / (1.0 + R))))
    return d1, d2


# --- R-sweep: circle-to-circle cost across the full ratio space -------------


def dv_over_v1(R: np.ndarray) -> np.ndarray:
    """Total dv in units of the inner circular speed, as a pure function of
    the radius ratio R = r2/r1 (r1 < r2)."""
    return (
        (np.sqrt(2.0 * R / (1.0 + R)) - 1.0)
        + (1.0 / np.sqrt(R)) * (1.0 - np.sqrt(2.0 / (1.0 + R)))
    )


def peak_of_cost_curve() -> dict:
    """Locate the interior maximum of dv_total/v1 vs R (dense grid + parabolic
    refinement). The literature value is R* ~ 15.58, dv/v1 ~ 0.5363."""
    Rs = np.logspace(0.0, 3.0, 2000001)
    cost = dv_over_v1(Rs)
    i = int(np.argmax(cost))
    # Parabolic refinement on three neighbours of the grid maximum:
    # fit y = a x^2 + b x + c in x = ln R, optimum at x = -b/(2a).
    x = np.log(Rs[i - 1 : i + 2])
    y = cost[i - 1 : i + 2]
    a_fit = ((y[2] - y[1]) / ((x[2] - x[1]) * (x[2] - x[0]))
             - (y[1] - y[0]) / ((x[1] - x[0]) * (x[2] - x[0])))
    b_fit = (y[1] - y[0]) / (x[1] - x[0]) - a_fit * (x[1] + x[0])
    x_opt = -b_fit / (2.0 * a_fit)
    R_star = float(np.exp(x_opt))
    val = float(a_fit * x_opt**2 + b_fit * x_opt
                + (y[0] - a_fit * x[0] ** 2 - b_fit * x[0]))
    return {"R_star": R_star, "dv_over_v1_peak": val,
            "grid_argmax_R": float(Rs[i]), "grid_dv_over_v1": float(cost[i])}


def ratio_sweep() -> dict:
    """Sweep R = r2/r1 across the full range, with asymptote checks."""
    R_grid = np.array([
        1.000001, 1.00001, 1.0001, 1.001, 1.01, 1.05, 1.1, 1.2, 1.5, 2.0,
        3.0, 4.0, 6.409676, 8.0, 10.0, 11.94, 15.581885187539546, 20.0,
        30.0, 50.0, 1e2, 1e3, 1e6, 1e9, 1e12,
    ])
    v1 = 1.0
    cells = []
    for R in R_grid:
        r1 = 1.0
        r2 = R * r1
        dv1_t, dv2_t = hohmann_dv_stable(r1, r2)
        dv1 = hohmann_dv1(r1, r2)
        dv2 = hohmann_dv2(r1, r2)
        cells.append({
            "R": float(R),
            "dv1_over_v1": dv1 / v1,
            "dv2_over_v1": dv2 / v1,
            "dv_total_over_v1": (dv1 + dv2) / v1,
            "dv_stable_over_v1": (dv1_t + dv2_t) / v1,
            "rel_diff_stable_vs_textbook": abs(
                (dv1_t + dv2_t) - (dv1 + dv2)) / (dv1_t + dv2_t),
            "transfer_time": hohmann_transfer_time(r1, r2),
            "e_transfer": transfer_elements(r1, r2)["e"],
        })
    asym_small = 0.5 * (R_grid - 1.0)
    cells_arr = np.array([c["dv_total_over_v1"] for c in cells])
    R_star_far = 1e12
    asym_far = np.sqrt(2.0) - 1.0
    peak = peak_of_cost_curve()
    return {
        "mu": MU,
        "R_grid": R_grid.tolist(),
        "cells": cells,
        "asymptote_small_R": "dv/v1 ~ (R-1)/2",
        "small_R_ratio_at_1p0001": float(
            cells_arr[2] / asym_small[2]),
        "small_R_ratio_at_1p001": float(cells_arr[3] / asym_small[3]),
        "asymptote_large_R": float(asym_far),
        "large_R_dv_over_v1_at_1e12": float(
            dv_over_v1(np.array([R_star_far]))[0]),
        "large_R_rel_gap": float(abs(dv_over_v1(np.array([R_star_far]))[0]
                                    - asym_far) / asym_far),
        "peak": peak,
        "outward_inward_symmetry_max_rel_diff": _symmetry_check(cells),
    }


def _symmetry_check(cells: list[dict]) -> float:
    """dv_total for transferring from radius r1 to R*r1 must equal the cost of
    the inward transfer at the same ratio (pure function of R)."""
    worst = 0.0
    for c in cells:
        R = c["R"]
        r1 = 1.0
        out = hohmann_dv_total(r1, R * r1)
        inn = hohmann_dv_total(R * r1, r1)  # inward at the same ratio
        worst = max(worst, abs(out - inn) / out)
    return float(worst)


# --- RK4 trajectory validation of the transfer -------------------------------


def transfer_case_steps(e: float, base_per_orbit: int = 512) -> int:
    """Steps for the half-orbit transfer, from the periapsis-resolution law
    of Experiment 002: steps/orbit ~ base * (1-e)^(-3/2)."""
    return max(64, int(np.ceil(0.5 * base_per_orbit / (1.0 - e) ** 1.5)))


def validate_transfer_rk4(
    r1: float, r2: float, mu: float = MU, steps: int | None = None
) -> dict:
    """Propagate the complete Hohmann transfer: burn1, half-orbit coast, burn2,
    and one orbit of the resulting circular orbit. Validates every closed-form
    quantity against the verified RK4 propagator of Experiment 002."""
    split = hohmann_split(r1, r2, mu)
    R = split["R"]
    dv1 = split["dv1"]
    dv2 = split["dv2"]
    t_tr = split["transfer_time"]
    v_r1 = split["v_at_r1"]
    v_r2 = split["v_at_r2"]
    v2 = split["v2"]
    a_t = split["transfer"]["a"]
    e_t = split["transfer"]["e"]
    h_t = split["transfer"]["h"]

    # --- Burn 1 at r1: to the transfer-ellipse speed at r1 (prograde for
    #     r2 > r1, retrograde for r2 < r1; the magnitude |dv1| is the burn).
    v1_circ = np.sqrt(mu / r1)
    assert abs(abs(v_r1 - v1_circ) - dv1) < 1e-9 * v_r1  # closed-form check
    r0 = np.array([r1, 0.0])
    v0 = np.array([0.0, v_r1])

    if steps is None:
        steps = transfer_case_steps(e_t)
    t = np.linspace(0.0, t_tr, steps + 1)
    states = propagate_rk4(r0, v0, mu, t)

    # --- Arrival at r2: apoapsis for r2 > r1, periapsis for r2 < r1. The
    # final propagated point must be the extremal radius (the apside).
    r_fin = np.hypot(states[-1, 0], states[-1, 1])
    v_fin = states[-1, 2:]
    v_fin_r = float(v_fin[0])  # radial component (x-axis at the apside)
    v_fin_mag = float(np.hypot(*v_fin))
    r_along = np.hypot(states[:, 0], states[:, 1])
    r_max_idx = int(np.argmax(r_along))
    r_min_idx = int(np.argmin(r_along))
    apsis_at_final = r_max_idx == len(states) - 1 or r_min_idx == len(states) - 1

    # Analytic reference for the same flight: the closed-form transfer orbit.
    # kepler_solution starts at periapsis; the outward flight (r2 > r1) does
    # too, but the inward flight (r2 < r1) starts at apoapsis, i.e. the same
    # ellipse one half-period later. Shift the reference times accordingly so
    # the analytic endpoint is the arrival apside (r2, v_r2), not the
    # departure one (avoids a misleading 100 %/50 % mismatch for r2 < r1).
    t_analytic = t + (t_tr if R < 1.0 else 0.0)
    ana = kepler_solution(a_t, e_t, mu, t_analytic)
    ana_r_fin = np.hypot(ana[-1, 0], ana[-1, 1])
    ana_v_fin = np.hypot(ana[-1, 2], ana[-1, 3])

    # Energy and angular momentum drift over the half orbit (RK4).
    dE, dH = 0.0, 0.0
    eps0 = -mu / (2.0 * a_t)
    for i in range(0, len(states), 4):
        els = orbital_elements(states[i, :2], states[i, 2:], mu)
        dE = max(dE, abs(els["energy"] - eps0) / abs(eps0))
        dH = max(dH, abs(els["angular_momentum"] - h_t) / abs(h_t))

    # --- Burn 2: slow to v2 at apoapsis; then one full circular orbit.
    v_after = np.array([0.0, -v2])  # (-r2, 0), prograde = -y direction
    dv2_measured = float(np.hypot(*(v_after - v_fin)))
    t_circ = np.linspace(0.0, 2.0 * np.pi * np.sqrt(r2**3 / mu), 1025)
    circ = propagate_rk4(np.array([-r2, 0.0]), v_after, mu, t_circ)
    r_circ = np.hypot(circ[:, 0], circ[:, 1])

    result = {
        "case": {"r1": r1, "r2": r2, "mu": mu, "R": R,
                 "e_transfer": e_t, "a_transfer": a_t,
                 "transfer_time": t_tr, "steps_half_orbit": steps},
        "burn1": {"dv1": dv1, "v_after_burn": v_r1, "v1_circ": v1_circ},
        "arrival_rk4": {
            "r_final": float(r_fin),
            "r2": r2,
            "rel_r_error": abs(r_fin - r2) / r2,
            "v_final_mag": v_fin_mag,
            "v_r2": v_r2,
            "rel_v_error": abs(v_fin_mag - v_r2) / v_r2,
            "v_radial": v_fin_r,
            "v_radial_over_v": abs(v_fin_r) / v_fin_mag,
            "apsis_at_final": apsis_at_final,
        },
        "arrival_analytic": {
            "r_final": float(ana_r_fin),
            "v_final": float(ana_v_fin),
            "rel_r_error": abs(ana_r_fin - r2) / r2,
            "rel_v_error": abs(ana_v_fin - v_r2) / v_r2,
        },
        "max_rel_drift": {"energy": dE, "angular_momentum": dH},
        "burn2": {
            "dv2": dv2,
            "dv2_measured_against_v2": dv2_measured,
            "rel_dv2_error": abs(dv2_measured - dv2) / dv2,
            "v2_circ": v2,
        },
        "post_burn_circular_orbit": {
            "radius_max_rel_variation": float(
                np.max(np.abs(r_circ - r2)) / r2),
            "speed_rel_error_vs_sqrt_mu_r2": float(
                max(abs(np.hypot(circ[i, 2], circ[i, 3]) - v2) / v2
                    for i in range(0, len(circ), 8))),
        },
    }
    return result


# --- Optimality: two-impulse family scan ------------------------------------


def two_impulse_cost(r1: float, r2: float, mu: float, r_p: float, r_a: float) -> float:
    """Total dv of a coplanar two-impulse transfer whose transfer ellipse has
    periapsis r_p <= r1 and apoapsis r_a >= r2, with both impulses executed at
    r1 and r2 on the outward leg (no tangency assumed). The cost is the sum of
    the vector velocity mismatches against the circular states."""
    e = (r_a - r_p) / (r_a + r_p)
    p = 2.0 * r_p * r_a / (r_p + r_a)
    h = np.sqrt(mu * p)
    nu1 = np.arccos(np.clip((p / r1 - 1.0) / e, -1.0, 1.0))
    nu2 = np.arccos(np.clip((p / r2 - 1.0) / e, -1.0, 1.0))
    v_r1 = (mu / h) * e * np.sin(nu1)
    v_t1 = h / r1
    v_r2 = (mu / h) * e * np.sin(nu2)
    v_t2 = h / r2
    v1c = np.sqrt(mu / r1)
    v2c = np.sqrt(mu / r2)
    dv1 = np.hypot(v_t1 - v1c, v_r1)
    dv2 = np.hypot(v2c - v_t2, v_r2)
    return dv1 + dv2


def optimality_scan(r1: float, r2: float, mu: float = MU) -> dict:
    """Grid over (r_p, r_a): the two-impulse cost family, minimised by the
    Hohmann ellipse (r_p = r1, r_a = r2)."""
    n_p, n_a = 121, 131
    r_p_grid = np.linspace(0.02 * r1, r1, n_p)
    r_a_grid = np.logspace(np.log10(r2), np.log10(200.0 * r2), n_a)
    Rp, Ra = np.meshgrid(r_p_grid, r_a_grid, indexing="ij")
    cost = np.empty_like(Rp)
    for i in range(n_p):
        for j in range(n_a):
            cost[i, j] = two_impulse_cost(r1, r2, mu, Rp[i, j], Ra[i, j])
    imin = np.unravel_index(int(np.argmin(cost)), cost.shape)
    dv_hohmann = hohmann_dv_total(r1, r2, mu)
    # 1D family: tangential departure (r_p = r1), r_a >= r2.
    r_a_curve = np.logspace(np.log10(r2), np.log10(200.0 * r2), 400)
    cost_curve = np.array(
        [two_impulse_cost(r1, r2, mu, r1, ra) for ra in r_a_curve])
    # 1D family: r_a = r2, r_p <= r1 (advantages of lower periapsis = none).
    r_p_curve = np.linspace(0.02 * r1, r1, 400)
    cost_rp = np.array(
        [two_impulse_cost(r1, r2, mu, rp, r2) for rp in r_p_curve])
    cost_edge = cost_curve[0]
    return {
        "r1": r1,
        "r2": r2,
        "R": r2 / r1,
        "grid": {
            "n_rp": n_p,
            "n_ra": n_a,
            "cost_min": float(np.min(cost)),
            "cost_min_at_hohmann_corner": float(
                two_impulse_cost(r1, r2, mu, r1, r2)),
            "dv_hohmann": dv_hohmann,
            "rel_gap_min_vs_hohmann": abs(float(np.min(cost)) - dv_hohmann)
            / dv_hohmann,
            "argmin_rp_over_r1": float(Rp[imin] / r1),
            "argmin_ra_over_r2": float(Ra[imin] / r2),
            "argmin_is_hohmann_corner": bool(
                abs(Rp[imin] / r1 - 1.0) < 1e-6
                and abs(Ra[imin] / r2 - 1.0) < 1e-6),
        },
        "curve_tangent_departure": {
            "min_rel_excess_over_hohmann": float(
                (np.min(cost_curve) - dv_hohmann) / dv_hohmann),
            "min_at_r_a_over_r2": float(
                r_a_curve[int(np.argmin(cost_curve))] / r2),
            "monotone_increasing_after_min": bool(
                np.all(np.diff(cost_curve[int(np.argmin(cost_curve)):]) > -1e-12)),
        },
        "curve_apoapsis_fixed": {
            "min_rel_excess_over_hohmann": float(
                (np.min(cost_rp) - dv_hohmann) / dv_hohmann),
            "min_at_r_p_over_r1": float(
                r_p_curve[int(np.argmin(cost_rp))] / r1),
            "monotone_increasing_toward_hohmann": bool(
                np.all(np.diff(cost_rp) < 1e-12 * (1.0 + np.abs(cost_rp[:-1])))
                and cost_rp[-1] == np.min(cost_rp)),
        },
        "sample_grid": {
            "r_p_axis_over_r1": (r_p_grid / r1)[::12].tolist(),
            "r_a_axis_over_r2": (r_a_grid / r2)[::16].tolist(),
            "cost_row_rp_eq_r1": (cost[-1, ::16]).tolist(),
        },
        "curve_selected": {
            "r_a_over_r2": (r_a_curve / r2)[::20].tolist(),
            "cost_minus_dv_hohmann_over_v1": (
                (cost_curve - dv_hohmann) / np.sqrt(mu / r1))[::20].tolist(),
            "r_p_over_r1": (r_p_curve / r1)[::20].tolist(),
            "cost_rp_minus_dv_hohmann_over_v1": (
                (cost_rp - dv_hohmann) / np.sqrt(mu / r1))[::20].tolist(),
        },
    }


# --- Real-system anchors ----------------------------------------------------


def hyperbolic_excess_requires(vinf: float, v_circ: float) -> float:
    """Delta-v to reach hyperbolic excess vinf from a circular orbit: the
    departure speed is sqrt(2 v_circ^2 + vinf^2) (Oberth / vis-viva)."""
    return np.sqrt(2.0 * v_circ**2 + vinf**2) - v_circ


def tsiolkovsky_fraction(dv: float, isp: float) -> float:
    """Propellant mass fraction 1 - exp(-dv / (Isp * g0)) of a single stage.

    ``dv`` is in km/s (the experiment's unit); the exhaust velocity Isp*g0 is
    in m/s, so the exponent is dv*1000/(Isp*g0).
    """
    return 1.0 - np.exp(-dv * 1000.0 / (isp * G0))


def real_cases() -> dict:
    """Reality anchors built from IAU nominal constants and JPL mean orbits:
    LEO -> GEO, Earth -> Mars (heliocentric), Earth -> Venus, the
    trans-Mars injection from LEO, and Tsiolkovsky propellant budgets."""
    # --- LEO 200 km -> GEO (Earth) ---
    r_leo = R_EARTH_KM + LEO_ALT_KM
    r_geo = R_EARTH_KM + GEO_ALT_KM
    leo_geo = hohmann_split(r_leo, r_geo, MU_EARTH_KM3S2)
    leo_geo["transfer_time_hours"] = leo_geo["transfer_time"] / 3600.0
    leo_geo["fuel_frac_isp300"] = tsiolkovsky_fraction(
        leo_geo["dv_total"], 300.0)
    leo_geo["fuel_frac_isp450"] = tsiolkovsky_fraction(
        leo_geo["dv_total"], 450.0)
    for k in ("dv1", "dv2", "dv_total", "v1", "v2"):
        leo_geo[k] = leo_geo[k]  # already km/s

    # --- Earth -> Mars (heliocentric, circular mean orbits) ---
    a_earth = AU_KM
    a_mars = MARS_A_AU * AU_KM
    e_m = hohmann_split(a_earth, a_mars, MU_SUN_KM3S2)
    e_m["transfer_time_days"] = e_m["transfer_time"] / 86400.0
    # v_infinity at departure = dv1 (both measured against the same circular
    # speed); arrival v_infinity = dv2 by construction of the Hohmann orbit.
    e_m["vinf_departure_km_s"] = e_m["dv1"]
    e_m["vinf_arrival_km_s"] = e_m["dv2"]
    # From 200 km LEO: escape speed, then to v_infinity.
    v_leo = np.sqrt(MU_EARTH_KM3S2 / r_leo)
    e_m["tmi_from_200km_leo_km_s"] = hyperbolic_excess_requires(
        e_m["vinf_departure_km_s"], v_leo)
    e_m["fuel_frac_leo_to_tmi_isp300"] = tsiolkovsky_fraction(
        e_m["tmi_from_200km_leo_km_s"], 300.0)
    # Mars orbit insertion at 300 km: capture from v_infinity.
    r_mars_300 = 3.3895e3 + 3.0e2
    v_mars_circ = np.sqrt(4.2828e4 / r_mars_300)
    e_m["mars_orbit_insertion_from_vinf_km_s"] = (
        hyperbolic_excess_requires(e_m["vinf_arrival_km_s"], v_mars_circ))
    e_m["total_leo_to_mars_orbit_km_s"] = (
        e_m["tmi_from_200km_leo_km_s"]
        + e_m["mars_orbit_insertion_from_vinf_km_s"])

    # --- Earth -> Venus (inward heliocentric: r2 < r1) ---
    a_venus = VENUS_A_AU * AU_KM
    e_v = hohmann_split(a_earth, a_venus, MU_SUN_KM3S2)  # inward
    e_v["transfer_time_days"] = e_v["transfer_time"] / 86400.0
    e_v["vinf_departure_km_s"] = e_v["dv1"]
    e_v["vinf_arrival_km_s"] = e_v["dv2"]
    e_v["tmi_from_200km_leo_km_s"] = hyperbolic_excess_requires(
        e_v["vinf_departure_km_s"], v_leo)

    return {"leo_geo": leo_geo, "earth_mars": e_m, "earth_venus": e_v}


# --- Figures ----------------------------------------------------------------


def make_figures(sweep: dict, rk4: list[dict], opt: list[dict]) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (RESULTS_DIR / "figures").mkdir(parents=True, exist_ok=True)
    paths = []

    # 1. Cost curve vs radius ratio with asymptotes and peak.
    Rs = np.logspace(np.log10(1.0001), 3.0, 400)
    cost = dv_over_v1(Rs)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.semilogx(Rs, cost, "C0-", lw=1.6, label="dv_total / v1")
    Rs_s = np.linspace(1.0, 1.3, 60)
    ax.semilogx(Rs_s, 0.5 * (Rs_s - 1.0), "C1--", lw=1.0,
                label="asymptote (R-1)/2 as R -> 1")
    ax.axhline(np.sqrt(2.0) - 1.0, color="C2", ls=":", lw=1.2,
               label="sqrt(2)-1 (escape burn, R -> inf)")
    pk = sweep["peak"]
    ax.plot(pk["R_star"], pk["dv_over_v1_peak"], "ks", ms=6,
            label=f"peak R* = {pk['R_star']:.2f}, dv/v1 = {pk['dv_over_v1_peak']:.4f}")
    ax.set_xlabel("radius ratio R = r2/r1")
    ax.set_ylabel("total delta-v / v1 (inner circular speed)")
    ax.set_title("Hohmann transfer cost across the radius-ratio space")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "dv_vs_radius_ratio.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 2. LEO -> GEO transfer geometry (real units, analytic trajectory).
    r_leo = R_EARTH_KM + LEO_ALT_KM
    r_geo = R_EARTH_KM + GEO_ALT_KM
    a_t = 0.5 * (r_leo + r_geo)
    e_t = (r_geo - r_leo) / (r_geo + r_leo)
    mu_e = MU_EARTH_KM3S2
    t_tr = hohmann_transfer_time(r_leo, r_geo, mu_e)
    t = np.linspace(0.0, t_tr, 900)
    tr = kepler_solution(a_t, e_t, mu_e, t)
    th = np.linspace(0.0, 2.0 * np.pi, 600)
    fig, ax = plt.subplots(figsize=(6.8, 6.8))
    ax.plot(r_leo * np.cos(th), r_leo * np.sin(th), "C0-", lw=1.0,
            label="LEO 200 km")
    ax.plot(r_geo * np.cos(th), r_geo * np.sin(th), "C2-", lw=1.0,
            label="GEO (35 786 km)")
    ax.plot(tr[:, 0], tr[:, 1], "C1-", lw=1.8,
            label="Hohmann transfer ellipse")
    ax.plot(tr[0, 0], tr[0, 1], "r^", ms=9, label="burn 1 (dv1 = 2.45 km/s)")
    ax.plot(tr[-1, 0], tr[-1, 1], "rv", ms=9, label="burn 2 (dv2 = 1.48 km/s)")
    ax.plot(0, 0, "k*", ms=12, label="Earth (nominal R = 6378.1 km)")
    ax.set_aspect("equal")
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_title("LEO -> GEO Hohmann transfer (IAU nominal constants)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "leo_geo_transfer.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 3. Earth -> Mars heliocentric geometry.
    a_t = 0.5 * AU_KM * (1.0 + MARS_A_AU)
    e_t = (MARS_A_AU - 1.0) / (MARS_A_AU + 1.0)
    t_tr = hohmann_transfer_time(AU_KM, MARS_A_AU * AU_KM, MU_SUN_KM3S2)
    tr = kepler_solution(a_t, e_t, MU_SUN_KM3S2,
                         np.linspace(0.0, t_tr, 700))
    th = np.linspace(0.0, 2.0 * np.pi, 600)
    fig, ax = plt.subplots(figsize=(6.8, 6.8))
    ax.plot(AU_KM * np.cos(th), AU_KM * np.sin(th), "C0-", lw=1.0,
            label="Earth orbit (1 au)")
    ax.plot(MARS_A_AU * AU_KM * np.cos(th), MARS_A_AU * AU_KM * np.sin(th),
            "C2-", lw=1.0, label="Mars orbit (1.523679 au)")
    ax.plot(tr[:, 0], tr[:, 1], "C1-", lw=1.8,
            label="Hohmann transfer (259 d)")
    ax.plot(AU_KM, 0.0, "C0o", ms=5)
    ax.plot(-MARS_A_AU * AU_KM, 0.0, "C2o", ms=5)
    ax.plot(0, 0, "k*", ms=12, label="Sun")
    ax.set_aspect("equal")
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_title("Earth -> Mars Hohmann transfer (heliocentric)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "earth_mars_transfer.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 4. Optimality: cost excess over Hohmann for the two-impulse family.
    r1 = opt[0]["r1"]
    r2 = opt[0]["r2"]
    n_p, n_a = 121, 131
    r_p_grid = np.linspace(0.02 * r1, r1, n_p)
    r_a_grid = np.logspace(np.log10(r2), np.log10(200.0 * r2), n_a)
    cost = np.empty((n_p, n_a))
    for i in range(n_p):
        for j in range(n_a):
            cost[i, j] = two_impulse_cost(r1, r2, 1.0, r_p_grid[i], r_a_grid[j])
    dv_h = hohmann_dv_total(r1, r2, 1.0)
    X, Y = np.meshgrid(r_p_grid / r1, np.log10(r_a_grid / r2), indexing="ij")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4))
    cf = ax1.contourf(X, Y, np.log10(cost - dv_h + 1e-12), levels=30, cmap="viridis")
    ax1.plot(1.0, 0.0, "r*", ms=14, label="Hohmann (r_p=r1, r_a=r2)")
    ax1.set_xlabel("transfer periapsis r_p / r1")
    ax1.set_ylabel("log10(r_a / r2)")
    ax1.set_title("2-impulse cost excess over Hohmann (R = "
                  f"{r2 / r1:.2g}, log10 scale)")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    fig.colorbar(cf, ax=ax1)
    # 1D cut: tangential departure, r_a swept.
    r_a_curve = np.logspace(np.log10(r2), np.log10(200.0 * r2), 400)
    cost_curve = np.array(
        [two_impulse_cost(r1, r2, 1.0, r1, ra) for ra in r_a_curve])
    ax2.semilogx(r_a_curve / r2, (cost_curve - dv_h) / np.sqrt(1.0 / r1),
                 "C0-", lw=1.5)
    ax2.axvline(1.0, color="r", ls="--", lw=1.0)
    ax2.set_xlabel("transfer apoapsis r_a / r2")
    ax2.set_ylabel("(dv - dv_Hohmann) / v1")
    ax2.set_title("Tangent-departure family: minimum at r_a = r2")
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "optimality_surface.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    return paths


# --- Main -------------------------------------------------------------------


def main() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Closed-form self-consistency (also asserted in the test suite).
    split1 = hohmann_split(1.0, 2.0)
    split_mars = hohmann_split(AU_KM, MARS_A_AU * AU_KM, MU_SUN_KM3S2)

    sweep = ratio_sweep()
    rk4_cases = [
        validate_transfer_rk4(1.0, 1.5),
        validate_transfer_rk4(1.0, 6.409676),
        validate_transfer_rk4(1.0, 20.0),
        validate_transfer_rk4(1.0 / 1.3825, 1.0),  # Venus-like ratio (R = 1.3825)
        validate_transfer_rk4(2.0, 1.0),  # inward (R = 0.5, r2 < r1)
    ]
    opt_cases = [
        optimality_scan(1.0, 2.0),
        optimality_scan(1.0, 6.409676),
        optimality_scan(1.0, 20.0),
    ]
    real = real_cases()
    figures = make_figures(sweep, rk4_cases, opt_cases)

    print("=== Closed-form budget, r1=1, r2=2 (canonical units) ===")
    print(
        f"dv1 = {split1['dv1']:.8f}  dv2 = {split1['dv2']:.8f}  "
        f"total = {split1['dv_total']:.8f}  (dv/v1 = {split1['dv_total_over_v1']:.6f})"
    )
    a_t = split1["transfer"]["a"]
    print(f"transfer time = {split1['transfer_time']:.6f} "
          f"(half of period 2*pi*sqrt(a_t^3/mu) = "
          f"{2.0 * np.pi * np.sqrt(a_t**3 / MU):.6f})")
    print("=== R-sweep (dv_total/v1; asymptotes; peak) ===")
    print(
        f"R=1.0001: {sweep['cells'][2]['dv_total_over_v1']:.8f} vs (R-1)/2 = "
        f"{0.5 * 1e-4:.3e} (ratio {sweep['small_R_ratio_at_1p0001']:.6f})"
    )
    print(
        f"R=1e12:  {sweep['large_R_dv_over_v1_at_1e12']:.6f} vs sqrt(2)-1 = "
        f"{sweep['asymptote_large_R']:.6f} (rel gap {sweep['large_R_rel_gap']:.2e})"
    )
    pk = sweep["peak"]
    print(
        f"peak: R* = {pk['R_star']:.5f}, dv/v1 = {pk['dv_over_v1_peak']:.6f}"
        f"  (grid argmax {pk['grid_argmax_R']:.5f})"
    )
    print(
        f"inward/outward symmetry: max rel diff "
        f"{sweep['outward_inward_symmetry_max_rel_diff']:.2e}"
    )
    print("=== RK4 transfer validation ===")
    for case in rk4_cases:
        c = case["case"]
        a1, a2 = case["arrival_rk4"], case["arrival_analytic"]
        print(
            f"R={c['R']:.4g} e_t={c['e_transfer']:.4f} steps={c['steps_half_orbit']:>6}: "
            f"rk4 rel r-err {a1['rel_r_error']:.2e} v-err {a1['rel_v_error']:.2e} "
            f"radial/v {a1['v_radial_over_v']:.2e} | analytic {a2['rel_r_error']:.1e}/"
            f"{a2['rel_v_error']:.1e} | burn2 rel err {case['burn2']['rel_dv2_error']:.2e} "
            f"| circular orbit r-variation "
            f"{case['post_burn_circular_orbit']['radius_max_rel_variation']:.2e}"
        )
    print("=== Optimality (2-impulse family) ===")
    for o in opt_cases:
        g = o["grid"]
        print(
            f"R={o['R']:.4g}: min cost gap vs Hohmann {g['rel_gap_min_vs_hohmann']:.2e}"
            f"  argmin (r_p/r1={g['argmin_rp_over_r1']:.4f}, "
            f"r_a/r2={g['argmin_ra_over_r2']:.4f})  corner={g['argmin_is_hohmann_corner']}"
        )
    print("=== Real anchors ===")
    lg = real["leo_geo"]
    print(
        f"LEO(200 km)->GEO: dv1={lg['dv1']:.4f} dv2={lg['dv2']:.4f} "
        f"total={lg['dv_total']:.4f} km/s, t_tr={lg['transfer_time_hours']:.2f} h, "
        f"fuel(Isp300)={lg['fuel_frac_isp300']:.3f}"
    )
    em = real["earth_mars"]
    print(
        f"Earth->Mars: t_tr={em['transfer_time_days']:.2f} d, "
        f"v_inf dep={em['vinf_departure_km_s']:.4f} km/s, "
        f"v_inf arr={em['vinf_arrival_km_s']:.4f} km/s, "
        f"TMI from LEO={em['tmi_from_200km_leo_km_s']:.4f} km/s, "
        f"LEO->Mars orbit total={em['total_leo_to_mars_orbit_km_s']:.4f} km/s"
    )
    ev = real["earth_venus"]
    print(
        f"Earth->Venus (inward): t_tr={ev['transfer_time_days']:.2f} d, "
        f"v_inf dep={ev['vinf_departure_km_s']:.4f} km/s, "
        f"v_inf arr={ev['vinf_arrival_km_s']:.4f} km/s, "
        f"TMI from LEO={ev['tmi_from_200km_leo_km_s']:.4f} km/s"
    )

    result = {
        "closed_form_budget_r1_1_r2_2": split1,
        "earth_mars_budget_heliocentric": split_mars,
        "ratio_sweep": sweep,
        "rk4_transfer_validation": rk4_cases,
        "optimality_scan": opt_cases,
        "real_anchors": real,
        "constants": {
            "mu_earth_km3_s2": MU_EARTH_KM3S2,
            "r_earth_km": R_EARTH_KM,
            "geo_alt_km": GEO_ALT_KM,
            "leo_alt_km": LEO_ALT_KM,
            "au_km": AU_KM,
            "mu_sun_km3_s2": MU_SUN_KM3S2,
            "mars_a_au": MARS_A_AU,
            "venus_a_au": VENUS_A_AU,
        },
        "figures": [Path(p).name for p in figures],
    }
    path = save_json_result(
        RESULTS_DIR / "results.json",
        result,
        name="hohmann_transfer",
        description=(
            "Hohmann two-impulse transfer between coplanar circular orbits: "
            "closed-form dvs and transfer time verified from first principles, "
            "full radius-ratio sweep with asymptotes and the interior maximum, "
            "RK4 trajectory validation of the complete transfer, optimality "
            "scan over the two-impulse family, and real-system anchors "
            "(LEO->GEO, Earth->Mars, Earth->Venus) from IAU nominal constants."
        ),
    )
    print(f"\nSaved results -> {path}")
    return result


if __name__ == "__main__":
    main()