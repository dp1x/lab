"""Kepler orbit validation: Newtonian two-body gravity reproduces Kepler's laws.

Validates a deterministic fixed-step RK4 propagator for the two-body problem
against the closed-form Keplerian solution:

  - Kepler I:   bound orbits are ellipses (pointwise match of the propagated
                state with the analytic conic solution; conic equation holds
                along the trajectory).
  - Kepler II:  equal areas are swept in equal times (areal velocity dA/dt =
                h/2 where h is the specific angular momentum).
  - Kepler III: T^2 proportional to a^3, with constant 4*pi^2/mu.

Plus: conservation of specific energy and specific angular momentum, a
propagator convergence study, and a real-units Earth orbit whose period is
compared against the IAU-defined (GM)_Sun and the astronomical unit.

References: R. R. Bate, D. D. Mueller, J. E. White, "Fundamentals of
Astrodynamics", Dover, 1971; H. D. Curtis, "Orbital Mechanics for Engineering
Students", 4th ed., Elsevier, 2021; C. D. Murray, S. F. Dermott, "Solar System
Dynamics", Cambridge UP, 1999.

Canonical units: mu = 1, a = 1 gives period T = 2*pi.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lab_utils.metrics import convergence_rate, max_abs_error
from lab_utils.results import save_json_result

# --- Physical constants (verified against primary sources, 2026-08-13) -----

# (GM)_Sun IAU 2015 Resolution B3 nominal heliocentric gravitational constant.
MU_SUN_KM3_S2 = 1.3271244e11  # km^3 / s^2  [iau.org, IAUGA2015-Resolution-B3]
# Astronomical unit, IAU 2012 Resolution B2: exactly 149597870.7 km.
AU_KM = 1.495978707e8  # km
# Earth mean orbital eccentricity, JPL SSD (Standish & Williams 1992),
# mean ecliptic J2000, EM barycenter.
EARTH_E = 0.01671123
# Sidereal year, Astronomical Almanac for the Year 2025 (USNO/HMNAO), p. C2.
SIDEREAL_YEAR_DAYS = 365.256363

# --- Experiment parameters -------------------------------------------------

MU = 1.0  # canonical gravitational parameter
A_DEFAULT = 1.0  # canonical semi-major axis
ORBITS_K1 = 5  # orbits propagated for the Kepler-I pointwise check
STEPS_PER_ORBIT = 512  # stepsize h = T / STEPS_PER_ORBIT
K2_E = 0.6  # eccentricity for the equal-areas study
K2_INTERVALS = 12  # equal-time intervals partitioning one orbit
K2_STEPS_PER_INTERVAL = 64  # propagation steps per interval
CONS_ORBITS = 10  # conservation-law horizon (orbits)
# Kepler-III sweep grid (mean anomaly and eccentricity).
K3_A = [0.5, 1.0, 2.0, 4.0, 8.0]
K3_E = [0.10, 0.30, 0.60, 0.85]
K3_ORBITS = 2.05  # just over two orbits: two periapsis passages are measured

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"


# --- Two-body dynamics -----------------------------------------------------


def acceleration(r: np.ndarray, mu: float = MU) -> np.ndarray:
    """Central-force acceleration a = -mu * r / |r|^3 (planar)."""
    rr = np.hypot(r[0], r[1])
    return -mu * r / rr**3


def propagate_rk4(
    r0: np.ndarray,
    v0: np.ndarray,
    mu: float,
    t: np.ndarray,
) -> np.ndarray:
    """Classical RK4 propagation of the planar two-body problem.

    Returns states of shape (len(t), 4) ordered [x, y, vx, vy]. RK4 is not
    symplectic; for this central-force problem the measured relative drift of
    h = |r x v| is O(1e-11) and of the specific energy O(1e-9) over 10 orbits
    (discretization-order smallness, not round-off-level exactness).
    """
    h = t[1] - t[0]
    state = np.array([r0[0], r0[1], v0[0], v0[1]], dtype=float)
    states = np.empty((len(t), 4))
    states[0] = state
    for i in range(1, len(t)):
        r, v = state[:2], state[2:]
        k1v = acceleration(r, mu)
        k2v = acceleration(r + 0.5 * h * v, mu)
        k3v = acceleration(r + 0.5 * h * v + 0.25 * h * h * k1v, mu)
        k4v = acceleration(r + h * v + 0.5 * h * h * k2v, mu)
        v_new = v + (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
        r_new = r + h * v + (h * h / 6.0) * (k1v + k2v + k3v)
        state = np.concatenate((r_new, v_new))
        states[i] = state
    return states


# --- Orbital elements and the analytic solution ----------------------------


def orbital_elements(r: np.ndarray, v: np.ndarray, mu: float = MU) -> dict[str, float]:
    """Classical two-body elements from a state vector (planar, perifocal).

    Returns specific energy, specific angular momentum, semi-major axis,
    eccentricity (from the eccentricity vector) and semi-latus rectum.
    """
    rr = np.hypot(r[0], r[1])
    vv = np.hypot(v[0], v[1])
    eps = 0.5 * vv**2 - mu / rr
    h_vec = np.array([0.0, 0.0, r[0] * v[1] - r[1] * v[0]])
    h = h_vec[2]
    a = -mu / (2.0 * eps)
    r3 = np.array([r[0], r[1], 0.0])
    v3 = np.array([v[0], v[1], 0.0])
    e_vec = np.cross(v3, h_vec) / mu - r3 / rr
    e = float(np.hypot(e_vec[0], e_vec[1]))
    return {
        "energy": float(eps),
        "angular_momentum": float(h),
        "semi_major_axis": float(a),
        "eccentricity": e,
        "semi_latus_rectum": float(h**2 / mu),
        "eccentricity_vector": e_vec.tolist(),
    }


def solve_kepler(M: np.ndarray, e: float, tol: float = 1e-14) -> np.ndarray:
    """Solve Kepler's equation M = E - e sin E for the eccentric anomaly.

    Deterministic Newton iteration started from E0 = M + e sin M (converges
    for all e in [0, 1)); fixed iteration cap of 100.
    """
    E = M + e * np.sin(M)
    for _ in range(100):
        f = E - e * np.sin(E) - M
        E = E - f / (1.0 - e * np.cos(E))
        if np.max(np.abs(f)) < tol:
            break
    if np.max(np.abs(f)) >= tol:
        raise ValueError(f"Kepler's equation did not converge (e={e}, max|f|={np.max(np.abs(f)):.2e})")
    return E


def kepler_solution(a: float, e: float, mu: float, t: np.ndarray) -> np.ndarray:
    """Closed-form Keplerian state [x, y, vx, vy] at times t.

    Periapsis at t = 0, orbit in the xy-plane (perifocal frame). From
    M = n t, E via Kepler's equation: r = a(cosE - e), a*sqrt(1-e^2) sinE.
    """
    n = np.sqrt(mu / a**3)
    M = n * t
    E = solve_kepler(M, e)
    cosE, sinE = np.cos(E), np.sin(E)
    denom = 1.0 - e * cosE
    x = a * (cosE - e)
    y = a * np.sqrt(1.0 - e**2) * sinE
    vx = -a * n * sinE / denom
    vy = a * np.sqrt(1.0 - e**2) * n * cosE / denom
    return np.column_stack([x, y, vx, vy])


def initial_state(a: float, e: float, mu: float = MU) -> tuple[np.ndarray, np.ndarray]:
    """Perifocal initial conditions: periapsis r_p = a(1-e), tangential speed
    from the vis-viva equation v_p = sqrt(mu (1+e) / (a (1-e)))."""
    r_p = a * (1.0 - e)
    v_p = np.sqrt(mu * (1.0 + e) / (a * (1.0 - e)))
    return np.array([r_p, 0.0]), np.array([0.0, v_p])


def true_anomaly(r: np.ndarray, e_vec: np.ndarray) -> float:
    """True anomaly from position and the eccentricity vector (planar)."""
    e = np.hypot(e_vec[0], e_vec[1])
    rr = np.hypot(r[0], r[1])
    if e < 1e-14:
        return 0.0
    cos_nu = float(np.dot(e_vec[:2], r) / (e * rr))
    sin_nu = float(e_vec[0] * r[1] - e_vec[1] * r[0]) / (e * rr)
    return float(np.arctan2(sin_nu, cos_nu))


# --- Kepler's law checks ---------------------------------------------------


def case_steps_per_orbit(e: float, base: int = STEPS_PER_ORBIT) -> int:
    """Steps per orbit that resolve the periapsis passage for eccentricity e.

    The periapsis pass occupies a fraction ~ (1-e)^(3/2) of the period, so
    uniform-time RK4 needs steps ~ (1-e)^(-3/2) to keep per-step error bounded
    near periapsis. Measured: at e = 0.85, 512 steps/orbit give ~100% position
    error over 5 orbits; 4096 give ~8e-5.
    """
    return max(base, int(np.ceil(base / (1.0 - e) ** 1.5)))


def run_case(
    a: float,
    e: float,
    mu: float,
    orbits: float,
    steps_per_orbit: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate and return (t, states) for `orbits` orbits of case (a, e).

    steps_per_orbit defaults to the periapsis-resolved value.
    """
    period = 2.0 * np.pi * np.sqrt(a**3 / mu)
    if steps_per_orbit is None:
        steps_per_orbit = case_steps_per_orbit(e)
    n_steps = max(2, int(round(orbits * steps_per_orbit)))
    t = np.linspace(0.0, orbits * period, n_steps + 1)
    r0, v0 = initial_state(a, e, mu)
    states = propagate_rk4(r0, v0, mu, t)
    return t, states


def kepler1_shape_check(a: float, e: float, mu: float) -> dict:
    """Kepler I: the propagator reproduces the analytic conic.

    (a) Pointwise: max relative position error between the RK4 trajectory and
        the closed-form Kepler solution over ORBITS_K1 orbits.
    (b) Conic:   with elements from the t=0 state, |r| must satisfy the conic
        equation r = p / (1 + e cos nu) along the entire trajectory.
    """
    t, states = run_case(a, e, mu, ORBITS_K1)
    ana = kepler_solution(a, e, mu, t)
    pos_err = np.hypot(states[:, 0] - ana[:, 0], states[:, 1] - ana[:, 1])
    ana_r = np.hypot(ana[:, 0], ana[:, 1])
    max_rel_err = float(np.max(pos_err / ana_r))

    r0, v0 = initial_state(a, e, mu)
    els = orbital_elements(r0, v0, mu)
    e_vec = np.array(els["eccentricity_vector"])
    p, ee = els["semi_latus_rectum"], els["eccentricity"]
    r_num = np.hypot(states[:, 0], states[:, 1])
    conic_errors = np.empty(len(states))
    for i in range(len(states)):
        nu = true_anomaly(states[i, :2], e_vec)
        conic_errors[i] = abs(r_num[i] - p / (1.0 + ee * np.cos(nu))) / r_num[i]
    return {
        "a": a,
        "e_input": e,
        "e_from_eccentricity_vector": els["eccentricity"],
        "orbits": ORBITS_K1,
        "steps_per_orbit": case_steps_per_orbit(e),
        "total_steps": len(states) - 1,
        "max_rel_pointwise_error": max_rel_err,
        "max_rel_conic_error": float(np.max(conic_errors)),
    }


def circular_orbit_check(a: float, mu: float) -> dict:
    """Kepler I for the degenerate e = 0 case: constant radius, e_vec ~ 0."""
    t, states = run_case(a, 0.0, mu, ORBITS_K1, STEPS_PER_ORBIT)
    r = np.hypot(states[:, 0], states[:, 1])
    r0, v0 = initial_state(a, 0.0, mu)
    els = orbital_elements(r0, v0, mu)
    return {
        "a": a,
        "e_input": 0.0,
        "max_rel_radius_variation": float(np.max(np.abs(r - a)) / a),
        "measured_eccentricity": els["eccentricity"],
        "orbits": ORBITS_K1,
        "steps_per_orbit": STEPS_PER_ORBIT,
    }


def kepler2_equal_areas(a: float, e: float, mu: float) -> dict:
    """Kepler II: equal areas in equal times, and dA/dt = h/2.

    One orbit is split into K2_INTERVALS equal-time intervals; the sector area
    of each interval is the triangle sum 1/2 * |r_i x r_{i+1}|, which for a
    central-force trajectory converges to the exact sector area.
    """
    period = 2.0 * np.pi * np.sqrt(a**3 / mu)
    interval_t = period / K2_INTERVALS
    h = interval_t / K2_STEPS_PER_INTERVAL
    n_steps = K2_INTERVALS * K2_STEPS_PER_INTERVAL
    t = np.linspace(0.0, period, n_steps + 1)
    r0, v0 = initial_state(a, e, mu)
    states = propagate_rk4(r0, v0, mu, t)

    areas = []
    for k in range(K2_INTERVALS):
        lo = k * K2_STEPS_PER_INTERVAL
        hi = (k + 1) * K2_STEPS_PER_INTERVAL
        seg = 0.0
        for i in range(lo, hi):
            seg += 0.5 * abs(
                states[i, 0] * states[i + 1, 1] - states[i, 1] * states[i + 1, 0]
            )
        areas.append(seg)

    mean_area = float(np.mean(areas))
    max_rel_dev = float(np.max(np.abs(np.array(areas) - mean_area)) / mean_area)

    # Full-orbit sector area must be pi * a * b = pi * a^2 * sqrt(1 - e^2).
    full_area = float(sum(areas))
    ellipse_area = np.pi * a**2 * np.sqrt(1.0 - e**2)

    # Areal velocity: cumulative swept area grows linearly with time.
    cum = np.cumsum(np.array(areas))
    t_mid = (np.arange(K2_INTERVALS) + 0.5) * interval_t
    slope = float(np.polyfit(t_mid, cum, 1)[0])
    h_ref = orbital_elements(r0, v0, mu)["angular_momentum"]
    return {
        "a": a,
        "e": e,
        "intervals": K2_INTERVALS,
        "steps_per_interval": K2_STEPS_PER_INTERVAL,
        "interval_areas": areas,
        "mean_interval_area": mean_area,
        "max_rel_interval_deviation": max_rel_dev,
        "full_orbit_sector_area": full_area,
        "ellipse_area_pi_ab": float(ellipse_area),
        "full_orbit_area_rel_error": abs(full_area - ellipse_area) / ellipse_area,
        "areal_velocity_slope": slope,
        "expected_areal_velocity_h_over_2": h_ref / 2.0,
        "areal_velocity_rel_error": abs(slope - h_ref / 2.0) / (h_ref / 2.0),
    }


def measure_period(t: np.ndarray, states: np.ndarray, e: float) -> float:
    """Orbital period from successive periapsis passages (e > 0), or from the
    unwrapped polar angle (e = 0), with sub-step interpolation."""
    r = np.hypot(states[:, 0], states[:, 1])
    if e > 0.0:
        minima = [
            i
            for i in range(1, len(r) - 1)
            if r[i] < r[i - 1] and r[i] < r[i + 1]
        ]
        if len(minima) < 2:
            raise ValueError("fewer than two periapsis passages detected")
        refined = []
        for i in minima:
            denom = r[i - 1] - 2.0 * r[i] + r[i + 1]
            if abs(denom) < 1e-300:
                refined.append(t[i])
            else:
                refined.append(
                    t[i] + (t[1] - t[0]) * (r[i - 1] - r[i + 1]) / (2.0 * denom)
                )
        return refined[1] - refined[0]
    theta = np.unwrap(np.arctan2(states[:, 1], states[:, 0]))
    hstep = t[1] - t[0]
    crossing = []
    for k in (1, 2):
        target = 2.0 * np.pi * k
        idx = int(np.searchsorted(theta, target))
        frac = (target - theta[idx - 1]) / (theta[idx] - theta[idx - 1])
        crossing.append(t[idx - 1] + frac * hstep)
    return crossing[1] - crossing[0]


def kepler3_harmonic_law(mu: float) -> dict:
    """Kepler III: T^2 = (4 pi^2 / mu) a^3 across a sweep of (a, e)."""
    cells = []
    for a in K3_A:
        for e in K3_E:
            period_theory = 2.0 * np.pi * np.sqrt(a**3 / mu)
            t, states = run_case(a, e, mu, K3_ORBITS)
            t_meas = measure_period(t, states, e)
            cells.append(
                {
                    "a": a,
                    "e": e,
                    "period_theory": period_theory,
                    "period_measured": t_meas,
                    "period_rel_error": abs(t_meas - period_theory) / period_theory,
                    "T2_over_a3": t_meas**2 / a**3,
                }
            )
    log_a = np.array([c["a"] for c in cells])
    log_T = np.array([c["period_measured"] for c in cells])
    slope, intercept = np.polyfit(np.log10(log_a), np.log10(log_T), 1)
    t2a3 = np.array([c["T2_over_a3"] for c in cells])
    expected = 4.0 * np.pi**2 / mu
    return {
        "mu": mu,
        "cells": cells,
        "loglog_slope": float(slope),
        "fit_intercept": float(intercept),
        "expected_4pi2_over_mu": float(expected),
        "T2_over_a3_min": float(np.min(t2a3)),
        "T2_over_a3_max": float(np.max(t2a3)),
        "T2_over_a3_max_rel_err": float(
            np.max(np.abs(t2a3 - expected)) / expected
        ),
    }


def conservation_study(a: float, e: float, mu: float) -> dict:
    """Specific energy and angular momentum drift over CONS_ORBITS orbits."""
    t, states = run_case(a, e, mu, CONS_ORBITS)
    r0, v0 = initial_state(a, e, mu)
    els0 = orbital_elements(r0, v0, mu)
    eps0, h0 = els0["energy"], els0["angular_momentum"]
    eps_max, h_max = 0.0, 0.0
    for i in range(len(states)):
        r, v = states[i, :2], states[i, 2:]
        els = orbital_elements(r, v, mu)
        eps_max = max(eps_max, abs(els["energy"] - eps0) / abs(eps0))
        h_max = max(h_max, abs(els["angular_momentum"] - h0) / abs(h0))
    return {
        "a": a,
        "e": e,
        "orbits": CONS_ORBITS,
        "steps_per_orbit": case_steps_per_orbit(e),
        "total_steps": len(states) - 1,
        "energy_analytic": eps0,
        "angular_momentum_analytic": h0,
        "max_rel_energy_drift": eps_max,
        "max_rel_angular_momentum_drift": h_max,
    }


def propagator_convergence(a: float, e: float, mu: float) -> dict:
    """Global position error of RK4 vs the analytic solution vs stepsize."""
    steps_list = [128, 256, 512, 1024, 2048]
    errors = []
    for spp in steps_list:
        t, states = run_case(a, e, mu, 1.0, spp)
        ana = kepler_solution(a, e, mu, t)
        errors.append(
            max_abs_error(states[:, :2], ana[:, :2])
        )
    stepsizes = np.array([1.0 / s for s in steps_list])
    return {
        "a": a,
        "e": e,
        "steps_per_orbit": steps_list,
        "stepsizes": stepsizes.tolist(),
        "max_errors": errors,
        "measured_order": convergence_rate(
            np.array(errors), stepsizes
        ).tolist(),
        "theoretical_order": 4.0,
    }


def earth_case() -> dict:
    """Real-units anchor: Earth orbit about the Sun with IAU constants.

    Two-body idealization (Sun + Earth, GM from IAU 2015 B3 nominal), a = 1 au
    (IAU 2012 B2), e from JPL SSD. The predicted period 2*pi*sqrt(a^3/mu) =
    365.256898 d is compared against the sidereal year (Astronomical Almanac
    2025, 365.256363 d); the +1.5e-6 residual is physical (Earth-Moon system
    mass and mean-vs-nominal a), not numerical.
    """
    a, e, mu = AU_KM, EARTH_E, MU_SUN_KM3_S2
    period_theory = 2.0 * np.pi * np.sqrt(a**3 / mu)
    orbit = 3.05
    n_steps = int(round(orbit * 1024))
    t = np.linspace(0.0, orbit * period_theory, n_steps + 1)
    r0, v0 = initial_state(a, e, mu)
    states = propagate_rk4(r0, v0, mu, t)
    t_meas = measure_period(t, states, e)
    e_meas = orbital_elements(r0, v0, mu)["eccentricity"]
    return {
        "mu_km3_s2": mu,
        "a_km": a,
        "e_input": e,
        "e_measured": e_meas,
        "period_theory_s": period_theory,
        "period_theory_days": period_theory / 86400.0,
        "period_measured_days": t_meas / 86400.0,
        "period_rel_error": abs(t_meas - period_theory) / period_theory,
        "sidereal_year_days": SIDEREAL_YEAR_DAYS,
        "rel_diff_vs_sidereal_year": (period_theory / 86400.0 - SIDEREAL_YEAR_DAYS)
        / SIDEREAL_YEAR_DAYS,
    }


# --- Figures ---------------------------------------------------------------


def make_figures(k2: dict, k3: dict, conv: dict) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (RESULTS_DIR / "figures").mkdir(parents=True, exist_ok=True)
    paths = []

    # 1. Orbit shapes across eccentricities.
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    for e, color in zip([0.0, 0.3, 0.6, 0.85], ["C0", "C1", "C2", "C3"]):
        t, states = run_case(1.0, e, MU, 1.2, 512)
        if e == 0.0:
            ax.plot(states[:, 0], states[:, 1], color=color, label=f"e = {e}")
        else:
            ax.plot(states[:, 0], states[:, 1], color=color, label=f"e = {e}")
    ax.plot(0, 0, "k*", ms=12, label="focus (central body)")
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Kepler I: two-body bound orbits are ellipses (a = 1)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    p = RESULTS_DIR / "figures" / "orbit_shapes.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 2. Equal-areas fan + linear swept area.
    a, e, mu = 1.0, K2_E, MU
    period = 2.0 * np.pi * np.sqrt(a**3 / mu)
    n_steps = K2_INTERVALS * K2_STEPS_PER_INTERVAL
    t, states = run_case(a, e, mu, 1.0, n_steps)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))
    ax1.plot(states[:, 0], states[:, 1], "C0", lw=1)
    for k in range(K2_INTERVALS + 1):
        i = k * K2_STEPS_PER_INTERVAL
        ax1.plot([0, states[i, 0]], [0, states[i, 1]], "k-", lw=0.4, alpha=0.5)
    ax1.plot(0, 0, "k*", ms=10)
    ax1.set_aspect("equal")
    ax1.set_title(f"Kepler II: one orbit split into {K2_INTERVALS} equal-time sectors (e = {e})")
    ax1.grid(True, alpha=0.3)
    interval_t = period / K2_INTERVALS
    t_mid = (np.arange(K2_INTERVALS) + 0.5) * interval_t
    cum = np.cumsum(np.array(k2["interval_areas"]))
    ax2.plot(t_mid, cum, "o-", label="cumulative swept area")
    ax2.plot(t_mid, k2["expected_areal_velocity_h_over_2"] * t_mid, "r--",
             label="(h/2) t, h from state")
    ax2.set_xlabel("time t")
    ax2.set_ylabel("swept area")
    ax2.set_title("Swept area grows linearly: dA/dt = h/2")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "equal_areas.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 3. Kepler III: log T vs log a, and T^2/a^3 vs 4 pi^2/mu.
    cells = k3["cells"]
    log_a = np.log10(np.array([c["a"] for c in cells]))
    log_T = np.log10(np.array([c["period_measured"] for c in cells]))
    slope = k3["loglog_slope"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))
    ax1.plot(log_a, log_T, "o", ms=4)
    fit_a = np.linspace(log_a.min(), log_a.max(), 50)
    ax1.plot(fit_a, slope * fit_a + k3["fit_intercept"], "r--",
             label=f"fit slope = {slope:.4f} (theory 1.5)")
    ax1.set_xlabel("log10 a")
    ax1.set_ylabel("log10 T_measured")
    ax1.set_title("Kepler III: T ~ a^{3/2}")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    t2a3 = np.array([c["T2_over_a3"] for c in cells])
    ax2.semilogy(range(len(cells)), t2a3, "o", ms=4)
    ax2.axhline(k3["expected_4pi2_over_mu"], color="r", ls="--",
                label="4 pi^2 / mu")
    ax2.set_xlabel("sweep cell (a, e)")
    ax2.set_ylabel("T^2 / a^3")
    ax2.set_title("Kepler III: T^2/a^3 = 4 pi^2 / mu for all (a, e)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "kepler3.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 4. Conservation of energy and angular momentum.
    t, states = run_case(1.0, 0.6, MU, CONS_ORBITS, STEPS_PER_ORBIT)
    r0, v0 = initial_state(1.0, 0.6, MU)
    eps0 = orbital_elements(r0, v0, MU)["energy"]
    h0 = orbital_elements(r0, v0, MU)["angular_momentum"]
    dE = np.empty(len(t))
    dH = np.empty(len(t))
    for i in range(len(states)):
        els = orbital_elements(states[i, :2], states[i, 2:], MU)
        dE[i] = abs(els["energy"] - eps0) / abs(eps0)
        dH[i] = abs(els["angular_momentum"] - h0) / abs(h0)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.semilogy(t / (2 * np.pi), dE, label="specific energy |dE/E|")
    ax.semilogy(t / (2 * np.pi), dH, label="angular momentum |dh/h|")
    ax.set_xlabel("orbits elapsed")
    ax.set_ylabel("relative deviation")
    ax.set_title("Conservation of invariants over 10 orbits (e = 0.6)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "conservation.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 5. Propagator convergence.
    fig, ax = plt.subplots(figsize=(7.5, 5))
    stepsizes = np.asarray(conv["stepsizes"])
    ax.loglog(stepsizes, conv["max_errors"], "o-",
              label="RK4 vs analytic solution (one orbit, e=0.6)")
    h_ref = stepsizes[1]
    serr = conv["max_errors"][1]
    ax.loglog(stepsizes, serr * (stepsizes / h_ref) ** 4, "r--",
              label="slope 4 reference")
    ax.set_xlabel("stepsize h (fraction of period)")
    ax.set_ylabel("max position error")
    ax.set_title("RK4 propagator convergence: O(h^4)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "propagator_convergence.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    return paths


# --- Main ------------------------------------------------------------------


def main() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    k1_cases = [
        kepler1_shape_check(1.0, e, MU) for e in (0.3, 0.6, 0.85)
    ]
    circle = circular_orbit_check(1.0, MU)
    k2 = kepler2_equal_areas(1.0, K2_E, MU)
    k3 = kepler3_harmonic_law(MU)
    cons = conservation_study(1.0, 0.6, MU)
    conv = propagator_convergence(1.0, 0.6, MU)
    earth = earth_case()
    figures = make_figures(k2, k3, conv)

    print("=== Kepler I: propagate vs analytic conic (a=1, 5 orbits) ===")
    for case in k1_cases:
        print(
            f"e={case['e_input']:<5} pointwise rel err {case['max_rel_pointwise_error']:.2e}"
            f"  conic rel err {case['max_rel_conic_error']:.2e}"
            f"  e_meas={case['e_from_eccentricity_vector']:.6f}"
        )
    print(
        f"circular: max radius variation {circle['max_rel_radius_variation']:.2e}"
        f"  e_meas={circle['measured_eccentricity']:.2e}"
    )
    print("=== Kepler II (a=1, e=0.6) ===")
    print(
        f"intervals={k2['intervals']}  max rel area deviation {k2['max_rel_interval_deviation']:.2e}"
        f"  full-orbit area vs pi*a*b {k2['full_orbit_area_rel_error']:.2e}"
        f"  dA/dt vs h/2 {k2['areal_velocity_rel_error']:.2e}"
    )
    print("=== Kepler III (mu=1) ===")
    print(
        f"log-log slope {k3['loglog_slope']:.5f} (theory 1.5)"
        f"  T^2/a^3 rel err vs 4 pi^2/mu {k3['T2_over_a3_max_rel_err']:.2e}"
    )
    print("=== Conservation (e=0.6, 10 orbits) ===")
    print(
        f"max rel energy drift {cons['max_rel_energy_drift']:.2e}"
        f"  max rel h drift {cons['max_rel_angular_momentum_drift']:.2e}"
    )
    print("=== Propagator convergence ===")
    avg_order = float(np.mean(conv["measured_order"]))
    print(f"measured order {avg_order:.3f} (theory 4)")
    print("=== Earth orbit (real units) ===")
    print(
        f"T_pred {earth['period_theory_days']:.6f} d, "
        f"T_meas {earth['period_measured_days']:.6f} d "
        f"(rel err {earth['period_rel_error']:.2e}), "
        f"vs sidereal year {earth['sidereal_year_days']} d "
        f"(rel diff {earth['rel_diff_vs_sidereal_year']:.2e}), "
        f"e_meas {earth['e_measured']:.6f}"
    )

    result = {
        "kepler1": {"cases": k1_cases, "circular": circle},
        "kepler2": k2,
        "kepler3": k3,
        "conservation": cons,
        "propagator_convergence": conv,
        "earth_anchor": earth,
        "figures": [Path(p).name for p in figures],
    }
    path = save_json_result(
        RESULTS_DIR / "results.json",
        result,
        name="kepler_orbit_validation",
        description=(
            "Validation of a two-body RK4 propagator against Kepler's three "
            "laws, conservation invariants, and an IAU-real-units Earth anchor."
        ),
    )
    print(f"\nSaved results -> {path}")
    return result


if __name__ == "__main__":
    main()