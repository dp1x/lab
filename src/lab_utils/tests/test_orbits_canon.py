"""Equivalence pins for lab_utils.orbits (canonical copies of donor machinery).

Each canon function is a verbatim transcription of a donor experiment's version
(Exp 008 groundtracks / Exp 009 j2Precession). These tests load the donors via
importlib and assert agreement on probe grids -- the transitional dual-source-of-
truth is honest because it is machine-checked. Donors remain frozen and untouched.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np

from lab_utils import orbits

_EXPERIMENTS_DIR = Path(__file__).resolve().parents[3] / "research" / "orbital-mechanics" / "experiments"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gt008 = _load("gt008_for_canon", _EXPERIMENTS_DIR / "groundtracks" / "experiment.py")
j2009 = _load("j2009_for_canon", _EXPERIMENTS_DIR / "j2Precession" / "experiment.py")


# --------------------------------------------------------------------------- #
# L2 -- pinned equivalence vs donors
# --------------------------------------------------------------------------- #
def test_constants_match_donors():
    assert orbits.MU_EARTH_KM3S2 == gt008.MU_EARTH_KM3S2 == j2009.MU_EARTH_KM3S2
    assert orbits.R_EARTH_KM == gt008.R_EARTH_KM
    assert orbits.OMEGA_EARTH_RAD_S == gt008.OMEGA_EARTH_RAD_S
    assert orbits.J2_EARTH == j2009.J2_EARTH


def test_solve_kepler_matches_donor_on_grid():
    M_grid = np.linspace(-np.pi, 3 * np.pi, 97)
    for e in (0.0, 0.3, 0.7, 0.99):
        mine = orbits.solve_kepler(M_grid.copy(), e)
        theirs = gt008.solve_kepler(M_grid.copy(), e)
        assert np.array_equal(np.asarray(mine), np.asarray(theirs)), f"e={e} not bit-equal"


def test_true_anomaly_from_E_matches_donor():
    E_grid = np.linspace(-np.pi, np.pi, 65)
    for e in (0.0, 0.4, 0.95):
        assert np.array_equal(
            np.asarray(orbits.true_anomaly_from_E(E_grid.copy(), e)),
            np.asarray(gt008.true_anomaly_from_E(E_grid.copy(), e)),
        )


def test_period_mean_motion_rotation_matrix_match_donor():
    for a in (6778.0, 26562.0):
        assert orbits.orbital_period(a) == gt008.orbital_period(a)
        assert orbits.mean_motion(a) == gt008.mean_motion(a)
    Q1 = orbits.rotation_matrix_313(0.3, 0.7, 1.1)
    Q2 = gt008.rotation_matrix_313(0.3, 0.7, 1.1)
    assert np.array_equal(Q1, Q2)
    assert np.max(np.abs(Q1.T @ Q1 - np.eye(3))) < 1e-15


def test_coe_to_rv_and_seed_state_match_donors():
    args = (8000.0, 0.2, 0.6, 1.0, 0.4, 0.9)
    r1, v1 = orbits.coe_to_rv_eci(*args[:5], args[5])
    r2, v2 = gt008.coe_to_rv_eci(*args[:5], args[5])
    assert np.array_equal(r1, r2) and np.array_equal(v1, v2)
    s1 = orbits.seed_state(*args[:6], mu=orbits.MU_EARTH_KM3S2)
    s2 = j2009.seed_state(*args[:6], mu=j2009.MU_EARTH_KM3S2)
    for a_, b_ in zip(s1, s2):
        assert np.array_equal(np.asarray(a_), np.asarray(b_))


def test_rv_to_coe_eci_matches_donor_with_guards():
    rng_states = [
        ((7000.0, 0.0, 500.0), (0.0, 7.5, 0.0)),  # circular equatorial: node/ecc guards active
        ((8000.0, 1000.0, 300.0), (-1.0, 6.0, 0.5)),
    ]
    for r, v in rng_states:
        mine = orbits.rv_to_coe_eci(np.array(r), np.array(v))
        theirs = j2009.rv_to_coe_eci(np.array(r), np.array(v))
        for key in ("a", "e", "inc", "Omega", "omega", "nu"):
            a_val, b_val = mine[key], theirs[key]
            if np.isnan(a_val) or np.isnan(b_val):
                assert np.isnan(a_val) and np.isnan(b_val), f"{key} guard mismatch"
            else:
                assert a_val == b_val, f"{key}: {a_val} != {b_val}"


def test_steps_per_orbit_matches_donor():
    for e in (0.0, 0.2, 0.5, 0.8, 0.95):
        assert orbits.steps_per_orbit(e) == j2009.steps_per_orbit(e)


# --------------------------------------------------------------------------- #
# Inline physics checks (theory duplicated here on purpose)
# --------------------------------------------------------------------------- #
def test_kepler_solver_round_trip_physics():
    """M(E) must reproduce M for extreme eccentricities."""
    for e in (0.0, 0.5, 0.98):
        E = np.linspace(-3.0, 3.0, 41)
        M = E - e * np.sin(E)
        E_back = np.asarray(orbits.solve_kepler(M, e))
        assert np.max(np.abs(E_back - E)) < 1e-10


def test_orbital_period_value_anchor():
    # T(7000 km) ~ 5755 s; closed form inline
    expected = 2 * math.pi * math.sqrt(7000.0**3 / 398600.4418)
    assert abs(orbits.orbital_period(7000.0) - expected) < 1e-9
    assert 5820 < expected < 5835


# --------------------------------------------------------------------------- #
# J2 RHS graduation (donor: Exp 009 inline accel inside propagate_3d_rk4_j2)
# --------------------------------------------------------------------------- #
def test_j2_rhs_trajectory_matches_donor_j2_on():
    """Full-force J2 Cowell trajectory: rk4_propagate(j2_rhs) vs donor loop."""
    from lab_utils.integrators import rk4_propagate

    mu = orbits.MU_EARTH_KM3S2
    a, e, inc = orbits.R_EARTH_KM + 550.0, 0.05, np.radians(51.6)
    r0, v0, _ = orbits.seed_state(a, e, inc, 0.3, 1.1, 0.7, mu)
    T = orbits.orbital_period(a, mu)
    t = np.linspace(0.0, 3.0 * T, 3 * 512 + 1)
    mine = rk4_propagate(orbits.j2_rhs(mu, orbits.J2_EARTH), t,
                         np.concatenate([r0, v0]))
    theirs = j2009.propagate_3d_rk4_j2(r0, v0, mu, t, orbits.J2_EARTH)
    rel = np.max(np.abs(mine - theirs)) / np.max(np.abs(theirs))
    assert rel < 1e-12, f"J2-on donor equivalence broken: {rel:.3e}"


def test_j2_rhs_bit_exact_vs_donor_j2_off():
    """j2 == 0 path must reproduce the Kepler-only donor loop bit-for-bit."""
    from lab_utils.integrators import rk4_propagate

    mu = orbits.MU_EARTH_KM3S2
    a, e, inc = orbits.R_EARTH_KM + 420.0, 0.01, np.radians(97.8)
    r0, v0, _ = orbits.seed_state(a, e, inc, 0.0, 0.0, 0.0, mu)
    T = orbits.orbital_period(a, mu)
    t = np.linspace(0.0, 2.0 * T, 2 * 512 + 1)
    x0 = np.concatenate([r0, v0])
    mine = rk4_propagate(orbits.j2_rhs(mu, 0.0), t, x0)
    theirs = j2009.propagate_3d_rk4_j2(r0, v0, mu, t, 0.0)
    assert np.array_equal(mine, theirs), "J2=0 path not bit-equal to donor"


def test_j2_rhs_acceleration_sign_structure():
    """Inline gradient sign check: bulge pulls equatorially, thins polar radius."""
    from lab_utils.integrators import rk4_step

    mu, j2 = orbits.MU_EARTH_KM3S2, orbits.J2_EARTH
    rhs = orbits.j2_rhs(mu, j2)
    v = np.array([0.0, 7.0, 0.0])
    # Equatorial point: J2 adds inward radial pull (|a| > mu/r^2).
    x_eq = np.array([orbits.R_EARTH_KM + 500.0, 0.0, 0.0])
    a_eq = rhs(0.0, np.concatenate([x_eq, v]))[3:]
    rm = np.linalg.norm(x_eq)
    assert np.linalg.norm(a_eq) > mu / rm**2
    # Polar point: z-pull reduced by the exact bulge factor (1 - 3 J2 (R/r)^2).
    x_pol = np.array([0.0, 0.0, orbits.R_EARTH_KM + 500.0])
    a_pol = rhs(0.0, np.concatenate([x_pol, v]))[3:]
    a_kep_z = -mu / (rm**2)
    expected_z = a_kep_z * (1.0 - 3.0 * j2 * (orbits.R_EARTH_KM / rm) ** 2)
    assert abs(a_pol[2] - expected_z) <= 1e-14 * abs(expected_z)


# --------------------------------------------------------------------------- #
# SSO inclination graduation (donor: Exp 012 orbitClasses.solve_sso_inclination)
# --------------------------------------------------------------------------- #
def test_sso_inclination_anchors_match_orbit_classes():
    """i_SSO at canonical altitudes matches the orbitClasses donor literals.

    Pinned by Exp 012 results.json and orbit-classes.md knowledge note.
    The lab_utils version is the same closed form with the no-silent-clip
    contract enforced via ValueError (no `np.clip` of cos_i).
    """
    oc012 = _load("oc012_for_sso", _EXPERIMENTS_DIR / "orbitClasses" / "experiment.py")
    for h_km, expected_deg in [
        (500.0, 97.401785943095),
        (600.0, 97.787646791197),
        (800.0, 98.603085267154),
    ]:
        a = orbits.R_EARTH_KM + h_km
        mine_deg = np.degrees(orbits.sso_inclination_rad(a, 0.0))
        donor = oc012.solve_sso_inclination(a, 0.0)
        assert donor["status"] == "OK"
        donor_deg = np.degrees(donor["incl_rad"])
        # The lab_utils implementation uses np.arccos(-ratio) directly while
        # the orbitClasses donor uses the same np.arccos; they should agree
        # to machine precision.
        assert abs(mine_deg - donor_deg) < 5e-5, (
            f"h={h_km} km: lab_utils {mine_deg} vs donor {donor_deg} deg"
        )
        # And the canonical pinned literal
        assert abs(mine_deg - expected_deg) < 5e-5, (
            f"h={h_km} km: {mine_deg} vs pinned {expected_deg} deg"
        )


def test_sso_inclination_retrograde_branch_strictly_above_90():
    for h_km in (500.0, 600.0, 700.0, 800.0):
        a = orbits.R_EARTH_KM + h_km
        i = orbits.sso_inclination_rad(a, 0.0)
        assert i > np.pi / 2, f"h={h_km} km: i = {np.degrees(i)} deg <= 90"


def test_sso_inclination_no_clip_raises_above_a_max():
    """Above a_max the cos_i exceeds [-1, 1] -> no real SSO solution exists.
    The lab_utils version raises ValueError; the orbitClasses version returns
    a typed `NO_REAL_SOLUTION` sentinel. Both reject the impossible case.
    """
    a_max = orbits.sso_existence_max_sma(0.0)
    # a_max itself gives cos i = -1 (i = 180 deg) -- the boundary
    # (within float64 ULP at the arccos singularity, ~3.4e-6 deg).
    i_at_max = orbits.sso_inclination_rad(a_max, 0.0)
    assert abs(np.degrees(i_at_max) - 180.0) < 1e-4
    # a > a_max must raise
    import pytest
    with pytest.raises(ValueError, match="no real SSO solution"):
        orbits.sso_inclination_rad(a_max * 1.001, 0.0)


def test_sso_existence_max_sma_value_pin():
    """a_max at e=0 is the lab's pinned 12352.505076 km (Exp 012 headline)."""
    a_max = orbits.sso_existence_max_sma(0.0)
    assert abs(a_max - 12352.505076) < 1e-3, f"a_max = {a_max} km"


def test_sso_target_deg_day_pinned_literal():
    """SSO_TARGET_DEG_DAY is the mean-solar-year rate 360/365.2422.

    Documented convention firewall: sidereal year (365.256) gives
    0.98560912 deg/day, Julian year (365.25) gives 0.98564685 deg/day;
    these would shift i_SSO at 600 km by 3.0e-4 and 1.7e-4 deg, caught
    at 5e-5 tolerance. Tropical year (365.24219) is 2.1e-7 below behavior
    discrimination; pinned by literal.
    """
    expected = 360.0 / 365.2422
    assert abs(orbits.SSO_TARGET_DEG_DAY - expected) < 1e-15
    # and the textbook value
    assert abs(orbits.SSO_TARGET_DEG_DAY - 0.985647332099) < 1e-9
