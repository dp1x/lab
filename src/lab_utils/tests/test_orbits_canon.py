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
