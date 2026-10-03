"""Equivalence + regression tests for ``lab_utils.ephemeris``.

Every function here was *moved*, not rewritten. These tests prove the move was
faithful by comparing against the donor module (Mission 3's ``mission_experiment.py``)
on the real byte-pinned DE441 snapshots, and they lock down the two historical
defect classes that lived in this code:

1. the IAU-1976 precession ``_rot3`` transpose (audit-019), and
2. the DE441 interpolation sign-reversed weights (2026-09), which moved the
   Mission-2 cross coefficient by ~5x and manufactured a spurious 3.9-sigma
   "detection".
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from lab_utils import ephemeris as EPH  # noqa: E402

_DONOR_PATH = (ROOT / "research/orbital-mechanics/missions"
              / "mission_j2_precession_modulated_lunisolar_term"
              / "mission_experiment.py")

SUN_SHA = "f2c4f04824b07d02638b3c6cac1e72385bedf56c7575709ab1158391f92889db"
MOON_SHA = "aee8509932c1ea169df1f518feebd0169cf1fd8ebd07f153f6fe80f4ee7d8c59"


def _donor():
    if "donor_me_eph" not in sys.modules:
        spec = importlib.util.spec_from_file_location("donor_me_eph", _DONOR_PATH)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["donor_me_eph"] = mod
        spec.loader.exec_module(mod)
    return sys.modules["donor_me_eph"]


def _snapshots():
    sun_p, moon_p = EPH.de441_snapshot_paths(ROOT)
    return EPH.load_snapshot(sun_p), EPH.load_snapshot(moon_p)


# --- provenance ------------------------------------------------------------ #

def test_snapshot_paths_resolve_inside_the_repo():
    sun_p, moon_p = EPH.de441_snapshot_paths(ROOT)
    assert sun_p.is_file() and moon_p.is_file()
    assert ROOT in sun_p.parents and ROOT in moon_p.parents, \
        "authoritative data must live in the repo, never on scratch"


def test_snapshot_checksums_match_the_committed_provenance():
    sun, moon = _snapshots()
    assert sun["sha256"] == SUN_SHA
    assert moon["sha256"] == MOON_SHA


def test_snapshot_load_matches_donor():
    d = _donor()
    sun, moon = _snapshots()
    ds = d.load_snapshot(d.SUN_SNAPSHOT)
    assert sun["n_points"] == ds["n_points"]
    assert np.array_equal(sun["t_s"], ds["t_s"])
    assert np.array_equal(sun["r_eci_km"], ds["r_eci_km"])


# --- precession ------------------------------------------------------------ #

def test_precession_matrix_is_orthonormal_and_right_handed():
    for t_s in (0.0, 8.0e8, 1.2e9):
        R = EPH.precession_j2000_to_mod(t_s)
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-14)
        assert abs(np.linalg.det(R) - 1.0) < 1e-14


def test_precession_matches_donor_and_is_identity_at_j2000():
    d = _donor()
    for t_s in (0.0, 3.15e8, 9.2e8):
        assert np.array_equal(EPH.precession_j2000_to_mod(t_s),
                              d.precession_j2000_to_mod(t_s))
    assert np.allclose(EPH.precession_j2000_to_mod(0.0), np.eye(3), atol=1e-15)


def test_precession_actually_rotates_at_later_epochs():
    """Guards the audit-019 transpose: a transpose is still orthonormal, so
    orthonormality alone cannot catch it. The composed convention
    ``_rot3(-z) @ _rot2(theta) @ _rot3(-zeta)`` yields a NEGATIVE rotation
    angle; a transposed _rot3 flips this sign, which is the discriminator.
    """
    R = EPH.precession_j2000_to_mod(9.2e8)          # ~2029, T ~ 0.29 centuries
    ang = math.atan2(R[1, 0], R[0, 0])
    assert ang < 0.0, "precession angle sign flipped -> _rot3 transpose regression"
    # General precession in longitude is ~50.3"/yr; 29 yr -> ~0.4 deg ~ 7e-3 rad.
    assert -0.02 < ang < -1e-3, f"precession angle {math.degrees(ang):.4f} deg out of band"


# --- interpolation --------------------------------------------------------- #

def test_interp_matches_donor_bitwise_on_many_epochs():
    d = _donor()
    sun, _ = _snapshots()
    donor_sun = d.load_snapshot(d.SUN_SNAPSHOT)   # loaded once, not per sample
    for t in np.linspace(sun["t_s"][0], sun["t_s"][-1], 37):
        assert np.array_equal(EPH.interp_snapshot(float(t), sun, True),
                              d.interp_snapshot(float(t), donor_sun, True))


def test_interp_clamps_outside_the_snapshot_span():
    sun, _ = _snapshots()
    assert np.array_equal(EPH.interp_snapshot(sun["t_s"][0] - 1e6, sun, False),
                          sun["r_eci_km"][0])
    assert np.array_equal(EPH.interp_snapshot(sun["t_s"][-1] + 1e6, sun, False),
                          sun["r_eci_km"][-1])


def test_interp_is_linear_between_samples():
    sun, _ = _snapshots()
    i = 100
    t_lo, t_hi = sun["t_s"][i], sun["t_s"][i + 1]
    mid = EPH.interp_snapshot(0.5 * (t_lo + t_hi), sun, False)
    expect = 0.5 * (sun["r_eci_km"][i] + sun["r_eci_km"][i + 1])
    assert np.allclose(mid, expect, rtol=0, atol=1e-9)


def test_buggy_interpolator_still_differs_from_the_fixed_one():
    """REGRESSION GUARD. If these ever agree in the interior, someone has silently
    reverted the fix that changed a11 from -7.85e-4 (SNR 6.89) to -1.49e-4
    (SNR 1.79). Endpoints are excluded because both branches clamp identically
    outside the tabulated range, which is correct behaviour, not a fix failure.
    """
    sun, _ = _snapshots()
    lo, hi = sun["t_s"][0], sun["t_s"][-1]
    interior = np.linspace(lo + 1.0, hi - 1.0, 51)
    diffs = []
    for t in interior:
        good = EPH.interp_snapshot(float(t), sun, False)
        bad = EPH.interp_snapshot_buggy(float(t), sun, False)
        diffs.append(float(np.linalg.norm(good - bad)))
    assert min(diffs) > 0.0, "buggy twin is identical to the fixed interpolator"
    assert max(diffs) > 1e3, "expected the sign-reversed weights to be a large effect"


# --- force model ----------------------------------------------------------- #

def test_third_body_accel_matches_donor_bitwise():
    d = _donor()
    sun, moon = _snapshots()
    ds = d.load_snapshot(d.SUN_SNAPSHOT)
    dm = d.load_snapshot(d.MOON_SNAPSHOT)
    r = np.array([7000.0, 1000.0, 400.0])
    for t in (8.3e8, 9.0e8, 1.1e9):
        assert np.array_equal(EPH.third_body_accel(r, t, sun, moon),
                              d.third_body_accel(r, t, ds, dm))


def test_third_body_accel_respects_switches_and_lambda():
    sun, moon = _snapshots()
    r = np.array([7000.0, 1000.0, 400.0])
    t = 9.0e8
    both = EPH.third_body_accel(r, t, sun, moon)
    assert np.allclose(EPH.third_body_accel(r, t, sun, moon, include_moon=False)
                       + EPH.third_body_accel(r, t, sun, moon, include_sun=False),
                       both, rtol=1e-14, atol=1e-18)
    assert np.allclose(EPH.third_body_accel(r, t, sun, moon, lambda_3body=0.0),
                       np.zeros(3), atol=1e-18)
    assert np.allclose(EPH.third_body_accel(r, t, sun, moon, lambda_3body=2.0),
                       2.0 * both, rtol=1e-14)


def test_third_body_accel_is_small_next_to_point_mass():
    """Sanity scale: the lunisolar term must be a small perturbation of Kepler."""
    sun, moon = _snapshots()
    from lab_utils import MU_EARTH_KM3S2
    r = np.array([6978.137, 0.0, 0.0])
    gm_over_r2 = MU_EARTH_KM3S2 / r[0] ** 2
    third = float(np.linalg.norm(EPH.third_body_accel(r, 9.0e8, sun, moon)))
    assert third / gm_over_r2 < 1e-4


# --- initial conditions ---------------------------------------------------- #

def test_circular_ic_matches_donor():
    d = _donor()
    for inc in (97.7876, 90.0, 30.0):
        x_mine = EPH.circular_ic(6978.137, np.radians(inc), 820476800.0)
        x_theirs = d.circular_ic(6978.137, np.radians(inc), 820476800.0)
        assert np.array_equal(x_mine, x_theirs)


def test_circular_ic_has_expected_energy_and_inclination():
    from lab_utils import MU_EARTH_KM3S2, R_EARTH_KM
    for inc in (97.7876, 30.0):
        x = EPH.circular_ic(R_EARTH_KM + 600.0, np.radians(inc), 820476800.0)
        r, v = x[:3], x[3:]
        energy = np.linalg.norm(v) ** 2 / 2 - MU_EARTH_KM3S2 / np.linalg.norm(r)
        assert abs(energy - (-MU_EARTH_KM3S2 / (2 * (R_EARTH_KM + 600.0)))) < 1e-6
        h = np.cross(r, v)
        assert abs(np.degrees(np.arccos(h[2] / np.linalg.norm(h))) - inc) < 1e-9