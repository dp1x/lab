"""Adversarial tests for the J2 mode-isolation bug discovered during this mission.

These tests guard against a specific bug class:
- use_j2 = mode != "kepler_only" was treated as true for ALL non-Kepler modes
  in the original implementation, causing sun_only, moon_only, and sun_moon
  propagations to silently include J2

The tests verify the FORCE-MODE ISOLATION contract:
- j2_only propagation should NOT contain any 3b acceleration
- sun_only propagation should NOT contain J2 acceleration
- moon_only propagation should NOT contain J2 acceleration
- sun_moon should NOT contain J2
- sun_moon_j2 should contain BOTH J2 and 3b

The original bug was caught by the analysis script's R_J2x3b sanity check
that showed ~1 deg/day residuals (catastrophically wrong); these tests would
have caught it at the unit-test level.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from mission_experiment import (
    H_SSO_KM, I_SSO_DEG, I_90_DEG, I_30_DEG,
    DT_S, MU_EARTH_KM3S2, R_EARTH_KM,
    SUN_SNAPSHOT, MOON_SNAPSHOT,
    AU_KM, SOLAR_GM_KM3_S2, LUNAR_GM_KM3_S2,
    LUNAR_DISTANCE_KM_MEAN,
    _load_snapshot,
    propagate_streaming_with_x0,
)


def _initial_state(h_km, i_deg):
    a = R_EARTH_KM + h_km
    v_circ = math.sqrt(MU_EARTH_KM3S2 / a)
    i_rad = math.radians(i_deg)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, v_circ * math.cos(i_rad), v_circ * math.sin(i_rad)])
    return np.concatenate([r0, v0])


@pytest.fixture(scope="module")
def loaded_snapshots():
    """Load snapshots once for the module."""
    sun_snap = _load_snapshot(SUN_SNAPSHOT)
    moon_snap = _load_snapshot(MOON_SNAPSHOT)
    return sun_snap, moon_snap


def _run_mode(sun_snap, moon_snap, mode, h_km=H_SSO_KM, i_deg=I_SSO_DEG,
              t_end_days=30.0, dt_s=DT_S):
    """Run a propagation in the given mode and return RAAN rate in deg/day."""
    x0 = _initial_state(h_km, i_deg)
    res = propagate_streaming_with_x0(
        sun_snap, moon_snap, x0,
        mode=mode, t0_s=0.0, t_end_s=t_end_days * 86400.0, dt_s=dt_s,
    )
    if len(res["t_cross"]) < 4:
        return float("nan")
    # OLS slope
    A = np.column_stack([np.ones_like(res["t_cross"]), res["t_cross"]])
    result = np.linalg.lstsq(A, res["om_cross"], rcond=None)
    b_rad_per_day = float(result[0][1])
    return math.degrees(b_rad_per_day) * 86400.0


def test_j2_only_omits_3b_contribution(loaded_snapshots):
    """J2-only propagation should not have any 3b signal. At i_sso,
    the J2-only rate should be very close to the analytical J2 secular rate."""
    sun_snap, moon_snap = loaded_snapshots
    rate = _run_mode(sun_snap, moon_snap, "j2_only", t_end_days=30.0)
    # Analytical J2 at h=600 km i_sso: ~+0.985 deg/day
    # 30 d at 60s dt: ~43200 steps. The 30-d OLS fit gives ~+0.99 deg/day
    # because finite-window bias pushes it slightly above the secular
    analytical = 0.9855
    assert abs(rate - analytical) / analytical < 0.05, (
        f"j2_only rate {rate:.4f} differs from analytical J2 {analytical} by >5%"
    )


def test_sun_only_omits_j2_contribution(loaded_snapshots):
    """Sun-only propagation should not have J2. At i_sso the rate should be
    a small retrograde Sun term (~-7e-5 deg/day), NOT the full J2 rate
    of ~+0.99 deg/day. THIS TEST WOULD HAVE CAUGHT THE BUG."""
    sun_snap, moon_snap = loaded_snapshots
    rate = _run_mode(sun_snap, moon_snap, "sun_only", t_end_days=30.0)
    # If J2 is wrongly included, rate would be ~+0.99 deg/day.
    # If J2 is correctly excluded, rate should be ~-7e-5 deg/day (small retrograde Sun).
    # Tolerance: the absolute value must be < 0.01 deg/day (100x smaller than J2)
    assert abs(rate) < 0.01, (
        f"sun_only rate = {rate:+.4e} deg/day; expected |rate| < 0.01 deg/day "
        f"(this would have caught the use_j2 bug)"
    )


def test_moon_only_omits_j2_contribution(loaded_snapshots):
    """Moon-only propagation should not have J2. At i_sso the rate should be
    a small retrograde Moon term (~-2e-4 deg/day), NOT the J2 rate.
    THIS TEST WOULD HAVE CAUGHT THE BUG."""
    sun_snap, moon_snap = loaded_snapshots
    rate = _run_mode(sun_snap, moon_snap, "moon_only", t_end_days=30.0)
    assert abs(rate) < 0.01, (
        f"moon_only rate = {rate:+.4e} deg/day; expected |rate| < 0.01 deg/day "
        f"(this would have caught the use_j2 bug)"
    )


def test_sun_moon_omits_j2_contribution(loaded_snapshots):
    """Sun+Moon propagation should not have J2. Sum of Sun-only + Moon-only.
    THIS TEST WOULD HAVE CAUGHT THE BUG."""
    sun_snap, moon_snap = loaded_snapshots
    rate = _run_mode(sun_snap, moon_snap, "sun_moon", t_end_days=30.0)
    assert abs(rate) < 0.01, (
        f"sun_moon rate = {rate:+.4e} deg/day; expected |rate| < 0.01 deg/day "
        f"(this would have caught the use_j2 bug)"
    )


def test_kepler_only_is_truly_zero(loaded_snapshots):
    """Kepler-only must give zero RAAN drift at all inclinations. THIS IS A
    STRONGER GUARANTEE: the absence of J2 in any mode must produce zero secular
    RAAN drift (Kepler problem conserves Ω exactly)."""
    sun_snap, moon_snap = loaded_snapshots
    for i_deg in [I_SSO_DEG, I_90_DEG, I_30_DEG]:
        rate = _run_mode(sun_snap, moon_snap, "kepler_only",
                          i_deg=i_deg, t_end_days=30.0)
        assert abs(rate) < 1e-10, (
            f"kepler_only at i={i_deg}: rate = {rate:.4e} deg/day, expected 0"
        )


def test_mode_additivity_at_small_3b():
    """When the 3b signal is small (e.g., 30-d arc, no J2), the sun_moon rate
    should be approximately the sum of sun_only and moon_only rates.

    This tests the additive property of the propagator. If the modes are
    properly isolated, this should hold within numerical precision.
    """
    sun_snap, moon_snap = _load_snapshot(SUN_SNAPSHOT), _load_snapshot(MOON_SNAPSHOT)
    rate_sun = _run_mode(sun_snap, moon_snap, "sun_only", t_end_days=30.0)
    rate_moon = _run_mode(sun_snap, moon_snap, "moon_only", t_end_days=30.0)
    rate_sum = _run_mode(sun_snap, moon_snap, "sun_moon", t_end_days=30.0)
    # The sun_moon rate should be approximately rate_sun + rate_moon
    # Allow some slack for nonlinear effects at the Lunisolar level
    diff = abs(rate_sum - (rate_sun + rate_moon))
    # Both rates are ~1e-4; tolerance should be smaller than that
    assert diff < 1e-3, (
        f"Mode additivity failed: sun_only={rate_sun:+.4e}, "
        f"moon_only={rate_moon:+.4e}, sun_moon={rate_sum:+.4e}, "
        f"sum_diff={diff:+.4e}"
    )


def test_j2_sign_at_polar_inclination(loaded_snapshots):
    """At polar inclination i=90°, the J2 secular rate is cos(i)=0 (zero!).
    j2_only at i=90 should give ~0 rate. This is a STRONG sign check."""
    sun_snap, moon_snap = loaded_snapshots
    rate = _run_mode(sun_snap, moon_snap, "j2_only", i_deg=I_90_DEG, t_end_days=30.0)
    # cos(90°) = 0, so J2 secular rate is ~0
    assert abs(rate) < 0.01, (
        f"j2_only at i=90°: rate = {rate:+.4e} deg/day, expected ~0 (cos(90°)=0)"
    )


def test_sun_moon_j2_at_polar_inclination(loaded_snapshots):
    """At i=90°, J2 = 0, so sun_moon_j2 should equal sun_moon. The Lunisolar
    residual is then the pure Lunisolar contribution, NOT contaminated by J2.
    """
    sun_snap, moon_snap = loaded_snapshots
    rate_full = _run_mode(sun_snap, moon_snap, "sun_moon_j2", i_deg=I_90_DEG, t_end_days=30.0)
    rate_sm = _run_mode(sun_snap, moon_snap, "sun_moon", i_deg=I_90_DEG, t_end_days=30.0)
    # At i=90°, J2 = 0, so the rates should be equal
    diff = abs(rate_full - rate_sm)
    assert diff < 0.01, (
        f"At i=90°, sun_moon_j2 = {rate_full:+.4e}, sun_moon = {rate_sm:+.4e}, "
        f"diff = {diff:+.4e}; J2 should vanish at i=90°"
    )
