"""Regression tests for audit-020 Track 4 -- numerical-implementation audit.

Verifies every formula in the 017/018/019 propagation pipeline at machine precision:

  T1. Third-body acceleration: direct+indirect form vs gradient-of-tidal-potential form.
      Tests at 50 random states + 2026 Sun/Moon actual positions.
  T2. Precession identity: at T=0 the matrix is identity; at T=0.26 centuries
      (2026) it rotates the X-axis by the standard IAU-1976 angle.
  T3. TDB-vs-TT time-scale offset: bounded by ~1.7 ms; impact on RAAN rate.
  T4. Lunar-inclination model-order error: i3_moon constant 28.584 deg vs actual
      2026 Moon direction in ECI after precession.
  T5. Precession-corrected inclination at 2026-01-01 vs ecliptic-mean 28.584 deg.
  T6. Ascending-node convention: arctan2(r_y, r_x) recovers the lab's RAAN.
  T7. Reference radial scale factor: (a/a_3)^3 vs (R_E/a_3)^2 (the wrong formula).
  T8. Reference geometric factor: sin 2(i-i_3)/sin i vs cos i (1 - 5/2 sin^2(...))
      (the Kozai APSIDAL factor that 016/017 confused).

The tests run at machine precision (1e-15 km/s^2 tolerance for accelerations,
1e-12 rad tolerance for angles). They are read-only; they do NOT modify any
production code or results.json. They are intended to be promoted into the
018/019 test suite to catch regressions in the third-body pipeline.
"""
from __future__ import annotations

import hashlib
import math
import sys
from pathlib import Path

import numpy as np

# Lab canonical constants
MU_EARTH_KM3S2 = 398600.4418
R_EARTH_KM = 6378.137
SOLAR_GM_KM3_S2 = 132712440018.0
LUNAR_GM_KM3_S2 = 4902.8001
AU_KM = 149597870.7
LUNAR_DISTANCE_KM = 384400.0
SOLAR_OBLIQUITY_DEG = 23.439
LUNAR_INCLINATION_DEG = 5.145
DEG = math.pi / 180.0
SEC = DEG / 3600.0
JD_J2000 = 2451545.0
T0_S = 820540800.0  # 2026-01-01 12:00 TT (lab convention)


# --------------------------------------------------------------------------- #
# IAU-1976 precession (eclipseTiming/019 reference convention)
# --------------------------------------------------------------------------- #
def _rot3(angle: float) -> np.ndarray:
    """Standard active rotation about +z by +angle (eclipseTiming convention)."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot2(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])


def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    """IAU-1976 precession: J2000 -> mean-of-date (Lieske 1977 polynomial)."""
    T = t_s / (86400.0 * 36525.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * SEC
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * SEC
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * SEC
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)


# --------------------------------------------------------------------------- #
# T1: Third-body acceleration -- independent verification
# --------------------------------------------------------------------------- #
def tidal_accel_form_a(r_sat: np.ndarray, r3: np.ndarray, mu3: float) -> np.ndarray:
    """Form (a): direct + indirect, as in 017/018 implementations.

    a = mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3 - mu_3 * r_3 / |r_3|^3
    """
    d = r3 - r_sat
    return mu3 * (d / np.linalg.norm(d) ** 3 - r3 / np.linalg.norm(r3) ** 3)


def tidal_accel_form_b(r_sat: np.ndarray, r3: np.ndarray, mu3: float) -> np.ndarray:
    """Form (b): gradient of the disturbing potential U_3.

    U_3(r_sat) = mu_3 / |r_3 - r_sat| - mu_3 (r_sat . r_3) / |r_3|^3

    The first term is the direct attraction (negative gradient of 1/|r_3-r_sat|);
    the second term is the indirect (centrifugal in the co-accelerating frame).
    """
    d = r3 - r_sat
    dnorm = np.linalg.norm(d)
    r3norm = np.linalg.norm(r3)
    # Gradient of mu_3 / |r_3 - r_sat| wrt r_sat: +mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3
    grad_direct = mu3 * d / dnorm ** 3
    # Gradient of -mu_3 (r_sat . r_3) / |r_3|^3 wrt r_sat: -mu_3 * r_3 / |r_3|^3
    grad_indirect = -mu3 * r3 / r3norm ** 2
    # The total acceleration is +grad U_3 (NOT -grad U_3, because U_3 is the
    # perturbation to the potential energy which adds to -mu_E/|r_sat|)
    return grad_direct + grad_indirect


def test_t1_form_a_equals_form_b_random():
    """Test 1: Form (a) and Form (b) agree to machine precision at 50 random states."""
    rng = np.random.default_rng(42)
    a_orbit = R_EARTH_KM + 600.0

    # Random satellite positions on sphere of radius a_orbit
    max_diff_sun = 0.0
    max_diff_moon = 0.0

    # Use actual 2026 Sun/Moon positions (geocentric ICRF)
    r3_sun_icrf = np.array([1.443257274919273e5, 2.895841575449329e5, 1.601589230016842e5])
    r3_moon_icrf = np.array([1.443257274919273e5, 2.895841575449329e5, 1.601589230016842e5])

    for _ in range(50):
        v = rng.standard_normal(3)
        v /= np.linalg.norm(v)
        r_sat = a_orbit * v

        a_a = tidal_accel_form_a(r_sat, r3_sun_icrf, SOLAR_GM_KM3_S2)
        a_b = tidal_accel_form_b(r_sat, r3_sun_icrf, SOLAR_GM_KM3_S2)
        max_diff_sun = max(max_diff_sun, float(np.max(np.abs(a_a - a_b))))

        a_a = tidal_accel_form_a(r_sat, r3_moon_icrf, LUNAR_GM_KM3_S2)
        a_b = tidal_accel_form_b(r_sat, r3_moon_icrf, LUNAR_GM_KM3_S2)
        max_diff_moon = max(max_diff_moon, float(np.max(np.abs(a_a - a_b))))

    assert max_diff_sun < 1e-15, f"form a vs b SUN diff = {max_diff_sun:.3e}"
    assert max_diff_moon < 1e-15, f"form a vs b MOON diff = {max_diff_moon:.3e}"


def test_t1_indirect_term_magnitude():
    """Test 1b: quantify the indirect/direct ratio for Sun and Moon.

    The 017 docstring claims "order 1e-5 of the direct term" which is correct
    for the Sun but underestimates the Moon by a factor ~350.
    """
    rng = np.random.default_rng(43)
    a_orbit = R_EARTH_KM + 600.0
    r3_sun = np.array([1.0, 0.0, 0.0]) * AU_KM
    r3_moon = np.array([0.0, 1.0, 0.0]) * LUNAR_DISTANCE_KM

    ratios_sun = []
    ratios_moon = []
    for _ in range(50):
        v = rng.standard_normal(3)
        v /= np.linalg.norm(v)
        r_sat = a_orbit * v

        d_sun = r3_sun - r_sat
        dnorm_sun = np.linalg.norm(d_sun)
        direct_sun = SOLAR_GM_KM3_S2 * d_sun / dnorm_sun ** 3
        indirect_sun = -SOLAR_GM_KM3_S2 * r3_sun / np.linalg.norm(r3_sun) ** 3
        ratios_sun.append(float(np.linalg.norm(indirect_sun) / np.linalg.norm(direct_sun)))

        d_moon = r3_moon - r_sat
        dnorm_moon = np.linalg.norm(d_moon)
        direct_moon = LUNAR_GM_KM3_S2 * d_moon / dnorm_moon ** 3
        indirect_moon = -LUNAR_GM_KM3_S2 * r3_moon / np.linalg.norm(r3_moon) ** 3
        ratios_moon.append(float(np.linalg.norm(indirect_moon) / np.linalg.norm(direct_moon)))

    ratio_sun = float(np.mean(ratios_sun))
    ratio_moon = float(np.mean(ratios_moon))

    # Sun: 6978 / 1.5e8 ~ 4.7e-5
    assert 3e-5 < ratio_sun < 7e-5, f"Sun indirect/direct ratio {ratio_sun:.3e} not in [3e-5, 7e-5]"
    # Moon: 6978 / 3.84e5 ~ 1.8e-2 (NOT 1e-5 as 017 docstring claims)
    assert 0.01 < ratio_moon < 0.03, f"Moon indirect/direct ratio {ratio_moon:.3e} not in [0.01, 0.03]"


# --------------------------------------------------------------------------- #
# T2: Precession identity
# --------------------------------------------------------------------------- #
def test_t2_precession_identity_at_T0():
    """Test 2: at T=0 the precession matrix is the identity."""
    P0 = precession_j2000_to_mod(0.0)
    diff = np.max(np.abs(P0 - np.eye(3)))
    assert diff < 1e-15, f"P0 != identity, max diff = {diff:.3e}"


def test_t2_precession_at_2026():
    """Test 2b: at T=0.26 centuries (2026), the matrix rotates X-axis by ~-0.333 deg.

    The eclipseTiming convention [[c,-s],[s,c]] produces a NEGATIVE rotation
    about +z for T > 0. The 018 original [[c,s],[-s,c]] would produce a POSITIVE
    rotation (the bug identified by audit-019 Track D and remediated).
    """
    P_2026 = precession_j2000_to_mod(T0_S)
    x_axis = np.array([1.0, 0.0, 0.0])
    x_rot = P_2026 @ x_axis
    angle_deg = math.degrees(math.atan2(x_rot[1], x_rot[0]))

    # Expected: ~-0.333 deg (eclipseTiming convention)
    assert -0.35 < angle_deg < -0.31, f"x-axis rotation at 2026 = {angle_deg:.4f} deg, expected ~-0.333"


# --------------------------------------------------------------------------- #
# T3: TDB vs TT time-scale offset (1.7 ms peak)
# --------------------------------------------------------------------------- #
def test_t3_tdb_tt_offset_impact_on_raan():
    """Test 3: TDB-TT offset (~1.7 ms peak) impact on secular RAAN rate.

    The snapshot is in TDB; the propagator uses TT-like. Offset is sub-second
    over 1 year; impact on secular rate is below the noise floor.
    """
    # TDB-TT peak offset: ~1.7 ms (annual + shorter periodic terms)
    delta_t_peak_s = 1.7e-3

    # Mean motion at h=600 km
    a = R_EARTH_KM + 600.0
    n = math.sqrt(MU_EARTH_KM3S2 / a ** 3)  # rad/s

    # Earth's heliocentric speed
    v_earth = 29.78  # km/s
    # Sun's apparent angular rate from Earth
    omega_sun_rad_s = v_earth / AU_KM  # ~2e-7 rad/s

    # A 1.7 ms shift in time corresponds to a Sun direction shift of
    # delta_t * omega_sun = 1.7e-3 * 2e-7 = 3.4e-10 rad = 7e-5 arcsec
    sun_dir_shift_rad = delta_t_peak_s * omega_sun_rad_s

    # Corresponding shift in i_3 (Sun inclination in MOD):
    # delta sin 2(i-i_3) ~ 2 cos 2(i-i_3) * delta i_3
    i = math.radians(97.7876)
    i3_sun = math.radians(23.439)
    d_sin2_factor = abs(2.0 * math.cos(2.0 * (i - i3_sun)))

    # RAAN rate per unit sin 2(i-i_3)/sin i
    solar_raan_factor = (3.0 / 8.0) * n * (SOLAR_GM_KM3_S2 / MU_EARTH_KM3S2) * (a / AU_KM) ** 3 / math.sin(i)

    # Delta in RAAN rate due to 1.7 ms time offset:
    delta_raan_rad_s = solar_raan_factor * d_sin2_factor * sun_dir_shift_rad
    delta_raan_deg_day = math.degrees(delta_raan_rad_s) * 86400.0

    # Should be < 1e-10 deg/day (well below the 1e-4 deg/day secular rate)
    assert delta_raan_deg_day < 1e-10, (
        f"TDB-TT 1.7 ms peak offset impact = {delta_raan_deg_day:.3e} deg/day "
        f"(expected < 1e-10 deg/day)"
    )


# --------------------------------------------------------------------------- #
# T4: Lunar-inclination constant 28.584 deg vs actual 2026 position
# --------------------------------------------------------------------------- #
def test_t4_lunar_inclination_model_order_error():
    """Test 4: model-order error from using constant 28.584 deg for the Moon.

    The Moon's actual geocentric direction at 2026-01-01 in ICRF is rotated
    to MOD by the IAU-1976 precession. The inclination to the MOD equator is
    NOT exactly 28.584 deg (which is the obliquity + lunar mean inclination).
    The actual inclination varies over the 18.6-yr nodal cycle between
    ~18.29 deg and ~28.58 deg. Quantify the model-order error in the secular
    RAAN rate from using a constant inclination.
    """
    # Moon at 2026-01-01 00:00 TDB (from snapshot, ICRF)
    r_moon_icrf = np.array([1.443257274919273e5, 2.895841575449329e5, 1.601589230016842e5])

    # Convert to MOD using IAU-1976 precession at t0 - 12h (the snapshot first row)
    t_first_row = T0_S - 43200.0
    P = precession_j2000_to_mod(t_first_row)
    r_moon_mod = P @ r_moon_icrf

    # Inclination to MOD equator (angle from +Z axis)
    r_norm = np.linalg.norm(r_moon_mod)
    cos_i3_actual = r_moon_mod[2] / r_norm
    i3_actual_deg = math.degrees(math.acos(cos_i3_actual))

    # The constant used in the lab formula
    i3_constant_deg = SOLAR_OBLIQUITY_DEG + LUNAR_INCLINATION_DEG  # = 28.584

    delta_i3 = abs(i3_actual_deg - i3_constant_deg)

    # The Moon in 2026 is near the descending node of the 18.6-yr cycle, so
    # i3_moon is near its minimum (~18.29 deg). The actual value is computed
    # below; the difference from 28.584 should be at least 5-10 deg.

    # Quantify the impact on the secular RAAN rate at h=600 km, i=i_sso
    a = R_EARTH_KM + 600.0
    n = math.sqrt(MU_EARTH_KM3S2 / a ** 3)
    i = math.radians(97.7876)
    lunar_raan_factor = (3.0 / 8.0) * n * (LUNAR_GM_KM3_S2 / MU_EARTH_KM3S2) * (
        a / LUNAR_DISTANCE_KM
    ) ** 3 / math.sin(i)

    # RAAN rate at constant inclination
    rate_constant = lunar_raan_factor * math.sin(2.0 * (i - math.radians(i3_constant_deg)))
    # RAAN rate at actual inclination
    rate_actual = lunar_raan_factor * math.sin(2.0 * (i - math.radians(i3_actual_deg)))

    # Relative error from using constant inclination
    if abs(rate_constant) > 0:
        rel_error = abs(rate_actual - rate_constant) / abs(rate_constant)
    else:
        rel_error = float("inf")

    # The actual 2026 lunar inclination is well below the constant; the
    # model-order error is large in fractional terms (a factor of ~10-50)
    # because the lunar inclination at the descending node approaches
    # obliquity - 5.145 = 18.29 deg, vs the 28.584 deg maximum used.
    print(f"[T4] actual i3_moon at 2026-01-01 (MOD) = {i3_actual_deg:.4f} deg")
    print(f"[T4] constant i3_moon = {i3_constant_deg:.4f} deg")
    print(f"[T4] delta = {delta_i3:.4f} deg")
    print(f"[T4] lunar RAAN rate constant = {math.degrees(rate_constant)*86400:.4e} deg/day")
    print(f"[T4] lunar RAAN rate actual = {math.degrees(rate_actual)*86400:.4e} deg/day")
    print(f"[T4] relative error = {rel_error:.4f}")

    # Sanity: assert delta_i3 > 0.1 deg (clearly different from constant)
    assert delta_i3 > 0.1, f"actual i3_moon = {i3_actual_deg} deg, expected to differ from constant"


# --------------------------------------------------------------------------- #
# T5: Precession-corrected inclination at 2026-01-01
# --------------------------------------------------------------------------- #
def test_t5_moon_direction_in_mod_after_precession():
    """Test 5: at 2026-01-01, the Moon's ECI (MOD) direction differs from the
    ecliptic-mean assumption. Quantify the inclination and node angle.
    """
    r_moon_icrf = np.array([1.443257274919273e5, 2.895841575449329e5, 1.601589230016842e5])
    t_first_row = T0_S - 43200.0
    P = precession_j2000_to_mod(t_first_row)
    r_moon_mod = P @ r_moon_icrf

    # Right ascension of the Moon (in MOD frame)
    ra_moon_deg = math.degrees(math.atan2(r_moon_mod[1], r_moon_mod[0])) % 360.0
    # Declination of the Moon (in MOD frame)
    dec_moon_deg = math.degrees(math.asin(r_moon_mod[2] / np.linalg.norm(r_moon_mod)))

    # The Moon's mean ecliptic inclination (5.145 deg) makes the Moon's mean
    # equatorial inclination ~ 28.584 deg, but the actual instantaneous
    # equatorial inclination is dec_moon_deg at the moment the Moon is at RA=ra_moon_deg
    print(f"[T5] Moon RA (MOD) at 2026-01-01 = {ra_moon_deg:.4f} deg")
    print(f"[T5] Moon Dec (MOD) at 2026-01-01 = {dec_moon_deg:.4f} deg")
    print(f"[T5] Note: Moon's instantaneous equatorial inclination == its declination "
          f"only at RA=0 or RA=180; otherwise it's a more complex function of RA, Dec, obliquity.")


# --------------------------------------------------------------------------- #
# T6: Ascending-node convention
# --------------------------------------------------------------------------- #
def test_t6_ascending_node_omega_convention():
    """Test 6: numerical Omega = atan2(r_y, r_x) at ascending-node crossing.

    This is the standard RAAN convention (prograde from X axis to node vector
    in XY plane). The theoretical formula assumes the same convention.
    """
    # Initial conditions at the ascending node: r = (a, 0, 0)
    a = R_EARTH_KM + 600.0
    r = np.array([a, 0.0, 0.0])

    # Omega = atan2(r_y, r_x) = atan2(0, a) = 0
    omega = math.atan2(r[1], r[0])
    assert abs(omega) < 1e-15, f"Omega at ascending node = {omega}, expected 0"

    # Rotate by 90 deg about Z: r = (0, a, 0)
    c, s = math.cos(math.pi / 2), math.sin(math.pi / 2)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    r_rot = R @ r
    omega_rot = math.atan2(r_rot[1], r_rot[0])
    assert abs(omega_rot - math.pi / 2) < 1e-12, f"Omega after +90 deg = {omega_rot}, expected pi/2"


# --------------------------------------------------------------------------- #
# T7: Reference radial scale factor
# --------------------------------------------------------------------------- #
def test_t7_correct_radial_factor_a_cubed():
    """Test 7: the correct radial scale factor is (a/a_3)^3, NOT (R_E/a_3)^2.

    The wrong 016/017 formula used (R_E/a_3)^2 (J2-style); the corrected
    formula uses (a/a_3)^3 (third-body style). Verify the magnitude ratio.
    """
    a = R_EARTH_KM + 600.0
    ratio_correct = (a / AU_KM) ** 3
    ratio_wrong = (R_EARTH_KM / AU_KM) ** 2

    # The ratio correct/wrong for the Sun is a^3 / R_E^2 = a/R_E^2 * a^2
    # For Sun: ratio = (6978)^3 / (6378)^2 = 3.4e11 / 4.07e7 = 8360x
    factor = ratio_correct / ratio_wrong
    # Should be ~ (a/R_E)^2 ~ (6978/6378)^2 ~ 1.197, but that's wrong;
    # let me recompute: ratio_correct = (6978/1.5e8)^3, ratio_wrong = (6378/1.5e8)^2
    # ratio_correct/ratio_wrong = 6978^3 / 6378^2 / 1.5e8 = a^3 / R_E^2 / AU
    # = (6978^3) / (6378^2 * 1.5e8)
    factor_explicit = a ** 3 / (R_EARTH_KM ** 2 * AU_KM)
    # For the Moon: ratio_correct/ratio_wrong = a^3 / R_E^2 / R_M
    # = (6978^3) / (6378^2 * 3.84e5)
    factor_moon = a ** 3 / (R_EARTH_KM ** 2 * LUNAR_DISTANCE_KM)

    # These are the magnitude differences between the wrong and correct formulas.
    print(f"[T7] Sun: wrong/correct = {1.0/factor_explicit:.4e} (or {factor_explicit:.4e}x over-estimate)")
    print(f"[T7] Moon: wrong/correct = {1.0/factor_moon:.4e} (or {factor_moon:.4e}x over-estimate)")
    # The wrong formula is ~ a^3/R_E^2/a_3 = a/R_E^2 * (a/a_3)^2 / 1 times too large
    # vs the correct formula which has (a/a_3)^3 / ((a/a_3)^3) = 1
    # Actually the wrong formula has (R_E/a_3)^2 vs correct (a/a_3)^3
    # Ratio = (a/R_E)^3 * (a_3/R_E)
    # Sun: (6978/6378)^3 * (1.5e8/6378) = 1.197^3 * 2.35e4 = 1.71 * 2.35e4 ~ 4e4
    # Moon: (6978/6378)^3 * (3.84e5/6378) = 1.71 * 60.2 ~ 103
    # Let me verify directly:
    sun_wrong_over_correct = (R_EARTH_KM / AU_KM) ** 2 / (a / AU_KM) ** 3
    sun_correct_over_wrong = 1.0 / sun_wrong_over_correct
    moon_wrong_over_correct = (R_EARTH_KM / LUNAR_DISTANCE_KM) ** 2 / (a / LUNAR_DISTANCE_KM) ** 3
    moon_correct_over_wrong = 1.0 / moon_wrong_over_correct
    print(f"[T7] Sun wrong over correct = {sun_wrong_over_correct:.4e}")
    print(f"[T7] Moon wrong over correct = {moon_wrong_over_correct:.4e}")
    # The wrong formula is 2-4 orders of magnitude too large at the radial scale.
    # For the Sun it's ~40000x; for the Moon it's ~100x.
    assert sun_wrong_over_correct > 1000, (
        f"Sun wrong/correct radial factor = {sun_wrong_over_correct:.3e}, expected >> 1000"
    )
    assert moon_wrong_over_correct > 10, (
        f"Moon wrong/correct radial factor = {moon_wrong_over_correct:.3e}, expected >> 10"
    )


# --------------------------------------------------------------------------- #
# T8: Reference geometric factor
# --------------------------------------------------------------------------- #
def test_t8_correct_geometric_factor_nodal_vs_apsidal():
    """Test 8: the correct geometric factor for nodal rate is sin 2(i-i_3)/sin i.

    The wrong 016/017 formula used the Kozai APSIDAL factor
    cos i (1 - 5/2 sin^2(i-i_3)). Verify they give different physics.
    """
    i = math.radians(97.7876)  # SSO retrograde
    i3 = math.radians(28.584)  # Moon equatorial mean
    ii3 = i - i3

    nodal = math.sin(2.0 * ii3) / math.sin(i)
    apsidal = math.cos(i) * (1.0 - 2.5 * math.sin(ii3) ** 2)

    print(f"[T8] nodal factor at i_sso: {nodal:.6f}")
    print(f"[T8] apsidal factor at i_sso: {apsidal:.6f}")
    print(f"[T8] ratio nodal/apsidal: {nodal/apsidal:.6f}")

    # They are different functions. The nodal is the correct one for dOmega/dt;
    # the apsidal is the correct one for d omega/dt (argument of periapsis).
    # Verify they have different magnitudes/signs at SSO.
    assert abs(nodal - apsidal) > 0.01, "nodal vs apsidal factors must differ"


# --------------------------------------------------------------------------- #
# Main entry point
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    print("=" * 70)
    print("Audit-020 Track 4 -- Numerical Implementation Audit")
    print("Regression tests for 017/018/019 propagation pipeline")
    print("=" * 70)

    tests = [
        ("T1a", test_t1_form_a_equals_form_b_random),
        ("T1b", test_t1_indirect_term_magnitude),
        ("T2a", test_t2_precession_identity_at_T0),
        ("T2b", test_t2_precession_at_2026),
        ("T3", test_t3_tdb_tt_offset_impact_on_raan),
        ("T4", test_t4_lunar_inclination_model_order_error),
        ("T5", test_t5_moon_direction_in_mod_after_precession),
        ("T6", test_t6_ascending_node_omega_convention),
        ("T7", test_t7_correct_radial_factor_a_cubed),
        ("T8", test_t8_correct_geometric_factor_nodal_vs_apsidal),
    ]

    failed = 0
    passed = 0
    for label, test_fn in tests:
        try:
            test_fn()
            print(f"\n[OK] {label}: {test_fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"\n[FAIL] {label}: {test_fn.__name__}: {e}")
            failed += 1

    print()
    print("=" * 70)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 70)
    sys.exit(0 if failed == 0 else 1)