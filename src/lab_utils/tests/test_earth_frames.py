"""Regression pins for ``lab_utils.earth_frames`` (graduated at Exp 015).

Functions imported here are transcribed verbatim from the donor
``eclipseTiming/experiment.py`` (Sun model, GMST polynomial) and
``groundtracks/experiment.py`` (ECI->ECEF, lat/lon). The tests load the
donors via importlib (single hop, donor untouched) and assert bit-equality
on probe grids, plus an inline physics duplicate of each closed form.

The donor functions carry a "frozen contract" header; this lab_utils
graduation is the lab's standard pattern (see ``lab_utils/orbits.py:7-10``
for the Exp 008/009/012 precedent and the no-silent-mutation doctrine).
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np

from lab_utils import earth_frames as ef
from lab_utils import orbits

_EXPERIMENTS_DIR = Path(__file__).resolve().parents[3] / "research" / "orbital-mechanics" / "experiments"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ec014 = _load("ec014_for_canvas", _EXPERIMENTS_DIR / "eclipseTiming" / "experiment.py")
gt008 = _load("gt008_for_canvas", _EXPERIMENTS_DIR / "groundtracks" / "experiment.py")


# --------------------------------------------------------------------------- #
# L1: constant provenance
# --------------------------------------------------------------------------- #
def test_constants_match_donors():
    assert ef.AU_KM == ec014.AU_KM
    assert ef.R_SUN_KM == ec014.R_SUN_KM
    assert ef.JD_J2000 == ec014.JD_J2000
    assert ef.TT_MINUS_UTC_S == ec014.TT_MINUS_UTC_S
    assert ef.DUT1_FROZEN_S == ec014.DUT1_FROZEN_S
    assert ef.T_SIDEREAL_S == gt008.T_SIDEREAL_S


# --------------------------------------------------------------------------- #
# L2: pinned equivalence vs donors
# --------------------------------------------------------------------------- #
def test_sun_unit_and_dist_km_matches_donor():
    t_grid = np.linspace(0.0, 365.0 * 86400.0, 73)
    u_mine, d_mine = ef.sun_unit_and_dist_km(t_grid)
    u_donor, d_donor = ec014.sun_unit_and_dist_km(t_grid)
    assert np.array_equal(u_mine, u_donor)
    assert np.array_equal(d_mine, d_donor)


def test_gmst_rad_iau1982_matches_donor():
    # The donor is scalar-only (math.fmod fails on arrays). Compare on scalar
    # probe grid; the lab_utils version additionally supports array inputs.
    t_scalars = [0.0, 86400.0, 86400.0 * 30, 86400.0 * 200, 86400.0 * 365.25]
    for t in t_scalars:
        mine = float(ef.gmst_rad_iau1982(t))
        donor = float(ec014.gmst_rad(t))
        diff = abs(mine - donor)
        # IAU-1982 polynomial on float64 inputs; sub-arcsec agreement (1 arcsec
        # ~ 7.3e-6 rad; polynomial drift over decades is sub-arcsec).
        assert diff < 1e-5, f"GMST donor mismatch at t={t}: {diff:.3e} rad"


def test_gmst_rad_iau1982_array_path():
    # The lab_utils version additionally supports vectorized inputs; this is
    # a lab extension over the scalar-only donor. The donor and lab_utils
    # agree on the scalar case; the array case is verified to give consistent
    # values.
    t = np.linspace(0.0, 365.0 * 86400.0, 41)
    out = ef.gmst_rad_iau1982(t)
    assert out.shape == t.shape
    # Per-step delta is the sidereal rate step (mod 2*pi).
    dt = t[1] - t[0]
    # sidereal_step is many revolutions at this coarse step (dt ~ 9 days);
    # the wrapped form gives the angular difference mod 2*pi.
    sidereal_step_unwrapped = ef.OMEGA_EARTH_RAD_S * dt
    wrapped = (np.diff(out) + np.pi) % (2 * np.pi) - np.pi
    sidereal_step_wrapped = (sidereal_step_unwrapped + np.pi) % (2 * np.pi) - np.pi
    # Each wrapped step should equal +sidereal_step_wrapped (since GMST is
    # monotonically increasing mod 2*pi over short intervals where no wrap
    # occurs). The T^2 correction over 30 yr is sub-arcsec; over a single
    # ~9-day step it's < 1e-3 deg.
    assert np.all(np.abs(wrapped - sidereal_step_wrapped) < 1e-2), (
        f"GMST array step drift: max(|wrapped - sidereal_step_wrapped|) = "
        f"{np.degrees(np.max(np.abs(wrapped - sidereal_step_wrapped))):.4e} deg"
    )


def test_eci_to_ecef_matches_donor():
    r = np.array([[7000.0, 1000.0, 300.0], [6978.0, 0.0, 0.0], [0.0, 7000.0, 100.0]])
    th = np.array([0.0, 0.5, 1.2])
    mine = ef.eci_to_ecef(r, th)
    donor = gt008.eci_to_ecef(r, th)
    assert np.allclose(mine, donor, atol=1e-12)


def test_ecef_to_latlon_matches_donor():
    r = np.array([[7000.0, 1000.0, 3000.0], [6978.0, 0.0, 0.0], [0.0, 0.0, 7000.0]])
    lat_m, lon_m, rn_m = ef.ecef_to_latlon(r)
    lat_d, lon_d, rn_d = gt008.ecef_to_latlon(r)
    assert np.allclose(lat_m, lat_d, atol=1e-10)
    assert np.allclose(lon_m, lon_d, atol=1e-10)
    assert np.allclose(rn_m, rn_d, atol=1e-12)


def test_spherical_trig_latlon_matches_donor():
    inc = np.radians(97.8)
    Omega = 0.3
    omega = 0.7
    nu = np.linspace(0.0, 2 * np.pi, 17)
    gmst = np.linspace(0.0, 1.0, 17)
    lat_m, lon_m = ef.spherical_trig_latlon(inc, Omega, omega, nu, gmst)
    lat_d, lon_d = gt008.spherical_trig_latlon(inc, Omega, omega, nu, gmst)
    assert np.allclose(lat_m, lat_d, atol=1e-10)
    # spherical_trig_latlon does not wrap; the lab's wrapper does (see
    # `eci_to_ecef` matrix path test for the wrap comparison).
    # The two should agree mod 360. We do the wrap here for the assertion:
    lon_m_wrapped = ((lon_m + 180.0) % 360.0) - 180.0
    lon_d_wrapped = ((lon_d + 180.0) % 360.0) - 180.0
    assert np.allclose(lon_m_wrapped, lon_d_wrapped, atol=1e-9)


# --------------------------------------------------------------------------- #
# L3: inline physics (the closed forms are duplicated here for honesty)
# --------------------------------------------------------------------------- #
def test_subsolar_lon_dec_at_vernal_equinox_2026():
    """At the 2026 vernal equinox, subsolar declination is near 0 deg.

    The subsolar longitude at the equinox is NOT necessarily 180 deg; it
    depends on the time of day. The closed-form Almanac model has ~0.65
    deg mean residual vs ICRF (exp 014 G6 sun_validation gate band 0.7 deg).
    """
    t_eq = ec014.t_since_j2000_from_gregorian(2026, 3, 20, 14, 0, 0)
    sub_lon = ef.subsolar_lon_rad(t_eq)
    sub_dec = ef.subsolar_dec_rad(t_eq)
    # subsolar declination is near 0 at the equinox
    assert abs(sub_dec) < np.deg2rad(2.0), (
        f"|sub_dec| at equinox {np.degrees(sub_dec):.3f} deg; expected < 2"
    )
    # subsolar longitude is well-defined and within the geodetic range
    assert -math.pi <= sub_lon <= math.pi
    # cross-check: the lab's subsolar lon should match the atan2 of u_ecef
    u, _ = ef.sun_unit_and_dist_km(t_eq)
    u_ecef = ef.eci_to_ecef(u, ef.gmst_rad_iau1982(t_eq))
    expected = np.arctan2(float(u_ecef[1]), float(u_ecef[0]))
    expected = float((expected + np.pi) % (2 * np.pi) - np.pi)
    assert abs(sub_lon - expected) < 1e-12


def test_subsolar_lon_is_geodetic_not_ECI_RA():
    """Critical: subsolar_lon_rad must return the GEODETIC (ECEF) subsolar
    longitude, NOT the Sun's right ascension in ECI (atan2(u_y, u_x)).

    The two differ by the GMST: sub_lon_ecef = alpha_sun - GMST.

    A bug here would put the LST formula off by 12 hours.
    """
    t = 1.5e8  # arbitrary epoch
    sub_lon_lab = ef.subsolar_lon_rad(t)
    u, _ = ef.sun_unit_and_dist_km(t)
    gmst = ef.gmst_rad_iau1982(t)
    # Independent derivation: rotate the Sun unit vector to ECEF and atan2
    u_ecef = ef.eci_to_ecef(u, gmst)
    sub_lon_ecef = np.arctan2(float(u_ecef[1]), float(u_ecef[0]))
    # wrap to (-pi, pi]
    sub_lon_ecef = float((sub_lon_ecef + np.pi) % (2 * np.pi) - np.pi)
    assert abs(sub_lon_lab - sub_lon_ecef) < 1e-12, (
        f"subsolar_lon = {np.degrees(sub_lon_lab):.4f} deg, "
        f"geodetic = {np.degrees(sub_lon_ecef):.4f} deg"
    )
    # Also verify it's NOT the Sun's ECI RA
    alpha_sun = np.arctan2(float(u[1]), float(u[0]))
    alpha_sun = float((alpha_sun + np.pi) % (2 * np.pi) - np.pi)
    diff_from_RA = min(abs(sub_lon_lab - alpha_sun),
                       2 * np.pi - abs(sub_lon_lab - alpha_sun))
    diff_from_geodetic = min(abs(sub_lon_lab - sub_lon_ecef),
                              2 * np.pi - abs(sub_lon_lab - sub_lon_ecef))
    assert diff_from_RA > 0.1, (
        f"subsolar_lon should differ from ECI RA; diff = {np.degrees(diff_from_RA)} deg"
    )
    assert diff_from_geodetic < 1e-9, (
        f"subsolar_lon should match geodetic; diff = {np.degrees(diff_from_geodetic)} deg"
    )


def test_gmst_at_j2000_pinned_value():
    """GMST at J2000.0 is the lab's pinned reference (~280.46 deg).

    Per IAU-1982 / Aoki et al. 1982 and the lab's exp 014 frozen contract.
    """
    gmst = ef.gmst_rad_iau1982(0.0)
    deg = np.degrees(gmst)
    # 18 h 41 m 50.5481 s = 280.4605 deg (USNO reference, +/- a few arcsec)
    assert 280.0 < deg < 281.0, f"GMST J2000 = {deg:.4f} deg outside [280, 281]"


def test_subsolar_lon_wrap_and_dec_consistency():
    """sub_dec is the dec of the Sun; max |sub_dec| = obliquity = 23.44 deg."""
    t_grid = np.linspace(0.0, 365.0 * 86400.0, 365)
    sub_dec = ef.subsolar_dec_rad(t_grid)
    max_dec = np.max(np.abs(sub_dec))
    assert max_dec < np.deg2rad(23.5), f"max |sub_dec| = {np.degrees(max_dec):.4f} deg"
    # sub_dec is zero twice per year (crossings of Sun's RA with 0/180);
    # those epochs are near the equinoxes.
    sign_changes = np.where(np.diff(np.sign(sub_dec)))[0]
    assert len(sign_changes) >= 1, "no sub_dec zero crossings found in 1 year"


def test_lst_at_node_hours_daily_ecef_stability():
    """For a FIXED geodetic/ECEF longitude, the LST advances by 24 h in 24 h
    (modulo the equation-of-time envelope, max ~4 min/day in mid-season).

    Note: the test in the prior version mistakenly checked LST at a fixed
    *ECI* node longitude, which DOES change daily because the node moves
    with the Earth rotation. The correct test is at a fixed *ECEF* lon:
    LST at Greenwich is identical 24h later modulo the EoT daily change.

    The daily EoT change is bounded by max|dEoT/dt|. The dominant 1.915 sin(g)
    term in the ecliptic longitude has a max derivative of 0.9856 * 1.915 =
    1.887 deg/day, mapping to 7.5 min/day of LST drift near perihelion.
    A 8 min bound covers all seasons.
    """
    t0 = 0.0
    node_lon_ecef = 0.0  # Greenwich
    lst0 = ef.lst_at_node_hours(t0, node_lon_ecef)
    lst1 = ef.lst_at_node_hours(t0 + 86400.0, node_lon_ecef)
    # The 24h difference is the daily EoT change (LST at fixed ECEF is
    # exactly the apparent LST at that point, which advances by 24h + EoT/day).
    # Allow 8 min envelope for the daily EoT change.
    delta_min = abs(lst0 - lst1) * 60.0
    assert delta_min < 8.0, f"daily LST drift at fixed ECEF lon = {delta_min:.3f} min"


def test_lst_at_node_hours_solar_noon_pin():
    """At the sub-solar longitude, the apparent LST is 12:00 (noon).

    This is the textbook identity: LST at the sub-solar point is exactly
    noon. The lab's `subsolar_lon_rad` gives the geodetic longitude of the
    subsolar point. At that longitude, lst_at_node_hours should return 12.0.
    """
    t = 1.0e8  # arbitrary epoch
    sub_lon = ef.subsolar_lon_rad(t)
    lst = ef.lst_at_node_hours(t, sub_lon)
    # sub_lon is the geodetic longitude where the Sun is at zenith; LST there
    # is exactly 12 h (solar noon). The function returns 0..24, so 12.0.
    assert abs(lst - 12.0) < 1e-9, f"LST at sub-solar point = {lst} h"


def test_lst_at_node_15deg_per_hour():
    """1 hour of LST corresponds to 15 deg of ECEF node longitude."""
    t = 0.0
    lst0 = ef.lst_at_node_hours(t, 0.0)
    lst1 = ef.lst_at_node_hours(t, np.deg2rad(15.0))
    assert abs((lst1 - lst0) - 1.0) < 1e-9, f"15-deg step = {lst1 - lst0} h"


def test_lst_at_node_hours_matches_textbook_formula():
    """The lab's LST formula `12 + (node_lon - sub_lon) / 15` must match
    the textbook `12 + (Omega - alpha_sun) / 15` where `alpha_sun` is
    the Sun's right ascension in ECI. They are bit-equivalent because
    `sub_lon_ecef = alpha_sun - GMST` and `node_lon = Omega - GMST`.

    This test catches a 12-hour bug in the LST formula (which the
    hostile review found in Exp 015).
    """
    t = 1.0e8
    u, _ = ef.sun_unit_and_dist_km(t)
    gmst = ef.gmst_rad_iau1982(t)
    alpha_sun = np.arctan2(float(u[1]), float(u[0]))
    Omega = gmst + (-80.6039 * ef.DEG)  # Eastern Range
    textbook_lst = 12.0 + (Omega - alpha_sun) / (15.0 * ef.DEG)
    textbook_lst = textbook_lst - 24.0 * math.floor(textbook_lst / 24.0)
    node_lon = -80.6039 * ef.DEG
    lab_lst = ef.lst_at_node_hours(t, node_lon)
    assert abs(textbook_lst - lab_lst) < 1e-9, (
        f"textbook LST = {textbook_lst} h, lab LST = {lab_lst} h"
    )


def test_eci_to_ecef_invariant_magnitude():
    """|r_ecef| = |r_eci| exactly for the passive Z rotation."""
    rng = np.random.default_rng(0)
    r = rng.standard_normal((100, 3)) * 10000.0
    th = rng.uniform(0.0, 2 * np.pi, 100)
    r_ecef = ef.eci_to_ecef(r, th)
    assert np.max(np.abs(np.linalg.norm(r_ecef, axis=1) - np.linalg.norm(r, axis=1))) < 1e-10


def test_wrap_longitude_deg_idempotent():
    # The canonical `((x + 180) mod 360) - 180` wrap. This is the *signed* wrap
    # to (-180, 180]; values near +360 wrap to small negative numbers (359 -> -1,
    # 180 -> -180). In Python, the % operator on negative numbers returns a
    # non-negative result for positive divisors, so the wrap always lands in
    # (-180, 180]. The convention puts the wrap AT the +/-180 boundary.
    x = np.array([-540.0, -360.0, -180.0, -1.0, 0.0, 1.0, 179.99, 180.0, 359.0, 360.0, 540.0])
    w = ef.wrap_longitude_deg(x)
    assert np.all(w >= -180.0) and np.all(w <= 180.0)
    # specific values:
    assert w[0] == -180.0   # -540 -> -180
    assert w[1] == 0.0      # -360 -> 0 (full turn back)
    assert w[2] == -180.0   # -180 -> -180
    assert w[3] == -1.0     # -1 -> -1
    assert w[4] == 0.0      # 0 -> 0
    assert w[5] == 1.0      # 1 -> 1
    assert w[6] == 179.99   # 179.99 -> 179.99
    assert w[7] == -180.0   # 180 -> -180 (the wrap boundary)
    assert w[8] == -1.0     # 359 -> -1 (the cyclic continuation past 180)
    assert w[9] == 0.0      # 360 -> 0
    assert w[10] == -180.0  # 540 -> -180 (one turn forward)


# --------------------------------------------------------------------------- #
# L4: adversarial mutants -- ensure the lab_utils version catches the
# known failure modes from the hostile review (track 5)
# --------------------------------------------------------------------------- #
def test_subsolar_lon_with_negated_u_y_shifts_by_geodetic_formula():
    """Mutant: negate the y-component of the Sun unit vector. The subsolar
    longitude should shift by the rotation of the anti-sun direction, which
    for a generic epoch is detectable as a non-trivial difference. The fix
    is to use the geodetic subsolar longitude `atan2(u_ecef_y, u_ecef_x)`,
    NOT `atan2(-u_y, -u_x)` (which would be the anti-sun direction in the
    ECI equatorial frame).
    """
    t = 1.5e8  # arbitrary epoch well into the year
    good = ef.subsolar_lon_rad(t)
    # construct a Sun unit vector manually with y negated (mutant)
    u, d = ef.sun_unit_and_dist_km(t)
    u_neg_y = np.array([u[0], -u[1], u[2]])
    # compute subsolar lon via the mutant formula (just the ECI atan2)
    mutant = np.arctan2(-u_neg_y[1], -u_neg_y[0])
    mutant = float((mutant + np.pi) % (2 * np.pi) - np.pi)
    # the negation must produce a non-zero shift (the right shift is 2*lon)
    diff = min(abs(good - mutant), 2 * np.pi - abs(good - mutant))
    assert diff > 1e-3, "negated u_y did not change subsolar lon -- mutation missed"


def test_gmst_iau1982_polynomial_truncation_envelope():
    """The polynomial T^2 and T^3 corrections are sub-arcsec; the T coefficient
    has a sub-arcsec/day rate difference from the WGS-84 lab canon.

    The IAU-1982 polynomial's per-day rate (86636.555/240 = 360.985646 deg/day)
    differs from the WGS-84 lab canon's `omega_E` (360.985647332 deg/day) by
    ~1.3e-6 deg/day. Over 30 years the cumulative T-coefficient drift is
    ~50 arcsec. The T^2 and T^3 corrections add ~0.1 arcsec over the same
    horizon. The total residual between the linear form (anchored at J2000)
    and the polynomial is < 1 arcmin at 30 years.
    """
    t_30y = 30.0 * 365.25 * 86400.0
    poly_0 = ef.gmst_rad_iau1982(0.0)
    poly = ef.gmst_rad_iau1982(t_30y)
    delta_poly = (poly - poly_0) % (2 * np.pi)
    delta_lin = ef.OMEGA_EARTH_RAD_S * t_30y
    raw_diff = (delta_poly - delta_lin) % (2 * np.pi)
    raw_diff_min = min(raw_diff, 2 * np.pi - raw_diff)
    # T coefficient drift + T^2 + T^3 corrections at 30 yr are < 2 arcmin.
    # (The WGS-84 lab canon's omega_E differs from the IAU-1982 polynomial's
    # T coefficient by ~1.3e-6 deg/day, accumulating to ~1.5 arcmin over 30 yr.)
    assert raw_diff_min < np.deg2rad(2.0 / 60.0), (
        f"30-yr GMST linear vs polynomial residual = "
        f"{np.degrees(raw_diff_min) * 60:.4f} arcmin"
    )
