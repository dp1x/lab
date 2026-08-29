"""Earth-frame helpers: Sun model, GMST polynomial, ECI<->ECEF, lat/lon, LST.

Graduated at Exp 015 (dawn-dusk SSO) for the second real consumer:
- Exp 014 (`eclipseTiming/experiment.py`) is the original donor
- Exp 015 (`dawnDuskSSO/experiment.py`) is the second consumer

Functions are verbatim transcriptions of the donor implementations where
possible, with a thin common interface (the lab's `t_s` is TT-like seconds
since J2000 throughout, see Exp 014 frozen contract).

Frame conventions (pinned by Exp 014):
- Earth: spherical, R_E = R_EARTH_KM (WGS-84 equatorial), uniform rotation at
  omega_E = OMEGA_EARTH_RAD_S (sidereal rate).
- Lab ECI: pseudo-inertial, J2000-anchored, but Sun direction is GEOMETRIC in
  mean equator/equinox of date (Astronomical Almanac low-precision formulas).
  The of-date vs J2000 precession bias is ~50.3"/yr -> at LEO event-rate
  floor of ~1 s/yr; named exclusion per Exp 014 contract.
- UTC: TT - 69.184 s (pinned Bulletin C era), UT1 := UTC, DUT1 = 0 frozen
  (+/- 0.9 s envelope disclosed).
- GMST: Aoki et al. 1982 IAU-1982 polynomial on UT1 (matches the Exp 014
  convention; the groundtracks linear `gmst0 + omega_E * t` is the simpler
  uniform-rotation approximation retained for invariant tests).

Covered by tests in ``src/lab_utils/tests/test_earth_frames.py`` with
regression pins against the donor at ``eclipseTiming/experiment.py`` and
``groundtracks/experiment.py`` (importlib-loaded for bit-equality).
"""
from __future__ import annotations

import math

import numpy as np

from lab_utils.orbits import (
    J2_EARTH,
    MU_EARTH_KM3S2,
    OMEGA_EARTH_RAD_S,
    R_EARTH_KM,
)

__all__ = [
    "AU_KM",
    "R_SUN_KM",
    "JD_J2000",
    "TT_MINUS_UTC_S",
    "DUT1_FROZEN_S",
    "DEG",
    "T_SIDEREAL_S",
    "sun_unit_and_dist_km",
    "subsolar_lon_rad",
    "subsolar_dec_rad",
    "gmst_rad_iau1982",
    "eci_to_ecef",
    "ecef_to_latlon",
    "spherical_trig_latlon",
    "wrap_longitude_deg",
    "lst_at_node_hours",
    "node_lon_from_raan_gmst",
]

# --------------------------------------------------------------------------- #
# Earth/Sun constants (provenance pinned by Exp 014)
# --------------------------------------------------------------------------- #
AU_KM = 149597870.7  # IAU 2012 Resolution B2 (exact)
R_SUN_KM = 695700.0  # IAU 2015 Resolution B3 nominal solar radius
JD_J2000 = 2451545.0  # Julian Date at J2000.0 TDB
TT_MINUS_UTC_S = 69.184  # IERS Bulletin C era (Exp 013/014 doctrine)
DUT1_FROZEN_S = 0.0  # declared; envelope +/- 0.9 s disclosed
DEG = math.pi / 180.0
T_SIDEREAL_S = 2.0 * math.pi / OMEGA_EARTH_RAD_S  # 86164.0905 s (Exp 008/014 canon)


# --------------------------------------------------------------------------- #
# Sun model: Astronomical Almanac low-precision (geometric, mean of date)
# --------------------------------------------------------------------------- #
def sun_unit_and_dist_km(t_s):
    """Geocentric geometric Sun unit vector (mean-of-date) and distance (km).

    Vectorized; t_s in lab TT-like seconds since J2000. Mean-of-date
    convention is consistent with Exp 014's frozen contract; the ICRF
    vs of-date precession bias is named-excluded per `localdocs/roadmap.md`
    and validated at ~0.65 deg against the byte-pinned 2026 Horizons Sun
    snapshot (Exp 014 G6 sun_validation gate, band 0.7 deg).
    """
    n = np.asarray(t_s, dtype=float) / 86400.0
    L = np.mod(280.460 + 0.9856474 * n, 360.0)
    g = np.mod(357.528 + 0.9856003 * n, 360.0)
    lam = np.deg2rad(L + 1.915 * np.sin(np.deg2rad(g)) + 0.020 * np.sin(np.deg2rad(2.0 * g)))
    eps = np.deg2rad(23.439 - 0.0000004 * n)
    u = np.stack(
        [np.cos(lam), np.cos(eps) * np.sin(lam), np.sin(eps) * np.sin(lam)], axis=-1
    )
    R_AU = 1.00014 - 0.01671 * np.cos(np.deg2rad(g)) - 0.00014 * np.cos(np.deg2rad(2.0 * g))
    return u, R_AU * AU_KM


def subsolar_lon_rad(t_s):
    """Subsolar-point geodetic longitude in the lab ECI mean-of-date frame (rad).

    Returns the *geodetic* (ECEF) longitude of the subsolar point: the
    longitude on Earth where the Sun is directly overhead. The Sun's
    apparent direction in ECI mean-of-date is `u = (u_x, u_y, u_z)`;
    after the ECI->ECEF rotation by the GMST, the geodetic coordinates
    of the subsolar point are at `atan2(u_ecef_y, u_ecef_x)`. Result in
    (-pi, pi] (west-positive, matching the lab's lat/lon wrap convention
    from Exp 008).

    This is *distinct* from the Sun's right ascension in ECI
    (which is `atan2(u_y, u_x)`). The two differ by the GMST:
    `subsolar_lon = alpha_sun - GMST` (mod 2*pi). The LST formula uses
    the geodetic subsolar longitude directly:
    `LST_at_node_lon = 12 + (node_lon - subsolar_lon) / 15 deg/h`.
    """
    u, _ = sun_unit_and_dist_km(t_s)
    u_ecef = eci_to_ecef(u, gmst_rad_iau1982(t_s))
    u_x_ecef = float(u_ecef[..., 0])
    u_y_ecef = float(u_ecef[..., 1])
    lon = np.arctan2(u_y_ecef, u_x_ecef)
    # wrap to (-pi, pi]
    lon = (lon + np.pi) % (2 * np.pi) - np.pi
    return float(lon) if np.ndim(lon) == 0 else lon


def subsolar_dec_rad(t_s):
    """Subsolar-point geocentric declination in the lab mean-of-date frame (rad).

    Companion to `subsolar_lon_rad`. Declination is `arcsin(u_z)` of the Sun
    unit vector. Vectorized.
    """
    u, _ = sun_unit_and_dist_km(t_s)
    u_z = np.asarray(u[..., 2], dtype=float)
    u_z = np.clip(u_z, -1.0, 1.0)
    dec = np.arcsin(u_z)
    return float(dec) if np.ndim(dec) == 0 else dec


# --------------------------------------------------------------------------- #
# GMST: Aoki et al. 1982 (IAU-1982) polynomial, matches eclipseTiming.frozen
# --------------------------------------------------------------------------- #
def gmst_rad_iau1982(t_s):
    """GMST (rad) on the lab TT-like time scale (Aoki et al. 1982 IAU-1982).

    Per the Exp 014 frozen contract: UT1 := UTC = TT - TT_MINUS_UTC_S - DUT1_FROZEN_S.
    Returns radians; the standard `gmst_sec/240 * DEG` form. Not wrapped to
    [0, 2 pi); callers wrap as needed (subtractive differences are unaffected).
    """
    jd_ut1 = JD_J2000 + (np.asarray(t_s, dtype=float) - TT_MINUS_UTC_S - DUT1_FROZEN_S) / 86400.0
    Tu = (jd_ut1 - JD_J2000) / 36525.0
    gmst_sec = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * Tu
        + 0.093104 * Tu**2
        - 6.2e-6 * Tu**3
    )
    if np.ndim(gmst_sec) == 0:
        gmst_sec = float(gmst_sec)
        gmst_sec = math.fmod(gmst_sec, 86400.0)
        if gmst_sec < 0.0:
            gmst_sec += 86400.0
        return gmst_sec / 240.0 * DEG
    gs = np.mod(np.asarray(gmst_sec, dtype=float), 86400.0)
    return gs / 240.0 * DEG


# --------------------------------------------------------------------------- #
# ECI <-> ECEF: passive Z rotation, west-positive lon drift (Exp 008 canon)
# --------------------------------------------------------------------------- #
def eci_to_ecef(r_eci, gmst):
    """Rotate ECI vectors to ECEF via passive Z rotation R(theta_G).

    R = [[c s 0], [-s c 0], [0 0 1]]  with c=cos(theta), s=sin(theta).
    Then lon_ecef = atan2(y_ecef, x_ecef) = lon_eci - theta_G (west drift).
    Preserves |r| to machine precision. Vectorized over leading time axis.
    r_eci: (N, 3) or (3,). gmst: scalar or (N,). Returns same shape as input.
    """
    r = np.asarray(r_eci, dtype=float)
    th = np.asarray(gmst, dtype=float)
    single = r.ndim == 1 and th.ndim == 0
    if r.ndim == 1:
        r = r[None, :]
    if th.ndim == 0:
        th = np.full(r.shape[0], float(th), dtype=float)
    else:
        th = th.reshape(-1)
        assert th.shape[0] == r.shape[0], "gmst and r_eci length mismatch"
    c = np.cos(th)
    s = np.sin(th)
    out = np.empty_like(r)
    out[:, 0] = c * r[:, 0] + s * r[:, 1]
    out[:, 1] = -s * r[:, 0] + c * r[:, 1]
    out[:, 2] = r[:, 2]
    if single:
        return out[0]
    return out


def ecef_to_latlon(r_ecef):
    """Geocentric lat/lon from ECEF vectors (spherical Earth).

    phi = arcsin(z / r)  in [-90, +90] deg
    lon = atan2(y, x)    in (-180, +180] deg wrapped (west-positive)

    Returns (lat_deg, lon_deg, r_norm_km) arrays shape (N,) or scalars.
    At the exact pole (hypot(x, y) ~ 0) lon is set to 0 (undefined; the wrap
    is unstable there). This matches the Exp 008 / groundtracks convention.
    """
    r = np.asarray(r_ecef, dtype=float)
    single = r.ndim == 1
    if single:
        r = r[None, :]
    x, y, z = r[:, 0], r[:, 1], r[:, 2]
    rn = np.sqrt(x * x + y * y + z * z)
    sin_phi = np.clip(z / rn, -1.0, 1.0)
    lat_rad = np.arcsin(sin_phi)
    lon_rad = np.arctan2(y, x)
    eps = 1e-12
    pole_mask = np.hypot(x, y) < eps * rn
    lon_rad = np.where(pole_mask, 0.0, lon_rad)
    lat_deg = np.degrees(lat_rad)
    lon_deg = np.degrees(lon_rad)
    lon_deg = ((lon_deg + 180.0) % 360.0) - 180.0
    if single:
        return float(lat_deg[0]), float(lon_deg[0]), float(rn[0])
    return lat_deg, lon_deg, rn


def spherical_trig_latlon(inc, Omega, omega, nu, gmst):
    """Independent spherical-trig lat/lon (different algebra, same physics).

    u = omega + nu (argument of latitude)
    sin phi = sin i * sin u  -> phi = arcsin(clip)
    lon_eci = Omega + atan2(cos i * sin u, cos u)
    lon_ecef = lon_eci - gmst (wrapped to (-pi, pi] then deg)

    inc, Omega, omega in rad; nu, gmst arrays (rad). Returns (lat_deg, lon_deg).
    Vectorized. Independent of the matrix `eci_to_ecef` path; cross-check
    against `ecef_to_latlon` is a known-good invariant (Exp 008 test suite).
    """
    u = np.asarray(omega, dtype=float) + np.asarray(nu, dtype=float)
    gmst_arr = np.asarray(gmst, dtype=float)
    sg = np.sin(inc) * np.sin(u)
    sg = np.clip(sg, -1.0, 1.0)
    lat_rad = np.arcsin(sg)
    lon_eci = Omega + np.arctan2(np.cos(inc) * np.sin(u), np.cos(u))
    lon_rad = lon_eci - gmst_arr
    lon_rad = (lon_rad + np.pi) % (2 * np.pi) - np.pi
    lat_deg = np.degrees(lat_rad)
    lon_deg = np.degrees(lon_rad)
    lon_deg = ((lon_deg + 180.0) % 360.0) - 180.0
    return lat_deg, lon_deg


def wrap_longitude_deg(lon_deg):
    """Wrap degrees to (-180, 180]."""
    return (np.asarray(lon_deg, dtype=float) + 180.0) % 360.0 - 180.0


# --------------------------------------------------------------------------- #
# LST: local solar time at a given geodetic/ECEF longitude
# --------------------------------------------------------------------------- #
def node_lon_from_raan_gmst(raan_rad, gmst_rad_val):
    """Geodetic longitude (rad) of the orbit-plane node at a given GMST.

    The inertial-frame RAAN (in lab ECI) lies on the longitude
    `raan - theta_G` in ECEF (the satellite crosses the equator at this
    geodetic longitude at the moment its GMST matches). Returns rad,
    wrapped to (-pi, pi].
    """
    lon = np.asarray(raan_rad, dtype=float) - np.asarray(gmst_rad_val, dtype=float)
    lon = (lon + np.pi) % (2 * np.pi) - np.pi
    return float(lon) if np.ndim(lon) == 0 else lon


def lst_at_node_hours(t_s, node_lon_rad, *, subsolar_lon_rad_fn=None):
    """Local solar time at the geodetic longitude `node_lon_rad` (hours, 0-24).

    The LST at a geodetic longitude `node_lon_rad` is
    `LST = 12 h + (node_lon_rad - subsolar_lon) / (15 deg/h)` (mod 24),
    where `subsolar_lon` is the geodetic subsolar longitude from
    `subsolar_lon_rad(t_s)`. This is the apparent LST (no EoT correction
    needed because the subsolar point is by definition the apparent Sun
    direction projected to the geoid).

    Two independent paths for the LST at a given point on Earth:
    - Path 1: `LST = 12 + (node_lon - subsolar_lon) / 15` (this function).
    - Path 2: `LST = 12 + (GMST + node_lon - alpha_sun) / 15` where
      `alpha_sun = atan2(u_y, u_x)` is the Sun's right ascension in ECI.
    The two are bit-equivalent because `subsolar_lon = alpha_sun - GMST`.

    `node_lon_rad` is the geodetic (ECEF) longitude of the point of
    interest (e.g., the orbit-plane ascending node, which is at
    `Omega - GMST` in ECI). Returns a value in [0, 24).
    """
    if subsolar_lon_rad_fn is None:
        sub_lon = subsolar_lon_rad(t_s)
    else:
        sub_lon = np.asarray(subsolar_lon_rad_fn(t_s), dtype=float)
    delta_h = (np.asarray(node_lon_rad, dtype=float) - sub_lon) / (15.0 * DEG)
    lst = 12.0 + delta_h
    lst = lst - 24.0 * np.floor(lst / 24.0)
    return float(lst) if np.ndim(lst) == 0 else lst
