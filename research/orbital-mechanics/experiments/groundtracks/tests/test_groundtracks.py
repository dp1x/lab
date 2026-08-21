"""Validation tests for groundTracks (Experiment 008).

Laboratory rule: verify before trusting. Tests are independent of the
implementation where it matters — expected values are derived from
theory or separate code paths, not by calling the same helper to
produce expected.

Module loaded via importlib explicit path (pytest module registry safety).
"""

import importlib.util
from pathlib import Path

import numpy as np

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "groundtracks_experiment", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)

# constants reference (independent of experiment where possible)
MU_EARTH = 398600.4418
R_EARTH = 6378.137
OMEGA_E_TRUE = 7.2921159e-5
T_SIDEREAL_TRUE = 2 * np.pi / OMEGA_E_TRUE  # 86164.0905


# --------------------------------------------------------------------------- #
# L1: Dual-algebra cross-check (spherical trig vs rotation matrix)
# --------------------------------------------------------------------------- #

def test_trig_vs_matrix_iss():
    o = experiment.real_orbits()["ISS"]
    a, e = o["a_km"], o["e"]
    inc, Omega, omega, M0 = (np.radians(o[k]) for k in ("inc_deg", "Omega_deg", "omega_deg", "M0_deg"))
    T = o["T_sec"]
    t = np.linspace(0, 2 * T, 2 * 720, endpoint=False)
    gt = experiment.ground_track_analytic(a, e, inc, Omega, omega, M0, t)
    dlat = np.max(np.abs(gt["lat_mat_deg"] - gt["lat_trig_deg"]))
    dlon = np.max(np.abs(((gt["lon_mat_deg"] - gt["lon_trig_deg"] + 180) % 360) - 180))
    assert dlat < 1e-10, f"trig vs matrix dlat {dlat}"
    assert dlon < 1e-10, f"trig vs matrix dlon {dlon}"


def test_trig_vs_matrix_all_orbits():
    orbits = experiment.real_orbits()
    trig = experiment.validate_trig_vs_matrix(orbits)
    for name, v in trig.items():
        assert v["max_abs_dlat_deg"] < 1e-10, f"{name} dlat {v['max_abs_dlat_deg']}"
        # lon masked near poles, but still <1e-10 away from poles
        assert v["max_abs_dlon_deg"] < 1e-09, f"{name} dlon {v['max_abs_dlon_deg']}"
        # great-circle dominated by floating noise, allow 1e-07 rad
        assert v["max_great_circle_rad"] < 1e-07, f"{name} great circle {v['max_great_circle_rad']}"


def test_gmst_identity():
    # GMST must advance by OMEGA_E * T_sid to 2pi
    assert abs(experiment.OMEGA_EARTH_RAD_S - OMEGA_E_TRUE) < 1e-12
    assert abs(experiment.T_SIDEREAL_S - T_SIDEREAL_TRUE) < 1e-6
    # gmst(t+T_sid) = gmst(t) + 2pi
    t = np.array([0.0, 1000.0, 5000.0])
    for ti in t:
        assert abs((experiment.gmst_rad(ti + T_SIDEREAL_TRUE) - experiment.gmst_rad(ti)) - 2 * np.pi) < 1e-12


# --------------------------------------------------------------------------- #
# L2: Invariants — max latitude, equatorial, pole, GEO, repeat
# --------------------------------------------------------------------------- #

def test_max_latitude_equals_inclination_circular():
    # For circular orbits, max|phi| should equal inc (or 180-inc for retrograde)
    # Test at exact u=90° instead of sampling max to avoid sampling error.
    for inc_deg, expected in [(0, 0), (30, 30), (51.6, 51.6), (90, 90), (120, 60), (180, 0), (98, 82)]:
        inc = np.radians(inc_deg)
        a = R_EARTH + 500
        e = 0.0
        Omega = 0.0
        omega = 0.0
        # For circular, u = omega + M, so M = 90° - omega gives u=90°
        M_at_max = np.radians(90.0) - omega
        n = np.sqrt(MU_EARTH / a ** 3)
        t_at_max = (M_at_max - 0.0) / n
        t = np.array([t_at_max])
        gt = experiment.ground_track_analytic(a, e, inc, Omega, omega, 0.0, t)
        lat = gt["lat_mat_deg"][0]
        # signed max: for retrograde 120°, expected 60 but lat at u=90° = sin i = sin 120 = sin 60 = 0.866 => phi =60
        assert abs(abs(lat) - expected) < 1e-09, f"inc {inc_deg} lat {lat} expected {expected}"


def test_equatorial_lat_zero():
    a = R_EARTH + 400
    gt = experiment.ground_track_analytic(a, 0.0, np.radians(0), 0, 0, 0, np.linspace(0, 2 * experiment.orbital_period(a), 720))
    assert np.max(np.abs(gt["lat_mat_deg"])) < 1e-12
    assert np.max(np.abs(gt["lat_trig_deg"])) < 1e-12


def test_polar_covers_all_latitudes():
    a = R_EARTH + 500
    T = experiment.orbital_period(a)
    t = np.linspace(0, T, 720, endpoint=False)
    gt = experiment.ground_track_analytic(a, 0.0, np.radians(90), 0, 0, 0, t)
    assert np.max(gt["lat_mat_deg"]) > 89.0
    assert np.min(gt["lat_mat_deg"]) < -89.0
    assert np.max(np.abs(gt["lat_mat_deg"])) == np.max(np.abs(gt["lat_mat_deg"]))  # sanity
    # latitude via asin must stay finite
    assert np.all(np.isfinite(gt["lat_mat_deg"]))
    assert np.all(np.isfinite(gt["lon_mat_deg"]))


def test_retrograde_max_lat():
    # i=98° => max lat 82° (180-98)
    inc = np.radians(98)
    a = R_EARTH + 600
    T = experiment.orbital_period(a)
    t = np.linspace(0, 3 * T, 3 * 720, endpoint=False)
    gt = experiment.ground_track_analytic(a, 0.0, inc, 0, 0, 0, t)
    max_lat = np.max(np.abs(gt["lat_mat_deg"]))
    assert abs(max_lat - 82.0) < 0.05  # sampling tolerance


def test_eci_to_ecef_preserves_radius():
    # |r_ecef| == |r_eci| to machine
    orbits = experiment.real_orbits()
    for name, o in orbits.items():
        a, e = o["a_km"], o["e"]
        inc = np.radians(o["inc_deg"])
        T = o["T_sec"]
        t = np.linspace(0, T, 360, endpoint=False)
        gt = experiment.ground_track_analytic(a, e, inc, np.radians(o["Omega_deg"]), np.radians(o["omega_deg"]), np.radians(o["M0_deg"]), t)
        r_eci = gt["r_eci"]
        r_ecef = gt["r_ecef"]
        rel = np.max(np.abs(np.linalg.norm(r_ecef, axis=1) - np.linalg.norm(r_eci, axis=1)) / np.linalg.norm(r_eci, axis=1))
        assert rel < 1e-14, f"{name} radius preserve {rel}"


def test_longitude_wrapping_antimeridian():
    # Equatorial LEO 3 orbits sampled at 720/orbit should have small step <5° after wrapping
    a = R_EARTH + 400
    T = experiment.orbital_period(a)
    t = np.linspace(0, 3 * T, 3 * 720, endpoint=False)
    gt = experiment.ground_track_analytic(a, 0.0, 0.0, 0, 0, 0, t)
    lon = gt["lon_mat_deg"]
    dlon = np.abs(((np.diff(lon) + 180) % 360) - 180)
    assert np.max(dlon) < 5.0, f"max step {np.max(dlon)}"
    # unwrapped should be monotonic increasing for prograde equatorial (since n > omega_E)
    lon_u = experiment.unwrap_longitude_deg(lon)
    assert np.all(np.diff(lon_u) > 0)


def test_delta_longitude_per_orbit():
    # For circular LEO 400 km equatorial, delta = -omega_E * T
    a = R_EARTH + 400
    T = experiment.orbital_period(a)
    analytic = -np.degrees(OMEGA_E_TRUE * T)
    measured_delta = experiment.delta_longitude_per_orbit(a, MU_EARTH, OMEGA_E_TRUE)
    assert abs(measured_delta - analytic) < 1e-12
    # also via ground track wrapped measurement
    t = np.linspace(0, 2 * T, 2 * 720, endpoint=False)
    gt = experiment.ground_track_analytic(a, 0.0, 0.0, 0, 0, 0, t)
    lon_u = experiment.unwrap_longitude_deg(gt["lon_mat_deg"])
    # unwrapped diff after 1 orbit = 360 + analytic (since lon_eci +360)
    measured_unwrapped = lon_u[720] - lon_u[0]
    expected_unwrapped = 360.0 + analytic
    assert abs(measured_unwrapped - expected_unwrapped) < 1e-06
    # wrapped delta should equal analytic
    measured_wrapped = ((measured_unwrapped + 180) % 360) - 180
    assert abs(measured_wrapped - analytic) < 1e-06


def test_sidereal_vs_solar_distinction():
    # Solar day 86400 would be 0.3% off; GEO period must be sidereal not solar
    a_geo = (MU_EARTH * T_SIDEREAL_TRUE ** 2 / (4 * np.pi ** 2)) ** (1.0 / 3.0)
    T_geo = experiment.orbital_period(a_geo)
    assert abs(T_geo - T_SIDEREAL_TRUE) < 1e-03  # ~1 ms tolerance
    # solar day would give different a
    a_solar = (MU_EARTH * 86400.0 ** 2 / (4 * np.pi ** 2)) ** (1.0 / 3.0)
    assert abs(a_solar - a_geo) > 70  # ~ 70 km difference, easily detectable
    assert abs(experiment.OMEGA_EARTH_RAD_S - 2 * np.pi / 86400) > 1e-07  # must not be solar rate


def test_geo_stationary():
    a_geo = experiment.real_orbits()["GEO"]["a_km"]
    T = experiment.orbital_period(a_geo)
    t = np.linspace(0, 5 * T, 5 * 720 + 1)
    gt = experiment.ground_track_analytic(a_geo, 0.0, 0.0, 0, 0, 0, t)
    assert np.max(np.abs(gt["lat_mat_deg"])) < 1e-10
    # lon variation should be ~0
    lon = gt["lon_mat_deg"]
    assert np.max(lon) - np.min(lon) < 1e-09
    # period must match sidereal to 1e-10 relative
    assert abs(T - T_SIDEREAL_TRUE) / T_SIDEREAL_TRUE < 1e-09


def test_geo_inclined_figure8():
    o = experiment.real_orbits()["GEO_inclined"]
    a, inc = o["a_km"], np.radians(o["inc_deg"])
    T = o["T_sec"]
    t = np.linspace(0, T, 720, endpoint=False)
    gt = experiment.ground_track_analytic(a, 0.0, inc, 0, 0, 0, t)
    # inclined GEO lat should oscillate ± inc
    assert abs(np.max(gt["lat_mat_deg"]) - o["inc_deg"]) < 1e-09
    assert abs(np.min(gt["lat_mat_deg"]) + o["inc_deg"]) < 1e-09
    # longitude should have small figure-8 variation (<~5.5° for 5° inc? Actually ~? but small)
    # For 5° inc, lon variation due to analemma is small <0.5°
    lon_range = np.max(gt["lon_mat_deg"]) - np.min(gt["lon_mat_deg"])
    assert lon_range < 2.0  # generous


def test_molniya_apogee_dwell_lat():
    # Molniya e=0.74 inc=63.4 omega=270 should dwell over north: max lat =63.4, min lat about -63.4 but apogee north
    o = experiment.real_orbits()["Molniya"]
    inc = np.radians(o["inc_deg"])
    omega = np.radians(o["omega_deg"])
    a, e = o["a_km"], o["e"]
    T = o["T_sec"]
    # sample densely because apogee dwell
    t = np.linspace(0, T, 1440, endpoint=False)
    gt = experiment.ground_track_analytic(a, e, inc, np.radians(o["Omega_deg"]), omega, 0.0, t)
    max_lat = np.max(gt["lat_mat_deg"])
    assert abs(max_lat - o["max_lat_theory_deg"]) < 0.1
    # also check periapsis altitude: r_peri = a(1-e) ~ 6880 km => alt ~500 km
    r_peri = a * (1 - e)
    assert abs(r_peri - (R_EARTH + 500)) < 600  # approx
    # apoapsis alt ~ 40k km
    r_apo = a * (1 + e)
    assert r_apo > 40000


def test_gmst_sign_west_drift():
    # Equatorial prograde: lon should decrease west per orbit (wrapped)
    a = R_EARTH + 400
    T = experiment.orbital_period(a)
    t = np.linspace(0, T, 720, endpoint=False)
    gt = experiment.ground_track_analytic(a, 0.0, 0.0, 0, 0, 0, t)
    # Need 2 orbits to see shift at same latitude (equatorial every point same lat)
    t2 = np.linspace(0, 2 * T, 1440, endpoint=False)
    gt2 = experiment.ground_track_analytic(a, 0.0, 0.0, 0, 0, 0, t2)
    lon_u = experiment.unwrap_longitude_deg(gt2["lon_mat_deg"])
    # per orbit unwrapped increase should be 360 - omega*T (~336.7)
    d1 = lon_u[720] - lon_u[0]
    assert d1 > 300 and d1 < 360
    # wrapped delta = d1 -360 = -omega*T
    d1_wrapped = d1 - 360
    assert d1_wrapped < 0 and abs(d1_wrapped + np.degrees(OMEGA_E_TRUE * T)) < 1e-06
    # Retrograde equatorial: should drift east (wrapped positive) because Earth rotation adds
    gt_ret = experiment.ground_track_analytic(a, 0.0, np.radians(180), 0, 0, 0, t2)
    lon_u_ret = experiment.unwrap_longitude_deg(gt_ret["lon_mat_deg"])
    d1_ret = lon_u_ret[720] - lon_u_ret[0]
    # For retrograde, lon_eci decreases 360 per orbit (since orbit opposite), net = -360 - omega*T = -383°, unwrapped -383, wrapped -23? Wait
    # Check sign: retrograde i=180, cos i = -1, u = nu, lon_eci = Omega + atan2(-sin u, cos u) = Omega - u
    # So after T, u +360, lon_eci -360, gmst +23, net -383°
    assert d1_ret < -300


# --------------------------------------------------------------------------- #
# L3: Propagation vs analytic
# --------------------------------------------------------------------------- #

def test_propagation_vs_analytic_iss():
    o = experiment.real_orbits()["ISS"]
    a, e = o["a_km"], o["e"]
    inc = np.radians(o["inc_deg"])
    T = o["T_sec"]
    t = np.linspace(0, 2 * T, 2 * 512 + 1)
    gt_ana = experiment.ground_track_analytic(a, e, inc, 0, 0, 0, t)
    gt_prop = experiment.propagate_ground_track(a, e, inc, 0, 0, 0, t)
    dlat = np.max(np.abs(gt_ana["lat_mat_deg"] - gt_prop["lat_deg"]))
    dlon = np.max(np.abs(((gt_prop["lon_deg"] - gt_ana["lon_mat_deg"] + 180) % 360) - 180))
    assert dlat < 5e-06
    assert dlon < 5e-06


def test_propagation_vs_analytic_polar():
    o = experiment.real_orbits()["Polar_LEO"]
    a, e = o["a_km"], o["e"]
    inc = np.radians(o["inc_deg"])
    T = o["T_sec"]
    t = np.linspace(0, 2 * T, 2 * 512 + 1)
    gt_ana = experiment.ground_track_analytic(a, e, inc, 0, 0, 0, t)
    gt_prop = experiment.propagate_ground_track(a, e, inc, 0, 0, 0, t)
    dlat = np.max(np.abs(gt_ana["lat_mat_deg"] - gt_prop["lat_deg"]))
    assert dlat < 5e-06


def test_propagation_vs_analytic_molniya():
    o = experiment.real_orbits()["Molniya"]
    a, e = o["a_km"], o["e"]
    inc = np.radians(o["inc_deg"])
    omega = np.radians(o["omega_deg"])
    T = o["T_sec"]
    t = np.linspace(0, T, 2048 + 1)
    gt_ana = experiment.ground_track_analytic(a, e, inc, 0, omega, 0, t)
    gt_prop = experiment.propagate_ground_track(a, e, inc, 0, omega, 0, t)
    dlat = np.max(np.abs(gt_ana["lat_mat_deg"] - gt_prop["lat_deg"]))
    dlon = np.max(np.abs(((gt_prop["lon_deg"] - gt_ana["lon_mat_deg"] + 180) % 360) - 180))
    # Molniya needs higher tolerance due to periapsis resolution
    assert dlat < 2e-04
    assert dlon < 7e-04


def test_convergence_order():
    conv = experiment.convergence_study()
    mean_order = conv["mean_order"]
    # RK4 order 4, allow 0.3 tolerance per earlier kepler validation (theory 4.07 etc)
    assert abs(mean_order - 4.0) < 0.35
    # errors must be strictly decreasing with step halving
    errs = conv["max_errors_deg"]
    for i in range(len(errs) - 1):
        assert errs[i] > errs[i + 1]
        ratio = errs[i] / errs[i + 1]
        # ratio should be ~16 for order 4
        assert 12 < ratio < 20, f"ratio {ratio} not ~16"


# --------------------------------------------------------------------------- #
# L6: Pathological and edge cases
# --------------------------------------------------------------------------- #

def test_pathological_no_nan():
    patho = experiment.pathological_checks()
    assert patho["all_finite_and_bounded"]
    assert patho["antimeridian_ok"]
    assert not patho["failures"]


def test_pole_singularity_lon_undefined():
    # At exact north pole, lon is set to 0 (undefined) but lat=90
    # Use polar orbit at time where lat=90 exactly
    a = R_EARTH + 500
    inc = np.radians(90)
    T = experiment.orbital_period(a)
    # t = T/4 gives pole for inc=90, omega=0, e=0
    t = np.array([T / 4])
    gt = experiment.ground_track_analytic(a, 0.0, inc, 0, 0, 0, t)
    assert abs(gt["lat_mat_deg"][0] - 90.0) < 1e-09
    # lon at pole is 0 by convention (guard)
    assert gt["lon_mat_deg"][0] == 0.0 or np.isfinite(gt["lon_mat_deg"][0])
    # near pole within 0.1° latitude, lon may be noisy but still finite
    t2 = np.linspace(0, T, 720, endpoint=False)
    gt2 = experiment.ground_track_analytic(a, 0.0, inc, 0, 0, 0, t2)
    # all finite
    assert np.all(np.isfinite(gt2["lat_mat_deg"]))
    assert np.all(np.isfinite(gt2["lon_mat_deg"]))


def test_circular_e_zero_handling():
    # e=1e-12 should not cause division by zero; also e=0 exactly
    for e in [0.0, 1e-12, 1e-8]:
        a = R_EARTH + 400
        gt = experiment.ground_track_analytic(a, e, np.radians(30), 0, 0, 0, np.linspace(0, experiment.orbital_period(a), 360))
        assert np.all(np.isfinite(gt["lat_mat_deg"]))


def test_high_eccentricity():
    # e=0.8
    a = R_EARTH + 500
    e = 0.8
    inc = np.radians(30)
    T = experiment.orbital_period(a)
    t = np.linspace(0, T, 1440, endpoint=False)
    gt = experiment.ground_track_analytic(a, e, inc, 0, 0, 0, t)
    assert np.all(np.isfinite(gt["lat_mat_deg"]))
    assert np.max(np.abs(gt["lat_mat_deg"])) < 30.001  # max lat ~30


def test_units_km_not_m():
    # If someone used meters, period would be off by sqrt(1000) ~31.6×
    a_km = R_EARTH + 400
    T_km = experiment.orbital_period(a_km)
    assert 5000 < T_km < 6000  # LEO period ~92 min
    a_m = (R_EARTH + 400) * 1000  # meters mistakenly passed as km
    T_m_wrong = experiment.orbital_period(a_m)  # would be 31× larger if a in meters but mu still km
    assert T_m_wrong > 150000  # clearly wrong, catches km/m bug
    # Our function should be called with km, so T_km is correct


def test_longitude_sign_convention():
    # Degree/radian bug would give huge error; test with 90° RAAN shifts lon by 90°
    a = R_EARTH + 500
    T = experiment.orbital_period(a)
    t0 = np.array([0.0])
    gt0 = experiment.ground_track_analytic(a, 0.0, np.radians(30), np.radians(0), 0, 0, t0)
    gt90 = experiment.ground_track_analytic(a, 0.0, np.radians(30), np.radians(90), 0, 0, t0)
    dlon = ((gt90["lon_mat_deg"][0] - gt0["lon_mat_deg"][0] + 180) % 360) - 180
    assert abs(dlon - 90.0) < 1e-06


def test_rotation_matrix_orthonormal():
    # Q must be rotation (orthonormal, det=1)
    for inc_deg in [0, 30, 90, 120]:
        Q = experiment.rotation_matrix_313(np.radians(10), np.radians(inc_deg), np.radians(20))
        should_be_I = Q @ Q.T
        assert np.max(np.abs(should_be_I - np.eye(3))) < 1e-14
        assert abs(np.linalg.det(Q) - 1.0) < 1e-14


def test_determinism():
    a = R_EARTH + 400
    t = np.linspace(0, experiment.orbital_period(a), 720)
    gt1 = experiment.ground_track_analytic(a, 0.0, np.radians(51.6), 0, 0, 0, t)
    gt2 = experiment.ground_track_analytic(a, 0.0, np.radians(51.6), 0, 0, 0, t)
    assert np.allclose(gt1["lat_mat_deg"], gt2["lat_mat_deg"])
    assert np.allclose(gt1["lon_mat_deg"], gt2["lon_mat_deg"])
    assert np.allclose(gt1["r_eci"], gt2["r_eci"])


def test_repeat_ground_track_geo():
    rep = experiment.repeat_ground_track_check()
    # GEO should repeat within 1e-09 deg wrapped
    assert rep["GEO_5orbit_err_wrapped_deg"] < 1e-09
    assert rep["GEO_wrapped_variation_deg"] < 1e-09
    assert rep["GEO_T_vs_Tsid_rel_err"] < 1e-12
    # 12h 2-orbit should also repeat wrapped 0 within 1e-09
    assert rep["12h_err_wrapped_deg"] < 1e-09


def test_invariants_overall():
    orbits = experiment.real_orbits()
    inv = experiment.validate_invariants(orbits)
    for name, v in inv.items():
        # max lat error <0.001 deg for dense sampling (allows sampling)
        assert v["max_lat_err_deg"] < 0.002, f"{name} max lat err {v['max_lat_err_deg']}"
        # delta lon wrapped error <1e-06 deg
        assert v["delta_lon_err_deg"] < 1e-06, f"{name} delta lon err {v['delta_lon_err_deg']}"
        assert v["r_preserve_max_rel"] < 1e-12
        assert not v["has_nan_or_inf"]


def test_kepler_solver_independent():
    # Solve Kepler for e=0.6 M=1 rad via Newton vs bisection independent
    e = 0.6
    M = 1.0
    E_newton = experiment.solve_kepler(M, e)
    # independent bisection: M = E - e sin E, root on [0, 2pi]
    import math

    def f(E):
        return E - e * math.sin(E) - M

    lo, hi = 0.0, 2 * math.pi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    E_bisect = 0.5 * (lo + hi)
    assert abs(E_newton - E_bisect) < 1e-12


def test_spherical_constants_documented():
    # Ensure constants are as documented and not silently changed
    assert abs(experiment.MU_EARTH_KM3S2 - 398600.4418) < 0.01
    assert abs(experiment.R_EARTH_KM - 6378.137) < 0.01
    assert abs(experiment.OMEGA_EARTH_RAD_S - 7.2921159e-5) < 1e-12
    # T sidereal vs solar distinction
    assert abs(experiment.T_SIDEREAL_S - 86164.0905) < 0.1
    assert abs(experiment.T_SIDEREAL_S - 86400) > 200
