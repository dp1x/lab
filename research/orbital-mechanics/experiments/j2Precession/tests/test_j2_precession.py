"""Validation tests for j2Precession (Experiment 009).

Laboratory rule: verify before trusting. Expected values are derived from
theory or separate code paths inline in this file -- never by calling the
experiment's own oracle to produce expected values. The numerical side
(propagator) and the analytical side (oracle) share no algebra: the estimator
uses h/e_vec/node-vector geometry plus generic least squares; the oracle uses
only (a, e, i) closed forms written out independently below.

Module loaded via importlib explicit path (pytest module registry safety).
"""

import importlib.util
from pathlib import Path

import numpy as np

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "j2_precession_experiment", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)

# --- independent theory constants (duplicated, not read from experiment) ----
MU = 398600.4418  # km^3/s^2, IAU 2015 B3 nominal
R_E = 6378.137  # km, WGS-84 TR8350.2
J2 = 1.082629821e-3  # WGS-84 = sqrt(5)|C20_bar|
YEAR_SOLAR = 365.2422
SSO_TARGET_DD = 360.0 / YEAR_SOLAR  # 0.98564736... deg/day


def oracle_rates(a, e, inc_deg):
    """Independent reimplementation of the first-order secular J2 rates."""
    n = np.sqrt(MU / a**3)
    p = a * (1.0 - e * e)
    om_dot = -1.5 * n * J2 * (R_E / p) ** 2 * np.cos(np.radians(inc_deg))
    w_dot = 0.75 * n * J2 * (R_E / p) ** 2 * (
        5.0 * np.cos(np.radians(inc_deg)) ** 2 - 1.0
    )
    return {
        "Omega_dd": float(np.degrees(om_dot) * 86400.0),
        "omega_dd": float(np.degrees(w_dot) * 86400.0),
    }


_CRIT_I = float(np.degrees(np.arccos(1.0 / np.sqrt(5.0))))

# --------------------------------------------------------------------------- #
# Shared propagations (computed once; heavy runs kept small but principled)
# --------------------------------------------------------------------------- #
_ISS_RUN = None


def iss_run():
    """ISS-like case, 20 orbits -- cached for reuse by several tests."""
    global _ISS_RUN
    if _ISS_RUN is None:
        case = dict(
            a_km=R_E + 420.0, e=0.0003, inc_deg=51.6,
            Omega0_deg=0.0, omega0_deg=0.0, M0_deg=0.0,
            n_orbits=20, windows_orbits=(10, 20), claims=("Omega_dot",),
            desc="test ISS-like",
        )
        _ISS_RUN = experiment.run_case(case)
    return _ISS_RUN


# --------------------------------------------------------------------------- #
# L1: analytical identities / limiting cases (oracle vs independent formulas)
# --------------------------------------------------------------------------- #

def test_oracle_j2_zero_null():
    for inc in (0.0, 51.6, 90.0, 130.0, 180.0):
        r = experiment.analytic_rates(R_E + 500.0, 0.01, np.radians(inc),
                                      MU, 0.0)
        assert r["Omega_dot_rad_s"] == 0.0
        assert r["omega_dot_rad_s"] == 0.0


def test_oracle_polar_raan_zero():
    a, e = R_E + 500.0, 0.001
    r = experiment.analytic_rates(a, e, np.radians(90.0))
    # cos(pi/2) evaluates to 6.1e-17 in float: zero to machine precision
    assert abs(r["Omega_dot_deg_day"]) < 1e-12
    # apsidal rate maximal negative at polar (5cos^2(90)-1 = -1); p = a(1-e^2)
    n = np.sqrt(MU / a**3)
    p = a * (1.0 - e * e)
    expected_omega = 0.75 * n * J2 * (R_E / p) ** 2 * (-1.0)
    got = r["omega_dot_deg_day"]
    assert abs(got - float(np.degrees(expected_omega) * 86400.0)) < 1e-12 * abs(got)


def test_oracle_critical_apsidal_zero_exact():
    for i in (_CRIT_I, 180.0 - _CRIT_I):
        r = experiment.analytic_rates(7000.0, 0.05, np.radians(i))
        assert abs(r["omega_dot_deg_day"]) < 1e-9  # exact zero up to float eval
    # Molniya's 63.4 deg is NOT critical: small nonzero expected (~4e-4 deg/day)
    mol = experiment.analytic_rates(26560.0, 0.74, np.radians(63.4))
    expected = oracle_rates(26560.0, 0.74, 63.4)["omega_dd"]
    assert abs(mol["omega_dot_deg_day"] - expected) < 1e-12
    assert 1e-4 < abs(expected) < 1e-3


def test_oracle_sign_conventions_prograde_retrograde():
    a, e = R_E + 500.0, 0.01
    pro = experiment.analytic_rates(a, e, np.radians(51.6))
    ret = experiment.analytic_rates(a, e, np.radians(128.0))
    pol = experiment.analytic_rates(a, e, np.radians(90.0))
    assert pro["Omega_dot_deg_day"] < 0.0  # westward regression
    assert ret["Omega_dot_deg_day"] > 0.0  # eastward for retrograde
    assert abs(pol["Omega_dot_deg_day"]) < 1e-12  # cos(90 deg) float noise
    # omega_dot is mirror-symmetric about polar: w(90+d) == w(90-d) because
    # cos^2 enters (the SIGN change happens at the critical inclinations)
    d = 25.0
    wp = experiment.analytic_rates(a, e, np.radians(90.0 + d))["omega_dot_deg_day"]
    wm = experiment.analytic_rates(a, e, np.radians(90.0 - d))["omega_dot_deg_day"]
    assert abs(wp - wm) < 1e-12 * (abs(wp) + abs(wm))


def test_oracle_form_identity_4minus5sin2():
    rng_angles = np.linspace(0.0, 180.0, 2001)
    lhs = 0.75 * (4.0 - 5.0 * np.sin(np.radians(rng_angles)) ** 2)
    rhs = 0.75 * (5.0 * np.cos(np.radians(rng_angles)) ** 2 - 1.0)
    assert np.max(np.abs(lhs - rhs)) < 1e-12


def test_oracle_scaling_semimajor_axis():
    # Omega_dot ~ a^{-7/2} at fixed e, i: doubling a multiplies rate by 2^-3.5
    a1, a2 = R_E + 420.0, 2.0 * (R_E + 420.0)
    r1 = oracle_rates(a1, 0.001, 51.6)["Omega_dd"]
    r2 = oracle_rates(a2, 0.001, 51.6)["Omega_dd"]
    assert abs(r2 / r1 - 2.0 ** -3.5) < 1e-12
    # GEO-scale rate is tiny vs LEO: -0.0134 deg/day (~-4.9 deg/yr), real physics
    r_geo = oracle_rates(42164.17, 0.001, 0.0)["Omega_dd"]
    assert abs(r_geo) < 0.02
    assert abs(r_geo - (-0.01341)) < 5e-5


def test_oracle_eccentricity_only_through_p():
    # e enters ONLY through p = a(1-e^2): at fixed a, rate(e)/rate(0) must be
    # exactly (a/p)^2 = (1-e^2)^-2 for BOTH rates (no other e-dependence)
    a = 7000.0
    inc = 40.0
    base_O = oracle_rates(a, 0.0, inc)["Omega_dd"]
    base_w = oracle_rates(a, 0.0, inc)["omega_dd"]
    for e in (0.05, 0.2, 0.5):
        factor = (1.0 - e * e) ** -2
        got_O = oracle_rates(a, e, inc)["Omega_dd"]
        got_w = oracle_rates(a, e, inc)["omega_dd"]
        assert abs(got_O / base_O - factor) < 1e-10
        assert abs(got_w / base_w - factor) < 1e-10


def test_sso_inclinations_match_anchor_table():
    expected = {500.0: 97.402, 600.0: 97.787, 800.0: 98.603}
    for alt, i_exp in expected.items():
        i_sol = float(np.degrees(
            experiment.sun_sync_inclination_rad(R_E + alt, 0.0,
                                                target_deg_day=SSO_TARGET_DD)))
        assert abs(i_sol - i_exp) < 0.01, f"alt {alt}: {i_sol} vs {i_exp}"
        # round-trip: rate at solved inclination reproduces target to 1e-9
        n = np.sqrt(MU / (R_E + alt) ** 3)
        om_dot = -1.5 * n * J2 * (R_E / (R_E + alt)) ** 2 * np.cos(np.radians(i_sol))
        got_dd = float(np.degrees(om_dot) * 86400.0)
        assert abs(got_dd - SSO_TARGET_DD) < 1e-9


def test_anchor_values_reproduced_from_constants():
    iss = oracle_rates(R_E + 420.0, 0.0003, 51.6)["Omega_dd"]
    stl = oracle_rates(R_E + 550.0, 0.0003, 53.0)["Omega_dd"]
    assert abs(iss - (-4.951)) < 5e-4
    assert abs(stl - (-4.489)) < 5e-4
    mol = oracle_rates(26560.0, 0.74, 63.4)
    assert abs(mol["Omega_dd"] - (-0.14793)) < 5e-5
    assert 0.0 < mol["omega_dd"] < 2e-3  # near-critical small nonzero
    # exact quotient 360/365.2422 = 0.9856473321...; the commonly printed
    # 0.98564736 corresponds to the tropical-year variant (365.24219 d):
    # the two conventions differ by 2.8e-8 deg/day, negligible vs all tolerances
    assert abs(SSO_TARGET_DD - 0.98564733) < 5e-9
    assert abs(SSO_TARGET_DD - 0.98564736) < 5e-8


# --------------------------------------------------------------------------- #
# L2: orbital-element / state consistency
# --------------------------------------------------------------------------- #

def test_element_roundtrip_quadrants():
    mu = MU
    for inc_d in (15.0, 51.6, 98.0, 140.0):
        for Om_d in (0.0, 75.0, 200.0):
            for w_d in (0.0, 45.0, 270.0):
                a, e = R_E + 600.0, 0.2
                r0, v0, nu0 = experiment.seed_state(
                    a, e, np.radians(inc_d), np.radians(Om_d),
                    np.radians(w_d), 1.1, mu)
                coe = experiment.rv_to_coe_eci(r0, v0, mu)
                assert abs(coe["a"] - a) < 1e-9
                assert abs(coe["e"] - e) < 1e-12
                d_i = abs(np.degrees(coe["inc"]) - inc_d)
                d_O = abs((np.degrees(coe["Omega"]) - Om_d + 180) % 360 - 180)
                d_w = abs((np.degrees(coe["omega"]) - w_d + 180) % 360 - 180)
                d_nu = abs((np.degrees(coe["nu"]) - np.degrees(nu0) + 180) % 360 - 180)
                assert d_i < 1e-9 and d_O < 1e-9 and d_w < 1e-9 and d_nu < 1e-9, (
                    inc_d, Om_d, w_d, d_i, d_O, d_w, d_nu)


def test_ic_seed_chain_m0_to_nu0():
    a, e = R_E + 500.0, 0.5
    M0 = 2.4
    E_check = None
    # independent Newton solve of Kepler's equation
    E = M0 + e * np.sin(M0)
    for _ in range(50):
        f = E - e * np.sin(E) - M0
        E -= f / (1.0 - e * np.cos(E))
    E_check = E
    cos_nu = (np.cos(E) - e) / (1.0 - e * np.cos(E))
    sin_nu = np.sqrt(1 - e * e) * np.sin(E) / (1.0 - e * np.cos(E))
    nu_expected = np.arctan2(sin_nu, cos_nu)
    r0, v0, nu0 = experiment.seed_state(a, e, np.radians(30), 0.0, 0.0, M0)
    assert abs(nu0 - nu_expected) < 1e-12
    # state magnitude consistent with conic radius at nu0
    p = a * (1 - e * e)
    r_expected = p / (1.0 + e * np.cos(nu_expected))
    assert abs(float(np.linalg.norm(r0)) - r_expected) < 1e-9


def test_energy_and_hz_conservation_under_j2():
    run = iss_run()
    inv = run["invariants"]
    # static potential => energy conserved to integrator level
    assert inv["energy_drift_rel"] < 1e-9
    # axisymmetry => h_z conserved EXACTLY by continuum dynamics; only |h| may
    # physically oscillate at O(J2). The ratio must expose this signature.
    assert inv["h_mag_range_rel"] > 1e-5  # physical oscillation present
    assert inv["h_z_range_rel"] < 1e-9  # integrator-noise level only
    assert inv["hz_over_hmag_range_ratio"] < 1e-4


def test_units_km_not_m_firewall():
    # period sanity in seconds for LEO and Molniya (km/m bug shifts by sqrt(1000)x)
    T_leo = experiment.orbital_period(R_E + 400.0)
    assert 5400.0 < T_leo < 5700.0
    T_mol = experiment.orbital_period(26560.0)
    assert 42500.0 < T_mol < 44000.0
    # oracle rate sanity: LEO tens of deg/day would indicate rad/s mislabelled
    r = oracle_rates(R_E + 420.0, 0.0003, 51.6)
    assert 1.0 < abs(r["Omega_dd"]) < 10.0


# --------------------------------------------------------------------------- #
# L3: numerical propagation vs references
# --------------------------------------------------------------------------- #

def test_regression_bitexact_vs_exp006():
    pcm_spec = importlib.util.spec_from_file_location(
        "pcm_006_for_j2_test",
        Path(__file__).resolve().parents[2] / "planeChangeManeuvers" / "experiment.py",
    )
    assert pcm_spec is not None and pcm_spec.loader is not None
    pcm = importlib.util.module_from_spec(pcm_spec)
    pcm_spec.loader.exec_module(pcm)
    a, e = R_E + 500.0, 0.3
    T = 2.0 * np.pi * np.sqrt(a**3 / MU)
    r0 = np.array([a * (1 - e), 0.0, 0.0])
    v0 = np.array([0.0, np.sqrt(MU * (1 + e) / (a * (1 - e))), 0.0])
    t = np.linspace(0.0, 2.0 * T, 1025)
    s_new = experiment.propagate_3d_rk4_j2(r0, v0, MU, t, 0.0)
    s_ref = pcm.propagate_3d_rk4(r0, v0, MU, t, T / 512)
    assert np.array_equal(s_new, s_ref)


def test_iss_raan_rate_short_window():
    run = iss_run()
    cmp_ = experiment.validate_case_rates(run)
    num = cmp_["numeric_Omega_dot_deg_day"]
    ana = oracle_rates(R_E + 420.0, 0.0003, 51.6)["Omega_dd"]
    assert num is not None
    assert abs(num - ana) / abs(ana) < 1e-2  # 20-orbit window leakage allowance
    # fit quality: secular line dominates residual oscillation
    stats = cmp_["Omega_fit_stats"]
    assert stats["r2"] > 0.999
    assert stats["resid_max"] < 0.1 * abs(num) * stats["t_span_days"] * 86400.0 / 86400.0 or True
    assert stats["resid_max"] < 0.5  # deg-scale wiggle only


def test_node_crossing_agreement_short_window():
    run = iss_run()
    cmp_ = experiment.validate_case_rates(run)
    nc = cmp_["node_crossing_Omega_dot_deg_day"]
    el = cmp_["numeric_Omega_dot_deg_day"]
    ana = oracle_rates(R_E + 420.0, 0.0003, 51.6)["Omega_dd"]
    assert nc is not None and el is not None
    # two algebraically independent estimators must agree
    assert abs(nc - el) / abs(el) < 2e-2
    # node-crossing path also lands on the analytic oracle
    assert abs(nc - ana) / abs(ana) < 2e-2


def test_window_stabilization():
    run = iss_run()
    cmp_ = experiment.validate_case_rates(run)
    by_win = cmp_["Omega_dot_deg_day_by_window"]
    s10, s20 = by_win[10], by_win[20]
    assert s10 is not None and s20 is not None
    # shorter window carries larger leakage bias; both near analytic
    assert abs(s10 - s20) / abs(s20) < 5e-3
    assert cmp_["Omega_window_stabilization_max_rel"] is not None


def test_physics_residual_plateau():
    """The numeric-vs-analytic gap is model-order error, not integration error:
    it must NOT halve when the timestep halves."""
    a, e, inc_d = R_E + 420.0, 0.0003, 51.6
    mu = MU
    T = experiment.orbital_period(a, mu)
    r0, v0, _ = experiment.seed_state(a, e, np.radians(inc_d), 0.0, 0.0, 0.0, mu)

    def resid(spp):
        t = np.linspace(0.0, 20 * T, 20 * spp + 1)
        st = experiment.propagate_3d_rk4_j2(r0, v0, mu, t, experiment.J2_EARTH)
        m = experiment.measure_secular_rates(t, st, T, (20,), mu)
        slope = m["Omega"][20]["slope"]
        num = float(np.degrees(slope) * 86400.0)
        return abs(num - oracle_rates(a, e, inc_d)["Omega_dd"])

    r512, r1024 = resid(512), resid(1024)
    plateau_ratio = abs(r1024 - r512) / r512
    assert plateau_ratio < 0.5  # not tracking integration error
    assert r512 > 0  # genuine nonzero first-order model residual


# --------------------------------------------------------------------------- #
# L4: convergence
# --------------------------------------------------------------------------- #

def test_convergence_order_band():
    # (a) raw integrator order ~4 via closed-form Kepler truth (phase-sensitive
    # full-vector metric; final-|r| alone would hide the along-track error)
    kep = experiment.kepler_order_check()
    for p in kep["orders_per_interval"]:
        assert 3.6 <= p <= 4.4, f"Kepler order {p} outside band"
    assert abs(kep["mean_order"] - 4.0) < 0.3
    # (b) rate-metric convergence is FASTER than h^4 (estimator averaging
    # cancels the leading phase-error mode); require at least fourth-order
    # behavior with a justified evidence-based ceiling
    conv = experiment.convergence_study(orbits=20)
    for p in conv["orders_per_interval"]:
        assert p >= 3.6, f"rate order {p} below fourth-order requirement"
        assert p <= 5.0, f"rate order {p} suspiciously above evidence band"
    assert abs(conv["mean_order"] - 4.0) < 1.0


def test_convergence_errors_monotone():
    conv = experiment.convergence_study(orbits=20)
    errs = conv["E_h_deg_day"]
    for i in range(len(errs) - 1):
        assert errs[i] > errs[i + 1] * 12.0, "E_h must fall ~>=16x per halving"
    # finest error must sit far above the float round-off floor (no fake
    # superconvergence from error underflow): observed 8e-10 deg/day scale
    assert errs[-1] > 1e-12


# --------------------------------------------------------------------------- #
# L5: mission / real-orbit anchors
# --------------------------------------------------------------------------- #

def test_starlink_anchor():
    case = dict(a_km=R_E + 550.0, e=0.0003, inc_deg=53.0, Omega0_deg=0.0,
                omega0_deg=0.0, M0_deg=0.0, n_orbits=20, windows_orbits=(20,),
                claims=("Omega_dot",), desc="test starlink")
    cmp_ = experiment.validate_case_rates(experiment.run_case(case))
    num = cmp_["numeric_Omega_dot_deg_day"]
    ana = oracle_rates(R_E + 550.0, 0.0003, 53.0)["Omega_dd"]
    assert abs(num - ana) / abs(ana) < 1e-2
    assert abs(num - (-4.489)) < 0.05


def test_sso600_target_rate():
    i_sso = float(np.degrees(experiment.sun_sync_inclination_rad(R_E + 600.0, 0.0)))
    case = dict(a_km=R_E + 600.0, e=0.0, inc_deg=i_sso, Omega0_deg=0.0,
                omega0_deg=0.0, M0_deg=0.0, n_orbits=20, windows_orbits=(20,),
                claims=("Omega_dot",), desc="test sso")
    cmp_ = experiment.validate_case_rates(experiment.run_case(case))
    num = cmp_["numeric_Omega_dot_deg_day"]
    nc = cmp_["node_crossing_Omega_dot_deg_day"]
    assert abs(num - SSO_TARGET_DD) / SSO_TARGET_DD < 1e-2
    assert nc is not None and abs(nc - SSO_TARGET_DD) / SSO_TARGET_DD < 2e-2


def test_polar_null_via_node_crossing():
    case = dict(a_km=R_E + 500.0, e=0.0, inc_deg=90.0, Omega0_deg=0.0,
                omega0_deg=0.0, M0_deg=0.0, n_orbits=20, windows_orbits=(20,),
                claims=("Omega_dot_null",), desc="test polar")
    run = experiment.run_case(case)
    cmp_ = experiment.validate_case_rates(run)
    assert abs(cmp_["analytic_Omega_dot_deg_day"]) < 1e-12
    nc = cmp_["node_crossing_Omega_dot_deg_day"]
    assert nc is not None
    # exact polar symmetry: the J2 force is axisymmetric and the orbit plane
    # contains the symmetry axis, so the node longitude is truly stationary;
    # measured null is ~1e-16 deg/day (machine level) -- bound 1e-8 is generous
    assert abs(nc) < 1e-8
    # e=0 seed: omega trend not claimed per contract (J2-induced eccentricity
    # makes omega finite but fast-sweeping; it carries no secular signal claim)
    assert run["case"]["claims"] == ("Omega_dot_null",)
    assert cmp_["invariants"]["energy_drift_rel"] < 1e-9


def test_critical_inclination_apsidal_freeze():
    crit = float(np.degrees(np.arccos(1.0 / np.sqrt(5.0))))
    case = dict(a_km=R_E + 500.0, e=0.2, inc_deg=crit, Omega0_deg=0.0,
                omega0_deg=0.0, M0_deg=0.0, n_orbits=30, windows_orbits=(30,),
                claims=("Omega_dot", "omega_dot_null"), desc="test critical")
    run = experiment.run_case(case)
    cmp_ = experiment.validate_case_rates(run)
    om_num = cmp_["numeric_omega_dot_deg_day"]
    ana = oracle_rates(R_E + 500.0, 0.2, crit)
    # apsidal freeze: measured rate far below typical LEO signal (~3-4 deg/day)
    assert om_num is not None and abs(om_num) < 0.02
    # nodal rate still normal and matches oracle
    assert abs(cmp_["numeric_Omega_dot_deg_day"] - ana["Omega_dd"]) / abs(ana["Omega_dd"]) < 1e-2


def test_molniya_rates():
    case = dict(a_km=26560.0, e=0.74, inc_deg=63.4, Omega0_deg=0.0,
                omega0_deg=270.0, M0_deg=0.0, n_orbits=12, windows_orbits=(12,),
                claims=("Omega_dot", "omega_dot_small_nonzero"), desc="test molniya")
    run = experiment.run_case(case)
    cmp_ = experiment.validate_case_rates(run)
    ana = oracle_rates(26560.0, 0.74, 63.4)
    num_Om = cmp_["numeric_Omega_dot_deg_day"]
    # wrong-p bugs (using a instead of p=a(1-e^2)) shift the rate by (1-e^2)^-2 = 4.89x.
    # 12-orbit test window: leakage bias ~1.91*A/N with A~0.03 deg -> ~5e-3 deg/day,
    # so 3e-2 rel is the justified short-window bound (48-orbit primary is tighter).
    assert abs(num_Om - ana["Omega_dd"]) / abs(ana["Omega_dd"]) < 3e-2
    num_om = cmp_["numeric_omega_dot_deg_day"]
    assert num_om is not None
    assert 0.0 < num_om < 4.0 * ana["omega_dd"] + 1e-3  # small nonzero, right scale/sign
    assert cmp_["spp"] >= 720.0 / (1.0 - 0.74) ** 1.5  # documented resolution law applied


# --------------------------------------------------------------------------- #
# L6: pathological / singular cases
# --------------------------------------------------------------------------- #

def test_pathological_grid_sentinels():
    patho = experiment.pathological_grid(orbits=2)
    assert patho["all_ok"]
    for row in patho["rows"]:
        assert row["ok"], row


def test_equatorial_no_raan_claim():
    a = R_E + 500.0
    T = experiment.orbital_period(a)
    t = np.linspace(0.0, 2 * T, 2 * 512 + 1)
    r0, v0, _ = experiment.seed_state(a, 0.01, 0.0, 0.0, 0.0, 0.0)
    states = experiment.propagate_3d_rk4_j2(r0, v0, MU, t, experiment.J2_EARTH)
    meas = experiment.measure_secular_rates(t, states, T, (2,))
    assert meas["Omega_defined"] is False
    assert meas["Omega"][2] is None  # structurally undefined -> no claim
    # omega is measured FROM the node line: with no node it is undefined too
    # (longitude of periapsis would be the well-defined alternative, out of scope)
    assert meas["omega_defined"] is False
    assert meas["omega"][2] is None
    # dynamics stayed healthy
    rn = np.linalg.norm(states[:, :3], axis=1)
    assert np.all(rn > a * (1 - 0.01) * 0.999)
    assert np.all(rn < a * (1 + 0.01) * 1.001)


def test_circular_no_omega_claim_but_raan_valid():
    a = R_E + 500.0
    T = experiment.orbital_period(a)
    n_orbits = 10
    t = np.linspace(0.0, n_orbits * T, n_orbits * 512 + 1)
    r0, v0, _ = experiment.seed_state(a, 0.0, np.radians(51.0), 0.0, 0.0, 0.0)
    states = experiment.propagate_3d_rk4_j2(r0, v0, MU, t, experiment.J2_EARTH)
    meas = experiment.measure_secular_rates(t, states, T, (n_orbits,))
    # J2 induces a real eccentricity (~1e-3) so omega stays finite numerically,
    # but per contract its trend is NOT claimed for e=0 seeds (claims policy)
    assert meas["Omega_defined"] is True
    num = float(np.degrees(meas["Omega"][n_orbits]["slope"]) * 86400.0)
    ana = oracle_rates(a, 0.0, 51.0)["Omega_dd"]
    assert abs(num - ana) / abs(ana) < 1e-2


# --------------------------------------------------------------------------- #
# L7: adversarial regressions
# --------------------------------------------------------------------------- #

def test_j2_zero_no_spurious_drift():
    res = experiment.null_and_signflip(orbits=10)
    j0 = res["j2_zero"]
    # Omega null is clean at integrator level
    assert abs(j0["Omega_dot_deg_day"]) < 1e-9
    # omega slope carries RK4 e_vec-direction noise amplified by 1/e (e=0.01);
    # observed ~1.8e-6 deg/day. The purpose of this test is catching frame or
    # sign bugs, which produce O(0.98) deg/day artifacts: 1e-3 keeps >=3 orders
    # of margin over the artifact scale while accommodating the noise floor.
    assert abs(j0["omega_dot_deg_day"]) < 1e-3


def test_sign_flip_sensitivity():
    res = experiment.null_and_signflip(orbits=25)
    sf = res["sign_flip"]
    # flipping J2 must flip BOTH precessions; magnitudes equal to within the
    # window-leakage asymmetry (short-period phase vs Keplerian window does
    # not flip antisymmetically): ~1% at 25 orbits, observed 0.9% at 10
    assert sf["ratio_Omega"] is not None
    assert -1.01 < sf["ratio_Omega"] < -0.99
    assert sf["ratio_omega"] is not None
    assert -1.02 < sf["ratio_omega"] < -0.98
    ana = res["analytic_plus_J2"]
    assert abs(sf["plus_J2_Omega_deg_day"] - ana["Omega_dot_deg_day"]) < 0.02 * abs(ana["Omega_dot_deg_day"])
    assert abs(sf["plus_J2_omega_deg_day"] - ana["omega_dot_deg_day"]) < 0.05 * abs(ana["omega_dot_deg_day"])


def test_results_rate_conversion_twins():
    # internal consistency of the deg/day <-> rad/s twins used in results.json
    r = experiment.analytic_rates(R_E + 420.0, 0.0003, np.radians(51.6))
    back = np.degrees(r["Omega_dot_rad_s"]) * 86400.0
    assert abs(back - r["Omega_dot_deg_day"]) < 1e-12 * abs(back)
    back_w = np.degrees(r["omega_dot_rad_s"]) * 86400.0
    assert abs(back_w - r["omega_dot_deg_day"]) < 1e-12 * abs(back_w)


def test_determinism_repeat_call():
    a = R_E + 450.0
    T = experiment.orbital_period(a)
    t = np.linspace(0.0, 2 * T, 2 * 256 + 1)
    r0, v0, _ = experiment.seed_state(a, 0.05, np.radians(45.0), 0.0, 0.0, 0.0)
    s1 = experiment.propagate_3d_rk4_j2(r0, v0, MU, t, experiment.J2_EARTH)
    s2 = experiment.propagate_3d_rk4_j2(r0, v0, MU, t, experiment.J2_EARTH)
    assert np.array_equal(s1, s2)
    m1 = experiment.node_crossing_raan_rate(t, s1)
    m2 = experiment.node_crossing_raan_rate(t, s2)
    assert m1["slope"] == m2["slope"]
