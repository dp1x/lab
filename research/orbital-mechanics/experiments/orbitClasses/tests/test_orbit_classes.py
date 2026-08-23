"""Experiment 012 - Orbit Classes: validation suite.

Layer map
---------
L1  closed-form identities & convention firewalls (constants duplicated inline on
    purpose -- tests never import expected values from the module under test;
    anti-shared-algebra doctrine, Exp 009/011 precedent)
L2  numerical recovery: propagation-based confirmation of class constraints
L3  convergence, invariants & determinism
L4  adversarial mutant battery: every realistic wrong implementation listed by the
    adversarial track must be caught by a named discriminator; known-blind test
    families are pre-registered and their compensating checks enforced
L5  committed-artifact integrity (results.json / figures)

Doctrine notes
--------------
* omega_E is the MASTER sidereal constant; P_sidereal is derived, never hard-coded.
* The solar-vs-sidereal-year discrimination (3.0e-4 deg in i_SSO) and the binding
  solar-vs-Julian-year pair (1.67e-4 deg) are ANALYTIC-layer assertions: they are
  far below any honest propagation band. The tropical-year variant (2.1e-7 deg)
  is behaviorally indistinguishable -- documented blindness, pinned by constant
  literal instead (pre-registered survivor M02b).
* Near the critical inclination the FULL J2 problem carries short-period dynamics
  (measured: osculating-a excursions ~ +160 km, event-period excess ~ +325 s/orbit,
  energy-conserving and converged). Millisecond period-split detection is therefore
  NOT claimed; the claimable freeze evidence lives on the orbit-averaged
  element-regression path (Exp 009 doctrine).
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("orbit_classes_experiment", EXPERIMENT_DIR / "experiment.py")
assert _spec is not None and _spec.loader is not None
EXP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(EXP)

# --- inline-duplicated constants (provenance: same sources as lab_utils canon) -
MU = 398600.4418          # IAU 2015 Resolution B3 nominal GM_E [km^3/s^2]
RE = 6378.137             # WGS-84 TR8350.2 equatorial radius [km]
J2 = 1.082629821e-3       # WGS-84 sqrt(5)|C20_bar|, C20_bar = -0.484166774985e-3
WE = 7.2921159e-5         # WGS-84 / Vallado Table 3-1 [rad/s] (MASTER sidereal constant)
YEAR_SOLAR = 365.2422     # mean-solar/Gregorian-mean convention d
SUN_RATE = 360.0 / YEAR_SOLAR  # deg/day
PSID = 2.0 * np.pi / WE   # sidereal day [s], derived from omega_E


def _deg_day(rad_s: float) -> float:
    return float(np.degrees(rad_s)) * 86400.0


def _oracle_omega_dot(a: float, e: float, inc_rad: float) -> float:
    """First-order secular nodal rate [deg/day]; independent inline duplicate."""
    n = math.sqrt(MU / a ** 3)
    p = a * (1.0 - e * e)
    return _deg_day(-1.5 * n * J2 * (RE / p) ** 2 * math.cos(inc_rad))


def _oracle_sso_cos(a: float, e: float) -> float:
    """cos i_SSO (retrograde branch, negative); inline duplicate of the closed form."""
    n = math.sqrt(MU / a ** 3)
    p = a * (1.0 - e * e)
    tgt_rad_s = math.radians(SUN_RATE) / 86400.0
    return -tgt_rad_s * p * p / (1.5 * n * J2 * RE ** 2)


# --- cached heavier propagations ---------------------------------------------
_RUNS: dict = {}


def _sso600_run():
    if "sso600" not in _RUNS:
        a = RE + 600.0
        sol = EXP.solve_sso_inclination(a, 0.0)
        _RUNS["sso600"] = EXP.propagate_case(a, 0.0, sol["incl_rad"], 40, 512)
    return _RUNS["sso600"]


def _lock_run():
    if "lock" not in _RUNS:
        _RUNS["lock"] = EXP.propagate_case(
            EXP.A_SEMISYNC_KM, 0.74, np.radians(EXP.CRITICAL_INC_DEG), 12,
            EXP.steps_per_orbit(0.74))
    return _RUNS["lock"]


def _drift_run():
    if "drift" not in _RUNS:
        _RUNS["drift"] = EXP.propagate_case(
            EXP.A_SEMISYNC_KM, 0.74, np.radians(EXP.CRITICAL_INC_DEG), 13,
            EXP.steps_per_orbit(0.74))
    return _RUNS["drift"]


_BLOCK_CACHE: dict = {}


def _block(name: str, builder):
    """Memoize the heavy analysis blocks -- each runs at most once per session."""
    if name not in _BLOCK_CACHE:
        _BLOCK_CACHE[name] = builder()
    return _BLOCK_CACHE[name]


def _fam():
    return _block("fam", EXP.molniya_family)


def _per():
    return _block("per", EXP.molniya_periods_and_drift)


def _crit():
    return _block("crit", EXP.critical_inclination_sweep)


def _gto():
    return _block("gto", EXP.gto_analysis)


def _geo():
    return _block("geo", EXP.geo_anchors)


def _conv():
    return _block("conv", EXP.convergence_study)


# ============================================================================ #
# L1 -- closed-form identities & conventions
# ============================================================================ #
def test_sso_solver_anchors_and_branch():
    for h, expected in ((500.0, 97.401785943095), (600.0, 97.787646791197),
                        (800.0, 98.603085267154)):
        sol = EXP.solve_sso_inclination(RE + h, 0.0)
        assert sol["status"] == "OK"
        i_deg = math.degrees(sol["incl_rad"])
        assert abs(i_deg - expected) <= 5e-5, (h, i_deg)
        # retrograde branch selection: i strictly above 90 deg
        assert 90.0 < i_deg < 180.0
        # inline oracle agreement (anti-shared-algebra twin)
        assert abs(i_deg - math.degrees(math.acos(_oracle_sso_cos(RE + h, 0.0)))) < 1e-9


def test_sso_wrong_branch_gives_exact_negated_rate():
    a = RE + 600.0
    sol = EXP.solve_sso_inclination(a, 0.0)
    mirror = math.pi - sol["incl_rad"]
    got = _oracle_omega_dot(a, 0.0, mirror)
    want = _oracle_omega_dot(a, 0.0, sol["incl_rad"])
    assert abs(want - SUN_RATE) < 1e-9          # solution produces +sun rate
    assert abs(got + SUN_RATE) < 1e-9           # mirror produces exactly -sun rate
    assert got < 0.0 < want                     # signed discriminators, not |.| ones


def test_cos_identity_neg_apow72():
    a_max = EXP.sso_existence_max_sma(0.0)
    for h in (400.0, 800.0, 1200.0):
        a = RE + h
        sol = EXP.solve_sso_inclination(a, 0.0)
        lhs = -math.cos(sol["incl_rad"])
        rhs = (a / a_max) ** 3.5
        assert abs(lhs - rhs) <= 1e-12 * abs(rhs)


def test_sso_family_monotone_in_altitude():
    incs = []
    for h in range(400, 1401, 100):
        sol = EXP.solve_sso_inclination(RE + float(h), 0.0)
        assert sol["status"] == "OK"
        incs.append(math.degrees(sol["incl_rad"]))
    assert all(b > a for a, b in zip(incs[:-1], incs[1:]))


def test_eccentric_sso_shifts_inclination_up():
    rp = RE + 600.0
    e = 0.2
    sol = EXP.solve_sso_inclination(rp / (1.0 - e), e)
    assert sol["status"] == "OK"
    i_ecc = math.degrees(sol["incl_rad"])
    assert abs(i_ecc - 105.82454538603649) <= 1e-4
    # vs the circular SSO at the SAME PERIGEE radius: eccentricity (higher apogee,
    # weaker J2) demands a stronger tilt
    i_circ_same_rp = math.degrees(EXP.solve_sso_inclination(rp, 0.0)["incl_rad"])
    assert i_ecc > i_circ_same_rp
    # vs the circular SSO at the SAME semi-major axis: smaller p strengthens J2,
    # so the required tilt is REDUCED (both directions documented)
    a = rp / (1.0 - e)
    i_circ_same_a = math.degrees(EXP.solve_sso_inclination(a, 0.0)["incl_rad"])
    assert i_ecc < i_circ_same_a


def test_existence_limit_value_sentinel_and_no_silent_clip():
    a_max = EXP.sso_existence_max_sma(0.0)
    assert abs(a_max - 12352.505076188283) <= 1e-6
    below = EXP.solve_sso_inclination(a_max * 0.999999, 0.0)
    assert below["status"] == "OK"
    assert math.degrees(below["incl_rad"]) > 179.5       # approaches 180, not 90
    nearer = EXP.solve_sso_inclination(a_max * (1.0 - 1e-10), 0.0)
    assert math.degrees(nearer["incl_rad"]) > math.degrees(below["incl_rad"])  # monotone approach
    above = EXP.solve_sso_inclination(a_max * 1.000001, 0.0)
    assert above["status"] == "NO_REAL_SOLUTION"         # typed sentinel, never clipped
    assert math.isnan(above["incl_rad"])
    # eccentricity EXTENDS the limit: a_max(e) = a_max(0) (1-e^2)^(-4/7)
    e = 0.2
    assert abs(EXP.sso_existence_max_sma(e)
               - a_max * (1.0 - e * e) ** (-4.0 / 7.0)) <= 1e-9


def test_year_convention_discriminators_analytic_layer():
    a = RE + 600.0
    base = math.degrees(EXP.solve_sso_inclination(a, 0.0)["incl_rad"])
    sid = math.degrees(EXP.solve_sso_inclination(a, 0.0, 360.0 / 365.25636)["incl_rad"])
    jul = math.degrees(EXP.solve_sso_inclination(a, 0.0, 360.0 / 365.25)["incl_rad"])
    trop = math.degrees(EXP.solve_sso_inclination(a, 0.0, 360.0 / 365.24219)["incl_rad"])
    assert abs(base - sid) >= 2.0 * EXP.I_SSO_TOL_DEG     # sidereal year rejected by tolerance
    assert abs(base - jul) >= EXP.I_SSO_TOL_DEG           # Julian year rejected (binding pair)
    assert abs(base - trop) < 1e-6                        # tropical: documented blindness (M02b)
    assert abs(SUN_RATE - 0.9856473320990837) < 1e-15     # target constant literal pin


def test_earth_rotation_confusion_guard():
    # confusing omega_E (per-day 360.98565 deg/day) with the sun rate must trip the
    # domain guard, not return a number
    status = EXP.solve_sso_inclination(RE + 600.0, 0.0,
                                       float(np.degrees(WE)) * 86400.0)["status"]
    assert status == "NO_REAL_SOLUTION"


def test_critical_inclination_exact_supplement_and_identity():
    icrit = EXP.CRITICAL_INC_DEG
    assert abs(icrit - 63.43494882292201) <= 1e-12
    assert abs(EXP.CRITICAL_INC_SUPP_DEG - 116.56505117707799) <= 1e-12
    assert abs(math.cos(np.radians(icrit)) ** 2 - 0.2) <= 1e-15


def test_semisync_radius_resonance_identity():
    a = EXP.A_SEMISYNC_KM
    assert abs(a - 26561.762328167155) <= 1e-9
    n = math.sqrt(MU / a ** 3)
    assert abs(n - 2.0 * WE) <= 1e-15 * (2.0 * WE)                 # n = 2 omega_E
    assert abs(2.0 * np.pi / n - PSID / 2.0) <= 1e-9               # T = P_sidereal/2
    # canonical Molniya geometry at e = 0.74
    assert abs(a * 0.26 - RE - 527.921205323461) <= 1e-6
    assert abs(a * 1.74 - RE - 39839.329451010846) <= 1e-6


def test_geo_radius_period_identity_and_altitude_anchor():
    a = EXP.A_GEO_KM
    assert abs(a - 42164.169461861835) <= 1e-9
    T = 2.0 * math.pi * math.sqrt(a ** 3 / MU)
    assert abs(T - PSID) <= 1e-12 * PSID                            # construction identity
    assert abs((a - RE) - 35786.03246186183) <= 1e-6                # classic altitude anchor


def test_molniya_off_lock_residual_is_finite_not_zero():
    # 63.4 deg is NOT the exact lock: residual must be the documented +4.0e-4 scale
    ana = EXP.analytic_rates(EXP.A_SEMISYNC_KM, 0.74, np.radians(63.4))
    assert 3.0e-4 < ana["omega_dot_deg_day"] < 5.0e-4


def test_mean_anomaly_rate_bracket_distinct_from_apsidal():
    # M-dot bracket (3cos^2 i - 1) differs from omega_dot's (5cos^2 i - 1)
    a, inc = RE + 500.0, np.radians(30.0)
    n = math.sqrt(MU / a ** 3)
    x = J2 * (RE / a) ** 2
    mdot_minus_n = EXP.mean_anomaly_rate(a, 0.0, inc) - n
    assert abs(mdot_minus_n - 0.75 * n * x * (3.0 * math.cos(inc) ** 2 - 1.0)) <= 1e-18
    # sign structure: at i=0 both positive; at polar i=90 M-dot correction negative
    assert EXP.mean_anomaly_rate(a, 0.0, 0.0) > n
    assert EXP.mean_anomaly_rate(a, 0.0, np.radians(90.0)) < n
    # ECCENTRIC absolute pin (e > 0 is the ONLY regime where a p:=a substitution is
    # visible; at e = 0 p == a identically and every e=0 check is structurally blind):
    # at Molniya elements / exact lock, sqrt(1-e^2)*(3cos^2 i - 1) =
    # 0.672607*(-0.4); M_dot - n must equal the inline p-based value.
    am, em = EXP.A_SEMISYNC_KM, 0.74
    nm = math.sqrt(MU / am ** 3)
    pm = am * (1.0 - em * em)
    icrit_rad = math.acos(1.0 / math.sqrt(5.0))
    mdot_mol = EXP.mean_anomaly_rate(am, em, icrit_rad)
    expected = nm + 0.75 * nm * J2 * (RE / pm) ** 2 * math.sqrt(1.0 - em * em) * (3.0 * math.cos(icrit_rad) ** 2 - 1.0)
    assert abs(mdot_mol - expected) <= 1e-15
    # the p:=a mutant shifts this by a factor (1-e^2)^-2 = 4.888 -- far outside
    distorted = nm + 0.75 * nm * J2 * (RE / am) ** 2 * math.sqrt(1.0 - em * em) * (3.0 * math.cos(icrit_rad) ** 2 - 1.0)
    assert abs(mdot_mol - distorted) > 1e-9


def test_gto_budget_signed_burns_ordering_and_swap_trap():
    gto = _gto()
    row = gto["rows"][2]  # hp = 300 km canonical
    assert abs(row["dv1_km_s"] - 2.425732698207) <= 1e-9
    assert abs(row["dv2_km_s"] - 1.466824319545) <= 1e-9
    assert abs(row["dv_total_km_s"] - 3.892557017752) <= 1e-9
    assert row["dv1_km_s"] > row["dv2_km_s"]                       # raising burn is larger
    vp, va = row["vp_km_s"], row["va_km_s"]
    vcrp, vcra = row["vcirc_rp_km_s"], row["vcirc_ra_km_s"]
    assert vp > vcrp > vcra > va                                   # speed ordering chain
    rp, ra = RE + 300.0, EXP.A_GEO_KM
    assert abs(vp * rp / (va * ra) - 1.0) <= 1e-12                 # h conservation (dual form)
    # swap-mutant trap: total is invariant under dv1<->dv2 -- total-only asserts are
    # structurally blind; per-burn signed asserts above are the catch layer
    swapped = row["dv2_km_s"] + row["dv1_km_s"]
    assert abs(swapped - row["dv_total_km_s"]) <= 1e-12


def test_gto_exp004_anchor_continuity():
    anchor = _gto()["exp004_anchor"]
    assert abs(anchor["dv_total_km_s"] - 3.931859594283) <= 1e-9
    assert abs(anchor["dv_total_km_s"] - 3.9319) <= 1e-4           # Exp 004 recorded value


def test_gto_geometry_bugs_quantified():
    bugs = _gto()["bugs"]
    # bug X: GEO altitude used as radius -> budget short by ~122.6 m/s, final orbit
    # circularized at the wrong radius has nowhere near the synchronous period
    bx = bugs["bug_X_altitude_as_radius"]
    assert abs(bx["dv_total_km_s"] - 3.76999) <= 2e-3
    # circularizing at the WRONG apogee radius gives an 18.71 h orbit -- hours away
    # from the 23.93 h synchronous day the mission needs
    assert abs(bx["circular_period_at_wrong_apogee_hours"] - 18.71452493290484) <= 1e-6
    assert abs(bx["circular_period_at_wrong_apogee_hours"] - 23.9345) > 3.0
    # bug Y: GEO radius used as altitude -> budget overshoots by ~82 m/s
    assert abs(bugs["bug_Y_radius_as_altitude"]["dv_total_km_s"] - 3.97441) <= 2e-3


def test_circular_speed_on_ellipse_is_caught_by_vis_viva_identity():
    errs = _gto()["bugs"]["circular_v_mutant_errors_km_s"]
    # v = sqrt(mu/a) misestimates apsidal speeds by > 1 km/s at Molniya-class e
    assert abs(errs["at_molniya_perigee"]) > 1.0
    assert abs(errs["at_molniya_apogee"]) > 1.0
    # vis-viva identity holds at sampled anomalies of a true ellipse (the checker)
    a, e = EXP.A_SEMISYNC_KM, 0.74
    for nu in np.linspace(0.0, 2.0 * np.pi, 25):
        r = a * (1.0 - e * e) / (1.0 + e * math.cos(nu))
        v = math.sqrt(MU * (2.0 / r - 1.0 / a))
        assert abs(v - math.sqrt(MU * (1.0 + 2.0 * e * math.cos(nu) + e * e) / (a * (1.0 - e * e)))) \
            <= 1e-9


def test_dwell_closed_form_values_limits_monotonicity():
    for delta, expected in ((30.0, 0.6064438460305431), (60.0, 0.8421914604381311),
                            (90.0, 0.9236066186864406)):
        assert abs(EXP.dwell_fraction_closed_form(0.74, delta) - expected) <= 1e-12
    # limit: for a +/-90 deg window f -> 1/2 as e -> 0 (half the orbit either side);
    # narrower windows tend to their circular-geometry value (e.g. +/-60 deg -> 1/3)
    assert abs(EXP.dwell_fraction_closed_form(1e-12, 90.0) - 0.5) <= 1e-9
    assert abs(EXP.dwell_fraction_closed_form(1e-12, 60.0) - 1.0 / 3.0) <= 1e-9
    # full window (+/-180 deg) always covers the whole orbit
    for e in (0.3, 0.6, 0.74):
        assert abs(EXP.dwell_fraction_closed_form(e, 180.0) - 1.0) <= 1e-12
    # monotone growth of the apogee dwell with eccentricity
    fracs = [EXP.dwell_fraction_closed_form(e, 60.0) for e in (0.0, 0.3, 0.6, 0.74)]
    assert all(b > a for a, b in zip(fracs[:-1], fracs[1:]))


# ============================================================================ #
# L2 -- numerical recovery
# ============================================================================ #
def test_sso_numeric_closure_at_solved_inclination():
    run = _sso600_run()
    meas = EXP.measured_secular_trends(run, (20, 40))
    numA = meas["Omega_dot_deg_day"]
    numB = meas["node_crossing_Omega_dot_deg_day"]
    relA = abs(numA - SUN_RATE) / SUN_RATE
    assert relA <= 1.0e-2                    # model-order disclosure band (measured 0.6%)
    assert numA > 0.0                        # eastward: sign discriminator (branch mutants)
    assert meas["node_crossing_Omega_dot_deg_day"] > 0.0
    assert abs(numA - numB) / abs(numA) <= 1.0e-2   # dual-path (leakage-scale disclosed)
    inv = EXP.invariants_gate(run)
    assert inv["energy_drift_rel"] < 1e-8


def test_molniya_freeze_on_element_regression_path():
    run = _lock_run()
    meas = EXP.measured_secular_trends(run, (6, 12))
    assert abs(meas["omega_dot_deg_day"]) <= EXP.OMEGA_NULL_ABS_DEG_DAY  # null-consistent
    nc = meas["node_crossing_Omega_dot_deg_day"]
    ana = EXP.analytic_rates(EXP.A_SEMISYNC_KM, 0.74, np.radians(EXP.CRITICAL_INC_DEG))
    assert abs(nc - ana["Omega_dot_deg_day"]) / abs(ana["Omega_dot_deg_day"]) <= 2.0e-2
    inv = EXP.invariants_gate(run)
    assert inv["energy_drift_rel"] < 1e-8


def test_kepler_machinery_check_bit_tight():
    chk = _per()["kepler_machinery_check"]
    assert chk["T_peri_rel_err_vs_kepler"] <= 1e-9
    assert chk["T_node_rel_err_vs_kepler"] <= 1e-9


def test_j2on_kepler_excess_reported_with_plateau():
    near = _per()["near_critical_63p4"]
    ex = near["kepler_excess_peri_ms"]
    assert 250000.0 < ex < 400000.0          # measured excess ~ +323 s/orbit
    ratio = near["convergence_ratio_excess"]
    assert EXP.KEPLER_EXCESS_PLATEAU_RATIO_BAND[0] <= ratio <= EXP.KEPLER_EXCESS_PLATEAU_RATIO_BAND[1]
    # first-order split is DISCLOSED, never claimed as detected
    assert abs(near["split_ms_first_order_theory_disclosure"] - 24.059563336778) <= 1e-6
    assert "NOT claimed as detected" in near["split_detection_note"]


def test_critical_sweep_localizes_sign_change():
    sw = _crit()
    zc = sw["zero_crossing_localized_deg"]
    assert zc is not None
    assert abs(zc - EXP.CRITICAL_INC_DEG) <= 0.15
    # antisymmetry about the lock at +/-0.5 deg
    assert sw["antisym_ratio_pm05deg"] <= EXP.ANTISYM_RATIO_TOL
    # slope vs theory (finite-step side, +/-10%)
    theory = sw["slope_theory_deg_day_per_deg"]
    assert abs(sw["slope_fit_deg_day_per_deg"] - theory) <= 0.10 * abs(theory)
    # plateau under step halving at the lock
    pl = sw["plateau_at_lock_by_spp"]
    vals = [abs(v) for v in pl.values()]
    assert max(vals) <= 5.0 * max(min(vals), 1e-6)


def test_dwell_numeric_matches_closed_form():
    fam = _fam()
    for row in fam["rows"]:
        assert row["dwell_abs_err"] <= EXP.DWELL_NUM_TOL_ABS
    row74 = next(r for r in fam["rows"] if abs(r["e"] - 0.74) < 1e-12)
    assert abs(row74["dwell90_closed"] - 0.9236066186864406) <= 1e-9


def test_p_scaling_law_of_apsidal_rate():
    fam = _fam()
    # omega_dot prop-to (R/p)^2 => log-log slope vs (1-e^2)^-2 equals 1
    assert abs(fam["p_scaling_loglog_slope"] - 1.0) <= 0.02


def test_apogee_event_rate_identity():
    drift = EXP.apogee_longitude_drift(_drift_run()["t"], _drift_run()["states"])
    assert drift["n_apogees"] >= 12
    icrit = np.radians(EXP.CRITICAL_INC_DEG)
    ana = EXP.analytic_rates(EXP.A_SEMISYNC_KM, 0.74, icrit)
    inertial = ana["Omega_dot_deg_day"] + ana["omega_dot_deg_day"] - float(np.degrees(WE)) * 86400.0
    aps = EXP.apsis_crossings(_drift_run()["t"], _drift_run()["states"])
    T_orb = EXP.slope_vs_index(aps["peri"])
    pred = 360.0 * 86400.0 / T_orb + inertial
    assert abs(drift["event_rate_deg_day"] - pred) <= 0.05
    assert drift["fit_residual_rms_deg"] <= 0.1
    # mean-line disclosure value recorded for the record
    res = json.loads((EXPERIMENT_DIR / "results" / "results.json").read_text(encoding="utf-8"))
    assert abs(res["results"]["headline"]["apogee_mean_line_disclosure_deg_day"] + 1.133368081634) <= 1e-9


def test_repeat_lattice_keplerian_closure():
    lat = _per()["repeat_lattice_keplerian"]
    for entry in lat:
        assert abs(entry["closure_deg_after_m_orbits"]) <= 1e-6
    by_key = {e["m_per_k_days"]: e["a_res_km"] for e in lat}
    assert abs(by_key["1:1"] - EXP.A_GEO_KM) <= 1e-6              # GEO is the 1:1 member
    assert abs(by_key["2:1"] - EXP.A_SEMISYNC_KM) <= 1e-6         # Molniya is the 2:1 member


def test_repeat_corrected_radius_first_order_disclosure():
    per = _per()
    a_rep = per["repeat_corrected_radius_km"]
    assert abs(a_rep - 26553.420404899953) <= 1e-6
    assert a_rep < per["naive_radius_km"]                          # regression shrinks the orbit


def test_geo_stationarity_negative_control_and_inclined_nodal():
    geo = _geo()
    # nonzero stationarity residual kills any "everything vanishes at GEO" mutant
    assert abs(geo["stationarity_residual_deg_day_keplerian_a"] - 0.026828214162) <= 2e-6
    # individual rates are nonzero but tied: omega_dot = -2 Omega_dot at i = 0
    assert geo["relation_omegadot_eq_minus2Omegadot_i0_rel_err"] <= 1e-12
    # construction identity
    assert geo["period_match_rel_err"] <= 1e-12
    # inclined-GEO draconitic shift: first-order prediction vs measurement inside
    # the disclosed model-order band (short-period jitter floor ~ seconds)
    assert abs(geo["inclined_i5_T_node_pred_minus_sid_ms"] + 9507.2) <= 2.0
    assert abs(geo["inclined_i5_T_node_num_minus_sid_ms"]
               - geo["inclined_i5_T_node_pred_minus_sid_ms"]) <= 1000.0 * EXP.GEO_INCLINED_NODAL_TOL_S


def test_gto_rk4_flight_arrival():
    flight = _gto()["canonical_flight"]["arrival"]
    assert flight["kepler_only"]["arrival_rel_err"] <= 1e-6
    # J2-on coast deviates at model-order level (disclosed, not hidden)
    assert 0.0 < flight["j2_on"]["arrival_rel_err"] <= 1.0e-2


# ============================================================================ #
# L3 -- convergence, invariants & determinism
# ============================================================================ #
def test_state_space_convergence_order_band():
    conv = _conv()
    assert 3.6 <= conv["state_space"]["mean_order"] <= 4.4
    assert all(3.6 <= o <= 4.4 for o in conv["state_space"]["orders_per_interval"])


def test_rate_metric_convergence_rule():
    conv = _conv()
    orders = conv["rate_metric"]["orders_per_interval"]
    assert all(o >= 3.6 for o in orders)      # fourth-order-or-better per interval
    assert all(o <= 5.0 for o in orders)
    errs = conv["rate_metric"]["errors_deg_day"]
    assert errs[-1] < errs[0]                 # monotone decay overall


def test_period_ladder_report_only_floor():
    # near i_crit the event-jitter floor dominates (~0.5 s); recorded, not claimed
    conv = _conv()
    assert conv["period_metric"]["finest_gap_ms"] <= 2000.0
    assert len(conv["period_metric"]["errors_s"]) == 3


def test_pathological_grid_all_ok():
    patho = EXP.pathological_grid()
    assert patho["all_ok"]
    assert len(patho["rows"]) == 20


def test_determinism_repeat_call_bitwise():
    sol1 = EXP.solve_sso_inclination(RE + 600.0, 0.0)
    sol2 = EXP.solve_sso_inclination(RE + 600.0, 0.0)
    assert sol1 == sol2
    run_a = EXP.propagate_case(RE + 500.0, 0.001, np.radians(51.6), 2, 512)
    run_b = EXP.propagate_case(RE + 500.0, 0.001, np.radians(51.6), 2, 512)
    assert np.array_equal(run_a["states"], run_b["states"])
    d1 = EXP.dwell_fraction_numeric(run_a["t"], run_a["states"], 60.0, run_a["T_kepler_s"])
    d2 = EXP.dwell_fraction_numeric(run_b["t"], run_b["states"], 60.0, run_b["T_kepler_s"])
    assert d1 == d2


# ============================================================================ #
# L4 -- adversarial battery & firewalls
# ============================================================================ #
def test_adversarial_battery_records_detectable_symptoms():
    adv = EXP.adversarial_battery()
    flip = adv["j2_sign_flip_wrong_branch"]
    assert abs(flip["i_wrong_branch_deg"] - 82.21235320880295) <= 1e-6   # prograde twin
    # evaluated under TRUE J2 the flipped solution produces exactly the negated rate
    assert abs(flip["produced_Omega_dot_true_J2_deg_day"] + SUN_RATE) <= 1e-9
    assert flip["produced_Omega_dot_true_J2_deg_day"] < 0.0
    solar = adv["solar_halfday_resonance_mutant"]
    assert abs(solar["delta_km"] - 48.4) <= 1.0                          # absolute-a catches
    psub = adv["p_equals_a_mutant"]
    assert abs(psub["Omega_dot_distortion_factor"] - (1.0 / 0.4524 ** 2)) <= 1e-6
    dwellm = adv["dwell_linear_mutant"]
    assert abs(dwellm["error_pp"] - 42.36) <= 0.1


def test_mutant_blind_spots_have_compensating_absolute_asserts():
    # M12d survivor class: frozen-null tests cannot see p:=a (zeros map to zeros).
    # Compensating absolute-value asserts (this test) pin the Molniya rates.
    ana = EXP.analytic_rates(EXP.A_SEMISYNC_KM, 0.74, np.radians(EXP.CRITICAL_INC_DEG))
    assert abs(ana["Omega_dot_deg_day"] - (-0.1477185028649066)) <= 1e-6
    # naive p:=a would inflate |Omega_dot| by (1-e^2)^-2 = 4.888x -- rejected here
    distorted = -0.1477185028649066 / (0.4524 ** 2)
    assert abs(distorted - ana["Omega_dot_deg_day"]) > 0.5
    # M06a survivor class: total-dv asserts cannot see swap; signed per-burn pins
    row = _gto()["rows"][2]
    assert abs(row["dv1_km_s"] - 2.425732698207) <= 1e-9
    assert abs(row["dv2_km_s"] - 1.466824319545) <= 1e-9


def test_unit_firewalls_km_m_and_deg_rad():
    # partial km->m conversion in R_eq blows up the J2 term: solver must land outside
    # the family table (altitude-trend killer, adversarial M04a)
    bad = EXP.solve_sso_inclination(RE + 600.0, 0.0, R_eq_km=RE * 1000.0)
    if bad["status"] == "OK":
        assert math.degrees(bad["incl_rad"]) < 90.0 + 1e-2   # collapses onto 90 deg, off-table
    else:
        assert bad["status"] == "NO_REAL_SOLUTION"
    # degree/radian firewall: feeding the sun rate as a rad/s-magnitude number
    crazy = EXP.solve_sso_inclination(RE + 600.0, 0.0, target_deg_day=0.9856473320990837e-6)
    assert crazy["status"] == "NO_REAL_SOLUTION" or math.degrees(crazy["incl_rad"]) < 91.0
    # sidereal-vs-solar day confusion in the rad/s conversion shifts i_SSO by ~2.2e-2 deg,
    # three orders above the 5e-5 deg tolerance (caught analytically)
    a = RE + 600.0
    n = math.sqrt(MU / a ** 3)
    K = 1.5 * n * J2 * (RE / a) ** 2
    tgt_wrong_rad_s = math.radians(SUN_RATE * 86400.0 / 86164.0905) / 86400.0
    i_conv = math.degrees(math.acos(-tgt_wrong_rad_s / K))
    i_true = math.degrees(EXP.solve_sso_inclination(a, 0.0)["incl_rad"])
    assert abs(i_conv - i_true) >= 1.5e-4


def test_unwrap_guard_policy():
    clean = np.linspace(0.0, 10.0, 500)
    assert EXP.unwrap_guard_ok(clean)
    aliased = np.array([0.0, 0.1, 0.2, 3.0, 3.1])   # jump > pi/4 between samples
    assert not EXP.unwrap_guard_ok(aliased)


def test_constants_match_canon_and_provenance_complete():
    from lab_utils import orbits as canon
    assert EXP.MU_EARTH_KM3S2 == canon.MU_EARTH_KM3S2 == MU
    assert EXP.R_EARTH_KM == canon.R_EARTH_KM == RE
    assert EXP.J2_EARTH == canon.J2_EARTH == J2
    assert EXP.OMEGA_EARTH_RAD_S == canon.OMEGA_EARTH_RAD_S == WE
    res = json.loads((EXPERIMENT_DIR / "results" / "results.json").read_text(encoding="utf-8"))
    consts = res["results"]["constants"]
    for key in ("mu_provenance", "R_E_provenance", "J2_provenance", "omega_E_provenance",
                "sso_convention", "frame_convention", "units"):
        assert isinstance(consts[key], str) and len(consts[key]) > 20, key


# ============================================================================ #
# L5 -- committed-artifact integrity
# ============================================================================ #
def _results():
    return json.loads((EXPERIMENT_DIR / "results" / "results.json").read_text(encoding="utf-8"))


def test_results_json_headline_pins():
    h = _results()["results"]["headline"]
    pins = {
        "i_SSO_600km_deg": (97.787646791197, 5e-5),
        "sso_existence_a_max_km": (12352.505076188283, 1e-6),
        "critical_inclination_deg": (63.434948822922, 1e-9),
        "zero_crossing_localized_deg": (63.42989447034195, 1e-6),
        "a_semisync_km": (26561.762328167155, 1e-9),
        "molniya_kepler_excess_peri_ms": (323005.16084601986, 1.0),
        "apogee_event_rate_measured_deg_day": (355.460235485893, 0.05),
        "repeat_corrected_radius_km": (26553.420404899953, 1e-6),
        "a_geo_km": (42164.169461861835, 1e-9),
        "gto_dv_total_300km_km_s": (3.892557017752, 1e-9),
        "state_order_mean": (4.092491779905, 0.35),
        "rate_order_mean": (4.506, 0.5),
        "worst_sso_numeric_closure_rel": (0.006143428213, 1e-3),
    }
    for key, (expected, tol) in pins.items():
        assert abs(h[key] - expected) <= tol, (key, h[key], expected)


def test_results_figures_registered_and_present():
    res = _results()["results"]
    fig_dir = EXPERIMENT_DIR / "results" / "figures"
    for name in res["figures"]:
        assert (fig_dir / name).exists() and (fig_dir / name).stat().st_size > 1000
    assert len(res["figures"]) == 6


def test_results_contract_block_discloses_dependencies():
    contract = _results()["results"]["contract"]
    assert contract["rates_kind"].startswith("SECULAR")
    assert "model-order" in contract["model_residual_policy"]
    assert abs(contract["i_sso_tolerance_deg"] - 5e-5) <= 1e-12
    assert "Julian" in contract["i_sso_tolerance_justification"]
