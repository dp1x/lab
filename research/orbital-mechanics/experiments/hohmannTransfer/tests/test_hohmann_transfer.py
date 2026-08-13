"""Validation tests for the Hohmann transfer experiment.

These must pass BEFORE any results are trusted (laboratory rule: verify before
trust). The tests check the closed-form machinery from first principles
(vis-viva, energy, angular momentum, Kepler III), the radius-ratio sweep with
its asymptotes and interior maximum, the RK4 trajectory validation of the
complete transfer, the two-impulse optimality scan, the real-system anchors
(LEO->GEO, Earth->Mars, Earth->Venus, trans-Mars injection), and determinism.

The experiment module and its dependency (Experiment 002's verified machinery)
are loaded via importlib from explicit paths (see tools/new_experiment.py).
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "hohmann_transfer_experiment", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)

AU_KM = experiment.AU_KM
GEO_ALT_KM = experiment.GEO_ALT_KM
LEO_ALT_KM = experiment.LEO_ALT_KM
MARS_A_AU = experiment.MARS_A_AU
MU = experiment.MU
MU_EARTH_KM3S2 = experiment.MU_EARTH_KM3S2
MU_SUN_KM3S2 = experiment.MU_SUN_KM3S2
R_EARTH_KM = experiment.R_EARTH_KM
VENUS_A_AU = experiment.VENUS_A_AU
hohmann_dv1 = experiment.hohmann_dv1
hohmann_dv2 = experiment.hohmann_dv2
hohmann_dv_stable = experiment.hohmann_dv_stable
hohmann_dv_total = experiment.hohmann_dv_total
hohmann_split = experiment.hohmann_split
hohmann_transfer_time = experiment.hohmann_transfer_time
hyperbolic_excess_requires = experiment.hyperbolic_excess_requires
kov = experiment.kov
optimality_scan = experiment.optimality_scan
peak_of_cost_curve = experiment.peak_of_cost_curve
ratio_sweep = experiment.ratio_sweep
real_cases = experiment.real_cases
transfer_elements = experiment.transfer_elements
tsiolkovsky_fraction = experiment.tsiolkovsky_fraction
two_impulse_cost = experiment.two_impulse_cost
validate_transfer_rk4 = experiment.validate_transfer_rk4


# --- Closed-form machinery from first principles ---------------------------


def test_dv1_equals_vis_viva_periapsis_speed_minus_circular():
    for r1, r2 in ((1.0, 2.0), (1.0, 6.41), (1.0, 20.0), (0.7, 1.0)):
        v1 = np.sqrt(MU / r1)
        a_t = 0.5 * (r1 + r2)
        v_p = np.sqrt(MU * (2.0 / r1 - 1.0 / a_t))  # vis-viva
        assert hohmann_dv1(r1, r2) == pytest.approx(v_p - v1, rel=1e-12)


def test_dv2_equals_circular_minus_vis_viva_apoapsis_speed():
    for r1, r2 in ((1.0, 2.0), (1.0, 6.41), (1.0, 20.0), (0.7, 1.0)):
        v2 = np.sqrt(MU / r2)
        a_t = 0.5 * (r1 + r2)
        v_a = np.sqrt(MU * (2.0 / r2 - 1.0 / a_t))
        assert hohmann_dv2(r1, r2) == pytest.approx(v2 - v_a, rel=1e-12)


def test_transfer_time_is_half_the_transfer_ellipse_period():
    for r1, r2 in ((1.0, 2.0), (1.0, 6.41), (1.0, 20.0)):
        a_t = 0.5 * (r1 + r2)
        period = 2.0 * np.pi * np.sqrt(a_t**3 / MU)
        assert hohmann_transfer_time(r1, r2) == pytest.approx(
            0.5 * period, rel=1e-12)


def test_angular_momentum_constant_across_the_transfer():
    """r1*v(r1) must equal r2*v(r2) (apsidal h conservation) regardless of
    direction."""
    for r1, r2 in ((1.0, 2.0), (1.0, 6.41), (1.0, 20.0), (0.7, 1.0)):
        els = transfer_elements(r1, r2)
        assert r1 * els["v_at_r1"] == pytest.approx(r2 * els["v_at_r2"], rel=1e-12)
        assert els["h"] == pytest.approx(r1 * els["v_at_r1"], rel=1e-12)
        assert els["h"] == pytest.approx(r2 * els["v_at_r2"], rel=1e-12)


def test_transfer_energy_and_geometry():
    for r1, r2 in ((1.0, 2.0), (1.0, 6.41), (1.0, 20.0)):
        els = transfer_elements(r1, r2)
        assert els["a"] == pytest.approx(0.5 * (r1 + r2), rel=1e-12)
        assert els["e"] == pytest.approx((r2 - r1) / (r2 + r1), rel=1e-12)
        assert els["energy"] == pytest.approx(-MU / (2.0 * els["a"]), rel=1e-12)
        assert els["p"] == pytest.approx(
            els["a"] * (1.0 - els["e"] ** 2), rel=1e-12)


def test_transfer_elements_match_002_elements_machinery():
    """The transfer ellipse built from (r1, r2) equals the elements recovered
    from the propagated departure state by Experiment 002's machinery."""
    for r1, r2 in ((1.0, 2.0), (1.0, 6.41)):
        els = transfer_elements(r1, r2)
        r0 = np.array([r1, 0.0])
        v0 = np.array([0.0, els["v_at_r1"]])
        rec = kov.orbital_elements(r0, v0, MU)
        assert rec["semi_major_axis"] == pytest.approx(els["a"], rel=1e-12)
        assert rec["eccentricity"] == pytest.approx(els["e"], rel=1e-12)
        assert rec["angular_momentum"] == pytest.approx(els["h"], rel=1e-12)
        assert rec["energy"] == pytest.approx(els["energy"], rel=1e-12)


def test_inward_and_outward_burns_swap_magnitudes():
    """The same transfer ellipse in both directions: the burn magnitudes at
    each radius are identical (direction only swaps the order of dv1/dv2)."""
    for R in (1.01, 1.5, 6.409676, 20.0):
        out = hohmann_split(1.0, R)
        inn = hohmann_split(R, 1.0)
        assert inn["dv1"] == pytest.approx(out["dv2"], rel=1e-12), f"R={R}"
        assert inn["dv2"] == pytest.approx(out["dv1"], rel=1e-12), f"R={R}"


def test_degenerate_ratio_r2_equals_r1_is_zero_cost():
    split = hohmann_split(1.0, 1.0)
    assert split["dv_total"] == 0.0
    assert split["transfer"]["e"] == 0.0
    assert split["transfer"]["a"] == pytest.approx(1.0, rel=1e-12)


def test_first_burn_never_exceeds_the_escape_burn():
    """dv1 < (sqrt(2)-1)*v1 = the escape burn from r1, at any radius ratio."""
    v_esc = np.sqrt(2.0) - 1.0
    for R in (1.01, 1.5, 6.41, 20.0, 100.0, 1e4):
        dv1 = hohmann_dv1(1.0, R)
        assert dv1 < v_esc, f"R={R}: dv1={dv1} >= escape burn {v_esc}"


def test_stable_forms_agree_with_textbook_forms():
    for R in (1.0001, 1.001, 1.01, 2.0, 6.41, 20.0, 1000.0):
        d1, d2 = hohmann_dv_stable(1.0, R)
        assert d1 == pytest.approx(hohmann_dv1(1.0, R), rel=1e-9), f"R={R}"
        assert d2 == pytest.approx(hohmann_dv2(1.0, R), rel=1e-9), f"R={R}"


# --- R-sweep: asymptotes, interior maximum, symmetry ------------------------


def test_small_R_asymptote_half_R_minus_one():
    """dv/v1 ~ (R-1)/2 as R -> 1 (leading order, from the stable forms)."""
    for R in (1.0001, 1.001):
        dv = hohmann_dv1(1.0, R) + hohmann_dv2(1.0, R)
        ratio = dv / (0.5 * (R - 1.0))
        assert abs(ratio - 1.0) < 1e-3, f"R={R}: ratio {ratio}"


def test_large_R_asymptote_escape_burn():
    """dv/v1 -> sqrt(2)-1 as R -> infinity (the parabola limit: dv1 tends to
    the escape burn and dv2 to zero)."""
    sweep = ratio_sweep()
    assert sweep["large_R_rel_gap"] < 1e-4
    assert abs(sweep["asymptote_large_R"] - (np.sqrt(2.0) - 1.0)) < 1e-12


def test_peak_location_and_value():
    pk = peak_of_cost_curve()
    assert 15.5 < pk["R_star"] < 15.65, pk["R_star"]
    assert 0.535 < pk["dv_over_v1_peak"] < 0.538
    # Consistency with the dense-grid argmax.
    assert abs(pk["R_star"] - pk["grid_argmax_R"]) / pk["grid_argmax_R"] < 1e-3


def test_peak_sits_above_both_asymptotes():
    """Total dv peaks above sqrt(2)-1 and far above (R-1)/2: fuel cost is
    non-monotone in the radius ratio."""
    pk = peak_of_cost_curve()
    assert pk["dv_over_v1_peak"] > (np.sqrt(2.0) - 1.0)
    assert pk["dv_over_v1_peak"] > 0.53


def test_inward_outward_symmetry():
    sweep = ratio_sweep()
    assert sweep["outward_inward_symmetry_max_rel_diff"] < 1e-12


def test_dv1_increases_and_dv2_vanishes_at_both_limits():
    """dv1 rises monotonically with R; dv2 starts at 0 (R=1), peaks inside
    (max 0.190 at R ~ 5.88), and decays to 0 (R -> inf): the circularization
    burn is largest for intermediate-ratio targets."""
    cells = ratio_sweep()["cells"]
    dv1s = [c["dv1_over_v1"] for c in cells]
    dv2s = [c["dv2_over_v1"] for c in cells]
    assert all(b > a for a, b in zip(dv1s, dv1s[1:]))
    assert all(d > 0.0 for d in dv2s)
    assert dv2s[0] < 1e-5  # R = 1.000001
    assert dv2s[-1] < 1e-3  # R = 1e12
    assert max(dv2s) < 0.191  # interior maximum, measured peak ~0.1900 at R~5.88


def test_transfer_ellipse_eccentricity_tends_to_1():
    cells = ratio_sweep()["cells"]
    assert cells[-1]["e_transfer"] > 0.999999  # R = 1e12
    assert cells[0]["e_transfer"] < 1e-5  # R = 1.000001


# --- RK4 trajectory validation of the transfer ------------------------------


def test_rk4_arrival_matches_closed_form():
    for case in (dict(r1=1.0, r2=1.5), dict(r1=1.0, r2=6.409676),
                 dict(r1=1.0, r2=20.0)):
        r = validate_transfer_rk4(case["r1"], case["r2"])
        tol = 1e-4 if case["r2"] == 20.0 else 1e-6
        assert r["arrival_rk4"]["rel_r_error"] < tol, (
            f"R={case['r2']}: {r['arrival_rk4']['rel_r_error']}")
        assert r["arrival_rk4"]["rel_v_error"] < tol
        assert r["arrival_rk4"]["v_radial_over_v"] < tol


def test_inward_transfer_also_lands_on_target():
    """True inward transfers (r2 < r1) must land on r2 with speed v_r2."""
    for r1, r2 in ((1.0 / 1.3825, 1.0), (2.0, 1.0), (1.0, 0.5)):
        r = validate_transfer_rk4(r1, r2)
        assert r["arrival_rk4"]["rel_r_error"] < 1e-6
        assert r["arrival_rk4"]["rel_v_error"] < 1e-6
        assert r["arrival_rk4"]["v_radial_over_v"] < 1e-6


def test_analytic_arrival_is_exact():
    """The closed-form transfer orbit (kepler_solution of the ellipse) must
    hit r2 with speed v_a at t_tr to machine precision - this validates the
    transfer-ellipse algebra itself, independent of the integrator. Covered
    for both directions (the inward flight starts at apoapsis, i.e. the
    same ellipse a half-period later; the reference is phase-shifted)."""
    for r1, r2 in ((1.0, 1.5), (1.0, 6.409676), (1.0, 20.0),
                   (2.0, 1.0), (1.0, 0.5)):
        r = validate_transfer_rk4(r1, r2)
        assert r["arrival_analytic"]["rel_r_error"] < 1e-9
        assert r["arrival_analytic"]["rel_v_error"] < 1e-9


def test_apex_is_reached_exactly_at_transfer_time():
    """The final propagated point is the apside (apoapsis outward, periapsis
    inward), i.e. the extremal radius occurs at t = t_tr."""
    assert validate_transfer_rk4(1.0, 6.409676)["arrival_rk4"]["apsis_at_final"] is True
    assert validate_transfer_rk4(1.0 / 1.3825, 1.0)["arrival_rk4"]["apsis_at_final"] is True


def test_post_burn_orbit_is_circular_at_r2():
    r = validate_transfer_rk4(1.0, 6.409676)
    circ = r["post_burn_circular_orbit"]
    assert circ["radius_max_rel_variation"] < 1e-6
    assert circ["speed_rel_error_vs_sqrt_mu_r2"] < 1e-6


def test_invariants_drift_on_the_transfer_is_small():
    r = validate_transfer_rk4(1.0, 6.409676)
    assert r["max_rel_drift"]["energy"] < 1e-6
    assert r["max_rel_drift"]["angular_momentum"] < 1e-6


def test_burn2_velocity_difference_equals_closed_form():
    """The measured vector difference v2_circ - v_final must equal dv2."""
    for r2 in (1.5, 6.409676, 20.0):
        r = validate_transfer_rk4(1.0, r2)
        assert r["burn2"]["rel_dv2_error"] < 1e-4, f"R={r2}"


# --- Two-impulse optimality scan --------------------------------------------


def test_two_impulse_cost_reproduces_hohmann_at_corner():
    for r1, r2 in ((1.0, 2.0), (1.0, 6.409676), (1.0, 20.0)):
        corner = two_impulse_cost(r1, r2, MU, r1, r2)
        assert corner == pytest.approx(hohmann_dv_total(r1, r2), rel=1e-12)


def test_grid_minimum_is_the_hohmann_transfer():
    for R in (2.0, 6.409676, 20.0):
        o = optimality_scan(1.0, R)
        g = o["grid"]
        assert g["rel_gap_min_vs_hohmann"] < 1e-5, f"R={R}"
        assert g["argmin_is_hohmann_corner"] is True, f"R={R}"
        assert abs(g["argmin_rp_over_r1"] - 1.0) < 1e-6
        assert abs(g["argmin_ra_over_r2"] - 1.0) < 1e-6


def test_tangent_departure_family_minimum_at_r_a_equals_r2():
    o = optimality_scan(1.0, 6.409676)
    c = o["curve_tangent_departure"]
    assert c["min_rel_excess_over_hohmann"] < 1e-9
    assert abs(c["min_at_r_a_over_r2"] - 1.0) < 1e-6
    assert c["monotone_increasing_after_min"] is True


def test_fixed_apoapsis_family_minimum_at_r_p_equals_r1():
    o = optimality_scan(1.0, 6.409676)
    c = o["curve_apoapsis_fixed"]
    assert c["min_rel_excess_over_hohmann"] < 1e-9
    assert abs(c["min_at_r_p_over_r1"] - 1.0) < 1e-6
    assert c["monotone_increasing_toward_hohmann"] is True


# --- Real-system anchors ----------------------------------------------------


def test_leo_geo_budget_matches_canonical():
    """The canonical LEO(200 km) -> GEO budget is ~3.93 km/s total, ~5.3 h
    coast (Curtis Ch. 6; widely quoted as ~3.9-3.95 km/s)."""
    lg = real_cases()["leo_geo"]
    assert 3.90 < lg["dv_total"] < 3.99
    assert 2.40 < lg["dv1"] < 2.50
    assert 1.40 < lg["dv2"] < 1.55
    assert 5.0 < lg["transfer_time_hours"] < 5.6
    assert 0.0 < lg["fuel_frac_isp300"] < 1.0
    assert lg["fuel_frac_isp300"] > lg["fuel_frac_isp450"] > 0.5


def test_earth_mars_anchors_match_published_values():
    """Published canonical values: transfer ~259 days (NASA/marspedia), v_inf
    at departure ~2.94 km/s and arrival ~2.65 km/s, trans-Mars injection from
    LEO ~3.6 km/s (Wikipedia; public calculators)."""
    em = real_cases()["earth_mars"]
    assert 258.0 < em["transfer_time_days"] < 261.0
    assert 2.90 < em["vinf_departure_km_s"] < 3.00
    assert 2.60 < em["vinf_arrival_km_s"] < 2.72
    assert 3.50 < em["tmi_from_200km_leo_km_s"] < 3.75
    assert 5.5 < em["total_leo_to_mars_orbit_km_s"] < 6.0


def test_earth_venus_inward_transfer():
    """Inward transfer Earth -> Venus: ~146 d, departure excess ~2.5 km/s,
    arrival excess ~2.7 km/s (computed from IAU/JPL constants; the return
    transfer proportions follow from the burn-swap identity)."""
    ev = real_cases()["earth_venus"]
    assert 140.0 < ev["transfer_time_days"] < 152.0
    assert 2.35 < ev["vinf_departure_km_s"] < 2.65
    assert 2.55 < ev["vinf_arrival_km_s"] < 2.85
    assert 3.35 < ev["tmi_from_200km_leo_km_s"] < 3.65


def test_hyperbolic_excess_consistency():
    v = 7.784  # km/s, LEO-like
    assert hyperbolic_excess_requires(0.0, v) == pytest.approx(
        (np.sqrt(2.0) - 1.0) * v, rel=1e-12)
    # Reaching any v_infinity costs more than escaping alone.
    assert hyperbolic_excess_requires(2.94, v) > hyperbolic_excess_requires(0.0, v)


def test_tsiolkovsky_fraction_monotone():
    assert tsiolkovsky_fraction(0.0, 300.0) == pytest.approx(0.0, abs=1e-12)
    f1 = tsiolkovsky_fraction(1.0, 300.0)
    f2 = tsiolkovsky_fraction(4.0, 300.0)
    assert 0.0 < f1 < f2 < 1.0
    # Sanity against a hand-computed value: 3930 m/s at Isp 300 s
    # (exhaust 2941.995 m/s) -> 1 - exp(-1.33583) = 0.73721.
    assert tsiolkovsky_fraction(3.93, 300.0) == pytest.approx(0.7372, rel=1e-3)


def test_earth_mars_budget_matches_geometry():
    """The heliocentric budget must equal the canonical-units curve value at
    the same ratio: dv/v1(R=1.523679) ~ 0.1792 (dv1) + small dv2."""
    em = real_cases()["earth_mars"]
    R = MARS_A_AU
    split = hohmann_split(AU_KM, R * AU_KM, MU_SUN_KM3S2)
    assert em["dv_total"] == pytest.approx(split["dv_total"], rel=1e-12)


# --- Constants sanity -------------------------------------------------------


def test_constants_sanity():
    assert abs(MU_EARTH_KM3S2 - 3.986004e5) < 1.0
    assert abs(R_EARTH_KM - 6.3781e3) < 1.0
    assert abs(GEO_ALT_KM - 3.5786e4) < 1.0
    assert abs(AU_KM - 1.495978707e8) < 1.0
    assert 1.52 < MARS_A_AU < 1.53
    assert 0.72 < VENUS_A_AU < 0.73


# --- Determinism ------------------------------------------------------------

# (Deliberately commented-out pattern documentation: the sweep and the
#  optimality scan are pure functions; determinism is enforced by the absence
#  of RNG. The subprocess check below follows the 002/003 pattern.)


def test_determinism_across_processes():
    """A fresh interpreter must produce bit-identical numerical output."""
    import json
    import subprocess
    import sys

    exp_dir = Path(__file__).resolve().parents[1]
    script = (
        "import importlib.util, json, sys\n"
        f"spec = importlib.util.spec_from_file_location('exp4', "
        f"{str(exp_dir / 'experiment.py')!r})\n"
        "m = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(m)\n"
        "payload = {\n"
        "    'dv': m.hohmann_dv_total(1.0, 6.409676),\n"
        "    't': m.hohmann_transfer_time(1.0, 6.409676),\n"
        "    'rk4': m.validate_transfer_rk4(1.0, 6.409676)['arrival_rk4'],\n"
        "    'opt': m.optimality_scan(1.0, 6.409676)['grid']['cost_min'],\n"
        "}\n"
        "print(json.dumps(payload, sort_keys=True, default=float))\n"
    )
    out = subprocess.check_output([sys.executable, "-c", script], text=True).strip()
    here = {
        "dv": hohmann_dv_total(1.0, 6.409676),
        "t": hohmann_transfer_time(1.0, 6.409676),
        "rk4": validate_transfer_rk4(1.0, 6.409676)["arrival_rk4"],
        "opt": optimality_scan(1.0, 6.409676)["grid"]["cost_min"],
    }
    expected = json.dumps(here, sort_keys=True, default=float)
    assert out == expected, "results differ between interpreter processes"