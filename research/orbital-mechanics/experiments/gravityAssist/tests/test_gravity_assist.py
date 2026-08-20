"""Validation tests for Experiment 007: planetary gravity assist (patched conic).

Validation ladder:
  L1  analytic identities (turn-angle forms, impact parameter, ceiling)
  L2  independent vector construction (B-plane rotation, magnitude preservation)
  L3  independent 3D Cowell propagation with element recovery
  L4  conservation (planet-frame energy/speed)
  L5  canonical published anchors (Voyager 1/2 Jupiter, Voyager 1 Saturn)
  L6  pathological regimes (near-parabolic, ultra-weak, monotonicity)
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("gravity_assist_exp", _DIR / "experiment.py")
assert _spec is not None and _spec.loader is not None
ga = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ga)


# --------------------------------------------------------------------- #
# L1 — analytic identities
# --------------------------------------------------------------------- #
@pytest.mark.parametrize("x", [1e-6, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e6, 1e10])
def test_turn_angle_forms_agree(x):
    stable = ga.turn_angle_rad(x)
    naive = ga.turn_angle_naive_rad(x)
    assert abs(stable - naive) < 1e-12
    # arcsin(1/e) == atan2(1, sqrt(x(x+2))) identity
    assert abs(math.sin(stable / 2) - 1.0 / (1.0 + x)) < 1e-12


def test_turn_angle_near_pi_form():
    for x in (1e-12, 1e-9, 1e-6, 1e-3):
        assert abs((math.pi - ga.turn_angle_rad(x)) - ga.turn_angle_near_pi_rad(x)) < 1e-14


@pytest.mark.parametrize("x", [0.01, 1.0, 100.0])
def test_impact_parameter_consistency(x):
    """b = r_p sqrt(1+2/x) must equal (mu/v_inf^2) sqrt(e^2-1)."""
    mu, v_inf = 1000.0, 2.0
    r_p = x * mu / v_inf**2
    b1 = ga.impact_parameter_km(r_p, x)
    e = 1.0 + x
    b2 = mu / v_inf**2 * math.sqrt(e**2 - 1.0)
    assert abs(b1 - b2) / b2 < 1e-12


def test_periapsis_speed_vis_viva():
    mu, r_p, v_inf = 398600.435507, 7000.0, 3.0
    v_p = ga.periapsis_speed_kms(r_p, v_inf, mu)
    # energy: v_inf^2/2 == v_p^2/2 - mu/r_p
    assert abs(0.5 * v_inf**2 - (0.5 * v_p**2 - mu / r_p)) < 1e-9


def test_soi_values_match_published_scale():
    soi = {k: ga.soi_radius_km(v) for k, v in ga.PLANETS.items()}
    assert abs(soi["Earth"] - 9.2465e5) / 9.2465e5 < 1e-3
    assert abs(soi["Jupiter"] - 4.8206e7) / 4.8206e7 < 1e-3
    # Saturn's a_p varies ~1% across JPL products; SOI is a convention anyway
    assert abs(soi["Saturn"] - 5.4545e7) / 5.4545e7 < 1e-2
    assert soi["Mars"] < soi["Venus"] < soi["Earth"]


# --------------------------------------------------------------------- #
# L2 — independent vector construction
# --------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "delta_deg,beta", [(10.0, 0.0), (98.6051969, 0.7), (170.0, 2.1), (45.0, 5.0)]
)
def test_flyby_preserves_v_inf_magnitude(delta_deg, beta):
    v_in = np.array([3.0, 0.4, -1.0])
    v_out = ga.flyby_v_inf_out(v_in, math.radians(delta_deg), beta, b_km=1e5)
    assert abs(np.linalg.norm(v_out) - np.linalg.norm(v_in)) / np.linalg.norm(v_in) < 1e-12
    # rotation angle between in/out equals delta
    cosang = np.dot(v_in, v_out) / (np.linalg.norm(v_in) * np.linalg.norm(v_out))
    assert abs(math.acos(cosang) - math.radians(delta_deg)) < 1e-10


def test_energy_identity_d_eps_equals_V_dot_dv():
    """d_eps = 1/2(v_out^2 - v_in^2) must equal V_p . dv_helio exactly."""
    p = ga.PLANETS["Jupiter"]
    v_in = np.array([2.0, -1.0, 0.5])
    v_out = ga.flyby_v_inf_out(v_in, math.radians(98.6), 0.3, b_km=3e5)
    vi, vo = ga.assemble_heliocentric(p, v_in, v_out)
    d_eps = 0.5 * (np.dot(vo, vo) - np.dot(vi, vi))
    V = np.array([0.0, 0.0, ga.planet_speed_kms(p)])
    assert abs(d_eps - np.dot(V, vo - vi)) < 1e-9


def test_optimum_bend_parallel_to_V_p_and_value():
    """At alpha* = pi/2 + delta/2 (phi=0): dv parallel +V_p, d_eps = ceiling."""
    p = ga.PLANETS["Jupiter"]
    v_inf, rpf = 10.0, 1.02
    r_p = rpf * p["R_eq"]
    x = ga.x_parameter(r_p, v_inf, p["mu"])
    delta = ga.turn_angle_rad(x)
    h = delta / 2
    alpha_star = math.pi / 2 + h
    k = np.array([0.0, 0.0, 1.0])
    s = np.array([math.sin(alpha_star), 0.0, math.cos(alpha_star)])
    q_des = (k - math.cos(alpha_star) * s) / math.sin(alpha_star)
    beta = ga.beta_for_bend(q_des, s)
    v_out = ga.flyby_v_inf_out(v_inf * s, delta, beta, ga.impact_parameter_km(r_p, x))
    dv = v_out - v_inf * s
    parallel_err = np.linalg.norm(dv - np.dot(dv, k) * k) / np.linalg.norm(dv)
    assert parallel_err < 1e-10
    V = np.array([0.0, 0.0, ga.planet_speed_kms(p)])
    vi, vo = ga.assemble_heliocentric(p, v_inf * s, v_out)
    d_eps = 0.5 * (np.dot(vo, vo) - np.dot(vi, vi))
    ceiling = 2.0 * ga.planet_speed_kms(p) * v_inf * math.sin(h)
    assert abs(d_eps - ceiling) / ceiling < 1e-9


def test_orientation_grid_matches_analytic_max():
    """1-deg grid max of F must be within 0.1% of the analytic ceiling (F=1)."""
    for delta_deg in (30.0, 98.6051969, 170.0):
        h = math.radians(delta_deg / 2)
        f_grid, _, _ = ga.orientation_grid_max(h)
        assert f_grid <= 1.0 + 1e-12
        assert 1.0 - f_grid < 1e-3  # < 0.1% grid loss


def test_landscape_closed_form_against_direct_vectors():
    """d_eps(alpha, phi) closed form == direct vector evaluation."""
    p = ga.PLANETS["Earth"]
    v_inf, delta_deg, alpha_deg, phi_deg = 4.0, 80.0, 130.0, 37.0
    delta = math.radians(delta_deg)
    h = delta / 2
    alpha, phi = math.radians(alpha_deg), math.radians(phi_deg)
    k = np.array([0.0, 0.0, 1.0])
    # unit vector at exactly angle alpha from k (azimuth is irrelevant by symmetry)
    s = np.array([math.sin(alpha), 0.0, math.cos(alpha)])
    e1 = (k - np.dot(k, s) * s) / np.linalg.norm(k - np.dot(k, s) * s)
    e2 = np.cross(k, s) / np.linalg.norm(np.cross(k, s))
    q = math.cos(phi) * e1 + math.sin(phi) * e2
    v_out = v_inf * (math.cos(delta) * s + math.sin(delta) * q)
    vi, vo = ga.assemble_heliocentric(p, v_inf * s, v_out)
    d_eps_direct = 0.5 * (np.dot(vo, vo) - np.dot(vi, vi))
    A = 2.0 * ga.planet_speed_kms(p) * v_inf * math.sin(h)
    d_eps_formula = A * ga.orientation_factor(np.array([alpha]), np.array([phi]), h)[0]
    assert abs(d_eps_direct - d_eps_formula) / max(abs(d_eps_formula), 1e-12) < 1e-10


# --------------------------------------------------------------------- #
# L5 — canonical anchors
# --------------------------------------------------------------------- #
def test_anchor_voyager1_jupiter():
    a = ga.reconstruct_anchor(ga.ANCHORS[0])
    assert abs(a["v_inf_kms"] - 10.7691921) < 2e-4
    assert abs(a["r_p_km"] - 348435.3475) / 348435.3475 < 1e-6
    assert abs(a["delta_deg"] - 98.6051969) < 1e-4
    assert abs(a["dv_vector_kms"] - 16.3296256) < 2e-4
    assert abs(a["d_epsilon_km2s2"] - 200.832) / 200.832 < 1e-4
    assert abs(a["dv_scalar_kms"] - 10.9945) / 10.9945 < 1e-3
    # vector change vs scalar speed change are DIFFERENT observables
    assert a["dv_vector_kms"] > a["dv_scalar_kms"]


def test_anchor_voyager2_jupiter():
    a = ga.reconstruct_anchor(ga.ANCHORS[1])
    assert abs(a["v_inf_kms"] - 7.6159658) < 2e-4
    assert abs(a["r_p_km"] - 721375.5751) / 721375.5751 < 1e-6
    assert abs(a["delta_deg"] - 97.4795289) < 1e-4
    assert abs(a["dv_vector_kms"] - 11.4501783) < 2e-4
    assert abs(a["d_epsilon_km2s2"] - 151.759) / 151.759 < 1e-4


def test_anchor_voyager1_saturn():
    a = ga.reconstruct_anchor(ga.ANCHORS[2])
    assert abs(a["v_inf_kms"] - 15.1093416) < 2e-4
    assert abs(a["r_p_km"] - 184023.4753) / 184023.4753 < 1e-6
    assert abs(a["delta_deg"] - 56.6512718) < 1e-4
    assert abs(a["dv_vector_kms"] - 14.3382247) < 2e-4
    assert abs(a["d_epsilon_km2s2"] - 26.12146) / 26.12146 < 1e-3
    # the physical point of this anchor: big rotation, small energy gain
    assert a["dv_vector_kms"] > 10.0
    assert a["dv_scalar_kms"] < 2.0


# --------------------------------------------------------------------- #
# L3/L4 — independent propagation + conservation
# --------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "planet,rpf,v_inf",
    [("Earth", 1.10, 5.0), ("Jupiter", 1.02, 10.0), ("Jupiter", 1.05, 0.7)],
)
def test_L3_element_recovery_matches_closed_form(planet, rpf, v_inf):
    p = ga.PLANETS[planet]
    r_p = rpf * p["R_eq"]
    delta_an = ga.turn_angle_rad(ga.x_parameter(r_p, v_inf, p["mu"]))
    rec = ga.propagate_flyby_l3(r_p, v_inf, p["mu"])
    assert abs(rec["delta_rad"] - delta_an) / delta_an < 1e-6
    assert abs(rec["r_p_km"] - r_p) / r_p < 1e-6
    assert abs(rec["v_inf_kms"] - v_inf) / v_inf < 1e-8
    assert rec["eps_planet_conservation_rel"] < 1e-8


def test_L3_patch_radius_insensitivity():
    """Recovered delta must not depend on the patch radius (element recovery)."""
    p = ga.PLANETS["Jupiter"]
    deltas = [
        ga.propagate_flyby_l3(1.02 * p["R_eq"], 10.0, p["mu"], R0_factor=f)[
            "delta_rad"
        ]
        for f in (50.0, 100.0, 200.0)
    ]
    assert max(deltas) - min(deltas) < 1e-8


# --------------------------------------------------------------------- #
# L6 — pathological regimes
# --------------------------------------------------------------------- #
def test_pathological_checks_pass():
    res = ga.pathological_checks()
    assert res["all_finite"]
    assert res["delta_monotone_decreasing_in_x"]
    assert res["stable_vs_mpmath_max_rel_err"] < 1e-15
    assert res["near_pi_form_max_abs_err_below_x_1e-2"] < 1e-12


def test_extreme_regimes_no_nan_and_correct_limits():
    # near-parabolic: delta -> pi, |dv| -> 2 v_inf -> 0
    x = 1e-14
    assert abs(ga.turn_angle_rad(x) - math.pi) < 1e-5
    # ultra-weak: delta ~ 2/x (1 - 1/x + ...), relative error ~ 1/x
    for x in (1e6, 1e8, 1e10):
        approx = 2.0 / x
        assert abs(ga.turn_angle_rad(x) - approx) / approx < 2.0 / x
    # ceiling falls as ~1/v_inf at high speed: exact vs leading asymptote
    # differ by the O(1/x) correction, so require |exact/asymptotic - 1| < 1/x
    p = ga.PLANETS["Earth"]
    r_p = 1.02 * p["R_eq"]
    V_p = ga.planet_speed_kms(p)
    for v in (20.0, 30.0):
        x = ga.x_parameter(r_p, v, p["mu"])
        ceil_v = 2.0 * V_p * v / (1.0 + x)
        asym = 2.0 * V_p * p["mu"] / (r_p * v)
        assert abs(ceil_v / asym - 1.0) < 1.0 / x


def test_sweep_uses_physical_radii_not_factors():
    """Regression: the sweep grid must carry r_p in KM (factor x R_eq).

    With the factor mistaken for km, x collapses and the 'ceiling' degenerates
    to 2 V_p v_inf (sin(delta/2) -> 1).  Guard: the sweep value at
    (Jupiter, v_inf = 10 km/s, r_p = 1.02 R_eq) must equal the analytic
    ceiling 2 V_p v_inf/(1 + x) and be strictly below 2 V_p v_inf.
    """
    sweep = ga.run_sweep()
    s = sweep["Jupiter"]
    i_v = int(np.argmin(np.abs(ga.V_INF_KMS - 10.0)))
    v = float(ga.V_INF_KMS[i_v])
    r_p = 1.02 * ga.PLANETS["Jupiter"]["R_eq"]
    x = ga.x_parameter(r_p, v, ga.PLANETS["Jupiter"]["mu"])
    expected = 2.0 * s["V_p_kms"] * v / (1.0 + x)
    got = float(s["d_eps_max_grid_km2s2"][i_v, 0])
    assert got < 0.999 * 2.0 * s["V_p_kms"] * v  # sin factor strictly present
    assert abs(got - expected) / expected < 1e-12
    assert abs(x - 0.0576) < 0.01  # x ~ 0.058, NOT ~1e-6


def test_sweep_best_has_interior_optimum_for_earth():
    """d_eps_max(v) = 2 V_p v/(1 + r_p v^2/mu) peaks at v = sqrt(mu/r_p);
    for Earth at r_min that is ~7.8 km/s, inside the sweep, so the grid best
    must sit near it (not at the v = 30 edge)."""
    sweep = ga.run_sweep()
    s = sweep["Earth"]
    r_p = 1.02 * ga.PLANETS["Earth"]["R_eq"]
    v_star = float(np.sqrt(ga.PLANETS["Earth"]["mu"] / r_p))
    assert 6.0 < v_star < 10.0
    grid_at_rmin = s["d_eps_max_grid_km2s2"][:, 0]
    best_at_v = float(ga.V_INF_KMS[int(np.argmax(grid_at_rmin))])
    assert abs(best_at_v - v_star) / v_star < 0.15  # nearest log-grid point
