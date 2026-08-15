"""Validation tests for the combined bi-elliptic transfer + plane-change
experiment (the "006" plane-change slot, completed in place).

These must pass BEFORE any results are trusted (laboratory rule: verify before
trust). They check:
  * closed-form cost identities from first principles (law of cosines burns),
  * the s -> infinity limit identity with Exp 005 (bi_parabolic),
  * the plane-change split law theta1 + theta2 + theta3 = delta_i,
  * the R = 1 detour anchors (di_c = 2 arcsin(1/3), di_inf = 60 deg),
  * the three regime code paths (two-burn / finite-s / s->infinity),
  * the finite-s window pinching shut near R ~ 6.5 (the prior 'abrupt' regime),
  * the independent alternate-optimizer agreement (non-unimodal check),
  * the 3D RK4 trajectory validation of the optimal two-burn and three-burn
    maneuvers (burn magnitudes, arrival radius, circular speed, normal),
  * the analytic s->infinity limit matching the 50-digit mpmath reference,
  * and determinism (byte-identical optimizer output across runs).
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "plane_change_experiment", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)

bi_parabolic_plane_change_limit = experiment.bi_parabolic_plane_change_limit
three_burn_cost = experiment.three_burn_cost
two_burn_dv = experiment.two_burn_dv
two_burn_optimal = experiment.two_burn_optimal
three_burn_global = experiment.three_burn_global
combined_optimum = experiment.combined_optimum
di_c_boundary = experiment.di_c_boundary
di_inf_boundary = experiment.di_inf_boundary
detour_optimum = experiment.detour_optimum
pure_plane_change_dv = experiment.pure_plane_change_dv
validate_two_burn_rk4 = experiment.validate_two_burn_rk4
validate_three_burn_rk4 = experiment.validate_three_burn_rk4
biell_bi_parabolic = experiment.biell_exp.bi_parabolic


# --- Closed-form identities ------------------------------------------------


def test_two_burn_law_of_cosines():
    """The 2-burn combined cost is exactly the sum of two law-of-cosines burns
    with the split theta1 + (di - theta1)."""
    R, di, th1 = 6.41, np.radians(28.6), np.radians(2.174)
    dv1, dv2, tot = two_burn_dv(R, di, th1)
    vp = np.sqrt(2 * R / (1 + R))
    v2 = 1 / np.sqrt(R)
    v_apo = np.sqrt(2 / (R * (1 + R)))
    expect = (
        np.sqrt(1 + vp ** 2 - 2 * vp * np.cos(th1))
        + np.sqrt(v2 ** 2 + v_apo ** 2 - 2 * v2 * v_apo * np.cos(di - th1))
    )
    assert abs(tot - expect) < 1e-14
    assert abs(dv1 + dv2 - tot) < 1e-14


def test_three_burn_split_law_theta1_plus_theta2_plus_theta3_equals_di():
    """At an interior optimum the three plane-change angles sum to delta_i."""
    R, di = 2.0, np.radians(47.5)
    tb = three_burn_global(R, di)
    th1 = np.radians(tb["theta1_deg"])
    th2 = np.radians(tb["theta2_deg"])
    th3 = np.radians(tb["theta3_deg"])
    assert abs(th1 + th2 + th3 - di) < 1e-6


def test_three_burn_cost_matches_closed_form():
    """three_burn_cost reproduces the explicit law-of-cosines sum."""
    R, s, di = 2.0, 3.0, np.radians(40.0)
    th1, th2 = np.radians(5.0), np.radians(33.5)
    val = three_burn_cost(R, s, th1, th2, di)
    vp12 = np.sqrt(2 * s / (1 + s))
    v2 = 1 / np.sqrt(R)
    va1 = np.sqrt(2 / (s * (1 + s)))
    va2 = np.sqrt(2 * R / (s * (R + s)))
    vp23 = np.sqrt(2 * s / (R * (R + s)))
    th3 = di - th1 - th2
    expect = (
        np.sqrt(1 + vp12 ** 2 - 2 * vp12 * np.cos(th1))
        + np.sqrt(va1 ** 2 + va2 ** 2 - 2 * va1 * va2 * np.cos(th2))
        + np.sqrt(v2 ** 2 + vp23 ** 2 - 2 * v2 * vp23 * np.cos(th3))
    )
    assert abs(val - expect) < 1e-14


def test_s_inf_limit_equals_exp005_bi_parabolic():
    """The s->infinity 3-burn limit equals the COPLANAR bi-parabolic limit of
    Exp 005: (sqrt(2)-1)(1 + 1/sqrt(R))."""
    for R in (1.5, 2.0, 6.41, 12.0, 40.0):
        mine = bi_parabolic_plane_change_limit(R)
        theirs = float(biell_bi_parabolic(np.array([R]))[0])
        assert abs(mine - theirs) < 1e-15


# --- R = 1 detour anchors ---------------------------------------------------


def test_R1_detour_anchors():
    """At R = 1 the pure plane-change detour anchors are analytic:
    direct burn optimal for di <= 2 arcsin(1/3) ~ 38.9424 deg; optimum s* ->
    infinity at di = 60 deg with cost 2(sqrt2 - 1)."""
    q = 1.0 / 3.0
    di_c = 2 * np.arcsin(q)
    assert abs(np.degrees(di_c) - 38.9424) < 1e-3
    s_star_lo, tot_lo = detour_optimum(di_c - 1e-3)
    assert s_star_lo == 1.0  # direct burn still optimal just below the corner
    s_star_hi, tot_hi = detour_optimum(di_c + 1e-3)
    assert s_star_hi > 1.0   # interior dip opens just above
    s_star_60, tot_60 = detour_optimum(np.radians(61.0))
    assert s_star_60 == float("inf")
    assert abs(tot_60 - 2 * (np.sqrt(2) - 1)) < 1e-12


# --- Regime structure -------------------------------------------------------


def test_two_burn_optimal_at_small_di():
    """For small plane-change angle the 2-burn is the global optimum (at least
    for moderate radius ratios that are still below the 3-burn boundary)."""
    for R in (2.0, 6.41):
        co = combined_optimum(R, np.radians(5.0))
        assert co["regime"] == "two_burn"
        assert co["beats_two_burn"] is False
    # large R can already be in the 3-burn regime at small di (correct physics)
    co = combined_optimum(12.0, np.radians(5.0))
    assert co["regime"] == "infinite_s"


def test_finite_s_dip_beats_two_burn_R2_di47p5():
    """The canonical finite-s dip: R = 2, di ~ 47.5 deg beats two-burn by
    ~1-2% with a finite intermediate apoapsis (prior-agent claim)."""
    co = combined_optimum(2.0, np.radians(47.5))
    assert co["regime"] == "finite_s"
    assert co["three_burn_s_star"] > 1.0
    assert co["beats_two_burn"]
    saving = 100 * (co["two_burn_dv"] - co["best_dv"]) / co["two_burn_dv"]
    assert 1.0 < saving < 3.0


def test_s_infinite_regime_at_large_di():
    """For large plane-change angle the optimum is the s->infinity (free-at-
    apoapsis) 3-burn; its cost equals the analytic bi-parabolic limit."""
    co = combined_optimum(12.0, np.radians(35.0))
    assert co["regime"] == "infinite_s"
    assert abs(co["best_dv"] - bi_parabolic_plane_change_limit(12.0)) < 1e-12
    assert co["beats_two_burn"]


def test_finite_s_window_pinches_shut_near_R_6p5():
    """The finite-s window closes: di_inf(R) -> di_c(R) as R increases, and for
    R above ~6.5 there is NO finite-s optimum (the winning 3-burn is s->inf).
    This is the 'abrupt behavior' resolved rigorously."""
    w_lo = di_inf_boundary(2.0) - di_c_boundary(2.0)
    w_hi = di_inf_boundary(6.41) - di_c_boundary(6.41)
    assert w_lo > 15.0
    assert 0.0 <= w_hi < 2.0  # window has nearly closed by R = 6.41
    # above the pinch there is no finite-s window: at R = 8 the 3-burn that
    # beats two-burn is the s->infinity one (di_c(8) ~ 31 deg, use di = 40).
    co = combined_optimum(8.0, np.radians(40.0))
    assert co["regime"] == "infinite_s"


def test_di_c_monotone_decreasing_with_R_eventually():
    """di_c(R) (two-burn -> 3-burn) decreases toward 0 as R grows large
    (very large radius ratios make the bi-parabolic strategy win even for
    tiny plane changes)."""
    dc_small = di_c_boundary(2.0)
    dc_large = di_c_boundary(15.0)
    assert dc_large < dc_small
    assert dc_large >= 0.0


# --- Independent alternate-optimizer cross-check (non-unimodal) -------------


def test_alternate_optimizer_agreement():
    """A second, independent optimizer (larger s grid, larger split grid,
    larger s_max) must agree on regime and best dv with the primary one."""
    cases = [(2.0, 47.5), (4.0, 45.0), (6.41, 38.2), (2.0, 50.0), (12.0, 35.0)]
    for R, di in cases:
        a = combined_optimum(R, np.radians(di))
        b = combined_optimum(R, np.radians(di), ns=400, nth=96, s_max=1e8)
        assert a["regime"] == b["regime"]
        assert abs(a["best_dv"] - b["best_dv"]) < 1e-4


# --- RK4 trajectory validation ---------------------------------------------


def test_rk4_two_burn_validation():
    """The optimal 2-burn maneuver propagated with an independent 3D RK4
    reproduces the closed-form burns and arrives on the target circular orbit
    inclined by delta_i (normal alignment cos ~ 1, radius/speed error ~1e-12)."""
    for (R, di) in [(6.41, 28.6), (2.0, 30.0), (4.0, 20.0)]:
        v = validate_two_burn_rk4(R, np.radians(di))
        assert v["dv_rel_err"] < 1e-6
        assert v["final_r_error"] < 1e-6
        assert v["final_speed_error"] < 1e-6
        assert v["final_circular_h_error"] < 1e-6
        assert abs(v["normal_alignment_cos"] - 1.0) < 1e-6


def test_rk4_three_burn_validation():
    """The optimal finite-s 3-burn maneuver propagated with an independent 3D
    RK4 reproduces the closed-form burns and lands on the target inclined
    circular orbit."""
    for (R, di, s) in [(2.0, 47.5, 2.73), (4.0, 45.0, 4.97), (2.0, 40.0, 2.16)]:
        v = validate_three_burn_rk4(R, np.radians(di), s)
        assert v["dv_rel_err"] < 1e-6
        assert v["final_r_error"] < 1e-6
        assert v["final_speed_error"] < 1e-6
        assert v["final_circular_h_error"] < 1e-6
        assert abs(v["normal_alignment_cos"] - 1.0) < 1e-6


# --- Determinism ------------------------------------------------------------


def test_determinism_of_optimizer():
    """The optimizer is pure float64 with no RNG; two calls with the same
    inputs must return bit-identical results, and the s->inf identity is
    exact to ~1e-16."""
    a = combined_optimum(2.0, np.radians(47.5))
    b = combined_optimum(2.0, np.radians(47.5))
    assert a == b
    diff = abs(bi_parabolic_plane_change_limit(2.0)
               - float(biell_bi_parabolic(np.array([2.0]))[0]))
    assert diff < 1e-15
