"""Validation tests for the bi-elliptic vs Hohmann transfer experiment.

These must pass BEFORE any results are trusted (laboratory rule: verify before
trust). The tests check the closed-form three-burn machinery from first
principles (vis-viva, half-period legs), the s -> R degeneracy onto the
Hohmann transfer, the digit-safe forms near that corner, the two boundary
radius ratios (11.94 bi-parabolic crossover, 15.58 always-cheaper boundary
equal to the Hohmann cost maximum of Experiment 004), the corner identity
d/ds dv_biell|_{s=R} = d/dR dv_H, the classical crossover table (including
the 12 -> 815.82 value and the rejection of the copied 15.81 entry), the
three-regime claims on adversarial grids, the f(s) shape classification, the
fuel-saving peak and asymptote, the never-wins low family and the inward
time-reversal equivalence, the mpmath 50-digit cross-checks, the RK4
validation of the complete three-burn trajectory, real-system anchors
(including the Wikipedia R = 14 example), and determinism.

The experiment module and its dependency chain (Experiment 004, which in turn
loads Experiment 002's verified machinery) are loaded via importlib from
explicit paths (see tools/new_experiment.py).
"""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "bielliptic_vs_hohmann_experiment", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)

MU = experiment.MU
S2_MINUS_1 = experiment.S2_MINUS_1
bi_parabolic = experiment.bi_parabolic
bielliptic_burns = experiment.bielliptic_burns
bielliptic_dv_total = experiment.bielliptic_dv_total
bielliptic_transfer_time = experiment.bielliptic_transfer_time
burn1_stable = experiment.burn1_stable
burn2_stable = experiment.burn2_stable
burn3_stable = experiment.burn3_stable
corner_identity_check = experiment.corner_identity_check
crossover_s = experiment.crossover_s
d_dv_ho_dR = experiment.d_dv_ho_dR
dv_over_v1 = experiment.dv_over_v1
ellipse_speed = experiment.ellipse_speed
f_high = experiment.f_high
f_high_stable = experiment.f_high_stable
f_low = experiment.f_low
high_precision_verification = experiment.high_precision_verification
hohmann_dv_total = experiment.hohmann_dv_total
hohmann_transfer_time = experiment.hohmann_transfer_time
low_family_and_inward = experiment.low_family_and_inward
peak_of_cost_curve = experiment.peak_of_cost_curve
R_always_cheaper = experiment.R_always_cheaper
R_bp_crossover = experiment.R_bp_crossover
real_anchors = experiment.real_anchors
region_verification = experiment.region_verification
saving_and_time = experiment.saving_and_time
shape_diagnostics = experiment.shape_diagnostics
transfer_case_steps = experiment.transfer_case_steps
validate_bielliptic_rk4 = experiment.validate_bielliptic_rk4


# --- Closed-form machinery from first principles ---------------------------


def test_burn1_is_vis_viva_periapsis_speed_minus_circular():
    for r1, r2, r_b in ((1.0, 2.0, 5.0), (1.0, 6.41, 100.0),
                        (1.0, 14.0, 40.0), (1.0, 20.0, 20.5)):
        v1 = np.sqrt(MU / r1)
        a1 = 0.5 * (r1 + r_b)
        dv1, _, _ = bielliptic_burns(r1, r2, r_b)
        assert dv1 == pytest.approx(
            np.sqrt(MU * (2.0 / r1 - 1.0 / a1)) - v1, rel=1e-12)


def test_burn2_is_apoapsis_speed_difference_between_ellipses():
    for r1, r2, r_b in ((1.0, 2.0, 5.0), (1.0, 6.41, 100.0),
                        (1.0, 14.0, 40.0), (1.0, 20.0, 20.5)):
        a1 = 0.5 * (r1 + r_b)
        a2 = 0.5 * (r2 + r_b)
        _, dv2, _ = bielliptic_burns(r1, r2, r_b)
        assert dv2 == pytest.approx(
            abs(ellipse_speed(r_b, a2, MU) - ellipse_speed(r_b, a1, MU)),
            rel=1e-12)


def test_burn3_is_circular_minus_vis_viva_periapsis_speed():
    for r1, r2, r_b in ((1.0, 2.0, 5.0), (1.0, 6.41, 100.0),
                        (1.0, 14.0, 40.0), (1.0, 20.0, 20.5)):
        v2 = np.sqrt(MU / r2)
        a2 = 0.5 * (r2 + r_b)
        _, _, dv3 = bielliptic_burns(r1, r2, r_b)
        assert dv3 == pytest.approx(
            abs(v2 - np.sqrt(MU * (2.0 / r2 - 1.0 / a2))), rel=1e-12)


def test_total_cost_is_sum_of_burns():
    for r1, r2, r_b in ((1.0, 2.0, 5.0), (1.0, 6.41, 100.0),
                        (1.0, 14.0, 40.0), (1.0, 20.0, 20.5)):
        assert bielliptic_dv_total(r1, r2, r_b) == pytest.approx(
            sum(bielliptic_burns(r1, r2, r_b)), rel=1e-12)


def test_normalized_cost_equals_closed_form_f_high():
    """bielliptic_dv_total(1, R, s) / v1 must equal f_high(R, s) exactly."""
    for R, s in ((2.0, 5.0), (6.41, 100.0), (14.0, 40.0), (20.0, 20.5),
                 (12.0, 1000.0)):
        v1 = np.sqrt(MU / 1.0)
        assert bielliptic_dv_total(1.0, R, s) / v1 == pytest.approx(
            float(f_high(np.array([R]), np.array([s]))[0]), rel=1e-12)


def test_transfer_time_is_two_half_periods():
    for r1, r2, r_b in ((1.0, 2.0, 5.0), (1.0, 6.41, 100.0),
                        (1.0, 14.0, 40.0)):
        t = np.pi * (np.sqrt((0.5 * (r1 + r_b)) ** 3 / MU)
                     + np.sqrt((0.5 * (r2 + r_b)) ** 3 / MU))
        assert bielliptic_transfer_time(r1, r2, r_b) == pytest.approx(
            t, rel=1e-12)


def test_degenerate_s_equals_R_reduces_to_hohmann():
    """As r_b -> r2 the bi-elliptic degenerates into the Hohmann transfer:
    dv3 -> 0, dv1 -> the Hohmann departure burn and dv2 -> the Hohmann
    arrival burn."""
    for R in (1.5, 2.0, 6.41, 14.0, 20.0):
        dv1, dv2, dv3 = bielliptic_burns(1.0, R, R)
        assert dv3 == pytest.approx(0.0, abs=1e-12)
        dv1_ho = np.sqrt(2.0 * R / (1.0 + R)) - 1.0
        dv2_ho = (1.0 / np.sqrt(R)) * (1.0 - np.sqrt(2.0 / (1.0 + R)))
        assert dv1 == pytest.approx(dv1_ho, rel=1e-12)
        assert dv2 == pytest.approx(dv2_ho, rel=1e-12)
        assert dv1 + dv2 == pytest.approx(hohmann_dv_total(1.0, R), rel=1e-12)


def test_f_high_at_corner_equals_hohmann_cost():
    for R in (2.0, 6.41, 14.0, 20.0, 100.0):
        assert float(f_high(np.array([R]), np.array([R]))[0]) == pytest.approx(
            float(dv_over_v1(np.array([R]))[0]), rel=1e-9)


def test_bi_parabolic_limit_formula_and_approach():
    """s -> infinity: dv -> (sqrt(2)-1)(1 + 1/sqrt(R)). Leading correction
    is (sqrt(2)/2)(sqrt(R)-3)/s: approached from below for R < 9, from above
    for R > 9 (verified numerically, e.g. R = 14, s = 1e8 gives +5.24e-9)."""
    for R in (2.0, 6.41, 14.0, 20.0, 100.0):
        assert float(bi_parabolic(np.array([R]))[0]) == pytest.approx(
            S2_MINUS_1 * (1.0 + 1.0 / np.sqrt(R)), rel=1e-12)
        fs = float(f_high(np.array([R]), np.array([1e8]))[0])
        bp = float(bi_parabolic(np.array([R]))[0])
        assert abs(fs - bp) < 1e-6 * bp
        corr = (np.sqrt(2.0) / 2.0) * (np.sqrt(R) - 3.0) / 1e8
        assert abs((fs - bp) - corr) < 1e-8 * bp


def test_each_burn_bound_above_escape_or_deep_space():
    """dv1 < escape burn from r1; dv3 < (escape burn)/sqrt(R); dv2 -> 0 as
    s -> infinity and as s -> R."""
    for R in (2.0, 6.41, 14.0, 20.0):
            for s in (R * 1.001, 10.0 * R, 1e6 * R):
                dv1, dv2, dv3 = bielliptic_burns(1.0, R, s)
                assert dv1 < S2_MINUS_1, (R, s, dv1)
                assert dv3 < S2_MINUS_1 / np.sqrt(R), (R, s, dv3)
            assert dv2 < 1e-6  # deep-space burn vanishes at large s
            dv2_corner = bielliptic_burns(1.0, R, R)[1]
            dv2_ho = (1.0 / np.sqrt(R)) * (1.0 - np.sqrt(2.0 / (1.0 + R)))
            assert dv2_corner == pytest.approx(dv2_ho, rel=1e-6)  # Hohmann arrival burn


def test_stable_forms_agree_with_textbook_forms():
    """Where both forms are well conditioned they agree to ~1e-12; near the
    corner s = R the stable forms remain exact (they are the rearranged
    versions, free of difference-of-close-square-roots cancellation), and at
    the corner itself their sum reproduces dv_H to machine precision."""
    for R in (2.0, 6.41, 14.0, 20.0):
        # well away from the corner: agreement to 1e-12
        s = np.array([1.5 * R, 3.0 * R, 10.0 * R])
        assert np.allclose(f_high_stable(R * np.ones_like(s), s),
                           f_high(R * np.ones_like(s), s), rtol=1e-10,
                           atol=1e-14)
        # approaching the corner: both agree with each other and with dv_H
        # (the naive error is ~1e-16 absolute; the stable one is exact)
        s_near = np.array([R * (1.0 + 1e-10), R * (1.0 + 1e-12)])
        stable = f_high_stable(R * np.ones_like(s_near), s_near)
        naive = f_high(R * np.ones_like(s_near), s_near)
        assert np.allclose(stable, naive, rtol=1e-6, atol=1e-10)
        dh = float(dv_over_v1(np.array([R]))[0])
        assert np.all(np.abs(stable - dh) < 1e-8)
        # exact corner values in the stable forms (no cancellation at all)
        assert float(burn1_stable(np.array([1.0]))[0]) == 0.0
        assert float(burn3_stable(np.array([R]), np.array([R]))[0]) == 0.0
        # at s = R the stable sum equals the Hohmann cost to machine precision
        s_corner = np.array([R])
        corner = burn1_stable(s_corner) + burn2_stable(
            np.array([R]), s_corner) + burn3_stable(np.array([R]), s_corner)
        assert float(corner[0]) == pytest.approx(
            float(dv_over_v1(np.array([R]))[0]), rel=1e-12)


# --- The two boundary radius ratios ----------------------------------------


def test_R_bp_crossover_is_11_dot_9388():
    R_bp = R_bp_crossover()
    assert 11.93 < R_bp < 11.95, R_bp
    # defining property: the bi-parabolic limit ties the Hohmann cost
    assert float(bi_parabolic(np.array([R_bp]))[0]) == pytest.approx(
        float(dv_over_v1(np.array([R_bp]))[0]), rel=1e-10)


def test_R_star_is_hohmann_cost_maximum():
    R_star = R_always_cheaper()
    assert 15.5 < R_star < 15.65, R_star
    assert abs(d_dv_ho_dR(np.array([R_star]))[0]) < 1e-10
    # independent cross-check against Experiment 004's own bisection
    pk = peak_of_cost_curve()
    assert abs(R_star - pk["R_star"]) / pk["R_star"] < 1e-4


def test_corner_identity():
    """d/ds f_high at s = R equals d/dR dv_H (finite differences)."""
    res = corner_identity_check(np.array([1.5, 2.0, 6.41, 11.94, 15.0]))
    assert res["max_rel_diff"] < 1e-4
    # at the Hohmann maximum the corner slope is stationary (both ~ 0):
    # absolute check, since the relative one is ill-defined there
    for R in (15.5817,):
        h = 1e-6 * R
        dds = float(
            (
                (f_high_stable(np.array([R]), np.array([R + h]))
                 - f_high_stable(np.array([R]), np.array([R - h])))
                / (2.0 * h)
            )[0]
        )
        ddr = float(d_dv_ho_dR(np.array([R]))[0])
        assert abs(dds - ddr) < 1e-8


# --- Crossover curve and the classical table -------------------------------


def test_crossover_breaks_even_at_s_c():
    for R in (12.0, 13.0, 14.0, 15.0):
        s_c = crossover_s(R)
        g = f_high(R * np.ones(3), np.array([s_c * 0.99, s_c, s_c * 1.01]))
        g -= dv_over_v1(R * np.ones(3))
        assert g[0] > 0.0 and abs(g[1]) / 0.5 < 1e-6 and g[2] < 0.0, (R, s_c)


def test_classical_table_values():
    """Gobetz & Doll / Escobal values, computed here: 12 -> 815.82,
    13 -> 48.90, 14 -> 26.10, 15 -> 18.19, 15.58 -> 15.58."""
    assert 810.0 < crossover_s(12.0) < 830.0
    assert 48.0 < crossover_s(13.0) < 50.0
    assert 25.5 < crossover_s(14.0) < 26.7
    assert 17.9 < crossover_s(15.0) < 18.5
    assert abs(crossover_s(15.58) - 15.58) < 0.05


def test_wikipedia_12_entry_is_transcription_error():
    """The widely copied table entry '12 -> 15.81' cannot be right: the
    crossover for R = 12 is ~815.82 (a dropped digit)."""
    assert crossover_s(12.0) > 800.0


def test_crossover_curve_monotone_in_log_s():
    """s_c decreases smoothly from infinity (R -> R_bp+) to the corner
    (R -> R*-); near R* it still lies within 25% of R."""
    R_bp, R_star = R_bp_crossover(), R_always_cheaper()
    Rs = np.linspace(R_bp + 0.05, R_star - 0.05, 12)
    sc = np.array([crossover_s(float(R)) for R in Rs])
    assert np.all(np.diff(sc) < 0.0)
    assert sc[0] > 100.0  # diverges toward R_bp
    assert abs(sc[-1] - Rs[-1]) / Rs[-1] < 0.25  # approaches the corner


# --- Three-regime claims (adversarial) -------------------------------------


def test_region_verification_passes_and_records_structure():
    R_bp, R_star = R_bp_crossover(), R_always_cheaper()
    rv = region_verification(R_bp, R_star)
    lo, hi = rv["R_lo_region"], rv["R_hi_region"]
    # the function asserts internally on every grid line
    assert lo["n_R"] == 90 and hi["n_R"] == 90
    assert lo["worst_margin"] > 0.0
    assert hi["worst_margin"] > 0.0
    mid = rv["mid_crossing_counts"]
    assert mid["all_exactly_one"] is True
    assert mid["max_crossings_seen"] == 1 and mid["min_crossings_seen"] == 1


def test_hohmann_wins_everywhere_below_R_bp():
    for R in (1.01, 2.0, 6.41, 11.9):
        s = np.logspace(np.log10(R * 1.000001), np.log10(R * 1e6), 300)
        g = f_high_stable(R * np.ones_like(s), s) - dv_over_v1(
            R * np.ones_like(s))
        assert np.all(g > -1e-12), (R, g.min())


def test_bielliptic_wins_everywhere_above_R_star():
    for R in (15.7, 20.0, 50.0, 100.0, 1000.0):
        s = np.logspace(np.log10(R * 1.000001), np.log10(R * 1e6), 300)
        g = f_high_stable(R * np.ones_like(s), s) - dv_over_v1(
            R * np.ones_like(s))
        assert np.all(g < 1e-12), (R, g.max())


def test_unique_crossing_between_boundaries():
    for R in (12.0, 13.0, 14.0, 15.0, 15.4):
        s = np.logspace(np.log10(R * 1.000001), np.log10(R * 1e7), 500)
        g = f_high_stable(R * np.ones_like(s), s) - dv_over_v1(
            R * np.ones_like(s))
        ncross = int(np.sum(np.sign(g[:-1]) * np.sign(g[1:]) < 0))
        assert ncross == 1, (R, ncross)
        assert g[0] > 0.0 and g[-1] < 0.0, R


# --- Shape of f(s) ---------------------------------------------------------


def test_shape_classification_matches_regimes():
    """On the shape grid: monotone increasing below the hump onset (~9.5),
    a single hump between that onset and R*, monotone decreasing above R*.
    The hump onset lies BELOW R_bp: the hump exists for intermediate R yet
    never dips below dv_H until R > R_bp (that is what makes the crossover
    unique on (R_bp, R*))."""
    sd = shape_diagnostics()
    kinds = {s["shape"] for s in sd["shapes"]}
    assert kinds == {"monotone_increasing", "single_hump",
                     "monotone_decreasing"}
    by_R = {s["R"]: s["shape"] for s in sd["shapes"]}
    inc_Rs = [R for R, k in by_R.items() if k == "monotone_increasing"]
    hump_Rs = [R for R, k in by_R.items() if k == "single_hump"]
    dec_Rs = [R for R, k in by_R.items() if k == "monotone_decreasing"]
    R_bp, R_star = R_bp_crossover(), R_always_cheaper()
    assert max(inc_Rs) < R_bp
    assert min(hump_Rs) > max(inc_Rs)
    assert min(dec_Rs) > R_star
    assert max(hump_Rs) < R_star
    onset = sd["hump_onset_R_approx"]
    assert 8.0 < onset < 11.0  # hump appears below R_bp (~11.94)
    assert abs(min(hump_Rs) - onset) < 1.0


def test_hump_cases_decrease_monotonically_after_peak():
    """The monotone decrease after the hump is what makes the crossover
    unique (asserted inside shape_diagnostics; re-verified here)."""
    for R in (12.0, 13.0, 14.0, 15.0):
        s = np.logspace(np.log10(R * 1.000001), np.log10(R * 1e8), 4000)
        f = f_high(R * np.ones_like(s), s)
        i_hump = int(np.argmax(f))
        assert np.all(np.diff(f[i_hump:]) < 0.0), R


# --- Saving curve and time penalty -----------------------------------------


def test_saving_peak_near_R_50():
    st = saving_and_time()
    assert 0.035 < st["peak_saving_over_v1"] < 0.050
    assert 30.0 < st["peak_saving_at_R"] < 100.0


def test_large_R_asymptote_of_the_saving():
    st = saving_and_time()
    assert 0.95 < st["asymptote_ratio"] < 1.05
    assert abs(st["asymptote_pred_at_1e4"] - (2.0 - np.sqrt(2.0)) / 100.0
               ) < 1e-12


def test_time_penalty_grows_with_s():
    """The flight-time penalty is always > 1, grows with s for a fixed R, and
    is extreme for the large-s cases."""
    st = saving_and_time()
    ratios = [c["t_biell_over_t_hohmann"] for c in st["time_ratio_cases"]]
    assert all(r > 1.0 for r in ratios)
    assert ratios[0] < ratios[1]  # 1.9 < 19.5 (R=2/s=3 vs R=6.41/s=30)

    def ratio(R, s):
        t_ho = np.pi * np.sqrt((1.0 + R) ** 3 / 8.0)
        t_be = np.pi * (np.sqrt((1.0 + s) ** 3 / 8.0)
                        + np.sqrt((R + s) ** 3 / 8.0))
        return t_be / t_ho

    # fixed R = 2: ratio strictly increases with s
    ss = np.linspace(3.0, 200.0, 8)
    rr = np.array([ratio(2.0, float(s)) for s in ss])
    assert np.all(np.diff(rr) > 0.0)
    assert ratios[-1] > 10.0


# --- Low family and inward case --------------------------------------------


def test_low_family_never_wins():
    """Intermediates below both orbits: g_low > 0 for every R > 1."""
    for R in (1.1, 2.0, 11.94, 14.0, 20.0):
        lo_s = np.logspace(np.log10(R * 1e-6), np.log10(min(1.0, R) * 0.9999),
                           200)
        g = f_low(R * np.ones_like(lo_s), lo_s) - dv_over_v1(
            R * np.ones_like(lo_s))
        assert np.all(g > -1e-12), (R, g.min())


def test_inward_high_family_wins_above_threshold():
    """Outer/inner = 20 > R*: the inward bi-elliptic beats the inward
    Hohmann for every intermediate apoapsis."""
    li = low_family_and_inward(R_bp_crossover(), R_always_cheaper())
    cases = li["inward_high_family"]["cases"]
    assert len(cases) == 5
    assert all(c["wins"] for c in cases)
    assert li["time_reversal_burn_identity"] is True


def test_time_reversal_burn_identity():
    """Reversing the transfer swaps the first and third burns (the same two
    ellipses, burns applied in the opposite order): fwd[k] == bwd[2-k]."""
    for (a, b, c) in ((1.0, 14.0, 40.0), (1.0, 20.0, 100.0),
                      (14.0, 1.0, 0.025)):
        fwd = bielliptic_burns(a, b, c)
        bwd = bielliptic_burns(b, a, c)
        assert all(abs(fwd[k] - bwd[2 - k]) < 1e-12 for k in range(3))


# --- High-precision mpmath cross-checks ------------------------------------


def test_mpmath_boundaries_match_float64():
    hp = high_precision_verification()
    assert abs(float(hp["R_bp"]["mpmath"]) - hp["R_bp"]["float64"]) < 1e-9
    assert abs(float(hp["R_star"]["mpmath"]) - hp["R_star"]["float64"]) < 1e-9
    assert 11.938 < float(hp["R_bp"]["mpmath"]) < 11.940
    assert 15.581 < float(hp["R_star"]["mpmath"]) < 15.583


def test_mpmath_crossover_table_matches_float64():
    hp = high_precision_verification()["s_c_table"]
    for R, mpv, f64 in zip(hp["R"], hp["mpmath"], hp["float64"]):
        assert abs(float(mpv) - f64) / f64 < 1e-8, R


def test_mpmath_corner_identity_at_50_digits():
    """The identity holds to ~1e-29 at 50 digits (finite-difference balance:
    roundoff ~1e-50/h against truncation O(h^2), h = 1e-20)."""
    hp = high_precision_verification()
    assert float(hp["corner_identity_50_digits"]["rel_diff"]) < 1e-25


# --- RK4 trajectory validation (three burns) -------------------------------


def test_step_law_follows_002_periapsis_resolution():
    # e -> 0: ~ base steps per half orbit; e -> 1: blows up
    assert transfer_case_steps(0.0) == 256
    assert transfer_case_steps(0.9) > transfer_case_steps(0.5)
    assert transfer_case_steps(0.99) > transfer_case_steps(0.9)


def test_rk4_three_burn_transfer_closure():
    for (r1, r2, r_b, tol) in (
        (1.0, 2.0, 5.0, 1e-4),
        (1.0, 6.41, 100.0, 1e-4),
        (1.0, 14.0, 40.0, 1e-4),
        (1.0, 20.0, 20.5, 1e-4),
        (1.0, 20.0, 200.0, 1e-4),
    ):
        r = validate_bielliptic_rk4(r1, r2, r_b)
        assert r["arrival_apoapsis_rk4"]["rel_r_error"] < tol
        assert r["arrival_apoapsis_rk4"]["rel_v_error"] < tol
        assert r["arrival_apoapsis_rk4"]["apsis_at_final"] is True
        assert r["arrival_periapsis_rk4"]["rel_r_error"] < tol
        assert r["arrival_periapsis_rk4"]["rel_v_error"] < tol
        assert r["arrival_periapsis_rk4"]["apsis_at_final"] is True
        assert r["burn2"]["rel_dv2_error"] < 1e-4
        assert r["burn3"]["rel_dv3_error"] < 1e-4
        assert r["post_burn_circular_orbit"]["radius_max_rel_variation"] < 1e-5
        assert r["post_burn_circular_orbit"][
            "speed_rel_error_vs_sqrt_mu_r2"] < 1e-5
        assert r["max_rel_drift"]["energy"] < 1e-5
        assert r["max_rel_drift"]["angular_momentum"] < 1e-6


def test_rk4_analytic_references_are_exact():
    """The closed-form legs (kepler_solution, phase-shifted for the leg that
    starts at apoapsis) must match the geometry to machine precision."""
    for (r1, r2, r_b) in ((1.0, 2.0, 5.0), (1.0, 6.41, 100.0),
                          (1.0, 14.0, 40.0), (1.0, 20.0, 20.5)):
        r = validate_bielliptic_rk4(r1, r2, r_b)
        assert r["arrival_apoapsis_analytic"]["rel_r_error"] < 1e-9
        assert r["arrival_apoapsis_analytic"]["rel_v_error"] < 1e-9
        assert r["arrival_periapsis_analytic"]["rel_r_error"] < 1e-9
        assert r["arrival_periapsis_analytic"]["rel_v_error"] < 1e-9


# --- Real-system anchors ----------------------------------------------------


def test_r0_6700_14x_reproduces_wikipedia_example():
    """Wikipedia's worked R = 14, s = 40 example (r0 = 6700 km, mu =
    398600.4418 km3/s2): dv1 = 3061.04, dv2 = 608.825, dv3 = 447.662 m/s,
    total 4117.53 vs Hohmann 4133.72 m/s - reproduced exactly here."""
    a = real_anchors()["wiki_14x_r0_6700"]
    assert 4116.5 < 1000.0 * a["dv_bielliptic_km_s"] < 4118.5
    assert 4132.5 < 1000.0 * a["dv_hohmann_km_s"] < 4135.0
    assert 3059.0 < 1000.0 * a["dv1_km_s"] < 3063.0
    assert 607.5 < 1000.0 * a["dv2_km_s"] < 610.0
    assert 446.5 < 1000.0 * a["dv3_km_s"] < 449.0
    assert 15.0 < 1000.0 * a["saving_vs_hohmann_km_s"] < 17.0
    assert 0.30 < a["saving_percent_of_hohmann"] < 0.45


def test_leo_geo_hohmann_is_cheaper():
    """R = 6.41 < R_bp: Hohmann wins; the bi-elliptic saving is negative
    and loses by hundreds of m/s."""
    a = real_anchors()["leo_geo"]
    assert a["saving_vs_hohmann_km_s"] < 0.0
    assert abs(a["saving_vs_hohmann_km_s"]) > 0.1  # > 100 m/s worse


def test_earth_mars_hohmann_is_cheaper():
    a = real_anchors()["earth_mars"]
    assert a["saving_vs_hohmann_km_s"] < 0.0


def test_geo_15_58x_is_near_break_even():
    """R = 15.58 is just below R* ~ 15.5817: with s = 30 the bi-elliptic wins
    by ~0.6% of the Hohmann cost (the crossover s_c(15.58) ~ 15.58)."""
    a = real_anchors()["geo_15_58x"]
    assert abs(a["saving_percent_of_hohmann"]) < 1.0


def test_max_saving_anchor_close_to_bi_parabolic():
    a = real_anchors()["leo_50x_best_saving"]
    # ~0.041 v1 at R ~ 50: v1(LEO) ~ 7.78 km/s -> ~320 m/s = 0.32 km/s
    assert 0.25 < a["saving_vs_hohmann_km_s"] < 0.40
    # with s = 1e6 the cost is within a few m/s of the bi-parabolic
    assert abs(a["dv_bielliptic_km_s"] - a["dv_bi_parabolic_km_s"]) < 0.005


def test_time_ratio_real_budget():
    a = real_anchors()["wiki_14x_r0_6700"]
    assert a["t_bielliptic_days"] > a["t_hohmann_days"]
    assert a["t_ratio"] > 1.0


# --- Constants sanity -------------------------------------------------------


def test_constants_sanity():
    assert abs(experiment.MU_EARTH_KM3S2 - 3.986004e5) < 1.0
    assert abs(experiment.R_EARTH_KM - 6.3781e3) < 1.0
    assert abs(experiment.AU_KM - 1.495978707e8) < 1.0
    assert 1.52 < experiment.MARS_A_AU < 1.53


# --- Determinism ------------------------------------------------------------


def test_determinism_across_processes():
    """A fresh interpreter must produce bit-identical numerical output."""
    exp_dir = Path(__file__).resolve().parents[1]
    script = (
        "import importlib.util, json, sys\n"
        f"spec = importlib.util.spec_from_file_location('exp5', "
        f"{str(exp_dir / 'experiment.py')!r})\n"
        "m = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(m)\n"
        "payload = {\n"
        "    'R_bp': m.R_bp_crossover(),\n"
        "    'R_star': m.R_always_cheaper(),\n"
        "    'sc14': m.crossover_s(14.0),\n"
        "    'rk4': m.validate_bielliptic_rk4(1.0, 14.0, 40.0)"
        "['arrival_periapsis_rk4'],\n"
        "}\n"
        "print(json.dumps(payload, sort_keys=True, default=float))\n"
    )
    out = subprocess.check_output(
        [sys.executable, "-c", script], text=True).strip()
    here = {
        "R_bp": R_bp_crossover(),
        "R_star": R_always_cheaper(),
        "sc14": crossover_s(14.0),
        "rk4": validate_bielliptic_rk4(1.0, 14.0, 40.0)[
            "arrival_periapsis_rk4"],
    }
    expected = json.dumps(here, sort_keys=True, default=float)
    assert out == expected, "results differ between interpreter processes"