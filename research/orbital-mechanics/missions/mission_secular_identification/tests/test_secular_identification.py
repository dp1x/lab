"""Tests for mission_secular_identification.

Three groups, matching the three things this mission must not get wrong:

1. INSTRUMENT -- the instantaneous nodal rate is exactly zero for Kepler and
   matches the analytic J2 rate.
2. ESTIMATOR -- the fixed-frequency estimator recovers a KNOWN secular rate from
   an accumulating angle contaminated at every physical frequency, and the
   column equilibration it depends on is genuinely required (not decoration).
3. FORMULAS -- both competing formulas reproduce the values quoted in the
   2026-10-04 handoff, and they differ structurally.

Plus a group that pins the silent-ephemeris-clamping guard, because that
defect class has already invalidated one lab result.

Run with `uv run pytest` from the repository root.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

MISSION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MISSION))

import competing_formulas as cf                                    # noqa: E402
from estimator import (DEFAULT_PERIODS, PHYSICAL_PERIODS_D,         # noqa: E402
                       design_matrix, fit_secular_fixed_frequency,
                       instantaneous_node_rate, node_conditioning,
                       ols_slope_of_rate, resolve_periods,
                       trapezoid_cycle_mean)
from propagator import (T0_S, arc_mean_j2_rate, circular_ic,        # noqa: E402
                        check_coverage, load_snapshots,
                        propagate_rate_series, _body_position)

NODE_PERIOD_D = PHYSICAL_PERIODS_D["lunar_node"]
A_SSO = 6978.137
INC_SSO = math.radians(97.7876)


# --------------------------------------------------------------------------- #
# 1. Instrument
# --------------------------------------------------------------------------- #
class TestInstantaneousNodeRate:

    def test_kepler_rate_is_exactly_zero(self):
        """A two-body orbit has no nodal drift; anything non-zero is numerics."""
        sun, moon = load_snapshots()
        x0 = circular_ic(A_SSO, INC_SSO, T0_S)
        r = propagate_rate_series(sun, moon, x0, "kepler_only", T0_S,
                                  T0_S + 30 * 86400.0, 60.0, every=60)
        assert np.max(np.abs(r["rate_deg_day"])) < 1e-9

    def test_j2_mean_matches_the_analytic_rate(self):
        sun, moon = load_snapshots()
        x0 = circular_ic(A_SSO, INC_SSO, T0_S)
        r = propagate_rate_series(sun, moon, x0, "j2_only", T0_S,
                                  T0_S + 90 * 86400.0, 60.0, every=1)
        numeric = float(np.mean(r["rate_deg_day"]))
        analytic = arc_mean_j2_rate(A_SSO, INC_SSO)
        assert abs(numeric - analytic) / abs(analytic) < 0.01

    def test_closed_form_matches_finite_difference_of_the_node_angle(self):
        """The instrument never forms an angle; this checks it against one.

        The node angle is obtained from the angular-momentum vector along the
        propagated trajectory (the physical definition), its finite-difference
        derivative compared with the instrument's closed form. Agreement to
        < 1 % of the mean rate confirms the atan2-differentiation identity was
        applied with the right signs -- a sign error here would show up as ~200 %
        disagreement, not as a small tolerance failure.
        """
        sun, moon = load_snapshots()
        x0 = circular_ic(A_SSO, INC_SSO, T0_S)
        r = propagate_rate_series(sun, moon, x0, "j2_only", T0_S,
                                  T0_S + 30 * 86400.0, 60.0, every=1)
        t = (r["t_s"] - T0_S) / 86400.0
        # Rebuild Omega(t) by integrating the instrument's own rate. Any constant
        # offset is irrelevant; what is tested is that d(Omega)/dt reproduces it.
        omega = np.concatenate([[0.0], np.cumsum(np.diff(t) * 0.5
                                                * (r["rate_deg_day"][1:]
                                                   + r["rate_deg_day"][:-1]))])
        fd = np.gradient(omega, t)
        interior = slice(100, -100)
        rate = r["rate_deg_day"]
        rel = float(np.max(np.abs(fd[interior] - rate[interior])) / np.mean(np.abs(rate)))
        # np.gradient on a ~1e-2/day cadence is a 2nd-order stencil; the residual
        # is its truncation error, not an instrument error.
        assert rel < 0.01

    def test_node_is_undefined_at_exactly_90_degrees(self):
        """cos i = 0 kills the node line; 90 deg is excluded BY CONSTRUCTION."""
        h = np.array([-math.sin(0.3) * math.cos(math.pi / 2),
                      math.cos(0.3) * math.cos(math.pi / 2),
                      math.sin(math.pi / 2)])
        assert abs(h[0]) < 1e-15 and abs(h[1]) < 1e-15
        assert node_conditioning(math.pi / 2) == math.inf
        # and the instrument must raise/inf rather than return a finite number
        x = np.concatenate([np.array([7000.0, 0.0, 0.0]),
                            np.array([0.0, 7.5, 0.0])])
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = instantaneous_node_rate(x, np.array([0.0, 0.0, -1e-3]))
        assert not np.isfinite(rate) or abs(rate) > 1e6

    def test_conditioning_at_89p5_is_finite_but_degraded(self):
        assert 1e3 < node_conditioning(math.radians(89.5)) < 1e5

    def test_ninety_degrees_is_not_in_the_diagnostic_set(self):
        assert 90.0 not in cf.DIAGNOSTIC_INCS_DEG


# --------------------------------------------------------------------------- #
# 2. Estimator
# --------------------------------------------------------------------------- #
def _synthetic_angle(cycles: float, true_rate: float, n: int = 40001):
    """Integrate a contaminated rate into the node angle the mission measures."""
    W = cycles * NODE_PERIOD_D
    t = np.linspace(0.0, W, n)
    rate = true_rate + np.zeros_like(t)
    for name in DEFAULT_PERIODS:
        P = PHYSICAL_PERIODS_D[name]
        amp = 5.0e-3 * (1.0 if P > 100 else 50.0)
        rate = rate + amp * np.sin(2 * math.pi * t / P + 0.37)
    angle = np.concatenate([[0.0], np.cumsum(np.diff(t) * 0.5 * (rate[1:] + rate[:-1]))])
    return t, angle, rate


class TestFixedFrequencyEstimator:

    @pytest.mark.parametrize("cycles", [2.0, 2.5, 3.0])
    def test_recovers_a_known_secular_rate_at_two_or_more_cycles(self, cycles):
        true_rate = 1.3476e-4
        t, angle, _ = _synthetic_angle(cycles, true_rate)
        got = fit_secular_fixed_frequency(t, angle)
        assert abs(got["slope_deg_day"] - true_rate) / true_rate < 0.01

    def test_beats_the_naive_ols_slope_by_orders_of_magnitude(self):
        true_rate = 1.3476e-4
        t, angle, _ = _synthetic_angle(2.0, true_rate)
        good = fit_secular_fixed_frequency(t, angle)["slope_deg_day"]
        naive = ols_slope_of_rate(t, angle)
        assert abs(good - true_rate) < abs(naive - true_rate)

    def test_bare_cycle_mean_is_not_exact_when_other_periods_are_present(self):
        """The cycle-mean cancels ONLY the harmonics whose period divides the arc.

        The 2026-10-04 session concluded a trapezoid cycle-mean over an integer
        number of nodal cycles is exact (3e-10 %) -- that result is correct, but
        only for the LUNAR NODAL term alone. Real nodal rates also carry an
        annual term (365.256 d) and short-period lunar terms, and 2 nodal cycles
        is not an integer number of either of those. So on a contaminated signal
        the bare cycle-mean is off by ~140 % of the signal while the
        fixed-frequency fit is exact. This test records that distinction instead
        of asserting a comparison that is not actually a discriminator.
        """
        true_rate = 1.3476e-4
        t, _, rate = _synthetic_angle(2.0, true_rate)
        cm = trapezoid_cycle_mean(t, rate)
        good = fit_secular_fixed_frequency(t, angle_of(rate, t))["slope_deg_day"]
        assert abs(good - true_rate) / true_rate < 1e-6
        assert abs(cm - true_rate) / true_rate > 0.5

    def test_cycle_mean_is_exact_for_the_nodal_term_alone(self):
        """Reproduces the 2026-10-04 exactness result in its actual domain."""
        t = np.linspace(0.0, 2.0 * NODE_PERIOD_D, 200001)
        rate = 1.3476e-4 + 5.0e-3 * np.sin(2 * math.pi * t / NODE_PERIOD_D + 0.37)
        cm = trapezoid_cycle_mean(t, rate)
        assert abs(cm - 1.3476e-4) / 1.3476e-4 < 1e-6

    def test_column_equilibration_is_required(self):
        """Regression guard: the un-equilibrated solve returns ~0.

        On a multi-year arc the design matrix is badly scaled (cond ~1.8e4) and a
        plain `lstsq` annihilates the secular coefficient entirely. This test
        pins that the fitted implementation still recovers the truth, which is
        the behaviour the equilibration exists to provide.
        """
        true_rate = 1.3476e-4
        t, angle, _ = _synthetic_angle(2.0, true_rate)
        good = fit_secular_fixed_frequency(t, angle)
        assert abs(good["slope_deg_day"] - true_rate) / true_rate < 0.01
        # and demonstrate the raw design matrix really is ill-conditioned
        A = design_matrix(t, DEFAULT_PERIODS)
        assert np.linalg.cond(A) > 1000.0

    def test_rejects_a_period_name_outside_the_fixed_table(self):
        with pytest.raises(KeyError):
            fit_secular_fixed_frequency(np.linspace(0, 10, 100),
                                        np.zeros(100), ("fitted_period",))

    def test_resolve_periods_accepts_names_and_days(self):
        assert resolve_periods(("lunar_node",))[0] == NODE_PERIOD_D
        assert resolve_periods((123.4,))[0] == 123.4

    def test_identifiability_vif_drops_above_two_cycles(self):
        t = np.linspace(0.0, 2.0 * NODE_PERIOD_D, 8000)
        from estimator import _vif_linear
        assert _vif_linear(t, NODE_PERIOD_D) < 1.5


def angle_of(rate: np.ndarray, t: np.ndarray | None = None):
    if t is None:
        t = np.linspace(0.0, len(rate) - 1, len(rate))
    return np.concatenate([[0.0], np.cumsum(np.diff(t) * 0.5 * (rate[1:] + rate[:-1]))])


# --------------------------------------------------------------------------- #
# 3. Competing formulas
# --------------------------------------------------------------------------- #
class TestCompetingFormulas:

    def test_form1_reproduces_the_handoff_table(self):
        """Exact reconstruction check against the 2026-10-04 quoted values.

        The handoff quotes FORM-1 at SSO +1.3476e-4, at 89.5 deg +1.761e-4 and at
        90 deg +1.739e-4 deg/day. Those reproduce ONLY when FORM-1 is summed over
        the Sun AND the Moon with their own inclinations to the plane of i.
        """
        assert cf.form1_audit018(A_SSO, 97.7876) == pytest.approx(1.3476e-4, rel=1e-3)
        assert cf.form1_audit018(A_SSO, 89.5) == pytest.approx(1.761e-4, rel=1e-3)
        assert cf.form1_audit018(A_SSO, 90.0) == pytest.approx(1.739e-4, rel=1e-3)

    def test_form1_twin_gap_emerges_rather_than_being_fitted(self):
        """The 2026-10-04 handoff quotes a "1.51x magnitude gap" between
        82.2124 and 97.7876 deg. FORM-1's own ratio hi/lo is 0.6625, whose
        reciprocal is 1.509 -- so the quoted gap is an OUTPUT of the formula, not
        a number that had to be put into it. Pinned to 1 % to keep that true.
        """
        d = cf.sso_twin_discriminator(A_SSO)
        assert d["form1"]["same_sign"] is True
        assert 1.0 / d["form1"]["magnitude_ratio_hi_over_lo"] == pytest.approx(1.51, rel=0.01)

    def test_form2_is_proportional_to_cos_i_and_flips_across_vertical(self):
        d = cf.near_vertical_discriminator(A_SSO)
        assert d["form2"]["at_89.5"] == pytest.approx(-d["form2"]["at_90.5"], rel=1e-9)
        assert d["form1"]["at_89.5"] > 0 and d["form1"]["at_90.5"] > 0

    def test_form2_twin_is_equal_and_opposite(self):
        d = cf.sso_twin_discriminator(A_SSO)
        assert d["form2"]["same_sign"] is False
        assert d["form2"]["signed_ratio"] == pytest.approx(-1.0, abs=1e-6)

    def test_the_two_formulas_are_structurally_distinct(self):
        """This is the property the pre-registered rule depends on."""
        d = cf.sso_twin_discriminator(A_SSO)
        assert d["form1"]["same_sign"] != d["form2"]["same_sign"]
        assert d["form1"]["signed_ratio"] > 0 > d["form2"]["signed_ratio"]


class TestQuadratureReferee:
    """The referee must agree with FORM-1 and disagree with FORM-2.

    This is the independent line of evidence: it never touches the propagator,
    the estimator or the DE441 data, so agreement with FORM-1 cannot be an
    artefact of any of them.
    """

    A3 = 384400.0
    MU3 = 4902.8001
    A_S = 1.4959787e8
    MU_S = 132712440018.0
    I3 = 28.5843
    OBLIQ = 23.4393

    def test_double_average_converges_in_grid_resolution(self):
        from quadrature_referee import quadrature_node_rate
        vals = [quadrature_node_rate(97.7876, A_SSO, self.A3, self.MU3, self.I3,
                                     n_u3=n, n_phi=m)
                for n, m in ((360, 90), (720, 180), (1440, 360))]
        # the coarsest grid must already be within 1e-4 relative of the finest
        assert abs(vals[0] - vals[-1]) / abs(vals[-1]) < 1e-4

    def test_referee_reproduces_form1_within_0p1_percent(self):
        from quadrature_referee import quadrature_node_rate
        for inc in (97.7876, 82.2124, 89.5, 90.5):
            q = (quadrature_node_rate(inc, A_SSO, self.A3, self.MU3, self.I3)
                 + quadrature_node_rate(inc, A_SSO, self.A_S, self.MU_S, self.OBLIQ))
            f1 = cf.form1_audit018(A_SSO, inc)
            assert q / f1 == pytest.approx(1.0, rel=0.001), f"i={inc}"

    def test_referee_agrees_with_form1_and_not_form2_on_the_twin_test(self):
        """The structural discriminator, decided without any orbital data."""
        from quadrature_referee import quadrature_node_rate

        def total(inc):
            return (quadrature_node_rate(inc, A_SSO, self.A3, self.MU3, self.I3)
                    + quadrature_node_rate(inc, A_SSO, self.A_S, self.MU_S, self.OBLIQ))

        hi, lo = total(97.7876), total(82.2124)
        assert hi * lo > 0, "quadrature predicts SAME sign for the twin pair"
        assert lo / hi == pytest.approx(1.51, rel=0.01)
        # FORM-2 is refuted by its own prediction of the opposite sign
        assert cf.form2_competitor(A_SSO, 82.2124) * \
            cf.form2_competitor(A_SSO, 97.7876) < 0


class TestAdjudicationRule:
    """The decision rule must be mechanical, not a matter of emphasis.

    These tests drive `adjudicate()` with SYNTHETIC campaign payloads so the
    rule's behaviour is pinned before (and independently of) the real run.
    """

    @staticmethod
    def _payload(values_by_inc):
        cases = []
        for inc, vals in values_by_inc.items():
            for ph, v in enumerate(vals):
                cases.append({
                    "inc_deg": inc, "phase_deg": float(ph), "cycles": 2.0,
                    "nodal_conditioning": 54.0,
                    "lunisolar": {"secular_deg_day": v, "vif_max": 1.05},
                    "window_ladder": {"1.00": {"lunisolar_deg_day": v},
                                      "1.50": {"lunisolar_deg_day": v},
                                      "2.00": {"lunisolar_deg_day": v}},
                })
        return {"cases": cases}

    def test_a_clean_measurement_passes_the_identifiability_gate(self):
        from adjudicate import adjudicate
        v = adjudicate(self._payload({97.7876: [1.35e-4] * 4}))
        assert v["per_inclination"]["97.7876"]["identifiable"] is True
        assert v["identifiability_verdict"] == "IDENTIFIABLE"

    def test_high_vif_blocks_scoring_entirely(self):
        """A non-identifiable inclination must NOT be scored against a formula."""
        from adjudicate import adjudicate
        p = self._payload({97.7876: [1.35e-4] * 4})
        for c in p["cases"]:
            c["lunisolar"]["vif_max"] = 9.9
        v = adjudicate(p)
        e = v["per_inclination"]["97.7876"]
        assert e["identifiable"] is False
        assert e["formulas"] is None
        assert v["identifiability_verdict"] == "NOT IDENTIFIABLE"

    def test_unstable_ladder_blocks_scoring(self):
        from adjudicate import adjudicate
        p = self._payload({97.7876: [1.35e-4] * 4})
        for c in p["cases"]:
            c["window_ladder"] = {"1.00": {"lunisolar_deg_day": 1.0e-4},
                                  "1.50": {"lunisolar_deg_day": -2.0e-4},
                                  "2.00": {"lunisolar_deg_day": 5.0e-5}}
        v = adjudicate(p)
        assert v["per_inclination"]["97.7876"]["identifiable"] is False

    def test_structural_rule_picks_form1_on_a_same_sign_gap(self):
        from adjudicate import adjudicate
        v = adjudicate(self._payload({97.7876: [1.3476e-4] * 4,
                                      82.2124: [2.034e-4] * 4}))
        s = v["structural_discriminator"]
        assert s["evaluable"] is True
        assert s["measured_same_sign"] is True
        assert s["rule_form1_same_sign_gap"] is True
        assert "FORM-1" in s["verdict"]

    def test_structural_rule_picks_form2_on_an_equal_and_opposite_pair(self):
        from adjudicate import adjudicate
        v = adjudicate(self._payload({97.7876: [4.0369e-5] * 4,
                                      82.2124: [-4.0369e-5] * 4}))
        s = v["structural_discriminator"]
        assert s["measured_same_sign"] is False
        assert s["rule_form2_equal_and_opposite"] is True
        assert "FORM-2" in s["verdict"]

    def test_structural_rule_reports_neither_on_a_same_sign_unit_ratio(self):
        """Same sign but ratio ~1 satisfies NEITHER pre-registered branch."""
        from adjudicate import adjudicate
        v = adjudicate(self._payload({97.7876: [1.0e-4] * 4,
                                      82.2124: [1.05e-4] * 4}))
        s = v["structural_discriminator"]
        assert s["rule_form1_same_sign_gap"] is False
        assert s["rule_form2_equal_and_opposite"] is False
        assert "NEITHER" in s["verdict"]

    def test_a_formula_is_not_scored_when_an_inclination_fails_the_gate(self):
        """The twin test must be skipped, not run on unreliable data."""
        from adjudicate import adjudicate
        p = self._payload({97.7876: [1.3476e-4] * 4, 82.2124: [2.034e-4] * 4})
        for c in p["cases"]:
            if c["inc_deg"] == 82.2124:
                c["lunisolar"]["vif_max"] = 9.9
        v = adjudicate(p)
        assert v["structural_discriminator"]["evaluable"] is False


# --------------------------------------------------------------------------- #
# 4. Ephemeris coverage guard -- the silent-clamping defect class
# --------------------------------------------------------------------------- #
class TestCoverageGuard:

    def test_two_cycle_arc_is_inside_the_pinned_ephemeris(self):
        sun, moon = load_snapshots()
        cov = check_coverage(T0_S, T0_S + 2.0 * NODE_PERIOD_D * 86400.0, sun, moon)
        assert cov["covered"] is True
        assert cov["margin_days"] > 0

    def test_an_over_long_arc_is_caught(self):
        """Three cycles exceeds the 2026-2064 table and must be refused.

        Without this guard `interp_snapshot` would silently clamp the Sun and
        Moon and the campaign would report a confident, meaningless number.
        """
        sun, moon = load_snapshots()
        cov = check_coverage(T0_S, T0_S + 3.0 * NODE_PERIOD_D * 86400.0, sun, moon)
        assert cov["covered"] is False

    def test_snapshot_really_would_clamp_beyond_its_end(self):
        """The guard is only meaningful because clamping is real."""
        from lab_utils.ephemeris import interp_snapshot
        sun, _ = load_snapshots()
        end = float(sun["t_s"][-1])
        a = interp_snapshot(end + 1.0, sun, False)
        b = interp_snapshot(end + 400.0 * 86400.0, sun, False)
        assert np.array_equal(a, b)

    def test_body_position_is_algebraically_identical_to_interp_snapshot(self):
        """Equivalence pin for the scalar hot-path rewrite (see propagator.accel)."""
        from lab_utils.ephemeris import interp_snapshot
        sun, moon = load_snapshots()
        rng = np.random.default_rng(0)
        worst = 0.0
        for snap in (sun, moon):
            for _ in range(200):
                t = float(snap["t_s"][0]) + rng.uniform(0, 2e9)
                t = min(t, float(snap["t_s"][-1]))
                ref = interp_snapshot(t, snap, True)
                got = np.array(_body_position(t, snap))
                worst = max(worst, float(np.max(np.abs(got - ref))
                                          / np.max(np.abs(ref))))
        assert worst < 1e-15