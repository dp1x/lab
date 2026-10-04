"""Regression tests pinning the estimator-identifiability findings of 2026-10-04.

These are pure linear-algebra / quadrature properties with NO dependence on the
force model, so they are stable tests that will not need re-derivation.

What they pin:
  1. The committed 18.6-yr arc is 0.9993 lunar nodal cycles, NOT the 0.894
     claimed in localdocs/reports/mission-mean-element-leakage-2026-10-03.md.
  2. On that arc, a joint [1, t, cos(wt), sin(wt)] fit has a variance inflation
     factor of ~11 on the linear column for any period >= 6798 d, and the fitted
     secular slope is therefore not identifiable.
  3. The same fit on a 2-cycle arc has VIF ~ 1.01: identifiable.
  4. A trapezoid cycle-mean of a nodal-frequency rate is exact at integer cycles
     and has offset ~ 1/W away from them; np.mean carries an O(1/N) endpoint
     bias that must not be mistaken for physics.

Reconstructed after the R: scratch loss of 2026-10-04; see
localdocs/reports/mission-mean-element-leakage-2026-10-04.md.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

NODE_PERIOD_D = 6798.383          # lunar nodal regression period
COMMITTED_ARC_D = 18.6 * 365.25   # what `--years 18.6` actually produces


def _vif_linear(t: np.ndarray, period: float) -> float:
    """Variance inflation factor of the linear column against a harmonic pair."""
    w = 2.0 * math.pi / period
    A = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    coef, *_ = np.linalg.lstsq(A, t, rcond=None)
    r2 = 1.0 - (t - A @ coef).var() / t.var()
    return 1.0 / (1.0 - r2)


def _trapz_mean(f, n: int, window: float) -> float:
    t = np.linspace(0.0, window, n)
    fv = f(t)
    dt = window / (n - 1)
    return float(dt * (fv[0] / 2.0 + fv[1:-1].sum() + fv[-1] / 2.0) / window)


class TestCommittedArcLength:
    def test_committed_arc_is_one_nodal_cycle_not_0894(self):
        cycles = COMMITTED_ARC_D / NODE_PERIOD_D
        assert cycles == pytest.approx(0.9993, abs=2e-4)
        # the 2026-10-03 report states 0.894, which is inconsistent with 6793.6 d
        assert abs(6793.6 / NODE_PERIOD_D - 0.894) > 0.09

    def test_two_cycles_exceeds_the_pinned_snapshot(self):
        from lab_utils.ephemeris import de441_snapshot_paths, load_snapshot
        from pathlib import Path
        repo = Path(__file__).resolve().parents[5]
        sun_p, _ = de441_snapshot_paths(repo)
        snap = load_snapshot(sun_p)
        t0 = 820476800.0
        one_cycle_end = t0 + NODE_PERIOD_D * 86400.0
        two_cycle_end = t0 + 2.0 * NODE_PERIOD_D * 86400.0
        assert one_cycle_end < snap["t_s"][-1]      # 1 cycle fits
        assert two_cycle_end > snap["t_s"][-1]      # 2 cycles does not fit


class TestSnapshotClamping:
    def test_sun_vector_is_frozen_beyond_the_snapshot_end(self):
        from lab_utils.ephemeris import de441_snapshot_paths, interp_snapshot, load_snapshot
        from pathlib import Path
        repo = Path(__file__).resolve().parents[5]
        sun_p, _ = de441_snapshot_paths(repo)
        snap = load_snapshot(sun_p)
        end = float(snap["t_s"][-1])
        raw_after = interp_snapshot(end + 1.0, snap, False)
        raw_far = interp_snapshot(end + 400.0 * 86400.0, snap, False)
        assert np.array_equal(raw_after, raw_far)

    def test_clamping_is_invisible_after_precession(self):
        """The precessed comparison is non-zero even though the vector is frozen.

        Pinned because it caused a false 'not clamped' reading during the audit:
        the precession rotation depends on t, so a frozen J2000 vector still
        rotates. The raw (apply_precession=False) comparison is the correct test.
        """
        from lab_utils.ephemeris import de441_snapshot_paths, interp_snapshot, load_snapshot
        from pathlib import Path
        repo = Path(__file__).resolve().parents[5]
        sun_p, _ = de441_snapshot_paths(repo)
        snap = load_snapshot(sun_p)
        end = float(snap["t_s"][-1])
        pre_a = interp_snapshot(end + 1.0, snap, True)
        pre_b = interp_snapshot(end + 400.0 * 86400.0, snap, True)
        assert not np.array_equal(pre_a, pre_b)   # precessed -> looks unclamped


class TestIdentifiability:
    def test_vif_is_large_on_the_committed_arc(self):
        t = np.linspace(0.0, COMMITTED_ARC_D, 8000)
        assert _vif_linear(t, 9175.0) > 5.0
        assert _vif_linear(t, 9000.0) > 5.0

    def test_vif_is_small_on_a_two_cycle_arc(self):
        t = np.linspace(0.0, 2.0 * NODE_PERIOD_D, 8000)
        assert _vif_linear(t, 9175.0) < 1.5

    def test_slope_is_not_recoverable_from_a_degenerate_joint_fit(self):
        """Inject a known secular slope + known long-period term, then refit at a
        WRONG period: the fitted slope is wrong by >10x even though the fit
        explains >99.9% of the variance."""
        t = np.linspace(0.0, COMMITTED_ARC_D, 8000)
        s0 = -8.564733e-4
        series = s0 * t + 35.0 * np.sin(2.0 * math.pi * t / 9175.0 + 0.7)
        worst = 0.0
        for P in (6798.383, 8000.0, 9000.0, 9500.0, 12000.0):
            w = 2.0 * math.pi / P
            A = np.column_stack([np.ones_like(t), t, np.cos(w * t), np.sin(w * t)])
            coef, *_ = np.linalg.lstsq(A, series, rcond=None)
            resid = series - A @ coef
            frac_var = 1.0 - resid.var() / series.var()
            assert frac_var > 0.999          # "explains" essentially everything
            worst = max(worst, abs(coef[1] - s0) / abs(s0))
        assert worst > 10.0                 # ... yet the slope is unrecoverable


class TestCycleMeanEstimator:
    def test_trapezoid_cycle_mean_is_exact_at_integer_cycles(self):
        w = 2.0 * math.pi / NODE_PERIOD_D
        s0, amp, phase = 4.0e-7, 2.4e-2, 0.7
        rate = lambda t: s0 + amp * np.sin(w * t + phase)
        for k in (1, 2, 4):
            got = _trapz_mean(rate, 200001, k * NODE_PERIOD_D)
            assert abs(got - s0) / s0 < 1e-8

    def test_np_mean_carries_an_endpoint_bias_that_must_not_be_read_as_physics(self):
        w = 2.0 * math.pi / NODE_PERIOD_D
        s0, amp, phase = 4.0e-7, 2.4e-2, 0.7
        rate = lambda t: s0 + amp * np.sin(w * t + phase)
        W = NODE_PERIOD_D
        naive = float(np.mean(rate(np.linspace(0.0, W, 10001))))
        assert abs(naive - s0) / s0 > 1.0        # np.mean is badly biased here
        good = _trapz_mean(rate, 10001, W)
        assert abs(good - s0) / s0 < 1e-8         # trapezoid is exact

    def test_offset_grows_with_distance_from_an_integer_cycle(self):
        w = 2.0 * math.pi / NODE_PERIOD_D
        s0, amp, phase = 4.0e-7, 2.4e-2, 0.7
        rate = lambda t: s0 + amp * np.sin(w * t + phase)
        offsets = []
        for k in (0.99, 0.999, 0.9999):
            off = _trapz_mean(rate, 200001, k * NODE_PERIOD_D) - s0
            offsets.append(abs(off) / s0)
        assert offsets[0] > offsets[1] > offsets[2]   # 1/W law
        # the committed 18.6-yr arc is ~27x the signal away from an integer cycle
        off = _trapz_mean(rate, 200001, COMMITTED_ARC_D) - s0
        assert abs(off) / s0 > 10.0


class TestNodeConditioning:
    @pytest.mark.parametrize("inc_deg,expect_bad", [(30.0, False), (97.7876, False),
                                                   (89.5, True)])
    def test_node_longitude_conditioning_grows_as_inclination_approaches_90(self, inc_deg,
                                                                          expect_bad):
        i = math.radians(inc_deg)
        cond = 1.0 / math.cos(i) ** 2
        assert (cond > 1e4) is expect_bad

    def test_node_is_undefined_at_exactly_90_degrees(self):
        i = math.radians(90.0)
        h = np.array([-math.sin(0.3) * math.cos(i), math.cos(0.3) * math.cos(i),
                      math.sin(i)])
        assert h[0] == pytest.approx(0.0, abs=1e-15)
        assert h[1] == pytest.approx(0.0, abs=1e-15)
        # atan2(0, 0) is 0 in numpy but the quantity is undefined
        assert math.cos(i) == pytest.approx(0.0, abs=1e-15)

    def test_convention_is_raan_modulo_180(self):
        """atan2(-h_x, h_y) reproduces the RAAN for cos i > 0 and RAAN+180 for
        cos i < 0, so a constant offset cancels in any slope estimate.

        Derived from the general identity atan2(sin x, cos x) = x mod 2pi, applied
        to h_x = -sinO cos i, h_y = cosO cos i -> the argument is (O) when cos i>0
        and (O + pi) when cos i<0.
        """
        for inc_deg in (30.0, 60.0, 97.7876, 120.0):
            i, O = math.radians(inc_deg), math.radians(33.0)
            ci = math.cos(i)
            hx = -math.sin(O) * ci
            hy = math.cos(O) * ci
            got = math.atan2(-hx, hy)
            # atan2 scales by the sign of the common factor ci
            expect = O if ci > 0 else O + math.pi
            d = (got - expect + math.pi) % (2 * math.pi) - math.pi
            assert abs(d) < 1e-12
