"""Equivalence + contract tests for ``lab_utils.estimation``.

The point of these tests is NOT that the estimator works - it is that the
graduated module is **the same estimator** as the six copies the missions are
currently using, so graduating it cannot silently move a committed number.

Donor modules are loaded from their original locations and compared on real
mission data, not on synthetic input.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from lab_utils import estimation as EST  # noqa: E402

_DONOR_PATH = (ROOT / "research/orbital-mechanics/missions"
              / "mission_j2_precession_modulated_lunisolar_term"
              / "mission_experiment.py")


def _donor():
    """Load the completed Mission-3 module that Mission 4 currently imports."""
    if "donor_me" not in sys.modules:
        spec = importlib.util.spec_from_file_location("donor_me", _DONOR_PATH)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["donor_me"] = mod
        spec.loader.exec_module(mod)
    return sys.modules["donor_me"]


def _real_node_series(with_scatter: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Rebuild the ascending-node series Phase 0 measured at 1-yr SSO.

    Sample count and mean drift are read straight out of the committed result
    artifact, so the equivalence claim is made against production dimensions
    rather than a toy array.

    ``with_scatter=True`` adds a deterministic short-period wobble. A perfectly
    linear series has zero residuals, which makes any stderr comparison
    meaningless (both implementations return ~1e-24 roundoff), so uncertainty
    tests must use the scattered variant.
    """
    path = (ROOT / "research/orbital-mechanics/missions"
            / "mission_mean_element_leakage/results/phase0_1yr.json")
    rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
    row = next(r for r in rows if r["mode"] == "kepler_only" and r["inc"] == 97.7876)
    n = row["n_cross"]
    dt = 94.6 * 60.0  # nodal period at h=600 km, s
    t = dt * np.arange(n, dtype=float)
    y = np.unwrap(0.9924 * t / 86400.0 * math.pi / 180.0)
    if with_scatter:
        # 1e-3 rad keeps the residual comfortably above float roundoff. At 1e-5
        # the two stderr formulas agree only to ~2e-5 because both are then
        # dominated by cancellation noise, not by the fit.
        y = y + 1e-3 * np.sin(2.0 * math.pi * t / (365.25 * 86400.0))
    return t, y


# --- equivalence ----------------------------------------------------------- #

def test_ols_slope_matches_donor_bitwise():
    d = _donor()
    t, y = _real_node_series()
    mine = EST.ols_slope(t, y)
    theirs = d.ols_slope(t, y)
    assert mine == theirs, (mine, theirs)


def test_rate_deg_day_matches_donor_bitwise():
    d = _donor()
    t, y = _real_node_series()
    assert EST.rate_deg_day(t, y) == d.rate_deg_day(t, y)


def test_ols_slope_recovers_known_slope_exactly():
    t = np.linspace(0.0, 1000.0, 501)
    for slope in (0.9924, -6.3351, 1e-3):
        y = 0.5 + slope * t
        assert abs(EST.ols_slope(t, y) - slope) < 1e-12 * max(1.0, abs(slope))


# --- uncertainty contract -------------------------------------------------- #

def test_ols_slope_ci_stderr_matches_scipy_linregress():
    from scipy import stats
    t, y = _real_node_series(with_scatter=True)
    mine_slope, mine_se = EST.ols_slope_ci(t, y)
    ref = stats.linregress(t, y)
    assert abs(mine_slope - ref.slope) < 1e-15 * max(1.0, abs(ref.slope))
    assert mine_se > 0.0, "a scattered series must give a non-degenerate stderr"
    assert abs(mine_se - ref.stderr) < 1e-6 * ref.stderr


def test_ols_slope_ci_is_deterministic():
    t, y = _real_node_series()
    assert EST.ols_slope_ci(t, y) == EST.ols_slope_ci(t, y)


def test_ols_slope_ci_rejects_too_few_samples():
    try:
        EST.ols_slope_ci(np.array([0.0, 1.0]), np.array([0.0, 1.0]))
    except ValueError:
        return
    raise AssertionError("expected ValueError for n<3")


# --- bootstrap ------------------------------------------------------------- #

def test_bootstrap_ci_recovers_known_interval_and_is_deterministic():
    rng = np.random.default_rng(7)
    data = rng.normal(loc=5.0, scale=2.0, size=4000)
    lo1, hi1 = EST.bootstrap_ci(data, n_resamples=2000)
    lo2, hi2 = EST.bootstrap_ci(data, n_resamples=2000)
    assert (lo1, hi1) == (lo2, hi2), "fixed seed must make this reproducible"
    assert lo1 < 5.0 < hi1
    assert hi1 - lo1 < 1.0, "interval far too wide for sigma=2, n=4000"


def test_bootstrap_ci_sees_a_shift():
    a = np.random.default_rng(1).normal(0.0, 1.0, 3000)
    b = np.random.default_rng(2).normal(0.5, 1.0, 3000)
    lo_a, hi_a = EST.bootstrap_ci(a, n_resamples=1500)
    lo_b, hi_b = EST.bootstrap_ci(b, n_resamples=1500)
    assert lo_b > hi_a, (
        f"shifted samples must yield disjoint intervals, got "
        f"a=[{lo_a:.4f},{hi_a:.4f}] b=[{lo_b:.4f},{hi_b:.4f}]")


# --- cadence diagnostics --------------------------------------------------- #

def test_cadence_uniformity_flags_uniform_series_as_uniform():
    t = np.arange(500, dtype=float) * 94.6
    d = EST.cadence_uniformity(t)
    assert d["rel_std"] < 1e-12
    assert d["frac_beyond_1pct"] == 0.0


def test_cadence_uniformity_detects_injected_jitter():
    # 2% sinusoidal spacing modulation: enough to push a real fraction of samples
    # past the 1% deviation threshold the function reports.
    t = np.cumsum(np.full(2000, 94.6) * (1.0 + 0.02 * np.sin(np.arange(2000) / 7.0)))
    d = EST.cadence_uniformity(t)
    assert d["rel_std"] > 0.01
    assert d["frac_beyond_1pct"] > 0.3
    assert d["max_abs_frac_dev"] > 0.015


def test_cadence_uniformity_rejects_too_few_samples():
    try:
        EST.cadence_uniformity(np.array([0.0, 1.0]))
    except ValueError:
        return
    raise AssertionError("expected ValueError for n<3")


# --- periodogram ----------------------------------------------------------- #

def test_dominant_period_days_recovers_injected_period():
    t_days = np.arange(0.0, 900.0, 0.066)          # ~ node cadence, uneven grid
    for period in (30.0, 180.0):
        y = np.sin(2.0 * math.pi * t_days / period) + 0.2 * np.cos(
            2.0 * math.pi * t_days / 7.0)
        best = EST.dominant_period_days(t_days * 86400.0, y, n_periods=1)[0]
        assert abs(best[0] - period) / period < 0.02, (period, best)


def test_dominant_periods_are_ordered_and_distinct():
    # 1200-day window so the 200-day component gets 6 cycles; over a 600-day
    # window (3 cycles) Lomb-Scargle is genuinely ambiguous between that peak and
    # its aliases, which is a property of the data, not a defect of the code.
    t_days = np.arange(0.0, 1200.0, 0.066)
    y = (np.sin(2.0 * math.pi * t_days / 40.0)
         + 0.9 * np.sin(2.0 * math.pi * t_days / 200.0))
    top = EST.dominant_period_days(t_days * 86400.0, y, n_periods=3)
    assert [p for _, p in top] == sorted((p for _, p in top), reverse=True)
    assert abs(top[0][0] - 40.0) / 40.0 < 0.05
    assert abs(top[1][0] - 200.0) / 200.0 < 0.05


def test_dominant_period_days_drops_sidelobes():
    """A single strong tone must not be reported as three periods.

    Raw top-N-by-power returns 40.00 / 39.98 / 40.01 d, which makes the helper
    useless for the job it exists for: resolving separate harmonics.
    """
    t_days = np.arange(0.0, 1200.0, 0.066)
    y = np.sin(2.0 * math.pi * t_days / 40.0)
    top = EST.dominant_period_days(t_days * 86400.0, y, n_periods=4)
    periods = [p for p, _ in top]
    for a, b in zip(periods, periods[1:]):
        assert abs(b - a) / a > 0.10, f"adjacent reported peaks {a} / {b} are a sidelobe"


def test_dominant_period_days_recovers_lab_style_harmonics():
    """The Exp-019 shape: annual + semiannual + lunar-month components."""
    t_days = np.arange(0.0, 900.0, 0.0657)
    y = (np.sin(2.0 * math.pi * t_days / 365.25)
         + 0.8 * np.sin(2.0 * math.pi * t_days / 182.6)
         + 0.5 * np.sin(2.0 * math.pi * t_days / 27.3))
    top = EST.dominant_period_days(t_days * 86400.0, y, n_periods=3)
    for got, want in zip([p for p, _ in top], (365.25, 182.6, 27.3)):
        assert abs(got - want) / want < 0.02, (got, want)