"""Focused tests for mission_j2_precession_modulated_lunisolar_term (no empirical fit)."""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
HERE = Path(__file__).resolve().parents[1]
# Self-first insertion (same convention as sibling missions): when the full
# repository suite runs from a different CWD, sibling mission directories on
# sys.path can shadow THIS mission's modules with stale copies. Inserting
# HERE first pins imports to the module under test.
sys.path.insert(0, str(HERE))
# Evict any cached sibling mission module: mission_j2_lunisolar_coupling also
# ships mission_experiment.py. Whichever test file imports first wins
# sys.modules, and a path-insert alone cannot evict it -> AttributeError
# storms in the full suite. Re-import pinned to HERE after the eviction.
sys.modules.pop("mission_experiment", None)
import averaged_propagator as AV
import mission_experiment as ME
import theory_modulation as TH


def test_correct_interp_midpoint_is_average():
    import tempfile
    # synthetic snap: two points 0h and 24h
    snap = {"t_s": np.array([0.0, 86400.0]),
            "r_eci_km": np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])}
    mid = ME.interp_snapshot(43200.0, snap, apply_precession=False)
    assert np.allclose(mid, [0.5, 0.5, 0.0], atol=1e-12)
    buggy = ME.interp_snapshot_buggy(43200.0, snap, apply_precession=False)
    # buggy gives [1.5,-0.5,0]: must NOT equal average (guards the bug class)
    assert not np.allclose(buggy, [0.5, 0.5, 0.0], atol=1e-12)


def test_corrected_use_j2_isolation():
    # sun_only must not include J2: static check of source guard
    src = Path(ME.__file__).read_text(encoding="utf-8")
    assert 'mode in ("j2_only", "sun_moon_j2")' in src or "mode in ('j2_only', 'sun_moon_j2')" in src


def test_S0_pins_prior_values():
    a = 6378.137 + 600.0
    s_sso = TH.secular_total_deg_day(a, math.radians(97.7876))
    s90 = TH.secular_total_deg_day(a, math.radians(90.0))
    s30 = TH.secular_total_deg_day(a, math.radians(30.0))
    assert abs(s_sso - 1.348e-4) / 1.348e-4 < 0.05
    assert abs(s90 - 1.739e-4) / 1.739e-4 < 0.05
    assert abs(s30 - 4.547e-5) / 4.547e-5 < 0.05


def test_ols_bias_formula_vs_synthetic():
    # synthetic: s=0.5 deg/day + A cos(omega t), W=365, compare exact g vs lstsq
    A, om, W, phi, s = 2e-4, 2 * math.pi / 365.0, 365.0, 0.7, 0.5
    t = np.linspace(0, W, 2001)
    y = s * t + (A / om) * np.sin(om * t + phi)
    A_mat = np.column_stack([np.ones_like(t), t])
    coef, *_ = np.linalg.lstsq(A_mat, y, rcond=None)
    measured = float(coef[1]) - s
    predicted = TH.ols_bias_deg_day(A, om, W, phi)
    assert abs(measured - predicted) / max(abs(A), 1e-12) < 0.01
    # resonant limit
    assert abs(TH.ols_bias_deg_day(A, 0.0, W, phi) - A * math.cos(phi)) < 1e-15
    # no-secular masquerade: pure sinusoid can still give nonzero OLS slope (documents the trap)
    t2 = np.linspace(0, 90.0, 901)
    y2 = (A / om) * np.sin(om * t2 + 0.0)
    c2, *_ = np.linalg.lstsq(np.column_stack([np.ones_like(t2), t2]), y2, rcond=None)
    assert abs(float(c2[1])) > 0  # apparent slope without any secular


def test_octupole_negligible_vs_expected_R():
    a = 6378.137 + 600.0
    for inc in (97.7876, 90.0, 30.0):
        b = TH.octupole_bound_deg_day(a, math.radians(inc))
        assert b < 0.1 * 2.4e-4  # <10% of smallest prior |R|


def test_averaged_propagator_recovers_prescribed_secular():
    h = [(1, 0, 0.0, 0.0)]
    out = AV.propagate_averaged(6978.0, math.radians(97.0), W_day=365.0, dt_day=1.0,
                                oj2_deg_day=0.9, s0_deg_day=1e-4, harmonics=h,
                                n3_deg_day=0.0, o3dot_deg_day=0.0)
    s_hat = AV.ols_slope_deg_day(out["t_day"], out["Omega_deg"])
    assert abs(s_hat - (0.9 + 1e-4)) / (0.9 + 1e-4) < 1e-3


def test_frozen_plane_residual_survives_without_closed_loop_cancellation():
    frozen = json.loads((HERE / "results" / "frozen_365d.json").read_text(encoding="utf-8"))["summary"]
    modulated = json.loads((HERE / "results" / "headline_365d.json").read_text(encoding="utf-8"))["summary"]
    for inc in ("97.7876", "90.0", "30.0"):
        frozen_r = abs(float(frozen[inc]["R"]))
        modulated_r = abs(float(modulated[inc]["R"]))
        assert frozen_r > 1e-10, f"frozen-plane R at {inc} deg is numerically zero"
        if float(inc) == 90.0:
            assert abs(frozen_r - modulated_r) < 1e-6, "90 deg has no J2 precession to freeze"
        else:
            assert frozen_r < modulated_r, f"frozen-plane R at {inc} deg is not reduced"


def test_prescribed_rate_response_removes_commanded_rotation():
    import run_prescribed as RP
    rows = [
        {"omega_p": -1.0, "rate": -0.9999},
        {"omega_p": 0.0, "rate": -0.0001},
        {"omega_p": 1.0, "rate": 1.0002},
    ]
    corrected = [RP.response_minus_command(r) for r in rows]
    np.testing.assert_allclose(corrected[0]["response"], 0.0001, atol=1e-15)
    np.testing.assert_allclose(corrected[1]["response"], -0.0001, atol=1e-15)
    np.testing.assert_allclose(corrected[2]["response"], 0.0002, atol=1e-15)


def test_active_prescribed_control_precesses_force_free_plane():
    sun = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[1.5e8, 0, 0]] * 2),
           "sha256": "x", "n_points": 2}
    moon = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[3.8e5, 0, 0]] * 2),
            "sha256": "x", "n_points": 2}
    a, inc = 6978.0, math.radians(97.7876)
    x0 = ME.circular_ic(a, inc, 0.0)
    # Active rigid rotation: a force-free plane must precess at +1 deg/day.
    out = ME.propagate(x0=x0, sun=sun, moon=moon, mode="kepler_only",
                       t0_s=0.0, t_end_s=10 * 86400.0, dt_s=60.0,
                       omega_p_deg_day=1.0, subsample_every=5000)
    assert abs(ME.rate_deg_day(out["t_node"], np.array(out["omega_node"])) - 1.0) < 0.05
    # No control: force-free node stays flat.
    out0 = ME.propagate(x0=x0, sun=sun, moon=moon, mode="kepler_only",
                        t0_s=0.0, t_end_s=10 * 86400.0, dt_s=60.0,
                        subsample_every=5000)
    assert abs(ME.rate_deg_day(out0["t_node"], np.array(out0["omega_node"]))) < 0.05


def test_removed_passive_frame_kwargs_fail_loudly():
    sun = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[1.5e8, 0, 0]] * 2),
           "sha256": "x", "n_points": 2}
    moon = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[3.8e5, 0, 0]] * 2),
            "sha256": "x", "n_points": 2}
    a, inc = 6978.0, math.radians(97.7876)
    x0 = ME.circular_ic(a, inc, 0.0)
    import pytest
    with pytest.raises(TypeError):
        ME.propagate(x0=x0, sun=sun, moon=moon, mode="kepler_only",
                     t0_s=0.0, t_end_s=86400.0,
                     measure_in_rotating_frame=True)


def _node_series(out):
    return out["t_node"], np.array(out["omega_node"])


def test_freeze_is_common_rotation_bias_invariant():
    # Open-loop guard: biasing the freeze counter-rate by delta must shift
    # both raw J2 rates by ~-delta while the paired R stays put.
    # Proves the control is feed-forward subtraction, not closed-loop kill.
    # NOTE sign: apply_control rotates by -freeze_rate, so +delta in the
    # commanded rate shifts the measured node by -delta.
    sun = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[1.5e8, 0, 0]] * 2),
           "sha256": "x", "n_points": 2}
    moon = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[3.8e5, 0, 0]] * 2),
            "sha256": "x", "n_points": 2}
    a, inc = 6978.0, math.radians(97.7876)
    x0 = ME.circular_ic(a, inc, 0.0)
    kw = dict(sun=sun, moon=moon, x0=x0, t0_s=0.0, t_end_s=10 * 86400.0,
              dt_s=60.0, subsample_every=5000)
    r_j2 = ME.rate_deg_day(*_node_series(ME.propagate(mode="j2_only", freeze_plane=True, **kw)))
    r_full = ME.rate_deg_day(*_node_series(ME.propagate(mode="sun_moon_j2", freeze_plane=True, **kw)))
    delta = 0.05
    base = ME.propagate(mode="j2_only", freeze_plane=True, **kw)["freeze_plane_rate_deg_day"]
    kw2 = dict(kw, freeze_plane_rate_deg_day=base + delta)
    r_j2_b = ME.rate_deg_day(*_node_series(ME.propagate(mode="j2_only", freeze_plane=True, **kw2)))
    r_full_b = ME.rate_deg_day(*_node_series(ME.propagate(mode="sun_moon_j2", freeze_plane=True, **kw2)))
    assert abs((r_j2_b - r_j2) + delta) < 0.05
    assert abs((r_full_b - r_full) + delta) < 0.05


def test_freeze_and_kinematic_controls_hold():
    sun = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[1.5e8, 0, 0]] * 2),
           "sha256": "x", "n_points": 2}
    moon = {"t_s": np.array([0.0, 1e9]), "r_eci_km": np.array([[3.8e5, 0, 0]] * 2),
            "sha256": "x", "n_points": 2}
    a, inc = 6978.0, math.radians(97.7876)
    x0 = ME.circular_ic(a, inc, 0.0)
    # kinematic rotation without forces must precess at prescribed rate (10 d)
    out = ME.propagate(x0=x0, sun=sun, moon=moon, mode="kepler_only", t0_s=0.0,
                       t_end_s=10 * 86400.0, dt_s=60.0, omega_p_deg_day=1.0,
                       subsample_every=5000)
    s_hat = ME.rate_deg_day(out["t_node"], np.array(out["omega_node"]))
    assert abs(s_hat - 1.0) < 0.05
