"""Focused tests for mission_mean_element_leakage (Phase 0 baseline breathing)."""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.modules.pop("phase0_baseline_breathing", None)
import phase0_baseline_breathing as P0  # noqa: E402


def test_j2_sensitivity_matches_finite_difference():
    """Analytic dOmegaJ2/di must match a central difference of the same function."""
    a = 6978.137
    for inc_deg in (97.7876, 90.0, 30.0, 0.0):
        i0 = math.radians(inc_deg)
        d = 1e-6  # rad
        num = (P0.omega_dot_j2_deg_day(a, i0 + d) - P0.omega_dot_j2_deg_day(a, i0 - d)) / (2 * d)
        ana = 1.5 * (math.sqrt(P0.MU_EARTH_KM3S2 / a ** 3) * 86400.0 * 180.0 / math.pi) \
            * P0.J2_EARTH * (P0.R_EARTH_KM / a) ** 2 * math.sin(i0)
        assert abs(num - ana) / max(abs(ana), 1e-12) < 1e-6


def test_sensitivity_constant_matches_h600_value():
    """The headline sensitivity used in the card: ~0.1258 deg/day per degree at h=600."""
    a = P0.R_EARTH_KM + P0.H_KM
    i0 = math.radians(97.7876)
    d = math.pi / 180.0
    val = (P0.omega_dot_j2_deg_day(a, i0 + d) - P0.omega_dot_j2_deg_day(a, i0)) / 1.0
    assert abs(val - 0.12579) / 0.12579 < 0.01


def test_baseline_breathing_formula_closes_on_synthetic_offset():
    """A pure mean-inclination offset di must produce D = sensitivity * di."""
    a = P0.R_EARTH_KM + P0.H_KM
    i_j = math.radians(97.7876)
    di_deg = 0.05
    i_f = math.radians(97.7876 + di_deg)
    D = P0.omega_dot_j2_deg_day(a, i_f) - P0.omega_dot_j2_deg_day(a, i_j)
    expected = 0.12579 * di_deg
    assert abs(D - expected) / expected < 0.02  # curvature over 0.05 deg is ~1%


def test_propagator_matches_campaign_j2_only_rate():
    """The Phase 0 integrator must reproduce the campaign's j2_only node rate."""
    sun = P0.ME.load_snapshot(P0.ME.SUN_SNAPSHOT)
    moon = P0.ME.load_snapshot(P0.ME.MOON_SNAPSHOT)
    a0 = P0.R_EARTH_KM + P0.H_KM
    inc = 97.7876
    days = 20.0
    x0 = P0.ME.circular_ic(a0, math.radians(inc), P0.T0_S)
    rec = P0.propagate_record(sun, moon, x0, "j2_only", P0.T0_S,
                              P0.T0_S + days * 86400.0)
    mine = P0.ME.rate_deg_day(rec["t"], rec["om"])
    out = P0.ME.propagate(sun, moon, x0, mode="j2_only", t0_s=P0.T0_S,
                          t_end_s=P0.T0_S + days * 86400.0, dt_s=60.0)
    theirs = P0.ME.rate_deg_day(np.array(out["t_cross"]), np.array(out["om_cross"]))
    assert abs(mine - theirs) < 1e-9, (mine, theirs)


def test_phase0_result_pins_campaign_residuals():
    """If the Phase 0 result file exists, its R must match the campaign headline."""
    path = HERE / "results" / "phase0_1yr.json"
    if not path.exists():
        return  # pre-campaign; the pin activates once the run lands
    data = json.loads(path.read_text(encoding="utf-8"))
    pinned = {"97.7876": 1.2959476885023864e-3,
              "90.0": 6.703347307517931e-4,
              "30.0": -2.961594415091425e-4}
    for inc, R_ref in pinned.items():
        R = data["summary"][inc]["R_measured"]
        assert abs(R - R_ref) < 0.05 * abs(R_ref), (inc, R, R_ref)
    assert "wall_s" not in data  # no timing telemetry in artifacts
    assert data["provenance"]["sun_sha256"].startswith("f2c4f048")
    assert data["provenance"]["moon_sha256"].startswith("aee85099")
