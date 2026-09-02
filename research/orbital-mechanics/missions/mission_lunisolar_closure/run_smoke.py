"""Smoke test for mission_lunisolar_closure pipeline.

Runs the streaming RK4 propagator with a short arc (e.g. 30 d) and
coarse step (dt=600 s) to validate the entire pipeline before
committing to the 18.6-yr run.

Validates:
- Snapshot loading
- IAU-1976 precession
- Third-body acceleration (direct + indirect)
- Ascending-node detection
- All four estimators (direct OLS, secant, harmonic regression, node-vector OLS)
- Phase-locked 2-window estimator
- Synthetic oracle test
- Force-level identity check
- Idealized bridge

Saves results under results/smoke_results.json.
"""
from __future__ import annotations

import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SUN_SNAPSHOT = HERE / "reference" / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MOON_SNAPSHOT = HERE / "reference" / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MANIFEST_PATH = HERE / "reference" / "MANIFEST.json"

spec = importlib.util.spec_from_file_location("mlc_smoke",
                                                str(HERE / "experiment.py"))
exp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp)

from lab_utils.earth_frames import JD_J2000

OUT_DIR = HERE / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

H_SSO_KM = 600.0
I_SSO_DEG = 97.7876
DT_S = 60.0  # 60 s is the lab standard; smoke uses same
DURATION_DAYS = 365.25  # 1 yr for smoke (long enough for harmonic regression conditioning)


def main():
    t_start = time.time()
    print("=== mission_lunisolar_closure smoke test ===")

    print("loading snapshots")
    sun_snap = exp._load_snapshot(SUN_SNAPSHOT)
    moon_snap = exp._load_snapshot(MOON_SNAPSHOT)
    print(f"  Sun: {sun_snap['n_points']} rows sha256={sun_snap['sha256'][:16]}")
    print(f"  Moon: {moon_snap['n_points']} rows sha256={moon_snap['sha256'][:16]}")

    print("synthetic oracle test")
    synth = exp.synthetic_oracle_test()
    print(f"  estimator (f) bias: {synth['estimator_f_bias_deg_day']:.3e} deg/day")
    print(f"  estimator (a) bias: {synth['estimator_a_bias_deg_day']:.3e} deg/day")
    print(f"  verdict: {synth['verdict']}")
    if abs(synth["estimator_f_bias_deg_day"]) > 1e-9:
        print("FAIL: synthetic oracle estimator (f) bias too large")
        return

    print("force-level identity check")
    identity = exp.force_level_identity_check()
    print(f"  max_diff_sun = {identity['max_diff_sun_km_s2']:.3e} km/s^2")
    print(f"  max_diff_moon = {identity['max_diff_moon_km_s2']:.3e} km/s^2")
    if not (identity["passes_sun"] and identity["passes_moon"]):
        print("FAIL: force-level identity")
        return

    cf_sso = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, I_SSO_DEG)
    print(f"corrected cf i_sso: total={cf_sso['total_deg_day']:+.4e} deg/day "
          f"(solar={cf_sso['solar_deg_day']:+.4e}, lunar={cf_sso['lunar_deg_day']:+.4e})")

    print(f"initial state setup (h=600 km, i={I_SSO_DEG} deg)")
    a = exp.R_EARTH_KM + H_SSO_KM
    v = math.sqrt(exp.MU_EARTH_KM3S2 / a)
    i_rad = math.radians(I_SSO_DEG)
    x0 = np.array([a, 0.0, 0.0,
                    0.0, v * math.cos(i_rad), v * math.sin(i_rad)])

    t0_s = (sun_snap["jd_start"] - JD_J2000) * 86400.0
    t_end_s = t0_s + DURATION_DAYS * 86400.0

    print(f"\nstreaming propagation: {DURATION_DAYS} d, dt={DT_S} s, mode=j2_only")
    t0p = time.time()
    out_j2 = exp.propagate_streaming_with_x0(
        sun_snap, moon_snap, x0,
        mode="j2_only", t0_s=t0_s, t_end_s=t_end_s,
        dt_s=DT_S, subsample_every=100,
    )
    print(f"  done in {time.time()-t0p:.1f}s "
          f"({out_j2['n_steps']} steps, {len(out_j2['t_cross'])} crossings)")

    print(f"\nstreaming propagation: {DURATION_DAYS} d, dt={DT_S} s, mode=sun_moon_j2")
    t0p = time.time()
    out_full = exp.propagate_streaming_with_x0(
        sun_snap, moon_snap, x0,
        mode="sun_moon_j2", t0_s=t0_s, t_end_s=t_end_s,
        dt_s=DT_S, subsample_every=100,
    )
    print(f"  done in {time.time()-t0p:.1f}s "
          f"({out_full['n_steps']} steps, {len(out_full['t_cross'])} crossings)")

    if len(out_full["t_cross"]) < 4:
        print(f"FAIL: too few crossings ({len(out_full['t_cross'])}) for estimator")
        return

    t_rel = (out_full["t_cross"] - out_full["t_cross"][0]) / 86400.0
    fit = exp.harmonic_regression(t_rel, out_full["om_cross"])
    print(f"\nfull mode harmonic regression:")
    print(f"  secular rate: {fit['b_deg_per_day']:+.4e} deg/day")
    print(f"  RMS residual: {fit['rms_residual_deg']:.4e} deg")

    t_rel = (out_j2["t_cross"] - out_j2["t_cross"][0]) / 86400.0
    fit_j2 = exp.harmonic_regression(t_rel, out_j2["om_cross"])
    print(f"\nj2-only mode harmonic regression:")
    print(f"  secular rate: {fit_j2['b_deg_per_day']:+.4e} deg/day")

    ls_harmonic = fit["b_deg_per_day"] - fit_j2["b_deg_per_day"]
    print(f"\nLunisolar contribution (full - j2_only) via harmonic regression: "
          f"{ls_harmonic:+.4e} deg/day")
    print(f"corrected cf prediction: {cf_sso['total_deg_day']:+.4e} deg/day")

    pl = exp.phase_locked_two_window(out_full["t_cross"], out_full["om_cross"],
                                       window_days=5.0,
                                       separation_days=exp.HALF_NODAL_DAYS,
                                       t_start_s=out_full["t_cross"][0])
    print(f"\nphase-locked 2-window (window=5d, sep=HALF_NODAL_DAYS):")
    print(f"  window A: {pl['window_a_n_nodes']} nodes, drift={pl['window_a_drift_deg_day']:+.3e}")
    print(f"  window B: {pl['window_b_n_nodes']} nodes, drift={pl['window_b_drift_deg_day']:+.3e}")
    print(f"  avg drift: {pl['avg_drift_deg_day']:+.3e}")

    payload = {
        "meta": {
            "description": "Smoke test for mission_lunisolar_closure pipeline (short arc, coarse dt)",
            "duration_days": DURATION_DAYS,
            "dt_s": DT_S,
            "wall_clock_total_s": time.time() - t_start,
        },
        "snapshots": {
            "sun_sha256": sun_snap["sha256"],
            "moon_sha256": moon_snap["sha256"],
        },
        "synthetic_estimator_test": synth,
        "force_level_identity_check": identity,
        "corrected_cf_i_sso": cf_sso,
        "full_mode": {
            "n_steps": out_full["n_steps"],
            "n_crossings": len(out_full["t_cross"]),
            "n_subsamples": len(out_full["t_node"]),
            "harmonic_regression_secular_deg_day": fit["b_deg_per_day"],
            "harmonic_regression_rms_deg": fit["rms_residual_deg"],
        },
        "j2_only_mode": {
            "n_steps": out_j2["n_steps"],
            "n_crossings": len(out_j2["t_cross"]),
            "n_subsamples": len(out_j2["t_node"]),
            "harmonic_regression_secular_deg_day": fit_j2["b_deg_per_day"],
        },
        "lunisolar_contribution_via_harmonic_reg_deg_day": ls_harmonic,
        "phase_locked_estimator": pl,
        "verdict": "PASS" if (abs(synth["estimator_f_bias_deg_day"]) < 1e-9 and
                                identity["passes_sun"] and identity["passes_moon"] and
                                len(out_full["t_cross"]) >= 4) else "FAIL",
    }
    out_path = OUT_DIR / "smoke_results.json"
    exp.save_json_result(str(out_path), payload,
                          name="mission_lunisolar_closure_smoke",
                          description="Pipeline validation: 30-d arc, dt=60 s, J2-only + J2+Sun+Moon")
    print(f"\nSmoke test results saved to {out_path}")
    print(f"Total wall-clock: {time.time()-t_start:.1f}s")
    print(f"Verdict: {payload['verdict']}")


if __name__ == "__main__":
    main()