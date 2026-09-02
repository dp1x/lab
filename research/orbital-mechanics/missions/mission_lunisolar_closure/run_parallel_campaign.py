"""Parallel campaign runner for mission_lunisolar_closure.

Runs the 6 propagations (3 inclinations x 2 modes) in parallel across
all available CPU cores via multiprocessing.Pool. Each propagation is
independent (different mode + initial state) so parallelism is safe.

Designed for an 8-core CPU (Snapdragon X Elite); uses min(8, cpu_count()-1)
workers to leave 1 core for OS and other tasks.

For each propagation:
- Streams RK4 step-by-step (no full-trajectory storage)
- Collects ascending-node crossings (~100k over 18.6 yr)
- Subsamples node-vector samples every 100 steps (~98k samples)
- Total per propagation: ~9.8M RK4 steps, ~1.5 GB peak memory

After all propagations complete:
- Computes Lunisolar contribution as (full - j2_only) per inclination
- Applies 4 estimators: direct OLS, secant, harmonic regression, node-vector OLS
- Saves results.json
"""
from __future__ import annotations

import importlib.util
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUN_SNAPSHOT = HERE / "reference" / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MOON_SNAPSHOT = HERE / "reference" / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MANIFEST_PATH = HERE / "reference" / "MANIFEST.json"

OUT_DIR = HERE / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W_18P6_DAYS = 18.6 * 365.25
DT_S = 60.0
SUBSAMPLE_EVERY = 100
H_SSO_KM = 600.0
I_SSO_DEG = 97.7876

INCLINATIONS = {
    "i_sso": I_SSO_DEG,
    "i_90": 90.0,
    "i_30": 30.0,
}


# --------------------------------------------------------------------------- #
# Worker process: load module + run one propagation
# --------------------------------------------------------------------------- #
_WORKER_EXP = None
_WORKER_SUN = None
_WORKER_MOON = None
_WORKER_T0 = None


def _init_worker(exp_path, sun_snap_path, moon_snap_path, t0_s):
    """Per-process initializer: load module and snapshots once."""
    global _WORKER_EXP, _WORKER_SUN, _WORKER_MOON, _WORKER_T0
    # numpy is needed in worker for np.array in _run_one_propagation
    global np
    import numpy as np
    spec = importlib.util.spec_from_file_location("mlc_worker", str(exp_path))
    exp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exp)
    _WORKER_EXP = exp
    _WORKER_SUN = exp._load_snapshot(sun_snap_path)
    _WORKER_MOON = exp._load_snapshot(moon_snap_path)
    _WORKER_T0 = t0_s


def _run_one_propagation(args):
    """Run one (inclination, mode) propagation. Returns dict with
    the inputs + the streaming output + the analysis."""
    incl_name, i_deg, mode, label = args
    exp = _WORKER_EXP
    a = exp.R_EARTH_KM + H_SSO_KM
    v = math.sqrt(exp.MU_EARTH_KM3S2 / a)
    i_rad = math.radians(i_deg)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, v * math.cos(i_rad), v * math.sin(i_rad)])
    x0 = np.concatenate([r0, v0])
    t_end = _WORKER_T0 + W_18P6_DAYS * 86400.0
    t0p = time.time()
    out = exp.propagate_streaming_with_x0(
        _WORKER_SUN, _WORKER_MOON, x0,
        mode=mode, t0_s=_WORKER_T0, t_end_s=t_end,
        dt_s=DT_S, subsample_every=SUBSAMPLE_EVERY,
    )
    elapsed = time.time() - t0p

    # Analysis
    ana = {
        "n_nodes": int(len(out["t_cross"])),
        "n_subsamples": int(len(out["t_node"])),
        "wall_clock_s": float(elapsed),
    }
    if len(out["t_cross"]) >= 4:
        t_rel = (out["t_cross"] - out["t_cross"][0]) / 86400.0
        _, b_a = exp.ols_slope(t_rel, out["om_cross"])
        ana["direct_ols_deg_per_day"] = float(math.degrees(b_a))
        ana["secant_deg_per_day"] = float(math.degrees(
            (out["om_cross"][-1] - out["om_cross"][0]) / t_rel[-1]))
        fit_f = exp.harmonic_regression(t_rel, out["om_cross"])
        ana["harmonic_regression_deg_per_day"] = float(fit_f["b_deg_per_day"])
        ana["harmonic_regression_rms_residual_deg"] = float(fit_f["rms_residual_deg"])
        ana["harmonic_amplitudes_deg"] = fit_f["harmonic_amplitudes_deg"]
    if len(out["t_node"]) >= 10:
        t_rel = (out["t_node"] - out["t_node"][0]) / 86400.0
        _, b_n = exp.ols_slope(t_rel, out["omega_node"])
        ana["node_vector_deg_per_day"] = float(math.degrees(b_n))
    return {
        "label": label,
        "incl_name": incl_name,
        "i_deg": i_deg,
        "mode": mode,
        "analysis": ana,
    }


def main():
    t_start = time.time()

    # Determine worker count: cpu_count - 1 (leave 1 core for OS + lead)
    cpu_count = mp.cpu_count()
    n_workers = max(1, cpu_count - 1)
    print(f"=== mission_lunisolar_closure PARALLEL campaign ===")
    print(f"CPU cores: {cpu_count}; workers: {n_workers}")

    print("loading snapshots (lead process, for t0)")
    spec = importlib.util.spec_from_file_location("mlc_lead", str(HERE / "experiment.py"))
    exp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exp)
    from lab_utils.earth_frames import JD_J2000

    sun_snap = exp._load_snapshot(SUN_SNAPSHOT)
    moon_snap = exp._load_snapshot(MOON_SNAPSHOT)
    print(f"  Sun: {sun_snap['n_points']} rows sha256={sun_snap['sha256'][:16]}")
    print(f"  Moon: {moon_snap['n_points']} rows sha256={moon_snap['sha256'][:16]}")

    print("pre-flight: synthetic oracle")
    synth = exp.synthetic_oracle_test()
    print(f"  estimator (f) bias: {synth['estimator_f_bias_deg_day']:.3e} deg/day")
    print(f"  estimator (a) bias: {synth['estimator_a_bias_deg_day']:.3e} deg/day")
    print(f"  verdict: {synth['verdict']}")

    print("pre-flight: force-level identity check")
    identity = exp.force_level_identity_check()
    print(f"  max_diff_sun = {identity['max_diff_sun_km_s2']:.3e} km/s^2")
    print(f"  max_diff_moon = {identity['max_diff_moon_km_s2']:.3e} km/s^2")
    if not (identity["passes_sun"] and identity["passes_moon"]):
        print("FAIL")
        return

    cf_sso = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, I_SSO_DEG)
    cf_90 = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, 90.0)
    cf_30 = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, 30.0)
    print(f"corrected cf: i_sso={cf_sso['total_deg_day']:+.4e}, "
          f"i=90={cf_90['total_deg_day']:+.4e}, i=30={cf_30['total_deg_day']:+.4e}")

    bridge_sso = exp.idealized_circular_perturber_bridge(H_SSO_KM, I_SSO_DEG)
    print(f"idealized bridge (i_sso): ratio={bridge_sso['ratio']:.3f}")

    t0_s = (sun_snap["jd_start"] - JD_J2000) * 86400.0

    # Build the 6 tasks
    tasks = []
    for incl_name, incl_deg in INCLINATIONS.items():
        for mode in ("j2_only", "sun_moon_j2"):
            label = f"{incl_name}_{mode}_18p6yr"
            tasks.append((incl_name, incl_deg, mode, label))

    print(f"\nrunning {len(tasks)} propagations in parallel "
          f"(each {W_18P6_DAYS:.1f} d, dt={DT_S} s)")

    # The numpy import is needed in the worker for np.array in _run_one_propagation
    global np
    import numpy as np

    # Pool of workers; each loads its own snapshots
    ctx = mp.get_context("spawn")  # explicit start method for cross-platform
    with ctx.Pool(processes=n_workers,
                    initializer=_init_worker,
                    initargs=(HERE / "experiment.py",
                                SUN_SNAPSHOT,
                                MOON_SNAPSHOT,
                                t0_s)) as pool:
        results = pool.map(_run_one_propagation, tasks)

    # Aggregate
    all_analyses = {r["label"]: r["analysis"] for r in results}
    all_lunisolar = {incl_name: {} for incl_name in INCLINATIONS.keys()}
    for incl_name in INCLINATIONS.keys():
        full_key = f"{incl_name}_sun_moon_j2_18p6yr"
        j2_key = f"{incl_name}_j2_only_18p6yr"
        if full_key in all_analyses and j2_key in all_analyses:
            full = all_analyses[full_key]
            j2 = all_analyses[j2_key]
            for k in ("direct_ols_deg_per_day", "secant_deg_per_day",
                       "harmonic_regression_deg_per_day",
                       "node_vector_deg_per_day"):
                if k in full and k in j2:
                    all_lunisolar[incl_name][k] = full[k] - j2[k]

    print(f"\n=== Lunisolar contributions (full - j2_only) ===")
    for incl_name, ls in all_lunisolar.items():
        print(f"\n{incl_name}:")
        for k, v in ls.items():
            print(f"  {k}: {v:+.6e} deg/day")

    payload = {
        "meta": {
            "description": "mission_lunisolar_closure parallel campaign: 18.6-yr direct arc, "
                            "3 inclinations x 2 modes, single phase. "
                            "STREAMING RK4 (no full-trajectory storage). "
                            "Parallelized across all CPU cores via multiprocessing.Pool.",
            "name": "mission_lunisolar_closure_parallel",
            "n_workers": n_workers,
            "wall_clock_total_s": time.time() - t_start,
        },
        "code_sha256": exp.code_hashes(),
        "constants": {
            "R_E_km": exp.R_EARTH_KM,
            "J2": exp.J2_EARTH,
            "mu_E_km3_s2": exp.MU_EARTH_KM3S2,
            "mu_Sun_km3_s2": exp.SOLAR_GM_KM3_S2,
            "mu_Moon_km3_s2": exp.LUNAR_GM_KM3_S2,
            "AU_km": exp.AU_KM,
            "LUNAR_DISTANCE_KM_cf": exp.LUNAR_DISTANCE_KM_MEAN,
            "LUNAR_INCLINATION_DEG": exp.LUNAR_INCLINATION_DEG,
            "SOLAR_OBLIQUITY_DEG": exp.SOLAR_OBLIQUITY_DEG,
            "LUNAR_NODAL_PERIOD_DAYS": exp.LUNAR_NODAL_PERIOD_DAYS,
            "HALF_NODAL_PERIOD_DAYS": exp.HALF_NODAL_DAYS,
        },
        "contract": {
            "frame": "ECI mean-of-date; Sun and Moon rotated from ICRF/J2000 via FIXED IAU-1976 precession.",
            "units": "km, km^3/s^2, s since J2000 (TT-like); radians internal; degrees at I/O.",
            "horizons_arc_days": W_18P6_DAYS,
            "horizons_arc_years": 18.6,
            "phase": 0.0,
            "inclinations_deg": {"i_sso": I_SSO_DEG, "i_90": 90.0, "i_30": 30.0},
            "force_modes": ["j2_only", "sun_moon_j2"],
            "estimators": [
                "direct_OLS on ascending-node crossings",
                "secant on ascending-node crossings",
                "harmonic_regression (theory-driven OLS; Estimator f)",
                "node_vector (theory-INDEPENDENT kinematic observable; Estimator n)",
            ],
            "headline_estimator": "harmonic_regression at 18.6-yr (Estimator f)",
            "decision_rule": "Corrected formula is VERIFIED-WITH-LIMITATION if the 18.6-yr harmonic regression rate at i_sso agrees with the corrected formula within +/- 50%; the i=90 and i=30 controls must agree within +/- 100%.",
            "limitations": [
                "Single phase per inclination (lunar anomalistic zero); 18.6-yr direct fit over a full lunar nodal cycle averages over the nodal modulation of the secular rate.",
                "No J2 x Lunisolar coupling term in the corrected formula.",
                "No atmospheric drag.",
                "Point-mass Sun; no solar radiation pressure.",
            ],
        },
        "snapshots": {
            "sun_sha256": sun_snap["sha256"],
            "sun_n_points": sun_snap["n_points"],
            "sun_duration_days": sun_snap["duration_days"],
            "moon_sha256": moon_snap["sha256"],
            "moon_n_points": moon_snap["n_points"],
            "moon_duration_days": moon_snap["duration_days"],
            "manifest_sha256": exp._sha256(MANIFEST_PATH),
        },
        "synthetic_estimator_test": synth,
        "force_level_identity_check": identity,
        "corrected_closed_form": {
            "i_sso": cf_sso,
            "i_90": cf_90,
            "i_30": cf_30,
        },
        "idealized_circular_perturber_bridge_i_sso": bridge_sso,
        "headline_propagations": {
            incl_name: {
                f"{mode}_analysis": all_analyses.get(f"{incl_name}_{mode}_18p6yr", {})
                for mode in ("j2_only", "sun_moon_j2")
            }
            for incl_name in INCLINATIONS.keys()
        },
        "lunisolar_estimates": all_lunisolar,
        "comparison_with_corrected_formula": {
            "i_sso_cf_total_deg_day": cf_sso["total_deg_day"],
            "i_sso_harmonic_reg_lunisolar_deg_day": all_lunisolar["i_sso"].get("harmonic_regression_deg_per_day", float("nan")),
            "i_sso_ratio_harmonic_to_cf": (all_lunisolar["i_sso"].get("harmonic_regression_deg_per_day", 0.0) /
                                              cf_sso["total_deg_day"] if cf_sso["total_deg_day"] != 0 else float("nan")),
            "i_90_cf_total_deg_day": cf_90["total_deg_day"],
            "i_90_harmonic_reg_lunisolar_deg_day": all_lunisolar["i_90"].get("harmonic_regression_deg_per_day", float("nan")),
            "i_90_ratio_harmonic_to_cf": (all_lunisolar["i_90"].get("harmonic_regression_deg_per_day", 0.0) /
                                              cf_90["total_deg_day"] if cf_90["total_deg_day"] != 0 else float("nan")),
            "i_30_cf_total_deg_day": cf_30["total_deg_day"],
            "i_30_harmonic_reg_lunisolar_deg_day": all_lunisolar["i_30"].get("harmonic_regression_deg_per_day", float("nan")),
            "i_30_ratio_harmonic_to_cf": (all_lunisolar["i_30"].get("harmonic_regression_deg_per_day", 0.0) /
                                              cf_30["total_deg_day"] if cf_30["total_deg_day"] != 0 else float("nan")),
        },
    }
    out_path = OUT_DIR / "results.json"
    exp.save_json_result(str(out_path), payload,
                          name="mission_lunisolar_closure_parallel",
                          description="18.6-yr direct arc, 3 inclinations x 2 modes, parallelized")
    print(f"\nResults saved to {out_path}")
    print(f"Total wall-clock: {time.time()-t_start:.1f}s ({(time.time()-t_start)/60:.1f} min)")


if __name__ == "__main__":
    main()