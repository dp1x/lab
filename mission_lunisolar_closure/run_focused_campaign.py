"""Focused campaign runner using streaming propagation (memory-efficient).

Runs 6 propagations: 3 inclinations (i_sso, i=90, i=30) x 2 modes
(j2_only, sun_moon_j2). Each propagation is ~15 min on commodity
hardware (J2-only) or ~25-30 min (J2+Sun+Moon) because the J2-only
propagation has fewer RHS evaluations per step.

Each propagation:
- Streams RK4 step-by-step (no full-trajectory storage)
- Collects ascending-node crossings (n_crossings ~103k over 18.6 yr)
- Subsamples node-vector samples every 100 steps (~98k samples)
- Total per propagation: ~9.8M RK4 steps, ~1.5 GB peak memory

Applies 3 estimators: direct OLS, theory-driven harmonic regression,
theory-INDEPENDENT node-vector OLS. Computes Lunisolar contribution
as (full - j2_only) per inclination.
"""
from __future__ import annotations

import importlib.util
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LAB_ROOT = HERE.parents[1]
SUN_SNAPSHOT = HERE / "reference" / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MOON_SNAPSHOT = HERE / "reference" / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MANIFEST_PATH = HERE / "reference" / "MANIFEST.json"

spec = importlib.util.spec_from_file_location("mission_lunisolar_closure",
                                                str(HERE / "experiment.py"))
exp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp)

from lab_utils.earth_frames import JD_J2000

OUT_DIR = HERE / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W_18P6_DAYS = 18.6 * 365.25
W_1YR_DAYS = 365.25
HALF_NODAL_DAYS = exp.HALF_NODAL_DAYS

H_SSO_KM = 600.0
I_SSO_DEG = 97.7876

INCLINATIONS = {
    "i_sso": I_SSO_DEG,
    "i_90": 90.0,
    "i_30": 30.0,
}


def sha256(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_streaming(sun_snap, moon_snap, *, h_km: float, i_deg: float,
                     phase_d: float, mode: str, t0_s: float,
                     duration_days: float, dt_s: float = 60.0,
                     subsample_every: int = 100,
                     label: str = "") -> dict:
    a = exp.R_EARTH_KM + h_km
    v = math.sqrt(exp.MU_EARTH_KM3S2 / a)
    i_rad = math.radians(i_deg)
    r0_base = np.array([a, 0.0, 0.0])
    v0_base = np.array([0.0, v * math.cos(i_rad), v * math.sin(i_rad)])
    phase_rad = 2 * math.pi * phase_d / 27.5546
    cos_p, sin_p = math.cos(phase_rad), math.sin(phase_rad)
    r0 = np.array([
        r0_base[0] * cos_p - r0_base[1] * sin_p,
        r0_base[0] * sin_p + r0_base[1] * cos_p,
        r0_base[2],
    ])
    v0 = np.array([
        v0_base[0] * cos_p - v0_base[1] * sin_p,
        v0_base[0] * sin_p + v0_base[1] * cos_p,
        v0_base[2],
    ])
    x0 = np.concatenate([r0, v0])
    t_end = t0_s + duration_days * 86400.0
    print(f"[{label}] {duration_days:.1f} d, dt={dt_s}s, subsample={subsample_every} ... ", end="", flush=True)
    t0p = time.time()
    out = exp.propagate_streaming_with_x0(
        sun_snap, moon_snap, x0,
        mode=mode, t0_s=t0_s, t_end_s=t_end,
        dt_s=dt_s, subsample_every=subsample_every,
    )
    elapsed = time.time() - t0p
    print(f"done in {elapsed:.0f}s ({out['n_steps']} steps, "
          f"{len(out['t_cross'])} crossings, {len(out['t_node'])} nodes)")
    return out


def analyze_streaming(prop: dict) -> dict:
    """Apply direct OLS, harmonic regression, and node-vector OLS to the
    ascending-node-crossing and node-vector samples."""
    out = {
        "n_nodes": int(len(prop["t_cross"])),
        "n_subsamples": int(len(prop["t_node"])),
    }
    # Direct OLS on ascending-node crossings
    if len(prop["t_cross"]) >= 4:
        t_rel = (prop["t_cross"] - prop["t_cross"][0]) / 86400.0  # days
        _, b_a = exp.ols_slope(t_rel, prop["om_cross"])  # rad/day
        out["direct_ols_deg_per_day"] = math.degrees(b_a)
        # Secant
        out["secant_deg_per_day"] = math.degrees((prop["om_cross"][-1] - prop["om_cross"][0]) / t_rel[-1])
        # Harmonic regression
        fit_f = exp.harmonic_regression(t_rel, prop["om_cross"])
        out["harmonic_regression_deg_per_day"] = fit_f["b_deg_per_day"]
        out["harmonic_regression_rms_residual_deg"] = fit_f["rms_residual_deg"]
        out["harmonic_amplitudes_deg"] = fit_f["harmonic_amplitudes_deg"]
    # Node-vector OLS on subsampled samples
    if len(prop["t_node"]) >= 10:
        t_rel = (prop["t_node"] - prop["t_node"][0]) / 86400.0
        _, b_n = exp.ols_slope(t_rel, prop["omega_node"])
        out["node_vector_deg_per_day"] = math.degrees(b_n)
    return out


def main():
    t_start = time.time()
    print(f"=== Focused mission_lunisolar_closure campaign (STREAMING) ===")

    print(f"loading snapshots")
    sun_snap = exp._load_snapshot(SUN_SNAPSHOT)
    moon_snap = exp._load_snapshot(MOON_SNAPSHOT)
    print(f"  Sun: {sun_snap['n_points']} rows sha256={sun_snap['sha256'][:16]}")
    print(f"  Moon: {moon_snap['n_points']} rows sha256={moon_snap['sha256'][:16]}")

    print(f"pre-flight: synthetic oracle")
    synth = exp.synthetic_oracle_test()
    print(f"  estimator (f) bias: {synth['estimator_f_bias_deg_day']:.3e} deg/day")
    print(f"  estimator (a) bias: {synth['estimator_a_bias_deg_day']:.3e} deg/day")
    print(f"  verdict: {synth['verdict']}")

    print(f"pre-flight: force-level identity check")
    identity = exp.force_level_identity_check()
    print(f"  max_diff_sun = {identity['max_diff_sun_km_s2']:.3e} km/s^2")
    print(f"  max_diff_moon = {identity['max_diff_moon_km_s2']:.3e} km/s^2")
    if not (identity["passes_sun"] and identity["passes_moon"]):
        print(f"  FAIL")
        return

    cf_sso = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, I_SSO_DEG)
    cf_90 = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, 90.0)
    cf_30 = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, 30.0)
    print(f"corrected cf: i_sso={cf_sso['total_deg_day']:+.4e}, "
          f"i=90={cf_90['total_deg_day']:+.4e}, i=30={cf_30['total_deg_day']:+.4e}")

    bridge_sso = exp.idealized_circular_perturber_bridge(H_SSO_KM, I_SSO_DEG)
    print(f"idealized bridge (i_sso): ratio={bridge_sso['ratio']:.3f}")

    t0_s = (sun_snap["jd_start"] - JD_J2000) * 86400.0

    # Main campaign: 3 inclinations x 2 modes = 6 propagations
    all_propagations = {}
    all_analyses = {}
    all_lunisolar = {}

    for incl_name, incl_deg in INCLINATIONS.items():
        all_lunisolar[incl_name] = {}
        for mode in ("j2_only", "sun_moon_j2"):
            label = f"{incl_name}_{mode}_18p6yr"
            print(f"\n{label}:")
            prop = run_streaming(
                sun_snap, moon_snap,
                h_km=H_SSO_KM, i_deg=incl_deg, phase_d=0.0,
                mode=mode, t0_s=t0_s, duration_days=W_18P6_DAYS, dt_s=60.0,
                subsample_every=100, label=label,
            )
            all_propagations[label] = prop
            ana = analyze_streaming(prop)
            all_analyses[label] = ana
            for k, v in ana.items():
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    print(f"  {k}: {v:+.6e}" if isinstance(v, float) else f"  {k}: {v}")

    # Compute Lunisolar contribution per inclination
    for incl_name in INCLINATIONS.keys():
        full_key = f"{incl_name}_sun_moon_j2_18p6yr"
        j2_key = f"{incl_name}_j2_only_18p6yr"
        if full_key in all_analyses and j2_key in all_analyses:
            full = all_analyses[full_key]
            j2 = all_analyses[j2_key]
            ls = {}
            for k in ("direct_ols_deg_per_day", "secant_deg_per_day",
                       "harmonic_regression_deg_per_day", "node_vector_deg_per_day"):
                if k in full and k in j2:
                    ls[k] = full[k] - j2[k]
            all_lunisolar[incl_name] = ls
        # Also compute phase-locked estimator (would need windowed extraction;
        # omitted in streaming version since the trajectory isn't stored)

    print(f"\n=== Lunisolar contributions (full - j2_only) ===")
    for incl_name in INCLINATIONS.keys():
        print(f"\n{incl_name}:")
        for k, v in all_lunisolar[incl_name].items():
            print(f"  {k}: {v:+.6e} deg/day")

    payload = {
        "meta": {
            "description": "Mission 1 Lunisolar Capability Closure: 18.6-yr direct arc, "
                            "3 inclinations, J2-only and J2+Sun+Moon modes, single phase. "
                            "STREAMING RK4 (no full-trajectory storage).",
            "git_commit": "PENDING",
            "name": "mission_lunisolar_closure_focused_streaming",
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
                "Phase-locked 2-window estimator not implemented in the streaming version; the 18.6-yr harmonic regression with the lunar-nodal period in the basis subsumes its role.",
            ],
        },
        "snapshots": {
            "sun_sha256": sun_snap["sha256"],
            "sun_n_points": sun_snap["n_points"],
            "sun_duration_days": sun_snap["duration_days"],
            "moon_sha256": moon_snap["sha256"],
            "moon_n_points": moon_snap["n_points"],
            "moon_duration_days": moon_snap["duration_days"],
            "manifest_sha256": sha256(MANIFEST_PATH),
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
                f"{mode}_analysis": {k: v for k, v in ana.items() if not isinstance(v, np.ndarray)}
                for mode, ana in [
                    ("j2_only", all_analyses.get(f"{incl_name}_j2_only_18p6yr", {})),
                    ("sun_moon_j2", all_analyses.get(f"{incl_name}_sun_moon_j2_18p6yr", {})),
                ]
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
            "i_30_cf_total_deg_day": cf_30["total_deg_day"],
            "i_30_harmonic_reg_lunisolar_deg_day": all_lunisolar["i_30"].get("harmonic_regression_deg_per_day", float("nan")),
        },
    }
    out_path = OUT_DIR / "results.json"
    exp.save_json_result(payload, out_path)
    print(f"\nResults saved to {out_path}")
    print(f"Total wall-clock: {time.time()-t_start:.1f}s ({(time.time()-t_start)/60:.1f} min)")


if __name__ == "__main__":
    main()