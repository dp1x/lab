"""Focused minimum-sufficient-horizon campaign runner.

Runs the headline 18.6-yr direct arc at i_sso (J2+Sun+Moon and J2-only
controls) with the full estimator hierarchy. The phase-locked 2-window
estimator is extracted as a sub-trajectory analysis from the 18.6-yr
propagation. Includes the idealized bridge and synthetic oracle.

Cost: ~150-2 minutes (two 18.6-yr propagations at dt=60s).

Contract frozen before run (per constitution preamble).
"""
from __future__ import annotations

import sys
import importlib.util
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LAB_ROOT = HERE.parents[1]
SUN_SNAPSHOT = HERE / "reference" / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MOON_SNAPSHOT = HERE / "reference" / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
MANIFEST_PATH = HERE / "reference" / "MANIFEST.json"

# Import experiment.py as a module
spec = importlib.util.spec_from_file_location("mission_lunisolar_closure",
                                                str(HERE / "experiment.py"))
exp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp)

from lab_utils.earth_frames import JD_J2000

OUT_DIR = HERE / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W_18P6_DAYS = 18.6 * 365.25
W_1YR_DAYS = 365.25
HALF_NODAL_DAYS = exp.LUNAR_NODAL_PERIOD_DAYS / 2.0

H_SSO_KM = 600.0
I_SSO_DEG = 97.7876


def sha256(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_headline_propagation(sun_snap, moon_snap, *, h_km: float, i_deg: float,
                              phase_d: float, mode: str, t0_s: float,
                              duration_days: float, dt_s: float = 60.0,
                              label: str = "") -> dict:
    """Run one propagation. Returns dict with t_cross, om_cross, t_node, omega_node,
    and the full trajectory for sub-window analysis."""
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
    n_steps = int(math.ceil((t_end - t0_s) / dt_s))
    t_grid = np.linspace(t0_s, t_end, n_steps + 1)
    print(f"[{label}] {n_steps} steps over {duration_days:.2f} d ... ", end="", flush=True)
    t0_prop = time.time()
    f = exp.make_rhs(sun_snap, moon_snap, mode=mode, apply_precession=True)
    x_traj = exp.rk4_propagate(f, t_grid, x0)
    elapsed = time.time() - t0_prop
    t_cross, om_cross = exp.detect_ascending_nodes(t_grid, x_traj)
    t_node, omega_node = exp.node_vector_series(t_grid, x_traj)
    print(f"done in {elapsed:.1f}s ({n_steps} steps, {len(t_cross)} nodes)")
    return {
        "label": label,
        "t_grid": t_grid,
        "x_traj": x_traj,
        "t_cross": t_cross,
        "om_cross": om_cross,
        "t_node": t_node,
        "omega_node": omega_node,
        "wall_clock_s": elapsed,
        "n_steps": n_steps,
        "n_nodes": int(len(t_cross)),
    }


def analyze_propagation(prop: dict, *, t0_s: float) -> dict:
    """Apply all estimators to a propagation."""
    t_cross = prop["t_cross"]
    om_cross = prop["om_cross"]
    t_node = prop["t_node"]
    omega_node = prop["omega_node"]
    out = {"n_nodes": int(len(t_cross))}
    if len(t_cross) >= 4:
        t_rel = (t_cross - t_cross[0]) / 86400.0
        # Direct OLS over the full arc
        _, b_a = exp.ols_slope(t_rel, om_cross)
        out["direct_ols_deg_per_day"] = math.degrees(b_a) * 86400.0
        # Secant
        out["secant_deg_per_day"] = math.degrees((om_cross[-1] - om_cross[0]) / t_rel[-1]) * 86400.0
        # Harmonic regression
        fit_f = exp.harmonic_regression(t_rel, om_cross)
        out["harmonic_regression_deg_per_day"] = fit_f["b_deg_per_day"]
        out["harmonic_regression_rms_deg"] = fit_f["rms_residual_deg"]
        out["harmonic_amplitudes_deg"] = fit_f["harmonic_amplitudes_deg"]
        # Node-vector OLS
        if len(t_node) > 10:
            _, b_n = exp.ols_slope((t_node - t_node[0]) / 86400.0, omega_node)
            out["node_vector_deg_per_day"] = math.degrees(b_n) * 86400.0
        # Phase-locked 2-window estimator at half-nodal separation
        pl = exp.phase_locked_two_window(
            prop["t_grid"], prop["x_traj"],
            window_days=W_1YR_DAYS,
            separation_days=HALF_NODAL_DAYS,
            t_start_s=t0_s,
        )
        out["phase_locked_avg_deg_per_day"] = pl.get("phase_locked_avg_slope_deg_day", float("nan"))
        out["phase_locked_window_a_deg_per_day"] = pl.get("window_a_slope_deg_day", float("nan"))
        out["phase_locked_window_b_deg_per_day"] = pl.get("window_b_slope_deg_day", float("nan"))
    return out


def main():
    t_start = time.time()
    print(f"=== Focused mission_lunisolar_closure campaign ===")

    # Load snapshots
    print(f"loading snapshots")
    sun_snap = exp._load_snapshot(SUN_SNAPSHOT)
    moon_snap = exp._load_snapshot(MOON_SNAPSHOT)
    print(f"  Sun: {sun_snap['n_points']} rows sha256={sun_snap['sha256'][:16]}")
    print(f"  Moon: {moon_snap['n_points']} rows sha256={moon_snap['sha256'][:16]}")

    # Pre-flight checks
    print(f"pre-flight: synthetic oracle")
    synth = exp.synthetic_oracle_test()
    print(f"  estimator (f) bias on synthetic oracle: {synth['estimator_f_bias_deg_day']:.3e}")
    print(f"  estimator (a) bias on synthetic oracle: {synth['estimator_a_bias_deg_day']:.3e}")
    print(f"  verdict: {synth['verdict']}")
    if abs(synth["estimator_f_bias_deg_day"]) > 1e-9:
        print(f"  WARNING: estimator (f) bias exceeds 1e-9 on synthetic oracle; investigate")

    print(f"pre-flight: force-level identity check")
    identity = exp.force_level_identity_check()
    print(f"  max_diff_sun = {identity['max_diff_sun_km_s2']:.3e} km/s^2")
    print(f"  max_diff_moon = {identity['max_diff_moon_km_s2']:.3e} km/s^2")
    if not (identity["passes_sun"] and identity["passes_moon"]):
        print(f"  FAIL: force-level identity does not pass")
        return

    # Corrected cf at canonical inclinations
    cf_sso = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, I_SSO_DEG)
    cf_90 = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, 90.0)
    cf_30 = exp.corrected_secular_lunisolar_raan_rate_rad_s(H_SSO_KM, 30.0)
    print(f"corrected cf at i_sso: solar={cf_sso['solar_deg_day']:+.4e}, "
          f"lunar={cf_sso['lunar_deg_day']:+.4e}, total={cf_sso['total_deg_day']:+.4e}")
    print(f"corrected cf at i=90:  solar={cf_90['solar_deg_day']:+.4e}, "
          f"lunar={cf_90['lunar_deg_day']:+.4e}, total={cf_90['total_deg_day']:+.4e}")
    print(f"corrected cf at i=30:  solar={cf_30['solar_deg_day']:+.4e}, "
          f"lunar={cf_30['lunar_deg_day']:+.4e}, total={cf_30['total_deg_day']:+.4e}")

    # Idealized bridge
    bridge_sso = exp.idealized_circular_perturber_bridge(H_SSO_KM, I_SSO_DEG)
    print(f"idealized bridge (i_sso): orbit-averaged={bridge_sso['idealized_orbit_averaged_nodal_deg_day']:+.4e}, "
          f"cf_lunar={bridge_sso['cf_lunar_component_deg_day']:+.4e}, "
          f"ratio={bridge_sso['ratio']:.3f}")

    # Epoch for t=0
    t0_s = (sun_snap["jd_start"] - JD_J2000) * 86400.0

    # ---- Main campaign ----
    # Phase 1: 18.6-yr J2-only at i_sso
    print(f"\nPhase 1: 18.6-yr J2-only at i_sso, phase 0")
    prop_j2 = run_headline_propagation(
        sun_snap, moon_snap,
        h_km=H_SSO_KM, i_deg=I_SSO_DEG, phase_d=0.0,
        mode="j2_only", t0_s=t0_s, duration_days=W_18P6_DAYS, dt_s=60.0,
        label="i_sso_J2_only_18p6yr_phase0"
    )
    j2_analysis = analyze_propagation(prop_j2, t0_s=t0_s)
    print(f"  J2-only direct OLS: {j2_analysis['direct_ols_deg_per_day']:+.6f} deg/day")
    print(f"  J2-only harmonic reg: {j2_analysis['harmonic_regression_deg_per_day']:+.6f} deg/day")
    print(f"  J2-only node-vector:  {j2_analysis.get('node_vector_deg_per_day', float('nan')):+.6f} deg/day")

    # Phase 2: 18.6-yr J2+Sun+Moon at i_sso (THE HEADLINE)
    print(f"\nPhase 2: 18.6-yr J2+Sun+Moon at i_sso, phase 0 (HEADLINE)")
    prop_full = run_headline_propagation(
        sun_snap, moon_snap,
        h_km=H_SSO_KM, i_deg=I_SSO_DEG, phase_d=0.0,
        mode="sun_moon_j2", t0_s=t0_s, duration_days=W_18P6_DAYS, dt_s=60.0,
        label="i_sso_J2+LS_18p6yr_phase0"
    )
    full_analysis = analyze_propagation(prop_full, t0_s=t0_s)
    print(f"  J2+LS direct OLS: {full_analysis['direct_ols_deg_per_day']:+.6f} deg/day")
    print(f"  J2+LS harmonic reg: {full_analysis['harmonic_regression_deg_per_day']:+.6f} deg/day")
    print(f"  J2+LS node-vector:  {full_analysis.get('node_vector_deg_per_day', float('nan')):+.6f} deg/day")
    print(f"  J2+LS phase-locked avg: {full_analysis.get('phase_locked_avg_deg_per_day', float('nan')):+.6f} deg/day")
    print(f"  J2+LS phase-locked A: {full_analysis.get('phase_locked_window_a_deg_per_day', float('nan')):+.6f} deg/day")
    print(f"  J2+LS phase-locked B: {full_analysis.get('phase_locked_window_b_deg_per_day', float('nan')):+.6f} deg/day")

    # Lunisolar contribution = full - j2 (for each estimator)
    ls_results = {}
    for est_name in ("direct_ols_deg_per_day", "harmonic_regression_deg_per_day",
                       "node_vector_deg_per_day", "phase_locked_avg_deg_per_day",
                       "phase_locked_window_a_deg_per_day",
                       "phase_locked_window_b_deg_per_day"):
        if est_name in full_analysis and est_name in j2_analysis:
            ls_results[est_name] = full_analysis[est_name] - j2_analysis[est_name]
    print(f"\nLunisolar contribution (full - j2_only):")
    for k, v in ls_results.items():
        print(f"  {k}: {v:+.4e} deg/day")

    # Phase 3: 18.6-yr J2+Sun+Moon at i=90 deg (inclination structure)
    print(f"\nPhase 3: 18.6-yr J2+Sun+Moon at i=90 deg (inclination structure)")
    prop_90_full = run_headline_propagation(
        sun_snap, moon_snap,
        h_km=H_SSO_KM, i_deg=90.0, phase_d=0.0,
        mode="sun_moon_j2", t0_s=t0_s, duration_days=W_18P6_DAYS, dt_s=60.0,
        label="i_90_J2+LS_18p6yr_phase0"
    )
    full90_analysis = analyze_propagation(prop_90_full, t0_s=t0_s)
    prop_90_j2 = run_headline_propagation(
        sun_snap, moon_snap,
        h_km=H_SSO_KM, i_deg=90.0, phase_d=0.0,
        mode="j2_only", t0_s=t0_s, duration_days=W_18P6_DAYS, dt_s=60.0,
        label="i_90_J2_only_18p6yr_phase0"
    )
    j2_90_analysis = analyze_propagation(prop_90_j2, t0_s=t0_s)
    ls_90 = {}
    for est_name in ("direct_ols_deg_per_day", "harmonic_regression_deg_per_day",
                       "node_vector_deg_per_day", "phase_locked_avg_deg_per_day"):
        if est_name in full90_analysis and est_name in j2_90_analysis:
            ls_90[est_name] = full90_analysis[est_name] - j2_90_analysis[est_name]
    print(f"\nLunisolar at i=90 deg (full - j2_only):")
    for k, v in ls_90.items():
        print(f"  {k}: {v:+.4e} deg/day")

    # Phase 4: 18.6-yr J2+Sun+Moon at i=30 deg (inclination structure, prograde opposite-sign)
    print(f"\nPhase 4: 18.6-yr J2+Sun+Moon at i=30 deg (inclination structure)")
    prop_30_full = run_headline_propagation(
        sun_snap, moon_snap,
        h_km=H_SSO_KM, i_deg=30.0, phase_d=0.0,
        mode="sun_moon_j2", t0_s=t0_s, duration_days=W_18P6_DAYS, dt_s=60.0,
        label="i_30_J2+LS_18p6yr_phase0"
    )
    full30_analysis = analyze_propagation(prop_30_full, t0_s=t0_s)
    prop_30_j2 = run_headline_propagation(
        sun_snap, moon_snap,
        h_km=H_SSO_KM, i_deg=30.0, phase_d=0.0,
        mode="j2_only", t0_s=t0_s, duration_days=W_18P6_DAYS, dt_s=60.0,
        label="i_30_J2_only_18p6yr_phase0"
    )
    j2_30_analysis = analyze_propagation(prop_30_j2, t0_s=t0_s)
    ls_30 = {}
    for est_name in ("direct_ols_deg_per_day", "harmonic_regression_deg_per_day",
                       "node_vector_deg_per_day", "phase_locked_avg_deg_per_day"):
        if est_name in full30_analysis and est_name in j2_30_analysis:
            ls_30[est_name] = full30_analysis[est_name] - j2_30_analysis[est_name]
    print(f"\nLunisolar at i=30 deg (full - j2_only):")
    for k, v in ls_30.items():
        print(f"  {k}: {v:+.4e} deg/day")

    # ---- Save results ----
    payload = {
        "meta": {
            "description": "Focused minimum-sufficient-horizon mission_lunisolar_closure campaign: 18.6-yr direct arc at h=600 km, 3 inclinations (i_sso=97.79, i=90, i=30), J2-only and J2+Sun+Moon modes, single phase (lunar anomalistic zero).",
            "git_commit": "PENDING",
            "name": "mission_lunisolar_closure_focused",
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
            "frame": "ECI mean-of-date; Sun and Moon rotated from ICRF/J2000 via FIXED IAU-1976 precession (Track D 019 remediation).",
            "units": "km, km^3/s^2, s since J2000 (TT-like); radians internal; degrees at I/O.",
            "horizons_arc_days": W_18P6_DAYS,
            "horizons_arc_years": 18.6,
            "phase_locked_separation_days": HALF_NODAL_DAYS,
            "inclinations_deg": {"i_sso": I_SSO_DEG, "i_90": 90.0, "i_30": 30.0},
            "phases": [0.0],
            "force_modes": ["j2_only", "sun_moon_j2"],
            "estimators": [
                "direct_OLS over full 18.6-yr arc",
                "secant over full arc",
                "harmonic_regression (theory-driven basis)",
                "node_vector (theory-INDEPENDENT kinematic observable)",
                "phase_locked_two_window (9.3-yr separation)",
            ],
            "headline_estimator": "harmonic_regression at 18.6-yr (Estimator f)",
            "decision_rule": "The corrected formula is VERIFIED if the 18.6-yr harmonic regression rate at i_sso agrees with the corrected formula within +/- 50% (i.e., 0.5x to 2.0x ratio). The phase-locked 2-window estimator should agree with the harmonic regression within +/- 30% if the secular rate is approximately constant over 9.3 yr.",
            "limitations": [
                "Real DE441 Sun/Moon ephemeris; the corrected formula assumes idealized circular perturbers (bridge experiment quantifies this).",
                "Fixed-step RK4 at dt=60s; convergence ladder verified at sub-1% accuracy for J2-only.",
                "Single phase per inclination (lunar anomalistic zero); 4-phase ensemble in original Exp 020 contract was replaced by full 18.6-yr direct fit which captures the secular rate over a full lunar nodal cycle.",
                "Lunar nodal modulation of the secular rate (via i_3 oscillation) means the 18.6-yr harmonic regression recovers the cycle-mean rate; the phase-locked 2-window estimator gives an epoch-specific estimate that should agree within ~30%.",
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
            "i_sso_J2_only": {k: v for k, v in j2_analysis.items() if not isinstance(v, np.ndarray)},
            "i_sso_J2_plus_LS": {k: v for k, v in full_analysis.items() if not isinstance(v, np.ndarray)},
            "i_sso_lunisolar_contribution": ls_results,
            "i_90_J2_only": {k: v for k, v in j2_90_analysis.items() if not isinstance(v, np.ndarray)},
            "i_90_J2_plus_LS": {k: v for k, v in full90_analysis.items() if not isinstance(v, np.ndarray)},
            "i_90_lunisolar_contribution": ls_90,
            "i_30_J2_only": {k: v for k, v in j2_30_analysis.items() if not isinstance(v, np.ndarray)},
            "i_30_J2_plus_LS": {k: v for k, v in full30_analysis.items() if not isinstance(v, np.ndarray)},
            "i_30_lunisolar_contribution": ls_30,
        },
        "comparison_with_corrected_formula": {
            "i_sso_cf_total_deg_day": cf_sso["total_deg_day"],
            "i_sso_harmonic_reg_lunisolar_deg_day": ls_results.get("harmonic_regression_deg_per_day", float("nan")),
            "i_sso_phase_locked_avg_lunisolar_deg_day": ls_results.get("phase_locked_avg_deg_per_day", float("nan")),
            "i_sso_ratio_harmonic_to_cf": (ls_results.get("harmonic_regression_deg_per_day", 0.0) /
                                              cf_sso["total_deg_day"] if cf_sso["total_deg_day"] != 0 else float("nan")),
            "i_sso_ratio_phase_locked_to_cf": (ls_results.get("phase_locked_avg_deg_per_day", 0.0) /
                                                  cf_sso["total_deg_day"] if cf_sso["total_deg_day"] != 0 else float("nan")),
            "i_90_cf_total_deg_day": cf_90["total_deg_day"],
            "i_90_harmonic_reg_lunisolar_deg_day": ls_90.get("harmonic_regression_deg_per_day", float("nan")),
            "i_30_cf_total_deg_day": cf_30["total_deg_day"],
            "i_30_harmonic_reg_lunisolar_deg_day": ls_30.get("harmonic_regression_deg_per_day", float("nan")),
        },
    }
    out_path = OUT_DIR / "results.json"
    exp.save_json_result(payload, out_path)
    print(f"\nResults saved to {out_path}")
    print(f"Total wall-clock: {time.time()-t_start:.1f}s "
          f"({(time.time()-t_start)/60:.1f} min)")


if __name__ == "__main__":
    main()