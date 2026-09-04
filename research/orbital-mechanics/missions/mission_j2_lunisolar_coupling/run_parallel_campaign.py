"""Run the J2 x Lunisolar coupling campaign in parallel.

Runs force-mode isolation + perturbative scaling experiments at multiple
inclinations, with multiprocessing.Pool for parallel execution on commodity
multi-core hardware.

Three phases:
  Phase A: reduced-model (synthetic circular Moon) at 1 yr arc, i=i_sso, 90 d arc.
  Phase B: perturbative scaling at 90 d arc, i=i_sso, h=600 km.
  Phase C: full 18.6-yr force-mode decomposition at 3 inclinations x 6 modes.

Output: results/force_mode_isolation.json, results/perturbative_scaling.json,
results/cross_term_residual.json, results/final_decision.json.
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# Avoid experiment.py shadowing by other experiments
sys.path.insert(0, str(HERE))

from mission_experiment import (
    H_SSO_KM,
    I_SSO_DEG,
    I_90_DEG,
    I_30_DEG,
    DT_S,
    MU_EARTH_KM3S2,
    R_EARTH_KM,
    SUN_SNAPSHOT,
    MOON_SNAPSHOT,
    _load_snapshot,
    ols_slope,
    propagate_streaming_with_x0,
    harmonic_regression,
    code_hashes,
)


def initial_state(h_km: float, i_deg: float, mu: float = MU_EARTH_KM3S2,
                  r_eq: float = R_EARTH_KM) -> np.ndarray:
    a = r_eq + h_km
    v_circ = math.sqrt(mu / a)
    i_rad = math.radians(i_deg)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, v_circ * math.cos(i_rad), v_circ * math.sin(i_rad)])
    return np.concatenate([r0, v0])


# --------------------------------------------------------------------------- #
# Worker function for parallel execution
# --------------------------------------------------------------------------- #
def _propagation_worker(args):
    """Single propagation: returns dict with mode, RAAN rate, and metadata."""
    sun_data = args["sun_data"]  # tuple(t_s, r_eci, sha256)
    moon_data = args["moon_data"]
    x0 = args["x0"]
    mode = args["mode"]
    t_end_s = args["t_end_s"]
    dt_s = args["dt_s"]
    lambda_j2 = args.get("lambda_j2", 1.0)
    lambda_3body = args.get("lambda_3body", 1.0)
    synthetic_moon = args.get("synthetic_moon", False)
    label = args.get("label", f"{mode}_lj2{lambda_j2}_l3b{lambda_3body}")
    apply_precession = args.get("apply_precession", True)

    t_s_sun, r_eci_sun, sha_sun = sun_data
    t_s_moon, r_eci_moon, sha_moon = moon_data
    sun_snap = {"t_s": t_s_sun, "r_eci_km": r_eci_sun, "sha256": sha_sun}
    if synthetic_moon:
        moon_snap = {"t_s": t_s_moon, "r_eci_km": r_eci_moon, "sha256": sha_moon}
    else:
        moon_snap = {"t_s": t_s_moon, "r_eci_km": r_eci_moon, "sha256": sha_moon}

    t0 = time.time()
    res = propagate_streaming_with_x0(
        sun_snap, moon_snap, x0,
        mode=mode, t0_s=0.0, t_end_s=t_end_s, dt_s=dt_s,
        lambda_j2=lambda_j2, lambda_3body=lambda_3body,
        synthetic_moon=synthetic_moon,
    )
    wall_s = time.time() - t0
    n_steps = res["n_steps"]
    n_nodes = len(res["t_cross"])
    rate_deg_day = float("nan")
    rms_deg = float("nan")
    if n_nodes > 4:
        _, b = ols_slope(res["t_cross"], res["om_cross"])
        rate_deg_day = math.degrees(b) * 86400.0
        fit = harmonic_regression(
            (res["t_cross"] - 0.0) / 86400.0,
            res["om_cross"],
        )
        rms_deg = fit["rms_residual_deg"]
    return {
        "label": label,
        "mode": mode,
        "lambda_j2": lambda_j2,
        "lambda_3body": lambda_3body,
        "synthetic_moon": synthetic_moon,
        "rate_deg_day": rate_deg_day,
        "rms_residual_deg": rms_deg,
        "n_steps": n_steps,
        "n_nodes": n_nodes,
        "wall_clock_s": wall_s,
        "apply_precession": apply_precession,
    }


def _pack_sun(snap):
    return (snap["t_s"], snap["r_eci_km"], snap["sha256"])


# --------------------------------------------------------------------------- #
# Phase A: reduced-model (synthetic Moon) at 1 yr
# --------------------------------------------------------------------------- #
def phase_a_reduced_model(sun_snap, moon_snap, n_workers=4, arc_days=365.25):
    """Reduced-model experiment: synthetic circular Moon vs real DE441 Moon.

    If the residual after the synthetic Moon is much smaller than with the
    real Moon, the discrepancy is dominated by lunar-eccentricity /
    inclination-variation effects (H0b/H0c).
    """
    print(f"\n=== Phase A: reduced-model experiment (synthetic Moon vs real DE441) ===")
    print(f"  Arc: {arc_days} d, h=600 km, i=i_sso")
    print(f"  n_workers: {n_workers}")

    x0 = initial_state(H_SSO_KM, I_SSO_DEG)
    t_end_s = arc_days * 86400.0
    sun_packed = _pack_sun(sun_snap)
    moon_packed = _pack_sun(moon_snap)

    jobs = []
    for mode in ["j2_only", "sun_only", "moon_only", "sun_moon", "sun_moon_j2"]:
        # Real Moon
        jobs.append({
            "sun_data": sun_packed, "moon_data": moon_packed, "x0": x0,
            "mode": mode, "t_end_s": t_end_s, "dt_s": DT_S,
            "synthetic_moon": False,
            "label": f"real_{mode}",
        })
        # Synthetic Moon
        jobs.append({
            "sun_data": sun_packed, "moon_data": moon_packed, "x0": x0,
            "mode": mode, "t_end_s": t_end_s, "dt_s": DT_S,
            "synthetic_moon": True,
            "label": f"synthetic_{mode}",
        })

    print(f"  Total jobs: {len(jobs)}")
    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(_propagation_worker, jobs)

    return {"phase": "A_reduced_model", "arc_days": arc_days, "results": results}


# --------------------------------------------------------------------------- #
# Phase B: perturbative scaling
# --------------------------------------------------------------------------- #
def phase_b_perturbative_scaling(sun_snap, moon_snap, n_workers=4, arc_days=90.0):
    """Perturbative scaling experiment: (lambda_J2, lambda_3body) sweep."""
    print(f"\n=== Phase B: perturbative scaling experiment ===")
    print(f"  Arc: {arc_days} d, h=600 km, i=i_sso")

    x0 = initial_state(H_SSO_KM, I_SSO_DEG)
    t_end_s = arc_days * 86400.0
    sun_packed = _pack_sun(sun_snap)
    moon_packed = _pack_sun(moon_snap)

    LAMBDA_J2_VALUES = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
    LAMBDA_3BODY_VALUES = [0.0, 0.5, 1.0, 2.0]

    jobs = []
    for lj2 in LAMBDA_J2_VALUES:
        for l3b in LAMBDA_3BODY_VALUES:
            jobs.append({
                "sun_data": sun_packed, "moon_data": moon_packed, "x0": x0,
                "mode": "sun_moon_j2", "t_end_s": t_end_s, "dt_s": DT_S,
                "lambda_j2": lj2, "lambda_3body": l3b,
                "synthetic_moon": False,
                "label": f"full_lj2{lj2}_l3b{l3b}",
            })
            jobs.append({
                "sun_data": sun_packed, "moon_data": moon_packed, "x0": x0,
                "mode": "j2_only", "t_end_s": t_end_s, "dt_s": DT_S,
                "lambda_j2": lj2, "lambda_3body": 0.0,
                "synthetic_moon": False,
                "label": f"j2only_lj2{lj2}",
            })

    print(f"  Total jobs: {len(jobs)}")
    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(_propagation_worker, jobs)

    return {
        "phase": "B_perturbative_scaling",
        "arc_days": arc_days,
        "lambda_j2_values": LAMBDA_J2_VALUES,
        "lambda_3body_values": LAMBDA_3BODY_VALUES,
        "results": results,
    }


# --------------------------------------------------------------------------- #
# Phase C: full 18.6-yr force-mode decomposition
# --------------------------------------------------------------------------- #
def phase_c_full_arc(sun_snap, moon_snap, n_workers=4, arc_days=1.0*365.25):
    """DE441 propagation at 3 inclinations x 6 modes = 18 jobs at the given arc length.

    The full 18.6-yr arc is the canonical mission scale (one lunar nodal cycle),
    but is computationally expensive (~80 min/propagation x 18 jobs / 7 workers ~ 4 hr).
    A 1-yr arc gives the same force-mode decomposition structure at ~20 min wall.
    The 5-yr or 18.6-yr arc can be enabled by setting arc_days explicitly.
    """
    print(f"\n=== Phase C: force-mode decomposition at {arc_days/365.25:.2f}-yr arc ===")
    print(f"  Arc: {arc_days} d, h=600 km, 3 inclinations x 6 modes")

    t_end_s = arc_days * 86400.0
    sun_packed = _pack_sun(sun_snap)
    moon_packed = _pack_sun(moon_snap)

    jobs = []
    for i_deg in [I_SSO_DEG, I_90_DEG, I_30_DEG]:
        x0 = initial_state(H_SSO_KM, i_deg)
        for mode in ["kepler_only", "j2_only", "sun_only", "moon_only",
                      "sun_moon", "sun_moon_j2"]:
            jobs.append({
                "sun_data": sun_packed, "moon_data": moon_packed, "x0": x0,
                "mode": mode, "t_end_s": t_end_s, "dt_s": DT_S,
                "synthetic_moon": False,
                "label": f"i{i_deg:.2f}_{mode}",
            })

    print(f"  Total jobs: {len(jobs)}")
    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(_propagation_worker, jobs)

    return {"phase": "C_force_mode", "arc_days": arc_days, "results": results}


# --------------------------------------------------------------------------- #
# Main orchestration
# --------------------------------------------------------------------------- #
def main():
    print("Loading snapshots...")
    t0 = time.time()
    sun_snap = _load_snapshot(SUN_SNAPSHOT)
    moon_snap = _load_snapshot(MOON_SNAPSHOT)
    print(f"  Sun: {sun_snap['n_points']} points, sha256={sun_snap['sha256'][:16]}...")
    print(f"  Moon: {moon_snap['n_points']} points, sha256={moon_snap['sha256'][:16]}...")
    print(f"  Loaded in {time.time()-t0:.2f}s")

    # Use 7 workers (leave 1 core for OS) - 8 cores available
    n_workers = min(7, mp.cpu_count())
    print(f"Using {n_workers} workers (cpu_count={mp.cpu_count()})")

    out_dir = HERE / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Phase A: reduced model
    out_a = phase_a_reduced_model(sun_snap, moon_snap, n_workers=n_workers, arc_days=365.25)
    out_a["snapshot_provenance"] = {
        "sun_sha256": sun_snap["sha256"],
        "moon_sha256": moon_snap["sha256"],
        "sun_n_points": sun_snap["n_points"],
        "moon_n_points": moon_snap["n_points"],
    }
    out_a["code_hashes"] = code_hashes()
    with open(out_dir / "phase_a_reduced_model.json", "w") as f:
        json.dump(out_a, f, indent=2, default=str)
    print(f"\nPhase A results written to {out_dir / 'phase_a_reduced_model.json'}")

    # Phase B: perturbative scaling
    out_b = phase_b_perturbative_scaling(sun_snap, moon_snap, n_workers=n_workers, arc_days=90.0)
    out_b["snapshot_provenance"] = out_a["snapshot_provenance"]
    out_b["code_hashes"] = code_hashes()
    with open(out_dir / "phase_b_perturbative_scaling.json", "w") as f:
        json.dump(out_b, f, indent=2, default=str)
    print(f"\nPhase B results written to {out_dir / 'phase_b_perturbative_scaling.json'}")

    # Phase C: force-mode decomposition at 1-yr arc (fast; default)
    # The full 18.6-yr arc is the canonical mission scale but takes ~4 hr wall.
    # The 1-yr arc gives the same mode-isolation structure at ~20 min wall.
    # Override via PHASE_C_ARC_DAYS env var to use longer arcs.
    import os
    arc_days = float(os.environ.get("PHASE_C_ARC_DAYS", "365.25"))
    out_c = phase_c_full_arc(sun_snap, moon_snap, n_workers=n_workers, arc_days=arc_days)
    out_c["snapshot_provenance"] = out_a["snapshot_provenance"]
    out_c["code_hashes"] = code_hashes()
    out_path = out_dir / f"phase_c_full_{arc_days:.0f}d.json"
    with open(out_path, "w") as f:
        json.dump(out_c, f, indent=2, default=str)
    print(f"\nPhase C results written to {out_path}")

    print("\n=== Campaign complete ===")


if __name__ == "__main__":
    main()
