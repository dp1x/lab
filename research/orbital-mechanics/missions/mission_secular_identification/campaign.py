"""The >= 2-cycle campaign for mission_secular_identification.

For each inclination x phase, two propagations are run on the SAME grid:
``sun_moon_j2`` and ``j2_only``. The lunisolar rate is their difference, and
because both are estimated with the same estimator on the same grid the J2
common-mode cancels exactly. That is why this mission needs two force modes
per case rather than the five the coupling mission needed.

Every case is scored against BOTH competing formulas with the mission card's
pre-registered decision rule. Neither formula is fitted to anything.

Parallelism: `multiprocessing.Pool` over (inclination, phase) cases, one worker
per CPU. Determinism: fixed dt, fixed initial conditions, no RNG anywhere.

Writes `results/campaign.json`.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import competing_formulas as cf
from estimator import (DEFAULT_PERIODS, PHYSICAL_PERIODS_D, fit_secular_fixed_frequency,
                       node_conditioning, trapezoid_cycle_mean)
from propagator import (T0_S, check_coverage, circular_ic, load_snapshots,
                        propagate_rate_series)

HERE = Path(__file__).resolve().parent
NODE_PERIOD_D = PHYSICAL_PERIODS_D["lunar_node"]
A_SSO = 6978.137
MODES = ("j2_only", "sun_moon_j2")

#: Argument-of-latitude phases. Four phases per inclination, separated by a
#: quarter orbit, so no conclusion can rest on one favourable initial condition.
PHASES_DEG = (0.0, 90.0, 180.0, 270.0)


def _run_mode(sun, moon, inc_deg, phase_deg, mode, t_end_s, dt_s, every):
    x0 = circular_ic(A_SSO, math.radians(inc_deg), T0_S, math.radians(phase_deg))
    r = propagate_rate_series(sun, moon, x0, mode, T0_S, t_end_s, dt_s, every=every)
    t_days = (r["t_s"] - T0_S) / 86400.0
    return t_days, r["rate_deg_day"]


def _case(arg) -> dict:
    inc_deg, phase_deg, cycles, dt_s, every = arg
    sun, moon = load_snapshots()
    t_end_s = T0_S + cycles * NODE_PERIOD_D * 86400.0
    cov = check_coverage(T0_S, t_end_s, sun, moon)
    if not cov["covered"]:
        raise RuntimeError(f"arc not covered by pinned ephemeris: {cov}")

    out = {"inc_deg": inc_deg, "phase_deg": phase_deg, "cycles": cycles,
           "dt_s": dt_s, "nodal_conditioning": node_conditioning(math.radians(inc_deg)),
           "coverage": cov}
    per_mode = {}
    for mode in MODES:
        t, rate = _run_mode(sun, moon, inc_deg, phase_deg, mode, t_end_s, dt_s, every)
        fit = fit_secular_fixed_frequency(t, rate, DEFAULT_PERIODS)
        fit["cycle_mean_deg_day"] = trapezoid_cycle_mean(t, rate)
        fit["n_samples"] = int(t.size)
        per_mode[mode] = fit
    out["modes"] = per_mode

    # Lunisolar = difference of the two identically-estimated rates.
    f_full = per_mode["sun_moon_j2"]
    f_j2 = per_mode["j2_only"]
    luni = f_full["slope_deg_day"] - f_j2["slope_deg_day"]
    luni_se = math.hypot(f_full["slope_stderr_deg_day"], f_j2["slope_stderr_deg_day"])
    out["lunisolar"] = {
        "secular_deg_day": luni,
        "stderr_deg_day": luni_se,
        "abs_stderr_deg_day": abs(luni_se),
        "vif_max": max(f_full["vif_max"], f_j2["vif_max"]),
        "cycle_mean_deg_day": f_full["cycle_mean_deg_day"] - f_j2["cycle_mean_deg_day"],
        "max_abs_harmonic_deg_day": max(f_full["max_abs_harmonic_deg_day"],
                                        f_j2["max_abs_harmonic_deg_day"]),
    }
    f1 = cf.form1_audit018(A_SSO, inc_deg)
    f2 = cf.form2_competitor(A_SSO, inc_deg)
    out["formulas"] = {
        "form1_deg_day": f1, "form2_deg_day": f2,
        "measured_over_form1": (luni / f1) if f1 != 0 else math.inf,
        "measured_over_form2": (luni / f2) if f2 != 0 else math.inf,
        "same_sign_form1": bool(luni * f1 > 0),
        "same_sign_form2": bool(luni * f2 > 0),
        "ratio_in_0p5_to_2_form1": (0.5 <= abs(luni / f1) <= 2.0) if f1 != 0 else False,
        "ratio_in_0p5_to_2_form2": (0.5 <= abs(luni / f2) <= 2.0) if f2 != 0 else False,
    }
    return out


def _window_ladder(arg) -> dict:
    """Same propagation re-estimated on sub-windows of increasing cycle count.

    Identifiability means the estimate STOPS MOVING as the window grows. If it
    keeps drifting, the secular term is not separated from the forced content.
    """
    inc_deg, phase_deg, dt_s, every = arg
    sun, moon = load_snapshots()
    ladder = {}
    for cycles in (1.0, 1.5, 2.0, 2.5, 3.0):
        try:
            t_end_s = T0_S + cycles * NODE_PERIOD_D * 86400.0
            if not check_coverage(T0_S, t_end_s, sun, moon)["covered"]:
                continue
            vals, vifs = {}, []
            for mode in MODES:
                t, rate = _run_mode(sun, moon, inc_deg, phase_deg, mode,
                                    t_end_s, dt_s, every)
                f = fit_secular_fixed_frequency(t, rate, DEFAULT_PERIODS)
                vals[mode] = f["slope_deg_day"]
                vifs.append(f["vif_max"])
            ladder[f"{cycles:.1f}"] = {
                "lunisolar_deg_day": vals["sun_moon_j2"] - vals["j2_only"],
                "vif_max": max(vifs),
            }
        except RuntimeError:
            continue
    return {"inc_deg": inc_deg, "phase_deg": phase_deg, "ladder": ladder}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cycles", type=float, default=2.0,
                    help="Arc length in lunar nodal cycles.")
    ap.add_argument("--dt", type=float, default=30.0,
                    help="RK4 step (s). 30 s is the gate-validated default: it "
                         "agrees with dt=15 s to ~0.04 %% of the lunisolar "
                         "signal. dt=120 s is 1.37 %% off and dt=240 s is "
                         "catastrophically wrong (see results/instrument_gates.json).")
    ap.add_argument("--every", type=int, default=20, help="Sample every N steps.")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--incs", default=",".join(str(i) for i in cf.DIAGNOSTIC_INCS_DEG))
    ap.add_argument("--phases", default=",".join(str(p) for p in PHASES_DEG))
    ap.add_argument("--ladder", action="store_true",
                    help="Also run the window ladder (slow).")
    args = ap.parse_args()

    incs = [float(x) for x in args.incs.split(",")]
    phases = [float(x) for x in args.phases.split(",")]

    # Sample density: dt/every is the effective output cadence.
    out_cadence_s = args.dt * args.every
    (HERE / "results").mkdir(exist_ok=True)

    cases = [(i, p, args.cycles, args.dt, args.every) for i in incs for p in phases]
    t0 = time.time()
    with Pool(args.workers) as pool:
        rows = pool.map(_case, cases)
    wall = time.time() - t0

    ladders = []
    if args.ladder:
        lcases = [(i, p, args.dt, args.every) for i in incs for p in phases[:2]]
        with Pool(args.workers) as pool:
            ladders = pool.map(_window_ladder, lcases)

    payload = {
        "mission": "mission_secular_identification",
        "arc": {"t0_s": T0_S, "nodal_period_d": NODE_PERIOD_D,
                "cycles": args.cycles,
                "arc_days": args.cycles * NODE_PERIOD_D,
                "arc_years": args.cycles * NODE_PERIOD_D / 365.25},
        "integrator": {"scheme": "fixed-step RK4 on the full 6-vector",
                       "dt_s": args.dt, "sample_every": args.every,
                       "effective_cadence_s": out_cadence_s,
                       "n_samples_per_mode_per_case": int(
                           args.cycles * NODE_PERIOD_D * 86400.0 / out_cadence_s)},
        "estimator": {"periods_d": {k: PHYSICAL_PERIODS_D[k] for k in DEFAULT_PERIODS},
                      "frequencies": "FIXED at physical values; no empirical fitting"},
        "cases": rows,
        "window_ladder": ladders,
        "provenance": {
            "code_hash_scheme": "lf-normalized-v1",
            "code": {p.name: hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()[:16]
                for p in sorted(HERE.glob("*.py"))},
            "sun_sha256": load_snapshots()[0]["sha256"],
            "moon_sha256": load_snapshots()[1]["sha256"],
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
            "compute": "local",
            "rng_seed": None,
            "wall_clock_s": round(wall, 1),
            "n_workers": args.workers,
            "incs_deg": incs, "phases_deg": phases,
        },
    }
    name = f"campaign_{args.cycles:g}cyc.json" if not args.ladder else \
           f"campaign_{args.cycles:g}cyc_ladder.json"
    (HERE / "results" / name).write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # console summary
    print(f"arc: {args.cycles} nodal cycles = {args.cycles * NODE_PERIOD_D:.1f} d "
          f"({args.cycles * NODE_PERIOD_D / 365.25:.2f} yr), dt={args.dt}s, "
          f"cadence={out_cadence_s / 3600:.2f} h, wall={wall:.0f}s")
    print(f"{'inc':>8} {'phase':>6} {'luni deg/day':>16} {'|stderr|':>11} "
          f"{'VIF':>6} {'FORM-1':>13} {'m/f1':>9} {'FORM-2':>13} {'m/f2':>9}")
    for r in rows:
        L, F = r["lunisolar"], r["formulas"]
        print(f"{r['inc_deg']:>8.3f} {r['phase_deg']:>6.1f} "
              f"{L['secular_deg_day']:>16.6e} {L['abs_stderr_deg_day']:>11.2e} "
              f"{L['vif_max']:>6.2f} {F['form1_deg_day']:>13.5e} "
              f"{F['measured_over_form1']:>9.3f} {F['form2_deg_day']:>13.5e} "
              f"{F['measured_over_form2']:>9.3f}")
    if ladders:
        print("\nwindow ladder (lunisolar deg/day):")
        for lad in ladders:
            seq = "  ".join(f"{k}:{v['lunisolar_deg_day']:+.3e}"
                            for k, v in sorted(lad["ladder"].items(),
                                               key=lambda kv: float(kv[0])))
            print(f"  i={lad['inc_deg']:>8.3f} phase={lad['phase_deg']:>5.1f}  {seq}")
    print(f"\nwall = {wall:.0f}s on {args.workers} workers -> {name}")


if __name__ == "__main__":
    sys.exit(main())