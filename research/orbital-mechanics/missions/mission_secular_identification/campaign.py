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
    inc_deg, phase_deg, cycles, dt_s, every, ladder_spec = arg
    sun, moon = load_snapshots()
    t_end_s = T0_S + cycles * NODE_PERIOD_D * 86400.0
    cov = check_coverage(T0_S, t_end_s, sun, moon)
    if not cov["covered"]:
        raise RuntimeError(f"arc not covered by pinned ephemeris: {cov}")

    out = {"inc_deg": inc_deg, "phase_deg": phase_deg, "cycles": cycles,
           "dt_s": dt_s, "nodal_conditioning": node_conditioning(math.radians(inc_deg)),
           "coverage": cov}

    # Propagate ONCE per mode over the LONGEST arc, then derive every window by
    # TRUNCATION. This is exact, not an approximation: fixed-step RK4 advances a
    # state independently of the horizon, so the state sequence of a 3-cycle run
    # has a 2-cycle run as its bit-identical prefix. Verified:
    #     np.array_equal(short['rate_deg_day'], long['rate_deg_day'][:n]) -> True
    # This turns a 5x-propagation window ladder into a 1x cost.
    window_cycles = [c for c, _ in ladder_spec if c <= cycles]
    window_cycles.append(cycles)
    window_cycles = sorted(set(window_cycles))

    per_mode = {}
    ladder = {}
    for mode in MODES:
        x0 = circular_ic(A_SSO, math.radians(inc_deg), T0_S, math.radians(phase_deg))
        r = propagate_rate_series(sun, moon, x0, mode, T0_S, t_end_s, dt_s, every=every)
        t_all = (r["t_s"] - T0_S) / 86400.0
        rate_all = r["rate_deg_day"]
        # the headline window
        fit = fit_secular_fixed_frequency(t_all, rate_all, DEFAULT_PERIODS)
        fit["cycle_mean_deg_day"] = trapezoid_cycle_mean(t_all, rate_all)
        fit["n_samples"] = int(t_all.size)
        per_mode[mode] = fit
        # every sub-window
        for c in window_cycles:
            n = int(round(c * NODE_PERIOD_D * 86400.0 / (dt_s * every))) + 1
            n = min(n, t_all.size)
            tw, rw = t_all[:n], rate_all[:n]
            fw = fit_secular_fixed_frequency(tw, rw, DEFAULT_PERIODS)
            ladder.setdefault(f"{c:.2f}", {})[mode] = {
                "slope_deg_day": fw["slope_deg_day"],
                "vif_max": fw["vif_max"],
                "cycle_mean_deg_day": trapezoid_cycle_mean(tw, rw),
                "n_samples": int(tw.size),
            }

    out["modes"] = per_mode
    out["window_ladder"] = {
        k: {"lunisolar_deg_day": v["sun_moon_j2"]["slope_deg_day"]
               - v["j2_only"]["slope_deg_day"],
           "vif_max": max(v["sun_moon_j2"]["vif_max"], v["j2_only"]["vif_max"]),
           "cycle_mean_lunisolar_deg_day": v["sun_moon_j2"]["cycle_mean_deg_day"]
               - v["j2_only"]["cycle_mean_deg_day"],
           "n_samples": v["sun_moon_j2"]["n_samples"]}
        for k, v in sorted(ladder.items(), key=lambda kv: float(kv[0]))
    }

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


def _run_and_checkpoint(args) -> tuple:
    """Run ONE case, write it to its own file immediately, return its path.

    CHECKPOINTING (added 2026-10-05 after a ~64 CPU-hour run was lost).
    The first version of this campaign buffered every case through `pool.map`
    and wrote a single JSON only after all of them returned. Terminating it at
    10 hr produced NO recoverable state at all -- not even the 12 of 16 cases
    that had already finished. Each case is now written to its own file the
    moment it completes, so an interrupted run keeps everything it earned and a
    restart resumes rather than recomputes.
    """
    import tempfile

    arg, ckpt_dir = args
    inc_deg, phase_deg, cycles, dt_s, every, ladder_spec = arg
    tag = f"i{inc_deg:.4f}_ph{phase_deg:05.1f}".replace(".", "p").replace("-", "m")
    ckpt = Path(ckpt_dir) / f"case_{tag}.json"
    if ckpt.exists():
        return str(ckpt)                     # resume: already paid for this case
    t0 = time.time()
    row = _case(arg)
    row["wall_s"] = round(time.time() - t0, 1)
    ckpt.write_text(json.dumps(row, indent=2), encoding="utf-8")
    print(f"    [done] i={inc_deg:.3f} phase={phase_deg:g} "
          f"luni={row['lunisolar']['secular_deg_day']:+.4e} deg/day "
          f"({row['wall_s']:.0f}s)", flush=True)
    return str(ckpt)


#: Measured cost constants for this host (2026-10-05, calibrated against a real
#: 20-day benchmark: 19.8 s for 57 600 steps at dt=30 s in sun_moon_j2 mode).
#: The j2_only mode skips ephemeris interpolation and precession entirely and
#: measures ~0.35x the cost.
US_PER_STEP_SUN_MOON_J2 = 344.0
J2_ONLY_COST_RATIO = 0.35


def estimate_cost(cases: int, cycles: float, dt_s: float, every: int,
                  workers: int, modes: int = 2) -> dict:
    """Preflight cost estimate. REFUSES to launch a run that overruns budget.

    WHY THIS EXISTS. The 2026-10-04 session measured one case, assumed its cost
    scaled trivially, and launched a 32-propagation matrix without multiplying
    out. It overran by ~3x and was terminated at ~64 CPU-hours with no results
    written, because the run buffered everything through `pool.map`. This
    function is the guard against that class of error: it runs before a single
    CPU-second is spent, and it is calibrated from a MEASURED per-step cost
    rather than an assumed one.

    Calibration: 19.8 s / 57 600 steps = 344 us/step for `sun_moon_j2` on this
    host, measured after the scalar hot-path rewrite (which gave 1.42x). The
    first version of this function used 2.0 us/step and overstated cost by ~170x,
    which would have made every run look impossible; the number is pinned by
    `tests/test_secular_identification.py::test_cost_estimate_is_calibrated`.
    """
    steps = cycles * 6798.383 * 86400.0 / dt_s
    per_step = (US_PER_STEP_SUN_MOON_J2 + US_PER_STEP_SUN_MOON_J2 * J2_ONLY_COST_RATIO) / 1e6
    per_case_s = steps * per_step
    total_s = cases * per_case_s
    wall_s = total_s / workers
    return {
        "n_cases": cases, "modes_per_case": modes, "cycles": cycles, "dt_s": dt_s,
        "steps_per_propagation": int(steps),
        "us_per_step_sun_moon_j2": US_PER_STEP_SUN_MOON_J2,
        "est_seconds_per_case": round(per_case_s, 1),
        "est_total_cpu_hours": round(total_s / 3600.0, 2),
        "est_wall_hours": round(wall_s / 3600.0, 2),
        "workers": workers,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cycles", type=float, default=2.0,
                    help="Headline arc length in lunar nodal cycles.")
    ap.add_argument("--dt", type=float, default=30.0,
                    help="RK4 step (s). 30 s is the gate-validated default: it "
                         "agrees with dt=15 s to ~0.04 %% of the lunisolar "
                         "signal. dt=120 s is 1.37 %% off and dt=240 s is "
                         "catastrophically wrong (see results/instrument_gates.json).")
    ap.add_argument("--every", type=int, default=20, help="Sample every N steps.")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--incs", default=",".join(str(i) for i in cf.DIAGNOSTIC_INCS_DEG))
    ap.add_argument("--phases", default=",".join(str(p) for p in PHASES_DEG))
    ap.add_argument("--ladder-max-cycles", type=float, default=0.0)
    ap.add_argument("--budget-hours", type=float, default=10.0,
                    help="Hard wall-clock budget. The run REFUSES to start if "
                         "the preflight estimate exceeds this (LAB_CONSTITUTION "
                         "section 4.3). Use --force-budget to override knowingly.")
    ap.add_argument("--force-budget", action="store_true",
                    help="Launch anyway despite exceeding the budget.")
    args = ap.parse_args()

    incs = [float(x) for x in args.incs.split(",")]
    phases = [float(x) for x in args.phases.split(",")]

    prop_cycles = max(args.cycles, args.ladder_max_cycles)
    ladder_spec = [(c, c) for c in (1.0, 1.5, 2.0, 2.5, 3.0)]
    out_cadence_s = args.dt * args.every
    (HERE / "results").mkdir(exist_ok=True)

    cases = [(i, p, prop_cycles, args.dt, args.every, ladder_spec)
             for i in incs for p in phases]

    # ---- PREFLIGHT: refuse an over-budget run before spending any CPU -----
    est = estimate_cost(len(cases), prop_cycles, args.dt, args.every,
                        args.workers, modes=len(MODES))
    print(f"PREFLIGHT: {len(cases)} cases x {len(MODES)} modes, "
          f"{prop_cycles:g} nodal cycles ({prop_cycles * NODE_PERIOD_D:.0f} d), "
          f"dt={args.dt:g}s, {args.workers} workers")
    print(f"  estimate: {est['est_wall_hours']:.2f} h wall "
          f"({est['est_total_cpu_hours']:.1f} CPU-h), "
          f"{est['steps_per_propagation']:,} steps/propagation")
    print(f"  budget:   {args.budget_hours:.2f} h")
    if est["est_wall_hours"] > args.budget_hours and not args.force_budget:
        print("\nREFUSING TO LAUNCH: estimate exceeds budget by "
              f"{est['est_wall_hours'] - args.budget_hours:.2f} h.")
        print("Reduce --incs / --phases / --cycles, raise --dt (only within the "
              "dt gate), or pass --force-budget if you are deliberately "
              "accepting the overrun.")
        sys.exit(2)

    ckpt_dir = HERE / "results" / "checkpoints"
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    payload_args = [(c, str(ckpt_dir)) for c in cases]
    t0 = time.time()
    with Pool(args.workers) as pool:
        paths = pool.map(_run_and_checkpoint, payload_args)
    wall = time.time() - t0
    rows = [json.loads(Path(p).read_text(encoding="utf-8")) for p in paths]
    rows.sort(key=lambda r: (r["inc_deg"], r["phase_deg"]))

    payload = {
        "mission": "mission_secular_identification",
        "arc": {"t0_s": T0_S, "nodal_period_d": NODE_PERIOD_D,
                "headline_cycles": args.cycles,
                "propagated_cycles": prop_cycles,
                "headline_arc_days": args.cycles * NODE_PERIOD_D,
                "headline_arc_years": args.cycles * NODE_PERIOD_D / 365.25},
        "integrator": {"scheme": "fixed-step RK4 on the full 6-vector",
                       "dt_s": args.dt, "sample_every": args.every,
                       "effective_cadence_s": out_cadence_s,
                       "ladder_by_truncation": True},
        "estimator": {"periods_d": {k: PHYSICAL_PERIODS_D[k] for k in DEFAULT_PERIODS},
                      "frequencies": "FIXED at physical values; no empirical fitting"},
        "cases": rows,
        "cost_preflight": est,
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
            "checkpointed": True,
            "checkpoint_dir": "results/checkpoints",
        },
    }
    name = f"campaign_{args.cycles:g}cyc.json"
    (HERE / "results" / name).write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # console summary
    print(f"headline arc: {args.cycles} nodal cycles = {args.cycles * NODE_PERIOD_D:.1f} d "
          f"({args.cycles * NODE_PERIOD_D / 365.25:.2f} yr); propagated {prop_cycles} cyc "
          f"for the ladder; dt={args.dt}s, cadence={out_cadence_s / 3600:.2f} h, "
          f"wall={wall:.0f}s")
    print(f"{'inc':>8} {'phase':>6} {'luni deg/day':>16} {'|stderr|':>11} "
          f"{'VIF':>6} {'FORM-1':>13} {'m/f1':>9} {'FORM-2':>13} {'m/f2':>9}")
    for r in rows:
        L, F = r["lunisolar"], r["formulas"]
        print(f"{r['inc_deg']:>8.3f} {r['phase_deg']:>6.1f} "
              f"{L['secular_deg_day']:>16.6e} {L['abs_stderr_deg_day']:>11.2e} "
              f"{L['vif_max']:>6.2f} {F['form1_deg_day']:>13.5e} "
              f"{F['measured_over_form1']:>9.3f} {F['form2_deg_day']:>13.5e} "
              f"{F['measured_over_form2']:>9.3f}")
    if rows and rows[0].get("window_ladder"):
        print("\nwindow ladder (lunisolar deg/day, by truncation):")
        for r in rows:
            seq = "  ".join(f"{k}:{v['lunisolar_deg_day']:+.3e}"
                            for k, v in r["window_ladder"].items())
            print(f"  i={r['inc_deg']:>8.3f} phase={r['phase_deg']:>5.1f}  {seq}")
    print(f"\nwall = {wall:.0f}s on {args.workers} workers -> {name}")


if __name__ == "__main__":
    sys.exit(main())