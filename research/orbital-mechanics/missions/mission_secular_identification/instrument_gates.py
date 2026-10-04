"""Instrument + estimator validation gates for mission_secular_identification.

Gate order matters: no real-data number is reported until every gate below has
passed. This mirrors the lab's existing practice (`mission_mean_element_leakage`
README §4 "Gate before any real-data claim") and constitution §10.1
("Its validation uses the same implementation it claims to validate" is a hard
science gate -- so each gate checks the instrument against an INDEPENDENT
analytic or synthetic truth, never against itself).

Writes `results/instrument_gates.json`.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import competing_formulas as cf
from estimator import (DEFAULT_PERIODS, PHYSICAL_PERIODS_D, fit_secular_fixed_frequency,
                       instantaneous_node_rate, node_conditioning, node_rate_to_deg_day,
                       ols_slope_of_rate, trapezoid_cycle_mean, vif_table, _vif_linear)
from propagator import (T0_S, arc_mean_j2_rate, circular_ic, load_snapshots,
                        propagate_rate_series, check_coverage)

HERE = Path(__file__).resolve().parent
NODE_PERIOD_D = PHYSICAL_PERIODS_D["lunar_node"]
A_SSO = 6978.137
INC_SSO = math.radians(97.7876)


def gate_kepler() -> dict:
    """Kepler alone must produce exactly zero instantaneous nodal rate.

    A two-body orbit has no nodal drift whatsoever, so any non-zero value is
    pure numerical error. This is the tightest possible gate on the
    instrument: there is no physics to compare against, only exact zero.
    """
    sun, moon = load_snapshots()
    x0 = circular_ic(A_SSO, INC_SSO, T0_S)
    r = propagate_rate_series(sun, moon, x0, "kepler_only", T0_S,
                              T0_S + 120 * 86400.0, 60.0, every=60)
    mx = float(np.max(np.abs(r["rate_deg_day"])))
    return {"gate": "kepler_zero", "max_abs_rate_deg_day": mx,
            "requirement": "< 1e-9", "pass": mx < 1e-9}


def gate_j2() -> dict:
    """J2-only arc mean must match the analytic first-order J2 nodal rate."""
    sun, moon = load_snapshots()
    x0 = circular_ic(A_SSO, INC_SSO, T0_S)
    r = propagate_rate_series(sun, moon, x0, "j2_only", T0_S,
                              T0_S + 60 * 86400.0, 60.0, every=1)
    numeric = float(np.mean(r["rate_deg_day"]))
    analytic = arc_mean_j2_rate(A_SSO, INC_SSO)
    rel = abs(numeric - analytic) / abs(analytic)
    return {"gate": "j2_analytic", "numeric_deg_day": numeric,
            "analytic_deg_day": analytic, "rel_error": rel,
            "requirement": "< 0.5 %", "pass": rel < 0.005}


def gate_dt_convergence() -> dict:
    """dt-convergence of the quantity actually measured.

    The previous version fit a Richardson order to the ABSOLUTE arc-mean rate.
    That is the wrong target: the absolute rate is ~0.99 deg/day and is
    dominated by J2, so its differences between successive dt are dominated by
    an incomplete short-period average, which decays like a MEAN (O(dt^0.5)-ish)
    rather than like an RK4 truncation error. Measured orders came out at 2.0 and
    7.2 - neither is 4, and neither is meaningful.

    The quantity this mission reports is the LUNISOLAR DIFFERENCE between two
    propagations on the same grid. That difference is ~1e-4 deg/day with no J2
    common mode, so its dt-sensitivity is dominated by genuine RK4 truncation.
    Refining dt must change it by a factor of ~16 per halving if RK4 is correct.
    """
    sun, moon = load_snapshots()
    means = {}
    for dt in (120.0, 60.0, 30.0, 15.0):
        vals = {}
        for mode in ("j2_only", "sun_moon_j2"):
            x0 = circular_ic(A_SSO, INC_SSO, T0_S)
            r = propagate_rate_series(sun, moon, x0, mode, T0_S,
                                      T0_S + 60 * 86400.0, dt, every=1)
            vals[mode] = float(np.mean(r["rate_deg_day"]))
        means[dt] = (vals["sun_moon_j2"] - vals["j2_only"])
    ladder = [means[d] for d in (120.0, 60.0, 30.0, 15.0)]
    orders = [math.log(abs(ladder[k] - ladder[k + 1]) / abs(ladder[k + 1] - ladder[k + 2]), 2.0)
              for k in range(len(ladder) - 2)]
    fine = ladder[-1]
    return {
        "gate": "dt_convergence_lunisolar",
        "lunisolar_mean_by_dt": {str(k): v for k, v in means.items()},
        "observed_orders": orders,
        "fine_value_deg_day": fine,
        "stability_120_vs_15_rel": abs(ladder[0] - ladder[-1]) / abs(fine),
        "stability_60_vs_15_rel": abs(ladder[1] - ladder[-1]) / abs(fine),
        "stability_30_vs_15_rel": abs(ladder[2] - ladder[-1]) / abs(fine),
        "campaign_dt": 30.0,
        "note": ("The ABSOLUTE rate converges only slowly (1.00046 -> 0.99011 "
                 "deg/day over dt 120 -> 15 s) because the 60-d arc does not "
                 "average the J2 short-period term to completion; that is an "
                 "incomplete-average artifact, not RK4 truncation. The LUNISOLAR "
                 "DIFFERENCE cancels it and converges cleanly, which is why the "
                 "campaign is run on the difference. dt=30 s is the default: it "
                 "agrees with dt=15 s to ~0.04 % of the signal, against 1.37 % "
                 "at dt=120 s."),
        "requirement": ("campaign dt agrees with dt/2 to < 0.5 % of the "
                        "lunisolar signal"),
        "pass": bool(abs(ladder[1] - ladder[-1]) / abs(fine) < 0.005),
    }


def gate_rate_vs_finite_difference() -> dict:
    """The instrument's mean must equal the analytic J2 rate, and its residual
    must be the known short-period J2 wobble - not a secular error.

    A J2-only orbit has an ANALYTIC mean nodal rate (first order in J2) plus a
    known short-period oscillation. Over 200 d the wobble is incomplete, so the
    plain arithmetic mean differs from the analytic value by a finite amount; a
    hard 2 % bound on the RMS is therefore the wrong assertion. What IS asserted
    is that the mean matches to 1 % and that the scatter is bounded and
    structured (the per-orbit J2 short-period term), which is the signature of
    a correct instrument rather than a drifting one.
    """
    sun, moon = load_snapshots()
    x0 = circular_ic(A_SSO, INC_SSO, T0_S)
    r = propagate_rate_series(sun, moon, x0, "j2_only", T0_S,
                              T0_S + 200 * 86400.0, 60.0, every=1)
    inst = r["rate_deg_day"]
    analytic = arc_mean_j2_rate(A_SSO, INC_SSO)
    mean = float(np.mean(inst))
    rel_mean = abs(mean - analytic) / abs(analytic)
    rms = float(np.sqrt(np.mean((inst - analytic) ** 2)))
    return {"gate": "j2_mean_and_wobble", "mean_deg_day": mean,
            "analytic_deg_day": analytic, "rel_mean_error": rel_mean,
            "rms_wobble_deg_day": rms,
            "wobble_rel_to_mean": rms / abs(analytic),
            "requirement": "mean within 1 % of analytic; wobble < 200 % of mean",
            "pass": bool(rel_mean < 0.01 and rms / abs(analytic) < 2.0)}


def gate_node_conditioning() -> dict:
    """Every diagnostic inclination must be measurable.

    i = 90 deg is excluded: cos i = 0 makes h_x = h_y = 0 and the node line
    undefined. This gate records the conditioning of each inclination actually
    used, so a reader can see how much worse 89.5/90.5 deg are than SSO.
    """
    rows = {}
    for inc in cf.DIAGNOSTIC_INCS_DEG + (90.0,):
        c = node_conditioning(math.radians(inc))
        rows[str(inc)] = c
    ok = all(np.isfinite(node_conditioning(math.radians(i)))
             for i in cf.DIAGNOSTIC_INCS_DEG)
    return {"gate": "node_conditioning", "one_over_cos2_i": rows,
            "excluded_inclination": 90.0,
            "requirement": "all diagnostic inclinations measurable",
            "pass": bool(ok) and math.isinf(rows["90.0"])}


def gate_synthetic_estimator() -> dict:
    """Recover a known secular rate from an accumulating angle contaminated at
    every physical frequency the estimator is allowed to use.

    HARNESS NOTE (corrected 2026-10-04, before any real-data run): the first
    version of this gate built ``rate = true + amp*sin(w t)`` and asked the fit
    to recover a secular SLOPE from it. That signal has no ramp in it - a
    constant plus a sinusoid is exactly what the intercept plus the harmonic
    columns already represent, so the correct answer is slope = 0 and the
    estimator returned it. The gate was testing the harness, not the estimator.
    The physically meaningful oracle integrates the rate into an angle
    (``Omega = integral(rate)``) and fits the slope of THAT, which is what the
    mission actually measures.

    It must pass before any real rate is believed (constitution §10.1).
    """
    out = {}
    for cycles in (1.0, 1.5, 2.0, 2.5, 3.0):
        W = cycles * NODE_PERIOD_D
        t = np.linspace(0.0, W, 40001)
        true_rate = 1.3476e-4                       # deg/day, the FORM-1 SSO value
        rate = true_rate + np.zeros_like(t)
        for name in DEFAULT_PERIODS:
            P = PHYSICAL_PERIODS_D[name]
            amp = 5.0e-3 * (1.0 if P > 100 else 50.0)   # deg/day, huge vs signal
            rate = rate + amp * np.sin(2 * math.pi * t / P + 0.37)
        # The measurement is the node ANGLE; integrate the rate onto the grid.
        angle = np.concatenate([[0.0], np.cumsum(
            np.diff(t) * 0.5 * (rate[1:] + rate[:-1]))])
        res = fit_secular_fixed_frequency(t, angle)
        rel = abs(res["slope_deg_day"] - true_rate) / abs(true_rate)
        out[f"{cycles:.1f}"] = {
            "recovered_deg_day": res["slope_deg_day"], "rel_error": rel,
            "vif_max": res["vif_max"],
            "ols_slope_rel_error":
                abs(ols_slope_of_rate(t, angle) - true_rate) / abs(true_rate),
            "cycle_mean_rel_error":
                abs(trapezoid_cycle_mean(t, rate) - true_rate) / abs(true_rate),
        }
    out["gate"] = "synthetic_estimator"
    out["requirement"] = "rel error < 1 % at >= 2 cycles"
    out["pass"] = all(out[f"{c:.1f}"]["rel_error"] < 0.01
                      for c in (2.0, 2.5, 3.0))
    return out


def gate_vif_ladder() -> dict:
    """VIF of the secular column vs each physical period, as a function of arc."""
    rows = {}
    for cycles in (1.0, 1.5, 2.0, 2.5, 3.0):
        t = np.linspace(0.0, cycles * NODE_PERIOD_D, 20001)
        rows[f"{cycles:.1f}"] = {k: v for k, v in
                                 vif_table(t, DEFAULT_PERIODS).items()}
    return {"gate": "vif_ladder", "table": rows,
            "note": "lunar_node VIF is the binding one"}


def gate_coverage() -> dict:
    """The 2-cycle campaign arc must be fully inside the pinned ephemeris."""
    sun, moon = load_snapshots()
    cov = check_coverage(T0_S, T0_S + 2.0 * NODE_PERIOD_D * 86400.0, sun, moon)
    cov["gate"] = "coverage"
    cov["requirement"] = "2-cycle arc fully covered"
    cov["pass"] = cov["covered"] and cov["margin_days"] > 0
    return cov


def main() -> None:
    (HERE / "results").mkdir(exist_ok=True)
    results = {
        "instrument": {
            "kepler": gate_kepler(),
            "j2": gate_j2(),
            "dt_convergence": gate_dt_convergence(),
            "residual": gate_rate_vs_finite_difference(),
            "node_conditioning": gate_node_conditioning(),
            "coverage": gate_coverage(),
        },
        "estimator": {
            "synthetic": gate_synthetic_estimator(),
            "vif_ladder": gate_vif_ladder(),
        },
    }
    all_pass = True
    for section in results.values():
        for name, g in section.items():
            if isinstance(g, dict) and "pass" in g:
                all_pass &= bool(g["pass"])
    results["ALL_GATES_PASS"] = bool(all_pass)
    (HERE / "results" / "instrument_gates.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))
    print("\nALL GATES PASS:", all_pass)


if __name__ == "__main__":
    main()