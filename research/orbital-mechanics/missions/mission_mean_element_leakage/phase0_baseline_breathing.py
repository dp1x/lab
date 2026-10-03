"""Phase 0 - baseline-breathing pre-check for mission_mean_element_leakage.

Hypothesis H-baseline: the lunisolar chain's residual

    R = (sun_moon_j2 - j2_only) - [(sun_only - kepler) + (moon_only - kepler)]

is (at least partly) the error we make by subtracting a CONSTANT J2 node rate
measured from a J2-only run, while the real run's mean a and i drift under
Sun/Moon forcing.  Sensitivity at h=600 km, i=97.79 deg:

    dOmegaJ2/di ~ +0.126 deg/day per degree
    dOmegaJ2/da ~ -4.9e-4 deg/day per km
    (a 0.010 deg mean-tilt change already equals the 1-yr SSO R)

This script propagates the same force modes as the campaign with one shared
integrator, dt and snapshot, records osculating a, e, i at every ascending-node
crossing, and computes

    D = < OmegaJ2(a_full(t), i_full(t)) - OmegaJ2(a_j2(t), i_j2(t)) >

i.e. the spurious rate that the constant-baseline subtraction injects, next to
the measured node-drift residual R from the same runs.

Constants, DE441 snapshots and the interpolation-bug-corrected snapshot
interpolation are imported from the predecessor mission, so both missions share
exactly one physics source. Deterministic; no timing telemetry is written.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

# --- bootstrap ------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parents[4]          # repo root
HERE = Path(__file__).resolve().parent
PRED = HERE.parent / "mission_j2_precession_modulated_lunisolar_term"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(PRED))
_spec = importlib.util.spec_from_file_location("pred_mission_experiment", PRED / "mission_experiment.py")
ME = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = ME
_spec.loader.exec_module(ME)

from lab_utils import J2_EARTH, MU_EARTH_KM3S2, R_EARTH_KM  # noqa: E402
from lab_utils.integrators import rk4_step  # noqa: E402

T0_S = 820476800.0
H_KM = 600.0
INCS = [97.7876, 90.0, 30.0]
MODES = ["kepler_only", "j2_only", "sun_only", "moon_only", "sun_moon_j2"]
RE2 = R_EARTH_KM ** 2


def accel(r: np.ndarray, t: float, sun: dict, moon: dict, mode: str) -> np.ndarray:
    """Kepler + J2 (mode-dependent) + third body (shared source with campaign)."""
    rm = float(np.linalg.norm(r))
    a = -MU_EARTH_KM3S2 * r / rm ** 3
    if mode in ("j2_only", "sun_moon_j2"):
        z2 = (r[2] * r[2]) / (rm * rm)
        c = -1.5 * J2_EARTH * MU_EARTH_KM3S2 * RE2 / rm ** 5
        g = 1.0 - 5.0 * z2
        a = a + c * np.array([r[0] * g, r[1] * g, r[2] * (3.0 - 5.0 * z2)])
    if mode != "kepler_only":
        a = a + ME.third_body_accel(r, t, sun, moon,
                                    include_sun=mode in ("sun_only", "sun_moon", "sun_moon_j2"),
                                    include_moon=mode in ("moon_only", "sun_moon", "sun_moon_j2"))
    return a


def omega_dot_j2_deg_day(a_km: float, inc_rad: float) -> float:
    """Analytic J2 secular node rate (first order, circular)."""
    n = math.sqrt(MU_EARTH_KM3S2 / a_km ** 3) * 86400.0 * 180.0 / math.pi
    return -1.5 * n * J2_EARTH * (R_EARTH_KM / a_km) ** 2 * math.cos(inc_rad)


def elements(r: np.ndarray, v: np.ndarray) -> tuple:
    """Osculating (a_km, e, i_rad)."""
    rm = float(np.linalg.norm(r))
    vm2 = float(np.dot(v, v))
    a = 1.0 / (2.0 / rm - vm2 / MU_EARTH_KM3S2)
    h = np.cross(r, v)
    hn = float(np.linalg.norm(h))
    inc = math.acos(max(-1.0, min(1.0, h[2] / hn)))
    ev = np.cross(v, h) / MU_EARTH_KM3S2 - r / rm
    return a, float(np.linalg.norm(ev)), inc


def propagate_record(sun: dict, moon: dict, x0: np.ndarray, mode: str,
                     t0_s: float, t_end_s: float, dt_s: float = 60.0) -> dict:
    """Streaming RK4; record t, a, e, i, Omega at every ascending-node crossing."""
    t_cross, om_cross, a_s, e_s, i_s = [], [], [], [], []
    x, t = np.asarray(x0, dtype=float).copy(), t0_s
    z_prev = x[2]
    h0 = np.cross(x[:3], x[3:])
    om_prev = math.atan2(-h0[0], h0[1])
    while t < t_end_s:
        x_new = rk4_step(lambda tt, xx: np.concatenate(
            [xx[3:], accel(xx[:3], tt, sun, moon, mode)]), t, x, dt_s)
        z_curr = x_new[2]
        if z_prev <= 0.0 < z_curr:
            frac = -z_prev / (z_curr - z_prev)
            rc = x[:3] + frac * (x_new[:3] - x[:3])
            vc = x[3:] + frac * (x_new[3:] - x[3:])
            t_cross.append(t + frac * dt_s)
            h = np.cross(rc, vc)
            o = math.atan2(-h[0], h[1])
            while o < om_prev - math.pi:
                o += 2 * math.pi
            while o > om_prev + math.pi:
                o -= 2 * math.pi
            om_prev = o
            om_cross.append(o)
            a_osc, e_osc, i_osc = elements(rc, vc)
            a_s.append(a_osc)
            e_s.append(e_osc)
            i_s.append(i_osc)
        z_prev, x, t = z_curr, x_new, t + dt_s
    om = np.unwrap(np.array(om_cross)) if len(om_cross) > 1 else np.array(om_cross)
    return {"t": np.array(t_cross), "om": om, "a": np.array(a_s),
            "e": np.array(e_s), "i": np.array(i_s), "n_cross": len(t_cross)}


def rolling_mean(y: np.ndarray, win: int) -> np.ndarray:
    """Causal-free boxcar (centred, 'same' length) used to approximate mean elements."""
    if win < 3 or y.size < win:
        return y.copy()
    k = np.ones(win) / win
    pad = win // 2
    yp = np.concatenate([np.full(pad, y[0]), y, np.full(win - 1 - pad, y[-1])])
    return np.convolve(yp, k, mode="valid")


def job(arg) -> dict:
    mode, inc, years = arg
    sun = ME.load_snapshot(ME.SUN_SNAPSHOT)
    moon = ME.load_snapshot(ME.MOON_SNAPSHOT)
    x0 = ME.circular_ic(R_EARTH_KM + H_KM, math.radians(inc), T0_S)
    rec = propagate_record(sun, moon, x0, mode, T0_S, T0_S + years * 365.25 * 86400.0)
    rate = ME.rate_deg_day(rec["t"], rec["om"]) if rec["n_cross"] > 10 else float("nan")
    oj2 = [omega_dot_j2_deg_day(aa, ii) for aa, ii in zip(rec["a"], rec["i"])]
    # 30-day boxcar => proxy for mean elements (removes the ~1-rev osculating wobble)
    win = max(3, int(round(30.0 * max(1, rec["n_cross"]) / max(years * 365.25, 1.0))))
    a_sm = rolling_mean(rec["a"], win)
    i_sm = rolling_mean(rec["i"], win)
    oj2_sm = [omega_dot_j2_deg_day(aa, ii) for aa, ii in zip(a_sm, i_sm)]
    return {"mode": mode, "inc": inc, "years": years, "rate": rate,
            "n_cross": rec["n_cross"],
            "a_mean": float(np.mean(rec["a"])), "a_first_km": float(rec["a"][0]),
            "i_mean_deg": float(np.degrees(np.mean(rec["i"]))),
            "i_first_deg": float(np.degrees(rec["i"][0])),
            "e_mean": float(np.mean(rec["e"])),
            "oj2_mean": float(np.mean(oj2)),
            "oj2_first": float(oj2[0]),
            "oj2_smoothed_mean": float(np.mean(oj2_sm)),
            "a_smoothed_first": float(a_sm[0]), "a_smoothed_last": float(a_sm[-1]),
            "i_smoothed_first_deg": float(np.degrees(i_sm[0])),
            "i_smoothed_last_deg": float(np.degrees(i_sm[-1])),
            "smoothing_window_crossings": int(win)}


def summary_for(by: dict, inc: float, modes: list) -> dict:
    """Per-inclination summary from whatever modes are present.

    With all five modes this reproduces the 1-yr schema exactly. With a reduced
    mode set (the 18.6-yr arm runs j2_only + sun_moon_j2 only, ~45 min/propagation)
    the isolated/residual terms are OMItTED rather than guessed: R needs the
    kepler/sun/moon isolation, which the long arc does not re-run.
    """
    def g(m):
        return by[(m, inc)]
    j2, full = g("j2_only"), g("sun_moon_j2")
    combined = full["rate"] - j2["rate"]
    # Baseline breathing: J2 rate evaluated with each run's own elements.
    D = full["oj2_mean"] - j2["oj2_mean"]
    D_sm = full["oj2_smoothed_mean"] - j2["oj2_smoothed_mean"]
    s = {
        "combined": combined,
        "D_baseline_breathing": D,
        "D_smoothed": D_sm,
        "combined_corrected_smoothed": combined - D_sm,
        "delta_i_mean_deg": full["i_mean_deg"] - j2["i_mean_deg"],
        "delta_a_mean_km": full["a_mean"] - j2["a_mean"],
        "full_smoothed_di_deg_over_arc": (full["i_smoothed_last_deg"]
                                          - full["i_smoothed_first_deg"]),
        "full_smoothed_da_km_over_arc": full["a_smoothed_last"] - full["a_smoothed_first"],
        "sensitivity_check_deg_day": 0.12579 * (full["i_smoothed_last_deg"]
                                                - full["i_smoothed_first_deg"]),
    }
    if {"kepler_only", "sun_only", "moon_only"} <= set(modes):
        kep, su, mo = g("kepler_only"), g("sun_only"), g("moon_only")
        isolated = (su["rate"] - kep["rate"]) + (mo["rate"] - kep["rate"])
        R = combined - isolated
        s["R_measured"] = R
        s["isolated"] = isolated
        s["fraction_of_R_explained_raw"] = (D / R) if R not in (0.0, None) else None
        s["fraction_of_R_explained_smoothed"] = (D_sm / R) if R not in (0.0, None) else None
        s["implied_di_for_R_deg"] = (R / 0.12579) if R not in (0.0, None) else None
    return s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--modes", default=",".join(MODES),
                    help="comma-separated subset of: %s" % ",".join(MODES))
    args = ap.parse_args()
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    unknown = [m for m in modes if m not in MODES]
    if unknown:
        ap.error("unknown mode(s): %s" % ", ".join(unknown))
    missing = [m for m in ("j2_only", "sun_moon_j2") if m not in modes]
    if missing:
        ap.error("j2_only and sun_moon_j2 are required; they define D and combined")
    (HERE / "results").mkdir(exist_ok=True)
    tasks = [(m, inc, args.years) for inc in INCS for m in modes]
    with Pool(args.workers) as p:
        rows = p.map(job, tasks)
    by = {(r["mode"], r["inc"]): r for r in rows}
    out = {"years": args.years, "modes": modes, "rows": rows, "summary": {}}
    for inc in INCS:
        out["summary"][str(inc)] = summary_for(by, inc, modes)
    out["provenance"] = {
        "sun_sha256": ME.load_snapshot(ME.SUN_SNAPSHOT)["sha256"],
        "moon_sha256": ME.load_snapshot(ME.MOON_SNAPSHOT)["sha256"],
        "code": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                 for p in sorted(HERE.glob("*.py"))},
    }
    name = "phase0_%gyr.json" % args.years
    out["mode_set_complete"] = set(modes) == set(MODES)
    (HERE / "results" / name).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out["summary"], indent=2))


if __name__ == "__main__":
    main()
