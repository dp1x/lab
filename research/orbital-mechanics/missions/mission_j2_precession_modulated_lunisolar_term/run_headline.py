"""Headline 365-d force-mode campaign (correct interp) + frozen-plane + altitude.

Writes results/headline_365d.json, results/frozen_365d.json, results/altitude_365d.json.
Parallelized (7 workers). Streaming, no full-trajectory storage.
"""
import json
import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "src")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mission_experiment as ME

T0_S = 820476800.0
W_D = 365.25
MODES = ["kepler_only", "j2_only", "sun_only", "moon_only", "sun_moon", "sun_moon_j2"]
INCS = [97.7876, 90.0, 30.0]
ALTS = [500.0, 800.0]


def job(arg):
    mode, inc, h, frozen, label = arg
    sun = ME.load_snapshot(ME.SUN_SNAPSHOT)
    moon = ME.load_snapshot(ME.MOON_SNAPSHOT)
    a = 6378.137 + h
    x0 = ME.circular_ic(a, math.radians(inc), T0_S)
    out = ME.propagate(sun, moon, x0, mode=mode, t0_s=T0_S, t_end_s=T0_S + W_D * 86400.0,
                       dt_s=60.0, freeze_plane=frozen)
    import numpy as np
    r = ME.rate_deg_day(np.array(out["t_cross"]), np.array(out["om_cross"])) if len(out["t_cross"]) > 10 else float("nan")
    return {"mode": mode, "inc": inc, "h": h, "frozen": frozen, "rate": r, "n": len(out["t_cross"])}


def run_set(tasks):
    # Pool(4): host commit-charge headroom is limited (no pagefile margin);
    # 7 workers + OpenBLAS threading caused WinError 1455 on spawn.
    with Pool(4) as p:
        return p.map(job, tasks)


def R_of(rows):
    by = {(r["mode"]): r["rate"] for r in rows}
    kep = by["kepler_only"]
    combined = by["sun_moon_j2"] - by["j2_only"]
    isolated = (by["sun_only"] - kep) + (by["moon_only"] - kep)
    # Gate-(f) residual used by the frozen regression test: summary {combined,
    # isolated, R}. j2_only/full rows stay in the per-mode rows only.
    return combined, isolated, combined - isolated


def summary_with_provenance(extra=None):
    prov = {"sun_sha256": ME.load_snapshot(ME.SUN_SNAPSHOT)["sha256"],
            "moon_sha256": ME.load_snapshot(ME.MOON_SNAPSHOT)["sha256"],
            "code": ME.code_hashes()}
    if extra:
        prov.update(extra)
    return prov


def main():
    t0 = time.monotonic()
    # headline: 3 inc x 6 modes at h=600
    tasks = [(m, inc, 600.0, False, "head") for inc in INCS for m in MODES]
    rows = run_set(tasks)
    head = {}
    for inc in INCS:
        sub = [r for r in rows if r["inc"] == inc]
        c, iso, R = R_of(sub)
        head[str(inc)] = {"combined": c, "isolated": iso, "R": R, "rows": sub}
    (HERE / "results" / "headline_365d.json").write_text(json.dumps(
        {"summary": {k: {kk: v[kk] for kk in ("combined", "isolated", "R")} for k, v in head.items()},
         "rows": rows, "provenance": summary_with_provenance({"W_d": W_D, "h_km": 600.0})}, indent=2))
    print("HEADLINE", json.dumps({k: {kk: v[kk] for kk in ("combined", "isolated", "R")} for k, v in head.items()}, indent=2))
    # frozen at same grid
    tasks2 = [(m, inc, 600.0, True, "froz") for inc in INCS for m in MODES]
    rows2 = run_set(tasks2)
    froz = {}
    for inc in INCS:
        sub = [r for r in rows2 if r["inc"] == inc]
        c, iso, R = R_of(sub)
        froz[str(inc)] = {"combined": c, "isolated": iso, "R": R}
    (HERE / "results" / "frozen_365d.json").write_text(json.dumps(
        {"summary": froz, "rows": rows2,
         "provenance": summary_with_provenance({"W_d": W_D, "h_km": 600.0, "control": "freeze_plane open-loop"})}, indent=2))
    print("FROZEN", json.dumps(froz, indent=2))
    # altitude at i_sso
    tasks3 = [(m, 97.7876, h, False, "alt") for h in ALTS for m in MODES]
    rows3 = run_set(tasks3)
    alt = {}
    for h in ALTS:
        sub = [r for r in rows3 if r["h"] == h]
        c, iso, R = R_of(sub)
        alt[str(h)] = {"combined": c, "isolated": iso, "R": R}
    (HERE / "results" / "altitude_365d.json").write_text(json.dumps(
        {"summary": alt, "rows": rows3,
         "provenance": summary_with_provenance({"W_d": W_D, "inc_deg": 97.7876})}, indent=2))
    print("ALT", json.dumps(alt, indent=2))


if __name__ == "__main__":
    main()
