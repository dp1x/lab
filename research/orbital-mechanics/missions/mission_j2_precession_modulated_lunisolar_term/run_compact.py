"""Compact 90-d reproduction: 6 modes x {correct,buggy} interp at h=600 i_sso.

Quantifies the shared-lineage interpolation bug impact and reproduces R
before any headline interpretation. Writes results/compact_interp_check.json.
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

H_KM = 600.0
I_DEG = 97.7876
T0_S = (820476800.0)
W_D = 90.0
MODES = ["kepler_only", "j2_only", "sun_only", "moon_only", "sun_moon", "sun_moon_j2"]


def job(args):
    mode, buggy = args
    sun = ME.load_snapshot(ME.SUN_SNAPSHOT)
    moon = ME.load_snapshot(ME.MOON_SNAPSHOT)
    a = 6378.137 + H_KM
    x0 = ME.circular_ic(a, math.radians(I_DEG), T0_S)
    out = ME.propagate(sun, moon, x0, mode=mode, t0_s=T0_S, t_end_s=T0_S + W_D * 86400.0,
                       dt_s=60.0, use_buggy_interp=buggy)
    import numpy as np
    r = ME.rate_deg_day(np.array(out["t_cross"]), np.array(out["om_cross"])) if len(out["t_cross"]) > 10 else float("nan")
    return {"mode": mode, "buggy": buggy, "rate": r, "n_cross": len(out["t_cross"])}


def main():
    t0 = time.monotonic()
    tasks = [(m, b) for m in MODES for b in (False, True)]
    # Pool(4): commit-charge headroom (see run_headline note).
    with Pool(4) as p:
        rows = p.map(job, tasks)
    by = {(r["mode"], r["buggy"]): r["rate"] for r in rows}

    def R(buggy):
        f = by[("sun_moon_j2", buggy)] - by[("j2_only", buggy)]
        iso = by[("sun_only", buggy)] + by[("moon_only", buggy)] - by[("kepler_only", buggy)]
        # R = combined - isolated ; combined=(full-j2), isolated=(sun+moon-kepler - kepler?) use standard:
        # full-j2 = luni_combined; sun+moon-2*kepler? kepler rate ~0; use (sun-kepler)+(moon-kepler)
        kep = by[("kepler_only", buggy)]
        iso2 = (by[("sun_only", buggy)] - kep) + (by[("moon_only", buggy)] - kep)
        return f, iso2, f - iso2

    f_c, iso_c, r_c = R(False)
    f_b, iso_b, r_b = R(True)
    out = {"rows": rows, "correct": {"combined": f_c, "isolated": iso_c, "R": r_c},
           "buggy": {"combined": f_b, "isolated": iso_b, "R": r_b},
           "sun_sha": ME.load_snapshot(ME.SUN_SNAPSHOT)["sha256"],
           "moon_sha": ME.load_snapshot(ME.MOON_SNAPSHOT)["sha256"]}
    (HERE / "results" / "compact_interp_check.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
