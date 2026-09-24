"""Re-run frozen-plane set with open-loop correction (365-d, 3 inc x 6 modes)."""
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


def job(arg):
    mode, inc = arg
    sun = ME.load_snapshot(ME.SUN_SNAPSHOT)
    moon = ME.load_snapshot(ME.MOON_SNAPSHOT)
    a = 6378.137 + 600.0
    x0 = ME.circular_ic(a, math.radians(inc), T0_S)
    out = ME.propagate(sun, moon, x0, mode=mode, t0_s=T0_S, t_end_s=T0_S + W_D * 86400.0,
                       dt_s=60.0, freeze_plane=True)
    import numpy as np
    r = ME.rate_deg_day(np.array(out["t_cross"]), np.array(out["om_cross"])) if len(out["t_cross"]) > 10 else float("nan")
    # also raw j2/full for diagnostics
    return {"mode": mode, "inc": inc, "rate": r}


def main():
    t0 = time.monotonic()
    tasks = [(m, inc) for inc in INCS for m in MODES]
    # Pool(4): commit-charge headroom (see run_headline note).
    with Pool(4) as p:
        rows = p.map(job, tasks)
    out = {}
    for inc in INCS:
        by = {r["mode"]: r["rate"] for r in rows if r["inc"] == inc}
        kep = by["kepler_only"]
        c = by["sun_moon_j2"] - by["j2_only"]
        iso = (by["sun_only"] - kep) + (by["moon_only"] - kep)
        out[str(inc)] = {"combined": c, "isolated": iso, "R": c - iso, "j2_only": by["j2_only"],
                         "full": by["sun_moon_j2"]}
    (HERE / "results" / "frozen_365d.json").write_text(json.dumps({"rows": rows, "summary": out}, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
