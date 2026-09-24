"""Final 18.6-yr validation (correct interp, 2 modes x 3 inc, dt=60s).

Writes results/final_18yr.json. ~3.5 hr single-core, ~30 min on 7 workers.
Only run after headline/robustness decision rule is recorded.
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
W_D = 6798.4  # one lunar nodal cycle
MODES = ["j2_only", "sun_moon_j2"]
INCS = [97.7876, 90.0, 30.0]


def job(arg):
    mode, inc = arg
    sun = ME.load_snapshot(ME.SUN_SNAPSHOT)
    moon = ME.load_snapshot(ME.MOON_SNAPSHOT)
    a = 6378.137 + 600.0
    x0 = ME.circular_ic(a, math.radians(inc), T0_S)
    out = ME.propagate(sun, moon, x0, mode=mode, t0_s=T0_S,
                       t_end_s=T0_S + W_D * 86400.0, dt_s=60.0)
    import numpy as np
    r = ME.rate_deg_day(np.array(out["t_cross"]), np.array(out["om_cross"]))
    return {"mode": mode, "inc": inc, "rate": r, "n": len(out["t_cross"])}


def main():
    t0 = time.monotonic()
    tasks = [(m, inc) for inc in INCS for m in MODES]
    # Pool(6): 6 jobs fit in one wave; BLAS threads pinned to 1 so commit
    # headroom is ~6x66 MB (see run_headline note on WinError 1455).
    with Pool(6) as p:
        rows = p.map(job, tasks)
    by = {(r["mode"], r["inc"]): r["rate"] for r in rows}
    luni = {str(inc): by[("sun_moon_j2", inc)] - by[("j2_only", inc)] for inc in INCS}
    out = {"rows": rows, "lunisolar": luni,
           "sun_sha": ME.load_snapshot(ME.SUN_SNAPSHOT)["sha256"],
           "moon_sha": ME.load_snapshot(ME.MOON_SNAPSHOT)["sha256"],
           "code": ME.code_hashes()}
    (HERE / "results" / "final_18yr.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
