"""Perturbative scaling + symmetry + dt ladder (90-d, i_sso, h=600).

Repeats Mission-2 methodology with correct interp, adds small-scale
symmetry set and dt ladder for detection-limit + stability.
Writes results/scaling_90d.json.
"""
import json
import math
import sys
import time
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mission_experiment as ME

T0_S = 820476800.0
W_D = 90.0
LJ2 = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
L3B = [0.0, 0.5, 1.0, 2.0]
SYM = [-0.5, -0.25, 0.25, 0.5]


def job(arg):
    lj2, l3b, dt, tag = arg
    sun = ME.load_snapshot(ME.SUN_SNAPSHOT)
    moon = ME.load_snapshot(ME.MOON_SNAPSHOT)
    a = 6378.137 + 600.0
    x0 = ME.circular_ic(a, math.radians(97.7876), T0_S)
    # full with multipliers
    out = ME.propagate(sun, moon, x0, mode="sun_moon_j2", t0_s=T0_S,
                       t_end_s=T0_S + W_D * 86400.0, dt_s=dt,
                       lambda_j2=lj2, lambda_3body=l3b)
    r_full = ME.rate_deg_day(np.array(out["t_cross"]), np.array(out["om_cross"]))
    # j2-only at same lj2 for subtraction (lunisolar = full - j2only)
    out2 = ME.propagate(sun, moon, x0, mode="j2_only", t0_s=T0_S,
                        t_end_s=T0_S + W_D * 86400.0, dt_s=dt, lambda_j2=lj2)
    r_j2 = ME.rate_deg_day(np.array(out2["t_cross"]), np.array(out2["om_cross"]))
    return {"lj2": lj2, "l3b": l3b, "dt": dt, "tag": tag,
            "full": r_full, "j2": r_j2, "luni": r_full - r_j2}


def fit_a11(rows):
    # fit luni = a01*l3b + a11*lj2*l3b + a02*l3b^2 (+a10*lj2 absorbed? luni(0)=0 so no lj2-only)
    # use full 2D poly like prior: f = a10 lj2 + a01 l3b + a11 lj2 l3b + a20 lj2^2 + a02 l3b^2 on full rate
    X, y = [], []
    for r in rows:
        X.append([r["lj2"], r["l3b"], r["lj2"] * r["l3b"], r["lj2"] ** 2, r["l3b"] ** 2])
        y.append(r["full"])
    X = np.array(X)
    y = np.array(y)
    coef, res, _, _ = np.linalg.lstsq(np.column_stack([np.ones(len(y)), X]), y, rcond=None)
    # coef[0]=const, [1]=a10,[2]=a01,[3]=a11,[4]=a20,[5]=a02
    pred = np.column_stack([np.ones(len(y)), X]) @ coef
    dof = max(len(y) - 6, 1)
    sigma2 = float(np.sum((y - pred) ** 2) / dof)
    cov = sigma2 * np.linalg.inv(np.column_stack([np.ones(len(y)), X]).T @ np.column_stack([np.ones(len(y)), X]))
    se = np.sqrt(np.diag(cov))
    return {"a11": float(coef[3]), "se11": float(se[3]), "snr": float(abs(coef[3]) / se[3]),
            "coef": [float(c) for c in coef], "se": [float(s) for s in se]}


def main():
    t0 = time.monotonic()
    grid = [(lj, lb, 60.0, "grid") for lj in LJ2 for lb in L3B]
    sym = [(s, 0.5, 60.0, "sym") for s in SYM] + [(0.5, s, 60.0, "sym") for s in SYM]
    ladder = [(1.0, 1.0, dt, "ladder") for dt in (120.0, 60.0, 30.0)] + \
             [(1.0, 0.0, dt, "ladder") for dt in (120.0, 60.0, 30.0)]
    tasks = grid + sym + ladder
    # Pool(4): commit-charge headroom (see run_headline note).
    with Pool(4) as p:
        rows = p.map(job, tasks)
    grid_rows = [r for r in rows if r["tag"] == "grid"]
    fit = fit_a11(grid_rows)
    # detection limit from lj2=0 or l3b=0 scatter
    zeros = [r["luni"] for r in grid_rows if r["lj2"] == 0 or r["l3b"] == 0]
    det = float(np.std(zeros)) if len(zeros) > 2 else float("nan")
    out = {"rows": rows, "fit": fit, "detection_limit": det,
           "sun_sha256": ME.load_snapshot(ME.SUN_SNAPSHOT)["sha256"],
           "moon_sha256": ME.load_snapshot(ME.MOON_SNAPSHOT)["sha256"],
           "code": ME.code_hashes()}
    (HERE / "results" / "scaling_90d.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({"fit": fit, "det": det}, indent=2))


if __name__ == "__main__":
    main()
