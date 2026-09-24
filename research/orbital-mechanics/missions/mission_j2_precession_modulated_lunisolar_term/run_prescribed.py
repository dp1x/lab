"""Prescribed-precession sweep (90-d, i_sso, 3body-only + kinematic rotation).

Tests H-mod frequency dependence: R vs Omega_dot_p at fixed 3body.
Writes results/prescribed_90d.json. 7 rates x 1 mode (+J2-only baseline).
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
W_D = 90.0
RATES = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]  # x OJ2(600,SSO)=0.9856


def job(om_p):
    sun = ME.load_snapshot(ME.SUN_SNAPSHOT)
    moon = ME.load_snapshot(ME.MOON_SNAPSHOT)
    a = 6378.137 + 600.0
    x0 = ME.circular_ic(a, math.radians(97.7876), T0_S)
    # 3body-only with kinematic precession (no J2 force)
    out = ME.propagate(sun, moon, x0, mode="sun_moon", t0_s=T0_S,
                       t_end_s=T0_S + W_D * 86400.0, dt_s=60.0,
                       omega_p_deg_day=om_p)
    import numpy as np
    r = ME.rate_deg_day(np.array(out["t_cross"]), np.array(out["om_cross"]))
    # baseline stationary (om_p=0 already in sweep); also J2-only rate for reference
    return {"omega_p": om_p, "rate": r}


def response_minus_command(row):
    """Return the physical response after removing the imposed node-rate command."""
    response = row["rate"] - row["omega_p"]
    return {**row, "response": response}


def main():
    # Pool(4): commit-charge headroom (see run_headline note).
    with Pool(4) as p:
        rows = p.map(job, RATES)
    corrected = [response_minus_command(r) for r in rows]
    # The response is compared against the stationary third-body-only baseline.
    base = [r for r in corrected if r["omega_p"] == 0.0][0]["response"]
    for r in corrected:
        r["R_vs_stationary"] = r["response"] - base
    out = {"rows": corrected, "baseline_stationary": base,
           "sun_sha": ME.load_snapshot(ME.SUN_SNAPSHOT)["sha256"],
           "moon_sha": ME.load_snapshot(ME.MOON_SNAPSHOT)["sha256"],
           "code": ME.code_hashes()}
    (HERE / "results" / "prescribed_90d.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
