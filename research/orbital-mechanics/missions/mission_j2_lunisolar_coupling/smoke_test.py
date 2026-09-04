__test__ = False  # not a pytest test module
"""Quick smoke test for the streaming RK4 propagator with mode isolation."""
import math
import sys
import time

import numpy as np

sys.path.insert(0, '.')

from mission_experiment import (
    H_SSO_KM,
    I_SSO_DEG,
    DT_S,
    MU_EARTH_KM3S2,
    R_EARTH_KM,
    SUN_SNAPSHOT,
    MOON_SNAPSHOT,
    PARENT_MANIFEST,
    _load_snapshot,
    ols_slope,
    propagate_streaming_with_x0,
    synthetic_circular_moon_state,
)


def initial_state(h_km, i_deg, mu=MU_EARTH_KM3S2, r_eq=R_EARTH_KM):
    a = r_eq + h_km
    v_circ = math.sqrt(mu / a)
    i_rad = math.radians(i_deg)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, v_circ * math.cos(i_rad), v_circ * math.sin(i_rad)])
    return np.concatenate([r0, v0])


def main():
    print("Loading snapshots...")
    t0 = time.time()
    sun_snap = _load_snapshot(SUN_SNAPSHOT)
    moon_snap = _load_snapshot(MOON_SNAPSHOT)
    print(f"  Sun: {sun_snap['n_points']} points, sha256={sun_snap['sha256'][:16]}...")
    print(f"  Moon: {moon_snap['n_points']} points, sha256={moon_snap['sha256'][:16]}...")
    print(f"  Loaded in {time.time()-t0:.2f}s")

    x0 = initial_state(H_SSO_KM, I_SSO_DEG)
    t_end_s = 90.0 * 86400.0  # 90-day smoke test

    print(f"\nSmoke test: 90 d at h=600 km i=i_sso={I_SSO_DEG} deg, dt={DT_S}s")
    print(f"  Initial state: r=({x0[0]:.1f}, {x0[1]:.1f}, {x0[2]:.1f}) km, "
          f"v=({x0[3]:.3f}, {x0[4]:.3f}, {x0[5]:.3f}) km/s")

    for mode in ["kepler_only", "j2_only", "sun_only", "moon_only",
                  "sun_moon", "sun_moon_j2"]:
        t1 = time.time()
        res = propagate_streaming_with_x0(
            sun_snap, moon_snap, x0,
            mode=mode, t0_s=0.0, t_end_s=t_end_s, dt_s=DT_S,
        )
        n_steps = res["n_steps"]
        n_nodes = len(res["t_cross"])
        if n_nodes > 4:
            _, b = ols_slope(res["t_cross"], res["om_cross"])
            rate_deg_day = math.degrees(b) * 86400.0
            print(f"  {mode:15s}: {n_steps} steps, {n_nodes} nodes, "
                  f"slope={rate_deg_day:+.6e} deg/day  ({time.time()-t1:.1f}s)")
        else:
            print(f"  {mode:15s}: {n_steps} steps, {n_nodes} nodes  ({time.time()-t1:.1f}s)")

    print("\nDone.")


if __name__ == "__main__":
    main()
