"""Run the perturbative scaling experiment for mission_j2_lunisolar_coupling.

For (lambda_J2, lambda_3body) pairs, propagate for 90 d at h=600 km i=i_sso
and measure the Lunisolar RAAN rate (full minus J2-only, with both modes
scaled by their respective multipliers).

The 2-D polynomial response is fit to:
    f(lambda_J2, lambda_3body) = a10 * lambda_J2 + a01 * lambda_3body
                                + a11 * lambda_J2 * lambda_3body
                                + a20 * lambda_J2^2 + a02 * lambda_3body^2

The cross coefficient a11 is the perturbation-scaling discriminator.

A 90-d arc is used here to keep the compute modest; the scaling experiment
is about the algebraic structure of the response, not the long-term secular
behavior.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent.parent))  # for lab_utils

from mission_experiment import (
    H_SSO_KM,
    I_SSO_DEG,
    DT_S,
    MU_EARTH_KM3S2,
    R_EARTH_KM,
    SUN_SNAPSHOT,
    MOON_SNAPSHOT,
    _load_snapshot,
    harmonic_regression,
    ols_slope,
    propagate_streaming_with_x0,
    code_hashes,
)

DAYS_90_S = 90.0 * 86400.0
DAYS_365_S = 365.25 * 86400.0

# Perturbative scaling values (lambda_J2, lambda_3body)
LAMBDA_J2_VALUES = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
LAMBDA_3BODY_VALUES = [0.0, 0.5, 1.0, 2.0]


def initial_state(h_km: float, i_deg: float, mu: float = MU_EARTH_KM3S2,
                  r_eq: float = R_EARTH_KM) -> np.ndarray:
    """Circular-orbit initial state at given altitude and inclination."""
    a = r_eq + h_km
    v_circ = math.sqrt(mu / a)
    i_rad = math.radians(i_deg)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, v_circ * math.cos(i_rad), v_circ * math.sin(i_rad)])
    return np.concatenate([r0, v0])


def fit_polynomial_2d(lam_j2: np.ndarray, lam_3b: np.ndarray,
                       rates: np.ndarray) -> dict:
    """Fit 2-D polynomial response surface."""
    # Design matrix: [1, lj2, l3b, lj2*l3b, lj2^2, l3b^2]
    A = np.column_stack([
        np.ones_like(lam_j2),
        lam_j2,
        lam_3b,
        lam_j2 * lam_3b,
        lam_j2 ** 2,
        lam_3b ** 2,
    ])
    result = np.linalg.lstsq(A, rates, rcond=None)
    coeffs = result[0]
    residuals = rates - A @ coeffs
    rms = float(np.sqrt(np.mean(residuals ** 2)))
    # Uncertainties (assuming uncorrelated noise; use lstsq residuals to estimate)
    if len(rates) > len(coeffs):
        sigma2 = np.sum(residuals ** 2) / (len(rates) - len(coeffs))
        cov = sigma2 * np.linalg.inv(A.T @ A)
        sigma = np.sqrt(np.diag(cov))
    else:
        sigma = np.full_like(coeffs, float("nan"))
    return {
        "a_const": float(coeffs[0]),
        "a10_lambda_j2": float(coeffs[1]),
        "a01_lambda_3body": float(coeffs[2]),
        "a11_cross": float(coeffs[3]),
        "a20_lambda_j2_sq": float(coeffs[4]),
        "a02_lambda_3body_sq": float(coeffs[5]),
        "a_const_sigma": float(sigma[0]),
        "a10_sigma": float(sigma[1]),
        "a01_sigma": float(sigma[2]),
        "a11_sigma": float(sigma[3]),
        "a20_sigma": float(sigma[4]),
        "a02_sigma": float(sigma[5]),
        "rms_residual": rms,
        "n_points": len(rates),
        "snr_cross": abs(coeffs[3]) / sigma[3] if sigma[3] > 0 else float("nan"),
        "snr_lambda_j2_sq": abs(coeffs[4]) / sigma[4] if sigma[4] > 0 else float("nan"),
        "snr_lambda_3body_sq": abs(coeffs[5]) / sigma[5] if sigma[5] > 0 else float("nan"),
    }


def run_perturbative_scaling(arc_days: float = 90.0,
                              i_deg: float = I_SSO_DEG,
                              h_km: float = H_SSO_KM) -> dict:
    """Run the perturbative scaling experiment.

    Returns dict with all (lambda_J2, lambda_3body) RAAN rates and the
    polynomial fit coefficients.
    """
    print(f"Loading snapshots...")
    sun_snap = _load_snapshot(SUN_SNAPSHOT)
    moon_snap = _load_snapshot(MOON_SNAPSHOT)
    print(f"  Sun: {sun_snap['n_points']} points, sha256={sun_snap['sha256'][:16]}...")
    print(f"  Moon: {moon_snap['n_points']} points, sha256={moon_snap['sha256'][:16]}...")

    x0 = initial_state(h_km, i_deg)
    t_end_s = arc_days * 86400.0

    results_grid = []
    lam_j2_arr = []
    lam_3b_arr = []
    rate_full_arr = []
    rate_j2_arr = []

    print(f"\nRunning perturbative scaling sweep at h={h_km} km, i={i_deg} deg, arc={arc_days} d")
    print(f"  Lambda_J2 values: {LAMBDA_J2_VALUES}")
    print(f"  Lambda_3body values: {LAMBDA_3BODY_VALUES}")
    print(f"  Total propagations: {len(LAMBDA_J2_VALUES) * len(LAMBDA_3BODY_VALUES) * 2}")

    for lj2 in LAMBDA_J2_VALUES:
        for l3b in LAMBDA_3BODY_VALUES:
            t0 = time.time()
            # Full model (with multipliers)
            res_full = propagate_streaming_with_x0(
                sun_snap, moon_snap, x0,
                mode="sun_moon_j2",
                t0_s=0.0, t_end_s=t_end_s, dt_s=DT_S,
                lambda_j2=lj2, lambda_3body=l3b,
            )
            # J2-only (with same lambda_j2 but no 3b)
            res_j2 = propagate_streaming_with_x0(
                sun_snap, moon_snap, x0,
                mode="j2_only",
                t0_s=0.0, t_end_s=t_end_s, dt_s=DT_S,
                lambda_j2=lj2, lambda_3body=0.0,
            )
            t_day_full = (res_full["t_cross"] - 0.0) / 86400.0
            t_day_j2 = (res_j2["t_cross"] - 0.0) / 86400.0
            # OLS slope of the unwrapped RAAN at ascending-node crossings
            _, b_full = ols_slope(res_full["t_cross"], res_full["om_cross"])
            _, b_j2 = ols_slope(res_j2["t_cross"], res_j2["om_cross"])
            rate_full_deg_day = math.degrees(b_full) * 86400.0
            rate_j2_deg_day = math.degrees(b_j2) * 86400.0
            rate_luni = rate_full_deg_day - rate_j2_deg_day
            results_grid.append({
                "lambda_j2": lj2,
                "lambda_3body": l3b,
                "rate_full_deg_day": rate_full_deg_day,
                "rate_j2_deg_day": rate_j2_deg_day,
                "rate_lunisolar_deg_day": rate_luni,
                "n_nodes_full": len(res_full["t_cross"]),
                "n_nodes_j2": len(res_j2["t_cross"]),
                "wall_clock_s": time.time() - t0,
            })
            lam_j2_arr.append(lj2)
            lam_3b_arr.append(l3b)
            rate_full_arr.append(rate_full_deg_day)
            rate_j2_arr.append(rate_j2_deg_day)
            print(f"  lambda_j2={lj2:.2f}, lambda_3b={l3b:.2f}: "
                  f"full={rate_full_deg_day:+.6f}, j2={rate_j2_deg_day:+.6f}, "
                  f"luni={rate_luni:+.6e} deg/day  ({time.time()-t0:.1f}s)")

    # Fit the 2-D polynomial to the FULL rate (J2 + 3b response surface)
    lam_j2_arr = np.array(lam_j2_arr)
    lam_3b_arr = np.array(lam_3b_arr)
    rate_full_arr = np.array(rate_full_arr)
    rate_j2_arr = np.array(rate_j2_arr)

    fit_full = fit_polynomial_2d(lam_j2_arr, lam_3b_arr, rate_full_arr)
    fit_j2 = fit_polynomial_2d(lam_j2_arr, lam_3b_arr, rate_j2_arr)

    print(f"\nPolynomial fit to FULL (J2+3b) rate:")
    print(f"  a11 (cross)    = {fit_full['a11_cross']:+.4e} +/- {fit_full['a11_sigma']:.2e}  "
          f"SNR={fit_full['snr_cross']:.2f}")
    print(f"  a20 (J2^2)     = {fit_full['a20_lambda_j2_sq']:+.4e} +/- {fit_full['a20_sigma']:.2e}  "
          f"SNR={fit_full['snr_lambda_j2_sq']:.2f}")
    print(f"  a02 (3b^2)     = {fit_full['a02_lambda_3body_sq']:+.4e} +/- {fit_full['a02_sigma']:.2e}  "
          f"SNR={fit_full['snr_lambda_3body_sq']:.2f}")
    print(f"  RMS residual   = {fit_full['rms_residual']:.4e} deg/day")

    return {
        "config": {
            "arc_days": arc_days,
            "h_km": h_km,
            "i_deg": i_deg,
            "lambda_j2_values": LAMBDA_J2_VALUES,
            "lambda_3body_values": LAMBDA_3BODY_VALUES,
        },
        "grid_results": results_grid,
        "polynomial_fit_full": fit_full,
        "polynomial_fit_j2_only": fit_j2,
        "snapshot_provenance": {
            "sun_sha256": sun_snap["sha256"],
            "moon_sha256": moon_snap["sha256"],
            "sun_n_points": sun_snap["n_points"],
            "moon_n_points": moon_snap["n_points"],
        },
        "code_hashes": code_hashes(),
    }


if __name__ == "__main__":
    out = run_perturbative_scaling(arc_days=90.0)
    out_path = HERE / "results" / "perturbative_scaling.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nResults written to {out_path}")
