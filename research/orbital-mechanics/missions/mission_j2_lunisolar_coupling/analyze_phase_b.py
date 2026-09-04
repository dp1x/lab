__test__ = False  # not a pytest test module
"""Analyze Phase B perturbative scaling results: fit 2D polynomial, check for cross term."""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main():
    with open(HERE / "results" / "phase_b_perturbative_scaling.json") as f:
        data = json.load(f)

    results = data["results"]
    # Group: for each (lj2, l3b) pair, get the FULL rate and the J2-only rate
    # Each pair has one full and one j2_only job
    full_rates = {}  # key: (lj2, l3b) -> rate
    j2_rates = {}    # key: (lj2) -> rate
    for r in results:
        lj2 = r["lambda_j2"]
        l3b = r["lambda_3body"]
        if r["mode"] == "sun_moon_j2":
            full_rates[(lj2, l3b)] = r["rate_deg_day"]
        elif r["mode"] == "j2_only":
            j2_rates[lj2] = r["rate_deg_day"]

    print("=== Phase B perturbative scaling analysis ===")
    print(f"Lambda_J2 values: {data['lambda_j2_values']}")
    print(f"Lambda_3body values: {data['lambda_3body_values']}")
    print()

    # Compute Lunisolar rate = full - J2-only at each (lj2, l3b)
    # Note: the full propagations use lambda_j2=lj2, lambda_3body=l3b
    # The J2-only propagations use lambda_j2=lj2, lambda_3body=0
    # So Lunisolar = full(lj2, l3b) - j2_only(lj2)
    print("Lunisolar rate (full - J2-only) at each (lambda_J2, lambda_3body):")
    print("  lambda_J2 \\ lambda_3body | ", end="")
    for l3b in data["lambda_3body_values"]:
        print(f"  {l3b:.2f}    ", end="")
    print()
    for lj2 in data["lambda_j2_values"]:
        print(f"  {lj2:.2f}                    | ", end="")
        for l3b in data["lambda_3body_values"]:
            luni = full_rates[(lj2, l3b)] - j2_rates[lj2]
            print(f"  {luni:+.4e}", end="")
        print()

    print()
    print("=== Polynomial fit to Lunisolar rate ===")
    print("Fit: f(lj2, l3b) = a10*lj2 + a01*l3b + a11*lj2*l3b + a20*lj2^2 + a02*l3b^2")
    lj2_arr = []
    l3b_arr = []
    luni_arr = []
    for lj2 in data["lambda_j2_values"]:
        for l3b in data["lambda_3body_values"]:
            lj2_arr.append(lj2)
            l3b_arr.append(l3b)
            luni_arr.append(full_rates[(lj2, l3b)] - j2_rates[lj2])
    lj2_arr = np.array(lj2_arr)
    l3b_arr = np.array(l3b_arr)
    luni_arr = np.array(luni_arr)

    A = np.column_stack([
        np.ones_like(lj2_arr), lj2_arr, l3b_arr,
        lj2_arr * l3b_arr, lj2_arr ** 2, l3b_arr ** 2,
    ])
    coeffs, residuals, rank, sv = np.linalg.lstsq(A, luni_arr, rcond=None)
    rms = np.sqrt(np.mean((luni_arr - A @ coeffs) ** 2))

    # Uncertainty: assume noise in residuals
    n_data = len(luni_arr)
    n_coeffs = len(coeffs)
    if n_data > n_coeffs:
        sigma2 = np.sum((luni_arr - A @ coeffs) ** 2) / (n_data - n_coeffs)
        cov = sigma2 * np.linalg.inv(A.T @ A)
        sigma = np.sqrt(np.diag(cov))
    else:
        sigma = np.full_like(coeffs, float("nan"))

    names = ["a_const", "a10_lj2", "a01_l3b", "a11_cross", "a20_lj2^2", "a02_l3b^2"]
    print()
    print(f"  {'name':18s} {'value':>15s} {'sigma':>15s} {'SNR':>8s}")
    for i, name in enumerate(names):
        snr = abs(coeffs[i]) / sigma[i] if sigma[i] > 0 else float("nan")
        print(f"  {name:18s} {coeffs[i]:+15.6e} {sigma[i]:15.6e} {snr:8.2f}")
    print(f"  RMS residual: {rms:.4e} deg/day")

    print()
    print("=== Verdict on cross-term a11 ===")
    a11 = coeffs[3]
    a11_sigma = sigma[3]
    snr_a11 = abs(a11) / a11_sigma if a11_sigma > 0 else float("nan")
    if snr_a11 > 3.0:
        print(f"  Cross term a11 = {a11:+.4e} +/- {a11_sigma:.4e}, SNR = {snr_a11:.2f}")
        print(f"  -> Cross term IS statistically significant (>3σ)")
        print(f"  -> The J2 x Lunisolar coupling MAY be present at this arc length")
    else:
        print(f"  Cross term a11 = {a11:+.4e} +/- {a11_sigma:.4e}, SNR = {snr_a11:.2f}")
        print(f"  -> Cross term NOT significant at this arc length")
        print(f"  -> Note: 90-day arc may be too short to detect a small cross term")

    print()
    print("=== Polynomial fit to FULL (J2 + 3b) rate ===")
    full_arr = []
    for lj2 in data["lambda_j2_values"]:
        for l3b in data["lambda_3body_values"]:
            full_arr.append(full_rates[(lj2, l3b)])
    full_arr = np.array(full_arr)
    coeffs_full, residuals_full, rank_full, sv_full = np.linalg.lstsq(A, full_arr, rcond=None)
    rms_full = np.sqrt(np.mean((full_arr - A @ coeffs_full) ** 2))
    if n_data > n_coeffs:
        sigma2_full = np.sum((full_arr - A @ coeffs_full) ** 2) / (n_data - n_coeffs)
        cov_full = sigma2_full * np.linalg.inv(A.T @ A)
        sigma_full = np.sqrt(np.diag(cov_full))
    else:
        sigma_full = np.full_like(coeffs_full, float("nan"))
    print(f"  {'name':18s} {'value':>15s} {'sigma':>15s} {'SNR':>8s}")
    for i, name in enumerate(names):
        snr = abs(coeffs_full[i]) / sigma_full[i] if sigma_full[i] > 0 else float("nan")
        print(f"  {name:18s} {coeffs_full[i]:+15.6e} {sigma_full[i]:15.6e} {snr:8.2f}")
    print(f"  RMS residual: {rms_full:.4e} deg/day")

    print()
    print("=== Verdict on J2-only scaling ===")
    a10 = coeffs_full[1]
    a10_sigma = sigma_full[1]
    snr_a10 = abs(a10) / a10_sigma if a10_sigma > 0 else float("nan")
    print(f"  J2 coefficient a10 = {a10:+.6e} +/- {a10_sigma:.6e}, SNR = {snr_a10:.2f}")
    if snr_a10 > 100:
        print(f"  -> J2 response is LINEAR (a10 dominates, J2^2 a20 = {coeffs_full[4]:+.2e})")
    else:
        print(f"  -> J2 response shows some non-linearity")


if __name__ == "__main__":
    main()
