"""Synthetic estimator test for Experiment 020 / Track 3.

Self-contained READ-ONLY analysis script. It does NOT touch any
production code. It builds a synthetic signal

    y(t) = a_secular * t + sum_k A_k * cos(omega_k * t + phi_k) + small noise

with the actual physical frequencies and FFT amplitudes from the
019 results.json (fft_periodicity_i_sso top-5 dominant periods and
amplitudes), samples it at the 019 ascending-node cadence (~14.91/day),
and applies the eight estimators specified in the task:

    (a) Direct OLS slope over [0, W]
    (b) OLS slope over [0, W] with the secular removed first
        (fit residual to zero-mean intercept = same as (a) numerically)
    (c) Polynomial in 1/W fit to the W-dependent slope
    (d) Cycle-averaged estimator (K equal segments, each ~T/N)
    (e) Linear-fit slope with the dominant periodic term subtracted
    (f) Linear-fit slope with N terms subtracted (theory-driven harmonic regression)
    (g) Direct secant estimator: (y(T) - y(0)) / T (no fitting)
    (h) Medians-of-segments estimator (median of per-segment slopes)

For each estimator, we report:
    - recovered slope at W in {365, 730, 1825, 3650} d
    - bias (recovered - true secular)
    - convergence rate as a function of W

The 'true secular' is set to a_secular = 1.0e-3 deg/day (the order of the
018 numerical Lunisolar contribution at i_sso).

This script is deliberately simple and deterministic (numpy only, no
network). It is a test of the estimator's theory, not a benchmark of
the 019 numerical pipeline.
"""
from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Synthetic signal definition (FACT: from 019 results.json)
# --------------------------------------------------------------------------- #
A_SECULAR_DEG_DAY = 1.0e-3  # deg/day

N_PER_DAY = 14.91  # crossings per day
DT_DAYS = 1.0 / N_PER_DAY

# Physical harmonics (FACT: from 019 results.json fft_periodicity_i_sso
# top-5 dominant periods and amplitudes). The task spec also requires the
# named physical drivers (evection 27.55 d, variation 14.77 d, lunar nodal
# 18.6 yr). We use 019 FFT amplitudes for the visible harmonics and
# standard Kaula/Murray-Dermott orders for evection/variation
# (Track B sections 4.1, 4.2).
PHYSICAL_HARMONICS = [
    # (period_days, amplitude_deg, label)
    (365.2422, 0.103, "annual solar forcing"),         # 019 top-1
    (182.6, 0.025, "half-annual solar"),                # 019 top-2
    (121.675, 0.012, "evection alias at 121.7 d"),      # 019 top-3 (evection residual)
    (91.256, 0.007, "lunar annual modulation"),         # 019 top-4
    (73.005, 0.005, "tertiary beat"),                   # 019 top-5
    (27.5546, 0.004, "evection direct (27.55 d)"),      # Kaula amplitude
    (14.7653, 0.003, "variation direct (14.77 d)"),     # Kaula amplitude
    (6798.4, 0.002, "lunar nodal (18.6 yr)"),           # Kaula upper bound
]

PHI_K = 0.0  # rad; for all harmonics (worst-case aligned)
NOISE_STD_DEG = 1.0e-4
RNG_SEED = 42

# Window lengths to test (synthetic extrapolation)
W_DAYS_LIST = (30.0, 90.0, 180.0, 365.0, 730.0, 1825.0, 3650.0)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def build_signal(t_days: np.ndarray) -> tuple:
    """Return (y_true, y_noisy) at the query times t_days."""
    y = A_SECULAR_DEG_DAY * t_days
    for period_days, amp_deg, _label in PHYSICAL_HARMONICS:
        omega = 2.0 * math.pi / period_days
        y = y + amp_deg * np.cos(omega * t_days + PHI_K)
    rng = np.random.default_rng(RNG_SEED)
    noise = rng.normal(0.0, NOISE_STD_DEG, size=t_days.shape)
    return y, y + noise


def ols_slope(t: np.ndarray, y: np.ndarray) -> float:
    """OLS slope of y on t."""
    t_mean = np.mean(t)
    y_mean = np.mean(y)
    cov = np.mean((t - t_mean) * (y - y_mean))
    var = np.mean((t - t_mean) ** 2)
    if var == 0:
        return float("nan")
    return float(cov / var)


def ols_intercept(t: np.ndarray, y: np.ndarray) -> tuple:
    """Full OLS: returns (intercept, slope)."""
    A = np.column_stack([np.ones_like(t), t])
    res = np.linalg.lstsq(A, y, rcond=None)
    return float(res[0][0]), float(res[0][1])


# --------------------------------------------------------------------------- #
# Estimators
# --------------------------------------------------------------------------- #
def estimator_a_direct_ols(t: np.ndarray, y: np.ndarray) -> float:
    """(a) Direct OLS slope over [0, W]."""
    return ols_slope(t, y)


def estimator_b_ols_residual(t: np.ndarray, y: np.ndarray) -> float:
    """(b) OLS with secular removed first; numerically identical to (a)
    because the secular IS the linear fit's center of mass in t."""
    _, slope = ols_intercept(t, y)
    return slope


def estimator_c(W_arr_days: np.ndarray, slopes: np.ndarray) -> dict:
    """(c) Omega_dot_fit(W) = a + b/W + c/W^2 fit (matches 019 model)."""
    x1 = 1.0 / W_arr_days
    x2 = 1.0 / W_arr_days ** 2
    A = np.column_stack([np.ones_like(W_arr_days), x1, x2])
    res = np.linalg.lstsq(A, slopes, rcond=None)
    a_mean = float(res[0][0])
    b_lin = float(res[0][1])
    c_quad = float(res[0][2])
    pred = a_mean + b_lin * x1 + c_quad * x2
    rms = float(np.sqrt(np.mean((slopes - pred) ** 2)))
    return {
        "intercept_a": a_mean,
        "b_1_over_W": b_lin,
        "c_1_over_W_sq": c_quad,
        "rms_residual": rms,
    }


def estimator_d_cycle_averaged(t: np.ndarray, y: np.ndarray,
                               n_segments: int = 12) -> dict:
    """(d) Cycle-averaged estimator (K equal segments)."""
    t_lo, t_hi = t[0], t[-1]
    seg_len = (t_hi - t_lo) / n_segments
    slopes = []
    for k in range(n_segments):
        mask = (t >= t_lo + k * seg_len) & (t <= t_lo + (k + 1) * seg_len)
        if mask.sum() < 5:
            continue
        slopes.append(ols_slope(t[mask], y[mask]))
    if not slopes:
        return {"mean": float("nan"), "std": float("nan"), "slopes": []}
    return {
        "mean": float(np.mean(slopes)),
        "std": float(np.std(slopes)),
        "slopes": [float(s) for s in slopes],
    }


def estimator_e_dominant_removed(t: np.ndarray, y: np.ndarray,
                                 period_days: float,
                                 amp_deg: float) -> float:
    """(e) Linear-fit slope with the dominant periodic term subtracted."""
    omega = 2.0 * math.pi / period_days
    y_corr = y - amp_deg * np.cos(omega * t + PHI_K)
    return ols_slope(t, y_corr)


def estimator_f_harmonic_regression(t: np.ndarray, y: np.ndarray,
                                    harmonics: list) -> float:
    """(f) Theory-driven harmonic regression: fit y = a + b*t +
    sum_k [c_k cos(omega_k t) + s_k sin(omega_k t)] and report b."""
    cols = [np.ones_like(t), t]
    for period_days, _amp, _label in harmonics:
        omega = 2.0 * math.pi / period_days
        cols.append(np.cos(omega * t))
        cols.append(np.sin(omega * t))
    A = np.column_stack(cols)
    res = np.linalg.lstsq(A, y, rcond=None)
    return float(res[0][1])


def estimator_g_secant(t: np.ndarray, y: np.ndarray) -> float:
    """(g) Direct secant: (y(T) - y(0)) / T (no fitting)."""
    T = t[-1] - t[0]
    if T <= 0:
        return float("nan")
    return float((y[-1] - y[0]) / T)


def estimator_h_median_of_segments(t: np.ndarray, y: np.ndarray,
                                   n_segments: int = 12) -> dict:
    """(h) Medians-of-segments: split [0, T] into K segments, take median slope."""
    t_lo, t_hi = t[0], t[-1]
    seg_len = (t_hi - t_lo) / n_segments
    slopes = []
    for k in range(n_segments):
        mask = (t >= t_lo + k * seg_len) & (t <= t_lo + (k + 1) * seg_len)
        if mask.sum() < 5:
            continue
        slopes.append(ols_slope(t[mask], y[mask]))
    if not slopes:
        return {"median": float("nan"), "mad": float("nan")}
    return {
        "median": float(statistics.median(slopes)),
        "mad": float(statistics.median(
            [abs(s - statistics.median(slopes)) for s in slopes]
        )),
    }


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run_window(W_days: float) -> dict:
    n_samples = max(20, int(round(W_days * N_PER_DAY)))
    t = np.linspace(0.0, W_days, n_samples)
    _y_true, y_noisy = build_signal(t)

    a_direct = estimator_a_direct_ols(t, y_noisy)
    b_resid = estimator_b_ols_residual(t, y_noisy)
    d_cycle = estimator_d_cycle_averaged(t, y_noisy, n_segments=12)
    e_dom = estimator_e_dominant_removed(
        t, y_noisy, period_days=365.2422, amp_deg=0.103)
    f_harm = estimator_f_harmonic_regression(t, y_noisy, PHYSICAL_HARMONICS)
    g_sec = estimator_g_secant(t, y_noisy)
    h_med = estimator_h_median_of_segments(t, y_noisy, n_segments=12)

    return {
        "W_days": W_days,
        "n_samples": n_samples,
        "estimators": {
            "a_direct_ols": a_direct,
            "b_ols_residual": b_resid,
            "d_cycle_averaged": d_cycle["mean"],
            "d_cycle_averaged_std": d_cycle["std"],
            "e_dominant_removed": e_dom,
            "f_harmonic_regression": f_harm,
            "g_secant": g_sec,
            "h_median_segments": h_med["median"],
        },
        "bias": {
            "a_direct_ols": a_direct - A_SECULAR_DEG_DAY,
            "b_ols_residual": b_resid - A_SECULAR_DEG_DAY,
            "d_cycle_averaged": d_cycle["mean"] - A_SECULAR_DEG_DAY,
            "e_dominant_removed": e_dom - A_SECULAR_DEG_DAY,
            "f_harmonic_regression": f_harm - A_SECULAR_DEG_DAY,
            "g_secant": g_sec - A_SECULAR_DEG_DAY,
            "h_median_segments": h_med["median"] - A_SECULAR_DEG_DAY,
        },
    }


def main() -> dict:
    print(f"[track-3] true secular = {A_SECULAR_DEG_DAY:+.6e} deg/day")
    print(f"[track-3] n_harmonics = {len(PHYSICAL_HARMONICS)}")
    print(f"[track-3] sample cadence = {N_PER_DAY} / day")

    per_window = []
    for W in W_DAYS_LIST:
        rec = run_window(W)
        per_window.append(rec)
        print(f"\n[track-3] W = {W} d, n = {rec['n_samples']} samples:")
        for name, val in rec["estimators"].items():
            print(f"    {name:30s} = {val:+.6e}")

    W_arr = np.array([r["W_days"] for r in per_window])
    slope_a = np.array([r["estimators"]["a_direct_ols"] for r in per_window])
    slope_d = np.array([r["estimators"]["d_cycle_averaged"] for r in per_window])
    slope_e = np.array([r["estimators"]["e_dominant_removed"] for r in per_window])
    slope_f = np.array([r["estimators"]["f_harmonic_regression"] for r in per_window])
    slope_g = np.array([r["estimators"]["g_secant"] for r in per_window])

    extrap_a = estimator_c(W_arr, slope_a)
    extrap_d = estimator_c(W_arr, slope_d)
    extrap_e = estimator_c(W_arr, slope_e)
    extrap_f = estimator_c(W_arr, slope_f)
    extrap_g = estimator_c(W_arr, slope_g)

    print("\n[track-3] Window-length extrapolation (1/W + 1/W^2 model):")
    print(f"  (a) direct OLS       intercept = {extrap_a['intercept_a']:+.6e} "
          f"(bias = {extrap_a['intercept_a']-A_SECULAR_DEG_DAY:+.3e})")
    print(f"  (d) cycle-averaged   intercept = {extrap_d['intercept_a']:+.6e} "
          f"(bias = {extrap_d['intercept_a']-A_SECULAR_DEG_DAY:+.3e})")
    print(f"  (e) dominant removed intercept = {extrap_e['intercept_a']:+.6e} "
          f"(bias = {extrap_e['intercept_a']-A_SECULAR_DEG_DAY:+.3e})")
    print(f"  (f) harmonic regress intercept = {extrap_f['intercept_a']:+.6e} "
          f"(bias = {extrap_f['intercept_a']-A_SECULAR_DEG_DAY:+.3e})")
    print(f"  (g) secant           intercept = {extrap_g['intercept_a']:+.6e} "
          f"(bias = {extrap_g['intercept_a']-A_SECULAR_DEG_DAY:+.3e})")

    print("\n[track-3] Convergence rate log|bias| vs log(W):")
    for est_name in ["a_direct_ols", "d_cycle_averaged",
                     "e_dominant_removed", "f_harmonic_regression",
                     "g_secant", "h_median_segments"]:
        bias_arr = np.array([abs(r["bias"][est_name]) for r in per_window])
        valid = (bias_arr > 1e-20) & np.isfinite(bias_arr)
        if valid.sum() >= 2:
            log_w = np.log(W_arr[valid])
            log_b = np.log(bias_arr[valid])
            slope_fit = float(np.polyfit(log_w, log_b, 1)[0])
            print(f"    {est_name:30s}: slope = {slope_fit:+.3f}  "
                  f"(i.e., bias ~ W^{slope_fit:.2f})")

    out = {
        "true_secular_deg_day": A_SECULAR_DEG_DAY,
        "n_harmonics": len(PHYSICAL_HARMONICS),
        "harmonics": [
            {"period_days": p, "amplitude_deg": a, "label": l}
            for (p, a, l) in PHYSICAL_HARMONICS
        ],
        "per_window": per_window,
        "extrap_a": extrap_a,
        "extrap_d": extrap_d,
        "extrap_e": extrap_e,
        "extrap_f": extrap_f,
        "extrap_g": extrap_g,
    }
    out_path = (Path(__file__).resolve().parent
                / "audit-020-track-3-synthetic-results.json")
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n[track-3] results -> {out_path}")
    return out


if __name__ == "__main__":
    main()