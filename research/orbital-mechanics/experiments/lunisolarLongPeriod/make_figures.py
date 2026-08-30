"""Experiment 019 -- Lunisolar Long-Period Terms figures.

Deterministic Agg backend, byte-stable for given input.

Figures:
- F1: window-length sensitivity (slopes vs W; corrected cf overlay)
- F2: window-length extrapolation fit (1/W, 1/W^2)
- F3: cycle-averaged estimator (12 monthly segments; mean and std)
- F4: FFT periodicity (dominant periods; annual, evection, variation)
- F5: Lunisolar decomposition (solar vs lunar; corrected vs numerical)
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REFERENCE_DIR = Path(__file__).resolve().parent
RESULTS_PATH = REFERENCE_DIR / "results" / "results.json"
FIG_DIR = REFERENCE_DIR / "results" / "figures"


def _load():
    return json.loads(RESULTS_PATH.read_text(encoding="utf-8"))["results"]


def _save(fig, name):
    fig.tight_layout()
    p = FIG_DIR / name
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p.name


def figure_window_sensitivity(payload):
    """F1: window-length sensitivity at h=600 km, i_sso and i=90."""
    by_alt = payload["window_sweeps"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)

    for ax, i_label in zip(axes, ("i97.79", "i90.00")):
        modes_shown = []
        for mode in ("sun_moon_j2", "sun_moon", "moon_only", "sun_only"):
            key = f"{i_label}_{mode}"
            if key not in by_alt:
                continue
            modes_shown.append(mode)
            data = by_alt[key]
            W = []
            S = []
            for w_str, r in data.items():
                W.append(float(w_str))
                S.append(r["slope_deg_per_day"])
            order = np.argsort(W)
            W = np.array(W)[order]
            S = np.array(S)[order]
            ax.plot(W, S, "o-", label=mode)

        # Overlay corrected cf
        cf_by_incl = payload["corrected_closed_form_by_inclination"]
        i_deg = 97.79 if i_label == "i97.79" else 90.0
        cf_total = cf_by_incl[f"i_{i_deg:.2f}"]["total_cf_deg_day"]
        cf_solar = cf_by_incl[f"i_{i_deg:.2f}"]["solar_cf_deg_day"]
        cf_lunar = cf_by_incl[f"i_{i_deg:.2f}"]["lunar_cf_deg_day"]
        # Reference baseline at zero
        ax.axhline(0, color="gray", lw=0.5, alpha=0.5)
        ax.set_xlabel("window length W (days)")
        ax.set_ylabel("slope (deg/day)")
        ax.set_xscale("log")
        ax.set_title(f"Window sensitivity at i={i_deg}°")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        # Text annotation of corrected cf
        ax.text(0.05, 0.05,
                f"corrected cf:\n  total={cf_total:+.4e}\n  solar={cf_solar:+.4e}\n  lunar={cf_lunar:+.4e} deg/day",
                transform=ax.transAxes, fontsize=8,
                verticalalignment="bottom",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    fig.suptitle("019 F1: Window-length sensitivity (full model + isolation modes)", fontsize=10)
    return _save(fig, "f1_window_length_sensitivity.png")


def figure_extrapolation(payload):
    """F2: window-length extrapolation fit Omega_dot(W) = a + b/W + c/W^2."""
    by_alt = payload["window_sweeps"]
    extraps = payload["window_length_extrapolation"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for ax, i_label in zip(axes, ("i97.79", "i90.00")):
        for mode in ("sun_moon_j2",):
            key = f"{i_label}_{mode}"
            if key not in by_alt:
                continue
            data = by_alt[key]
            W = np.array(sorted(float(w) for w in data.keys()))
            S = np.array([data[f"{w:.0f}"]["slope_deg_per_day"] for w in W])
            ax.plot(W, S, "o", label=f"numerical ({mode})", markersize=8)

            extrap = extraps[key]
            # Quadratic fit curve
            W_smooth = np.logspace(np.log10(W.min() / 2), np.log10(5000), 100)
            pred = (extrap["extrapolated_secular_deg_day"]
                    + extrap["b_1_over_W"] / W_smooth
                    + extrap["c_1_over_W_squared"] / W_smooth ** 2)
            ax.plot(W_smooth, pred, "--", color="gray",
                    label=f"1/W+1/W^2 fit\na={extrap['extrapolated_secular_deg_day']:+.4e}")

            # Linear 1/W fit
            pred_lin = (extrap["linear_1_over_W_secular_deg_day"]
                        + extrap["linear_1_over_W_b"] / W_smooth)
            ax.plot(W_smooth, pred_lin, ":", color="black",
                    label=f"linear 1/W fit\na={extrap['linear_1_over_W_secular_deg_day']:+.4e}")

            # Corrected cf overlay
            cf_by_incl = payload["corrected_closed_form_by_inclination"]
            i_deg = 97.79 if i_label == "i97.79" else 90.0
            cf_total = cf_by_incl[f"i_{i_deg:.2f}"]["total_cf_deg_day"]
            ax.axhline(cf_total, color="red", lw=1,
                       label=f"corrected cf\n+{cf_total:+.4e} deg/day")

            ax.set_xscale("log")
            ax.set_xlabel("window length W (days)")
            ax.set_ylabel("slope (deg/day)")
            ax.set_title(f"i={i_deg}°: window-length extrapolation")
            ax.legend(fontsize=7, loc="best")
            ax.grid(True, alpha=0.3)

    fig.suptitle("019 F2: Window-length extrapolation (secular-limit at W → ∞)", fontsize=10)
    return _save(fig, "f2_window_length_extrapolation.png")


def figure_cycle_averaged(payload):
    """F3: cycle-averaged estimator (12 monthly segments)."""
    cyc = payload["cycle_averaged_estimator"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for ax, key in zip(axes, ("i97.79_sun_moon_j2", "i90.00_sun_moon_j2")):
        if key not in cyc:
            continue
        data = cyc[key]
        slopes = np.array(data["slopes_deg_day"])
        n = len(slopes)
        x = np.arange(n)
        ax.bar(x, slopes, color="steelblue", edgecolor="black")
        ax.axhline(data["mean_deg_day"], color="red", lw=1,
                   label=f"mean = {data['mean_deg_day']:+.4e} deg/day")
        ax.axhline(data["mean_deg_day"] + data["std_deg_day"], color="gray", lw=0.5, ls="--")
        ax.axhline(data["mean_deg_day"] - data["std_deg_day"], color="gray", lw=0.5, ls="--")
        ax.set_xlabel("monthly segment index (k)")
        ax.set_ylabel("slope (deg/day)")
        ax.set_title(f"Cycle-averaged estimator ({n} segments)\n{key}")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle("019 F3: Cycle-averaged estimator (12 monthly segments)", fontsize=10)
    return _save(fig, "f3_cycle_averaged.png")


def figure_fft(payload):
    """F4: FFT periodicity (top dominant periods)."""
    fft_sso = payload["fft_periodicity_i_sso"]
    fft_90 = payload["fft_periodicity_i_90"]

    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    for ax, fft, label in zip(axes, (fft_sso, fft_90),
                              ("i_sso (97.79°)", "i=90°")):
        periods = fft["dominant_periods_day"]
        amplitudes = fft["dominant_amplitudes_deg"]
        x = np.arange(len(periods))
        ax.bar(x, amplitudes, color="steelblue", edgecolor="black")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{p:.1f}" for p in periods], rotation=45)
        ax.set_xlabel("period (days)")
        ax.set_ylabel("amplitude (deg)")
        ax.set_title(f"FFT top-5 periods at h=600 km, {label} (n={fft['n_points']} samples)")
        ax.grid(True, alpha=0.3)
        # Reference markers
        for ref_period, ref_label in [
            (365.24, "annual"),
            (27.55, "evection"),
            (14.77, "variation"),
        ]:
            ax.axvline(ref_period, color="red", lw=0.5, ls="--", alpha=0.5)
            ax.text(ref_period, max(amplitudes) * 1.1, ref_label,
                    rotation=90, fontsize=7, color="red", alpha=0.7)
    fig.suptitle("019 F4: FFT periodicity of osculating Omega(t) at h=600 km", fontsize=10)
    return _save(fig, "f4_fft_periodicity.png")


def figure_lunisolar_decomposition(payload):
    """F5: Lunisolar decomposition (solar vs lunar at i_sso only)."""
    by_alt = payload["window_sweeps"]
    fig, axes = plt.subplots(1, 1, figsize=(8, 5))
    ax = axes

    i_label = "i97.79"
    sun_key = f"{i_label}_sun_only"
    moon_key = f"{i_label}_moon_only"
    full_key = f"{i_label}_sun_moon_j2"
    W = sorted(float(w) for w in by_alt[full_key].keys())
    sun_data = []
    moon_data = []
    full_data = []
    for w in W:
        sun = by_alt[sun_key][f"{w:.0f}"]["slope_deg_per_day"]
        moon = by_alt[moon_key][f"{w:.0f}"]["slope_deg_per_day"]
        full = by_alt[full_key][f"{w:.0f}"]["slope_deg_per_day"]
        # Lunisolar isolation
        sun_data.append(full - moon)
        moon_data.append(full - sun)
        full_data.append(full)

    ax.plot(W, sun_data, "o-", color="orange", label="solar (full - moon_only)")
    ax.plot(W, moon_data, "s-", color="blue", label="lunar (full - sun_only)")
    ax.plot(W, full_data, "^-", color="black", label="full model")

    cf_by_incl = payload["corrected_closed_form_by_inclination"]
    cf = cf_by_incl[f"i_97.79"]
    ax.axhline(cf["solar_cf_deg_day"], color="orange", lw=0.7, ls="--",
               label=f"cf solar={cf['solar_cf_deg_day']:+.4e}")
    ax.axhline(cf["lunar_cf_deg_day"], color="blue", lw=0.7, ls="--",
               label=f"cf lunar={cf['lunar_cf_deg_day']:+.4e}")
    ax.axhline(cf["total_cf_deg_day"], color="red", lw=0.7, ls="--",
               label=f"cf total={cf['total_cf_deg_day']:+.4e}")

    ax.set_xlabel("window length W (days)")
    ax.set_ylabel("slope (deg/day)")
    ax.set_xscale("log")
    ax.set_title(f"Lunisolar decomposition at h=600 km, i=97.79° (i_sso)")
    ax.legend(fontsize=7, loc="best")
    ax.grid(True, alpha=0.3)
    fig.suptitle("019 F5: Lunisolar decomposition (solar vs lunar; corrected cf overlay)", fontsize=10)
    return _save(fig, "f5_lunisolar_decomposition.png")


def main():
    payload = _load()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    paths.append(figure_window_sensitivity(payload))
    paths.append(figure_extrapolation(payload))
    paths.append(figure_cycle_averaged(payload))
    paths.append(figure_fft(payload))
    paths.append(figure_lunisolar_decomposition(payload))
    print(f"[019] {len(paths)} figures written:")
    for p in paths:
        print(f"  - {p}")


if __name__ == "__main__":
    main()