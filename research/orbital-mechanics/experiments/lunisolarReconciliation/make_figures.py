"""Generate figures for Exp 018 from results.json."""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXP = Path(__file__).resolve().parent
RESULTS_PATH = EXP / "results" / "results.json"
FIGDIR = EXP / "results" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)


def make_figures():
    results = json.loads(RESULTS_PATH.read_text(encoding="utf-8"))
    r = results["results"]

    # F1: corrected cf vs numerical by altitude (Lunisolar-only contribution)
    fig, ax = plt.subplots(figsize=(9, 5))
    alts = [500, 600, 700, 800]
    corrected = [r["corrected_closed_form_by_altitude"][str(h)]["total_cf_deg_day"] for h in alts]
    ax.semilogy(alts, np.abs(corrected), "o-", color="darkorange",
                label="Corrected secular (Track B)")
    # 017 numerical (pre-registration): h=500/600/700/800 have +0.001320, +0.001284, +0.001249, +0.001215
    numerical_017 = [0.001320, 0.001284, 0.001249, 0.001215]
    ax.semilogy(alts, numerical_017, "s-", color="steelblue",
                label="017 numerical (1-yr fit, J2-sub)")
    # Compute 018 numerical from force isolation
    j2_018 = r["force_isolation_h600"]["j2_only"]["slope_deg_per_day"]
    # Use the corrected cf slope as a baseline; the 018 sun_moon_j2 slope at h=600 is
    # the numerical. For other altitudes, we don't have direct measurements,
    # so use 017 numerical as the reference (consistent with the 018 finding)
    ax.set_title(
        "Exp 018 F1: Corrected secular Lunisolar RAAN vs 017/018 numerical\n"
        "(byte-pinned DE441 Sun+Moon, 1-year arc)"
    )
    ax.set_xlabel("altitude (km)")
    ax.set_ylabel("|Lunisolar RAAN drift| (deg/day, log scale)")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGDIR / "f1_corrected_cf_vs_numerical.png", dpi=150)
    plt.close(fig)

    # F2: inclination sweep (slope vs inclination)
    fig, ax = plt.subplots(figsize=(9, 5))
    incl_data = r["inclination_sweep_h600"]
    incl_keys = sorted(incl_data.keys(), key=lambda x: float(x))
    incl_x = [float(k) for k in incl_keys]
    incl_y = [incl_data[k]["slope_deg_per_day"] for k in incl_keys]
    ax.plot(incl_x, incl_y, "o-", color="darkblue", label="018 numerical (1-yr fit)")
    # Overlay J2 only prediction
    R_E = 6378.137
    J2 = 0.001082629821
    a = R_E + 600.0
    import math as _math
    mu_E = 398600.4418
    n = _math.sqrt(mu_E / a ** 3)
    j2_pred_x = np.linspace(0, 180, 100)
    j2_pred_y = -1.5 * n * J2 * (R_E / a) ** 2 * np.cos(np.radians(j2_pred_x)) * 86400.0 * 180.0 / _math.pi
    ax.plot(j2_pred_x, j2_pred_y, "--", color="gray", alpha=0.5, label="J2-only prediction")
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(97.7876, color="red", lw=0.5, ls=":", alpha=0.5, label="i_sso(600)")
    ax.axvline(82.2124, color="orange", lw=0.5, ls=":", alpha=0.5, label="180-i_sso")
    ax.set_title(
        "Exp 018 F2: Sun+Moon+J2 RAAN drift vs inclination (h=600 km)\n"
        "J2 dominates; at i=90 deg J2 cos(i)=0 and Lunisolar secular is visible"
    )
    ax.set_xlabel("inclination (deg)")
    ax.set_ylabel("total RAAN drift (deg/day)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGDIR / "f2_inclination_sweep.png", dpi=150)
    plt.close(fig)

    # F3: window-length sensitivity
    fig, ax = plt.subplots(figsize=(9, 5))
    win_data = r["window_sensitivity_h600"]
    win_keys = sorted(win_data.keys(), key=lambda x: float(x))
    win_x = [float(k) for k in win_keys]
    win_y = [win_data[k]["slope_deg_per_day"] for k in win_keys]
    ax.plot(win_x, win_y, "o-", color="seagreen", label="018 numerical (1-yr or shorter fit)")
    ax.axhline(0.992014, color="gray", lw=0.5, ls="--", label="J2-only slope (control)")
    # Overlay Lunisolar residual
    j2 = r["force_isolation_h600"]["j2_only"]["slope_deg_per_day"]
    win_residual = [s - j2 for s in win_y]
    ax2 = ax.twinx()
    ax2.plot(win_x, win_residual, "x-", color="red", alpha=0.7,
             label="Lunisolar residual (slope - J2)")
    ax.set_title(
        "Exp 018 F3: Window-length sensitivity at h=600 km i_sso\n"
        "(Lunisolar-only contribution; secular term is the long-window limit)"
    )
    ax.set_xlabel("fit window (days)")
    ax.set_ylabel("total slope (deg/day)")
    ax2.set_ylabel("Lunisolar residual (deg/day)", color="red")
    ax2.tick_params(axis="y", labelcolor="red")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGDIR / "f3_window_sensitivity.png", dpi=150)
    plt.close(fig)

    # F4: precession on/off
    fig, ax = plt.subplots(figsize=(8, 5))
    with_p = r["precession_comparison_h600"]["with_precession"]["sun_moon_j2"]["slope_deg_per_day"]
    without_p = r["precession_comparison_h600"]["without_precession"]["sun_moon_j2"]["slope_deg_per_day"]
    ax.bar(["without precession", "with precession"],
           [without_p, with_p], color=["lightgray", "seagreen"], edgecolor="black")
    ax.set_title(
        f"Exp 018 F4: Precession on/off at h=600 km i_sso (1-yr arc)\n"
        f"Frame-mismatch bias: {(with_p-without_p)*365.2422:.3f} deg/year"
    )
    ax.set_ylabel("Sun+Moon+J2 RAAN drift (deg/day)")
    fig.tight_layout()
    fig.savefig(FIGDIR / "f4_precession_comparison.png", dpi=150)
    plt.close(fig)

    # F5: convergence ladder
    fig, ax = plt.subplots(figsize=(9, 5))
    conv = r["convergence"]
    dt_s = conv["dt_s"]
    r_diff = conv["max_r_diff_km"]
    ax.loglog(dt_s, r_diff, "o-", color="steelblue", label="numerical vs dt_ref")
    if r_diff[0] > 0 and r_diff[-1] > 0:
        c = r_diff[-1] / (dt_s[-1] ** 4)
        ref_line = [c * d ** 4 for d in dt_s]
        ax.loglog(dt_s, ref_line, "--", color="gray", label="order-4 reference")
    ax.set_title(
        f"Exp 018 F5: dt convergence ladder at h=600 km (1-day arc)\n"
        f"fitted order p_r = {conv['p_r']:.2f}, p_v = {conv['p_v']:.2f}"
    )
    ax.set_xlabel("dt (s)")
    ax.set_ylabel("max |r(dt) - r(dt_ref)| (km)")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGDIR / "f5_convergence_ladder.png", dpi=150)
    plt.close(fig)

    # F6: force isolation decomposition
    fig, ax = plt.subplots(figsize=(9, 5))
    fi = r["force_isolation_h600"]
    j2_slope = fi["j2_only"]["slope_deg_per_day"]
    sun_contrib = fi["sun_only"]["slope_deg_per_day"] - j2_slope
    moon_contrib = fi["moon_only"]["slope_deg_per_day"] - j2_slope
    print(f"sun_contrib: {sun_contrib:+.4e}")
    print(f"moon_contrib: {moon_contrib:+.4e}")
    cf_600 = r["corrected_closed_form_by_altitude"]["600"]
    labels = ["Solar\n(numerical)", "Lunar\n(numerical)",
              "Solar\n(corrected cf)", "Lunar\n(corrected cf)"]
    values = [sun_contrib, moon_contrib,
              cf_600["solar_cf_deg_day"], cf_600["lunar_cf_deg_day"]]
    colors = ["steelblue", "darkblue", "darkorange", "firebrick"]
    bars = ax.bar(labels, values, color=colors, edgecolor="black")
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, v + (0.0001 if v > 0 else -0.0001),
                f"{v:+.4e}", ha="center", fontsize=8,
                va="bottom" if v > 0 else "top")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_title(
        "Exp 018 F6: Lunisolar RAAN decomposition at h=600 km i_sso\n"
        "Numerical (1-yr fit, J2-sub) vs Corrected Secular (Track B)"
    )
    ax.set_ylabel("Lunisolar RAAN contribution (deg/day)")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(FIGDIR / "f6_lunisolar_decomposition.png", dpi=150)
    plt.close(fig)

    return ["f1_corrected_cf_vs_numerical.png", "f2_inclination_sweep.png",
            "f3_window_sensitivity.png", "f4_precession_comparison.png",
            "f5_convergence_ladder.png", "f6_lunisolar_decomposition.png"]


if __name__ == "__main__":
    paths = make_figures()
    for p in paths:
        print(f"[018] figure: results/figures/{p}")
