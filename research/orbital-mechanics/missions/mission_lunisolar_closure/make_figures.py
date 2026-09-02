"""Generate figures for mission_lunisolar_closure from results.json.

Reads the actual schema produced by run_focused_campaign.py:
  payload["lunisolar_estimates"][incl_name]    -> per-inclination estimator values
  payload["comparison_with_corrected_formula"]-> cf values + i_sso ratio
  payload["synthetic_estimator_test"]          -> synthetic-oracle diagnostics
  payload["headline_propagations"]             -> per-mode per-inclination raw analyses
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "results.json"
OUT_DIR = HERE / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


INCLINATIONS = ["i_sso", "i_90", "i_30"]
INCL_LABELS = {"i_sso": "i=97.79 (SSO)", "i_90": "i=90 (J2-clean)", "i_30": "i=30 (prograde)"}


def _get_ls(payload, incl, key, default=float("nan")):
    return float(payload.get("lunisolar_estimates", {}).get(incl, {}).get(key, default))


def _get_cf(payload, incl, default=float("nan")):
    return float(payload.get("comparison_with_corrected_formula", {}).get(
        f"{incl}_cf_total_deg_day", default))


def main():
    if not RESULTS.exists():
        print(f"missing: {RESULTS}; nothing to plot")
        return
    top = json.loads(RESULTS.read_text(encoding="utf-8"))
    # save_json_result wraps the payload under "results"; unwrap if present.
    payload = top.get("results", top)

    # Figure 1: corrected formula vs numerical estimators across inclinations
    fig, ax = plt.subplots(figsize=(9, 5))
    cf_total = [_get_cf(payload, k) for k in INCLINATIONS]
    harmonic = [_get_ls(payload, k, "harmonic_regression_deg_per_day") for k in INCLINATIONS]
    direct = [_get_ls(payload, k, "direct_ols_deg_per_day") for k in INCLINATIONS]
    node_vec = [_get_ls(payload, k, "node_vector_deg_per_day") for k in INCLINATIONS]
    x = np.arange(len(INCLINATIONS))
    width = 0.18
    ax.bar(x - 1.5 * width, cf_total, width, label="corrected cf")
    ax.bar(x - 0.5 * width, harmonic, width, label="harmonic reg (18.6 yr)")
    ax.bar(x + 0.5 * width, direct, width, label="direct OLS (18.6 yr)")
    ax.bar(x + 1.5 * width, node_vec, width, label="node-vector (18.6 yr)")
    ax.set_xticks(x)
    ax.set_xticklabels([INCL_LABELS[k] for k in INCLINATIONS])
    ax.set_ylabel("Lunisolar secular RAAN rate (deg/day)")
    ax.set_title("mission_lunisolar_closure: 18.6-yr DE441 vs corrected doubly-averaged\n"
                  "h=600 km; J2+Sun+Moon minus J2-only control")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    ax.axhline(0, color="k", linewidth=0.5)
    plt.tight_layout()
    out1 = OUT_DIR / "fig1_corrected_vs_numerical_by_inclination.png"
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"  saved {out1}")

    # Figure 2: estimator hierarchy at i_sso
    fig, ax = plt.subplots(figsize=(9, 5))
    methods = ["direct OLS", "secant", "harmonic reg", "node-vector"]
    keys = ["direct_ols_deg_per_day", "secant_deg_per_day",
            "harmonic_regression_deg_per_day", "node_vector_deg_per_day"]
    vals = [_get_ls(payload, "i_sso", k) for k in keys]
    cf = _get_cf(payload, "i_sso")
    bars = ax.bar(methods, vals, label="Lunisolar contribution (full - J2-only)")
    for bar, v in zip(bars, vals):
        if math.isfinite(v):
            ax.text(bar.get_x() + bar.get_width() / 2, v, f"{v:+.3e}",
                     ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
    ax.axhline(cf, color="r", linestyle="--",
                label=f"corrected cf = {cf:+.3e}" if math.isfinite(cf) else "corrected cf = NaN")
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_ylabel("Lunisolar RAAN rate at i_sso (deg/day)")
    ax.set_title("mission_lunisolar_closure: estimator hierarchy at h=600 km i_sso\n"
                  "18.6-yr direct arc, J2-only subtracted")
    ax.legend(loc="best", fontsize=9)
    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()
    out2 = OUT_DIR / "fig2_estimator_hierarchy_i_sso.png"
    fig.savefig(out2, dpi=150)
    plt.close(fig)
    print(f"  saved {out2}")

    # Figure 3: numerical / corrected ratio
    fig, ax = plt.subplots(figsize=(9, 5))
    ratios = []
    for k in INCLINATIONS:
        cf = _get_cf(payload, k)
        num = _get_ls(payload, k, "harmonic_regression_deg_per_day")
        if math.isfinite(cf) and math.isfinite(num) and abs(cf) > 0:
            ratios.append(num / cf)
        else:
            ratios.append(float("nan"))
    ax.bar(INCLINATIONS, ratios, color="steelblue")
    ax.axhline(1.0, color="r", linestyle="--", label="ratio = 1 (perfect)")
    ax.axhspan(0.5, 2.0, alpha=0.15, color="green",
                label="0.5x-2x band (VERIFIED-WITH-LIMITATION gate)")
    ax.set_ylabel("numerical / corrected (harmonic regression)")
    ax.set_title("mission_lunisolar_closure: numerical/analytical ratio across inclinations\n"
                  "harmonic regression at 18.6-yr vs corrected doubly-averaged formula")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    out3 = OUT_DIR / "fig3_numerical_to_cf_ratio.png"
    fig.savefig(out3, dpi=150)
    plt.close(fig)
    print(f"  saved {out3}")

    # Figure 4: synthetic oracle
    synth = payload.get("synthetic_estimator_test", {})
    if "a_true_deg_day" in synth:
        fig, ax = plt.subplots(figsize=(8, 5))
        methods = ["direct OLS (a)", "harmonic reg (f)"]
        bias = [abs(synth.get("estimator_a_bias_deg_day", 0.0)),
                abs(synth.get("estimator_f_bias_deg_day", 0.0))]
        ax.bar(methods, bias)
        ax.set_yscale("log")
        ax.set_ylabel("|bias| (deg/day)")
        ax.set_title(f"Synthetic oracle: estimator bias on known secular = {synth['a_true_deg_day']:.0e}\n"
                      "18.6-yr oracle; harmonic reg (f) recovers to machine precision")
        plt.tight_layout()
        out4 = OUT_DIR / "fig4_synthetic_oracle.png"
        fig.savefig(out4, dpi=150)
        plt.close(fig)
        print(f"  saved {out4}")

    # Figure 5: per-mode raw RAAN drift for i_sso
    hl = payload.get("headline_propagations", {})
    if "i_sso" in hl:
        fig, ax = plt.subplots(figsize=(8, 5))
        for mode in ("j2_only", "sun_moon_j2"):
            ana = hl["i_sso"].get(f"{mode}_analysis", {})
            label = f"{mode}" + (" (J2-only control)" if mode == "j2_only" else " (full)")
            vals = [ana.get(k, float("nan")) for k in
                    ("direct_ols_deg_per_day", "secant_deg_per_day",
                     "harmonic_regression_deg_per_day", "node_vector_deg_per_day")]
            ax.plot(["direct OLS", "secant", "harmonic reg", "node-vector"], vals,
                     marker="o", label=label)
        ax.axhline(0, color="k", linewidth=0.5)
        ax.set_ylabel("RAAN drift rate (deg/day) at i_sso")
        ax.set_title("Per-mode raw estimator values at h=600 km i_sso\n"
                      "Lunisolar contribution = full - J2-only control (Figure 2)")
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3)
        plt.xticks(rotation=15, ha="right")
        plt.tight_layout()
        out5 = OUT_DIR / "fig5_per_mode_raw_i_sso.png"
        fig.savefig(out5, dpi=150)
        plt.close(fig)
        print(f"  saved {out5}")


if __name__ == "__main__":
    main()