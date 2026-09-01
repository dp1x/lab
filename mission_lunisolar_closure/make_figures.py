"""Generate figures for mission_lunisolar_closure from results.json."""
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


def main():
    if not RESULTS.exists():
        print(f"missing: {RESULTS}; nothing to plot")
        return
    payload = json.loads(RESULTS.read_text(encoding="utf-8"))
    cmp = payload.get("comparison_with_corrected_formula", {})
    hl = payload.get("headline_propagations", {})

    # Figure 1: corrected formula vs numerical across inclinations
    fig, ax = plt.subplots(figsize=(8, 5))
    inclinations = ["i_sso", "i_90", "i_30"]
    incl_labels = ["i=97.79 (SSO)", "i=90", "i=30"]
    cf_total = [cmp.get(f"{k}_cf_total_deg_day", 0.0) for k in inclinations]
    harmonic = [hl.get(f"{k}_lunisolar_contribution", {}).get(
                    "harmonic_regression_deg_per_day", 0.0) for k in inclinations]
    direct = [hl.get(f"{k}_lunisolar_contribution", {}).get(
                  "direct_ols_deg_per_day", 0.0) for k in inclinations]
    node_vec = [hl.get(f"{k}_lunisolar_contribution", {}).get(
                    "node_vector_deg_per_day", 0.0) for k in inclinations]
    x = np.arange(len(inclinations))
    width = 0.18
    ax.bar(x - 2*width, cf_total, width, label="corrected cf")
    ax.bar(x - width, harmonic, width, label="harmonic reg (18.6 yr)")
    ax.bar(x, direct, width, label="direct OLS (18.6 yr)")
    ax.bar(x + width, node_vec, width, label="node-vector (18.6 yr)")
    ax.set_xticks(x)
    ax.set_xticklabels(incl_labels)
    ax.set_ylabel("Lunisolar secular RAAN rate (deg/day)")
    ax.set_title("Mission 1: 18.6-yr DE441 vs corrected doubly-averaged formula\n"
                  "h=600 km; J2+Sun+Moon, J2-only subtracted")
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    plt.tight_layout()
    out1 = OUT_DIR / "fig1_corrected_vs_numerical_by_inclination.png"
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"  saved {out1}")

    # Figure 2: phase-locked estimator at i_sso (if present)
    pl_data = hl.get("i_sso_lunisolar_contribution", {})
    if "phase_locked_avg_deg_per_day" in pl_data:
        fig, ax = plt.subplots(figsize=(8, 5))
        methods = ["direct OLS", "secant", "harmonic reg", "node-vector", "phase-locked avg"]
        vals = [pl_data.get("direct_ols_deg_per_day", 0.0),
                pl_data.get("estimator_g_lunisolar_deg_per_day", 0.0),
                pl_data.get("harmonic_regression_deg_per_day", 0.0),
                pl_data.get("node_vector_deg_per_day", 0.0),
                pl_data.get("phase_locked_avg_deg_per_day", 0.0)]
        cf = cmp.get("i_sso_cf_total_deg_day", 0.0)
        ax.bar(methods, vals, label="Lunisolar contribution")
        ax.axhline(cf, color='r', linestyle='--', label=f"corrected cf = {cf:+.3e}")
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_ylabel("Lunisolar RAAN rate at i_sso (deg/day)")
        ax.set_title("Mission 1: estimator hierarchy at h=600 km i_sso\n"
                      "18.6-yr direct arc, J2-only subtracted")
        ax.legend()
        plt.xticks(rotation=20, ha='right')
        plt.tight_layout()
        out2 = OUT_DIR / "fig2_estimator_hierarchy_i_sso.png"
        fig.savefig(out2, dpi=150)
        plt.close(fig)
        print(f"  saved {out2}")

    # Figure 3: ratio of numerical to corrected formula
    fig, ax = plt.subplots(figsize=(8, 5))
    if "i_sso_ratio_harmonic_to_cf" in cmp:
        ratios = [cmp.get("i_sso_ratio_harmonic_to_cf", float("nan")),
                  cmp.get("i_90_cf_total_deg_day", 1.0) /
                  (hl.get("i_90_lunisolar_contribution", {}).get(
                       "harmonic_regression_deg_per_day", float("nan")) + 1e-30),
                  cmp.get("i_30_cf_total_deg_day", 1.0) /
                  (hl.get("i_30_lunisolar_contribution", {}).get(
                       "harmonic_regression_deg_per_day", float("nan")) + 1e-30)]
        # Filter NaN
        ratios_clean = [r if math.isfinite(r) else float("nan") for r in ratios]
        x = np.arange(len(inclinations))
        ax.bar(x, ratios_clean, color='steelblue')
        ax.axhline(1.0, color='r', linestyle='--', label="ratio = 1 (perfect agreement)")
        ax.axhspan(0.5, 2.0, alpha=0.15, color='green', label='+/- 2x (VERIFIED-WITH-LIMITATION)')
        ax.set_xticks(x)
        ax.set_xticklabels(incl_labels)
        ax.set_ylabel("numerical / corrected (harmonic reg)")
        ax.set_title("Mission 1: numerical/analytical ratio across inclinations\n"
                      "0.5x to 2x = VERIFIED-WITH-LIMITATION gate")
        ax.legend()
        ax.grid(True, axis='y', alpha=0.3)
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
        ax.set_yscale('log')
        ax.set_ylabel("|bias| (deg/day)")
        ax.set_title(f"Synthetic oracle: estimator bias on known secular = {synth['a_true_deg_day']:.0e}\n"
                      "18.6-yr oracle; harmonic reg (f) recovers to machine precision")
        plt.tight_layout()
        out4 = OUT_DIR / "fig4_synthetic_oracle.png"
        fig.savefig(out4, dpi=150)
        plt.close(fig)
        print(f"  saved {out4}")


if __name__ == "__main__":
    main()