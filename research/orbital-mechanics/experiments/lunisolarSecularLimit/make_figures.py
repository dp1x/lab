"""Generate figures for Experiment 020 -- Lunisolar Long-Arc Secular-Limit Validation.

Reads results.json and produces the figures under results/figures/.
Deterministic; uses matplotlib only (no scipy).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_results() -> dict:
    with open(RESULTS_DIR / "results.json", encoding="utf-8") as f:
        wrapper = json.load(f)
    return wrapper["results"]


def f1_estimator_comparison(payload: dict) -> None:
    """Bar chart of the 4 estimators (a, f, g, n) compared to corrected cf."""
    cf = payload["corrected_closed_form_at_i_sso"]["total_cf_deg_day"]
    headline = payload["headline_secular_estimate"]
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = ["Direct OLS (a)", "Harmonic regr. (f)", "Secant (g)", "Node-vector (n)"]
    means = [
        headline["estimator_a_direct_ols"]["mean_deg_day"],
        headline["estimator_f_harmonic_regression"]["mean_deg_day"],
        headline["estimator_g_secant"]["mean_deg_day"],
        headline["node_vector_estimator"]["mean_deg_day"],
    ]
    stds = [
        headline["estimator_a_direct_ols"]["std_deg_day"],
        headline["estimator_f_harmonic_regression"]["std_deg_day"],
        headline["estimator_g_secant"]["std_deg_day"],
        headline["node_vector_estimator"]["std_deg_day"],
    ]
    x = np.arange(len(labels))
    ax.bar(x, means, yerr=stds, capsize=5, color=["C0", "C1", "C2", "C3"])
    ax.axhline(cf, color="red", linestyle="--", label=f"Corrected cf = {cf:+.3e} deg/day")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Lunisolar RAAN rate (deg/day)")
    ax.set_title("Exp 020: Estimator comparison vs corrected 018 cf (h=600 km i_sso)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "f1_estimator_comparison.png", dpi=120)
    plt.close(fig)


def f2_ratio_to_corrected_cf(payload: dict) -> None:
    """Bar chart of the ratio to corrected cf for each estimator."""
    headline = payload["headline_secular_estimate"]
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = ["Direct OLS (a)", "Harmonic regr. (f)", "Secant (g)", "Node-vector (n)"]
    ratios = [
        headline["estimator_a_direct_ols"]["ratio_to_corrected_cf"],
        headline["estimator_f_harmonic_regression"]["ratio_to_corrected_cf"],
        headline["estimator_g_secant"]["ratio_to_corrected_cf"],
        headline["node_vector_estimator"]["ratio_to_corrected_cf"],
    ]
    x = np.arange(len(labels))
    ax.bar(x, ratios, color=["C0", "C1", "C2", "C3"])
    ax.axhline(1.0, color="red", linestyle="--", label="ratio = 1 (perfect match)")
    ax.axhline(9.78, color="orange", linestyle=":", label="018 reported ratio = 9.78")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Ratio (numerical / corrected cf)")
    ax.set_yscale("log")
    ax.set_title("Exp 020: Lunisolar secular-rate ratio (numerical/corrected cf)")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "f2_ratio_to_corrected_cf.png", dpi=120)
    plt.close(fig)


def f3_phase_dependence(payload: dict) -> None:
    """Lunisolar secular estimate vs phase offset."""
    cf = payload["corrected_closed_form_at_i_sso"]["total_cf_deg_day"]
    lunisolar = payload["lunisolar_estimates"]
    if not lunisolar:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    phases_d = []
    estimator_f = []
    estimator_a = []
    estimator_g = []
    estimator_n = []
    for k, v in sorted(lunisolar.items()):
        phases_d.append(v["phase_offset_d"])
        estimator_f.append(v["estimator_f_lunisolar_deg_per_day"])
        estimator_a.append(v["estimator_a_lunisolar_deg_per_day"])
        estimator_g.append(v["estimator_g_lunisolar_deg_per_day"])
        estimator_n.append(v["node_vector_lunisolar_deg_per_day"])
    ax.plot(phases_d, estimator_f, "o-", label="Harmonic regression (f)")
    ax.plot(phases_d, estimator_a, "s-", label="Direct OLS (a)")
    ax.plot(phases_d, estimator_g, "^-", label="Secant (g)")
    ax.plot(phases_d, estimator_n, "v-", label="Node-vector (n)")
    ax.axhline(cf, color="red", linestyle="--", label=f"Corrected cf = {cf:+.3e}")
    ax.set_xlabel("Phase offset (days in lunar anomalistic cycle)")
    ax.set_ylabel("Lunisolar RAAN rate (deg/day)")
    ax.set_title("Exp 020: Phase dependence at h=600 km i_sso (1-yr arc)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "f3_phase_dependence.png", dpi=120)
    plt.close(fig)


def f4_harmonic_amplitudes_recovered(payload: dict) -> None:
    """Bar chart of the harmonic amplitudes recovered by the harmonic
    regression estimator, compared to the theoretical expected amplitudes."""
    synth = payload["synthetic_estimator_test"]
    fig, ax = plt.subplots(figsize=(10, 5))
    expected = [
        ("365.2422", 0.103),
        ("182.6211", 0.025),
        ("121.7474", 0.012),
        ("91.3106", 0.007),
        ("73.0484", 0.005),
        ("27.5546", 0.004),
        ("14.7653", 0.003),
        ("6798.4", 0.002),
    ]
    recovered = synth["harmonic_amplitudes_deg_recovered"]
    x = np.arange(len(expected))
    exp_amp = [e[1] for e in expected]
    rec_amp = [recovered[e[0]]["amp_deg"] for e in expected]
    width = 0.4
    ax.bar(x - width / 2, exp_amp, width, label="Expected (019 FFT + Track B)")
    ax.bar(x + width / 2, rec_amp, width, label="Recovered by harmonic regression")
    ax.set_xticks(x)
    ax.set_xticklabels([e[0] for e in expected], rotation=45)
    ax.set_ylabel("Amplitude (deg)")
    ax.set_yscale("log")
    ax.set_title("Exp 020: Synthetic estimator test - harmonic amplitude recovery")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "f4_harmonic_amplitudes_synthetic.png", dpi=120)
    plt.close(fig)


def f5_estimator_bias_synthetic(payload: dict) -> None:
    """Bar chart of estimator bias on synthetic oracle (Track 3 calibration)."""
    synth = payload["synthetic_estimator_test"]
    fig, ax = plt.subplots(figsize=(7, 5))
    labels = ["Direct OLS (a)", "Harmonic regr. (f)"]
    biases = [
        abs(synth["estimator_a_bias_deg_day"]),
        abs(synth["estimator_f_bias_deg_day"]),
    ]
    ax.bar(labels, biases, color=["C0", "C1"])
    ax.set_ylabel("|Bias| (deg/day)")
    ax.set_yscale("log")
    ax.set_title("Exp 020: Synthetic estimator test - bias on known secular")
    ax.grid(True, alpha=0.3, which="both")
    for i, b in enumerate(biases):
        ax.text(i, b * 1.1, f"{b:.2e}", ha="center", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "f5_estimator_bias_synthetic.png", dpi=120)
    plt.close(fig)


def main() -> None:
    payload = load_results()
    f1_estimator_comparison(payload)
    f2_ratio_to_corrected_cf(payload)
    f3_phase_dependence(payload)
    f4_harmonic_amplitudes_recovered(payload)
    f5_estimator_bias_synthetic(payload)
    print("[020-figs] all figures written to", FIG_DIR)


if __name__ == "__main__":
    main()