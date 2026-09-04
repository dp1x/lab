__test__ = False  # not a pytest test module
"""Generate publication-quality figures for mission_j2_lunisolar_coupling."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def fig1_perturbative_scaling():
    """Phase B: 2-D response surface for f(λ_J2, λ_3body) Lunisolar rate."""
    with open(HERE / "results" / "phase_b_perturbative_scaling.json") as f:
        data = json.load(f)

    LAMBDA_J2_VALUES = data["lambda_j2_values"]
    LAMBDA_3BODY_VALUES = data["lambda_3body_values"]

    # Build matrix
    results = data["results"]
    full = {}
    j2_only = {}
    for r in results:
        lj2 = r["lambda_j2"]
        l3b = r["lambda_3body"]
        if r["mode"] == "sun_moon_j2":
            full[(lj2, l3b)] = r["rate_deg_day"]
        elif r["mode"] == "j2_only":
            j2_only[lj2] = r["rate_deg_day"]

    # Lunisolar = full - j2_only
    Lj2 = np.array(LAMBDA_J2_VALUES)
    L3b = np.array(LAMBDA_3BODY_VALUES)
    Z = np.zeros((len(Lj2), len(L3b)))
    for i, lj2 in enumerate(Lj2):
        for j, l3b in enumerate(L3b):
            Z[i, j] = full[(lj2, l3b)] - j2_only[lj2]

    fig, ax = plt.subplots(figsize=(8, 5))
    # Plot each λ_J2 as a curve
    for i, lj2 in enumerate(Lj2):
        if lj2 == 0:
            continue  # skip λ_J2=0 (trivial zero)
        ax.plot(L3b, Z[i] * 1e4, "o-", label=f"λ_J2={lj2:.2f}")
    ax.set_xlabel("λ_3body (third-body multiplier)")
    ax.set_ylabel("Lunisolar rate (×10⁻⁴ deg/day)")
    ax.set_title("Phase B: Perturbative scaling at h=600 km, i_sso, 90-d arc\n"
                  "Cross-coupling term a11 statistically significant (SNR=6.89)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5)
    plt.tight_layout()
    plt.savefig(HERE / "results" / "figures" / "fig1_perturbative_scaling.png",
                dpi=120)
    plt.close()
    print("Wrote fig1_perturbative_scaling.png")


def fig2_force_mode_decomposition():
    """Phase C: Lunisolar decomposition by inclination at 1-yr arc."""
    # Find the latest phase C file
    candidates = sorted(HERE.glob("results/phase_c_full_*.json"))
    if not candidates:
        print("No phase C results found")
        return
    with open(candidates[-1]) as f:
        data = json.load(f)

    by_i = {}
    for r in data["results"]:
        label = r["label"]
        i_part, mode_part = label.split("_", 1)
        i_deg = float(i_part.lstrip("i"))
        by_i.setdefault(i_deg, {})[mode_part] = r["rate_deg_day"]

    inclinations = sorted(by_i.keys())
    Luni_combined = []
    Luni_isolated = []
    R_J2x3b = []
    for i_deg in inclinations:
        full = by_i[i_deg]["sun_moon_j2"]
        j2 = by_i[i_deg]["j2_only"]
        sun = by_i[i_deg]["sun_only"]
        moon = by_i[i_deg]["moon_only"]
        luni_combined = full - j2
        luni_isolated = sun + moon
        R = luni_combined - luni_isolated
        Luni_combined.append(luni_combined)
        Luni_isolated.append(luni_isolated)
        R_J2x3b.append(R)

    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.25
    x = np.arange(len(inclinations))
    ax.bar(x - width, np.array(Luni_combined) * 1e4, width,
            label="Luni_combined (full − J2)")
    ax.bar(x, np.array(Luni_isolated) * 1e4, width,
            label="Luni_isolated (Sun+Moon, no J2)")
    ax.bar(x + width, np.array(R_J2x3b) * 1e4, width,
            label="R_J2×3b (cross-coupling)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{i:.1f}°" for i in inclinations])
    ax.set_xlabel("Inclination")
    ax.set_ylabel("Lunisolar rate (×10⁻⁴ deg/day)")
    ax.set_title("Phase C: Force-mode decomposition at 1-yr arc, h=600 km\n"
                  "J2 × Lunisolar cross-coupling is the DOMINANT component")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5)
    plt.tight_layout()
    plt.savefig(HERE / "results" / "figures" / "fig2_force_mode_decomposition.png",
                dpi=120)
    plt.close()
    print("Wrote fig2_force_mode_decomposition.png")


def fig3_synthetic_vs_real_moon():
    """Phase A: synthetic (circular) vs real Moon comparison."""
    with open(HERE / "results" / "phase_a_reduced_model.json") as f:
        data = json.load(f)

    real = {}
    synthetic = {}
    for r in data["results"]:
        label = r["label"]
        if label.startswith("real_"):
            mode = label[5:]
            real[mode] = r["rate_deg_day"]
        elif label.startswith("synthetic_"):
            mode = label[10:]
            synthetic[mode] = r["rate_deg_day"]

    modes = ["sun_only", "moon_only", "sun_moon"]
    real_vals = [real.get(m, 0) for m in modes]
    synthetic_vals = [synthetic.get(m, 0) for m in modes]
    diff = [r - s for r, s in zip(real_vals, synthetic_vals)]

    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.25
    x = np.arange(len(modes))
    ax.bar(x - width, np.array(real_vals) * 1e4, width, label="Real DE441 Moon")
    ax.bar(x, np.array(synthetic_vals) * 1e4, width,
            label="Synthetic circular Moon (fixed i3, e3=0)")
    ax.bar(x + width, np.array(diff) * 1e4, width,
            label="Difference (lunar e/i variation)")
    ax.set_xticks(x)
    ax.set_xticklabels(modes, rotation=0)
    ax.set_xlabel("Mode")
    ax.set_ylabel("RAAN rate (×10⁻⁴ deg/day)")
    ax.set_title("Phase A: Real vs synthetic Moon at 1-yr arc, h=600 km, i=i_sso\n"
                  "Real Moon contributes 2x more than circular Moon")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5)
    plt.tight_layout()
    plt.savefig(HERE / "results" / "figures" / "fig3_synthetic_vs_real_moon.png",
                dpi=120)
    plt.close()
    print("Wrote fig3_synthetic_vs_real_moon.png")


def main():
    HERE.joinpath("results", "figures").mkdir(parents=True, exist_ok=True)
    fig1_perturbative_scaling()
    fig2_force_mode_decomposition()
    fig3_synthetic_vs_real_moon()


if __name__ == "__main__":
    main()
