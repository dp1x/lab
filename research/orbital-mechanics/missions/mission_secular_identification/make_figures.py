"""Figures for mission_secular_identification.

Three panels, each answering one question:

1. WINDOW LADDER -- is the estimate converging, or still drifting with the arc?
   This is the identifiability question, plotted directly.
2. MEASURED vs BOTH FORMULAS -- the adjudication, per inclination.
3. IDENTIFIABILITY DIAGNOSTICS -- VIF and phase spread against the gate
   thresholds, so a reader sees which inclinations were scoreable at all.

Deterministic: no randomness, fixed input file. Run after campaign.py.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402

import competing_formulas as cf                                    # noqa: E402
from adjudicate import _stability                                 # noqa: E402

HERE = Path(__file__).resolve().parent
FIG = HERE / "results" / "figures"
A_SSO = 6978.137


def fig_window_ladder(campaign: dict) -> Path:
    cases = campaign["cases"]
    incs = sorted({c["inc_deg"] for c in cases})
    phases = sorted({c["phase_deg"] for c in cases})
    fig, axes = plt.subplots(1, len(incs), figsize=(4.2 * len(incs), 3.6),
                             sharey=True)
    axes = [axes] if len(incs) == 1 else list(axes)
    for ax, inc in zip(axes, incs):
        rows = [c for c in cases if c["inc_deg"] == inc]
        any_ladder = False
        for ph in phases:
            row = next((r for r in rows if r["phase_deg"] == ph), None)
            if not row or not row.get("window_ladder"):
                continue
            lad = row["window_ladder"]
            ks = sorted(lad, key=float)
            xs = [float(k) for k in ks]
            ys = [lad[k]["lunisolar_deg_day"] for k in ks]
            ax.plot(xs, ys, "o-", ms=3, lw=1, alpha=0.75, label=f"phase {ph:g}$\\degree$")
            any_ladder = True
        f1 = cf.form1_audit018(A_SSO, inc)
        ax.axhline(f1, color="C0", ls="--", lw=1.2, label="FORM-1")
        f2 = cf.form2_competitor(A_SSO, inc)
        ax.axhline(f2, color="C1", ls=":", lw=1.2, label="FORM-2")
        ax.axhline(0.0, color="k", lw=0.5, alpha=0.4)
        ax.set_xlabel("arc length (lunar nodal cycles)")
        ax.set_title(f"i = {inc:.3f}$\\degree$", fontsize=10)
        ax.grid(alpha=0.25)
        if not any_ladder:
            ax.text(0.5, 0.5, "ladder unavailable", ha="center", transform=ax.transAxes)
    axes[0].set_ylabel("lunisolar $\\dot\\Omega$ (deg/day)")
    axes[0].legend(fontsize=7, loc="best", ncol=2)
    fig.suptitle("Window ladder: does the secular rate converge with arc length?", y=1.0)
    p = FIG / "fig1_window_ladder.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_measured_vs_formulas(campaign: dict, verdict: dict) -> Path:
    cases = campaign["cases"]
    per = verdict["per_inclination"]
    incs = sorted({c["inc_deg"] for c in cases})
    phases = sorted({c["phase_deg"] for c in cases})
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = list(range(len(incs)))
    w = 0.2
    for k, ph in enumerate(phases):
        ys = []
        for inc in incs:
            row = next((r for r in cases
                        if r["inc_deg"] == inc and r["phase_deg"] == ph), None)
            ys.append(row["lunisolar"]["secular_deg_day"] if row else math.nan)
        ax.errorbar([xi + (k - 1.5) * w for xi in x], ys,
                    yerr=[0.0] * len(ys), fmt="o", ms=4, alpha=0.55,
                    label=f"measured (phase {ph:g}$\\degree$)")
    ax.plot(x, [cf.form1_audit018(A_SSO, i) for i in incs], "C0--", lw=2,
            marker="s", ms=7, label="FORM-1 (audit-018)")
    ax.plot(x, [cf.form2_competitor(A_SSO, i) for i in incs], "C1:", lw=2,
            marker="^", ms=7, label="FORM-2 (competitor)")
    for xi, inc in zip(x, incs):
        e = per.get(str(inc))
        if e and not e["identifiable"]:
            ax.annotate("not identifiable", (xi, 0), ha="center", va="bottom",
                        fontsize=7, color="C3", rotation=90)
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{i:.3f}" for i in incs])
    ax.set_xlabel("inclination (deg)")
    ax.set_ylabel("lunisolar $\\dot\\Omega$ (deg/day)")
    ax.set_title(f"Measured secular rate vs both formulas "
                 f"({campaign['arc']['headline_cycles']:g} nodal cycles)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.25)
    p = FIG / "fig2_measured_vs_formulas.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_identifiability(campaign: dict, verdict: dict) -> Path:
    per = verdict["per_inclination"]
    incs = sorted(float(k) for k in per)
    cases = campaign["cases"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    vifs = [per[str(i)]["vif_max"] for i in incs]
    ax1.bar([str(i) for i in incs], vifs, color=["C0" if per[str(i)]["vif_ok"] else "C3"
                                                for i in incs])
    ax1.axhline(1.20, color="C3", ls="--", lw=1.4, label="gate: VIF $\\leq$ 1.20")
    ax1.axhline(1.0, color="k", lw=0.6, alpha=0.4)
    ax1.set_yscale("log")
    ax1.set_xlabel("inclination (deg)")
    ax1.set_ylabel("$\\max$ VIF")
    ax1.set_title("Identifiability gate: variance inflation")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.25, axis="y")

    spreads = [per[str(i)]["phase_spread_rel"] * 100 for i in incs]
    ax2.bar([str(i) for i in incs], spreads,
            color=["C0" if per[str(i)]["sign_agreement"] else "C3" for i in incs])
    ax2.axhline(100.0, color="C3", ls="--", lw=1.4, label="gate: spread < 100 %")
    ax2.set_xlabel("inclination (deg)")
    ax2.set_ylabel("phase-to-phase spread (% of median)")
    ax2.set_title("Phase decorrelation")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.25, axis="y")
    fig.suptitle("Which inclinations were scoreable at all?", y=1.02)
    p = FIG / "fig3_identifiability.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_referee(verdict: dict) -> Path:
    ref = verdict.get("quadrature_referee") or {}
    incs = [float(k) for k in ref if k != "error"]
    if not incs:
        return None
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    x = list(range(len(incs)))
    ax.plot(x, [ref[str(i)]["referee_deg_day"] for i in incs], "o-", color="k",
            lw=2, ms=7, label="quadrature referee (independent)")
    ax.plot(x, [ref[str(i)]["form1_deg_day"] for i in incs], "C0--", lw=1.6,
            marker="s", ms=6, label="FORM-1")
    ax.plot(x, [ref[str(i)]["form2_deg_day"] for i in incs], "C1:", lw=1.6,
            marker="^", ms=6, label="FORM-2")
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{i:g}" for i in incs])
    ax.set_xlabel("inclination (deg)")
    ax.set_ylabel("secular $\\dot\\Omega$ (deg/day)")
    ax.set_title("Independent referee: exact-potential double average vs both formulas")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
    p = FIG / "fig4_referee.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    campaign = json.loads((HERE / "results" / "campaign_2cyc.json").read_text(encoding="utf-8"))
    verdict = json.loads((HERE / "results" / "adjudication.json").read_text(encoding="utf-8"))
    out = [fig_window_ladder(campaign), fig_measured_vs_formulas(campaign, verdict),
           fig_identifiability(campaign, verdict)]
    r = fig_referee(verdict)
    if r:
        out.append(r)
    for p in out:
        print("wrote", p)


if __name__ == "__main__":
    main()