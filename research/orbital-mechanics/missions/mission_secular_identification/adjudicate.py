"""Adjudication of the campaign against the pre-registered decision rule.

This module READS the campaign output and applies the mission card's §4 rule.
It performs no estimation and no propagation: every number it judges was
produced elsewhere. Its job is to make the verdict mechanical rather than a
matter of narrative emphasis.

The rule (README §4, fixed before any result existed):

  4.1 Identifiability gate -- must pass before ANY formula is scored:
       VIF_max <= 1.20, AND stability across the window ladder, AND sign
       agreement across phases.
  4.2 Per-inclination scoring: same sign and |measured/formula| in [0.5, 2.0].
  4.3 Structural discriminator (across the inclination set, decisive):
       equal-and-opposite twin pair within 15 %      -> FORM-2
       same sign with |lo/hi| in [1.2, 2.0]         -> FORM-1
       otherwise                                    -> neither
  4.4 Magnitude verdict, only if 4.3 picked a winner.

A formula is NOT scored at an inclination that fails 4.1. Failing 4.1 is recorded
as H-null at that inclination, never as evidence for a formula.
"""
from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
import competing_formulas as cf                                    # noqa: E402

A_SSO = 6978.137
VIF_MAX = 1.20
TWIN_HI, TWIN_LO = 97.7876, 82.2124
STABILITY_FRAC = 0.25          # ladder spread allowed, as a fraction of the value


def _stability(case: dict) -> dict:
    """Does the estimate stop moving as the window grows?

    Identifiability means the LAST TWO ladder rungs agree. The full spread is
    also reported, because a value that drifts slowly across the whole ladder is
    still not a converged secular rate even if its tail happens to be flat.
    """
    lad = case.get("window_ladder") or {}
    keys = sorted(lad, key=lambda k: float(k))
    vals = [lad[k]["lunisolar_deg_day"] for k in keys]
    if len(vals) < 2:
        return {"stable": False, "reason": "ladder too short"}
    tail = abs(vals[-1] - vals[-2])
    spread = max(vals) - min(vals)
    ref = max(abs(vals[-1]), 1e-30)
    return {
        "stable": bool(tail / ref < STABILITY_FRAC),
        "tail_drift_rel": tail / ref,
        "full_spread_rel": spread / ref,
        "sign_changes": sum(1 for a, b in zip(vals, vals[1:]) if a * b < 0),
        "ladder": {k: lad[k]["lunisolar_deg_day"] for k in keys},
    }


def adjudicate(campaign: dict) -> dict:
    cases = campaign["cases"]
    by_inc: dict[float, list] = {}
    for c in cases:
        by_inc.setdefault(c["inc_deg"], []).append(c)

    per_inc = {}
    for inc, rows in sorted(by_inc.items()):
        vals = [r["lunisolar"]["secular_deg_day"] for r in rows]
        vifs = [r["lunisolar"]["vif_max"] for r in rows]
        stab = [_stability(r) for r in rows]
        sign_agree = all(v > 0 for v in vals) or all(v < 0 for v in vals)
        median = statistics.median(vals)
        # phase spread relative to the median: the honest uncertainty
        spread = (max(vals) - min(vals)) / abs(median) if median else math.inf

        vif_ok = max(vifs) <= VIF_MAX
        stable_ok = all(s["stable"] for s in stab)
        ident = bool(vif_ok and stable_ok and sign_agree and spread < 1.0)

        entry = {
            "inc_deg": inc,
            "n_phases": len(rows),
            "values_deg_day": vals,
            "median_deg_day": median,
            "mean_deg_day": statistics.fmean(vals),
            "phase_spread_rel": spread,
            "sign_agreement": bool(sign_agree),
            "vif_max": max(vifs),
            "vif_ok": bool(vif_ok),
            "stability": {"all_stable": bool(stable_ok),
                          "per_phase": stab},
            "identifiable": ident,
            "nodal_conditioning": rows[0]["nodal_conditioning"],
        }
        if ident:
            f1 = cf.form1_audit018(A_SSO, inc)
            f2 = cf.form2_competitor(A_SSO, inc)
            entry["formulas"] = {
                "form1_deg_day": f1, "form2_deg_day": f2,
                "measured_over_form1": median / f1,
                "measured_over_form2": median / f2,
                "same_sign_form1": bool(median * f1 > 0),
                "same_sign_form2": bool(median * f2 > 0),
                "h1_supported": bool(median * f1 > 0 and 0.5 <= abs(median / f1) <= 2.0),
                "h2_supported": bool(median * f2 > 0 and 0.5 <= abs(median / f2) <= 2.0),
            }
        else:
            entry["formulas"] = None
            entry["scored"] = False
        entry["scored"] = ident
        per_inc[str(inc)] = entry

    # ---- 4.3 structural discriminator -------------------------------------
    struct = {"evaluable": False}
    hi = per_inc.get(str(TWIN_HI))
    lo = per_inc.get(str(TWIN_LO))
    if hi and lo and hi["scored"] and lo["scored"]:
        m_hi, m_lo = hi["median_deg_day"], lo["median_deg_day"]
        same_sign = (m_hi * m_lo) > 0
        ratio = abs(m_lo / m_hi) if m_hi else math.inf
        struct = {
            "evaluable": True,
            "measured_hi": m_hi, "measured_lo": m_lo,
            "measured_same_sign": bool(same_sign),
            "measured_ratio_lo_over_hi": ratio,
            "rule_form2_equal_and_opposite":
                (not same_sign) and abs(ratio - 1.0) <= 0.15,
            "rule_form1_same_sign_gap":
                same_sign and 1.2 <= ratio <= 2.0,
        }
        if struct["rule_form2_equal_and_opposite"]:
            struct["verdict"] = "FORM-2 supported; FORM-1 refuted"
        elif struct["rule_form1_same_sign_gap"]:
            struct["verdict"] = "FORM-1 supported; FORM-2 refuted"
        else:
            struct["verdict"] = "NEITHER: measured twin pair matches neither prediction"

    # ---- quadrature referee (independent of the orbital measurement) --------
    # Deliberately kept in the SAME verdict document. The referee is an
    # independent line of evidence about the same physical quantity, so a reader
    # must be able to see the orbital measurement and its independent check side
    # by side rather than in separate reports.
    ref = {}
    try:
        import quadrature_referee as Q
        for inc in cf.DIAGNOSTIC_INCS_DEG:
            moon = Q.quadrature_node_rate(inc, A_SSO, 384400.0, 4902.8001, 28.5843)
            sun = Q.quadrature_node_rate(inc, A_SSO, 1.4959787e8,
                                         132712440018.0, 23.4393)
            f1, f2 = cf.form1_audit018(A_SSO, inc), cf.form2_competitor(A_SSO, inc)
            ref[str(inc)] = {"referee_deg_day": moon + sun,
                             "form1_deg_day": f1, "form2_deg_day": f2,
                             "referee_over_form1": (moon + sun) / f1}
    except Exception as exc:                                    # pragma: no cover
        ref = {"error": repr(exc)}

    n_id = sum(1 for v in per_inc.values() if v["identifiable"])
    return {
        "per_inclination": per_inc,
        "structural_discriminator": struct,
        "n_identifiable_inclinations": n_id,
        "n_inclinations": len(per_inc),
        "identifiability_verdict": (
            "IDENTIFIABLE" if n_id == len(per_inc) and per_inc
            else ("PARTIALLY IDENTIFIABLE" if n_id else "NOT IDENTIFIABLE")),
        "quadrature_referee": ref,
        "evidence_tier": "E3 (cross-method: orbital measurement + independent quadrature)",
        "protocol_tag": "P2 (card + decision rule frozen before any result)",
    }


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", default=str(HERE / "results" / "campaign_2cyc.json"))
    ap.add_argument("--out", default=str(HERE / "results" / "adjudication.json"))
    args = ap.parse_args()

    campaign = json.loads(Path(args.campaign).read_text(encoding="utf-8"))
    verdict = adjudicate(campaign)
    Path(args.out).write_text(json.dumps(verdict, indent=2), encoding="utf-8")

    print(f"IDENTIFIABILITY: {verdict['identifiability_verdict']} "
          f"({verdict['n_identifiable_inclinations']}/{verdict['n_inclinations']} inclinations)")
    print(f"\n{'inc':>9} {'median luni':>15} {'VIF':>6} {'ident?':>6} "
          f"{'FORM-1':>12} {'m/f1':>8} {'FORM-2':>12} {'m/f2':>8}")
    for k, v in verdict["per_inclination"].items():
        if v["scored"]:
            f = v["formulas"]
            print(f"{v['inc_deg']:>9.3f} {v['median_deg_day']:>15.6e} "
                  f"{v['vif_max']:>6.2f} {'yes':>6} "
                  f"{f['form1_deg_day']:>12.4e} {f['measured_over_form1']:>8.3f} "
                  f"{f['form2_deg_day']:>12.4e} {f['measured_over_form2']:>8.3f}")
        else:
            print(f"{v['inc_deg']:>9.3f} {v['median_deg_day']:>15.6e} "
                  f"{v['vif_max']:>6.2f} {'NO':>6} {'-':>12} {'-':>8} "
                  f"{'-':>12} {'-':>8}")
    s = verdict["structural_discriminator"]
    print(f"\nSTRUCTURAL DISCRIMINATOR: {s.get('verdict', 'not evaluable')}")
    if s.get("evaluable"):
        print(f"  measured lo/hi ratio = {s['measured_ratio_lo_over_hi']:.4f} "
              f"(same sign: {s['measured_same_sign']})")
    print(f"\n-> {args.out}")


if __name__ == "__main__":
    main()