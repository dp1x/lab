__test__ = False  # not a pytest test module
"""Analyze Phase C results: 18.6-yr force-mode decomposition.

Computes the non-additive residual R_J2x3b = full - J2 - Sun - Moon + 2-body
at 3 inclinations, and tests whether it can be explained by:
  - J2 x Lunisolar coupling (H1)
  - Higher-order Lunisolar terms (H0b)
  - Lunar eccentricity / inclination variation (H0c)
  - Mean-vs-osculating bias (H0d)
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main():
    # Find the latest phase C results file
    candidates = sorted(HERE.glob("results/phase_c_full_*.json"))
    if not candidates:
        print("No phase C results found")
        return
    phase_c_path = candidates[-1]
    print(f"Analyzing: {phase_c_path}")
    with open(phase_c_path) as f:
        data = json.load(f)

    results = data["results"]
    # Group by inclination and mode
    by_i = {}
    for r in results:
        # Label format: "i{deg}_{mode}"
        label = r["label"]
        i_part, mode_part = label.split("_", 1)
        i_deg = float(i_part.lstrip("i"))
        mode = mode_part
        by_i.setdefault(i_deg, {})[mode] = r

    print("=== Phase C: 18.6-yr force-mode decomposition ===")
    print(f"Arc: {data['arc_days']:.2f} d, h=600 km, 3 inclinations")
    print()

    inclinations = sorted(by_i.keys())
    modes = ["kepler_only", "j2_only", "sun_only", "moon_only",
             "sun_moon", "sun_moon_j2"]

    print("Per-mode 18.6-yr OLS RAAN rate (deg/day):")
    print(f"{'i (deg)':>10s} | ", end="")
    for m in modes:
        print(f"{m:>15s}", end=" | ")
    print()
    for i_deg in inclinations:
        print(f"{i_deg:>10.2f} | ", end="")
        for m in modes:
            r = by_i[i_deg].get(m)
            if r:
                print(f"{r['rate_deg_day']:>+15.6e}", end=" | ")
            else:
                print(f"{'N/A':>15s}", end=" | ")
        print()

    print()
    print("=== Force-mode decomposition (Lunisolar component) ===")
    print("Luni_combined = sun_moon_j2 - j2_only   (the COMBINED Lunisolar component, includes J2 x 3b coupling)")
    print("Luni_isolated = sun_only + moon_only      (the SUM of individual Lunisolar contributions, NO J2)")
    print("Cross-coupling R = Luni_combined - Luni_isolated = J2 x 3b coupling signal")
    print()
    print(f"{'i (deg)':>10s} | {'Luni_combined':>15s} | {'Luni_isolated':>15s} | {'Cross_R':>15s} | {'R/Luni_combined%':>15s}")
    for i_deg in inclinations:
        full = by_i[i_deg]["sun_moon_j2"]["rate_deg_day"]
        j2 = by_i[i_deg]["j2_only"]["rate_deg_day"]
        sun = by_i[i_deg]["sun_only"]["rate_deg_day"]
        moon = by_i[i_deg]["moon_only"]["rate_deg_day"]
        kep = by_i[i_deg]["kepler_only"]["rate_deg_day"]
        luni_combined = full - j2
        luni_isolated = sun + moon
        cross_R = luni_combined - luni_isolated
        R_pct = 100 * cross_R / luni_combined if abs(luni_combined) > 1e-30 else float("nan")
        print(f"{i_deg:>10.2f} | {luni_combined:>+15.6e} | {luni_isolated:>+15.6e} | {cross_R:>+15.6e} | {R_pct:>+14.1f}%")

    print()
    print("=== Non-additive residual R_J2x3b ===")
    print("R = sun_moon_j2 - j2_only - sun_only - moon_only + kepler_only")
    print("(the cross-coupling residual; should be ~0 if perturbations are additive)")
    print()
    for i_deg in inclinations:
        full = by_i[i_deg]["sun_moon_j2"]["rate_deg_day"]
        j2 = by_i[i_deg]["j2_only"]["rate_deg_day"]
        sun = by_i[i_deg]["sun_only"]["rate_deg_day"]
        moon = by_i[i_deg]["moon_only"]["rate_deg_day"]
        kep = by_i[i_deg]["kepler_only"]["rate_deg_day"]
        luni_combined = full - j2
        luni_isolated = sun + moon
        R_J2x3b = luni_combined - luni_isolated
        R_pct = 100 * R_J2x3b / luni_combined if abs(luni_combined) > 1e-30 else float("nan")
        print(f"  i={i_deg:6.2f}: R_J2x3b = {R_J2x3b:+.4e} deg/day "
              f"({R_pct:+.1f}% of combined Lunisolar)")

    print()
    print("=== Verdict on J2 x Lunisolar coupling at this arc length ===")
    print("Decision rule conditions from README §4.1:")
    print("  (a) residual > 1e-4 deg/day")
    print("  (b) residual scales with lambda_J2 * lambda_3body (Phase B)")
    print("  (c) sign compatible with observed retrograde at i_sso")
    print("  (d) residual at i=90 smaller than at i_sso (J2-modulated)")
    print("  (e) residual survives convergence ladder")

    for i_deg in inclinations:
        full = by_i[i_deg]["sun_moon_j2"]["rate_deg_day"]
        j2 = by_i[i_deg]["j2_only"]["rate_deg_day"]
        sun = by_i[i_deg]["sun_only"]["rate_deg_day"]
        moon = by_i[i_deg]["moon_only"]["rate_deg_day"]
        luni_combined = full - j2
        luni_isolated = sun + moon
        R_J2x3b = luni_combined - luni_isolated
        R_pct = 100 * R_J2x3b / luni_combined if abs(luni_combined) > 1e-30 else float("nan")
        print(f"  i={i_deg:6.2f}: R_J2x3b = {R_J2x3b:+.4e} deg/day ({R_pct:+.1f}% of Luni_combined)")

    # Save the analysis to a JSON
    out = {
        "phase": "C_analysis",
        "arc_days": data["arc_days"],
        "by_inclination": {},
    }
    for i_deg in inclinations:
        full = by_i[i_deg]["sun_moon_j2"]["rate_deg_day"]
        j2 = by_i[i_deg]["j2_only"]["rate_deg_day"]
        sun = by_i[i_deg]["sun_only"]["rate_deg_day"]
        moon = by_i[i_deg]["moon_only"]["rate_deg_day"]
        kep = by_i[i_deg]["kepler_only"]["rate_deg_day"]
        sm = by_i[i_deg]["sun_moon"]["rate_deg_day"]
        luni_combined = full - j2
        luni_isolated = sun + moon
        R_J2x3b = luni_combined - luni_isolated
        out["by_inclination"][str(i_deg)] = {
            "kepler_rate_deg_day": kep,
            "j2_only_rate_deg_day": j2,
            "sun_only_rate_deg_day": sun,
            "moon_only_rate_deg_day": moon,
            "sun_moon_rate_deg_day": sm,
            "sun_moon_j2_rate_deg_day": full,
            "Luni_combined_deg_day": luni_combined,
            "Luni_isolated_sun_plus_moon_deg_day": luni_isolated,
            "R_J2x3b_cross_residual_deg_day": R_J2x3b,
            "R_J2x3b_pct_of_luni_combined": (
                100 * R_J2x3b / luni_combined if abs(luni_combined) > 1e-30 else None
            ),
        }
    with open(HERE / "results" / "phase_c_analysis.json", "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nAnalysis written to {HERE / 'results' / 'phase_c_analysis.json'}")


if __name__ == "__main__":
    main()
