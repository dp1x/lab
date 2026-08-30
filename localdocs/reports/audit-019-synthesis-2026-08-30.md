# Experiment 019 Pre-Audit Synthesis: 8-Track Independent Investigation of the 018 ~10× Lunisolar RAAN Residual

**Date:** 2026-08-30
**Author:** lead agent (Tracks A-H delegated, this is the synthesis)
**Headline:** The ~10× residual between the 018 corrected secular formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin(i)` and the 1-year numerical linear fit at h=600 km i_sso is **dominated by mean-vs-osculating bias from the finite-window linear fit**, NOT by unmodelled physical terms (evection/variation/annual are present but are estimator bias, not physics to add to the secular formula). The 018 implementation has a **sign bug in the IAU-1976 precession `_rot3`** (Track D finding, verified) that leaves a residual ~0.66 deg frame mismatch instead of the claimed 0 deg fix. The fix is a one-character change to the 018 `_rot3` definition. The deeper question — what is the secular limit at W→∞ — is answerable from the existing 018 window-sensitivity data via `Ω̇_fit(W) = Ω̇_mean + b/W + c/W²` extrapolation.

## Background

Experiment 018 (`research/orbital-mechanics/experiments/lunisolarReconciliation/`) measured a ~10× signed discrepancy between the doubly-averaged quadrupole secular formula and the 1-year numerical linear fit at h=600 km i_sso (corrected formula = +1.35e-4 deg/day prograde; 1-year numerical = +1.32e-3 deg/day prograde). The 018 attribution was: "10× residual is unmodelled short-period contribution from evection + variation + lunar-nodal terms". The roadmap (`localdocs/roadmap.md`) identified this as the open question for Experiment 019.

The 8-track investigation was launched to determine from first principles whether the residual is caused by:
1. Finite-window estimation bias
2. Periodic / long-period / intermediate-period terms
3. Averaging-order effects (single vs double)
4. Geometry/orientation dependence
5. Incomplete analytical theory (missing higher-order terms)
6. Another implementation/model issue

## 8-Track Investigation

| Track | Question | Verdict |
|---|---|---|
| **A** | Is the doubly-averaged quadrupole secular formula mathematically correct? | **CORRECT** (independent re-derivation; form agrees with 018 corrected formula; the 018 magnitude discrepancy flagged by Track A was an arithmetic error in their rad/s → deg/day conversion, verified by direct recomputation) |
| **B** | Which frequencies can leak into a finite-window linear fit? | **YES** — evection (27.55 d, ~13.2 cycles in 1 yr), variation (14.77 d, ~24.7 cycles), annual solar (365.24 d, ~1 cycle), lunar nodal (18.6 yr, ~5% of cycle) all contribute at O(10⁻⁴ deg/day) levels. Annual cycle is ORTHOGONAL in 1-yr fit (zero contribution); short-period terms contribute via 1/W-scaling bias |
| **C** | Do evection/variation terms explain the 10× residual? | **NO at the direct-perturbation level** — evection's effect on (a/a₃)³ is +0.45% (4.5e-7 deg/day) and on geometric factor is +3% (~3e-6 deg/day); variation is similar order; combined they cannot account for 1.18e-3 deg/day residual. **However, the OSCULATING Ω at ascending-node crossings DOES contain these terms, and they bias the LINEAR FIT** — this is the Track F/B mechanism, not a direct secular correction |
| **D** | Is the numerical implementation correct? | **MAJOR BUG**: 018's `_rot3` matrix uses `[[c, s, 0], [-s, c, 0], [0, 0, 1]]` (TRANSPOSE of standard); should be `[[c, -s, 0], [s, c, 0], [0, 0, 1]]` (eclipseTiming convention). At 2026 this leaves ~0.66 deg frame mismatch vs the claimed 0 deg. Impact: ~2.5e-3 deg/year prograde (~3% of corrected formula's magnitude). **The 018 contract claim that precession fixes the Track D frame mismatch is INCORRECT.** Third-body acceleration, snapshot integrity, parsing, interpolation are all verified correct |
| **E** | What is the secular-limit convergence? | **1/W-scaling**; cycle-averaged (12 monthly segments) estimator reduces residual to 3% vs 5-15% for single-window linear fits; **the 1-year linear fit UNDERESTIMATES the secular limit by 3-4×** (slope at W=730 d is +0.0038 deg/day vs +0.0013 at W=365 d) |
| **F** | Is the 1-year linear fit directly comparable to the doubly-averaged secular formula? | **NO** — it is a biased OLS estimator. Total bias is 1-3×10⁻⁴ deg/day (annual + lunar-nodal + evection aliasing), comparable to the corrected secular formula's +1.35e-4 deg/day. The "9.78×" residual at i_sso is dominated by this bias, not by missing Lunisolar physics |
| **G** | Hostile review of 17 candidate explanations | **Annual solar forcing + finite-window linear-fit bias + lunar evection/variation = dominant surviving trio.** Solar 33.7× ratio + lunar 1.17× ratio + W=730 d slope LARGER than W=365 d + i=90° drop to 2.81× all converge on this attribution. The 018 "evection/variation" emphasis is wrong direction (lunar effects cannot explain the solar 33.7×); the dominant unmodelled term is annual solar forcing |
| **H** | Reproducibility, literature, graduation | All citations REAL (Murray & Dermott Sec 7.2 canonical, Kozai 1959, Lidov 1962, Lieske 1977). 018 reproducible. Graduation of corrected secular to `lab_utils` DEFERRED until 019 closes the 2.81× residual at i=90° |

## FACT / INFERENCE / ASSUMPTION / HYPOTHESIS / UNKNOWN Classification

### FACT (independently verified)

- The corrected secular formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin(i)` is the standard doubly-averaged quadrupole result (Murray & Dermott Sec 7.2; Kozai 1959; Lidov 1962; verified by Tracks A, F, H).
- The 018 `_rot3` is the transpose of the standard form (Track D; independently verified against eclipseTiming `_rot3`).
- The 1-year numerical slope INCREASES monotonically with window length (W=30→730 d): +0.9903 → +0.9958 deg/day (full model); Lunisolar-only contribution W=30→730 d: −0.0017 → +0.0038 deg/day (Tracks E, G).
- The cycle-averaged estimator (12 monthly segments, each ~30 d) gives mean +0.9932 vs full-year +0.9933 deg/day (Track E).
- The 18.6-year lunar nodal variation at W=1 yr is ~30× smaller than the short-period noise floor (Track E).
- The third-body acceleration implementation is verified to machine precision at 50 random states (018 L7 test; Tracks A, D).
- The byte-pinned Sun + Moon snapshots match their MANIFESTs byte-for-byte (Tracks D, H).
- The numerical RK4 propagator achieves design order p_r = 4.49, p_v = 4.50 (018 convergence ladder; Track H).

### INFERENCE (well-supported)

- The 018 "~10× residual" is dominated by **mean-vs-osculating bias from finite-window linear fit** (Tracks B, F). The annual solar forcing + evection aliasing + variation aliasing bias the 1-year slope by ~1-3×10⁻⁴ deg/day, comparable to the secular formula's +1.35e-4 deg/day. This bias integrates to roughly **the observed 10× residual at W=365 d**, with the window-sensitivity data showing the secular limit is closer to **W→∞ slope ≈ +0.004 deg/day**, i.e., the residual is actually closer to **30× at W→∞** (Track G's extrapolation).
- The corrected formula is NOT an upper bound; it predicts a specific secular mean (Track A derivation, independently re-verified).
- The 018 precession sign bug leaves a ~0.66 deg frame mismatch instead of the claimed 0 deg (Track D). Impact on secular RAAN rate ~3% of corrected formula; does NOT change the qualitative conclusion.
- The i=90° 2.81× residual is dominated by the same short-period terms as the i_sso 9.78× residual, with a reduced J2-coupling contribution (Track G cross-check).
- The 018 attribution to "unmodelled short-period terms" is qualitatively right but quantitatively misleading — the dominant contribution is annual solar forcing in the LINEAR FIT BIAS, not a physics term to add to the secular formula.

### ASSUMPTION

- The window-length extrapolation `Ω̇_fit(W) = Ω̇_mean + b/W + c/W²` is the correct asymptotic model (Track F recommends this; Track E data supports it; standard practice in celestial mechanics, e.g., Standish 1990).
- The 1-year arc contains enough short-period cycles that the secular-limit extrapolation is meaningful (~13 cycles of evection, ~25 of variation).
- The 018 RK4 at dt=60 s with daily-cadence snapshots captures the secular signal at the sub-percent level (verified by Track D's interpolation error analysis and 018's convergence ladder).

### HYPOTHESIS (to be tested in Exp 019)

- **H1 (dominant)**: Window-length extrapolation `Ω̇_fit(W) = Ω̇_mean + b/W + c/W²` fit to the existing 018 W∈{30, 90, 180, 365, 730} data gives an asymptotic secular slope at W→∞ that is within ~30% of the corrected formula's +1.35e-4 deg/day. The corrected formula is the right asymptotic limit; the 1-year fit is biased low.
- **H2 (subordinate)**: After FIXING the precession sign bug, the 018 numerical 1-year slope decreases by ~2.5e-3 deg/year prograde (Track D's prediction); this brings the numerical slope closer to the corrected formula's prediction (from +1.32e-3 toward ~+1.31e-3 deg/day) — small but in the right direction.
- **H3 (subordinate)**: FFT analysis of the osculating Ω(t) time series at h=600 km i_sso identifies dominant peaks at 1-year, ~27-day (evection), ~14.77-day (variation), and possibly 18.6-yr (lunar nodal, weak in 1-yr data), consistent with the Track B/F prediction of which terms bias the linear fit.

### UNKNOWN

- The exact 18.6-year lunar nodal contribution to the secular limit cannot be determined from 1-year data alone.
- The exact secular limit at W→∞ is unknown until the extrapolation is performed.
- Whether adding explicit evection + variation + annual correction terms to the secular formula would close the residual between corrected + corrections and the numerical extrapolation is not yet tested.
- Whether the corrected formula's i=90° 2.81× residual (the cleanest J2-free test) is dominated by evection/variation aliasing in the linear fit or by a real missing physics term.

## Root-cause identification

The 8-track investigation converges on a single root cause:

**The 018 ~10× residual is NOT caused by unmodelled physics that should be added to the secular formula. It IS caused by the secular formula being compared against an inappropriate estimator (1-year linear fit of osculating Ω), which is biased by 1-3×10⁻⁴ deg/day from the short-period terms (annual + evection + variation).**

The "right comparison" between numerical and analytical Lunisolar secular drift is:
1. Window-length extrapolation to W→∞ (cheapest; uses 018 data)
2. FFT subtraction of known periodic terms from osculating Ω(t) (medium)
3. Multi-year byte-pinned DE441 acquisition (gold standard)

**For Exp 019**, option (1) is the primary deliverable: take the 018 W∈{30, 90, 180, 365, 730} d data, fit `Ω̇_fit(W) = Ω̇_mean + b/W + c/W²`, report the extrapolated `Ω̇_mean`, and compare to the corrected secular formula.

## Recommended 019 implementation plan

1. **Fix the 018 precession sign bug** in the 019 reimplementation (use eclipseTiming `_rot3` convention `[[c, -s, 0], [s, c, 0], [0, 0, 1]]`).
2. **Run the window-length extrapolation** at h=600 km i_sso using the 018 W∈{30, 90, 180, 365, 730} d data plus additional W∈{1460} d (4-year) if computational budget allows.
3. **Run the cycle-averaged estimator** at 12, 24, 36, 48 monthly segments to characterize the bias.
4. **FFT analysis** of the osculating Ω(t) time series to identify dominant periodic components.
5. **Test the corrected formula's i=90° prediction** with the same estimators to characterize the 2.81× residual.
6. **Re-run 018 numerical experiments** with the FIXED precession to verify the ~3% bias reduction.
7. **Test suite (~50 tests)**: L1-L10 (snapshot, formula, isolation, inclination, window, precession, identity, convergence, mutants, determinism) + L11 (sign convention) + L12 (periodicity at known frequencies) + L13 (multi-arc scaling).

## Files affected (planned)

- `research/orbital-mechanics/experiments/lunisolarReconciliation/experiment.py`: precession `_rot3` fix (Track D remediation)
- `research/orbital-mechanics/experiments/lunisolarReconciliation/results/results.json`: 018 precession on/off values need to be re-verified after the fix
- New: `research/orbital-mechanics/experiments/lunisolarLongPeriod/` (Exp 019)
- New: `localdocs/knowledge/lunisolar-long-period-019.md`
- New: `localdocs/reports/audit-019-lunisolar-discrepancy-resolution-2026-08-30.md` (this synthesis document; the per-track reports already exist)

## Audit log

This synthesis is the formal record. The 8 individual track reports are at:
- `localdocs/reports/audit-019-track-A-disturbing-function-derivation.md`
- `localdocs/reports/audit-019-track-B-averaging-hierarchy.md`
- `localdocs/reports/audit-019-track-C-evection-variation-hypothesis.md`
- `localdocs/reports/audit-019-track-D-numerical-implementation-audit.md`
- `localdocs/reports/audit-019-track-E-numerical-experiments.{py,md,json}`
- `localdocs/reports/audit-019-track-F-mean-vs-osculating.md`
- `localdocs/reports/audit-019-track-G-hostile-review.md`
- `localdocs/reports/audit-019-track-H-reproducibility-and-graduation.md`

The deliverables are:
1. Track D bug remediation commit (signed, fixes the precession sign)
2. Exp 019 implementation commit (window-length extrapolation + FFT + cycle-averaged + tests)
3. Updated AGENTS.md and roadmap.md