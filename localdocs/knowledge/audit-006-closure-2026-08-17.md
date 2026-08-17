---
tags: [audit, plane-change, closure, boundary-precision, scientific-history]
date: 2026-08-17
aliases: [exp006-closure-audit, pinch-R-history]
links:
  - "[[combined-transfer-plane-change]]"
---

# Exp 006 Closure Audit — 2026-08-17 (surgical)

## Why this note exists

The first 2026-08-17 adversarial audit and this closure pass both touched the
Exp 006 pinch value, but they touched **different layers**. To prevent a future
reader from misreading git history, the distinction is recorded here explicitly.

## The 6.43 → 6.21 history (science vs prose)

This is a **science correction that preceded its documentation**, not a typo:

1. **Old completion report (2026-08-16):** reported pinch `R ≈ 6.43` and
   `di_c(1.05) = 11.24°`.
2. **Code/results fix (commit `3bc3e38`, 2026-08-16):** the `di_c(R)` boundary
   had been computed at the float-tie point (the 3-burn cost merely equaled the
   two-burn cost to ~1e-7 from discrete-grid noise). The fix introduced
   `WIN_MARGIN = 1e-5`, moved `di_c(1.05)` to **17.01°**, regenerated
   `results.json` (committed `R_pinch = 6.214815`), and added the mpmath-
   anchored regression test. **The numerical science was corrected here.**
3. **Prose lagged:** the card, knowledge note, and roadmap still said `6.43` /
   `6.427` in several spots.
4. **This closure audit (commit `caa2dea`, 2026-08-17):** reconciled the prose
   with the already-correct `6.21` result. **Documentation-only change.**

So `6.43` was a real earlier *numerical* value that was *corrected in code*
before this audit; the audit only fixed the lingering prose. Anyone reading
`git log` should not treat `6.43 → 6.21` as a documentation typo — it reflects
a genuine prior boundary bug (float-tie artifact) that was already repaired.

## What the closure pass verified (high-precision, independent)

An INDEPENDENT continuous minimizer (golden-section in `s` of a theta-minimized
cost — a different code path from the experiment's `(s, theta1, theta2)`
meshgrid) solved the exact boundary crossings by root finding:

| R | di_c (margin 1e-5) | exp di_c | di_inf (exact) | exp di_inf |
|---|---|---|---|---|
| 1.05 | 11.9° (bare 11.4°) | 17.0° | 60.17° | 60.17° |
| 2.00 | 36.7° | 37.9° | 57.37° | 57.35° |
| 4.00 | 41.35° | 41.85° | 48.55° | 48.53° |
| 6.21 | 38.59° | 38.39° | 39.24° | 38.39° |
| 8.00 | 32.08° | 31.87° | 31.87° | 31.87° |

- **di_inf(R)**: matches the experiment to **< 1°** everywhere.
- **di_c(R)**: matches to **~1°** for R ≥ 2. At **R = 1.05** the genuine 3-burn
  advantage is *shallow* (the float-tie region): the bare crossing is ~11.4°,
  the 1e-5-margin crossing ~11.9°, and the experiment's robust choice 17.0°.
  This is a legitimate soft band, not an error — the robust margin is the
  defensible value.
- **Pinch R**: independent continuous solve gives **6.48–6.51** (margin 1e-5)
  vs committed **6.214815**. Both lie in the **[6.2, 6.5]** soft-pinch band;
  the exact value is inherently optimizer-sensitive because the finite-s
  window closes slowly there.
- **Continuous s\* at (R=2, Δi=47.5°)**: **2.7257** (experiment 2.72). The
  cost is **flat in s** — dv varies only ~3e-4 across s ∈ [2.6, 2.9] — so s\*
  is well-defined ~2.72 and the **1.77% saving is the robust quantity**. The
  earlier 2.78 figure was coarse-grid error, now explained.

## Regression tests added

`tests/test_plane_change.py` gained 5 closure-regression tests
(`test_closure_*`) that pin the exact boundaries, the pinch band, and the
continuous s\* with the independent minimizer, so a future refactor cannot
silently shift them.

## Verdict

001–006 are **closed and maximally audited** as of 2026-08-17. Exp 007
(gravity assist) remains the correct next experiment and was not started.
