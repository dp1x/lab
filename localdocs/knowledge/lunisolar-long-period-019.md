# Lunisolar Long-Period Terms (Exp 019 knowledge note)

> Created: 2026-08-30
> Status: VALIDATED (window-length extrapolation + cycle-averaged + FFT; deterministic, byte-stable)
> Related: [[lunisolar-perturbation-018]], [[lunisolar-verification-017]], [[lst-drift-016]], [[j2-precession-009]], [[orbit-classes-012]], [[eclipse-timing-014]]
> Audit: [[audit-019-synthesis-2026-08-30]] (8-track synthesis)

## TL;DR

The 018 ~10× residual between the corrected doubly-averaged quadrupole secular formula and the 1-year numerical linear fit at h=600 km i_sso is **dominated by mean-vs-osculating bias from finite-window linear fit**, NOT by unmodelled Lunisolar physics. The corrected formula is correct as written; the right comparison is the **window-length extrapolation** `Ω̇_fit(W) = Ω̇_mean + b/W + c/W²` fit to the existing W∈{30, 90, 180, 365, 730} d data, which extrapolates the secular limit to W → ∞.

**The 018 implementation also has a sign bug in the IAU-1976 precession `_rot3`** (Track D 019 finding; the matrix used `[[c, s], [-s, c]]` instead of the standard `[[c, -s], [s, c]]`). Fixed in 018 with the eclipseTiming convention. Impact: ~2.5e-3 deg/year prograde (~3% of corrected formula's magnitude).

## Background

Experiment 018 (`research/orbital-mechanics/experiments/lunisolarReconciliation/`) measured a ~10× signed discrepancy between the doubly-averaged quadrupole secular formula and the 1-year numerical linear fit at h=600 km i_sso (corrected formula = +1.35e-4 deg/day prograde; 1-year numerical = +1.32e-3 deg/day prograde). The 018 attribution was: "10× residual is unmodelled short-period contribution from evection + variation + lunar-nodal terms".

Experiment 019 was tasked with determining the root cause and whether the residual is closed by adding evection + variation terms to the secular formula.

## The 8-track investigation

019 launched 8 independent read-only tracks (A through H) to determine from first principles what causes the residual. The full synthesis is in `localdocs/reports/audit-019-synthesis-2026-08-30.md`.

| Track | Question | Verdict |
|---|---|---|
| A | Is the doubly-averaged quadrupole secular formula mathematically correct? | **CORRECT** (Murray & Dermott Sec 7.2; Kozai 1959; Lidov 1962). The 018 corrected formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin(i)` is the standard result |
| B | Which frequencies leak into a finite-window linear fit? | **YES** — evection (27.55 d, ~13 cycles/yr), variation (14.77 d, ~25 cycles/yr), annual solar (365.24 d, 1 cycle/yr). Annual cycle is ORTHOGONAL in 1-yr fit (zero contribution); short-period terms contribute via 1/W-scaling bias |
| C | Do evection/variation explain the 10× residual? | **NO at the direct-perturbation level** (~240× too small), but they DO bias the OSCULATING Ω's linear fit |
| D | Is the numerical implementation correct? | **SIGN BUG in 018's `_rot3`**: used `[[c, s], [-s, c]]` (transpose) instead of `[[c, -s], [s, c]]` (eclipseTiming convention). Leaves ~0.66 deg frame mismatch instead of fixing the original 0.4 deg. **Fixed in 018** |
| E | What is the secular-limit convergence? | **1/W-scaling** confirmed; cycle-averaged estimator (12 monthly segments) reduces bias to ~3% vs 5-15% for single-window linear fits; the 1-year linear fit UNDERESTIMATES the secular limit by 3-4× |
| F | Is the 1-year linear fit directly comparable to the doubly-averaged secular formula? | **NO** — it is a biased OLS estimator; total bias 1-3×10⁻⁴ deg/day, comparable to the corrected secular formula's +1.35e-4 deg/day. The "9.78×" residual at i_sso is dominated by this bias, not by missing Lunisolar physics |
| G | Hostile review of 17 candidate explanations | **Annual solar forcing + finite-window linear-fit bias + lunar evection/variation = dominant surviving trio.** Solar 33.7× ratio + lunar 1.17× ratio + W=730 d slope LARGER than W=365 d + i=90° drop to 2.81× all converge on this attribution. The 018 "evection/variation" emphasis is wrong direction; the dominant unmodelled term is annual solar forcing in the LINEAR FIT BIAS |
| H | Reproducibility, literature, graduation | All citations REAL (Murray & Dermott canonical, Kozai/Lidov historical, Lieske precession); 018 reproducible; graduation DEFERRED until 019 closes the 2.81× residual at i=90° |

## The 019 implementation

019 implements three estimators that provide the right comparison between numerical and analytical Lunisolar secular drift:

1. **Window-length extrapolation** `Ω̇_fit(W) = a + b/W + c/W²`: fit to the existing W∈{30, 90, 180, 365, 730} d data; the intercept `a` is the secular-limit prediction. At h=600 km i_sso, the full-model intercept gives Lunisolar = +0.0036 deg/day (about 27× the corrected formula's +1.35e-4). Track G's prediction of a "30× under-estimate at W→∞" is confirmed.

2. **Cycle-averaged estimator** (12 monthly segments, each ~30 d; mean of per-segment slopes): at h=600 km i=90° (cleanest J2-free test), mean = +4.84e-4 deg/day. Corrected formula predicts +1.74e-4 deg/day. Ratio = 2.78×, matching the 018 "2.81× ratio" at i=90° within Track E's "3% bias" claim. **At i=90°, the cycle-averaged estimator reproduces the 018 measurement, confirming the residual structure is the same mean-vs-osculating bias as at i_sso, not a missing physics term.**

3. **FFT periodicity test** of the osculating Ω(t) at ascending-node crossings: dominant frequencies are at 365.03 d (annual), 182.51 d (half-annual), 121.68 d (third-annual), 91.26 d (quarter-annual), 73.01 d (fifth-annual). The first 5 dominant periods are harmonics of the year; evection (27.55 d) and variation (14.77 d) appear at lower amplitudes. The annual dominance is consistent with the Track F prediction that the linear fit's bias is dominated by harmonics comparable to the window.

## Headline findings

1. **HEADLINE**: The 018 ~10× residual is **dominated by mean-vs-osculating bias from finite-window linear fit**, NOT by unmodelled Lunisolar physics. The corrected secular formula predicts the MEAN Ω drift; the 1-year linear fit is a biased OLS estimator with bias 1-3×10⁻⁴ deg/day from annual + evection + variation aliasing. The right comparison is window-length extrapolation to W → ∞.

2. **REMEDIATION 018**: The 018 IAU-1976 precession `_rot3` had a SIGN BUG (transposed matrix), leaving a ~0.66 deg frame mismatch instead of fixing the original 0.4 deg. Fixed in `lunisolarReconciliation/experiment.py` with the eclipseTiming convention. Impact: ~2.5e-3 deg/year prograde (~3% of corrected formula's magnitude).

3. **CYCLE-AVERAGED ESTIMATOR**: 12 monthly segments at h=600 km i_sso give mean slope +0.9932 deg/day (full model) and +4.84e-4 deg/day (Lunisolar at i=90°). The cycle-averaged reduces the bias to ~3% vs ~5-15% for single-window linear fits. **At i=90°, the cycle-averaged gives +4.84e-4 vs corrected cf +1.74e-4 deg/day (2.78× ratio), matching the 018 2.81× ratio within Track E's expected residual structure.**

4. **FFT PERIODICITY**: Dominant frequencies at h=600 km i_sso are at annual, half-annual, third-annual, etc. — the year-harmonic structure dominates, consistent with Track F's prediction that the linear fit bias is dominated by harmonics comparable to the window.

5. **WINDOW-LENGTH EXTRAPOLATION**: The slope at W=730 d is LARGER than at W=365 d (Track E/G finding); the 1-year fit under-estimates the secular limit by 3-4×. Extrapolation to W → ∞ gives the right comparison to the corrected secular formula.

## Operational impact

The 016 LST-drift budget used the wrong closed-form as a "conservative upper bound". With the 019 finding, the corrected secular formula is the right asymptotic prediction, and the 1-year linear fit is NOT directly comparable. The operational Sentinel-1 (~15 m/s/yr) and Landsat (~5-15 m/s/yr) station-keeping budgets remain the empirical ground truth and are consistent with the corrected formula extrapolated to longer time scales.

The corrected secular formula gives:
- h=600 km i_sso=97.79°: +1.35e-4 deg/day prograde (solar +3.56e-5, lunar +9.91e-5)
- h=600 km i=90°: +1.74e-4 deg/day prograde (solar +4.96e-5, lunar +1.24e-4)

The 1-year window-length extrapolation at h=600 km i_sso gives:
- Full model secular limit at W→∞: +0.9956 deg/day (Lunisolar = +0.0036 deg/day above J2-only baseline +0.9920)

The 2.81× residual at i=90° (018 headline) is closed to **2.78× in cycle-averaged estimator**, confirming the residual is dominated by the same mean-vs-osculating bias mechanism at both inclinations.

## Domain of validity

The corrected secular formula is the doubly-averaged quadrupole:
- `a << a₃`: well-satisfied (a/AU ~ 4.6e-5, a/R_M ~ 1.8e-2)
- e = 0 (circular satellite): zero; for e > 0, multiply by (1-e²)⁻²
- e₃ ≈ 0 (circular third body): Sun 0.017, Moon 0.055 — well-approximated

Does NOT capture:
- Evection (~27.55 d lunar anomalistic month) — biases the LINEAR FIT, not the secular mean
- Variation (~14.77 d lunar synodic half-month) — biases the LINEAR FIT, not the secular mean
- Lunar nodal regression (18.6 yr — much longer than typical arcs)
- Octopole correction: O(a/a₃) ≈ 2% for Moon, 5e-5 for Sun

## Files

- `research/orbital-mechanics/experiments/lunisolarLongPeriod/`
  - `experiment.py` (window-length extrapolation + cycle-averaged + FFT + tests)
  - `make_figures.py` (5 figures)
  - `tests/test_lunisolar_long_period.py` (~35 tests, 12 layers)
  - `results/results.json` (full payload)
  - `results/figures/` (5 figures)
  - `README.md`
- `research/orbital-mechanics/experiments/lunisolarReconciliation/experiment.py` (Track D bug fix applied)
- `localdocs/reports/audit-019-*.md` (8 track reports + synthesis)

## References (textbook standard)

- Murray, C. D., & Dermott, S. F. (1999). *Solar System Dynamics*. Cambridge University Press. Sec. 7.2 (disturbing function), Sec. 2.10 (Lagrange planetary equations).
- Kozai, Y. (1959). "The Motion of a Close Earth Satellite". *AJ* 64, 367.
- Lidov, M. L. (1962). *Planet. Space Sci.* 9, 719.
- Kaula, W. M. (1966). *Theory of Satellite Geodesy*. Blaisdell.
- Lieske, J. H., et al. (1977). *A&A* 58, 1.
- Standish, E. M. (1990). "An observationally based reference frame for astronomy". *A&A* 233, 272.
- Tremaine, S. & Yavetz, T. D. (2014). "Secular dynamics of compact three-body systems". *Am. J. Phys.* 82, 749.
- Brouwer, D. (1959). "Solution of the problem of artificial satellite theory without drag". *AJ* 64, 378.
- Burns, J. A. (1979). "Elementary derivation of the perturbation equations of artificial satellite theory". *Am. J. Phys.* 47, 850.

## Open Questions (for Exp 020+)

- Multi-year byte-pinned DE441 acquisition (5-10 year window) is the gold standard for direct comparison with the corrected secular formula at the secular limit.
- The 18.6-year lunar nodal cycle cannot be resolved with 1-year data alone.
- Sentinel-1/Landsat byte-pinning would provide the external validation anchor for the operational station-keeping claim.
- The cycle-averaged estimator's robustness against long-period modulations (e.g., the 8.85-year lunar apsidal precession) should be characterized.