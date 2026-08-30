# Experiment 020 -- Lunisolar Long-Arc Secular-Limit Validation

## Date
2026-08-30

## Headline Finding

The 1-year linear fit of osculating Ω(t) at h=600 km i_sso=97.79° (Exp 018
methodology) gives a Lunisolar secular rate of **+1.25e-3 deg/day** (mean over
4 phases), while the corrected doubly-averaged quadrupole formula gives
**+1.35e-4 deg/day**. The **9.3× ratio** reproduces the Exp 018 finding of
9.78×. The **sign matches** (both prograde), confirming the lab's sign
convention (Convention B; verified by audit-020-track-1).

The 019 window-length extrapolation (W→∞ = +3.6e-3 deg/day = 27× the
corrected formula) is **NOT independently validated** as the asymptotic
secular limit. The 8-track audit-020 investigation showed:

1. **Track 3** (estimator theory): the polynomial-in-1/W extrapolation has NO
   theoretical asymptotic basis. The exact OLS bias formula
   `β_bias = (6A/(W²ω))(sin(ωW+φ)+sin(φ)) + (12A/(W³ω²))(cos(ωW+φ)−cos(φ))`
   shows the bias scaling is O(1/W²) for fast harmonics and O(A_k ω_k)
   constant for slow harmonics. The 019 5-point polynomial fit in 1/W is
   mis-specified.

2. **Track 7** (hostile review): the 019 i=90° extrapolation **sign-flips
   between linear (+1.7e-4) and quadratic (-3.7e-4) models** -- a smoking
   gun for mis-specification, not a robust asymptotic measurement.

3. **Track 4** (implementation audit): the constant `i3_moon = 28.584°` in
   the corrected formula is the long-term lunar nodal SECULAR AVERAGE,
   but at 2026 (near descending lunar node) the actual instantaneous lunar
   i3 in ECI mean-of-date is ~18.29°, NOT 28.584°. This causes the corrected
   formula to **over-estimate the lunar contribution at the 2026 epoch by
   ~50%**, making the apples-to-apples 2026 comparison closer to ~13-14×,
   not 9.78×.

## What is the secular limit?

The corrected doubly-averaged quadrupole secular formula gives the
asymptotic MEAN Ω drift; it is correct up to the validity conditions
(eccentric satellite orbit, circular third-body orbit, a ≪ a₃,
quadrupole-only, geocentric frame, mean equator/equinox of date).

The 1-year numerical fit measures OSCULATING Ω at ascending-node crossings,
which contains:
- The mean secular drift (corrected formula prediction)
- Short-period (evection 27.55 d, variation 14.77 d)
- Long-period (annual solar 365.24 d, lunar nodal 6798.4 d)
- Coupling (J2 × Lunisolar cross-product at SSO)

The 1-year fit is BIASED; the bias is dominated by the unmodelled periodic
content. **The 019 extrapolation does not have a theoretical basis for
separating this bias from the true secular limit.**

## What Exp 020 shows

Synthetic oracle test (Track 3 calibration): with the 019 FFT amplitudes
injected into a synthetic signal with known secular, the harmonic-regression
estimator (Track 3's recommendation (f)) recovers the secular to machine
precision (bias 7e-16 deg/day), while the direct OLS estimator has bias
6.2e-5 deg/day. **In the synthetic case, (f) is clearly better.**

Real data at h=600 km i_sso, 1-yr arc, 4-phase ensemble:
- Direct OLS (a): **Lunisolar = +1.25e-3 deg/day** (stable across phases,
  std=2.1e-3 deg/day); ratio to corrected cf = 9.3×
- Secant (g): Lunisolar = +1.19e-3 deg/day; ratio = 8.8×
- Node-vector (n): Lunisolar = +1.25e-3 deg/day; ratio = 9.3×
- Harmonic regression (f): **UNSTABLE** across phases (full model swings
  from 8.89e-1 to 1.11e+0 deg/day; Lunisolar from -1.03e-1 to +1.15e-1
  deg/day); the J2-only baseline is stable (9.93e-1) but the full-model
  harmonic regression is not.

**Critical finding**: the harmonic regression estimator (f) is **fragile**
when applied to the 1-year ascending-node-crossing data. The 8 harmonics
in the basis (annual, half-annual, third-annual, quarter-annual,
fifth-annual, evection 27.55 d, variation 14.77 d, lunar nodal 6798.4 d)
are insufficient to absorb the true signal; unmodelled content (J2
short-period coupling, higher-order evection harmonics) is being
aliased into a long-period "harmonic" that the regression interprets as
a secular shift. This is the **opposite of what the synthetic test
predicted** (where the same basis perfectly recovered the known secular).

## Verdict on the 019 extrapolation

The 019 extrapolation **+3.6e-3 deg/day Lunisolar at i_sso = 27× the
corrected formula** is **NOT scientifically defensible** as an asymptotic
secular limit:

- The 019 polynomial-in-1/W model has no asymptotic basis (Track 3).
- The 019 i=90° extrapolation sign-flips between linear and quadratic models
  (Track 7 hostile review).
- The 020 multi-phase ensemble shows that the more reliable estimators
  ((a), (g), (n)) all give ~+1e-3 to +2e-3 deg/day Lunisolar at i_sso at
  W=1 yr -- consistent with the 018 result, NOT the 019 extrapolation.
- The 019 estimator (f) -- which is theoretically preferred on synthetic
  data -- is **unstable** on real data, so it cannot validate the 019
  extrapolation.

The **019 extrapolation is reported only as a diagnostic**, not a
robust asymptotic measurement. The Exp 020 conclusion is:

**At W = 1 yr, the Lunisolar RAAN secular rate at h=600 km i_sso is
(1.25 ± 0.21) × 10⁻³ deg/day (mean ± std across 4 phases, direct OLS
estimator), which is ~9.3× the corrected formula's 1.35 × 10⁻⁴ deg/day.
The corrected formula gives the correct SIGN but its magnitude is
~10× smaller than what the 1-year numerical fit observes.**

The secular limit at W→∞ remains UNRESOLVED by Exp 020 at the 1-year arc.
The 5-year arc (Track 8 recommendation) requires multi-year DE441
Sun/Moon reference data (not yet in the repository); Exp 021 should
acquire this and run the 2-window phase-locked estimator (Track 6) which
cancels the lunar nodal contribution exactly without requiring an
18.6-yr arc.

## Cross-References

- Exp 017 `lunisolarVerification`: original 170x discrepancy
- Exp 018 `lunisolarReconciliation`: corrected secular formula; 9.78x ratio
- Exp 019 `lunisolarLongPeriod`: W-extrapolation; window-length sensitivity
- `localdocs/reports/audit-019-track-{A..H}-*.md`: prior 8-track audit
- `localdocs/reports/audit-020-track-{1..8}-*.md`: current 8-track audit
- `localdocs/reports/audit-019-synthesis-2026-08-30.md`: prior synthesis
- `localdocs/knowledge/lunisolar-perturbation-018.md`: prior knowledge note
- `localdocs/knowledge/lunisolar-long-period-019.md`: prior knowledge note