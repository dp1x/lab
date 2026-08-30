# Experiment 020 -- Lunisolar Long-Arc Secular-Limit Validation

## Research Question

At h = 600 km i_sso = 97.7876 deg, does the doubly-averaged quadrupole
lunisolar secular RAAN rate from the corrected Exp 018 formula
`(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i`
predict the secular limit that a sufficiently long controlled numerical
experiment converges to? Or does the formula systematically
under-estimate the actual secular drift (as the 019 extrapolation
`Ω̇_fit(W → ∞) = +3.6×10⁻³ deg/day` claimed, 27× the corrected formula)?

## Background

Exp 018 (`lunisolarReconciliation`) established the corrected doubly-averaged
quadrupole secular formula and showed that at h = 600 km i_sso it gives
+1.35×10⁻⁴ deg/day, while a 1-year numerical linear-fit of osculating Ω gives
+1.32×10⁻³ deg/day — a 9.78× discrepancy. Exp 019
(`lunisolarLongPeriod`) attributed the discrepancy to mean-vs-osculating bias
from the 1-year linear fit and used a window-length extrapolation
`Ω̇_fit(W) = a + b/W + c/W²` to estimate the W → ∞ secular limit at
+3.6×10⁻³ deg/day (27× the corrected formula).

Audit-020 (8-track independent investigation, `localdocs/reports/audit-020-track-{1..8}`)
found:
- **Track 1** (Track A independent re-derivation): the lab's sign convention
  (Convention B) is correct and consistent with the numerical at i_sso; the
  prior Track A audit's 35× discrepancy came from using Convention A.
- **Track 2** (periodic terms + bias): the 019 1/W extrapolation has NO
  theoretical asymptotic basis; the actual OLS bias scaling is O(1/W²) for
  fast harmonics and O(A_k ω_k) constant for slow harmonics.
- **Track 3** (estimator theory): derived the exact OLS bias formula
  `β_bias = (6A/(W²ω))(sin(ωW+φ)+sin(φ)) + (12A/(W³ω²))(cos(ωW+φ)−cos(φ))`;
  identified theory-driven harmonic regression as the best estimator (bias
  ~3×10⁻⁶ deg/day at W=1 yr vs ~10⁻⁴ deg/day for direct OLS).
- **Track 4** (implementation audit): no model-order errors in 017/018/019;
  DISCOVERY: at 2026 (near descending lunar node), the actual lunar i₃
  is ~18.29°, not the secular mean of 28.584°, making the actual 9.78× residual
  closer to **13-14× at 2026 epoch**.
- **Track 5** (independent estimator): recommended mean-element integration
  via Lagrange planetary equations (theory-aligned headline) plus angular-
  momentum-vector geometry (theory-INDEPENDENT cross-check).
- **Track 6** (long-arc design): recommended a **5-year baseline arc** as the
  minimum scientifically justified duration; a **2-window phase-locked
  estimator** at 9.3-yr separation cancels the lunar nodal contribution
  exactly without requiring an 18.6-yr arc.
- **Track 7** (hostile review of 019): the 019 i=90° extrapolation **sign-flips
  between linear and quadratic models** (linear: +1.7×10⁻⁴, quadratic:
  −3.7×10⁻⁴); the i_sso extrapolated +0.0036 deg/day is ~7× the i=90°
  cycle-averaged 0.000484 — an asymmetry the corrected formula cannot
  explain. Verdict: 019 extrapolation is **an artifact of model choice**,
  not a robust asymptotic measurement.
- **Track 8** (compute feasibility): 5-yr arc ≈ 3.5 hr single-core;
  10-yr ≈ 7 hr; 18.6-yr ≈ 13 hr. Storage 2.5-9 MB. Acquisition ~120 s HTTP.

## References

- Murray & Dermott (1999), *Solar System Dynamics*, Cambridge UP, Ch. 6 (lunar
  theory) and Ch. 7 (third-body disturbing function).
- Kaula (1962), "Development of the lunar and solar disturbing functions",
  *AJ* 67, 300.
- Kozai (1959), "The motion of a close earth satellite", *AJ* 64, 367.
- Lidov (1962), "The evolution of orbits of artificial satellites of
  planets", *PSS* 9, 719.
- Brouwer & Clemence (1961), *Methods of Celestial Mechanics*, Academic Press.
- Standish (1990), "An observationally based reference frame for astronomy",
  *A&A* 233, 272 (window-length extrapolation precedent).
- Exp 014 `eclipseTiming` Sun snapshot acquisition doctrine.
- Exp 017 `lunisolarVerification` Moon snapshot acquisition doctrine.
- Exp 018 `lunisolarReconciliation` corrected secular formula.
- Exp 019 `lunisolarLongPeriod` window-length extrapolation.

## Methodology

### Headline observable: theory-driven harmonic regression

We fit
`Ω(t) = a + b·t + Σ_k [c_k cos(2π t/T_k) + s_k sin(2π t/T_k)]`
simultaneously via OLS. `b` is the secular rate; its bias is **exactly zero
for each harmonic in the basis** (projection theorem); residual bias comes
only from unmodelled content.

The basis `T_k` includes the named physical drivers and the integer-cycle
annual harmonics detected by Exp 019's FFT:
- 365.2422 d (annual solar forcing)
- 182.6211 d (half-annual)
- 121.7474 d (third-annual)
- 91.3106 d (quarter-annual)
- 73.0484 d (fifth-annual)
- 27.5546 d (evection / lunar anomalistic)
- 14.7653 d (variation / lunar synodic half-month)
- 6798.4 d (lunar nodal)

### Cross-check observable: angular-momentum-vector secular rate

The kinematic node vector `n = z × h = z × (r × v)` is computed at every
RK4 step. `arctan2(n_y, n_x)` gives `Ω` without invoking ascending-node
detection. The OLS slope of this time series is the secular rate, independent
of the Lagrange planetary equation theory.

### Reference observables (cross-comparison)

- Direct OLS linear-fit slope (Exp 018 baseline)
- Cycle-averaged estimator (Exp 019 pattern)
- Window-length extrapolation in `1/W` (Exp 019 pattern)
- Exp 018 corrected secular formula (`+1.3475×10⁻⁴ deg/day` at i_sso)

### Numerical design

- **Arc length**: 5-year baseline (Track 8 recommendation), with fallback to
  1-year if multi-year snapshots unavailable.
- **Phase ensemble**: 4 phases spaced ~6.89 d apart (quarters of the lunar
  anomalistic month) to characterize phase dependence.
- **Force modes**: `sun_moon_j2`, `sun_moon`, `moon_only`, `sun_only`,
  `j2_only` (5 modes × 4 phases × 5-yr arc = 20 propagations).
- **Integrator**: fixed-step RK4 at dt = 60 s (Exp 018/019 verified design).
- **Reference data**: byte-pinned DE441 Sun + Moon at daily cadence.

### Validation

- Force-level identity check at 50 random states (machine precision).
- Convergence ladder on a 30-day subset (dt = 30, 60, 120 s).
- Precession identity at T = 0; non-identity at T = 0.26 centuries.
- Synthetic estimator test: known secular + 8 harmonics at Exp 019 FFT
  amplitudes; verify which estimator (a) vs (f) recovers the secular.

## Results

See `results/results.json` and `results/figures/` for the full numerical
output. Headline finding below.

## Limitations

- 5-year arc does NOT resolve the 18.6-year lunar nodal cycle; only partial
  nodal modulation is captured. Track 6's 2-window phase-locked estimator
  would cancel the nodal contribution exactly but requires a longer
  ephemeris.
- Multi-year (5-yr) byte-pinned DE441 Sun/Moon snapshots are not yet in the
  repository; Exp 020 currently uses the 1-year fallback.
- The harmonic-regression basis assumes the dominant periodic content is at
  integer-cycle annual harmonics + evection + variation + lunar nodal.
  Unmodelled harmonics with significant amplitude at intermediate periods
  can bias the secular estimator.

## Next Question (021)

The 2-window phase-locked estimator at 9.3-yr separation (Track 6) would
cancel the lunar nodal contribution exactly. Exp 021 should acquire a 10+
year DE441 Sun/Moon reference and run the phase-locked estimator as the
gold-standard confirmation of the 020 finding.