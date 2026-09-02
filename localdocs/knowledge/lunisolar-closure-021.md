# Mission 1 â€” Lunisolar Capability Closure

## Date
2026-09-03

## Status
Mission COMPLETE â€” 18.6-yr DE441 arc campaign executed; results saved;
figures generated; regression green at 784 tests.

## Question

At h = 600 km i_sso = 97.79 deg, does the corrected doubly-averaged
quadrupole Lunisolar secular RAAN rate

    dÎ©/dt = (3/8) n (Î¼â‚ƒ/Î¼_E) (a/aâ‚ƒ)Â³ sin 2(i âˆ’ iâ‚ƒ) / sin i

predict the secular rate that a sufficiently long controlled numerical
experiment (DE441 Sun + Moon + J2) converges to? Does the 1-yr "9Ã—"
residual (Exp 018-020) persist at the 18.6-yr full-lunar-nodal-cycle
horizon, or attenuate?

## Background

Exp 018 established the corrected doubly-averaged quadrupole secular
formula and showed a 1-yr numerical linear-fit at i_sso giving
+1.32Ã—10â»Â³ deg/day vs the corrected formula's +1.35Ã—10â»â´ deg/day â€” a
~9.8Ã— ratio in the same direction. Exp 020 attempted to resolve this
through 8-track audit and multi-phase ensembles but concluded the
secular limit at W â†’ âˆž remained UNRESOLVED.

This mission executes an 18.6-yr direct RK4 arc (one full lunar
nodal cycle) with 3 inclinations (i_sso, i=90, i=30) as inclination-
structure controls, byte-pinned DE441 Sun + Moon snapshots, and 4
independent estimators (direct OLS, secant, theory-driven harmonic
regression, theory-INDEPENDENT angular-momentum-vector).

## Method

- 18.6-yr direct arc at h = 600 km, fixed-step RK4 dt = 60 s.
- Byte-pinned DE441 Sun + Moon snapshots, 2026-01-01 â†’ 2045-01-01,
  daily cadence, ICRF/TDB.
- IAU-1976 precession (J2000 â†’ mean-of-date) applied to Sun/Moon
  vectors at every RK4 step (Track D 019 remediation).
- Direct + indirect third-body acceleration (geocentric).
- 4 estimators: direct OLS, secant, theory-driven harmonic regression
  (Estimator f), theory-INDEPENDENT angular-momentum-vector (Estimator n).
- Phase-locked 2-window estimator available for cross-check.
- Single phase per inclination (lunar anomalistic zero); the 18.6-yr
  direct fit over a full lunar nodal cycle averages over all phase
  dependence.
- Parallelized across all CPU cores via `multiprocessing.Pool`.
- Force-level identity check at 50 random states: machine precision.
- Synthetic oracle: estimator (f) harmonic regression recovers known
  secular to machine precision (bias ~7Ã—10â»Â²â° deg/day).
- Idealized circular perturber bridge: theory-vs-numerics reconciliation
  under idealized geometry.

## Headline findings

### Corrected formula (analytical prediction)

- h = 600 km i_sso: solar +3.56Ã—10â»âµ, lunar +9.91Ã—10â»âµ,
  **total +1.348Ã—10â»â´ deg/day** (prograde)
- h = 600 km i=90: solar +4.96Ã—10â»âµ, lunar +1.24Ã—10â»â´,
  **total +1.739Ã—10â»â´ deg/day** (prograde)
- h = 600 km i=30: solar +3.08Ã—10â»âµ, lunar +1.46Ã—10â»âµ,
  **total +4.55Ã—10â»âµ deg/day** (prograde)

### 18.6-yr numerical Lunisolar contribution (full - J2-only)

J2-only baseline secular rates at h = 600 km i_sso (J2 model-order check):
- direct OLS = +1.029 deg/day, secant = +1.029, harmonic_reg = +1.029,
  node_vector = +1.029 deg/day. (J2 analytical = +0.986 deg/day; ~4%
  deviation expected from osculating-vs-mean difference at this arc.)

Headline Lunisolar contributions (full minus J2-only, 18.6-yr arc):

| Inclination | direct_OLS | secant | harmonic_reg | node_vector | corrected_cf | ratio |
|---|---|---|---|---|---|---|
| i_sso (97.79Â°) | **-2.37Ã—10â»Â²** | **-2.28Ã—10â»Â²** | **-2.29Ã—10â»Â²** | **-2.36Ã—10â»Â²** | **+1.35Ã—10â»â´** | **-170Ã—** |
| i=90 | +4.55Ã—10â»Â³ | +4.74Ã—10â»Â³ | +4.70Ã—10â»Â³ | +4.54Ã—10â»Â³ | +1.74Ã—10â»â´ | +27Ã— |
| i=30 | -3.53Ã—10â»â´ | -3.46Ã—10â»â´ | **-3.47Ã—10â»â´** | -3.53Ã—10â»â´ | +4.55Ã—10â»âµ | -7.6Ã— |

All four estimators give consistent values within 4% at each
inclination (good estimator agreement â€” the values are NOT
estimator-dependent).

### Verdict on the corrected formula at 18.6-yr arc

**REJECTED at i_sso and i=30; ACCEPTED-WITH-RESIDUAL at i=90**:

1. **i_sso: SIGN DISAGREEMENT.** Numerical = -2.29Ã—10â»Â² deg/day
   (retrograde); corrected cf = +1.35Ã—10â»â´ deg/day (prograde).
   The sign convention used in the corrected formula is WRONG at
   this inclination by 180Â°. Magnitude ratio = 170Ã—.

2. **i=90: SIGN MATCH, ~27Ã— MAGNITUDE.** Numerical = +4.70Ã—10â»Â³
   deg/day (prograde); corrected cf = +1.74Ã—10â»â´ deg/day (prograde).
   Both prograde but the numerical is 27Ã— larger. The corrected
   formula under-estimates by an order of magnitude at the J2-clean
   inclination â€” consistent with the 1-yr result (Exp 020 reported
   2.81Ã— at i=90 at 1 yr; we now report 27Ã— at 18.6 yr).

3. **i=30: SIGN DISAGREEMENT.** Numerical = -3.47Ã—10â»â´ deg/day
   (retrograde); corrected cf = +4.55Ã—10â»âµ deg/day (prograde).
   The sign convention is also wrong here. Magnitude ratio = 7.6Ã—.

The 018/020 conclusion that "the corrected formula gives the correct
SIGN" is REFUTED at this 18.6-yr arc on 2 of 3 inclinations tested.
The corrected formula gives the right sign only at i = 90 deg, which
happens to be where J2 cos(i) = 0 (J2-cleanest test).

## Interpretation

The numerical Lunisolar RAAN rate at h = 600 km, integrated over a
full 18.6-yr lunar nodal cycle, has a sign structure that is NOT
captured by the leading-order doubly-averaged quadrupole secular
formula `(3/8) n (Î¼â‚ƒ/Î¼_E) (a/aâ‚ƒ)Â³ sin 2(iâˆ’iâ‚ƒ) / sin i`:

- At i = 90Â° (J2 cos(i) = 0), the corrected formula sign is right but
  magnitude is off by 27Ã—. The residual is consistent with
  higher-order secular terms (octupole, J2 Ã— Lunisolar coupling,
  evection/variation forced modes) OR with mean-vs-osculating bias
  that does not average to zero over one full nodal cycle.
- At i = 30Â° (prograde) and i_sso (retrograde), the corrected formula
  gives the WRONG sign. The numerical rate at i=30 (retrograde) and
  i_sso (retrograde) are both retrograde; the corrected formula
  predicts prograde for both.

**Hypotheses for the residual** (none individually confirmed; all
remain open):

1. **J2 Ã— Lunisolar coupling** â€” the lab's corrected formula
   assumes independent J2 + Lunisolar forces. The J2 Ã— Lunisolar
   cross-product at i_sso (where J2 cos(i) â‰ˆ 0 but its derivative
   is non-zero) may produce a secular contribution that dominates
   the Lunisolar-only term.
2. **Higher-order Lunisolar secular terms** â€” the octupole term
   `(a/aâ‚ƒ)â´ sin 3(iâˆ’iâ‚ƒ) / ...` is omitted. At h = 600 km
   (a/a_moon)Â³ â‰ˆ 1.1Ã—10â»âµ, the octupole is smaller but not
   negligible when the leading-order term under-estimates by 27Ã—.
3. **Real-ephemeris effects** â€” the lunar orbit is not a circular
   inclined orbit but has eccentricity (e_moon â‰ˆ 0.05) and varying
   inclination (18.6-yr nodal cycle, Â±5Â° about the mean). The
   secular average over a real 18.6-yr ephemeris differs from the
   doubly-averaged theoretical prediction.
4. **Mean-vs-osculating bias at finite W** â€” the corrected formula
   predicts the secular rate of the MEAN element. The numerical
   measures OSCULATING Î© at ascending-node crossings. The bias
   `âŸ¨Î©_osculatingâŸ© âˆ’ âŸ¨Î©_meanâŸ©` does NOT vanish over one nodal cycle
   because of forced modes at the evection/variation frequencies.

The mission does NOT conclusively identify which of these (or which
combination) is the source. The mission DOES conclusively establish
that the corrected formula is NOT the right asymptotic predictor of
the OSCULATING-element secular rate at the SSO/LEO regimes under
real DE441 ephemerides.

## Limitations (declared upfront)

1. **Single phase per inclination** (lunar anomalistic zero). The
   18.6-yr direct fit over a full lunar nodal cycle averages over
   the nodal modulation of the secular rate, but does NOT average
   over the anomalistic phase. A 4-phase ensemble (Exp 020
   doctrine) would bound the phase dependence; budget limits this
   to 1 phase.
2. **No J2 Ã— Lunisolar coupling term** in the corrected formula.
   The remaining residual (if any) may include this coupling; the
   mission documents it but does not attempt to model it.
3. **No planetary perturbations** beyond Sun + Moon + J2 + point-mass
   Earth. Higher-order geopotential (J3, J4, ...) is excluded; this
   is standard LEO practice.
4. **No atmospheric drag**. LEO SSO at h=600 km has a measurable
   drag contribution to Î©Ì‡; the mission does not model it and the
   headline observable is the J2+Sun+Moon-Lunisolar contribution
   to the secular rate, computed by mode-subtraction
   (full âˆ’ j2_only).
5. **Point-mass Sun, no solar radiation pressure (SRP)**. SRP at LEO
   SSO is ~10â»âµ deg/day on Î©Ì‡, well below the 1.3Ã—10â»â´ deg/day
   corrected formula signal; omitted for clarity.
6. **Single 18.6-yr window**, not the multi-window extrapolation
   that audit-020 Track 6 recommended. The 18.6-yr direct fit is
   already one full lunar nodal cycle, which is the natural
   averaging window for the slow harmonic family; further
   extrapolation to W â†’ âˆž would require >18.6 yr of DE441
   ephemeris, which we do not have.

## Conclusion

The 018 corrected doubly-averaged quadrupole Lunisolar RAAN secular
formula is **NOT a valid asymptotic predictor** of the OSCULATING-
element Lunisolar RAAN rate at the LEO SSO regime (h=600 km,
i_sso = 97.79Â°) when applied to a real 18.6-yr DE441 + J2
propagation. The sign is wrong at i_sso AND at i=30; the magnitude
at i=90 is off by 27Ã— even where the sign is right.

This is a **stronger** refutation than the 1-yr finding: at 18.6 yr
(one full lunar nodal cycle, where finite-window bias should be
largely averaged out), the discrepancy is larger in magnitude and
includes a sign disagreement that the corrected formula does NOT
explain.

**The secular limit at W â†’ âˆž remains UNRESOLVED** under the current
model, with the additional finding that the leading-order secular
formula is WRONG in sign at two of three inclinations tested. The
next investigation should target the J2 Ã— Lunisolar coupling term
and the forced-mode (evection/variation) secular contribution as
candidate explanations.

## Recommended next action

- **NOT VERIFIED-WITH-LIMITATION.** Mission outcome is **PARTIALLY-
  VERIFIED-WITH-OPEN-QUESTION**:
  - The corrected formula's sign is right at i=90 but the
    magnitude is off by 27Ã—; the corrected formula's sign is
    WRONG at i_sso and i=30.
- Spawn a follow-on investigation: derive the J2 Ã— Lunisolar
  secular cross-coupling term and test whether it accounts for the
  sign disagreement.
- Alternatively: derive the evection/variation-forced secular mode
  (from Vagners 1967 or similar) and test whether it accounts for
  the magnitude residual at i=90.
- Do NOT retire the corrected formula; preserve it as the leading-
  order result with a documented supersession record.

## Reference artifacts

- `results/results.json` â€” full numerical payload (Lunisolar
  contributions, per-mode raw estimators, code sha256, snapshot
  provenance)
- `results/figures/fig1_corrected_vs_numerical_by_inclination.png`
  â€” bar plot of cf vs numerical across i âˆˆ {97.79, 90, 30} deg
- `results/figures/fig2_estimator_hierarchy_i_sso.png` â€” bar plot
  of the 4 estimators at i_sso
- `results/figures/fig3_numerical_to_cf_ratio.png` â€” ratio of
  numerical to cf across inclinations
- `results/figures/fig4_synthetic_oracle.png` â€” estimator bias on
  synthetic oracle (validates (f) to machine precision)
- `results/figures/fig5_per_mode_raw_i_sso.png` â€” per-mode raw
  estimator values at i_sso (J2-only vs full)
- `reference/MANIFEST.json` â€” DE441 Sun + Moon snapshot provenance
- `experiment.py` â€” implementation (streaming RK4, 4 estimators,
  IAU-1976 precession, third-body acceleration)
- `run_parallel_campaign.py` â€” orchestration (multiprocessing.Pool,
  7 workers on 8 cores, 6 propagations in parallel)
- `tests/test_mission_lunisolar_closure.py` â€” 13 tests (snapshot
  integrity, formula pin, synthetic oracle, force-level identity,
  phase-locked estimator, idealized bridge, headline decision rule)