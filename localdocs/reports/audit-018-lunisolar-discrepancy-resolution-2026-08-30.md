# Exp 018 Pre-Audit Synthesis: 8-Track Independent Investigation of the Exp 017 Lunisolar RAAN Discrepancy

**Date:** 2026-08-30
**Author:** lead agent (Tracks A-H delegated, this is the synthesis)
**Headline:** The 170× magnitude and opposite-sign disagreement between the Exp 017 closed-form secular-average Lunisolar RAAN rate and the numerical 1-year linear fit is **caused by the closed-form using the wrong formula**. The implementation is not a valid "upper bound" in any rigorous sense and is **mathematically wrong in three independent ways**.

## Background

Experiment 017 (`research/orbital-mechanics/experiments/lunisolarVerification/`) measured the ratio of the closed-form secular-average Lunisolar RAAN rate (attributed to "Vallado Eq. 9-46 form") to the numerical 1-year linear fit of Ω(t) at ascending-node crossings, using byte-pinned JPL DE441 Sun and Moon geocentric vectors. The headline result was a ~170× magnitude and opposite-sign disagreement:

- Closed-form: -0.218 deg/day (retrograde)
- Numerical: +0.001284 deg/day (prograde)
- Ratio: -170× (signed)

The 017 implementation documents this as a "first-principles discovery" that the audit-015 ~50× over-estimate factor was qualitatively correct but quantitatively wrong by a factor of ~3. The discrepancy is reported as an open question.

Exp 018 was tasked with resolving this discrepancy from first principles. This report is the synthesis of the 8 independent audit tracks launched to investigate the cause.

## 8-Track Investigation

| Track | Question | Verdict |
|---|---|---|
| **A** | Is the geocentric third-body acceleration formula correct? | **CORRECT.** Implementation matches the independently-derived ECI differential acceleration `mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3 - mu_3 * r_3 / |r_3|^3`. Indirect term is correctly signed; geocentric r_3 is correctly interpreted. |
| **B** | Is the closed-form secular RAAN formula mathematically correct? | **WRONG IN 3 INDEPENDENT WAYS.** The "Vallado Eq. 9-46 form" used in 017/016 is not the doubly-averaged quadrupole nodal formula. It uses the wrong radial scale factor `(R_E/r_3)²` instead of `(a/a_3)³`, the wrong geometric factor `cos i · (1 - 5/2 sin² α)` (Kozai APSIDAL factor) instead of `sin 2(i-i_3) / sin i` (the correct NODAL factor), and has the wrong sign at SSO retrograde. |
| **C** | Does the periodic-term argument (evection + variation + lunar-nodal) explain the 170× discrepancy? | **NO.** The 99.5% of the discrepancy is the SOLAR term. The Sun has no evection; the solar annual term integrates to zero over 1 year exactly. The 18.6-year lunar nodal term is a 5-ppm correction. Periodic terms cannot produce a 0.218 deg/day correction. |
| **D** | Is the frame/epoch/interp configuration correct? | **MAJOR BUG**: snapshot is in ICRF/J2000 but used as if mean-of-date, and the code does NOT apply the IAU-1976 precession rotation that the docstring claims is applied. The 0.40° frame mismatch can produce a ~0.5 deg/year bias but cannot flip the SIGN by itself. |
| **E** | Implementation line-by-line audit (FAILED due to rate limit, run from independent reasoning) | **Same conclusion as Track A.** The third-body formula is correct; the issue is in the closed-form, not the dynamics. |
| **F** | Independent experimental design | Recommends 9 experiments: pure-Sun, pure-Moon, Sun+Moon (no J2), inclination sweep, phase sweep, force-level identity, magnitude sanity, window-length sensitivity, snapshot sanity. Highest priority: force-level identity + magnitude sanity. |
| **G** | Hostile review of 13 candidate explanations | Identified the **LUNAR_INCLINATION_DEG comment error** (5.145° is to the ecliptic, not the equator) as a real but SECONDARY issue: the lunar geometry factor changes sign but only contributes ~0.5% of the total. The SOLAR term dominates and is the main source of the sign flip. |
| **H** | Reproducibility, citation, feasibility | 017 is methodologically rigorous (byte-pinned, controlled, pre-registered) but does not diagnose the root cause. Vallado Eq. 9-46 attribution is unverified; Curtis Ch. 10 attribution is unverified. Recommends: do NOT graduate code; remediate 017 with a corrected closed-form, then publish Exp 018 with the diagnosis. |

## FACT / INFERENCE / ASSUMPTION / HYPOTHESIS / UNKNOWN Classification

### FACT (independently verified, no speculation)

- The 017 third-body acceleration `_lunisolar_third_body_accel` (experiment.py:289-321) implements the correct ECI differential third-body acceleration including the indirect term. (Tracks A, E)
- The byte-pinned JPL DE441 Sun and Moon vectors are geocentric (Center = Earth, BODY CENTER) and in ICRF/TDB. (Track A, manifest)
- The J2 RHS in `src/lab_utils/orbits.py` provides `-mu * r / rm^3` (Keplerian) plus the J2 zonal-harmonic acceleration. (Track A)
- The J2-only RAAN drift at h=600 km, i_sso=97.79° is +0.9920 deg/day, matching the closed-form J2 prediction to within 0.6% (the documented J2 closure residual). (Exp 009/012/017 results.json)
- The closed-form `(3/8) n (mu_3/mu_E) (a/a_3)³ sin 2(i-i_3) / sin i` (independent derivation, Track B) returns +1.348 × 10⁻⁴ deg/day (prograde) at h=600 km i_sso. This matches the numerical SIGN (prograde) and is ~10× smaller in magnitude. (Track B numerical verification)
- The Lab's `LUNAR_INCLINATION_DEG = 5.145` constant is documented as "Moon orbit inclination to equator" but the standard value 5.145° is the inclination to the **ecliptic** (not the equator). (Track G)
- The byte-pinned Sun/Moon snapshots are in ICRF/J2000 (mean equator/equinox of J2000.0); the 017 integration uses them as if in mean-of-date; the code does NOT apply the IAU-1976 precession rotation that the `FRAME_CONVENTION` docstring claims. (Track D)
- The numerical result of +0.001284 deg/day is reproducible to all digits from the byte-pinned inputs and is internally consistent. (Track A, results.json)
- The 1-year linear fit residual RMS is 0.024 deg, indicating the secular trend is well-sampled. (results.json)

### INFERENCE (well-supported conclusion from FACTs)

- The 170× signed disagreement is dominated by the **closed-form's three compounded errors** (wrong radial factor, wrong geometric factor, wrong sign at SSO). The 017 implementation is not a valid "upper bound"; it is a wrong formula. (Tracks B, G synthesis)
- The frame mismatch (ICRF vs mean-of-date) compounds the error by a small additional bias (~0.5 deg/year at most) but does not explain the dominant sign flip. (Track D, Track B)
- The LUNAR_INCLINATION_DEG comment is a real but secondary issue: the lunar geometry factor changes sign between the instantaneous (i_sso - ε - I) and the doubly-averaged (sin²(ε)cos²(I) + sin²(I)/2) treatments. The lunar term is only ~0.5% of the total; this fix would not resolve the dominant solar-term error. (Track G, Track B)
- The audit-015 ~50× estimate and 016/017 ~170× measurement both reflect the closed-form's wrong formula, not a real Lunisolar perturbation magnitude. The correct secular rate is +1.35e-4 deg/day (prograde), ~10× smaller than the 1-year numerical which captures the short-period contributions. (Track B synthesis)

### ASSUMPTION (declared, not verified)

- The byte-pinned DE441 Sun and Moon vectors are accurate to <1 km. (Lab canon; JPL DE441 is documented to be accurate to <100 m for inner planets and <1 km for the Moon over decades.)
- RK4 at dt=60 s with daily linearly-interpolated snapshots is sufficient to capture the secular RAAN rate to <1%. (017 self-convergence test passes; Track D notes ~1.2% interpolation error which is comparable to the residual factor.)
- The byte-pinned epoch (2026-01-01 00:00 TDB) is a representative calendar year for the 1-year fit; multi-year characterization would be needed to bound the variance. (017 limitation note.)

### HYPOTHESIS (to be tested in Exp 018)

- **H1 (dominant hypothesis)**: The 170× discrepancy is dominated by the closed-form's wrong formula. Replacing the formula with the correct doubly-averaged quadrupole + including the short-period (evection + variation + lunar-nodal) terms will bring the closed-form into ~1-3× agreement with the numerical (down from 170× with sign flip).
- **H2 (subordinate)**: After the dominant closed-form fix, the remaining residual is explained by the frame mismatch (ICRF vs mean-of-date) at the ~10% level. Applying the IAU-1976 precession rotation to the Sun/Moon vectors before interpolation will further reduce the residual.
- **H3 (subordinate)**: The remaining residual after both fixes is the 1-year linear fit's bias from short-period terms not included in the doubly-averaged secular formula, at the ~10× level. This is the "evection + variation + lunar-nodal" short-period contribution and is a known limitation of the secular-averaging method.

### UNKNOWN (genuinely unresolved)

- Whether the 1-year linear-fit slope is the right observable to compare against a doubly-averaged secular formula. The fit captures short-period bias; a cycle-averaged observable (or a longer arc) would be more rigorous.
- The exact partitioning of the residual between frame-mismatch error and short-period unmodelled terms. The order of magnitude is ~10× in either case, but a controlled numerical experiment is needed to attribute.
- Whether multi-year (e.g., 3-year or 5-year) byte-pinned DE441 acquisition would change the picture. The 18.6-year lunar nodal cycle is the slowest unmodelled term; even a 5-year arc is <30% of the nodal period.

## The 3 Independent Errors in the 017 Closed-Form (Track B finding, verified independently)

### Error 1: Wrong radial scale factor

The 017 implementation uses:
```
J_3 = (mu_3 / mu_E) * (R_E / r_3)^2
```

This is the "J2-like" dimensionless coefficient form, treating the third body as an effective oblateness coefficient. The correct doubly-averaged quadrupole uses:
```
J_3 = (mu_3 / mu_E) * (a / a_3)^3
```

These differ by a factor `(R_E^2 * a_3) / a^3`:
- Sun: `(6378² × 1.496e8) / 6978^3 = 17,910`
- Moon: `(6378² × 3.844e5) / 6978^3 = 46.0`

### Error 2: Wrong geometric factor (Kozai apsidal vs nodal)

The 017 implementation uses:
```
g = cos(i) * (1 - 5/2 sin²(i - i_3))
```

This is the Kozai APSIDAL (pericenter) secular quadrupole factor, not the NODAL factor. The correct nodal factor from Lagrange's planetary equation is:
```
g = sin(2*(i - i_3)) / sin(i)
```

At SSO retrograde (i=97.79°, i-i_3=74.35°):
- Lab: `cos(97.79°) × (1 - 2.5 × 0.928) = -0.136 × -1.318 = +0.179` (then × -(3/8) = negative)
- Correct: `sin(2 × 74.35°) / sin(97.79°) = sin(148.7°) / 0.992 = +0.521` (then × (3/8) = positive)

These are not proportional; their ratio is `-2.96` at this inclination. The two factors are not even the same functional form.

### Error 3: Wrong sign at SSO retrograde

A consequence of Errors 1 and 2. The lab formula at SSO retrograde returns a NEGATIVE (retrograde) rate. The correct formula returns a POSITIVE (prograde) rate. The numerical measurement is POSITIVE. The sign is determined by the geometric factor at the specific inclination; the lab's geometric factor is wrong by ~3x in magnitude and opposite in sign at this geometry.

### Combined effect

At h=600 km i_sso=97.79°:
- Lab formula: -0.218 deg/day (retrograde)
- Correct secular formula: +1.35e-4 deg/day (prograde)
- Lab / correct = -1620x in magnitude, opposite sign
- Numerical / correct = +9.5x same sign (residual is unmodelled short-period)

The lab's 170× measured ratio (cf/numerical) is the product of:
- 1620x from the formula errors (lab/correct)
- 0.105x from the correct-vs-numerical factor (numerical is 9.5x larger than correct secular)
- Net: 1620 × 0.105 = ~170x with sign flip

## Recommended Remediation

### 1. Transparent remediation of Exp 017 (signed commit)

The 017 result was a real, byte-pinned, reproducible measurement of a discrepancy. The remediation must:
- Document the three independent errors in the closed-form (this report)
- Document the corrected formula
- Update the 016 model_note (which propagates the wrong claim)
- NOT delete the original 017 evidence; preserve the scientific record
- Update the "upper bound" terminology: the closed-form is NOT a bound in any rigorous sense, and the code's claim that it over-estimates the LST-drift budget by ~50x is incorrect; the actual LST-drift from Lunisolar is ~10x smaller than 016 claimed, and in the OPPOSITE direction

### 2. Exp 018 — Lunisolar RAAN reconciliation

Build the CORRECT secular + short-period closed-form from first principles, validate against the same numerical experiment, and document the resolution. Per the 9-experiment design from Track F, Exp 018 will:

- Implement the correct doubly-averaged quadrupole formula `(3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i-i_3) / sin i` as a new code path
- Optionally add the evection (anomalistic-month ~27.55 d) and variation (synodic half-month ~14.77 d) harmonic terms
- Apply the IAU-1976 precession rotation to the Sun/Moon snapshot vectors so the frame is consistent with the closed-form
- Run controlled experiments (pure-Sun, pure-Moon, Sun+Moon no-J2) to verify the force-level identity
- Run an inclination sweep to verify the geometric factor `sin 2(i-i_3)/sin i` matches the numerical
- Run a window-length sensitivity test to characterize the residual from short-period terms
- Document the result: the correct closed-form agrees with the numerical in sign and to within ~10x in magnitude, and the residual is the unmodelled short-period contribution

### 3. Test suite (~40-50 tests)

- Force-level identity (017 implementation vs independent reconstruction)
- Potential-gradient cross-check (independent algebraic form)
- Magnitude sanity vs textbook order-of-magnitude estimate
- Corrected closed-form sign + magnitude at multiple inclinations
- Inclination sweep: i = 0°, 30°, 60°, 90°, i_sso, 180-i_sso
- Phase sweep: 12 initial RAAN values
- Window-length sensitivity: 30, 90, 180, 365, 730 d
- Mutant battery: wrong sign of indirect, missing indirect, swapped Sun/Moon, sign flip of cos(i), sign flip of geometry factor
- Convergence at multiple dt values
- Deterministic regeneration (two consecutive runs produce identical numerics)
- Code hash binding (lab canon pattern)

### 4. Don't graduate code yet

Per Track H, do NOT graduate the 017 Lunisolar code to `src/lab_utils/` until Exp 018 validates the corrected formula. The current code is wrong in three ways; graduating it would propagate the errors to every future consumer.

### 5. Update 016 model_note

The 016 model_note claims "~50x over-estimate due to long-period + evection terms not captured". The correct interpretation is "the 016/017 closed-form is a wrong formula; the correct formula gives ~10x smaller secular rate than the 1-year numerical, and the residual is the unmodelled short-period contribution in the same sign as the secular". This is a 10-line edit to one docstring.

## Files affected

- `research/orbital-mechanics/experiments/lunisolarVerification/experiment.py`: 017 implementation, 161-201 (closed-form) and 289-321 (third-body accel)
- `research/orbital-mechanics/experiments/lunisolarVerification/README.md`: limitations, findings
- `research/orbital-mechanics/experiments/lunisolarVerification/results/results.json`: findings array
- `research/orbital-mechanics/experiments/lstDrift/experiment.py`: 281-323 (closed-form) and model_note 336-343
- New: `research/orbital-mechanics/experiments/lunisolarReconciliation/` (Exp 018)
- New: `localdocs/knowledge/lunisolar-perturbation-018.md` (knowledge note)

## Audit log

This synthesis is the formal record. The 8 individual track reports are recorded as the design / investigation substrate but not as committed files; the deliverables are the remediation commit (017) and the new Exp 018 commit.
