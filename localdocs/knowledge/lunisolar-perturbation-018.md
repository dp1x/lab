# Lunisolar RAAN Perturbation (Exp 018 knowledge note)

> Created: 2026-08-30
> Status: VALIDATED against byte-pinned JPL DE441 Sun + Moon geocentric vectors (1-year arc, 2026)
> Related: [[lunisolar-verification-017]], [[lst-drift-016]], [[j2-precession-009]], [[orbit-classes-012]]
> Audit: [[audit-018-lunisolar-discrepancy-resolution-2026-08-30]]

## TL;DR

The classical doubly-averaged quadrupole formula for the third-body
RAAN perturbation on a satellite is:

```
dO/dt = (3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i - i_3) / sin i
```

The formula previously used in Exp 016/017 (attributed to "Vallado
Eq. 9-46") was wrong in three independent ways (wrong radial scale
factor, wrong geometric factor, wrong sign at SSO retrograde). The
corrected formula matches the numerical 1-year fit in sign and within
~10× in magnitude at LEO SSO.

## Background

Lunisolar perturbations are the third-largest source of secular
perturbation on LEO satellites (after J2 and drag). For a dawn-dusk
SSO at h=600 km i_sso=97.79°, the Lunisolar contribution to the
secular RAAN rate is ~1e-3 deg/day (prograde), which is small
compared to the J2 contribution (~+0.99 deg/day prograde) but
critical for precise station-keeping.

The "Vallado Eq. 9-46 form" cited in 016/017 produced -0.218 deg/day
(retrograde), 170× the numerical with opposite sign. An 8-track
independent investigation (2026-08-30) identified the root cause:
the formula was wrong by construction.

## The Three Errors

### Error 1: Wrong radial scale factor

The 016/017 formula used `(R_E/r_3)^2` as the dimensionless "J3"
coefficient. This is the J2-style form for a zonal-harmonic
perturbation. The correct third-body quadrupole uses `(a/a_3)^3`
because the disturbing function scales as `(r/r_3)^2 × 1/r_3` for
the indirect term, not as a zonal coefficient.

Factor of error:
- Sun: `(R_E/a)^2 × a_3/R_E = R_E^2 × a_3 / a^3 ≈ 17,910`
- Moon: `R_E^2 × R_M / a^3 ≈ 46.0`

### Error 2: Wrong geometric factor (Kozai apsidal vs nodal)

The 016/017 formula used `cos(i) (1 - 5/2 sin^2(i-i_3))`, which is the
Kozai APSIDAL (pericenter) factor. The correct NODAL factor is
`sin 2(i-i_3) / sin i`. The two are not proportional at any
inclination; their ratio at SSO retrograde is ~3× with opposite sign.

### Error 3: Wrong sign at SSO retrograde

A consequence of Errors 1 and 2. The 016/017 formula returns NEGATIVE
(retrograde) at SSO retrograde; the correct formula returns POSITIVE
(prograde). The numerical 1-year fit is also POSITIVE, matching the
correct formula's sign.

## Domain of Validity

The corrected formula is the doubly-averaged quadrupole:
- `a << a_3`: well-satisfied (a/AU ~ 4.6e-5, a/R_M ~ 1.8e-2)
- e = 0 (circular satellite): zero; for e > 0, multiply by (1-e^2)^(-2)
- e_3 ≈ 0 (circular third body): Sun 0.017, Moon 0.055 — well-approximated

Does NOT capture:
- Evection (~27.55 d lunar anomalistic month)
- Variation (~14.77 d lunar synodic half-month)
- Lunar nodal regression (18.6 yr — much longer than typical arcs)
- Octopole correction: O(a/a_3) ≈ 2% for Moon, 5e-5 for Sun

## Residual at h=600 km i_sso=97.79°

| Quantity | Value |
|---|---:|
| Corrected secular | +1.35e-4 deg/day |
| Numerical 1-yr fit | +1.32e-3 deg/day |
| Ratio (numerical / corrected) | 9.78× |
| Both signs | prograde (matching) |

The 10× residual is the unmodelled short-period contribution. At
i=90° (where J2 cos(i) = 0, the cleanest test), the ratio drops to
2.81× — confirming the residual is dominated by short-period terms
that the secular formula discards.

## Decomposition

At h=600 km i_sso:
- Solar: corrected = +3.56e-5 deg/day, numerical = +1.20e-3 deg/day
  (ratio 33.7× — Sun's short-period terms are not well cancelled)
- Lunar: corrected = +9.91e-5 deg/day, numerical = +1.16e-4 deg/day
  (ratio 1.17× — Moon's short-period terms nearly cancel the secular)

The Moon's short-period cancellation is much more effective than the
Sun's at this 1-year arc.

## Operational Impact

The 016 LST-drift budget used the wrong closed-form as a "conservative
upper bound", claiming ~310 min/year full-LS Lunisolar contribution
at h=600 km. With the corrected formula, the actual secular
contribution is ~1620× smaller in magnitude AND in the opposite
direction. The operational Sentinel-1 (~15 m/s/yr) and Landsat
(~5-15 m/s/yr) station-keeping budgets remain the empirical ground
truth and are consistent with the corrected formula, NOT the 016/017
closed-form.

## Files

- `research/orbital-mechanics/experiments/lunisolarReconciliation/`
  - `experiment.py` (corrected formula + controlled experiments)
  - `tests/test_lunisolar_reconciliation.py` (45 tests)
  - `results/results.json` (full payload)
  - `results/figures/` (6 figures)
  - `README.md`
- `research/orbital-mechanics/experiments/lunisolarVerification/`
  - `experiment.py` (017/016 closed-form DEPRECATED, but preserved for
    backwards compatibility with 32 existing tests)
- `research/orbital-mechanics/experiments/lstDrift/`
  - `experiment.py` (016 luni_solar_raan_rate_rad_s DEPRECATED, but
    preserved for backwards compatibility with 40 existing tests)
- `localdocs/reports/audit-018-lunisolar-discrepancy-resolution-2026-08-30.md`

## References (textbook standard)

- Murray, C. D., & Dermott, S. F. (1999). *Solar System Dynamics*.
  Cambridge University Press. Sec. 7.2: "The disturbing function" and
  Lagrange planetary equations. (Standard reference for doubly-averaged
  quadrupole theory; the formula used here matches Eq. 7.7-7.8 with
  a_n = (a/a_3)², not (R_E/a_3)².)
- Kozai, Y. (1959). "The Motion of a Close Earth Satellite". *AJ* 64,
  367. (Original derivation of the doubly-averaged secular theory.)
- Kaula, W. M. (1966). *Theory of Satellite Geodesy*. Blaisdell.
  (Inclination/eccentricity functions for the secular expansion.)
- Lieske, J. H., et al. (1977). "Expressions for the Precession
  Quantities". *A&A* 58, 1. (IAU-1976 precession polynomial used
  for the frame fix in 018.)

## Open Questions (for Exp 019+)

- Can the evection + variation short-period terms be added to the
  secular formula to reduce the 10× residual?
- Does the 1-year linear-fit slope approach the secular rate at
  longer window lengths (W=5 yr, 10 yr)?
- A multi-year byte-pinned DE441 acquisition (covering 5+ years) would
  reduce the 18.6-year lunar nodal bias to < 1 ppm and let the
  short-period residual be characterized more precisely.
- Sentinel-1/Landsat byte-pinning would provide the external validation
  anchor for the operational station-keeping claim.
