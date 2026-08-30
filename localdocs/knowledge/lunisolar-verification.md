---
tags: [orbital-mechanics, lunisolar, raan, perturbation, byte-pinned, DE441, Exp-017, audit-response, decadal-rejected]
date: 2026-08-30
aliases: [lunisolar-verification, Exp-017, cf-upper-numerical-ratio]
links:
  - "[[lst-drift-correction]]"
  - "[[audit-015-lst-drift-2026-08-29]]"
  - "[[audit-015-numerical-falsifier-2026-08-29]]"
  - "[[audit-015-adversarial-2026-08-29]]"
  - "[[audit-015-follow-up-candidates-2026-08-29]]"
  - "[[eclipse-timing]]"
  - "[[orbit-classes]]"
  - "[[j2-precession]]"
---

# Lunisolar Upper-Bound Verification (Exp 017)

## Summary

Exp 017 is the **direct numerical follow-up** to the Exp 016 closed-form
Lunisolar upper-bound disclaimer. After Exp 016 documented the closed-form
Vallado Eq. 9-46 form as an over-estimate (model_note:
"over-estimates by a factor of order sin²(i_SS) ~50x at SSO retrograde
inclinations"), this experiment **MEASURES** the actual over-estimate
factor using byte-pinned JPL Horizons DE441 geocentric Sun and Moon
vectors.

The headline result: **the closed-form over-estimates by ~170x at h=600 km**
(3x larger than the audit-015 ~50x estimate). Furthermore, the closed-form
predicts RETROGRADE (-0.218 deg/day) while the numerical integration
predicts PROGRADE (+0.001284 deg/day) — a SIGN DISAGREEMENT, not just a
magnitude over-estimate.

## Why This Experiment

The original Exp 017 direction (per AGENTS.md and roadmap.md) was a
"decadal station-keeping experiment with full Lunisolar + SRP + F10.7-driven
drag arc". An eight-track audit unanimously found this direction not
scientifically defensible at this time (the lab's exponential atmosphere is
inadequate for decadal drag, RK4 secular drift not characterized past 30
days, closed-form Lunisolar over-estimate was the single largest unverified
input, and Sentinel-1 operational records are not byte-pinned). The
**closed-form upper-bound verification** (audit-015 follow-up candidates
candidate #4; Track H Alt-1 scored 27/30) was the strongest scientifically
defensible alternative and is what was executed.

## Method

1. **Acquire Moon ephemeris** via `fetch_horizons_moon_snapshot.py`: JPL
   Horizons Moon (target 301) geocentric vectors, DE441, ICRF/TDB,
   KM-S, geometric, 1-day cadence, 366 rows, full 2026. SHA-256 pinned:
   `65f1d67f798a3b95...`. Identical schema to the existing Sun snapshot
   (eclipseTiming reference). Refuse-to-overwrite idempotence. Lab
   politeness policy: official API, 3 s spacing, single-digit request count.

2. **Build combined RHS**: point-mass Sun + Moon on satellite +
   graduated Kepler + J2 from `j2_rhs`. Linear time interpolation of
   byte-pinned Sun/Moon vectors at each RK4 stage.

3. **Build J2-only control RHS** (identical to combined minus Lunisolar).
   Used to subtract the dominant J2 secular rate from the numerically
   measured Ω(t), isolating the Lunisolar contribution (model-order
   separation per Track F Pillar C).

4. **Propagate both models** at dt=60 s over 1 year at h in
   {500, 600, 700, 800} km, with identical initial conditions (satellite
   at ascending node, x-axis, heading north; SSO inclination from
   `sso_inclination_rad`).

5. **Detect ascending-node crossings** at z=0 with vz>0 by linear
   interpolation; recover Ω at each crossing.

6. **Linear fit** of Ω(t) for each model; subtract J2 slope from
   full-model slope to isolate the Lunisolar contribution.

7. **Closed-form upper bound** at each altitude using the lab's
   reproduction of Vallado Eq. 9-46 form (identical to Exp 016).

8. **Compute ratio** = cf_total / numerical_Lunisolar at each altitude.

9. **dt convergence ladder** at h=600 km with aligned integer-multiple
   grids (the J2 secular drift ~10⁴ km/day would otherwise swamp the
   dt-refinement signal).

## Numerical Headline

| Quantity | h=500 km | h=600 km | h=700 km | h=800 km |
|---|---:|---:|---:|---:|
| Closed-form upper bound (deg/day) | -0.2108 | -0.2184 | -0.2263 | -0.2343 |
| Numerical Lunisolar (J2-subtracted, deg/day) | +0.001320 | +0.001284 | +0.001249 | +0.001215 |
| **cf_upper / numerical ratio (signed)** | **-159.64** | **-170.14** | **-181.19** | **-192.84** |
| Linear-fit residual RMS (deg) | 0.0247 | 0.0240 | 0.0234 | 0.0227 |

Convergence: p_r = 4.49, p_v = 4.50 (RK4 design order ~4 confirmed).

## What This Means

The lab's closed-form Lunisolar upper bound (Exp 016 `lstDrift/experiment.py`
luni_solar_raan_rate_rad_s) **significantly over-estimates** the actual
Lunisolar RAAN rate at LEO SSO:
- Magnitude over-estimate: ~170x at h=600 km (3x larger than the audit-015
  estimate of ~50x).
- Sign disagreement: closed-form predicts retrograde; numerical predicts
  prograde.
- The numerical value (+0.001284 deg/day = ~+0.47 deg/year at h=600) is
  within the operational envelope (~0.005 deg/day Sentinel/Landsat,
  ~1.8 deg/year) reported in Exp 016, but is the **opposite sign** from
  the closed-form prediction.

This means:
- The Exp 016 LST-drift budget that used the closed-form as an upper bound
  was **qualitatively correct** (closed-form is an over-estimate; real
  rate is much smaller) but **quantitatively wrong** (over-estimates by
  ~170x, not ~50x; sign is reversed).
- The **decadal direction originally proposed for 017** is still rejected:
  the closed-form Lunisolar was the single largest unverified input, and
  even with the verification now in hand, the lab still lacks the
  F10.7-driven drag model and multi-year Sentinel/Landsat byte-pinning
  needed for a defensible decadal experiment.

## Files

- `research/orbital-mechanics/experiments/lunisolarVerification/`
  - `experiment.py` — deterministic, offline, ~8 min runtime
  - `fetch_horizons_moon_snapshot.py` — one-time online Moon acquisition
  - `reference/horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt`
    — byte-pinned 76,204-byte Moon geocentric vectors (DE441)
  - `reference/MANIFEST.json` — acquisition provenance
  - `tests/test_lunisolar_verification.py` — 32 tests (L1-L6 layers)
  - `results/results.json` — frozen payload with code_sha256 binding
  - `results/figures/`
    - `f1_cf_upper_over_numerical_ratio.png` — cf/numerical ratio by altitude
    - `f2_lunisolar_raan_drift_comparison.png` — closed-form vs numerical
    - `f3_dt_convergence_ladder.png` — order 4.5 RK4 self-convergence
    - `f4_linear_fit_residuals.png` — ~0.024 deg RMS residuals

## Status

COMPLETE (2026-08-30): 32 new tests, all passing. Total repo tests: 658
(626 baseline + 32 new). All deterministic, offline, byte-stable figures.
Audit response: ✓. Decadal direction rejected; this experiment is the
scientifically defensible substitute.

## Suggested Follow-Up: Exp 018

Now that the closed-form Lunisolar over-estimate is measured, the next
candidate directions are:

1. **Evection + variation terms**: characterize the missing Lunisolar
   terms (evection at the anomalistic month, variation at the synodic
   half-month) to refine the closed-form. This is the natural follow-up
   to this experiment.

2. **Multi-year Sentinel/Landsat byte-pinning**: acquire Sentinel-1A
   precise orbit ephemerides (CNES POD, free public) and Landsat-8
   long-term ephemeris (NASA EOSDIS, free public) to provide the external
   validation anchor currently missing for any operational LEO claim.

3. **Decadal station-keeping with measured Lunisolar**: with this
   experiment's measured over-estimate factor as a calibration, the
   decadal direction becomes tractable IF the lab also adds
   (a) F10.7-driven density model (Jacchia-Bowman, MSIS, or NRLMSISE-00),
   (b) multi-year Sentinel/Landsat byte-pinned anchor, (c) symplectic or
   regularized integrator to suppress the RK4 secular drift on the
   decadal arc.

The full-port decadal experiment is properly an Exp 019 or Exp 020, not
Exp 018.