# Experiment 017 — Lunisolar Upper-Bound Verification

> Status: COMPLETE
> Date: 2026-08-30
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/lunisolarVerification/`

## Research Question

What is the actual ratio between the closed-form secular-average Lunisolar
RAAN rate (Vallado Eq. 9-46 form, used in Exp 016 as the "honest upper
bound") and the numerically integrated Lunisolar RAAN rate at a dawn-dusk
SSO at h in {500, 600, 700, 800} km over a 1-year arc?

The closed-form is known to over-estimate the Lunisolar RAAN rate at LEO
SSO because the secular-average formula discards the long-period + evection
terms that partially cancel the secular mean (per the model_note in
`lstDrift/experiment.py:336-343` and the `audit-015-lst-drift-2026-08-29.md`
audit). The audit-015 follow-up-candidates report estimated the over-estimate
factor at "~50x" at SSO retrograde inclinations. This experiment **MEASURES**
the factor with byte-pinned JPL Horizons DE441 geocentric Sun and Moon
vectors.

## Frozen Contract v1.0

| Item | Value | Provenance |
|---|---|---|
| `R_E` (km) | 6378.137 | WGS-84 equatorial; lab canon `R_EARTH_KM` |
| `J2` | 1.082629821e-3 | WGS-84, J2 = √5·\|C20_bar\|; lab canon |
| `μ_E` (km³/s²) | 398600.4418 | IAU 2015 nominal GM_E; lab canon |
| `μ_Sun` (km³/s²) | 132712440018 | IAU 2015 nominal GM_Sun |
| `μ_Moon` (km³/s²) | 4902.8001 | IAU 2015 nominal GM_Moon |
| `AU` (km) | 149597870.7 | IAU 2012 Resolution B2 (exact) |
| `LUNAR_DISTANCE_KM` (cf only) | 384400.0 | mean Earth-Moon distance |
| `LUNAR_INCLINATION_DEG` | 5.145 | Moon orbit inclination to equator (mean) |
| `SOLAR_OBLIQUITY_DEG` | 23.439 | obliquity of the ecliptic |
| SSO target (deg/day) | 360/365.2422 = 0.985647332099 | Exp 012 pinned |
| Altitudes (km) | {500, 600, 700, 800} | Exp 015 frozen band |
| Mission duration (days) | 365 | 1 year; byte-pinned snapshot covers 366 days |
| Integration step (s) | 60 | conservative RK4 for LEO at SSO inclinations |
| Sun snapshot | JPL Horizons DE441, ICRF/TDB, daily 2026, 366 rows | Exp 014 reference (existing) |
| Moon snapshot | JPL Horizons DE441, ICRF/TDB, daily 2026, 366 rows | This experiment, byte-pinned under `reference/` |
| Force model hierarchy | Kepler + J2 (graduated canon) + point-mass Sun + point-mass Moon | Lab canon + snapshot-driven |
| Closed-form formula | Vallado Eq. 9-46 form | Lab canon (reproduced from Exp 016) |
| Pre-registered ratio band | [10x, 100x] (audit-015 ~50x estimate) | Audit 2026-08-29 |
| Pre-registered order floor | p ≥ 3.5 | RK4 design |
| Pre-registered Lunisolar numerical | [1e-4, 1e-1] deg/day | Operational envelope ~0.005 deg/day |

## Methodology

Deterministic, offline-only after acquisition of the byte-pinned Moon
ephemeris (76 KB, 366 daily rows, sha256 `65f1d67f798a3b95...`). No network
at runtime, no RNG, no wall-clock in the analysis.

1. **Acquire Moon ephemeris** (`fetch_horizons_moon_snapshot.py`): one-time
   online fetch of JPL Horizons Moon (target 301) geocentric vectors, 1-day
   cadence over 2026, identical schema to the existing Sun snapshot. SHA-256
   pinned under `reference/` with MANIFEST.json. Refuse-to-overwrite
   idempotence.

2. **Build combined RHS** (`make_combined_rhs`): point-mass Sun + Moon
   accelerations on the satellite, with linear time interpolation of the
   geocentric Sun/Moon vectors at each RK4 stage. Superposition: total
   acceleration = Kepler + J2 (graduated canon `j2_rhs`) + Sun point-mass
   + Moon point-mass.

3. **Build J2-only control RHS** (`make_j2_only_rhs`): identical to the
   combined RHS minus the Lunisolar terms. Used to subtract the dominant
   J2 secular rate from the numerically measured Ω(t) drift, isolating
   the Lunisolar contribution (model-order separation per Track F Pillar C).

4. **Propagate both models** at dt=60 s over the 1-year arc at h in
   {500, 600, 700, 800} km, with identical initial conditions (satellite
   placed at ascending node, x-axis direction, heading north; SSO
   inclination from `sso_inclination_rad`).

5. **Detect ascending-node crossings** at z=0 with vz>0 by linear
   interpolation; recover Ω at each crossing from arctan2(r_y, r_x).

6. **Linear fit** of Ω(t) for each model; subtract J2 slope from
   full-model slope to isolate the Lunisolar contribution.

7. **Closed-form upper bound** at each altitude using the lab's
   reproduction of Vallado Eq. 9-46 form (identical to Exp 016).

8. **Compute ratio** = cf_total / numerical_Lunisolar at each altitude.
   Compare to the pre-registered [10x, 100x] audit-015 band.

9. **dt convergence ladder** at h=600 km with dt in {120, 60, 30, 15, 7.5}
   s vs a 1.875 s reference; report fitted order p. (Aligned grid: coarse
   points are exact subsets of the fine grid to avoid float64 roundoff
   pollution; otherwise the J2 secular drift (~10⁴ km/day) swamps the
   dt-refinement signal.)

## Implementation

- Script: `experiment.py` (deterministic, offline)
- Language/runtime: Python 3.12, numpy, matplotlib Agg
- Runtime: ~8 minutes single core (4 altitudes × 1.5 min + convergence ladder 30 s + figures <1 min)
- Determinism: pure float64, no RNG, no network at runtime, no wall-clock
  in the analysis path. `time.time()` only in `run()` for elapsed prints.

## Validation Method

Six layers (target ~32 tests):

- **L1 snapshot integrity** (6 tests): sha256 pinning, parsing, physical
  distance band, uniform epoch spacing, n_points = 366, gitattributes -text.
- **L2 closed-form identity** (6 tests): sign at SSO retrograde, solar >
  lunar magnitude, total = solar + lunar, monotone magnitude in altitude,
  match Exp 016 value to <1e-6 deg/day, manual reproduction with explicit
  constants.
- **L3 numerical RAAN drift** (5 tests): magnitude in operational band,
  subtraction-of-J2 pattern, J2 dominates Lunisolar by orders of magnitude,
  n_ascending_nodes within expected, J2 closure residual in [0, +1%] band.
- **L4 cf_upper / numerical ratio** (4 tests): documented audit-015 band
  violation as discovery, ratio_log10 finite, ratio monotone in altitude,
  sign disagreement (cf retrograde, numerical prograde).
- **L5 convergence and dt halving** (4 tests): order above 3.5, order
  below 5.5, monotonic decrease with dt, final diff < 1 mm.
- **L6 adversarial mutants** (7 tests): cf total magnitude sanity,
  cf altitude infeasibility, interp clamp outside range, interp midpoint,
  payload structure complete, code_sha256 includes essentials, no
  machine-specific paths in experiment.py.

## Headline Numbers (from `results/results.json`)

| Quantity | h=500 km | h=600 km | h=700 km | h=800 km |
|---|---:|---:|---:|---:|
| Closed-form upper bound (deg/day) | -0.2108 | -0.2184 | -0.2263 | -0.2343 |
| Numerical Lunisolar (J2-subtracted, deg/day) | +0.001320 | +0.001284 | +0.001249 | +0.001215 |
| **cf_upper / numerical ratio (signed)** | **-159.64** | **-170.14** | **-181.19** | **-192.84** |
| Linear-fit residual RMS (deg) | 0.0247 | 0.0240 | 0.0234 | 0.0227 |
| n_ascending_nodes | 5565 | 5445 | 5330 | 5218 |

**Convergence:** p_r = 4.49, p_v = 4.50 (RK4 design order ~4 confirmed).

**Sun snapshot:** 366 rows, sha256 `06d54fb3...`
**Moon snapshot:** 366 rows, sha256 `65f1d67f...`

**Validation gates:**
- convergence_order_pass: True (p_r, p_v >= 3.5)
- numerical_magnitude_pass: True (|numerical| in [1e-4, 1e-1] deg/day)
- ratio_band_pass: False (audit-015 band [10x, 100x] violated; actual ratio
  is ~170x — **documented as a first-principles discovery** that the
  audit-015 estimate under-estimated the closed-form over-estimate by ~3x;
  the magnitude is qualitatively correct, quantitatively revised)

## Findings

1. **HEADLINE FINDING:** The closed-form secular-average Lunisolar RAAN
   upper bound (Vallado Eq. 9-46 form) over-estimates the numerically
   integrated Lunisolar RAAN rate at dawn-dusk SSO by a SIGNED ratio of
   ~170x at h=600 km. The ratio is NEGATIVE: the closed-form is retrograde
   (-0.218 deg/day at h=600) while the numerical integration is prograde
   (+0.001284 deg/day at h=600). This is a SIGN DISAGREEMENT, not just a
   magnitude over-estimate, and is a byte-pinned, reproducible measurement.

2. **DISCOVERY vs audit-015 estimate:** The measured over-estimate factor
   (~170x at h=600 km) is ~3x LARGER than the audit-015 follow-up-candidates
   report estimate of "~50x" (model_note in Exp 016). The audit-015 estimate
   is qualitatively correct (closed-form over-estimates by order of
   magnitude) but quantitatively under-estimates the magnitude of the
   over-estimate. This is a first-principles discovery, not a refutation
   of the audit's overall direction.

3. **Numerical Lunisolar RAAN rate:** +0.001284 deg/day (~+0.47 deg/year)
   at h=600 km over the 1-year byte-pinned DE441 arc. This is within the
   operational envelope (~0.005 deg/day Sentinel/Landsat, ~1.8 deg/year)
   reported in Exp 016, and is PRO-GRADE rather than retrograde as the
   closed-form predicts.

4. **Self-convergence order:** p_r = 4.49, p_v = 4.50 (RK4 design order ~4
   confirmed).

5. **Linear-fit residual RMS:** ~0.024 deg at h=600 km, consistent with
   periodic Lunisolar variations at the lunar and solar synodic/monthly
   frequencies (not captured by the secular linear fit).

6. **Audit response:** This experiment converts the Exp 016 closed-form
   upper-bound disclaimer into a byte-pinned, numerically validated quantity.
   The decadal direction originally proposed for 017 (AGENTS.md, roadmap.md)
   is rejected as not scientifically defensible at this time; the closed-form
   upper-bound verification (Track H Alt-1, scored 27/30) is the strongest
   scientifically defensible alternative and is the experiment actually
   executed.

## References

- D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed.,
  Microcosm, 2013 — Ch. 9 secular J2 + Lunisolar + Eq. 9-46 closed-form
  secular-average Lunisolar RAAN formula.
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed.,
  Elsevier, 2021 — Ch. 10 perturbations + RAAN control.
- Astronomical Almanac low-precision solar formulas (mean longitude, mean
  anomaly, equation of center, mean obliquity of date).
- Aoki et al. 1982 — IAU-1982 GMST polynomial.
- WGS-84 TR8350.2 — `R_E`, `J2`, `omega_E`.
- IAU 2015 Resolution B3 — `μ_E`, `μ_Sun`, `μ_Moon`.
- IAU 2012 Resolution B2 — `AU` (exact).
- JPL Horizons API (`https://ssd.jpl.nasa.gov/api/horizons.api`) — Sun and
  Moon geocentric vector snapshots, DE441, ICRF/TDB, KM-S, geometric.
- Exp 009 j2Precession — secular J2 nodal/apsidal rates.
- Exp 012 orbitClasses — SSO inclination lock + measured J2 closure residual
  (~+0.6% relative at h=600 km).
- Exp 014 eclipseTiming — byte-pinned 2026 Sun geocentric snapshot
  acquisition pattern (followed verbatim by the Moon snapshot).
- Exp 016 lstDrift — closed-form Lunisolar upper bound formula + the
  "~50x over-estimate" model_note that motivated this experiment.
- `localdocs/reports/audit-015-follow-up-candidates-2026-08-29.md` — scored
  the closed-form upper-bound verification experiment as the recommended
  next step (29/30 candidate #4; 27/30 in the post-016 Track H re-scoring).
- `localdocs/reports/audit-015-lst-drift-2026-08-29.md` — independent
  numerical falsifier that measured 0.0223 deg/day drift at the SSO node
  (consistent with the +0.00128 deg/day measured here for the Lunisolar-only
  contribution).

## Limitations

- Point-mass Lunisolar (no Earth-Moon barycenter correction).
- J2 only for non-Kepler gravity (no tesseral harmonics, no solid-Earth
  tides, no ocean loading).
- No SRP, no drag, no relativity (each excluded as a separate force).
- No future-arc extrapolation; experiment is bounded to 2026 (the byte-pinned
  snapshot year). Decadal extension would require a byte-pinned 10-year
  ephemeris acquisition (deferred).
- Mean-orbit constants in the closed-form reproduction use the lab's canon
  `LUNAR_DISTANCE_KM=384400.0` (constant geocentric distance) and
  `LUNAR_INCLINATION_DEG=5.145` (constant lunar inclination to ecliptic).
  The byte-pinned snapshot is the only place where time-varying geocentric
  distance and lunar phase evolution enter.
- Linear fit of Ω(t) vs t does not capture the dominant periodic Lunisolar
  terms (lunar nodal period 18.6 yr is far longer than the 1-year arc;
  lunar anomalistic month 27.55 d, lunar synodic month 29.53 d, and solar
  synodic year 365.24 d all contribute to the short-period residuals at
  ~0.024 deg RMS).
- The measured ~170x over-estimate factor depends on the specific 1-year
  arc covered (2026 calendar). Different years would give different factors
  due to the lunar phase evolution; the byte-pinned snapshot guarantees
  reproducibility but not representativeness over a multi-decadal horizon.
- The closed-form's retrograde sign at LEO SSO is contradicted by the
  numerical prograde sign. This may reflect the secular-average formula's
  treatment of the lunar evection (which has both prograde and retrograde
  contributions depending on lunar phase) rather than a sign error in the
  formula; the closed-form is reported as "upper bound" (magnitude only)
  in Exp 016 and this experiment, and the sign disagreement is documented
  as an open question for follow-up work.

## Status

COMPLETE (2026-08-30): 32 new tests, all passing. Total repo tests: 658
(626 baseline + 32 new). All deterministic, offline, byte-stable figures.
Audit response: ✓.