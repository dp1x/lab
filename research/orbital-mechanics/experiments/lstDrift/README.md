# Experiment 016 — SSO LST-Drift Correction: First-Principles EoT Envelope + Multi-Year Perturbation Budget

> Status: in progress
> Date: 2026-08-30
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/lstDrift/`

## Research Question

What is the actual local-solar-time (LST) drift at the ascending node of a
true dawn-dusk SSO at h in {500, 600, 700, 800} km over a 1-5 year mission,
and what is the corresponding station-keeping Δv budget?

The question decomposes into four quantitatively testable sub-questions:

1. **First-order J2 + mean Sun lock**: does `dLST/dt = (dΩ/dt − dα_sun/dt)/15`
   vanish to leading order for an orbit with `dΩ/dt = SSO_TARGET_DEG_DAY`,
   and what is the residual once the lab's mean-of-date Almanac Sun model
   is used (not the idealized mean-sun-rate Sun)?
2. **Equation-of-Time envelope**: what is the peak-to-peak LST variation
   driven by the obliquity + eccentricity of Earth's orbit (the EoT
   envelope, ~±16 min) at an SSO node, and how does it compare with the
   byte-pinned 2026 Horizons Sun snapshot from Exp 014?
3. **J2 closure residual**: what is the `~0.6%` first-order vs higher-order
   J2 bias documented by Exp 009 / Exp 012, and what is the corresponding
   LST drift at an SSO node (over 1-5 year arcs)?
4. **Luni-solar + SRP + drag perturbations**: what is the LST drift
   driven by the lunar+solar third-body perturbation, solar radiation
   pressure, and atmospheric drag (at h in {500, 600, 700, 800} km),
   decomposed into secular, long-period, and short-period contributions?

The expected dominant secular LST drift is the J2 closure residual
(~2.2 deg/year = ~9 min/year) plus Lunisolar (~few minutes/year) plus
drag (altitude-dependent). The total should be of order 10-30 min/year
(consistent with Sentinel-1/Landsat operational envelopes of
~5-15 m/s/year station-keeping).

## Frozen Contract v1.0

| Item | Value | Provenance |
|---|---|---|
| `R_E` (km) | 6378.137 | WGS-84 equatorial; lab canon `R_EARTH_KM` |
| `J2` | 1.082629821e-3 | WGS-84, J2 = √5·\|C20_bar\|; lab canon |
| `μ_E` (km³/s²) | 398600.4418 | IAU 2015 nominal GM_E; lab canon |
| `ω_E` (rad/s) | 7.2921159e-5 | WGS-84 / Vallado Table 3-1 |
| SS0 target (deg/day) | 360/365.2422 = 0.985647332099 | Exp 012 pinned |
| LST target (h) | 18.0 (dusk-ascending terminator) | Exp 015 frozen |
| Altitudes (km) | {500, 600, 700, 800} | Exp 015 frozen band |
| Mission durations (years) | {1, 3, 5} | Pre-registered |
| EoT formula | Geometric Almanac mean-of-date (Exp 014 frozen) | lab canon |
| Sun reference | 2026 Horizons snapshot byte-pinned (`eclipseTiming/reference/`) | Exp 014 |
| Validation gate | EoT peak-to-peak vs Horizons snapshot | Pre-registered ±0.7 deg |
| Lunisolar | Point-mass third-body on circular lunar/solar orbits | Standard |
| SRP | Cannon-ball, A/m = 0.01 m²/kg (default), Earth-shadow cylindrical | Pre-registered |
| Drag | Exponential atmosphere, scale height H=8 km, base ρ = exp(-(h-500)/H*ln(10)) (fiducial) | Pre-registered |
| Reference Δv anchor | Sentinel-1 ~15 m/s/year, Landsat-7/8 ~5-15 m/s/year | Public-domain literature |

## Methodology

Deterministic, offline-only after acquisition of the 2026 Horizons Sun
snapshot (already in the repo from Exp 014). No network at runtime, no
RNG, no wall-clock in the analysis.

1. **EoT envelope (first-principles)**: evaluate `lst_at_orbit_node_at_t`
   at 10000 ascending-node crossings over 1 year; compute the peak-to-peak
   range; compare to the byte-pinned 2026 Horizons Sun snapshot
   (apparent-vs-mean residual = EoT).
2. **J2 closure residual**: propagate the closed-form mean-element Ω_dot
   (Exp 009 formula) at h in {500, 600, 700, 800} km; compare to
   SSO_TARGET_DEG_DAY; report the residual and the corresponding LST
   drift at 1, 3, 5 year arcs.
3. **Lunisolar perturbation budget**: analytical third-body formulas
   (Vallado Ch. 9); compute the secular and long-period RAAN perturbation
   for a 5-year arc at each altitude; report the corresponding LST drift.
4. **SRP perturbation budget**: cannon-ball model; compute the secular
   RAAN perturbation; report the LST drift at each altitude.
5. **Drag perturbation budget**: exponential atmosphere; compute the
   altitude-dependent RAAN perturbation (drag-induced Ω change coupled
   to J2); report the LST drift at each altitude.
6. **Total secular LST drift budget**: sum all contributions with
   conservative sign (worst-case envelope); compare to operational
   Sentinel-1/Landsat envelopes (~5-15 m/s/year).
7. **Station-keeping Δv**: closed-form RAAN-control Δv at each altitude
   for each 1, 3, 5 year arc; report as `Δv_5yr ~ 5*15 = 75 m/s/altitude`,
   compared to the 5-15 m/s/year operational anchor.

## Implementation

- Script: `experiment.py` (deterministic, offline)
- Language/runtime: Python 3.12, numpy, matplotlib Agg
- Runtime: < 5 min single core (no iterative sweeps; analytical formulas)
- Determinism: pure float64, no RNG, no network at runtime, no wall-clock
  in the analysis. `time.time()` only in `run()` for elapsed prints.

## Validation Method

Six layers (target ~30-50 new tests):

- **L1 closed-form identities**: LST formula at orbit-plane node matches
  `12 + (Ω - α_sun)/15` to 1e-12; GMST cancellation bit-equivalence.
- **L2 EoT envelope**: peak-to-peak range vs Horizons snapshot to ±0.7 deg
  band; signed vs unsigned; 2026 case pinned.
- **L3 J2 closure**: first-order J2 vs Sun rate residual at h in {500,600,
  700,800}; reproduce Exp 012 headline to 1e-3.
- **L4 Lunisolar + SRP + drag**: secular rate formula vs Vallado/Curtis
  reference values; altitude scaling law; drag-exponential consistency.
- **L5 station-keeping Δv**: closed-form RAAN-control Δv formula;
  orbit-plane maneuver magnitude; sign convention.
- **L6 adversarial mutants**: negated lunar node, swapped lunar/solar,
  EoT cleared, sign-flip Δv, signed vs unsigned envelope.

## Headline numbers (to be pinned in `results.json`)

| Quantity | Expected | Pre-registered band |
|---|---|---|
| EoT peak-to-peak | ~24 min | ±5 min around textbook |
| J2 closure residual | ~2.2 deg/year | within 1% of Exp 012 |
| Lunisolar LST drift (5 yr) | ~few min/year | ±50% of textbook |
| SRP LST drift (5 yr) | <1 min/year | within 2x of textbook |
| Drag LST drift (5 yr, h=600) | ~few min/year (downward) | depends on density |
| Total secular LST drift (5 yr, h=600) | ~10-30 min/year | within operational envelope |
| Station-keeping Δv (5 yr, h=600) | ~50-100 m/s | within 2x of operational 5-15 m/s/year |

## References

- D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed.,
  Microcosm, 2013 — Ch. 9 secular J2 + Lunisolar + SRP + drag.
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier,
  2021 — Ch. 10 perturbations + RAAN control.
- Bate/Mueller/White, *Fundamentals of Astrodynamics*, 1971 — Ch. 9 perturbations.
- Astronomical Almanac low-precision solar formulas (mean-of-date geometric Sun).
- Aoki et al. 1982 — IAU-1982 GMST polynomial.
- WGS-84 TR8350.2 — `R_E`, `J2`, `ω_E`.
- IAU 2015 Resolution B3 — `μ_E`.
- Sentinel-1 / Landsat flight dynamics reports (ESA Copernicus / NASA EOSDIS;
  reference anchor only).
- Exp 009 (j2Precession) — secular J2 rates + orbit-specific anchors.
- Exp 012 (orbitClasses) — SSO closure residual, finite-existence a_max.
- Exp 014 (eclipseTiming) — conical shadow model, byte-pinned Sun snapshot.
- Exp 015 (dawnDuskSSO) — multi-constraint mission analysis + corrected
  LST-drift narrative (2026-08-29 remediation).

## Limitations

- Spherical Earth, no tesseral harmonics.
- Point-mass Lunisolar (no Earth-Moon barycenter correction; out of scope
  for LST drift magnitude).
- Exponential atmosphere (no F10.7 / Jacchia-Bowman; out of scope for
  LST drift magnitude).
- No third-body Sun tidal torque on Earth (luni-solar solid-Earth tides).
- No relativity.
- Station-keeping Δv assumes impulsive burns at line of nodes (Vallado 8.5);
  finite-burn and plane-of-sky maneuvers are out of scope.

## Status

In progress; see `results.json` for the frozen output once complete.