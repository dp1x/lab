---
tags: [orbital-mechanics, orbit-classes, sun-synchronous, molniya, gto, geo, classification, j2, resonances]
date: 2026-08-23
aliases: [orbit-classes, constraint-defined-families, sso-existence-limit]
links:
  - "[[j2-precession]]"
  - "[[ground-tracks]]"
  - "[[hohmann-transfer]]"
  - "[[kepler-orbit-validation]]"
  - "[[combined-transfer-plane-change]]"
---

# Orbit Classes as Constraint-Defined Families

## Summary

Exp 012 reframed the classic "orbit catalogue" as a **classification problem**: each
operationally important Earth-orbit class is the solution set of coupled dynamical
constraints under a declared model (two-body + first-order secular J2 + spherical
Earth rotating at omega_E). Three mechanism classes emerge:

1. **Analytic locks** — an element is pinned by a secular-rate condition:
   SSO `Omega_dot(a,e,i) = +sun rate` gives `cos i_SSO = -(a/a_max)^(7/2)`; Molniya
   `omega_dot = 0` gives `cos^2 i = 1/5`; GEO `n = omega_E` gives a = 42164.169 km.
2. **Resonances** — Diophantine clock conditions `m T = k P_sidereal` produce a
   discrete lattice of repeat radii (GEO 1:1 = 42164 km, Molniya 2:1 = 26561.762 km,
   3:1 = 20270 km, ...).
3. **Existence boundaries** — demand can exceed supply: SSO exists only for
   `a <= a_max = (1.5 J2 sqrt(mu) R^2 / lambda)^(2/7) = 12352.505076 km`
   (h_max = 5974.368 km at i -> 180 deg, NOT 90 deg); eccentricity extends the
   limit via `(1-e^2)^(-4/7)`.

## Content

### The convention minefield (all machine-verified)

- The sun-tracking target is the MEAN-SOLAR-year rate `360/365.2422 =
  0.985647332099 deg/day`. Sidereal year (0.98560912) and Julian year (0.98564685)
  are wrong targets separated in i_SSO(600 km) by 3.0e-4 and 1.7e-4 deg — caught
  analytically with a 5e-5 deg tolerance. The tropical variant (365.24219 d) differs
  by only 2.1e-7 deg: behaviorally indistinguishable, pinned by literal instead.
- Confusing Earth's rotation (360.98565 deg/day) with the sun's apparent rate
  (0.98565 deg/day) drives cos i to -49.6: structurally outside [-1,1] — the solver
  must return a typed no-solution sentinel, never clip.
- omega_E is the MASTER sidereal constant; every resonance condition is written
  against `P_sidereal = 2*pi/omega_E`, never against a hard-coded day length.

### What full-force propagation adds near the critical inclination

The J2-on Kepler period excess of a Molniya-class orbit (e = 0.74, i = i_crit) is
**+323.0 s/orbit** beyond the Kepler value — plateaued under step refinement,
energy-conserving to the integrator floor, with osculating-a excursions ~+160 km.
This is small-divisor-amplified short-period dynamics that NO first-order secular
clock predicts (~+4 s). Consequences worth remembering:

- Millisecond-scale draconitic/anomalistic splits (+24.06 ms first-order at
  i = 63.4 deg; the naive `2pi/(n -/+ omega_dot)` form giving 48 ms is factor-2
  wrong because M-dot carries its own J2 term) are invisible to event timing.
- The claimable freeze evidence lives on orbit-averaged ELEMENT-regression paths
  ([[j2-precession]] doctrine): omega-dot at the lock measured -5.97e-05 deg/day.
- Event longitudes track the OSCULATING apsis, which sweeps ~360 deg/orbit; its
  measurable identity `Delta-lambda/orbit = 360 + (Omega_dot + omega_dot - omega_E)*T_days`
  held to 0.001 deg/orbit (355.460 vs 355.465 deg/day over a 12.5-day arc).

### Class numbers worth citing (pinned constants, machine precision)

| Quantity | Value |
|---|---|
| i_SSO @ 500/600/800 km (e=0) | 97.401786 / 97.787647 / 98.603085 deg |
| SSO existence limit | a_max = 12352.505076 km (h_max 5974.368 km) |
| Critical inclination | 63.43494882 deg (supplement 116.56505118 deg) |
| Semi-synchronous radius | 26561.762328 km (n = 2 omega_E) |
| Molniya apsides (e=0.74) | h_p/h_a = 527.92 / 39839.33 km |
| Apogee dwell (+/-90 deg, e=0.74) | f = (pi - E1 + e sin E1)/pi = 0.923607 |
| GEO radius | 42164.169462 km (h = 35786.032 km); omega_dot = -2 Omega_dot at i=0 |
| GEO stationarity residual (Keplerian radius) | +0.02683 deg/day NONZERO (negative control) |
| GTO budget (300 km -> GEO, coplanar) | dv1 = 2.42573, dv2 = 1.46682, total 3.89256 km/s, tof 5.275 h |
| Repeat-corrected Molniya radius (first-order) | 26553.420405 km (M_dot + omega_dot = 2(omega_E - Omega_dot)) |

### Validation-shape lessons

- Survivor-class mutants must be pre-registered with paired discriminators:
  omega-tests are blind to branch flips (5cos^2 i - 1 even in cos i); total-dv asserts
  are blind to burn swaps; frozen-null tests are blind to p:=a substitution; NO numeric
  threshold catches the tropical-year swap (documented blindness + constant pinning).
- Coarse grids quantize parabolic-refined crossing times to the sample step (observed
  lock artifact at GEO/512 spp): event-timing claims need resolution + multi-orbit
  averaging, or honest report-only floors.

## See Also

- [[j2-precession]] — secular rates and the estimator/oracle contract consumed here.
- [[ground-tracks]] — repeat-track closure machinery this experiment generalizes to m:n.
- [[hohmann-transfer]] — vis-viva budgets reused for the GTO class (3.9319 km/s anchor).
- [[combined-transfer-plane-change]] — where the GTO plane-change coupling lives.

## Status

Experiment 012 COMPLETE (2026-08-23): 43 tests, six figures, deterministic
double-run verified. Shared machinery graduated: `j2_rhs` in `src/lab_utils/orbits.py`.
