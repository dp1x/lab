---
tags: [orbital-mechanics, drag, orbit-decay, non-conservative, validation-doctrine]
date: 2026-08-22
aliases: [atmospheric drag, orbit decay]
links: ["[[j2-precession]]", "[[ground-tracks]]", "[[gravity-assist]]"]
---

# Orbit Decay / Atmospheric Drag

## Summary

Experiment 010 added the lab's first **non-conservative force**. Energy conservation
is replaced by **energy monotonicity + pointwise dissipation accounting**:
`d(eps_total)/dt = a_drag . v = -(1/2) kappa rho(h) |v_rel| (v_rel . v) <= 0`, with
`v_rel = v - omega_atm x r`. Circular decay obeys `da/dt = -kappa rho sqrt(mu a)`
with an exact erfi closed form per exponential layer; eccentric orbits follow
Gauss/King-Hele per-revolution deltas (apogee absorbs the decay; King-Hele few-%
forms only for e <= 0.1). Everything was validated against five independent oracle
paths and adversarial mutants.

## Content

### The convention trap that almost shipped

The handoff contract wrote `B = C_D A/m` but swept `{50..400 kg/m^2}` — magnitudes
that only exist for the INVERSE convention `beta = m/(C_D A)`. A real ISS-scale
spacecraft has `m/(C_D A) ~ 190 kg/m^2` but `C_D A/m ~ 0.005 m^2/kg`. Freezing the
wrong one is a 4x10^4 rate error that still "looks like decay". Lesson: **unit-check
the sweep magnitudes against physical spacecraft, not just the formula.**

### Estimator lesson: raw osculating elements lie under J2

Switching J2 on does NOT leave the measured decay rate unchanged — and that is not a
bug. Seeding a Kepler state relaxes mean elements by `(2a^2/mu)<U_J2>` (~5.7 km at
i=51.6°, because the inclination-dependent mean of the J2 potential is nonzero), so
decay proceeds in denser air: settled energy rates +12.3%. Raw osculating-a slopes
additionally suffer short-period ripple aliasing (block slopes swung -78..+71 m/day).
The honest comparison uses **energy-bookkeeping tail rates** with the altitude-shift
factor as prediction: measured/predicted = 1.1225/1.1229 (residual -0.04%).
Generalizes the [[j2-precession]] mean-vs-osculating doctrine to non-conservative forces.

### Superconvergence of crossing-type observables

The time-to-fall observable converged at orders 3.8 → 5.0 (above design order 4)
before hitting the floor: leading RK4 error terms cancel in phase/crossing errors.
Raw position error stayed >= order 4.4. Lesson: verify integrator order on a plain
observable too; a "too good" convergence slope is a property of the observable, not
the integrator.

### Co-rotation magnitude (handoff correction)

The handoff claimed atmosphere co-rotation shifts LEO drag by "~2-7%". Verified
physics: wind factor `(1 -/+ w)^2` with `w = omega_E r/v = 0.064` at 400 km equator
→ **±12-13% in drag magnitude, ±18-21% in dissipation rate** (prograde less,
retrograde more); for inclined circular orbits the first-order correction is
`w cos i`. Measured equatorial twin ratio 1.303 vs exact 1.295.

### Doctrine extensions (non-conservative validation)

- Pointwise identity residuals must be RELATIVE (an absolute bound inherited from
  conservative experiments was 17x larger than the signal).
- Monotonicity needs TOTAL energy when conservative perturbations are active
  (J2 trades Kepler energy back and forth; total including U_J2 is strictly decreasing).
- Plateau separation still works: law swap moves transit times ~103,000 s (flat in
  dt to 0.73%) while refinement contributes 2,781 s → 37x separation.
- Mutant harnesses (sign flip, kg->km^3, B inversion) prove the DETECTORS fire;
  each mutant produces a plausible-looking decaying orbit.

## Source Experiments

- [[j2-precession]] — generalized-force propagator pattern cloned here; bit-exact regressions
- ground tracks ([[ground-tracks]]) — Kepler machinery and constants
- ode-integration-basics — RK4 order discipline

## Key Takeaways

1. First non-conservative force: conservation contracts become monotonicity +
   dissipation accounting; identities must be evaluated on total energy.
2. Ballistic coefficient conventions are a 4x10^4 trap; freeze one and test the sweep
   magnitudes against real spacecraft.
3. erfi closed form exists for single exponential layers (h0 >= 120 km in float64);
   layer-joint-aligned quadrature is the general primary oracle.
4. Osculating elements under J2 carry a constant inclination-dependent offset — never
   fit secular rates to raw osculating series across a model switch.
5. Decade-wide benchmark bands only: quiet-time ISS decay (~2 km/month) binds;
   storm-time numbers are context, not validation.

## See Also

- [[j2-precession]] (perturbation doctrine predecessor)
- [[ground-tracks]]
- roadmap row 011 (Lagrange points) — next experiment
