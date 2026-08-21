---
tags: [orbital-mechanics, j2, precession, secular-perturbation, mean-vs-osculating, rk4, validation]
date: 2026-08-22
aliases: [j2-precession, nodal-precession, apsidal-precession, secular-rates]
links:
  - "[[ground-tracks]]"
  - "[[kepler-orbit-validation]]"
  - "[[combined-transfer-plane-change]]"
  - "[[ode-integration-basics]]"
---

# J2 Secular Precession: Making a Propagator Rediscover Theory

## Summary

Earth's oblateness (J2 = 1.082629821e-3, WGS-84 via sqrt(5)|C20_bar|) causes first-order
secular nodal and apsidal drift

```
Omega_dot = -(3/2) n J2 (R_E/p)^2 cos i          (westward for prograde)
omega_dot =  (3/4) n J2 (R_E/p)^2 (5cos^2 i - 1) (zero at i = arccos(1/sqrt(5)) = 63.4349 deg)
n = sqrt(mu/a^3), p = a(1-e^2)
```

Exp 009 established these analytically AND forced a full-force Cowell RK4 propagator
(`a = -mu r/r^3 + a_J2`) to rediscover them from raw Cartesian states: osculating elements
per sample -> unwrap -> integer-orbit-window least squares. Anchors reproduce to
0.42-0.68% (ISS -4.972 vs -4.951 deg/day; SSO600 at solved i_SSO=97.787647 deg hits the
solar-rate target +0.98565; polar null exact to machine precision via node crossings).
The residual is NOT integration error: it is stable across 20/50/100-orbit windows,
does not shrink when the timestep halves, and has one uniform sign — it is the
first-order model-order gap (mean-vs-osculating element offset delta-a/a ~ -1.2e-3 plus
second-order small-divisor terms near critical inclination). Raw RK4 order was proven
separately (4.09 via closed-form Kepler truth); the rate metric converges FASTER than
h^4 (orders 4.5-4.7) because orbit-averaging cancels the integrator's leading
phase-error mode.

## Content

### The estimator/oracle independence contract

The portfolio failure mode is "validation" where both sides share algebra. Exp 009 keeps
them disjoint by construction:

* Oracle input: `(a, e, i)` only. Estimator input: propagated Cartesian states only.
* Estimator geometry: `h = r x v`, node `= zhat x h`, `e_vec = (v x h)/mu - r/r`;
  generic closed-form OLS. No orbital-rate formula appears on the numerical side.
* Third path: regress inertial longitude AT ascending-node crossings (parabolic time
  refinement). First-order Omega short-period terms are ~cos(u) and vanish at the nodes,
  so this estimator gives CLEAN secular nulls — polar orbit measures Omega_dot ~ 1e-16
  deg/day, machine zero. It agrees with the element-fit estimator to ~8e-4 rel.
* Symmetry breakers that catch shared bugs: J2=0 null (frame bugs leak omega_E ~ 0.986
  deg/day artifacts), J2->-J2 flip (both rates must reverse), wrong-p discriminator
  (Molniya: using a instead of p shifts Omega_dot by (1-e^2)^{-2} = 4.89x).

### Sign/convention traps (each with a one-line test)

* Prograde Omega_dot < 0 (node regresses west); retrograde > 0; polar exactly 0
  (float-eval zero: cos(pi/2) = 6.1e-17, so test |x| < 1e-12, never == 0.0).
* omega_dot is mirror-symmetric about 90 deg — EQUAL at 90+d and 90-d (cos^2 enters),
  sign flips happen at 63.43/116.57 deg, NOT at polar.
* e enters ONLY through p: rate(e)/rate(0) = (1-e^2)^{-2} exactly. Same-p/different-a
  orbits have DIFFERENT rates (n depends on a) — a tempting but wrong "invariance" test.
* GEO nodal drift is -0.0134 deg/day (~-4.9 deg/yr): tiny but not negligible;
  don't set thresholds below real physics.

### Mean vs osculating is THE distinction

The double-averaged formulas consume MEAN elements; any state-based measurement yields
OSCULATING elements. Short-period oscillations average out over integer-orbit windows
(leakage bias ~ 1.91*A/N worst-case, A = short-period amplitude), but a constant offset
remains: evaluating the oracle at initial osculating elements differs from true mean
elements by O(J2(R/p)^2) in a — measured here as a systematic +0.42..0.68% relative rate
offset with uniform sign across all seven cases (implied delta-a/a ~ -1.2e-3 via
delta-Omega_dot/Omega_dot = -3.5 delta_a/a). Near critical inclination the second-order
secular terms carry small divisors (1-5cos^2 i) and amplify (Molniya 0.91%, CRITICAL 0.68%).
The plateau proof separates this from integration error: halving dt changes the residual
<50% while error-vs-reference metrics shrink 16x+ per halving.

### What "circular" means under J2

A seeded e=0 LEO develops a REAL induced eccentricity: radial excursion measured
+/-9.60 km at a=6878 km — exactly a*(3/2)J2(R_E/p)^2, from the effective-mu deficit of
the equatorial J2 term. Consequences: e_vec sits ~1e-3 above zero (far above any
reasonable guard), so omega stays numerically finite but sweeps once per orbit with no
secular content; the recovered "omega_dot" equals the mean motion (~5589 deg/day for ISS).
Claims policy needed: omega_dot claimed only for seed e >= 0.01; results.json carries an
explicit note field wherever the trap applies. At i = 0/180 RAAN is structurally
undefined (node vector vanishes) and omega — measured FROM the node — dies with it:
NaN sentinels, absolute-not-relative residuals against structural zeros.

### Reuse pattern

Local generalized RK4 whose j2==0 branch executes the identical float-op sequence as the
verified Exp 006 propagator -> regression test asserts np.array_equal bit-for-bit. Prior
experiment code stays untouched; verification transfers through equivalence. Element
helpers import single-hop from Exp 008; convergence_rate eps-guarded with
np.maximum(err, 1e-16) (lab_utils raises on nonpositive errors).

## Source Experiments

* `research/orbital-mechanics/experiments/j2Precession/` — runnable card, 32 tests,
  deterministic 53 s run, figures regenerate from results.json.
* [[kepler-orbit-validation]] — RK4 foundation, periapsis-resolution step law reused.
* [[combined-transfer-plane-change]] — verified 3-D Cowell propagator (bit-exactness donor).
* [[ground-tracks]] — constants/frame conventions and importlib reuse pattern; its
  "Limitations" deferred J2 to this experiment.

## Key Takeaways

* Separate the two questions always: (A) does the numerics converge (vs NUMERICAL
  reference, order ~4), and (B) how far is the converged answer from the analytical model
  (model-order residual, reported, never called error). Comparing finite-step numerics
  directly against analytics conflates them and produces fake convergence failures or fake
  agreement.
* A superconvergent derived metric is normal, not suspicious: averaging estimators cancel
  the integrator's leading error mode. Prove integrator order on a phase-sensitive
  full-vector metric vs closed-form truth; final-|r| alone hides along-track phase error
  and decays an order too fast.
* Node-crossing regression is a free, algebra-independent third estimator with exact nulls
  for symmetric configurations; use it as the symmetry breaker against shared-algebra bugs.
* Even-in-J2 physics breaks naive symmetry tests: flipping J2 flips rates by -1.0088, not
  -1.0000 (second-order secular terms). Tolerance bands must leave room for real physics,
  decided from evidence arithmetic, never post-hoc fitted to pass.
* Structural zeros need absolute residuals; structurally undefined elements need sentinels
  plus claims policy, not fabricated numbers. Store deg/day twins for small rad/s values
  (JSON rounding band hazard in save_json_result for values in [1e-10, 1e-7)).
* Anchors computed FROM pinned constants (SSO inclinations solved live) beat hard-coded
  table values: they stay consistent if constants ever change and make the provenance auditable.

## See Also

* [[ground-tracks]] — ECI/ECEF machinery, sidereal conventions, repeat tracks.
* [[kepler-orbit-validation]] — Kepler solution machinery, RK4 order methodology.
* [[combined-transfer-plane-change]] — 3-D Cowell propagation, rotation matrices.
* [[ode-integration-basics]] — integrator orders, symplectic vs non-symplectic behavior.
* [Vallado, Fundamentals of Astrodynamics and Applications, 4th ed., Ch.9]
* [Curtis, Orbital Mechanics for Engineering Students, 4th ed., Ch.10]
* [Bate, Mueller, White, Fundamentals of Astrodynamics, Ch.9]
* [Murray & Dermott, Solar System Dynamics, Ch.6]
