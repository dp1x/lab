---
tags: [orbital-mechanics, manuevers, transfers, hohmann, delta-v]
date: 2026-08-13
aliases: [hohmann-transfer, least-fuel-transfer]
links:
  - "[[kepler-orbit-validation]]"
  - "[[kepler-equation-solvers]]"
  - "[[ode-integration-basics]]"
---

# Hohmann Transfer: Least-Fuel Two-Impulse Transfer Between Circular Orbits

## Summary

Between two coplanar circular orbits of radii r1 < r2 around mass μ, the
minimum-delta-v two-impulse transfer is the Hohmann transfer: two tangential
burns at r1 and r2 on the half-ellipse with a_t = (r1+r2)/2, e_t = (r2−r1)/(r2+r1).
Measured (μ = r1 = 1): Δv_total = (√(2R/(1+R)) − 1 + √(1/R)(1 − √(2/(1+R)))) · v1
with R = r2/r1. On the verified RK4 machinery (Experiment 002), the complete
burn → coast → burn transfer closes the target circle to relative radius error
4.5e-10–4.2e-9 and 4.5e-11 in circularity across R ∈ {1.5, 6.41, 20}, a
Venus-like ratio (R = 1.3825) and a true inward transfer (R = 0.5).

## Key quantitative facts (verified in Experiment 004)

- **Closed form**: Δv1 = v1(√(2R/(1+R)) − 1), Δv2 = v2(1 − √(2/(1+R))),
  t_tr = π√((r1+r2)³/(8μ)) — re-derived from vis-viva and Kepler III, matched
  to 1e-9 by propagation.
- **Cost is NOT monotone**: interior max at R* = 15.5817, Δv_total/v1 = 0.536258.
  Past R* targets get cheaper in fuel but strictly more expensive in time on
  Hohmann paths (time grows ∝ R^{3/2}).
- **Asymptotes**: Δv_total/v1 → (R−1)/2 as R → 1 (measured 0.9999× at R = 1.0001);
  → √2−1 = 0.414214 (the escape burn) as R → ∞ (measured rel gap 2.4e-6 at
  R = 1e12). Reaching "infinity" costs exactly the escape burn. Δv1 has no burn
  above the escape value (√2 − 1)v1 anywhere; the circularization burn Δv2 peaks
  at 0.1900·v1 near R = 5.88 and vanishes at both ends.
- **Inward transfers**: exact burn-magnitude symmetry —
  Δv_total(r1, r2) = Δv_total(r2, r1) (numbers identical in the experiment); the
  textbook formula set assumes r2 > r1, the inward case swaps burn order. When
  validating a true inward flight (r2 < r1) against the closed-form ellipse,
  remember kepler-style solutions start at periapsis while the flight starts
  at apoapsis: compare against the same ellipse phase-shifted by a half-period
  (the leg actually flown), or the analytic arrival metric is meaningless.
- **Optimality**: within the two-impulse class with burns at r1, r2 (transfer
  ellipse r_p ≤ r1, r_a ≥ r2, tangency NOT assumed), the Hohmann corner is the
  strict minimum (121×131 grid gap ≤ 7.8e-16 at R ∈ {2, 6.41, 20}). Caveat:
  three-impulse bi-elliptic transfers beat Hohmann for R ≳ 11.94 when flight
  time is unconstrained — Experiment 005.
- **Real numbers** (IAU 2015 B3 nominal + JPL mean orbits; matches public
  canonical values): LEO(200 km)→GEO = 3.9319 km/s in 5.26 h (published ≈ 3.9 km/s);
  Earth→Mars = 258.87 d, v∞ 2.945 / 2.649 km/s, TMI from LEO = 3.6114 km/s
  (published ≈ 3.6, "0.4 km/s over escape"); Venus inward = 146.08 d,
  TMI 3.5036 km/s. Single-stage Isp 300 s LEO→GEO propellant fraction 0.737.

## Connections

- **[[kepler-orbit-validation]]** — the burn states, arrival checks and
  post-burn orbit verification reuse its validated RK4 propagator and
  periapsis-resolved step-count law ((1−e)^{−3/2}); transfer ellipse elements
  recovered via its Kepler-solution machinery.
- **[[kepler-equation-solvers]]** — inverse position along the transfer ellipse
  (position-on-ellipse at burn/arrival points) routes through Kepler solvers.
- **[[ode-integration-basics]]** — energy conservation/h integration bounds
  the propagation errors seen here (vector magnitude drift < 1e-9 over a
  full ellipse).
- Next in chain: **bi-elliptic vs Hohmann** (roadmap 005) — crossover ratio
  11.94/15.58 located geometrically; the same cost machinery extends to plane
  changes (006) and gravity assists (007).

## Reusable formula card

- Digit-safe near R = 1 (avoids catastrophic cancellation):
  √(2R/(1+R)) − 1 = (R−1)/((1+R)(1 + √(2R/(1+R)))); likewise for 1 − √(2/(1+R)).
- Hyperbolic excess for injection from low orbit (r_c, v_circ):
  v_inj = √(v_circ² + v∞²) — used to get 3.6114 km/s (Mars) from 2.9447 km/s v∞.
- Fuel fraction: m_prop/m_0 = 1 − exp(−Δv/(Isp·g0)) (units: km/s → m/s factor 1000).

## Open questions

- Where exactly does the bi-elliptic family cross Hohmann as a function of
  r1/r2 and the intermediate radius? (005)
- Combined plane-change + transfer minimum (cock from the two-impulse cost
  surface computed here) (006)
- Pork-chop launch-window ground truth vs mean-orbit anchor (013)