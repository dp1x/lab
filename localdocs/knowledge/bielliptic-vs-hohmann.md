---
tags: [orbital-mechanics, manuevers, transfers, bi-elliptic, hohmann, delta-v]
date: 2026-08-13
aliases: [bi-elliptic-transfer, bielliptic-vs-hohmann]
links:
  - "[[hohmann-transfer]]"
  - "[[kepler-orbit-validation]]"
---

# Bi-elliptic vs Hohmann: When Three Burns Beat Two

## Summary

Between two coplanar circular orbits (r1 = 1, v1 = sqrt(mu/r1) units, R = r2/r1 > 1)
the three-burn bi-elliptic transfer via an intermediate apoapsis s > 1 costs
less delta-v than the two-burn Hohmann transfer **iff R > 11.9387654726 and s
exceeds a crossover value s_c(R)**; for **R > 15.5817187388 it wins for every s**
(the "always cheaper" regime). Both boundaries were reproduced to 10 decimals in
Experiment 005 and tied to the Hohmann cost maximum of Experiment 004 by the
corner identity:

    d/ds f_high(R, s)|_{s=R}  =  d/dR dv_H(R)        (verified to 1e-29 at 50 digits)

The bi-parabolic limit s -> infinity costs f_bp(R) = (sqrt(2)-1)(1 + 1/sqrt(R))
with leading correction (sqrt(2)/2)(sqrt(R)-3)/s (from below for R < 9, above for
R > 9 — the "approached from below" folklore is only half right).

## Key quantitative facts (verified in Experiment 005)

- **Boundaries**: R_bp = 11.9387654726 (Hohmann wins for all s below this),
  R* = 15.5817187388 (bi-elliptic wins for all s above this). R* equals the
  Hohmann cost maximum R_star = 15.5817187369 (Exp 004) to 1e-10 relative —
  the corner identity makes the two problems coincide.
- **Crossover curve** s_c(R): diverges as R -> R_bp+ (s_c(12) = 815.82),
  falls steeply (s_c(13) = 48.90, s_c(14) = 26.10, s_c(15) = 18.19), and meets
  the corner s = R at R*. Exactly one crossing for every R in (R_bp, R*) —
  verified on a 90 x 400 sign grid (worst margin 1.7e-8 in the low region).
- **Why 11.94 and not the hump onset**: f_high(s) - dv_H develops a hump as soon
  as R ~ 9.53, but the hump stays above dv_H until R_bp; only after R_bp does it
  dip below. Shape classes: monotone increasing for R <= 8.93, single hump for
  9.53 <= R <= 15.00, monotone decreasing for R >= 16.00 (401-point s-grid).
- **Max saving**: 4.09% of v1 at R = 50.1 in the bi-parabolic limit; near-parabolic
  (s = 1e6) LEO-to-50x case saves 0.318 km/s (7.96%) but takes 5.5 million
  Hohmann times. Real budget: Wikipedia's worked example (r0 = 6700 km,
  R = 14, s = 40) reproduced exactly: 3.061043 + 0.608825 + 0.447662 =
  4.117530 km/s vs Hohmann 4.133716 km/s (saving 0.39%, t_ratio 11.3).
- **Inward low family** (r_b < r1) never wins: all three burns are strictly
  positive (the deep-space burn at r_b is a *boost*), and the s -> 1 corner
  ties Hohmann — verified across R > 1, s in (0, 1).
- **Time penalty**: t_biell/t_H = 3.7 (R = 2, s = 3) to 24.2 (R = 20, s = 100),
  unbounded as s grows; savings always cost time (see the trade-off curve).
- **Time-reversal symmetry**: reversing the transfer swaps the first and third
  burns (fwd[k] == bwd[2-k]) — both for the high family and, with the burn
  order swapped, the inward case.
- **Trajectory validation**: the full three-burn RK4 trajectory (Exp 002
  propagator) closes the target circle to radius rel error 2e-10..4e-8 and
  burn delta-v rel error 1e-9..1e-7 across six (R, s) cases; the deep-space
  burn at apoapsis is *perpendicular* to the circular velocity there (energy
  transfer without direction change at apogee).

## Connections

- Hohmann cost curve and its maximum: [[hohmann-transfer]]
- RK4 orbital propagator the trajectory validation uses: [[kepler-orbit-validation]]
- mpmath 50-digit cross-checks of the boundary ratios and the corner identity
  pin the float64 results (agreement 1e-29 at the corner).
