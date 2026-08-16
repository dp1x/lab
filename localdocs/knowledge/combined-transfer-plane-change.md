---
tags: [orbital-mechanics, maneuvers, transfers, plane-change, delta-v, bi-elliptic]
date: 2026-08-16
aliases: [combined-transfer-plane-change, super-synchronous-maneuver]
links:
  - "[[bielliptic-vs-hohmann]]"
  - "[[hohmann-transfer]]"
  - "[[kepler-orbit-validation]]"
---

# Combined Bi-elliptic Transfer + Plane Change: the Global Optimum Has Three Regimes

## Summary

The global minimum delta-v for a simultaneous radius change (r1 -> r2, ratio
R = r2/r1) and inclination change (delta_i) between two circular orbits has
**three distinct regimes**, not two. Below a boundary delta_i = di_c(R) the
ordinary **two-burn** combined Hohmann + plane change wins. Between di_c(R)
and a second boundary di_inf(R) a **finite intermediate-apoapsis three-burn**
(super-synchronous) dip beats it. Above di_inf(R) the optimal three-burn is the
**s -> infinity** (bi-parabolic) limit, where the entire plane change becomes
free at the near-rest apoapsis and the cost is exactly the coplanar
bi-parabolic value (sqrt(2)-1)(1 + 1/sqrt(R)), independent of delta_i.

The finite-s window **pinches shut at R_pinch ~ 6.43**: for R > 6.43 the only
three-burn regime that beats two-burn is the s->infinity one — no finite
intermediate apoapsis is ever optimal. This resolves the prior investigation's
"abrupt behavior near R ~ 6.41".

## Content

### Cost expressions (units of v1 = sqrt(mu/r1), r1 = 1, r2 = R)

Two-burn (split theta1 + (delta_i - theta1)):

    dv1 = sqrt(1 + vp^2 - 2 vp cos theta1)
    dv2 = sqrt(v2^2 + v_apo^2 - 2 v2 v_apo cos(delta_i - theta1))
    vp = sqrt(2R/(1+R)),  v2 = 1/sqrt(R),  v_apo = sqrt(2/(R(1+R)))

Three-burn via apoapsis s (split theta1 + theta2 + theta3 = delta_i):

    dv1 = sqrt(1 + vp12^2 - 2 vp12 cos theta1)
    dv2 = sqrt(va1^2 + va2^2 - 2 va1 va2 cos theta2)
    dv3 = sqrt(v2^2 + vp23^2 - 2 v2 vp23 cos theta3)
    vp12 = sqrt(2s/(1+s)),  va1 = sqrt(2/(s(1+s))),
    va2  = sqrt(2R/(s(R+s))),  vp23 = sqrt(2s/(R(R+s)))

s -> infinity limit: f_bp(R) = (sqrt(2)-1)(1 + 1/sqrt(R)) — identical to the
coplanar bi-parabolic limit of Exp 005 (verified to ~1e-16).

### Regime boundaries (selected, deg)

| R     | di_c   | di_inf  | finite-s window |
|-------|--------|---------|-----------------|
| 1.05  | 17.01  | 60.17   | 43.15           |
| 2.00  | 37.88  | 57.35   | 19.47           |
| 4.00  | 41.85  | 48.53   | 6.68            |
| 6.21* | 38.74  | 38.74   | 0.00 (pinch)    |
| 8.00  | 31.87  | 31.87   | 0               |
| 12.00 | 0      | 0       | 0               |

> **Boundary-resolution correction (audit 2026-08-16).** The earlier
> published table (di_c(1.05) = 11.24 deg, pinch at R ~ 6.427) was wrong at
> the float-noise tie. di_c is now defined at a robust 1e-5 win margin
> (above the ~1e-7 grid-optimization noise), giving di_c(1.05) = 17.01 deg
> (high-precision mpmath confirms no 3-burn win 11-15 deg, clear win at 18
> deg) and pinch at R ~ 6.21. The three-regime structure, the s->infinity
> identity, the R=2/47.5 deg dip (1.77%), and the SES-8 5.21% anchor are
> unchanged.

(*) Pinch at R ~ 6.427: window width -> 0. di_inf(R) decreases monotonically;
di_c(R) rises to ~40.9 deg near R = 3.8 then falls to 0 as R grows.

### Methodological lesson

The three-burn cost is **non-unimodal in s**, so a local (golden-section) s
optimizer is unreliable — it can miss the finite dip or lock onto s->infinity.
A dense global (s, theta1, theta2) grid is required, with the analytic
s->infinity limit folded in as a candidate. The two-burn split, by contrast,
is unimodal in theta1 and golden-section-safe.

### R = 1 detour anchors

Pure inclination change on one orbit: direct burn optimal for delta_i <=
2 arcsin(1/3) ~ 38.9424 deg; 3-burn detour cheaper for 38.9424 < delta_i < 60;
optimum s* -> infinity at delta_i >= 60 deg, cost -> 2(sqrt(2)-1).

### Real anchors (km/s)

- LEO(200km)->GEO, 28.6 deg and Curtis 300km->GEO, 28.6 deg: two_burn optimal
  (4.27 / 4.27 km/s; splitting the plane change ~2.2 deg at LEO saves ~0.6%
  over all-at-apogee).
- Super-synchronous 90 000 km apogee, 30 deg (SES-8 / Thaicom-6 class):
  **infinite_s regime wins by 5.21%** (4.096 vs 4.321 km/s) — the only anchor
  where a 3-burn beats the combined Hohmann + plane change.

## Source Experiments

- `research/orbital-mechanics/experiments/planeChangeManeuvers/` — the global
  optimum, regime boundaries, pinch point, RK4 validation, mpmath cross-checks.
- `research/orbital-mechanics/experiments/biellipticVsHohmann/` (Exp 005) —
  supplies the coplanar bi-parabolic limit f_bp(R) used as the s->infinity
  asymptote (identity verified to 1e-16).
- `research/orbital-mechanics/experiments/hohmannTransfer/` (Exp 004) — the
  two-burn Hohmann reference and IAU constants.

## Key Takeaways

- The combined transfer + plane-change problem has THREE regimes, with two
  boundaries di_c(R) and di_inf(R); the prior "two regimes" framing (all-at-r_b
  vs split) is insufficient — a finite-s intermediate optimum genuinely beats
  two-burn for moderate R and mid-range delta_i.
- The s->infinity (super-synchronous) strategy is the winning 3-burn for large
  delta_i at all R, and for ALL 3-burn-winning cases once R > 6.43.
- Non-unimodality in s is the key numerical trap; global grid search + analytic
  limit is the correct solver. RK4 trajectory validation (independent 3D
  Cowell) reproduces the closed-form burns and final inclined circular orbit to
  1e-11..1e-14.

## See Also

- [[bielliptic-vs-hohmann]]
- [[hohmann-transfer]]
- [[kepler-orbit-validation]]
