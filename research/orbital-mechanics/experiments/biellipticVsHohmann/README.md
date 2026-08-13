# Experiment Card: Bi-elliptic vs Hohmann Transfer Crossover

> Status: complete
> Date: 2026-08-13
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/biellipticVsHohmann/`

## Research Question

For a two-impulse Hohmann transfer between circular orbits `r1 -> r2` (R = r2/r1),
when does the three-burn bi-elliptic transfer (via an intermediate radius `r_b`, s =
r_b/r1) cost less in total delta-v, and what are the exact boundaries of the
crossover region? Concretely: (1) locate and verify the classical ratios R_bp =
11.94 (bi-parabolic crossover) and R* = 15.58 (always-cheaper boundary); (2) prove
numerically that R* coincides with the Hohmann cost maximum; (3) map the crossover
curve s_c(R) and the saving/time trade-off; (4) validate the closed-form cost model
against an RK4 three-burn trajectory and against published worked examples.

## Background Theory

Two circular orbits (r1 = 1, v1 = sqrt(mu/r1) units; r2 = R > 1). Hohmann:

    dv_H(R) = sqrt(2R/(1+R)) - 1 + (1/sqrt(R))(1 - sqrt(2/(1+R)))

Bi-elliptic via intermediate radius s > 1 (high family), three burns:

    f_high(R, s) = [sqrt(2s/(1+s)) - 1]
                   + (1/sqrt(s))[sqrt(2R/(R+s)) - sqrt(2/(1+s))]
                   + (1/sqrt(R))[1 - sqrt(2s/(R+s))]

Each term is the vis-viva speed difference at a burn point (departure r1, deep-space
s, arrival r2). Bi-parabolic limit s -> infinity:

    f_bp(R) = (sqrt(2)-1)(1 + 1/sqrt(R))

with leading correction (sqrt(2)/2)(sqrt(R)-3)/s (verified numerically: the limit is
approached from below for R < 9, from above for R > 9).

Classical results (Hoelker & Silber 1959; Vallado; Wikipedia): the bi-elliptic wins
for some s iff R > 11.94 (R_bp, where f_bp = dv_H), and wins for ALL s iff R > 15.58
(R*). At R* the crossover curve s_c(R) meets the corner s = R. The corner identity

    d/ds f_high(R, s)|_{s=R} = d/dR dv_H(R)

makes R* coincide with the maximum of the Hohmann cost curve (Exp 004: peak at
15.5817187369).

Low family (r_b < r1, inward excursion) is never cheaper: its three burns are all
positive costs (the deep-space burn at r_b is a boost), so it cannot beat the
two-impulse Hohmann optimum.

## References

- Hoelker, R. F. & Silber, R., *The Bi-Elliptical Transfer Between Circular
  Coplanar Orbits*, Army Ballistic Missile Agency, 1959.
- Vallado, D. A., *Fundamentals of Astrodynamics and Applications*, 4th ed., 2013.
- Wikipedia, "Bi-elliptic transfer" (worked example, R = 14, s = 40 with
  r0 = 6700 km, mu = 398600.4418 km^3/s^2), retrieved 2026.
- Exp 004 card: `research/orbital-mechanics/experiments/hohmannTransfer/`.

## Assumptions

- Circular coplanar starting/target orbits; impulsive burns; two-body problem with
  point masses (no oblateness, atmosphere, third bodies, finite burn times).
  [idealization]
- All radii normalized to r1 and speeds to the circular speed at r1 (model scale-free
  in the two-body problem). [verified]
- RK4 integrator step for trajectory validation only; closed-form results do not
  depend on it. [verified]
- Earth constants (r_E = 6378.1 km, mu = 398600.4418 km^3/s^2, LEO 200 km, GEO
  35786 km) for the real-anchor section. [published values]

## Methodology

- Closed-form normalized costs f_high(R, s) (high family), f_low (low family),
  bi-parabolic limit f_bp, flight times (half-period sums), all vectorized.
- R_bp solved from f_bp(R) = dv_H(R); R* solved from the corner identity
  d/ds f_high = d/dR dv_H; each cross-checked with a 50-digit mpmath recomputation
  (finite-difference h = 1e-20, giving ~1e-29 agreement).
- Crossover curve s_c(R) for R in (R_bp, R*) by bisection on the signed difference
  g = f_high - dv_H, 61 R-values.
- Adversarial region verification: 90 R-values per region x 400 s-values
  (log-spaced s/R in [1.000001, 1e6]) asserting the sign of g everywhere.
- Shape diagnostics classify g(s) per R as monotone increasing / single-hump /
  monotone decreasing (401-point s-grid) to justify the unique-crossing claim.
- RK4 validation of the full three-burn trajectory (Exp 002 propagator): burn
  placement, apoapsis/periapsis radii, circularization, energy/angular-momentum
  drift, for six (R, s) cases.
- Real-anchor budgets in km/s and days: LEO-GEO, GEO-lunar, GEO-15.58x, Earth-Mars,
  the Wikipedia 14x example (r0 = 6700 km), and a near-bi-parabolic 50x case.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib, mpmath
- Runtime: `uv run python experiment.py` (writes `results/results.json` + 6 figures)
- Determinism: pure float64 + fixed mpmath precision, no RNG; repeated runs are
  byte-identical apart from the timestamp (verified).
- Dependencies: numpy, matplotlib, mpmath (added via `uv add`, locked in uv.lock).

## Validation Method

- `tests/test_bielliptic_vs_hohmann.py` (44 tests): closed-form burn identities,
  corner = Hohmann degeneracy, bi-parabolic limit with leading correction,
  time-reversal symmetry (burns swap under reversal), time-penalty monotonicity,
  corner identity (finite differences), mpmath 50-digit tables, region sign,
  crossover-curve structure, shape classes, RK4 transfer validation, real anchors.
- All 175 lab tests pass (`uv run pytest`).
- External anchors: Wikipedia's worked R = 14, s = 40 example reproduced exactly
  (dv = 3.061043 + 0.608825 + 0.447662 = 4.117530 km/s vs Hohmann 4.133716 km/s);
  literature boundary ratios 11.94 / 15.58 reproduced to 10 decimals.

## Results

- R_bp = 11.9387654726 (literature 11.94); R* = 15.5817187388 (literature 15.58)
  = Hohmann cost maximum of Exp 004 (15.5817187369, cross-check 1e-10 rel).
- Corner identity: relative agreement 4.3e-9 by float64 finite differences;
  1.08e-29 at 50-digit precision.
- Crossover curve: s_c(R_bp) = infinity (diverges), s_c(12) = 815.82,
  s_c(13) = 48.90, s_c(14) = 26.10, s_c(15) = 18.19 (float64 vs mpmath agree to
  9+ digits), s_c(R*) = R* (15.5817): the curve spans a factor 5.4e4 in s.
- Region verification: g > 0 everywhere for R < R_bp (worst margin 1.72e-8);
  g < 0 everywhere for R > R* (worst margin 1.6e-2); exactly one crossing per R
  for R_bp < R < R* (61/61).
- Shape classes: monotone increasing for R <= 8.93; single hump for
  9.53 <= R <= 15.00; monotone decreasing for R >= 16.00. The hump onset
  (R ~ 9.53) lies BELOW R_bp: the hump exists for intermediate R but never dips
  under dv_H until R > R_bp - this resolves why the literature threshold is
  11.94, not the hump onset.
- Max saving: 4.09% of v1 at R = 50.1 (bi-parabolic limit, s -> infinity); time
  ratio t_biell/t_H reaches ~3.7 (R=2, s=3) to ~24 (R=20, s=100), unbounded as
  s grows.
- RK4: apoapsis/periapsis radii rel error ~1e-8..1e-10; burn deltas-v rel error
  ~1e-7..1e-9; energy drift ~1e-9; circular-orbit radius variation ~4.5e-11.
- Real anchors (km/s): LEO-GEO (R=6.41): 4.4709 vs 3.9319, -13.7% (Hohmann wins);
  GEO-lunar (R=11.47, s=27): -4.1%; Earth-Mars (R=4.6): -173% (Hohmann wins);
  GEO-15.58x (s=30): +0.58%; Wikipedia 14x (s=40): +0.39% saving, t_ratio 11.3;
  near-bi-parabolic 50x (s=1e6): +7.96% saving, t_ratio 5.5e6.
- Low family (r_b < 1): never cheaper than Hohmann for any R > 1, s < 1
  (all three burns are strictly positive costs).

Figures (`results/figures/`): cost_curves, crossover_map, saving_curve,
shape_per_R, time_penalty, trajectory_geometry.

## Limitations

- Two-body impulsive model: no J2, drag, gravity assists, finite burn arcs, or
  phasing; results are idealized delta-v comparisons, not mission delta-v.
- "Always cheaper" (R > R*) is proven numerically on a 90x400 grid, not
  analytically; the corner-identity + shape structure makes it compelling but the
  formal proof is outside the numerical experiment.
- The crossover curve near R_bp is ill-conditioned (s_c diverges); the recorded
  endpoint is symbolic (infinity) and float64/mpmath tables agree to 9+ digits
  elsewhere.
- Flight times assume half-ellipse arcs with no waiting/phasing constraints.

## Future Improvements

- Derivations/analytic proofs of the shape classes and the leading asymptote
  (sqrt(2)/2)(sqrt(R)-3)/s are documented but not formalized.
- Non-coplanar / non-circular generalizations (elliptic-to-elliptic bi-elliptic).
- Finite-burn (low-thrust) extension is a separate experiment.

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions (numpy, matplotlib, mpmath).
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
