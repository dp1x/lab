# Experiment Card: Combined Bi-elliptic Transfer + Plane Change (global optimum)

> Status: complete
> Date: 2026-08-16
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/planeChangeManeuvers/`
> Note: this is the "006 — Plane-change maneuvers" slot, executed here in
> place (the active work from the previous agent). Experiment 005 (coplanar
> bi-elliptic vs Hohmann) is a separate completed experiment; the combined
> transfer + plane-change question continues it and is closed by this card.

## Research Question

For a change of both **orbital radius** (r1 -> r2, radius ratio R = r2/r1) and
**inclination** (delta_i) between two circular orbits, what is the **global
minimum delta-v** when the maneuver may use any of:

1. **two burns** — the combined Hohmann + plane change (a 2-impulse transfer
   with the plane change split between the two burns);
2. **three burns** — raise apoapsis to an intermediate radius s, do the plane
   change (split across all three burns), then lower periapsis to r2 (the
   super-synchronous / bi-elliptic-with-plane-change strategy);
3. the **s -> infinity limit** of (2), where the apoapsis velocity vanishes and
   the entire plane change becomes "free" — cost independent of delta_i.

Does the true global optimum have **distinct regimes** (ordinary two-burn;
finite intermediate-apoapsis three-burn; asymptotic s -> infinity) and where
(in the (R, delta_i) plane) do the **boundaries between them** occur?

## Background Theory

Two circular orbits (r1 = 1, v1 = sqrt(mu/r1) units; r2 = R > 1; an inclination
change delta_i about the common node line). All costs are in units of v1.

**Two-burn combined transfer** (Part B). Burn 1 at r1 puts the spacecraft on
the Hohmann ellipse and rotates the velocity by theta1 about the node axis;
burn 2 at r2 circularizes and tilts by the remaining (delta_i - theta1):

    dv1 = sqrt(1 + vp^2 - 2 vp cos theta1)
    dv2 = sqrt(v2^2 + v_apo^2 - 2 v2 v_apo cos(delta_i - theta1))
    vp = sqrt(2R/(1+R)),  v2 = 1/sqrt(R),  v_apo = sqrt(2/(R(1+R)))

The cost is unimodal in the split theta1 in [0, delta_i], so golden section
finds the optimum reliably.

**Three-burn super-synchronous** (Part C). With an intermediate apoapsis
s >= R (r_b = s r1), the three impulse magnitudes are single law-of-cosines
burns with the split theta1 + theta2 + theta3 = delta_i distributed across the
departure, the deep-space apoapsis, and the arrival:

    dv1 = sqrt(1 + vp12^2 - 2 vp12 cos theta1)
    dv2 = sqrt(va1^2 + va2^2 - 2 va1 va2 cos theta2)
    dv3 = sqrt(v2^2 + vp23^2 - 2 v2 vp23 cos theta3)
    vp12 = sqrt(2s/(1+s)),  va1 = sqrt(2/(s(1+s))),
    va2  = sqrt(2R/(s(R+s))),  vp23 = sqrt(2s/(R(R+s)))

The cost surface over (s, theta1, theta2) is **NOT unimodal in s** — a local
optimizer (golden section in s) is therefore unreliable, which is why a dense
global grid is used.

**s -> infinity limit** (the "free plane change at infinity"). As s -> inf,
vp12 -> sqrt(2)-1, va1, va2, vp23 -> 0, and the deep-space burn vanishes; the
cost tends to the **coplanar bi-parabolic limit of Experiment 005**:

    f_bp(R) = (sqrt(2)-1)(1 + 1/sqrt(R))       [independent of delta_i]

This is a genuine limiting regime (verified to ~1e-16 against Exp 005).

**R = 1 detour anchors** (pure inclination change on one orbit). The classic
result: a single burn is optimal for delta_i <= 2 arcsin(1/3) ~ 38.9424 deg;
for 38.9424 < delta_i < 60 deg a 3-burn detour (raise apoapsis, plane-change at
the slow apoapsis, lower back) is cheaper; at delta_i >= 60 deg the optimum
diverges to infinity and the cost tends to 2(sqrt(2)-1) (plane change free at
the vanishing apoapsis speed).

## References

- Curtis, H. D., *Orbital Mechanics for Engineering Students*, 4th ed.,
  Elsevier 2021 — combined Hohmann + plane change (law of cosines); the
  300 km LEO, 28.6 deg -> GEO worked example (optimal split ~2.5 deg at LEO).
- Gonzalez, *Orbital Mechanics & Astrodynamics* (orbital-mechanics.space) —
  Plane-Change Maneuver example (same worked case).
- Wikipedia, "Geostationary transfer orbit" (law of cosines with cos(delta_i)
  at apogee) and "Supleynchronous orbit" (SES-8, Thaicom-6 apogee 90 000 km).
- Wakker, B., *Optimal Impulsive Orbit Transfers*, Springer 2015.
- Hoelker & Silber 1959; Exp 004 (Hohmann) and Exp 005 (coplanar bi-elliptic
  crossover / bi-parabolic limit) for the transfer machinery reused here.

## Assumptions

- Circular initial/target orbits; impulsive burns; two-body point-mass gravity
  (no oblateness, drag, finite burn arcs, phasing). [idealization, standard]
- Radius change and inclination change share the **common node line**; only the
  magnitude delta_i matters (combined inc + RAAN reduces to the total angle
  between plane normals). [verified]
- All radii normalized to r1 and speeds to v1; the two-body problem is
  scale-free so physical (km, km/s) values are v1 = sqrt(mu/r1) times the
  normalized numbers. [verified]
- The s -> infinity limit is taken as the analytic bi-parabolic value
  (delta_i-independent); the global optimizer compares finite dips against it.

## Methodology

- **Optimizer (non-unimodal)**: for each (R, delta_i) the 3-burn cost is
  minimized over a **dense log-spaced s grid (240 points, s in [1.0001 R, 1e6])**
  x a **2D plane-change split grid (64 x 64 in (theta1, theta2))** with two
  narrowing refinement passes and a final nested golden-section polish of the
  split. The global grid minimum is kept (NOT a local optimizer). The analytic
  s -> infinity limit is folded in so the true global minimum is chosen among
  {two-burn, finite-s 3-burn, s->infinity}.
- **Boundary mapping**: di_c(R) = smallest delta_i where 3-burn beats two-burn;
  di_inf(R) = delta_i above which the optimal 3-burn is the s->infinity regime
  (finite-s window upper edge). Both via bisection (tol 1e-4 deg) on the
  global optimizer, over R in ~[1.05, 50] (74 points, log-spaced plus a dense
  cluster 5.8-7.2 to resolve the pinch).
- **Independent validation**: (a) a second brute-force optimizer (ns=400,
  nth=96, s_max=1e8) cross-checked at 5 representative points; (b) a full 3D
  RK4 (Cowell) propagation of the optimal two-burn and finite-s three-burn
  maneuvers, checking burn magnitudes, arrival radius, circular speed, angular
  momentum, and the rotated orbital normal; (c) 50-digit mpmath recomputation
  of the s->infinity identity, the R=1 detour corners, and the R=2, delta_i=47.5
  finite-s dip; (d) real-system anchors (LEO->GEO 28.6 deg, Curtis 300 km /
  28.6 deg, GTO->GEO 5 deg, SES-8 super-synchronous 90 000 km / 30 deg).
- **Search-boundary check**: the s->infinity (s_max = 1e6) grid is confirmed
  against an enlarged s_max = 1e8 optimizer (no missed finite dips).
- **Determinism**: pure float64 + fixed mpmath precision (dps=50), no RNG.
  Repeated runs are byte-identical apart from the timestamp.

## Implementation

- Script: `experiment.py` (Parts A-E: pure plane change, two-burn, three-burn
  global optimizer, regime boundaries, RK4 validation, real anchors, mpmath
  cross-checks, figures).
- Language/runtime: Python 3.12, numpy 2.5.1, matplotlib 3.11.1, mpmath 1.4.1.
- Runtime: `uv run python experiment.py` (writes `results/results.json` + 3
  figures). Compute: ~8-10 min on a single core for the boundary + sweep grids.
- Determinism: documented above; the optimizer has no stochastic state.
- Dependencies: numpy, matplotlib, mpmath (already in `uv.lock` from Exp 005).

## Validation Method

- `tests/test_plane_change.py` (14 tests): closed-form burn identities, the
  s->infinity = Exp 005 bi-parabolic identity (to 1e-15), the split law
  theta1 + theta2 + theta3 = delta_i, the R = 1 detour corners, regime code
  paths, the finite-s window pinch near R ~ 6.4, alternate-optimizer agreement
  (non-unimodal check), 3D RK4 trajectory validation of two-burn and three-burn
  maneuvers, and optimizer determinism.
- Full lab suite (`uv run pytest`) passes (001-005 + this experiment).
- RK4 trajectory checks: burn dv rel error ~1e-11..1e-14, arrival radius error
  ~1e-11..1e-13, circular speed error ~1e-14, angular-momentum error
  ~1e-11..1e-13, normal-alignment cos = 1.0 to 1e-6.
- mpmath 50-digit: s->inf identity diff ~1.6e-16 at R = 2; the R = 2,
  delta_i = 47.5 finite-s dip beats two-burn by 1.765% (mpmath-agreed).

## Results

**Three distinct regimes exist**, separated by two boundaries di_c(R) and
di_inf(R):

- **two_burn** for delta_i < di_c(R).
- **finite_s** (finite intermediate apoapsis s* > R) for di_c(R) < delta_i <
  di_inf(R) — a genuine non-unimodal dip that beats two-burn.
- **infinite_s** (s -> infinity, plane change free at the near-rest apoapsis)
  for delta_i > di_inf(R).

Boundary highlights (R, di_c in deg, di_inf in deg):

| R     | di_c   | di_inf  | finite-s window (deg) |
|-------|--------|---------|-----------------------|
| 1.05  | 11.24  | 60.17   | 48.93                 |
| 2.00  | 36.17  | 57.35   | 21.18                 |
| 4.00  | 40.81  | 47.89   | 7.08                  |
| 6.41  | 37.88  | 38.42   | 0.54                  |
| 6.43* | 37.85  | 37.85   | 0.00 (pinch)          |
| 8.00  | ~31    | ~31     | 0 (no finite-s window)|
| 12.00 | 8.75   | 8.75    | 0                     |
| 20.00 | 0      | 0       | 0                     |

(*) **The finite-s window pinches shut at R_pinch ~ 6.43** (window width
-> 0 between R = 6.42 and R = 6.43). This is the "abrupt behavior near
R ~ 6.41" the prior investigation flagged: for R > 6.43 the only 3-burn
regime that beats two-burn is the s->infinity one — no finite intermediate
apoapsis is ever optimal. For small R the finite-s window is wide (21 deg at
R = 2); di_inf(R) decreases monotonically and di_c(R) first rises (to ~40.9
deg near R = 3.8) then falls toward 0 as R grows.

**Canonical finite-s dip**: R = 2, delta_i = 47.5 deg -> finite s* = 2.72,
dv = 0.6501 (v1 units) vs two-burn 0.6618 — a **1.77% saving** (the prior
"1-2%" claim, confirmed at 50-digit precision).

**s->infinity identity**: bi_parabolic_plane_change_limit(R) equals the
Exp 005 coplanar bi-parabolic limit to ~1e-16 — the combined problem's
asymptote is exactly the coplanar bi-parabolic cost, independent of delta_i.

**R = 1 detour**: di_c = 38.9424 deg (2 arcsin(1/3)), di_inf = 60 deg;
optimum s* -> infinity at delta_i >= 60 deg, cost -> 2(sqrt(2)-1).

**Real anchors** (km/s, v1 = sqrt(mu/r1)):

| Case | regime | best dv | two-burn dv | saving |
|------|--------|---------|-------------|--------|
| LEO(200km)->GEO, 28.6 deg | two_burn | 4.269 | 4.269 | 0% (Hohmann+split optimal) |
| Curtis 300km->GEO, 28.6 deg | two_burn | 4.266 | 4.266 | 0% (theta1=2.2 deg at LEO) |
| GTO->GEO, 5 deg | two_burn | 3.943 | 3.943 | 0% |
| Super-synchronous 90 000 km, 30 deg | infinite_s | 4.096 | 4.321 | **5.21%** |

The super-synchronous SES-8 / Thaicom-6 class is the real regime where the
3-burn (s->infinity) strategy wins by several percent — the only anchor here
where a 3-burn beats the combined Hohmann + plane change.

Figures (`results/figures/`): `regime_map.png` (the (R, delta_i) regime
diagram with both boundaries), `cost_vs_di.png` (cost vs delta_i at R = 1.5,
2, 6.41 with the three candidate curves), `s_star_vs_di.png` (finite-s s*
vs delta_i, showing the dip and divergence).

## Limitations

- Two-body impulsive model: no J2, drag, gravity assists, finite-burn arcs, or
  phasing constraints (the super-synchronous 5% saving is a delta-v optimum,
  not a mission delta-v; the time penalty is unbounded as s -> infinity).
- The regime boundaries are located numerically (bisection, tol 1e-4 deg) over
  a discrete R grid; the pinch R ~ 6.43 is bracketed to ~1e-3 in R, not proved
  analytically.
- The 3D RK4 validation uses a custom Cowell integrator (independent of Exp
  002's planar propagator); both are first-principles but the 3D one is new
  code added here and validated only against the closed-form burns it was
  built to reproduce.
- Node line assumed common; combined inclination + RAAN with a different node
  geometry is out of scope (reduces to the total normal angle anyway).

## Future Improvements

- Analytic proof of the window-pinch R and the di_c(R)/di_inf(R) curves.
- Elliptic-to-elliptic or non-coplanar bi-elliptic generalizations.
- Finite-burn (low-thrust) extension and a time-penalty trade-off surface
  (the s->infinity savings are unbounded in flight time).
- Direct comparison against a patched-conic / optimization tool (poliastro,
  GMAT) for the super-synchronous anchors.

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions (numpy 2.5.1, matplotlib 3.11.1,
  mpmath 1.4.1, Python 3.12).
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
- Determinism verified: two optimizer runs return bit-identical results (no RNG).
- Optimal (R, s) split satisfies theta1 + theta2 + theta3 = delta_i to 1e-6.
