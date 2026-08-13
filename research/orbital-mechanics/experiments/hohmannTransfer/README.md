# Hohmann Transfer: Least-Fuel Orbit-to-Orbit Transfer and Its Delta-v Budget

> Status: complete
> Date: 2026-08-13
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/hohmannTransfer`

## Research Question

Does the Hohmann two-impulse transfer — as implemented on the verified RK4
machinery of Experiment 002 — reproduce the closed-form delta-v equations,
transfer time and orbit geometry across the full radius-ratio space
R = r2/r1 from 1+ε to infinity, is it truly the minimum-delta-v transfer
within the two-impulse class, and does it reproduce the canonical real-world
budgets (LEO→GEO ≈ 3.93 km/s; Earth→Mars ≈ 259 d with v_infinity ≈ 2.94 /
2.65 km/s; trans-Mars injection ≈ 3.6 km/s)?

## Background Theory

The Hohmann transfer (Hohmann 1925) connects two coplanar circular orbits of
radii r1 < r2 around a central body μ with two tangential impulses on a
half-ellipse tangent to both circles. The transfer ellipse has semi-major
axis a_t = (r1+r2)/2 and eccentricity e_t = (r2−r1)/(r2+r1); the flight time
is half its period,

    t_tr = π √((r1+r2)³ / (8μ))            (from Kepler III)

and the burns follow from the vis-viva equation v² = μ(2/r − 1/a):

    Δv1 = v1 (√(2R/(1+R)) − 1),   R = r2/r1,   v1 = √(μ/r1)   (at r1)
    Δv2 = v2 (1 − √(2/(1+R))),                            v2 = √(μ/r2)   (at r2)

Direction note: the textbook formulas assume r2 > r1. For inward transfers
(r2 < r1) the same ellipse is traversed in the opposite direction; the two
burns swap order but keep their magnitudes per radius
(Δv@r_small = v_small(√(2R/(1+R)) − 1), Δv@r_big = v_big(1 − √(2/(1+R))),
R = r_big/r_small), so Δv_total(r1, r2) = Δv_total(r2, r1).

Structure of the cost curve (all in units of the inner circular speed v1):

- R → 1⁺: Δv_total/v1 ~ (R−1)/2 (leading order; burns vanish).
- Interior maximum at R* ≈ 15.58, Δv_total/v1 ≈ 0.5363 — cost is NOT
  monotone in destination distance.
- R → ∞: Δv1 → the escape burn (√2 − 1)v1 and Δv2 → 0, so
  Δv_total/v1 → √2 − 1 ≈ 0.41421: reaching "infinity" costs exactly the
  escape burn in the limit. Digit-safe forms for R → 1 (no cancellation):

    √(2R/(1+R)) − 1 = (R−1) / ((1+R)(1 + √(2R/(1+R))))
    1 − √(2/(1+R))  = (R−1) / ((1+R)(1 + √(2/(1+R))))

Optimality: among two-impulse transfers with burns at r1 and r2 (transfer
ellipse with periapsis r_p ≤ r1, apoapsis r_a ≥ r2, no tangency assumed),
the total delta-v is minimized at the Hohmann ellipse (r_p = r1, r_a = r2).
(Beyond the two-impulse class, three-impulse bi-elliptic transfers win for
R ≳ 11.94 — roadmap experiment 005.)

## References

- W. Hohmann, *Die Erreichbarkeit der Himmelskörper*, Oldenbourg, 1925 —
  the minimum-energy transfer between circular orbits.
- R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*,
  Dover, 1971, Ch. 6 — orbital maneuvers, Hohmann transfer equations.
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed.,
  Elsevier, 2021, Ch. 6 — Hohmann transfer, trans-Mars injection worked
  values (LEO→GEO ≈ 3.93 km/s family of examples).
- D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed.,
  Microcosm, 2013, Ch. 6 — transfer geometries, R* ≈ 15.58 peak discussions.
- IAU 2015 Resolution B3 nominal constants (Mamajek et al., arXiv:1510.07674):
  (GM)^N_Sun = 1.3271244e11 km³/s², (GM)^N_E = 3.986004e5 km³/s²,
  R^N_eE = 6378.1 km; IAU 2012 Resolution B2 au = 149597870.7 km (verified
  against the iau.org resolution text, 2026-08-13).
- JPL SSD, "Approximate Positions of the Planets" (Standish & Williams 1992)
  — Mars a = 1.523679 au, Venus a = 0.723332 au (mean ecliptic J2000,
  EM barycenter; same source family as Exp 002).
- Canonical public anchors verified 2026-08-13: Earth→Mars Hohmann ≈ 259 d
  with departure/arrival excess 2.94/2.64–2.65 km/s and total ≈ 5.6 km/s
  (NASA-adjacent and textbook-aligned public sources, e.g. marspedia.org),
  trans-Mars injection from LEO ≈ 3.6 km/s ("only about 0.4 km/s more than
  escape", Wikipedia), LEO→GEO ≈ 3.9 km/s (public engineering sources).

## Assumptions

- Two circular, coplanar, same-body orbits; impulsive burns (zero burn
  time); two-body problem (no third bodies, no oblateness) — the textbook
  idealization (idealization, standard for Hohmann analysis).
- Canonical units (μ = 1) for the core sweep; real units (IAU nominal
  constants, JPL mean orbits) for the anchors (verified).
- The real-system anchors idealize planetary orbits as circular mean orbits
  (Mars e ≈ 0.093 ignored); real launch windows add eccentricity and
  inclination corrections (plausible; standard first-order modeling).
- Tsiolkovsky fuel fractions assume a single rocket stage with constant
  Isp, no gravity losses (illustration, not a mission design).
- Floating-point evaluation of Δv formulas is digit-safe near R = 1 via
  rearranged forms (verified: masked cancellation).
- The 002 propagator (RK4, periapsis-resolved step counts) is reused as the
  verified integration layer; step counts follow the (1−e)^(−3/2) law of
  Experiment 002 (verified in this experiment's arrival errors).

## Methodology

All runs in `experiment.py`, deterministic (no RNG):

1. **Closed forms from first principles** — Δv1, Δv2 via vis-viva; transfer
   time via Kepler III; h conservation r_small·v(r_small) = r_big·v(r_big);
   energy −μ/(2a_t); the transfer ellipse's elements recovered by Exp 002's
   machinery from the burned state.
2. **Fully covers the ratio space** — R ∈ {1.000001 … 1e12} (25 cells,
   log-spaced), reporting Δv1, Δv2, Δv_total/v1, transfer time, transfer
   eccentricity; asymptotic checks (R−1)/2 at R → 1 and √2−1 at R → ∞;
   interior peak located by dense grid + parabolic refinement
   (R* = 15.5817, value 0.536258); digit-safe vs textbook forms.
3. **Complete transfer propagated with RK4** (Exp 002 machinery) for
   R ∈ {1.5, 6.41, 20} plus a Venus-like ratio case (R = 1.3825) and a true
   inward transfer (R = 0.5, r2 < r1): burn 1 at r1 → half-orbit coast →
   verify arrival radius, speed, radial velocity and apside timing → burn 2
   → one full orbit of the target circle verified. The analytic same-ellipse
   reference is phase-shifted by a half-period for inward flights (which
   start at apoapsis), so it always compares against the arrival apside.
4. **Optimality scan (two-impulse family)** at R ∈ {2, 6.41, 20}: cost
   Δv1 + Δv2 with vector velocity mismatches (no tangency assumed) over a
   121 × 131 grid of transfer ellipses (r_p ≤ r1, r_a ≥ r2), plus 1D
   families (tangent departure r_p = r1; fixed apoapsis r_a = r2).
5. **Real anchors**: LEO(200 km)→GEO with IAU nominal Earth; Earth→Mars and
   Earth→Venus heliocentric with IAU Sun + JPL mean orbits; trans-Mars and
   trans-Venus injection from 200 km LEO (hyperbolic excess via
   √(2v_circ² + v_inf²)); Tsiolkovsky propellant fractions at Isp 300/450 s.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib (no scipy)
- Runtime: `uv run python experiment.py`
- Determinism: no RNG; pure functions; subprocess bit-identical check in the
  test suite.
- Dependencies: numpy, matplotlib (already in pyproject; no new deps).
- Reuses `lab_utils.results` (save_json_result) and — deliberately — the
  verified propagator/Kepler machinery of Experiment 002
  (`keplerOrbitValidation/experiment.py`, loaded by explicit-path importlib):
  single source of truth for propagation, no scaffolding rebuilt.
- Resource footprint: seconds of CPU, < 1 MB outputs — no scratch needed.

## Validation Method

Independent checks in `tests/test_hohmann_transfer.py` (36 tests, all green):

- Closed forms re-derived in-test: vis-viva Δv1/Δv2, half-period flight
  time, h-conservation across the ellipse, energy = −μ/(2a_t), ellipse
  elements recovered by 002's machinery; degenerate R = 1 gives zero cost;
  Δv1 < escape burn everywhere.
- Direction symmetry: Δv_total(r1, r2) = Δv_total(r2, r1) and
  burn magnitudes swap per radius (|dv1|in = dv2out) to 1e-12.
- Asymptotes: (R−1)/2 at R = 1.001/1.0001 within 1e-3; √2−1 at R = 1e12
  within 1e-4 relative; peak R* ∈ (15.5, 15.65), value ∈ (0.535, 0.538);
  Δv2 vanishes at both R → 1 and R → ∞ (interior max 0.190 at R ≈ 5.88).
- RK4 trajectory: arrival radius/speed errors < 1e-6 (R ≤ 6.41) and < 1e-4
  (R = 20) — measured 4.5e-10…4.2e-9; analytic(same ellipse) arrival exact
  to < 1e-9; apside reached exactly at t_tr; burn 2 vector difference
  matches Δv2 to < 1e-4 (measured ≤ 5e-9); post-burn orbit circular to
  < 1e-6 (measured 4.5e-11); energy/h drift < 1e-6.
- Optimality: the 2D grid minimum coincides with the Hohmann corner
  (gap measured 0 to 7.8e-16); 1D families min at r_a = r2 / r_p = r1.
- Real anchors in published bands: LEO→GEO total 3.93 ± 0.05 km/s, coast
  5.3 ± 0.3 h; Earth→Mars 259 ± 1 d, v∞ 2.94 ± 0.03 / 2.65 ± 0.05 km/s,
  TMI 3.6 ± 0.1 km/s; Venus inward ~146 d.
- Determinism: bit-identical JSON in a fresh interpreter (subprocess).

An independent numeric cross-check (this session, mirrored in tests): the
burn-magnitude swap identity, the R→1 digit-safe forms, h-conservation and
the hand-derived escape-burn asymptote were re-derived by hand algebra and
verified numerically against the code on multiple grids.

## Results

All results in `results/results.json`; figures in `results/figures/`.

**Closed-form budget** (r1 = 1, r2 = 2, μ = 1): Δv1 = 0.15470, Δv2 = 0.12976,
total 0.284457 = 0.284457·v1; t_tr = 5.77147 = half of 2π√(a_t³/μ) =
11.54295.

**Cost curve across R** (Δv_total/v1):

| R | Δv1/v1 | Δv2/v1 | total/v1 | note |
|---|--------|--------|----------|------|
| 1.0001 | 2.5000e-5 | 2.4996e-5 | 4.99963e-5 | ≈ (R−1)/2 (0.9999× asymptote) |
| 1.5 | 0.09545 | 0.08620 | 0.181645 | |
| 2.0 | 0.15470 | 0.12976 | 0.284457 | |
| 6.41 | 0.31533 | 0.18978 | 0.505103 | LEO→GEO-like ratio |
| 11.94 | 0.35847 | 0.17562 | 0.534095 | first bi-elliptic crossover |
| 15.58 | 0.37090 | 0.16536 | **0.536258** | interior max |
| 20.0 | 0.38013 | 0.15460 | 0.534731 | |
| 100 | 0.40862 | 0.08450 | 0.493123 | |
| 1e12 | 0.41421 | 3.1e-7 | 0.414215 | → √2−1 (rel gap 2.4e-6) |

Interior maximum located at R* = 15.5817, Δv_total/v1 = 0.536258 — so there
is a radius ratio past which farther targets cost LESS fuel (but strictly
more time); the far-asymptote of 0.41421 · v1 equals the escape burn, and
Δv1 alone never exceeds it. Outward/inward symmetry is exact (0.0 relative
difference over the grid). Δv2 has its own interior maximum of 0.1900 at
R ≈ 5.88 (the circularization burn peaks for mid-ratio targets).

**RK4 validation of the complete transfer** (burn → coast → burn):

| case | e_t | steps | arrival rel r-err | v-err | radial/v | burn2 err | orbit circ. |
|------|-----|-------|-------------------|-------|----------|-----------|-------------|
| R = 1.5 | 0.200 | 358 | 4.5e-10 | 4.5e-10 | 1.5e-9 | 3.8e-9 | 4.5e-11 |
| R = 6.41 | 0.730 | 1826 | 1.3e-9 | 1.3e-9 | 5.1e-9 | 1.4e-9 | 4.5e-11 |
| R = 20 | 0.905 | 8711 | 4.2e-9 | 4.2e-9 | 2.4e-8 | 1.9e-9 | 4.5e-11 |
| Venus-like 1.383 | 0.161 | 333 | 4.7e-10 | 4.6e-10 | 1.5e-9 | 5.0e-9 | 4.5e-11 |
| inward R = 0.5 | 0.333 | 471 | 2.3e-10 | 2.2e-10 | 3.0e-10 | 1.7e-9 | 4.5e-11 |

The analytic same-ellipse reference hits r2/v(r2) exactly (0 to 1.4e-16) in
both directions (for r2 < r1 the reference is the same ellipse one
half-period later — the leg the inward flight actually flies).
The apside (apoapsis outward, periapsis inward) is reached precisely at
t_tr; the post-circularization orbit holds r2 to 4.5e-11 over one orbit.

**Optimality within two-impulse transfers**: the 121×131-grid minimum over
transfer ellipses (r_p ≤ r1, r_a ≥ r2, burns at r1 and r2, no tangency)
coincides with the Hohmann corner — measured gap 0 (R = 6.41, 20) and
7.8e-16 (R = 2), argmin exactly (r_p = r1, r_a = r2). The tangent-departure
family decreases to its minimum at r_a = r2; the fixed-apoapsis family
decreases to its minimum at r_p = r1.

**Real anchors** (IAU 2015 B3 nominal + JPL mean orbits):

- LEO (200 km, r = 6578.1 km) → GEO (35 786 km, r = 42164.1 km), Earth:
  Δv1 = 2.4546 km/s, Δv2 = 1.4773 km/s, **total 3.9319 km/s**, coast
  5.26 h; Tsiolkovsky propellant fraction 0.737 (Isp 300 s) / 0.590
  (Isp 450 s).
- Earth → Mars (helio, Mars at 1.523679 au): transfer 258.87 d,
  v∞ departure 2.9447 km/s, v∞ arrival 2.6489 km/s, heliocentric total
  5.5936 km/s; trans-Mars injection from 200 km LEO 3.6114 km/s
  (published 3.6); LEO → low-Mars-orbit total ≈ 5.70 km/s (adds the
  300 km insertion burn 2.0914 km/s).
- Earth → Venus (inward): transfer 146.08 d, v∞ departure 2.4954 km/s,
  arrival 2.7066 km/s; injection from LEO 3.5036 km/s — the inward
  burn-at-Earth is 0.45 km/s smaller than Mars's because lowering an orbit
  takes less than raising one at the same ratio asymmetry.

## Limitations

- Two-body, circular, coplanar, impulsive idealization: no eccentricity or
  inclination of real planets, no phasing/launch-window analysis, no
  non-impulsive (finite-thrust) burns, no gravity losses.
- Optimality proved only within the two-impulse class with burns at r1 and
  r2; the three-impulse bi-elliptic family (which beats Hohmann beyond
  R ≈ 11.94 when time is free) is quantified in roadmap experiment 005.
- The R → ∞ asymptotic statement is verified numerically at R = 1e12
  (2.4e-6 relative); the exact limit follows from the closed form.
- Fuel fractions are single-stage idealizations; real missions split
  staging, and real Isp/thrust regimes change the numbers.
- Real anchors use mean circular orbits; JPL Horizons ground truth for a
  real body is deferred to roadmap experiment 013.

## Future Improvements

- Bi-elliptic vs Hohmann crossover (roadmap 005): locate the 11.94/15.58
  crossover structure with the same burn machinery.
- Plane-change and combined maneuvers (roadmap 006) using the 2-impulse
  cost machinery (vector burns already computed here).
- Lambert-problem generalization for fixed transfer angles (the cost
  surface here is the 180° slice).
- Gravity-assist connection (roadmap 007): escape-burn asymptote is the
  natural zero-cost-satellite boundary.
- Guard-band real mission numbers: JPL Horizons pork-chop-style launch
  window analysis (roadmap 013).

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
- Numerical results are deterministic; re-running rewrites `results.json`
  with a fresh timestamp and generation commit (provenance), while the
  scientific output stays identical.

### Change Log

- 2026-08-13 (audit fix): `validate_transfer_rk4`'s analytic-arrival metric
  was misleading for genuine inward transfers (r2 < r1): the closed-form
  reference starts at periapsis while the inward flight starts at
  apoapsis, so it compared against the departure apside. The reference is
  now phase-shifted by a half-period for inward flights; a true inward
  case (R = 0.5) was added to the committed RK4 validation set and to the
  test suite (previously the "inward" case R = 1.3825 was actually an
  outward-branch call). No closed-form results changed.