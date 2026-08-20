# Experiment Card: Planetary Gravity Assist — Two-Body Hyperbolic Flyby

> Status: complete
> Date: 2026-08-21
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/gravityAssist/`

## Research Question

What is the **maximum heliocentric velocity/energy change** that an *unpowered*
planetary flyby (patched-conic two-body hyperbola) can produce, as a function of
planet, incoming hyperbolic-excess speed v∞, periapsis radius r_p, and the 3-D
orientation of the incoming asymptote relative to the planet's heliocentric
velocity — and does the closed-form patched-conic model reproduce canonical
published mission flybys?

Framing (deliberate): a pure gravity assist involves **no propulsive Δv**. The
planet-frame excess-speed magnitude is conserved (|v∞,out| = |v∞,in|); only its
*direction* rotates by the turn angle δ. The heliocentric velocity-vector change
Δ**v**_helio = **v**_out − **v**_in is reported as a **vector change**, and the
heliocentric specific-energy change as Δε = ½(v_out² − v_in²) = **V**_p · Δ**v**_helio.
These two observables must not be conflated (Voyager 1 @ Jupiter: |Δ**v**| =
16.33 km/s but Δ|v| = +10.99 km/s).

## Background Theory

Planet-centered hyperbola with periapsis r_p and excess speed v∞. With the
nondimensional parameter x = r_p v∞²/μ_p:

- eccentricity: e = 1 + x
- turn angle (cancellation-safe form): **δ = 2·atan2(1, √(x(x+2)))** ≡ 2·arcsin(1/e);
  near-π form π − δ = 2·atan2(√(x(x+2)), 1)
- impact parameter: b = r_p·√(1 + 2/x) = (μ_p/v∞²)·√(e²−1)
  (the simple form δ = 2·arctan(μ_p/(b v∞²)) uses **b, never r_p**)
- periapsis speed (vis-viva): v_p = √(v∞² + 2μ_p/r_p)

3-D geometry (B-plane, deterministic sign convention): with reference pole k̂0
(J2000 ecliptic north), ŝ = v∞,in/v∞, t̂ = k̂0×ŝ/|k̂0×ŝ|, r̂ = ŝ×t̂,
**B** = b(cosβ t̂ + sinβ r̂), bend direction q̂ = −**B**/b (gravity bends the
trajectory toward the planet), plane normal n̂ = ŝ×q̂, and Rodrigues rotation
**v**_∞,out = R_n̂(δ) **v**_∞,in = **v**_∞,in cosδ + (n̂×**v**_∞,in) sinδ.

Heliocentric assembly: **v**_in = **V**_p + **v**_∞,in, **v**_out = **V**_p +
**v**_∞,out, hence Δε = **V**_p·Δ**v**_helio (planet exchanges orbital energy
with the spacecraft; nothing is "created").

Exact orientation landscape (α = angle between ŝ and **V**_p; φ = flyby-plane
orientation about ŝ; h = δ/2):

    Δε(α,φ) = 2 V_p v∞ sin h · [ −sin h cos α + cos h sin α cos φ ]

Global maximum at **α* = π/2 + δ/2, φ* = 0** (the optimal incoming asymptote is
*more than 90°* from the planet velocity — "prograde" is the wrong intuition);
global minimum at α = π/2 − δ/2, φ = π; neutral locus cos φ = tan(δ/2)·cot α.
The extremal bend vector is **exactly parallel to +V_p** (Cauchy–Schwarz
equality condition), giving the closed-form ceiling

    |Δε|_max = 2 V_p v∞ sin(δ/2) = 2 V_p v∞ / (1 + r_p v∞²/μ_p)

Asymptotics: δ ~ 2/x as x→∞ (|Δ**v**| ~ 2μ_p/(r_p v∞) — energy exchange vanishes
as 1/v∞ at high speed); π − δ ~ 2√(2x) as x→0 while |Δ**v**| → 2v∞ → 0 (a
180° turn of a vanishing velocity is still a vanishing change).

## References

- C. D. Murray & S. F. Dermott, *Solar System Dynamics*, Cambridge Univ. Press,
  1999 (patched-conic / sphere-of-influence formulation).
- R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*,
  Dover, 1971 (patched-conic flyby; cited as the basis of NASA's Trajectory
  Browser model).
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier,
  2021, Ch. 8 (hyperbolic flyby, turn angle, v∞ vector addition).
- JPL Publication 82-43, *Spacecraft Trajectory Design and B-Plane targeting*
  (B-plane definition, impact parameter, r_p conversion), NASA NTRS.
- NASA/JPL, Voyager mission hyperbolic encounter elements (Voyager 1 Jupiter
  1979-03-05: a = −1,092,356 km, e = 1.318976; Voyager 2 Jupiter 1979-07-09:
  a = −2,184,140 km, e = 1.330279; Voyager 1 Saturn 1980-11-12: a = −166,152 km,
  e = 2.107561) and heliocentric pre/post-encounter semimajor axes (NASA Science
  Voyager planetary-elements table; provenance attributed to JPL trajectory
  engineering).
- JPL Solar System Dynamics, astrodynamic/physical parameter tables (DE440
  solar GM 132712440041.279 km³/s²; planet-only GMs — Jupiter 126,686,531.9
  km³/s² from JUP365, Saturn 37,931,206.23 km³/s² — *not* the DE440
  planetary-system GMs, which include the major moons).
- NASA Trajectory Browser user guide (patched two-body flyby model, periapsis
  feasibility floors, documented fidelity limitations).

## Assumptions

- Two-body planet-centered gravity during the encounter (patched conic);
  solar gravity during the encounter neglected — **idealization** (this is the
  model; roadmap 011/013 extend the fidelity).
- Planet on a circular heliocentric orbit, **V**_p = √(μ_☉/a_p) — **idealization**
  (Mars e = 0.0934 implies ±~10% slow variation of the true V_p).
- Impulsive (asymptotic) patched-conic encounter; no finite-duration effects —
  **idealization**.
- Planet-only GMs (Jupiter/Saturn): moons excluded from the central body —
  **convention, documented** (DE440 *system* GMs deliberately not used).
- Atmospheric/operational flyby floors: sweep uses r_p ≥ 1.02 R_eq as a
  **screening convention**, not a physical law (real missions use
  atmosphere/heat-rate constrained corridors; NASA Trajectory Browser rejects
  periapses below the atmosphere).
- Anchor elements are a published navigation product without released
  covariance — source-to-source spread (~0.13% Voyager 1 Jupiter closest
  approach across NASA products) is reported as **spread, not statistical
  uncertainty**.

## Methodology

1. **Closed-form core (L1/L2):** implement the numerical contract above
   (x-parameterization, cancellation-safe δ, B-plane geometry, Rodrigues
   rotation). Cross-check δ against 2·arcsin(1/e) and, at extreme x, against
   50-digit mpmath.
2. **Orientation landscape:** evaluate the closed-form Δε(α,φ) on a 1°×1° grid
   (181×360; grid-resolution loss at the optimum ≈ 0.008%, far below the 0.1%
   requirement) for representative cases; verify grid max ≈ analytic ceiling
   and that the analytic optimum (α*, 0) reproduces Δ**v**_helio ∥ +**V**_p.
3. **Independent numerical propagation (L3/L4):** for representative flybys,
   construct the exact conic state at radius R₀ on the *inbound* leg, integrate
   the raw two-body EOM through periapsis with the verified Exp 006 3D Cowell
   RK4 (deterministic angular step rule Δt = θ·r/v, θ = 0.02), then **recover
   the hyperbola from the final state** (ε, **h**, **e**) and compute the
   asymptotic directions analytically — never by comparing finite-radius
   velocity directions (which carries O(10⁻²) rad truncation at 5–10 SOI).
   Check recovered δ, r_p, v∞ and planet-frame energy conservation; repeat from
   R₀ ∈ {50, 100, 200}·r_p to show the recovered δ is patch-radius-insensitive.
4. **Canonical anchors (L5):** reproduce Voyager 1 Jupiter 1979, Voyager 2
   Jupiter 1979, Voyager 1 Saturn 1980 from the published (a, e) elements with
   planet-only GMs; compare v∞, r_p, δ, |Δ**v**|, and Δε (from the published
   heliocentric a_in/a_out) against the published/derived values.
5. **Parameter sweep:** 5 planets × 32 log-spaced v∞ ∈ [0.5, 30] km/s × 64
   log-spaced r_p ∈ [1.02 R_eq, 50 R_eq] = 10,240 cases; report the analytic
   ceiling |Δε|_max per case and the optimal geometry; the orientation grid
   serves as a numerical cross-check at representative cases (the analytic
   maximum makes a brute-force 6.6×10⁸-point orientation search unnecessary).
6. **Pathological regimes (L6):** x from ~10⁻¹⁸ to ~10¹² (near-parabolic →
   ultra-weak), δ monotonicity, no NaN/Inf, stable-vs-unstable formula
   agreement, mpmath verification at the extremes.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib, mpmath (verification only);
  reuses the verified 3D Cowell RK4 from Experiment 006 by explicit-path import
  (single source of truth; no scaffolding rebuilt).
- Runtime: `uv run python experiment.py` (~1–2 min single core).
- Determinism: pure float64 closed forms; no RNG anywhere; fixed grids; the L3
  integrator uses a deterministic angular step rule (no adaptive randomness).
- Dependencies: numpy, matplotlib, mpmath (all already in the lockfile).

## Validation Method

- Unit tests: `tests/test_gravity_assist.py` (pytest).
- Analytic identities: δ forms agree to < 1e-12 rad (and to ~1e-30 vs mpmath at
  extremes); |v∞,out| = |v∞,in| to machine precision; Δε = ½(v_out²−v_in²) =
  **V**_p·Δ**v**; grid max vs analytic ceiling within 0.1%.
- Independent integration: recovered δ/r_p/v∞ vs closed form; patch-radius
  insensitivity (R₀ ∈ {50,100,200}·r_p); planet-frame energy conservation.
- Published data: three Voyager anchors (see Results).

## Results

See `results/results.json` and figures. Headline numbers:

- **Orientation optimum:** α* = 90° + δ/2 (e.g. Voyager-1-like Jupiter flyby,
  δ = 98.6° → α* = 139.30°, i.e. 40.7° short of antiparallel); ceiling
  Δε_max = 2V_p·v∞/(1+x).
- **Anchors (derived from published elements, planet-only GMs):**
  Voyager 1 Jupiter: v∞ = 10.7692 km/s, r_p = 348,435 km (published closest
  approach 348,890 km / ~350,000 km — 0.13%/0.45% source spread), δ = 98.605°,
  |Δ**v**| = 16.330 km/s, Δε = +200.83 km²/s² (Δ|v| = +10.99 km/s).
  Voyager 2 Jupiter: v∞ = 7.6160 km/s, r_p = 721,376 km, δ = 97.480°,
  |Δ**v**| = 11.450 km/s, Δε = +151.76 km²/s².
  Voyager 1 Saturn: v∞ = 15.109 km/s, r_p = 184,023 km, δ = 56.651°,
  |Δ**v**| = 14.338 km/s, Δε = +26.12 km²/s² — a 14.3 km/s vector rotation
  buying only ~1.25 km/s of heliocentric speed: orientation dominates.
- **Ceiling examples (r_p = 1.02 R_eq):** Jupiter at v∞ = 10 km/s →
  |Δε|_max ≈ 247 km²/s²; Earth at v∞ = 3 km/s → ≈ 156 km²/s²; the ceiling
  falls as ~1/v∞ at high v∞ for every planet.
- **Interior optimum (sweep finding):** for fixed r_p the ceiling
  2V_p v∞/(1 + r_p v∞²/μ_p) is maximized at **v∞* = √(μ_p/r_p)** —
  7.8 km/s (Earth), 7.3 (Venus), 3.5 (Mars) at r_p = 1.02 R_eq, all inside
  the practical flyby regime; Jupiter/Saturn (v* = 41/25 km/s) peak at or
  beyond the v∞ = 30 km/s sweep edge. This is the quantitative reason
  Earth/Venus gravity assists dominate real mission design.
- **L3:** recovered δ agrees with the closed form to ≲ 1e-9 rad and is
  insensitive to the patch radius (element-recovery method), while direct
  finite-radius direction comparison at the same radii is wrong by 1e-3–1e-1
  rad — the documented reason the recovery method is mandatory.

## Limitations

- Patched conic only: no solar third-body perturbation during the encounter, no
  planetary oblateness (J2), no moons, no finite-burn/loiter effects, no
  ephemeris-level reconstruction (roadmap 011 → 013 extend fidelity).
- Circular heliocentric planet orbits; real V_p varies (Mars ±~10%).
- Anchor comparison inherits published-element provenance; no covariance was
  released, so agreement is quoted against *derived* values with source spread
  stated separately — no invented statistical uncertainties.
- Operational flyby floors (atmosphere, radiation, heat rate) are not modeled;
  r_p ≥ 1.02 R_eq is a screening convention only.
- The "sphere of influence" is a convention, not a physical boundary; the L3
  integrator validates the two-body model itself, which is exact inside the
  patched-conic idealization regardless of patch radius.

## Future Improvements

- Restricted three-body encounter dynamics (roadmap 011, Lagrange points first).
- Ephemeris-level reconstruction of a historical encounter (roadmap 013, JPL
  Horizons ground truth).
- B-plane targeting under launch-dispersion covariances (needs real navigation
  data).

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
