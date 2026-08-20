# Research Roadmap

The experiment sequence the laboratory follows. Any agent picking up a task should
read this first, then `charter.md`, then continue from the current active state.

## North Star

**Orbital mechanics is the flagship domain.** Compact physics + closed-form
verification + real-world data (NASA/JPL Horizons) as the answer key. Numerics is
the foundation that supplies verified numerical methods; energy systems is the
second pillar.

## Domain naming (fixed)

- `numerics/` — the substrate: verified numerical methods (was `physics/`; too
  broad — everything is physics). Exp 001 lives here as the seed.
- `orbital-mechanics/` — the flagship (was `aerospace/`; implies airplanes/CFD,
  which is supercomputer-only and NOT our domain).
- `energy`, `computer-architecture`, `cybersecurity` — roadmap text only, no
  empty folders until real content exists.

## Sequence

### Phase 1 — Numerics (foundation)

| # | Experiment | Status | Validation |
|---|-----------|--------|-----------|
| 001 | Numerical integrator study (Euler, RK2, RK4, symplectic Euler, velocity Verlet) | **COMPLETE** | analytic solution, convergence order, energy invariants |

Continue foundation work only where it directly serves the flagship.

### Phase 2 — Orbital Mechanics (FLAGSHIP)

| # | Experiment | Question | Validation |
|---|-----------|----------|-----------|
| 002 | Kepler orbit validation | Does Newtonian gravity reproduce elliptical orbits and Kepler's laws? | **COMPLETE** (2026-08-13): analytic conic pointwise (≤3e-6), equal areas (≤7.4e-5), T²/a³ = 4π²/μ (≤8.2e-9), invariants ≤1.2e-9 over 10 orbits, IAU Earth anchor 365.256898 d vs sidereal year |
| 003 | Kepler's equation solvers | Newton vs bisection vs series; convergence study | **COMPLETE** (2026-08-13): Newton order 2 (≤9 iters at e ≤ 0.99), bisection halving 0.49992 (49–50 iters), series q_meas matches Watson q(e) ≤ 0.14%, cross-solver agreement ≤ 1.8e-13 |
| 004 | Hohmann transfer | Least-fuel orbit-to-orbit transfer; Δv budget | **COMPLETE** (2026-08-13): closed forms vs RK4 (r/v err ≤ 4.2e-9), (R−1)/2 and √2−1 asymptotes (≤0.9999×, ≤2.4e-6), peak R* 15.5817 / 0.536258, inward symmetry 0.0, 2-impulse optimality grid gap ≤ 7.8e-16, LEO→GEO 3.9319 km/s, E→Mars 258.87 d / v∞ 2.945/2.649 / TMI 3.6114 km/s |
| 005 | Bi-elliptic vs Hohmann | Crossover radius law (β > 15.58) | **COMPLETE** (2026-08-13): R_bp 11.9387654726 & R* 15.5817187388 (= Exp 004 Hohmann peak 15.5817187369 via corner identity, 1e-29 @ 50 digits), crossover curve s_c(R) 12→815.8 … R*→R*, region signs on 90×400 grid (worst margin 1.7e-8), hump onset 9.53 < R_bp, max saving 4.09% v1 @ R=50.1, RK4 3-burn validation (≤4e-8), Wikipedia 14× example exact (4.117530 vs 4.133716 km/s) |
| 006 | Plane-change maneuvers (+ combined transfer+plane change) | Cost of changing inclination & the global optimum vs bi-elliptic super-synchronous | **COMPLETE** (2026-08-16, adversarial audit 2026-08-17): three regimes (two-burn / finite-s 3-burn / s→∞), boundaries di_c(R), di_inf(R); finite-s window pinches shut at R≈6.21 (re-audited); s→∞ limit = Exp 005 bi-parabolic (1e-16); R=2,Δi=47.5° finite dip beats two-burn 1.77%; SES-8 super-sync anchor wins 5.21%; 3D RK4 validation (≤1e-11); 14 tests |
| 007 | Gravity assist / slingshot | Velocity boost from a flyby | **COMPLETE** (2026-08-21): patched-conic 2-body hyperbolic flyby; exact 3-D orientation landscape Δε(α,φ) with global max at α\*=90°+δ/2 (bend ∥ V_p, Cauchy–Schwarz); ceiling Δε_max = 2V_p·v∞/(1+r_p·v∞²/μ_p) with interior optimum v∞\* = √(μ_p/r_p); cancellation-safe δ = 2·atan2(1,√(x(x+2))); B-plane sign convention; L3 Cowell + element recovery (δ to 3e-11 rel, patch-radius-insensitive, RK4 order 4 verified); Voyager 1/2 Jupiter + Voyager 1 Saturn anchors reproduced (Δε +200.83/+151.76/+26.12 km²/s²); 33 tests |
| 008 | Ground tracks | Path a satellite traces over Earth | Spherical geometry, real orbit params |
| 009 | J2 precession | Node drift from Earth's bulge | Analytic secular rates |
| 010 | Orbit decay | Drag over time → re-entry timeline | Published decay benchmarks |
| 011 | Lagrange points | 3-body stability zones | Known L1–L5 positions |
| 012 | Orbit classes | Sun-synchronous / Molniya / GTO specifics | Closed-form + real data |
| 013 | JPL ephemeris validation | Full propagator vs NASA's published positions | Real Horizons data as ground truth |
| 014+ | Eclipse timing, launch windows, trajectory optimization, solar sails, … | each seeds the next | known physics + JPL |

### Phase 3 — Energy Systems (second pillar)

Solar forecasting → battery degradation modelling → power flow on IEEE test grids
→ economic dispatch. Stands alone; does not block the flagship.

### Phase 4 — Computer Architecture

CPU pipeline simulator, cache simulator, scheduling algorithms.

### Phase 5 — Cybersecurity

Cryptographic analysis, secure protocol modelling, vulnerability-testing frameworks.

## Sweep methodology

One experiment = one research question, but the implementation may sweep hundreds
of parameter combinations. Big CSVs live in `data/` (gitignored); a small
sampled/summary JSON is committed to `results/` for reproducibility. This is how
we get 50–100 experiments of real value instead of 5.

## Hooks

- **JPL verification layer**: NASA/JPL Horizons API is free and returns real
  positions via plain `curl`/HTTP. Use it as the ultimate validation source.
  Quirk: START_TIME/STOP_TIME must wrap a range (start < stop).
- **Reuse**: experiments build on `src/lab_utils/` and templates — never rebuild
  scaffolding.
- **Consolidation**: every ~5 experiments, a synthesis report under
  `localdocs/reports/`.
- **Deterministic only**: fixed seeds, no time-dependent nondeterminism. Reality
  is the verification layer.