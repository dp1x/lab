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
| 002 | Kepler orbit validation | Does Newtonian gravity reproduce elliptical orbits and Kepler's laws? | Kepler's three laws; energy + angular-momentum conservation |
| 003 | Kepler's equation solvers | Newton vs bisection vs series; convergence study | Closed-form roots, known test orbits |
| 004 | Hohmann transfer | Least-fuel orbit-to-orbit transfer; Δv budget | Closed-form Δv equations |
| 005 | Bi-elliptic vs Hohmann | Crossover radius law (β > 15.58) | Closed-form comparison |
| 006 | Plane-change maneuvers | Cost of changing inclination | Closed-form Δv |
| 007 | Gravity assist / slingshot | Velocity boost from a flyby | Patched-conic / known flyby numbers |
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