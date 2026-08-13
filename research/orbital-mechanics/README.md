# Orbital Mechanics — Experiments (FLAGSHIP)

The laboratory's flagship domain: how objects move under gravity in space —
orbits, transfers, perturbations, trajectories. Compact physics with closed-form
checks and real-world data (NASA/JPL Horizons) as ultimate validation.

## Why this domain

1. **Verification goldmine:** NASA/JPL Horizons API is free, reachable via plain
   `curl`, returns real ephemeris data as ground truth.
2. **Pro-grade by construction:** the same equations mission planners use.
3. **Deep enough for 50–100+ experiments:** orbit types, solvers, transfers,
   perturbations, Lagrange points, ground tracks, eclipses, sweeps.

## Experiments

| # | Experiment | Status | Question |
|---|-----------|--------|----------|
| 002 | Kepler orbit validation | **complete** | Does Newtonian gravity reproduce elliptical orbits and Kepler's laws? |
| 003 | Kepler's equation solvers | **complete** | Newton vs bisection vs series; convergence study |
| 004 | Hohmann transfer | **complete** | Least-fuel orbit-to-orbit transfer; Δv budget |
| 005 | Bi-elliptic vs Hohmann | **complete** | When do three burns beat two? (R_bp = 11.94, R* = 15.58) |
| 006 | Plane-change maneuvers | planned | Cost of changing inclination |
| 007 | Gravity assist / slingshot | planned | Free velocity boost from a flyby |
| 008 | Ground tracks | planned | The path a satellite traces over Earth |
| 009 | J2 precession | planned | Orbit-node drift from Earth's bulge |
| 010 | Orbit decay | planned | Drag over time → re-entry timeline |
| 011 | Lagrange points | planned | 3-body stability zones |
| 012 | Orbit classes | planned | Sun-synchronous / Molniya / GTO |
| 013 | JPL ephemeris validation | planned | Propagator vs NASA real positions |
| 014+ | Eclipse timing, launch windows, trajectory optimization, … | planned | each seeds the next |

## Conventions

- Deterministic propagation (fixed steps/params); sweep methodology applies
  (CSVs → `data/`, summary JSON → `results/`).
- Validation layers: closed-form (Kepler/Hohmann equations) → invariants
  (energy, angular momentum) → real data (JPL Horizons).
- See `localdocs/roadmap.md` for sequence + hooks.