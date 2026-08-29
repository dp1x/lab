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
| 006 | Plane-change maneuvers | **complete** (audited 2026-08-17) | Cost of changing inclination; three regimes, boundaries + pinch |
| 007 | Gravity assist / slingshot | **complete** | Max heliocentric energy change from an unpowered flyby (patched conic) |
| 008 | Ground tracks | **complete** (2026-08-21) | Spherical-Earth ground track: dual-algebra lat/lon, invariants + RK4 + real LEO/GEO/Molniya anchors |
| 009 | J2 precession | **complete** (2026-08-22) | Secular Ω̇/ω̇ from Earth's oblateness: full-force RK4 rediscovers first-order rates via independent estimator; model-order residual separated from integration error; LEO/SSO/Molniya/critical-i anchors |
| 010 | Orbit decay | **complete** (2026-08-22) | Atmospheric drag (first non-conservative force): dissipation accounting + monotonicity doctrine; decay law vs erfi/quadrature oracles (3.6 m / 500 revs); scalings, co-rotation twins, J2 mean-element transient, reentry timing; 39 tests |
| 011 | Lagrange points | **complete** (2026-08-22) | Rotating-frame CR3BP: L1–L5 equilibria, Jacobi integral, Routh stability boundary, nonlinear perturbation signatures, dimensional cross-check, adversarial mutant battery; first rotating-frame experiment; shared integrators/orbits graduated to `src/lab_utils` |
| 012 | Orbit classes | **complete** (2026-08-23) | Constraint-defined families: SSO inclination lock + finite existence limit (a_max = 12352.5 km); Molniya apsidal freeze + semi-synchronous resonance + dwell geometry (+323 s/orbit short-period Kepler-excess finding near lock); GEO 1:1 fixed point with nonzero-rate negative control; GTO vis-viva budgets anchored to Exp 004; adversarial convention battery; `j2_rhs` graduated to `src/lab_utils` |
| 013 | JPL ephemeris validation | **complete** (2026-08-24) | Propagator vs pinned NASA/JPL Horizons ISS states (ICRF/TDB, sha256-pinned snapshot, offline deterministic): J2 removes 99.33% RMS; declared drag tier worsens at primary β=100 (reported verbatim), β-band zero-crossing at band edge documented without tuning; decomposition budget separates integration/interpolation/time/init/constants from a jointly-attributed remainder |
| 014 | Eclipse timing & launch windows | **complete** (2026-08-28) | Conical + cylindrical shadow geometry via dual formulations (apparent-angle vs shadow-axis, event-time agreement sub-second); closed-form g on analytic Kepler states decouples event error from integration step; pinned-ISS first-4 event epochs agree with real NASA trajectory to 5.5-13.5 s; ISS/GEO three-tier duration pins; cone−cylinder GEO boundary shift 64.3 s; Sun model agrees with byte-pinned 2026 Horizons snapshot to 0.65 deg; 40 new tests |
| 015 | Dawn-dusk SSO launch-window targeting | **complete** (2026-08-29) | Multi-constraint mission analysis: SSO lock + LST-at-node + J2 nodal drift + eclipse-free. Year-long feasible launch-time search for h in {500,600,700,800} km from Eastern Range; 266-295 components/altitude (monotone in h); held-out equinox weeks dominate feasibility (36.7 vs 11.9/day); LST drifts through 24 h/year at sidereal-solar differential (4 min/day); shared machinery graduated: `sso_inclination_rad` + `gmst_rad_iau1982` + `sun_unit_and_dist_km` + ECI→ECEF/lat-lon + `lst_at_node_hours` to `src/lab_utils`; 34 new tests |
| 016+ | Eclipse-aware station-keeping, ground-track targeting under J2 mean-vs-osculating, … | planned | each seeds the next |

## Conventions

- Deterministic propagation (fixed steps/params); sweep methodology applies
  (CSVs → `data/`, summary JSON → `results/`).
- Validation layers: closed-form (Kepler/Hohmann equations) → invariants
  (energy, angular momentum) → real data (JPL Horizons).
- See `localdocs/roadmap.md` for sequence + hooks.