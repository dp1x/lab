# Numerics — Experiments

Foundation domain: verified numerical methods (integrators, solvers, error
metrics) that all applied domains consume. Results are validated against
analytic solutions and invariants.

## Experiments

| # | Experiment | Status | Question | Results |
|---|-----------|--------|----------|---------|
| 001 | [ODE integrator comparison](experiments/odeIntegratorStudy/) | complete | Convergence order + energy preservation of Euler/RK2/RK4/symplectic-Euler/Verlet on the harmonic oscillator | Orders confirmed (1,2,4,1,2); symplectic bounded energy error |

## Next candidates

- Nonlinear (Duffing) oscillator: phase error, frequency drift, Poincaré sections
  (deferred — orbital mechanics launched first; see `localdocs/roadmap.md`).
- Scientific computing utilities grow here as experiments need them
  (`src/lab_utils/`).