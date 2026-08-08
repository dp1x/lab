---
tags: [numerics, physics, odes, symplectic-integration]
date: 2026-08-05
aliases: [ode-integrators, symplectic-integrators]
links:
  - "[[harmonic-oscillator]]"
---

# Numerical Integration of ODEs: Order and Symplecticity

## Summary

Five fixed-step integrators were compared on the simple harmonic oscillator
(x'' + ω²x = 0) — forward Euler, midpoint/RK2, RK4, symplectic Euler, and velocity
Verlet. All matched their theoretical convergence order (1, 2, 4, 1, 2), and the
energy study confirmed **symplectic integrators preserve energy in a bounded,
non-accumulating way over long horizons** while non-symplectic methods show secular
energy drift.

## Content

### Convergence order

Global error E(h) ≍ C·h^p where p is the method order. Measured (one period, ω=1):

| method           | order | err(h=0.1) |
|------------------|-------|------------|
| euler            | 1.07  | 3.68e-01   |
| rk2_midpoint     | 2.01  | 8.14e-03   |
| rk4              | 4.01  | 4.08e-06   |
| symplectic_euler | 1.02  | 5.20e-02   |
| velocity_verlet  | 2.00  | 2.01e-03   |

RK4 dominates per-step accuracy; interestingly symplectic Euler is ~5× more
accurate than forward Euler at the same order thanks to its implicit nature.

### Energy preservation (geometric integration)

Over t ∈ [0, 200π] at h = 0.05 (E0 = 0.5): non-symplectic methods accumulate energy
error proportional to integration time (forward Euler diverges, 2.1e13; RK4 grows
linearly, 1.4e-6·t). Symplectic methods (Verlet, symplectic Euler) oscillate in a
bounded energy shell: velocity Verlet final |ΔE| = 6.9e-7 with max 3.1e-4 over the
whole horizon.

### Method of measurement

Central lesson: **always evaluate the reference on the integrator's time grid**
`t_i = i·h`, not an independently generated grid. `np.linspace(0, T, n+1)` with
non-dividing T misaligns sample points and swamps high-order convergence — first
measured orders came out 0.65–0.87 until the grid was aligned.

## Source Experiments

- `research/numerics/experiments/odeIntegratorStudy/` — full experiment note, code,
  tests, results. Runnable: `uv run python research/numerics/experiments/odeIntegratorStudy/experiment.py`

## Key Takeaways

- Convergence order is a necessary but insufficient validation; always also check
  long-horizon invariants (energy, momentum, phase).
- Symplectic methods are preferred for Hamiltonian systems over long times even when
  lower order presents more per-step error.
- Reproducible numeric experiments need a rigorous reference grid.

## See Also

- [[Oscillator]]
- [Hairer, Lubich & Wanner, *Geometric Numerical Integration* (Springer 2006)]
- [Butcher, *Numerical Methods for ODEs* (Wiley 2016)]