# Comparative Numerical Integration of ODEs: Convergence Order and Energy Preservation

> Status: complete
> Date: 2026-08-05 (figures + results regenerated 2026-08-09)
> Domain: numerics
> Experiment dir: `research/numerics/experiments/odeIntegratorStudy`

## Research Question

How do five first-order ODE integrators — forward Euler, RK2 (midpoint), RK4,
symplectic (semi-implicit) Euler, and velocity Verlet — compare in (a) global
convergence order and (b) long-term preservation of total mechanical energy, when
applied to the linear harmonic oscillator?

## Background Theory

The simple harmonic oscillator `x'' + omega^2 x = 0` has the closed-form solution

    x(t) = x0 cos(omega t) + (v0 / omega) sin(omega t)

with total mechanical energy (unit mass)

    E = 1/2 v^2 + 1/2 omega^2 x^2 = const.

Classical Runge–Kutta theory predicts global error O(h^p) where p is the method's
order: Euler p=1, midpoint/RK2 p=2, RK4 p=4. Geometric/symplectic integration theory
(Hairer, Lubich, Wanner 2006) predicts that symplectic methods (symplectic Euler,
velocity Verlet) produce *bounded* energy error over long horizons, while non-symplectic
methods (Euler, RK2, RK4) exhibit secular energy drift.

## References

- J. C. Butcher, *Numerical Methods for Ordinary Differential Equations*, 3rd ed., Wiley, 2016.
- E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, 2nd ed., Springer, 2006.

## Assumptions

- Unit mass, omega = 1 (verified), x0 = 1, v0 = 0 (verified).
- Float64 arithmetic; double-precision round-off negligible at the stepsizes used.
- Energy study horizon 200π = 100 periods is "long-term" for this linear problem (plausible).

## Methodology

For each method and each requested stepsize h in {0.1, 0.05, 0.025, 0.0125}, integrate
the oscillator as a first-order system from t=0 to t=2pi (one period). The grid uses
`n = round(2pi/h)` steps of effective size `h_eff = 2pi/n` so the final grid point is
exactly 2pi (effective stepsizes ~0.0997, ~0.0499, ~0.0250, ~0.0125; recorded in
`results.json`). The reported error is the maximum absolute error over all grid
points in [0, T] (the standard definition of global error). The convergence
order is the log-log slope of error vs effective stepsize across the four refinements.

Energy behaviour is measured over the long horizon t in [0, 200pi] with h = 0.05,
recording max deviation and final deviation of E from its initial value.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib
- Runtime: `uv run python experiment.py`
- Determinism: no RNG involved; integrators are purely deterministic.

## Validation Method

- Exact analytic solution used as reference.
- `uv run pytest tests/test_experiment.py`:
  - analytic solution satisfies initial conditions and conserves energy;
  - measured convergence order matches theoretical order within 0.5;
  - symplectic methods preserve energy far better than non-symplectic ones;
  - runs are bit-reproducible (determinism test).

## Results

Measured global convergence orders (log-log slope over h = 0.1 … 0.0125) and
measured energy behaviour (E0 = 0.5, h = 0.05, horizon t ∈ [0, 200π]):

| method          | theoretical order | measured order | energy: max \|ΔE\| | energy: final \|ΔE\| |
|-----------------|-------------------|----------------|---------------------|------------------------|
| euler           | 1                 | 1.06           | 2.119e+13           | 2.119e+13             |
| rk2_midpoint    | 2                 | 2.01           | 9.915e-03           | 9.915e-03             |
| rk4             | 4                 | 4.01           | 1.363e-06           | 1.363e-06             |
| symplectic_euler| 1                 | 1.02           | 1.282e-02           | 1.630e-03             |
| velocity_verlet | 2                 | 2.00           | 3.125e-04           | 1.338e-06             |

All five methods match theoretical convergence order (Euler 1, RK2 2, RK4 4,
symplectic Euler 1, Verlet 2), confirming the implementations. At the coarsest grid
(h_eff ~ 0.0997) RK4 achieves max error 4.04e-06 vs 3.66e-01 for Euler over one period.

The energy study confirms geometric integration theory: non-symplectic methods
(Euler, RK2, RK4) show secular, integrating-time-proportional energy error (Euler
blows up to 2.1e13; RK4's 1.36e-06 is small but grows linearly with t), while
symplectic methods have bounded, oscillating energy error — velocity Verlet's max
deviation 3.125e-04 with final deviation 1.338e-06 (error does not accumulate),
and symplectic Euler oscillates (1.282e-02 max, 1.630e-03 final).

See `results/results.json` for exact numbers and `results/figures/` for the two
figures (convergence.png: error vs stepsize; energy_deviation.png: |E(t)-E(0)| over
the long horizon).

## Limitations

- Linear problem only; order and energy results do not directly transfer to
  nonlinearity or stiff systems.
- No adaptive stepping; fixed-step integration only.
- Energy preservation is one measure of geometric fidelity; phase-space area
  preservation and long-term frequency error are not measured here.

## Future Improvements

- Extend to nonlinear (Duffing) oscillator and n-body problems.
- Measure phase error / frequency drift and Poincare maps for symplectic methods.
- Add variable-stepsize adaptive integrators (Dormand–Prince) for comparison.
- Scale study: measure error growth with integration time for each method.

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command to reproduce (from the repo root):
  `uv sync && uv run pytest && uv run python research/numerics/experiments/odeIntegratorStudy/experiment.py`
- Numerical results are deterministic; re-running rewrites `results.json` with a
  fresh timestamp and the generation commit (provenance). This is intended: the
  metadata identifies the actual run, while the scientific output stays identical.