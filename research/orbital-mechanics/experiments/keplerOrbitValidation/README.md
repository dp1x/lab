# Kepler Orbit Validation: Newtonian Two-Body Gravity Reproduces Kepler's Laws

> Status: complete
> Date: 2026-08-13
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/keplerOrbitValidation`

## Research Question

Does Newtonian two-body gravity — as integrated by a deterministic fixed-step
RK4 propagator — reproduce elliptical orbits and Kepler's three laws, and do
the resulting orbits agree with the closed-form Keplerian solution and with
real Earth-system constants?

## Background Theory

The two-body problem (central force a = −μ r / r³) has the closed-form solution
of an ellipse in the orbital plane, parametrized by Kepler's equation

    M = E − e sin E,    with M = n t (t measured from periapsis), n = √(μ/a³)

position in the perifocal frame: x = a(cos E − e), y = a√(1−e²) sin E, where
r = a(1 − e cos E). The three laws take closed forms:

- Kepler I:   bound orbits are ellipses, r = p/(1 + e cos ν), p = h²/μ = a(1−e²).
- Kepler II:  areal velocity dA/dt = h/2 (h = |r × v|); a full orbit sweeps πab.
- Kepler III: T = 2π√(a³/μ), i.e. T²/a³ = 4π²/μ independent of eccentricity.

Conserved invariants: specific energy ε = ½v² − μ/r = −μ/(2a) (vis-viva:
v² = μ(2/r − 1/a)) and specific angular momentum h = √(μ a (1−e²)). The
eccentricity vector e_vec = (v × h_vec)/μ − r̂ recovers both magnitude and
orientation (points at periapsis).

Real-world anchor: with the IAU 2015 nominal solar mass parameter and the IAU
2012 definition of the astronomical unit, T(1 au) = 2π√(au³/GM☉) =
365.256898 d, which sits +1.47e-6 above the sidereal year (365.256363 d) — the
residual is physical (Earth–Moon system mass, mean-vs-nominal a), not numeric.

## References

- R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*,
  Dover, 1971 — two-body problem, Kepler's equation, vis-viva (Ch. 2).
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier,
  2021 — two-body problem and orbital position as a function of time (Ch. 2–3).
- C. D. Murray, S. F. Dermott, *Solar System Dynamics*, Cambridge UP, 1999 —
  two-body problem and areal velocity (Ch. 2).
- IAU 2015 Resolution B3, "Recommended nominal conversion constants for
  selected solar and planetary properties" (GM☉ = 1.3271244e11 km³/s²,
  exact nominal). iau.org resolution text; Mamajek et al., arXiv:1510.07674.
- IAU 2012 Resolution B2, "Astronomical unit" (au = 149597870.7 km exactly).
- JPL SSD, "Approximate Positions of the Planets" (Standish & Williams 1992),
  mean ecliptic J2000 EM-barycenter eccentricity e = 0.01671123.
- *Astronomical Almanac for the Year 2025*, USNO & HM Nautical Almanac Office,
  2024 — sidereal year 365.256363 d (pp. C2, L9).
- Values verified against primary sources on 2026-08-13 (see Knowledge note).

## Assumptions

- Two-body idealization: Sun + point Earth; no other planets, no Moon mass, no
  oblateness (verified as the governing model — the Earth anchor's 1.5e-6
  residual is this physics, not numerical error).
- Canonical units (μ = 1, a = 1, T = 2π) for the core sweep; real units for the
  Earth case only (verified).
- Deterministic fixed-step RK4; per-case stepsize chosen so the periapsis
  passage is resolved: steps/orbit ~ (1−e)^(−3/2) (plausible, measured in
  `case_steps_per_orbit`: at e = 0.85 uniform 512 steps/orbit give ~100%
  position error over 5 orbits; 4096 give ~8e-5).
- Float64 arithmetic; Newton's method for Kepler's equation converges ~6
  iterations at e ≤ 0.85 (verified tol 1e-14).

## Methodology

All runs start at periapsis r_p = a(1−e) with the vis-viva tangential speed
v_p = √(μ(1+e)/(a(1−e))), propagated by classical RK4 over the planar two-body
system. Validations:

1. **Convergence**: one orbit, e = 0.6, h ∈ T/128 … T/2048; global position
   error vs the analytic solution, log-log slopes.
2. **Kepler I**: pointwise max relative position error of the RK4 trajectory
   vs the closed-form Kepler solution over 5 orbits (e = 0.3, 0.6, 0.85), plus
   the conic equation r = p/(1 + e cos ν) along the trajectory and the
   eccentricity vector recovering e_input; circular case (e = 0) checked for
   constant radius.
3. **Kepler II**: one orbit split into 12 equal-time intervals of 64 steps;
   each interval's swept sector area via triangle sums ½|r_i × r_{i+1}|;
   full-orbit area vs πab; cumulative area slope vs h/2.
4. **Kepler III**: sweep a ∈ {0.5, 1, 2, 4, 8} × e ∈ {0.10, 0.30, 0.60, 0.85}
   (20 cells); measured period from two successive periapsis passages with
   parabolic vertex refinement; T²/a³ vs 4π²/μ and log-log slope vs 1.5.
5. **Conservation**: 10 orbits at e = 0.6; max relative drift of ε and h.
6. **Earth anchor**: a = au, μ = GM☉, e = 0.01671123; measured period vs
   2π√(a³/μ) and vs the sidereal year; measured eccentricity vs the JPL value.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib (no scipy)
- Runtime: `uv run python experiment.py`
- Determinism: no RNG; fixed Newton caps and tolerances; subprocess
  cross-process bit-identical check in the test suite.
- Dependencies: numpy, matplotlib (already in pyproject; no new deps).
- Reuses `lab_utils.metrics` (max_abs_error, convergence_rate) and
  `lab_utils.results` (save_json_result).

## Validation Method

Independent checks in `tests/test_kepler_orbit_validation.py` (20 tests,
all green):

- `kepler_solution` satisfies the ODE numerically (central differences) and
  the periapsis initial conditions; energy −μ/2a constant on the analytic
  solution; Kepler solver round-trips E → M → E.
- State-to-elements recovery (a, e, ε, h from initial conditions) for
  e ∈ {0, 0.3, 0.6, 0.85}; vis-viva ↔ energy consistency.
- RK4 vs analytic solution over 5 orbits (~1e-7 relative at e = 0.6);
  conic equation along the trajectory (< 1e-5); circular case r ≡ a.
- Equal areas: interval deviations < 1e-3; full-orbit area vs πab < 1e-3;
  areal velocity vs h/2 < 1e-3.
- Kepler III: log-log slope 1.5 ± 0.01; T²/a³ = 4π²/μ within 1e-3 on every
  cell; per-cell period vs theory < 1e-3.
- Conservation: relative drift < 1e-6 (ε) and < 1e-8 (h) over 10 orbits.
- Convergence: slopes → 4 (coarse-slope tolerance ±0.5, finest ±0.2).
- Earth: measured period vs closed form < 1e-4; e_meas vs 0.01671123 < 1e-3;
  closed form sits in the band (1e-6, 2e-6) above the sidereal year.
- Determinism: identical output in a fresh interpreter (importlib-loaded
  module, explicit path).

An independent adversarial review (subagent, 2026-08-13) re-derived the RK4
tableau, Kepler solver, elements, true-anomaly and period-measurement algebra,
re-ran the suite, and reproduced every recorded number with a second
integrator family (velocity Verlet), a bisection Kepler solver, and its own
periapsis detector: no numeric discrepancies, two documentation items fixed
(docstring wording, field naming).

## Results

All results in `results/results.json`; figures in `results/figures/`.

**Kepler I** (a = 1, 5 orbits, periapsis-resolved steps):

| e_input | steps/orbit | max rel pointwise error | max rel conic error | e_measured |
|---------|-------------|------------------------|---------------------|------------|
| 0.0     | 512         | 9.3e-10 (radius drift) | –                   | 0.0        |
| 0.3     | 875         | 2.58e-08               | 4.01e-09            | 0.300000   |
| 0.6     | 2024        | 1.77e-07               | 6.04e-09            | 0.600000   |
| 0.85    | 8814        | 3.04e-06               | 1.14e-08            | 0.850000   |

**Kepler II** (a = 1, e = 0.6): 12 equal-time intervals, max relative area
deviation 7.4e-5; full-orbit sector area 2.513219 vs πab = 2.513274 (rel
2.2e-5); areal-velocity slope 0.3999964 vs h/2 = 0.4 (rel 9.1e-6).

**Kepler III** (20 cells, a ∈ [0.5, 8], e ∈ [0.10, 0.85]): log-log slope
1.50000 (theory 1.5); T²/a³ ∈ [39.47841728, 39.47841763] vs 4π² = 39.47841760
(max rel error 8.2e-9); every measured period matches 2π√(a³/μ) to ≤ 4.1e-9.

**Conservation** (e = 0.6, 10 orbits): max relative drift 1.19e-9 (ε),
6.96e-11 (h). RK4 is not symplectic; these drifts are small discretization
artifacts, not round-off-level exactness (verified: h = 7e-11 ≫ 2.2e-16).

**Propagator convergence** (e = 0.6): per-orbit max position error
1.31e-3 → 1.12e-8 for h = T/128 → T/2048; log-log slopes 4.39, 4.24, 4.13,
4.07 → 4 (pre-asymptotic O(h⁵) contributions at coarse h).

**Earth anchor** (real units): T_pred = T_meas = 365.256898 d (rel error
2.5e-10); vs sidereal year 365.256363 d: +1.47e-6 (+46 s) — the expected
two-body idealization residual; e_measured = 0.01671123 vs JPL 0.01671123.

## Limitations

- Two-body, planar, point-mass only: no J2 oblateness, no third bodies, no
  Moon mass, no relativity — the Earth-anchor residual (1.5e-6) is exactly
  this missing physics.
- Fixed-step integration with periapsis-resolved step counts; no adaptive
  stepping, no long-horizon (centuries) propagation.
- Mean-anomaly-frame analytic reference assumes periapsis at t = 0; no
  argument-of-periapsis / inclination variation explored (planar only).
- RK4 phase error grows with integrated time (linear in orbits); 5-orbit and
  10-orbit horizons reported, not secular horizons.
- Period measurement uses successive periapsis passages; alternative
  definition (e.g., node crossings) not compared.

## Future Improvements

- Kepler solver internals now measured in detail (roadmap 003, complete):
  Newton order 2, bisection halving 1/2, Fourier-Bessel series decay q(e);
  Newton with the M + e sin M starter is the recommended path.
- Add eccentric-anomaly-based uniform-time sampling to decouple step count
  from eccentricity (regularization in eccentric anomaly).
- Hohmann transfer economics (roadmap 004), reusing this propagator.
- Extend to 3D classical elements (Ω, i, ω) and ground tracks (roadmap 008).
- JPL Horizons comparison for a real body (roadmap 013).

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
- Numerical results are deterministic; re-running rewrites `results.json` with
  a fresh timestamp and generation commit (provenance), while the scientific
  output stays identical.