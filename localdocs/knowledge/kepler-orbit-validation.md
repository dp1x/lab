---
tags: [orbital-mechanics, two-body, kepler, validation]
date: 2026-08-13
aliases: [kepler-orbit-validation, kepler-laws]
links:
  - "[[ode-integration-basics]]"
---

# Kepler Orbit Validation: Newtonian Gravity Reproduces Kepler's Laws

## Summary

A deterministic fixed-step RK4 two-body propagator was validated against the
closed-form Keplerian solution and Kepler's three laws. Pointwise trajectory
agreement with the analytic conic reaches 3e-8–3e-6 relative over five orbits
(e = 0.3 → 0.85); equal-areas law holds to 7.4e-5; T²/a³ = 4π²/μ across a
20-cell (a, e) sweep to 8.2e-9; specific energy and angular momentum drift
≤ 1.2e-9 over ten orbits. A real-units Earth orbit with IAU 2015/2012
constants reproduces T = 365.256898 d (measured == theory to 2.5e-10), which
sits +1.47e-6 above the sidereal year — the physical two-body residual
(Earth–Moon system mass).

## Content

### The closed-form machinery

For μ = GM, the two-body solution from periapsis at t = 0:

- Kepler's equation: M = n t = E − e sin E, n = √(μ/a³) — solves by Newton
  from E₀ = M + e sin M (≈6 iterations at e ≤ 0.85, tol 1e-14; non-convergence
  raises).
- Position: x = a(cos E − e), y = a√(1−e²) sin E; velocity:
  v = (∂/∂t) of position with dE/dt = n/(1 − e cos E).
- Elements from any state: ε = ½v² − μ/r = −μ/(2a); h = |r × v| = √(μa(1−e²));
  e_vec = (v × h_vec)/μ − r̂ (points to periapsis); p = h²/μ.
- Kepler III: T = 2π√(a³/μ) ⇒ T²/a³ = 4π²/μ ≈ 39.4784176 (μ = 1).

### The key measurement lessons

1. **Uniform-time RK4 must resolve the periapsis passage.** At e = 0.85, 512
   uniform steps/orbit produce ~100% position error after 5 orbits; the
   passage occupies a fraction ~(1−e)^{3/2} of the period, so scale
   steps/orbit ∝ (1−e)^{−3/2} (875 @ 0.3, 2024 @ 0.6, 8814 @ 0.85 → errors
   2.6e-8, 1.8e-7, 3.0e-6).
2. **Sector-area sums are O(h³)-accurate per step.** Triangle area
   ½|r_i × r_{i+1}| has no O(h²) correction because the O(h²) displacement
   term is ∥ r and r × a = 0 for a central force. The same trick (time t vs
   swept area) gives Kepler II from the polygonal trajectory directly.
3. **RK4 is not symplectic** — h-drift 7e-11, ε-drift 1.2e-9 over 10 orbits
   per this problem. Small by discretization order, but not round-off-level.
   (Contrast with Exp 001: symplectic Verlet keeps energy bounded indefinitely.)
4. **Reality anchor arithmetic**: with au = 149597870.7 km exactly and
   GM☉ = 1.3271244e11 km³/s² (both exact IAU-defined values), the two-body
   period at 1 au is 365.256898 d. The Astronomical Almanac sidereal year is
   365.256363 d. The +46 s gap (1.47e-6) is real physics — Earth+Moon mass at
   the Sun, mean-vs-nominal a — NOT numerical error; measured-vs-closed-form
   is 2.5e-10. This makes "1 au two-body orbit" a clean anchor for future
   propagation work (Exp 013 Horizons).

## Source Experiments

- `research/orbital-mechanics/experiments/keplerOrbitValidation/` — full card,
  code, 20 tests, results.json, figures. Runnable:
  `uv run python research/orbital-mechanics/experiments/keplerOrbitValidation/experiment.py`
- Independent adversarial review (2026-08-13): all closed forms and measured
  numbers reproduced with a second integrator family (velocity Verlet), a
  bisection Kepler solver, and an independent periapsis detector — no numeric
  discrepancies found.

## Key Takeaways

- Newtonian central gravity demonstrably produces the Kepler conics; the
  hand-derived closed forms and the integrator agree to 1e-8–1e-6 at
  periapsis-resolved step counts.
- Eccentricity is the resolution driver for time-uniform integration
  (steps ∝ (1−e)^{−3/2}); convergence studies must report it.
- Sector-area (areal-velocity) measurement is the cheapest rigorous Kepler II
  check and double-checks h via dA/dt = h/2.
- IAU-defined constants (au, GM☉) turn a textbook reproduction into a
  quantitative reality anchor: two-body Earth period = 365.256898 d with a
  physically explained 46 s offset from the sidereal year.
- Cross-experiment import hygiene: `tests/__init__.py` + duplicate test
  basenames break the suite once a second experiment exists; tests now load
  `experiment.py` via importlib from explicit paths (see
  `tools/new_experiment.py`).

## See Also

- [[ode-integration-basics]] — integrator order vs geometric fidelity
  (Exp 001), why symplectic methods matter for long horizons.
- [Bate, Mueller, White, *Fundamentals of Astrodynamics* (Dover 1971)]
- [Curtis, *Orbital Mechanics for Engineering Students*, 4th ed. (Elsevier 2021)]
- [Murray & Dermott, *Solar System Dynamics* (Cambridge UP 1999)]
- [IAU 2015 Resolution B3 (nominal GM☉, arXiv:1510.07674)]
- [IAU 2012 Resolution B2 (au = 149597870.7 km)]
- [JPL SSD approximate planetary positions (Standish & Williams 1992)]