# Experiment Card: 011 — Lagrange Points / Restricted Three-Body Dynamics

> Status: complete
> Date: 2026-08-22
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/lagrangePoints/`

## Research Question

Does an independently derived, rotating-frame CR3BP model rediscover the five
classical Lagrange equilibrium positions, conserve the Jacobi integral under RK4
propagation at the documented convergence order, reproduce the linear stability
classification (collinear unstable for every mass ratio; triangular stable iff
27μ(1−μ) < 1) together with its nonlinear perturbation signatures, and survive a
dimensional/nondimensional cross-check plus a pre-registered adversarial mutant
battery — while transitioning the laboratory from single-primary dynamics into
rotating-frame multi-body machinery?

This is the sequence's first rotating-frame / non-inertial experiment: every sign,
convention, and oracle was re-derived rather than inherited.

## Background Theory

**Frozen contract v1.0.** Nondimensional barycentric rotating frame; mass ratio
μ = m₂/(m₁+m₂) ∈ (0, ½]; +x from barycenter toward m₂; ω = +ẑ (prograde);
primaries fixed at (−μ, 0, 0) and (1−μ, 0, 0); length unit = separation a,
time unit = 1/n (n = √(G(m₁+m₂)/a³) ≡ 1), mass unit = m₁+m₂.

Effective potential (positive-style):

```
Ω(x,y,z) = (1−μ)/r₁ + μ/r₂ + (x²+y²)/2 ,  r₁ to m₁ (larger primary)
```

Equations of motion: `x″ − 2y′ = Ω_x`, `y″ + 2x′ = Ω_y`, `z″ = Ω_z`
(vector form r″ + 2ẑ×r′ = ∇Ω). The Coriolis term −2ω×v does no work; the
centrifugal +ω²r⊥ is absorbed in Ω.

Jacobi integral `C = 2Ω − |v_rot|²` is conserved along solutions. Exact scalings:
`C_dim = n²a²·C_nondim` (homogeneous), and the inertial bridge
`C = 2(n·h_z − E_I)` ⟺ `E_I + C/2 = L_z` (machine-identity). E_I alone is NOT
conserved — the moving primaries do net work at rate
dE_I/dt = μ(1−μ)·y_R·(r₁⁻³ − r₂⁻³).

Equilibria (∇Ω = 0): L4/L5 = (½−μ, ±√3/2) exactly (equilateral triangles;
completeness: exactly five equilibria exist). Collinear L1/L2/L3 are roots of the
axis scalar f(x) = x − (1−μ)(x+μ)/|x+μ|³ − μ(x−1+μ)/|x−1+μ|³, whose derivative
f′ = 1 + 2(1−μ)/|x+μ|³ + 2μ/|x−1+μ|³ ≥ 1 proves exactly three roots with proven
brackets [½−μ, 1−μ], [1−μ, 2−μ], [−1−μ, −1].

Stability: planar Jacobian has antisymmetric Coriolis coupling (+2/−2). At
collinear points A := (1−μ)/r₁³ + μ/r₂³ > 1 gives Ω_xx = 1+2A, Ω_yy = 1−A < 0,
Ω_zz = −A; spectrum {±σ, ±iν, ±i√A} with σ = √[((A−2)+√(A(9A−8)))/2] > 0 —
unstable for every μ. Triangular points: λ⁴ + λ² + (27/4)μ(1−μ) = 0; purely
imaginary iff 27μ(1−μ) < 1, i.e. μ < μ_Routh = (9−√69)/18 = 0.0385208965045514…
(Routh criterion); vertical mode ω_z ≡ 1 there. Above threshold the quartet
±α±iβ appears with α = √((√γ−½)/2), β = √((√γ+½)/2).

## References

- V. Szebehely, *Theory of Orbits*, Academic Press, 1967.
- C. Murray & S. Dermott, *Solar System Dynamics*, Cambridge UP, 1999, ch. 3.
- W. Koon, M. Lo, J. Marsden & S. Ross, *Dynamical Systems, the Three-Body Problem
  and Space Mission Design*, 2011.
- E. Routh, Proc. London Math. Soc. 6 (1875) — triangular-point linear criterion.
- G. Gascheau, C. R. Acad. Sci. 16 (1843) — necessary condition 27μ(1−μ) < 1.
- V. Arnold, Russ. Math. Surveys 18(6) (1963); A. Deprit & A. Deprit-Bartholomé,
  Astron. J. 72 (1967) — nonlinear stability of L4/L5 (cited, not re-derived).
- NASA public mission pages (JWST orbit page; SOHO/DSCOVR at SEL1; ARTEMIS at
  Earth–Moon libration points; Jupiter Trojans) — percent-level anchors only.

## Assumptions

- Classical circular CR3BP: primaries on exact circular rails (**idealization**;
  real lunar orbit e = 0.055, solar-perturbed).
- No ephemeris fidelity: anchors validate scale+geometry at percent level only
  (**verified**: computed values land within 0.03% of published distances).
- μ derived from declared GM constants (portfolio canon of Exp 008/009/010);
  literature value 0.012150585609624 differs by 1.26e-7 relative (**documented**,
  one convention pinned so the unit cross-check is internally exact).
- Fixed-step RK4 adequate away from grazing encounters (**verified**, with
  Jacobi-drift guards where trajectories approach the secondaries).

## Methodology

Six-track research panel (theory / equilibria / stability / validation design /
adversarial review / software) produced independent derivations that were then
cross-checked; disagreements were resolved by rederivation (three corrections to
the planning brief itself: the inertial Jacobi form, the collinear discriminant
A(9A−8), and the Routh closed form (9−√69)/18).

Production path: bracketed bisection + safeguarded Newton on the axis scalar
(proven brackets, analytic f′); propagation via the newly graduated generic
`lab_utils.integrators.rk4_propagate`; eigenvalues via `np.linalg.eigvals` on the
4×4 planar Jacobian plus decoupled vertical pair; nonlinear perturbation via
eigenvector-projected amplitudes (growth window [2ε, 100ε]) and normal-mode ICs
u₀ = X_L4 + amp·Re(w).

Cases: Earth–Moon (μ = GM-derived 0.012150584078), Sun–(Earth+Moon)
(3.040423452e-6), Routh boundary (0.0385208965…), above-threshold (0.05),
equal mass (0.5), singular limits (1e-3, 1e-6).

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib (Agg), mpmath (40–50 dps anchors)
- Runtime: `.venv\Scripts\python.exe experiment.py` (~100 s; deterministic)
- Determinism: no RNG; two full runs produce bit-identical JSON content
  (only the recorded timestamp differs)
- Dependencies: numpy, matplotlib, mpmath, lab_utils (integrators/metrics/results)

## Validation Method

Layered per the anti-shared-algebra doctrine (tests duplicate all theory inline):
- **L1 frames**: round trips, handedness, angular velocity, Coriolis isolation
  ([0, −2u, 0]), centrifugal direction, L4/L5 exactness.
- **L2 equilibria**: vector residual ‖∇Ω(L_i)‖∞ ≤ 7.3e-16 across all cases;
  TWO algebraically independent quintic families as polynomial oracles;
  mpmath 40-dps anchors agree to 0.0; ordering + bracket containment +
  γ₂ > γ₁; MU-firewall (bit-exact GM derivation, L4-x identity).
- **L3 propagation vs reference**: inertial-consistency round trip (rotating-frame
  propagation mapped to the inertial frame satisfies Newton's law with moving
  primaries) — clean residual 2.9e-13; kinematic closed-form round trip at
  measured order 4.00; mirror field law f(Mw) = −M f(w) exact bitwise.
- **L4 conservation**: Jacobi drift ladders over three trajectory classes — orders
  4.71–4.95 above floor, quantization floor plateau (1.33e-15, rung-independent);
  bridge identity ≤ 1.33e-15; E_I+C/2 = L_z ≤ 5.6e-16; drift-law mismatch
  9.2e-6 → 5.7e-7 under dt/4 (ratio 16.0 ≈ 4⁴); spatial-vs-planar-frozen
  evaluator discrimination.
- **L5 stability**: closed-form vs numeric rates ≤ 9.1e-16 relative; Routh
  criterion across the threshold grid; boundary degeneracy demonstrated in 50-dps
  arithmetic; nonlinear growth-rate recovery (rel err 2.9e-3 at ε=1e-4, bias
  ∝ ε with ratio 100.5–101.1) and long-period frequency recovery (7.5e-7 at
  amp=1e-4 after modal IC construction).
- **L6 units/determinism**: dimensional↔nondimensional equilibria ≤ 5.9e-11 km,
  C-scaling ≤ 1.5e-16 rel, 90-day trajectory correspondence ≤ 1.5e-14 rel,
  Jacobi-along-path ≤ 3.0e-16 rel; double-run bit-identity.
- **L7 adversarial mutants**: coriolis-flip, centrifugal-drop, mapping-sign-flip,
  μ-convention, body-swap domain rejection, planar-frozen Jacobi — ALL caught by
  their designated discriminators (spectra provably blind to Coriolis flip:
  max eigenvalue shift 1.4e-14 — documented why trajectory-level tests are the
  only reliable killer).
- Full suite: `tests/test_lagrange_points.py`, 46 tests, banners L1..L7.

## Results

Headline values (Earth–Moon unless noted):

| Quantity | Value |
|---|---|
| x_L1 / x_L2 / x_L3 | 0.836915133309 / 1.155682159554 / −1.005062645172 |
| γ₁ (EML1, from Moon) | 0.15093·384400 km = **58019 km** (published ~58000) |
| γ₂ (EML2, from Moon) | 0.16783·384400 km = **64515 km** (published ~64500) |
| γ₃ (EML3, from Earth) | **381675 km** (published ~381700) |
| SEL1 / SEL2 offset | 1.49762e6 / 1.50768e3·10³ km (published ~1.5e6 km) |
| C₁/C₂/C₃/C₄=C₅ | 3.188341 / 3.172160 / 3.012147 / 2.987997 |
| σ at L1 (EM) | 2.932055933643 (closed form agrees to 9.1e-16 rel) |
| ν_long / ν_short at L4 | 0.298208173056 / 0.954500863170; ω_z = 1 exactly |
| μ_Routh | 0.0385208965045514 (γ(μ_R)−¼ = 4.7e-51 at 50 dps) |
| Quartet at μ=0.05 | ±0.181986 ± 0.730150i (closed form to 12 digits) |
| Max equilibrium residual | 7.2e-16 (all cases, all points) |
| Jacobi drift orders | 4.71–4.95 (capture class); floor plateau 1.33e-15 |
| Inertial-consistency residual | 2.9e-13 clean vs O(10⁻¹⁺) for mutants |

Figures: `results/figures/f1_geometry_em.png` (rotating-frame geometry),
`f2_zvc_levels_em.png` (zero-velocity structure at critical levels),
`f3_eigen_spectrum.png` (stability spectra incl. unstable quartet),
`f4_jacobi_convergence.png` (drift ladders + floor).

## Limitations

- Circular, patchless model: real libration-point orbits require eccentricity and
  solar-perturbation corrections (percent-level anchor tolerance is honest).
- Collinear points are linearly unstable — station keeping and halo/Lyapunov
  orbit families are out of scope (natural follow-up experiment).
- Linear stability ≠ nonlinear stability in general; the L4/L5 nonlinear claim
  relies on cited KAM refinements (Arnold 1963; Deprit & Deprit-Bartholomé 1967).
- Fixed-step RK4 requires drift guards near grazing encounters; escape through an
  open neck is phase-dependent (temporary capture band measured at
  1.02–1.05 × v_crit within a 60-tu horizon), so "critical escape speed" is a
  necessary condition, not a sufficient one.

## Future Improvements

- Halo/Lyapunov periodic-orbit continuation from the L1/L2 collinear points
  (third-order analytical generators + differential correction).
- Elliptical restricted three-body problem (ERTBP) for the Moon's eccentricity.
- Bi-circular four-body (Sun perturbation) for long-term station-keeping budgets.
- Feed equilibrium/eigenvalue machinery into Exp 014+ trajectory optimization.

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
  (this session used the pinned `.venv` directly; results verified bit-reproducible).
