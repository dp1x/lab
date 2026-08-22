---
tags: [orbital-mechanics, cr3bp, lagrange-points, rotating-frame, jacobi-constant, stability, validation]
date: 2026-08-22
aliases: [lagrange-points, cr3bp, restricted-three-body, libration-points]
links:
  - "[[j2-precession]]"
  - "[[orbit-decay]]"
  - "[[ground-tracks]]"
  - "[[ode-integration-basics]]"
---

# CR3BP / Lagrange Points: Entering Rotating-Frame Dynamics

## Summary

Exp 011 transitions the lab from single-primary dynamics to the rotating frame.
The circular restricted three-body problem (CR3BP) freezes two primaries on the
x-axis of a frame rotating at n = sqrt(G(m1+m2)/a^3): effective potential
Omega = (1-mu)/r1 + mu/r2 + (x^2+y^2)/2 with EOM r'' + 2 zhat x r' = grad(Omega),
Jacobi C = 2*Omega - v_rot^2 conserved. The five equilibria: L4/L5 exact
equilateral points (1/2-mu, +-sqrt(3)/2); collinear L1/L2/L3 from a scalar axis
equation whose derivative is >= 1 everywhere (exactly three roots, proven
brackets). Collinear points are unstable for EVERY mass ratio; L4/L5 stable iff
27*mu*(1-mu) < 1 (Routh threshold mu_R = (9-sqrt69)/18 = 0.0385208965).

## Content

### Convention traps that survive into every later experiment

* **Coriolis sign is trajectory-only observable.** Flipping the Coriolis term is
  exact time reversal: equilibria identical, characteristic polynomial even in the
  coupling => eigenvalue spectra IDENTICAL to solver noise (1.4e-14), and Jacobi
  conservation unaffected. Only an inertial-consistency round trip (propagate in
  rotating frame -> map -> compare against Newton's law with moving primaries)
  kills it (clean 2.9e-13 vs mutant O(10^-1+)). Any future rotating-frame code
  MUST carry this discriminator; spectra/energy checks give false confidence.
* **Inertial mechanical energy is NOT conserved** in CR3BP — moving primaries do
  net work at rate dE/dt = mu(1-mu) y_R (r1^-3 - r2^-3). The correct inertial
  bridge is C = 2(n h_z - E_I) <=> E_I + C/2 = L_z (kinematic identity, exact).
  A "conserved E_I" in a CR3BP code means primaries were frozen by mistake.
* **Primary distances in inertial-frame formulas need MOVING primaries**
  p_i(th) = R(th) p_i0 — using the static rotating-axis offsets silently breaks
  identities away from th=0 (this bit us during implementation).
* **Inertial validators must advance stage-time-dependent terms through RK4
  stages**: freezing moving primaries inside stages secretly reduces RK4 to
  order 1.
* **Velocity-flip retrace does not work in the rotating frame** (Coriolis is odd
  in v). Time reversal = negative-step flow integration, or compose with the
  mirror M = flip(y, vx, vz) which satisfies the anti-equivariance f(Mw) = -M f(w)
  exactly (bitwise).
* **mu conventions differ**: m2/(m1+m2) vs m2/m1 (ratio ~1.23% for Earth-Moon);
  percent-level mission anchors CANNOT catch the wrong one (gamma shifts only
  0.4%) — pin mu by bit-exact re-derivation from declared GMs.

### Numbers worth keeping

* Earth-Moon (GM-derived mu = 0.012150584078): x_L1 = 0.836915133309,
  x_L2 = 1.155682159554, x_L3 = -1.005062645172; gamma1*a = 58019 km,
  gamma2*a = 64515 km (published ~58k/64.5k); SEL1/SEL2 offsets 1.4976e6 /
  1.5077e6 km (the "~1.5 million km" JWST/SOHO region).
* Asymptotics: gamma_L1,L2 = alpha(1 -/+ alpha/3 - alpha^2/9), alpha=(mu/3)^(1/3);
  x_L3 = -(1 + 5mu/12); gamma3 = 1 - 7mu/12 (deficit from unity!).
* Eigenvalues at collinear points from single quantity A = (1-mu)/r1^3 + mu/r2^3:
  sigma = sqrt(((A-2)+sqrt(A(9A-8)))/2) — discriminant A(9A-8), NOT 9A^2 (a
  seducing factorization slip (1+2A)(1-A) = 1+A-2A^2 produces wrong "clean" forms).
  EM values: sigma(L1)=2.93205593, sigma(L2)=2.15867432, sigma(L3)=0.17787536.
* Triangular: lambda^4 + lambda^2 + (27/4)mu(1-mu) = 0; nu_short -> 1 as mu->0
  (NOT sqrt(27/4)); vertical omega_z = 1 EXACTLY at L4/L5 (frame frequency).
* At mu_R the spectrum is defective (+-i/sqrt2 double root): float64 splits it by
  ~7e-8 (sqrt(eps)-scaling) — expected, do not "fix".
* Nonlinear perturbation doctrine: growth rate measured on the eigenvector-
  projected amplitude over window [2eps, 100eps]; bias ~ 0.69*c_top scales
  linearly with eps (ratios 100.5-101.1 measured across eps ladder). Normal-mode
  ICs: u0 = X + amp*Re(w) with w the FIRST-ORDER eigenvector — velocity slots
  already carry the nu-scaling; scaling them again double-counts and contaminates
  frequency fits by ~1000x.

### Zero-velocity topology is quantitative, not decorative

Neck openings at critical C values are testable: at C = C2 + delta the axis wall
sits strictly left of L2 (bisection on 2*Omega(x,0) = C), particles released at
rest can never pass; at C = C2 - delta escape happens — but through an OPEN neck
escape is phase-dependent: temporary capture held speeds up to 1.05*v_crit bound
past 60 time units. Treat v_crit as necessary condition only; guard every ZVC
classification run with max|dC| <= 1e-6 because grazing flybys corrupt C at
percent levels under fixed-step RK4.

### Shared machinery graduated

`src/lab_utils/integrators.py` (generic rk4_step/rk4_propagate, f(t,x) signature
for non-autonomous RHS) and `src/lab_utils/orbits.py` (element/Kepler canon with
2-5 existing consumers each) now exist; equivalence pinned against donor
experiments. Exp 011 consumes integrators immediately; CR3BP-specific code stays
experiment-local until a second consumer appears (Exp 014+).

### Validation ladder that worked

Scalar root-finding cross-checked against vector residuals AND two algebraically
independent quintic families AND mpmath 40-dps anchors; eigenvalues against closed
forms; C against the bridge identity; frames against inertial consistency; units
against a full dimensional twin pipeline; six adversarial mutants against their
designated killers. Every asserted bound was recorded in results["tolerances"].

## Status

Experiment 011 COMPLETE (2026-08-22). 46 tests green (L1..L7 banners); full-suite
regression maintained; results bit-reproducible across runs.

## Next

Halo/Lyapunov orbit families around L1/L2 (periodic-orbit continuation) would be
the natural 011b; roadmap continues with 012 orbit classes, then 013 JPL
ephemeris validation where Lagrange anchors become sanity checks.
