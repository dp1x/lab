---
tags: [synthesis, orbital-mechanics, numerics, knowledge-base, closure]
date: 2026-08-17
aliases: [om-001-006-synthesis, lab-synthesis-1]
links:
  - "[[ode-integration-basics]]"
  - "[[kepler-orbit-validation]]"
  - "[[kepler-equation-solvers]]"
  - "[[hohmann-transfer]]"
  - "[[bielliptic-vs-hohmann]]"
  - "[[combined-transfer-plane-change]]"
  - "[[audit-006-closure-2026-08-17]]"
---

# Synthesis: Orbital-Mechanics Experiments 001–006

> **Scope note (authoritative):** This is a **strictly orbital-mechanics
> 001–006 synthesis**. Experiment 001 is the Numerics foundation that *feeds*
> the Orbital Mechanics flagship; 002–006 are the first five experiments of
> the Orbital Mechanics phase (Phase 2 of the roadmap). This boundary exists
> because the lab's **synthesis cadence (~5 experiments, `roadmap.md`) happened
> to trigger here** — it is an **"all validated experiments to date /
> cadence-triggered synthesis boundary," NOT a formally defined scientific
> phase**. The roadmap does not establish 001–006 as a phase, and Phase 2
> itself is still continuing.
>
> **Three concepts kept distinct (permanent lab convention, earned in the
> Exp 006 audit):**
> **experiment sequence ≠ scientific phase ≠ synthesis-cadence boundary.**
>
> No 13-part laboratory timeline is referenced here, because the repository
> contains no evidence for one (see `[[audit-006-closure-2026-08-17]]` / A0
> scope audit).

## 1. The chain and what each experiment contributed

| # | Experiment | Primary contribution | Reused by later experiments |
|---|-----------|----------------------|------------------------------|
| 001 | Numerical integrators (ODE) | Verified RK4 convergence (order 4.01) + symplectic energy behavior; rigorous reference-grid method | All later RK4 propagation; the grid-alignment lesson (use `T/n` steps) |
| 002 | Kepler orbit validation | Verified two-body RK4 propagator + Kepler solver + elements extraction; closed-form anchors (Kepler I/II/III), IAU Earth year | 004, 005 (loaded by path); 006 built a *new* 3D Cowell propagator on the same principles |
| 003 | Kepler's-equation solvers | Newton/bisection/Fourier–Bessel solvers; Watson q(e) asymptotics; Miller recurrence for J_n | Time-to-position in transfers; cross-checked in 002 |
| 004 | Hohmann transfer | Closed-form Δv, R* = 15.5817187388 peak, (R−1)/2 and √2−1 asymptotes, 2-impulse optimality, real anchors (LEO→GEO 3.9319 km/s, E→Mars 258.87 d) | 005, 006 (IAU constants, Hohmann cost surface) |
| 005 | Bi-elliptic vs Hohmann | R_bp = 11.9387654726, corner identity R* = Exp 004 peak (1e-29 @ 50 digit), crossover curve s_c(R), hump-onset resolution | 006 (bi-parabolic limit f_bp(R)) |
| 006 | Combined transfer + plane change | Three regimes (two-burn / finite-s / s→∞); boundaries di_c(R), di_inf(R); s→∞ = bi-parabolic identity (1e-16); R=2/47.5° dip 1.77%; SES-8 5.21% | — (frontier of the closed chain) |

**Reuse discipline honored:** no scaffolding was rebuilt; 004 imports 002's
propagator, 005 imports 004, 006 imports 002/004/005 by explicit-path importlib.
006's 3D Cowell integrator is new code but built on the same first-principles
principles and validated only against its own closed-form burns.

## 2. Knowledge hierarchy (claim classification)

### Established analytical results (derived + closed-form, high confidence)
- **R\*** = 15.5817187388, the Hohmann cost maximum (independently re-derived
  in the 2026-08-17 audit; dv/v1 = 0.536258).
- **R_bp** = 11.9387654726, the bi-parabolic crossover (bisection, 1e-14).
- Hohmann asymptotic behavior: Δv/v1 → (R−1)/2 as R→1; → √2−1 (escape burn)
  as R→∞.
- Bi-elliptic crossover structure: s_c(R) diverges at R_bp, meets the corner
  at R*; exactly one crossing per R ∈ (R_bp, R*).
- Exp 006 three-regime structure (two-burn / finite-s 3-burn / s→∞).
- **s→∞ identity:** the combined-problem asymptotic equals the coplanar
  bi-parabolic limit f_bp(R) = (√2−1)(1 + 1/√R), **independent of Δi** (verified
  to ~1e-13 mpmath, ~1e-16 float64) — the plane change becomes "free" at the
  near-rest apoapsis.
- R=1 detour corners: di_c = 2·arcsin(1/3) = 38.9424°, di_inf = 60°.

### Numerically established results (resolved by computation; confidence tied to method)
- Optimizer maps: di_c(R), di_inf(R) regime boundaries; the finite-s window
  **pinch** where the window closes.
- **Representative minima:** R=2, Δi=47.5° finite-s dip beats two-burn by
  **1.77%** (s* ≈ 2.72, dv ≈ 0.6501); SES-8 super-synchronous (R≈13.7, Δi=30°)
  wins **5.21%** (4.096 vs 4.321 km/s).
- Real-system anchors: LEO→GEO 28.6° (two-burn optimal, 0% over split),
  Curtis 300 km→GEO, GTO→GEO 5° — all two-burn optimal.
- 004 real anchors: LEO→GEO 3.9319 km/s / 5.26 h; Earth→Mars 258.87 d,
  TMI 3.6114 km/s; Venus inward 146.08 d.

### Model-dependent results (idealizations; limits explicitly noted)
- All transfers assume **two-body, circular, coplanar (or common-node plane
  change), impulsive burns, point-mass gravity**; no J2, drag, third bodies,
  finite burn arcs, or phasing.
- Exp 006's 5% super-synchronous saving is a **Δv optimum, not a mission
  Δv** (unbounded time penalty as s→∞).
- Patched-conic-type approximations (used in 005/006 limiting identities) are
  model choices, not exact physics.
- Neglect of ephemeris / real planetary states means anchors are validated
  against *mean* orbits and published worked examples, not live Horizons data
  (that is roadmap experiment 013).

## 3. The boundary-definition principle (from the Exp 006 audit)

When any experiment reports a regime/phase boundary, distinguish and **name**
which of these is being quoted:

1. **Mathematical optimum / root** — the exact crossing of continuous cost
   curves.
2. **Numerically resolved boundary** — what a given optimizer + resolution
   actually resolves.
3. **Operational classification threshold** — a deliberate robustness margin
   (e.g. Exp 006's `WIN_MARGIN = 1e-5`, requiring the 3-burn to beat two-burn
   by ≥ 1e-5 before a regime change is declared).

Concrete example (Exp 006, R=1.05): the **mathematical** di_c crossing is
≈ 11.4–11.9°; the **operational/robust** boundary actually reported is
**17.01°**; `R_pinch = 6.214815` is the **operational** pinch, while the
independent continuous mathematical pinch is ≈ 6.48–6.51. Both soft bands are
legitimate; the reported values are *operational*, not exact roots. Regression
tests assert **tolerance bands**, not false exactness. This distinction
prevented the final audit from mistaking a robustness convention for a physical
law and must carry into all future experiments (e.g. Exp 007's flyby-geometry
optimality boundaries).

## 4. What is NOT yet done (remaining Orbital-Mechanics roadmap)

Per the authoritative `roadmap.md` (followed literally; **not revised**):

| # | Experiment | Note |
|---|-----------|------|
| 007 | Gravity Assist / slingshot | patched-conic / known flyby numbers — **next** |
| 008 | Ground Tracks | spherical geometry, real orbit params |
| 009 | J2 Precession | node drift from Earth's bulge (analytic secular rates) |
| 010 | Orbit Decay | drag over time → re-entry |
| 011 | Lagrange Points | **first restricted-three-body** experiment |
| 012 | Orbit Classes | Sun-synchronous / Molniya / GTO |
| 013 | JPL Ephemeris Validation | **first ephemeris** validation vs NASA Horizons |
| 014+ | Eclipse timing, launch windows, trajectory opt, solar sails | |

The conceptual *model-fidelity progression* 2-body → patched-conic →
restricted 3-body (011) → ephemeris (013) is noted here as a **research
insight**, but it is **not** the experiment numbering or immediate sequence;
the roadmap order above is authoritative.

## 5. Carry-forward rules for the lab
- Reuse verified machinery (002 propagator, 004 constants, 005 bi-parabolic
  limit); never rebuild scaffolding.
- Separate **established analytical / numerically established / model-dependent**
  claims in every card, knowledge note, and synthesis.
- Apply the **mathematical-root vs operational-threshold** distinction to every
  boundary reported.
- Validate before trusting; independent optimizer + mpmath + RK4 cross-checks
  are the established standard (proven effective in 004–006).

## Source authority
All claims above are drawn from the experiment cards/READMEs, knowledge notes,
`results.json` artifacts, the 2026-08-17 adversarial + closure audits, and
`localdocs/roadmap.md` / `charter.md`. No experiment implementation or result
was altered to produce this synthesis. Head of chain at time of writing:
`b864798` (origin/main, clean tree, 195 tests passing).
