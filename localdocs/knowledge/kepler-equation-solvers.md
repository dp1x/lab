---
tags: [orbital-mechanics, kepler, numerical-methods, series, bessel]
date: 2026-08-13
aliases: [kepler-equation-solvers, fourier-bessel-kepler]
links:
  - "[[kepler-orbit-validation]]"
  - "[[ode-integration-basics]]"
---

# Kepler's Equation Solvers: Newton, Bisection, Fourier-Bessel Series

## Summary

Kepler's equation M = E − e sin E was solved three ways plus an independent
reference (fixed point): Newton (order 2 measured, 4–16 iterations to 1e-14),
bisection (halving 0.49992, always 49–50 iterations), and the Fourier-Bessel
series E = M + Σ (2/n) J_n(n e) sin(nM). The series' measured per-term decay
ratio matches Watson's q(e) = e e^χ/(1+χ), χ = √(1−e²), to 0.13% on the
coefficient tail; q → 1 as e → 1, so the truncated series needs ~3400 terms
at e = 0.95 and ~10⁴ at e = 0.99 to reach machine precision — versus 5–9
Newton iterations. All four solvers agree to ≤ 1.8e-13 over a 13-point M grid
at five eccentricities.

## Content

### The Fourier-Bessel identity and its numerical trap

E − M = Σ_{n≥1} (2/n) J_n(n e) sin(nM) — Lagrangian from 1770, formalized
with Bessel functions (Borghi 2024; Philcox-Goodman-Slepian 2021). The
coefficient magnitudes decay geometrically: J_n(ne) ~ (2πχn)^{−1/2} q^n
(Watson §8.4), q(e) = e e^χ/(1+χ) → 1 as e → 1.

**Numerical trap (lesson):** evaluating J_n(ne) by the textbook alternating
power series is hopeless at z = ne ≈ n: the largest term is ~e^{0.46n} while
the result is ~e^{−0.03n} — hundreds of digits of cancellation. Log-space
summation cannot dodge cancellation. The stable route is Miller's backward
recurrence J_{k−1} = (2k/z)J_k − J_{k+1} from a large index, normalized by
the exact identity J_0 + 2Σ_{k≥1} J_{2k} = 1 — only even orders enter
because J_{−k} = (−1)^k J_k (a wrong J_0 + 2Σ all-k normalization yields
values off by exactly Σ_odd J_k, e.g. 1/1.4897 at z = 0.5). Raw sequence
peaks at e^{~1600} for z ~ 600: rescale the 2-term window AND the running
sum in-flight; capture the target r_n only at matching scale (rescaling a
pre-captured value is wrong — it underflows to 0 and must be returned as 0.0,
which is the correct rounding for true values below float64 range).

### Measured results (canonical, tol = 1e-14)

- Newton: order 1.949–2.035 (plateau-corrected local order), worst-case
  iterations 4→16 (e = 0.1→0.99, starter M) and 4→9 (starter M + e sin M) —
  the second starter wins at high e.
- Bisection: exactly 49–50 iterations (theory ceil(log2(π/1e-14)) = 49),
  halving 0.49992, 1 eval/iteration vs Newton's 2.
- q_measured vs q_theory: 0.636314/0.637034 … 0.987939/0.989271 — the
  −0.13% offset is the finite-n prefactor c_n ~ n^{−3/2}q^n, not a theory
  mismatch (mpmath independently reproduced −0.16%).
- Work comparison: series needs N ~ (1/q)·log(1e-14)-ish terms: ~3400 at
  e = 0.95, ~10⁴ at 0.99 — Newton for the win.
- Cross-checks: |bisection − Newton| ≤ 1.6e-14, |series − Newton| ≤ 3.3e-14,
  |fixed-point − Newton| ≤ 1.8e-13; E(0)=0, E(π)=π, E(2π)=2π exact; fixed
  point's linear rate is e·cos E* (0.0363 measured vs 0.0354 expected at
  e = 0.5, M = 1), not e.

### Operational rule for propagation work

For any code path needing E(M, e) many times (Exp 002, Exp 004+): use Newton
with starter M + e·sin M (e < 0.99: ≤ 8 iterations; at worst 9 for
e = 0.99). Fall back to bisection only for adversarial e/M or when a strictly
bracketed guarantee is required; the nested-function design of
`solve_newton`/`solve_bisection`/`solve_series` is reusable as-is.

## Source Experiments

- `research/orbital-mechanics/experiments/keplerEquationSolvers/` — full
  card, code, 46 tests, results.json, 3 figures. Runnable:
  `uv run python research/orbital-mechanics/experiments/keplerEquationSolvers/experiment.py`
- Independent adversarial review (2026-08-13): mpmath 60-digit re-derivation
  of Miller, the series identity, and the Debye asymptotics; J_n matched to
  ≤ 6.9e-15; series-vs-Newton ≤ 8e-15 on 25 random points; two fresh runs
  bit-identical.

## Key Takeaways

- Kepler's equation has no cheap closed form, but Newton (msin starter) is
  effectively it: ≤ 9 iterations to 1e-14 anywhere in e < 0.99.
- The Fourier-Bessel series is beautiful theory with terrible cost near the
  parabolic limit: cost ~ 1/√(1−e²)-per-digit; q(e) quantifies it exactly.
  Use it where uniform analytic structure matters (e.g., series reversion),
  not as a solver.
- J_n(z) at z ≈ n is a classic cancellation trap; the backward recurrence with
  an even-orders-only normalization identity is the robust float64 route —
  and the normalization identity (J_{−k} = (−1)^k J_k) must be gotten right.
- Convergence-order diagnostics must skip the round-off plateau: the last
  residual step is ~1e-16 noise, not part of the geometric tail.
- Fixed-point iteration converges linearly with rate e·cos E* ≤ e; useful as
  an independent cross-check, too slow as a solver at high e.

## See Also

- [[kepler-orbit-validation]] — Exp 002 used Newton-on-Kepler for the
  closed-form reference; this note explains the solver internals.
- [[ode-integration-basics]] — Exp 001; numerical-fidelity context for why
  solver accuracy to 1e-14 matters below integrator error.
- [Borghi, *Mathematics* 12(1):154, 2024 (arXiv:2312.01437)]
- [Philcox, Goodman & Slepian, *MNRAS* 506, 6111 (2021, arXiv:2103.15829)]
- [Watson, *Treatise on the Theory of Bessel Functions*, 2nd ed., §8.4]
- [Gil, Segura & Temme, *Numerical Methods for Special Functions*, SIAM 2007]