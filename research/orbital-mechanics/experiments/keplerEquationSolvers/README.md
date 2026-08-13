# Kepler's Equation Solvers: Newton, Bisection and the Fourier-Bessel Series

> Status: complete
> Date: 2026-08-13
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/keplerEquationSolvers`

## Research Question

How do the classical solvers for Kepler's equation M = E − e sin E compare —
Newton iteration, bisection, and the Fourier-Bessel (Lagrange/Bessel) series
E = M + Σ (2/n) J_n(n e) sin(nM)? Can the empirically observed series decay
ratio be quantitatively explained by the asymptotic theory q(e) (Watson,
Sec. 8.4), and how pathological is the series as e → 1?

## Background Theory

Kepler's equation M = E − e sin E (M mean anomaly, E eccentric anomaly,
0 ≤ e < 1) is the universal time-to-position gate of astrodynamics (Exp 002
used a Newton solver for it). It has no closed-form solution in elementary
functions; three classical families exist:

- **Newton iteration** on f(E) = E − e sin E − M: f is strictly increasing
  (f′ = 1 − e cos E > 0 for e < 1), so convergence from any starter on either
  side of the root; expected quadratic order p = 2.
- **Bisection** on the canonical bracket [M, π] (M ≤ π), extended to
  M ∈ (π, 2π) by E(2π − M) = 2π − E(M); exactly one function evaluation per
  iteration, guaranteed linear at rate 1/2.
- **Fourier-Bessel series** (Lagrange 1770; Bessel; see Borghi 2024;
  Philcox, Goodman & Slepian 2021):

      E − M = Σ_{n≥1} (2/n) J_n(n e) sin(nM)

  The coefficients decay geometrically with ratio (Watson 1944, Sec. 8.4;
  from J_n(n e) ~ (2πχn)^{−1/2} q^n, χ = √(1−e²)):

      q(e) = e·exp(√(1−e²)) / (1 + √(1−e²))  →  1 as e → 1

  so the series needs O(1/(1−q)) ≈ O(1/√(1−e²)) terms to reach a given
  accuracy — catastrophically slow near the parabolic limit. Evaluating the
  coefficients needs J_n(ne) at z ≈ n, where the alternating power series
  suffers catastrophic cancellation (largest term ~ e^{0.46 n} vs result
  ~ e^{−0.03 n}); instead the coefficients are evaluated with Miller's
  backward recurrence, which is stable for z < n (DLMF 10.74(iii)).

**Fixed-point iteration** E ← M + e sin E serves as a fourth, independent
reference solver (linear convergence, asymptotic rate e·cos E* ≤ e).

## References

- P. Colwell, *Solving Kepler's Equation over Three Centuries*, Willmann-Bell,
  1993 — history and methods.
- R. Borghi, "The Kepler equation: Fourier–Bessel expansions", *Mathematics*
  12(1):154, 2024 (arXiv:2312.01437) — the (2/n)J_n(ne) coefficient identity
  and its convergence.
- O. H. E. Philcox, J. Goodman, Z. Slepian, "Kepler's Goat Herd: An update to
  the series expansion method for Kepler's equation", *MNRAS* 506, 6111–6116,
  2021 (arXiv:2103.15829) — modern treatment of the series method.
- G. N. Watson, *A Treatise on the Theory of Bessel Functions*, 2nd ed.,
  Cambridge UP, 1944, Sec. 8.4 — J_n(ne) ~ (2πχn)^{−1/2} q^n asymptotics.
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed.,
  Elsevier, 2021, Ch. 3 — Newton and starter choices.
- DLMF §10.74(iii), Gil, Segura & Temme, *Numerical Methods for Special
  Functions*, SIAM, 2007 — stable backward-recurrence evaluation of J_n.

## Assumptions

- Kepler's equation on M ∈ [0, 2π) with 0 ≤ e < 1; all solvers target
  residual |f(E)| < 1e-14 (verified).
- Fourier-Bessel identities and q(e) asymptotics taken from the cited
  literature; verified numerically in this experiment (verified).
- Coefficients evaluated in float64 via Miller backward recurrence, judged
  against DLMF/Abramowitz-Stegun published J_n values (verified).
- Fixed-point reference solver included purely as an independent cross-check,
  not as a practical solver (plausible).

## Methodology

All studies in `experiment.py`, deterministic (no RNG):

1. **Newton order**: residual history per (e, starter) at M = 1; convergence
   order read from the residual steps clear of the round-off plateau
   (p = log(r_{k+1}/r_k)/log(r_k/r_{k−1})); starters M and M + e sin M.
2. **Bisection**: bracket-width history, measured halving factor (~0.5) and
   iteration counts; analytic count ceil(log2(π/1e-14)) = 48.4 ≈ 49.
3. **Eccentricity sweep**: worst-case iterations + function evaluations over
   the 13-point M grid for e ∈ {0.10 … 0.99}.
4. **Series decay**: max residual over the M grid for N ∈ {2 … 2048} terms at
   e ∈ {0.5, 0.7, 0.85, 0.9, 0.95}; measured decay ratio from the geometric
   mean of |c_{n+1}/c_n| over the coefficient tail (n ≥ 512) vs q_theory(e).
5. **Cross-solver agreement**: max |ΔE| over the M grid at
   e ∈ {0.1, 0.3, 0.6, 0.85, 0.95}, Newton as common reference.
6. **Special values**: E(0) = 0, E(π) = π, E(2π) = 2π for all e.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib (no scipy)
- Runtime: `uv run python experiment.py`
- Determinism: no RNG; J_n(ne) by Miller backward recurrence (exact
  deterministic arithmetic, in-flight rescaling); subprocess cross-process
  bit-identical check in the test suite.
- Dependencies: numpy, matplotlib (already in pyproject; no new deps).
- Reuses `lab_utils.results` (save_json_result).

## Validation Method

Independent checks in `tests/test_kepler_equation_solvers.py` (46 tests,
all green):

- Closed-form anchors: E(0) = 0, E(π) = π, E(2π) = 2π; all solvers reach
  |residual| < 1e-13 across 5 e × 6 M (parametrized).
- Newton order ~2 (1.95–2.04 measured, plateau-corrected); bisection
  halving = 0.5 ± 0.01 and 44–60 iterations; fixed-point rate = e·cos E*
  (e.g. 0.0363 measured vs 0.0354 expected at e = 0.5, M = 1).
- Series coefficients against *published* constants: c₁ = 2J₁(0.5) =
  0.4845369153, c₂ = J₂(1.0) = 0.1149034849, J₀(1) = 0.7651976866
  (DLMF 10.2), rel 1e-9; series vs Newton < 1e-9 at N = 2048 for
  e ∈ {0.3 … 0.9}; measured q vs q_theory within 0.15% at e ∈ {0.5, 0.7, 0.85};
  q(e) monotone increasing, q(0.99) > 0.998; residual at N = 64 is still
  ~1e-3 at e = 0.9 (genuinely slow), < 1e-9 at N = 2048.
- Cross-solver agreement over the full grid < 1e-12 (Newton vs bisection),
  < 1e-11 (fixed point), < 1e-7 (series at e = 0.95).
- Determinism: bit-identical JSON in a fresh interpreter.

An independent adversarial review (subagent, 2026-08-13) re-derived the
methods with mpmath at 60-digit precision: jn_miller matched independent J_n
to ≤ 6.9e-15 on every (n, z) the studies use; the Fourier-Bessel identity
held to 0.0 at N ≤ 40 and 8e-15 at N = 2048 over 25 random (e, M); the
Debye asymptotics J_n(ne)√(2πnχ)/q^n → 1 confirmed the q(e) formula; two
fresh runs were bit-identical. Three doc items and one latent out-of-range
edge (post-capture rescaling of the coefficient underflow) were fixed; no
published number changed.

## Results

All results in `results/results.json`; figures in `results/figures/`.

**Convergence orders/rates** (M = 1.0 rad):

| solver | e | iterations | measured order/rate |
|--------|---|-----------|----------------------|
| Newton (msin) | 0.3 | 4 | 2.00 |
| Newton (msin) | 0.6 | 5 | 2.00 |
| Newton (msin) | 0.9 | 5 | 2.00 |
| Newton (m)    | 0.9 | 7 | 1.97 |
| Bisection     | 0.3…0.9 | 49 | halving 0.49992 |
| Fixed point   | 0.5 (M=1) | ~17 | 0.0363 vs e cos E* = 0.0354 |

**Worst-case iterations over the M grid** (e → 1): Newton-M 4 → 16,
Newton-(M + e sin M) 4 → 9, bisection 50, as e: 0.1 → 0.99. Newton's
starter M + e sin M is a clear win at high e (16 → 9 iterations at e = 0.99).

**Series decay ratio vs Watson's theory**:

| e | q_measured | q_theory | rel diff |
|---|-----------|----------|----------|
| 0.5  | 0.636314 | 0.637034 | −1.1e-3 |
| 0.7  | 0.832936 | 0.834064 | −1.4e-3 |
| 0.85 | 0.941528 | 0.942802 | −1.4e-3 |
| 0.9  | 0.967919 | 0.969228 | −1.4e-3 |
| 0.95 | 0.987939 | 0.989271 | −1.3e-3 |

The systematic −0.13% offset is the finite-n prefactor c_n ~ n^{−3/2} q^n
(ratio q·(1−3/(2n)) at n ≈ 1000), independently reproduced with mpmath at
−0.16%. The q → 1 behavior quantifies the pathology: at e = 0.95 the series
still needs ~3400 terms for 1e-14, at e = 0.99 ~10⁴ terms — versus 5–9
Newton iterations.

**Cross-solver agreement** (M grid; Newton = reference): |bisection −
Newton| ≤ 1.6e-14, |series(N=2048) − Newton| ≤ 3.3e-14, |fixed-point −
Newton| ≤ 1.8e-13 across e ∈ {0.1 … 0.95}. Special values exact.

## Limitations

- Series method tested to N = 2048 and e ≤ 0.95; at e ≥ 0.97 the term count
  for machine precision exceeds the study's cap. The series is interesting
  mathematically and historically, but Newton with the M + e sin M starter
  dominates it operationally at every e tested (5–16 iterations vs 10²–10⁴
  terms).
- J_n(ne) coefficients underflow float64 for e ≲ 0.35 and n ≳ 350 (true
  values < 1e-308); handled by rounding to 0 — harmless to every sum here,
  but the coefficient array itself should not be read as exact below that
  scale.
- q_measured carries the systematic finite-n bias of the n^{−3/2} prefactor
  (−0.13% at n ≥ 512; the residual-based alternative fit is noisier because
  max-over-M residual decay oscillates with N).
- Fixed-point solver is a reference only; no acceleration (Aitken) tried.
- Starter study covers M + e sin M and M; Danby-method type starters not
  compared.

## Future Improvements

- Hohmann transfer economics (roadmap 004) will reuse these solvers.
- Compare against a starter-of-choice survey (Danby 1957/1991, Gooding &
  Odell, Markley) on the e → 1 frontier.
- Derive and verify the residual-envelope bound Σ_{n>N} c_n sin(nM) vs q^N
  beyond the coefficient-ratio measurement used here.
- Couple with Exp 002's periapsis-resolved propagation for E-based
  uniform-anomaly sampling (regularization use case).

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
- Numerical results are deterministic; re-running rewrites `results.json` with
  a fresh timestamp and generation commit (provenance), while the scientific
  output stays identical.