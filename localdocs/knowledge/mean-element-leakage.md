---
title: Mean-element leakage and the lunisolar RAAN residual
status: UNRESOLVED — the residual is not window-independent
date: 2026-10-03
mission: mission_mean_element_leakage
supersedes: "[[lunisolar-closure-021]]" residual interpretation (partially)
related: "[[j2-precession-modulated-lunisolar-term]]", "[[lunisolar-secular-limit-020]]", "[[scientific-stack]]"
---

# The SSO lunisolar RAAN residual is a windowing artifact of a long-period term

## Claim

The residual left after the "baseline-breathing" correction at sun-synchronous
inclination, h = 600 km — **−8.564733e-4 deg/day** over 18.6 years — **is not a
window-independent quantity and therefore is not evidence of a dynamical interaction.**

Measured on one propagation, with one code version, varying only the analysis window:

| window | residual (deg/day) |
|---|---|
| 1.00 yr | +1.3929e-4 |
| 4.65 yr | +4.8318e-4 |
| 9.30 yr | +1.6934e-3 |
| 12.00 yr | +1.9164e-3 |
| 18.60 yr (the committed value) | **−8.5792e-4** |

Sliding the start epoch of a 4.65-yr window spans **−6.574e-4 … +5.615e-4**
(peak-to-peak 1.219e-3, one sign change); a 9.30-yr window spans
**−1.899e-3 … +1.693e-3** (peak-to-peak 3.593e-3). The committed value is **smaller than
its own window-to-window standard deviation** in the 4.65-yr ladder (2.09σ) and sits at
0.73σ in the 9.30-yr ladder.

A secular rate cannot reverse sign when the window moves.

## Why

The detrended node-difference series `Ω(sun_moon_j2) − Ω(j2_only)` is **99.4%** explained
by a single sinusoid of ~34.8 deg amplitude at a fitted ~9175-day period
(RMS misfit 18.70 → 0.106 deg). All other periods are ≤ 0.26 deg (annual 0.256,
semiannual 0.093, lunar monthly ≤ 0.024).

The arc spans **6793.6 d = 0.894 lunar nodal cycles**. Over less than one cycle, the OLS
slope of `A·sin(ωt+φ)` is a **phase-dependent** quantity of order `A·ω` — a constant
independent of window length (audit-020 Track 3's slow-harmonic regime). Here
`A·ω/2π ≈ 3.3e-3 deg/day`, roughly **4× the residual**.

Fitting the secular and 6798-day terms jointly is **degenerate** over 0.894 cycles (the
fitted period runs to the search bound). That degeneracy *is* the diagnosis: on this arc
the secular term and the lunar-nodal term cannot be separated, and the least-squares
slope of one absorbs the other.

## What is falsified

| Hypothesis | Verdict | Decisive evidence |
|---|---|---|
| `a` decays 165 km over 18.6 yr | **FALSIFIED — RK4 step artifact** | two-body control (where `a` is exactly conserved) gives p = **5.001/5.000** at dt 60/30/15 s; extrapolated −143.0 km vs −164.8 km committed |
| the `a` drift contaminates the residual | **FALSIFIED — common-mode, cancels** | pure-`a` channel of `D` is 0.034% of `D` at SSO, **0.000%** at 90° where `a` still dropped 86 km (`cos 90° = 0`) |
| mean-element treatment is inadequate | **FALSIFIED** | see below |
| node-crossing sampling biases the rate | **EXONERATED** | node-detection-free estimator reproduces committed `combined` to **5.8e-9 deg/day** at 1 yr, all three inclinations |
| `ols_slope` conditioning artifact (+2.17e-3) | **DOES NOT REPRODUCE** | +4.7e-18 on the real grid; re-scoring moves the residual by 8.9e-16 |

## The mean-element question, settled against the lab's own prior fear

Three independent lines converge:

1. **First-principles.** The J2 short-period terms are
   `Ω_sp = +3J₂(R/a)²cos i·sin u`, `i_sp = −(3/2)J₂(R/a)²sin 2i·cos u`,
   `a_sp = +2J₂(R/a)e·sin u`. At the ascending node **`u ≡ 0`, so `Ω_sp` is exactly
   zero there** — the sample point is the node of the short-period node term's own
   argument. Confirmed numerically at dt = 7.5 s: the detrended short-period `Ω` is
   ≤ 3.2e-4 deg, four orders of magnitude too small to bias an 8.6e-4 deg/day rate.
2. **Aliasing.** A per-revolution harmonic sampled at one point per orbit contributes
   **exactly zero** to an OLS slope. There is no small-divisor channel between
   short-period content and nodal drift in this sampling scheme.
3. **Literature.** Vallado (2013) §9.6 p.654: *"for a first-order theory, it's immaterial
   whether we use mean or osculating elements on the right-hand sides … the resulting
   errors will be on the order of J²."* The lab's canon is first-order, so the
   mean/osculating error is **O(J₂²) by construction** — below its own truncation level.

The one genuine defect of the 30-day boxcar is the **Jensen gap** (averaging `a,i` then
evaluating the nonlinear `Ω̇_J2`), measured from the committed scalars at **0.17% of the
residual** at SSO — three orders of magnitude too small.

### Retracted prior claim

`audit-019-track-F` §4's *"~0.86 deg of short-period scatter at each crossing"* is **not
reproducible** and is retracted. It multiplies the short-period amplitude by
`T_snap/T_orb ≈ 15`, which is invalid. The lab's separate conclusion — that a full Brouwer
transform is the wrong tool here — is **confirmed and strengthened**: at the node there is
literally nothing for a Brouwer short-period subtraction to remove.

### Citation correction

"Brouwer (1958) *The artificial satellite orbits of Earth*" **does not exist**. The
correct source is **Brouwer, D. (1959), "Solution of the Problem of Artificial Satellite
Theory Without Drag," *Astron. J.* 64(1274), 378–397**, DOI 10.1086/107958 (verified via
Crossref and Vallado's reference list). Kozai (1959) is *Astron. J.* 64(1274), 367–377.

## A nuisance term that swamps the 30° column

`D = ⟨Ω̇_J2⟩_full − ⟨Ω̇_J2⟩_j2` differences two *mean analytic* rates, but the chain subtracts
`j2`'s *OLS slope*. The gap is a systematic offset:

| i | contaminant `⟨Ω̇_J2⟩_j2 − rate_j2` | residual | ratio |
|---|---|---|---|
| 97.7876° | +1.059212e-3 | −8.564733e-4 | 1.24× |
| 90.0° | 1.6e-18 | −4.243816e-4 | ~0 |
| 30.0° | **−1.237715e-2** | −1.951493e-4 | **63×** |

At 30° the estimand is 63× below its own nuisance floor: **that column carries no
information about H-baseline at either window.** The fix is to regress against the single
analytic baseline `full.rate − Ω̇_J2(⟨a⟩,⟨i⟩)` instead of differencing two baselines.

## Numerical note that matters for every future arc

`dt = 60 s` gives only **96.7 steps per orbit** and under-resolves the once-per-orbit J2
short-period content of `a` (±9.2 km at SSO). It biases the crossings-sampled `⟨a⟩` by
−13.7 km (−9.9 km of which is aliasing, the rest RK4 drift), which shifts the **absolute**
`Ω̇_J2` level by **0.50%**. That is common-mode and cancels in the difference, so the
residual is safe — but **any absolute rate claim needs dt ≤ 15 s**.

## Open

1. **The identity of the ~9000-day period.** Fitted 9175 d, +35% vs the 6798.4-d lunar
   nodal regression, degenerate with the linear term. Needs > 1 lunar nodal cycle.
2. **Whether any secular term survives.** Not measured.
3. **Multi-cycle DE441 span.** > 37 yr, or an explicit two-phase-locked window design.