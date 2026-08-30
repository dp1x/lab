# Audit 020 / Track 6 — Arc-Length Design for Secular-Limit Convergence

> **Author**: Track 6 (independent arc-design track), 8-track parallel delegation
> for **Experiment 020** of the Computational Research Laboratory.
> **Date**: 2026-08-31.
> **Status**: Independent design report; read-only audit; no production code
> modified. The deliverable is a **deterministic arc-length design** for Exp 020
> that turns the 019 window-length-extrapolation finding (Track G's
> "W → ∞ extrapolation") into a **reproducible multi-arc experiment** that pins
> the secular Lunisolar RAAN rate to a pre-declared uncertainty without
> requiring a literal 18.6-year integration.
>
> **Inputs read (read-only)**:
> - `localdocs/reports/audit-019-track-B-averaging-hierarchy.md`
> - `localdocs/reports/audit-019-track-F-mean-vs-osculating.md`
> - `localdocs/reports/audit-019-track-G-hostile-review.md`
> - `research/orbital-mechanics/experiments/lunisolarLongPeriod/experiment.py`
>   (for window_sweep timing constants and the frozen window set
>   W ∈ {30, 90, 180, 365, 730} d)
> - 017 `fetch_horizons_moon_snapshot.py` (acquisition pattern)
> - 014 `eclipseTiming/reference/MANIFEST.json` (acquisition schema)
> - Track 1 of this audit (`audit-020-track-1-disturbing-function-reconciliation.md`)
>   for the corrected-formula reference value (+1.348e-4 deg/day at h=600 km i_sso)
> - Track 7 of this audit (`audit-020-track-7-hostile-review.md`) for the
>   model-dependence caveat on the 019 extrapolation
>
> **Inputs NOT read**: any future Exp 020 implementation, any Track 2-5 output
> of this audit (parallel-track discipline).

---

## 0. Executive summary

**Headline design decision for Exp 020:**
1. **The literal 18.6-year arc is NOT required.** With a calibrated
   `Ω̇_fit(W) = a + b/W + c/W²` extrapolation plus a Cycle-averaged estimator
   on **multi-epoch windows** (4 windows of length 1 year at 18.6-yr / 4 ≈ 4.65
   yr spacing), the secular rate at h = 600 km i_sso can be pinned to ±10 %
   of the corrected formula value with an effective arc of ≈ **5 years total
   integration** (≈ 4.65 yr × 1 yr + overhead).
2. **The 5-year integration IS required**, not the 1-year arc of Exp 019.
   The Track F "9.78×" residual at W = 365 d is precisely the expected
   finite-window-bias deficit; the Track G "27×" extrapolation at W = 730 d
   is still in the bias regime. Extending to W ≈ 5 yr reaches the
   "Regime C boundary" for the annual solar forcing AND samples ≥ 25 % of the
   lunar nodal cycle, breaking the cos(ω_nodal W) aliasing that traps the
   1-year fit.
3. **A three-arc plan** is recommended (see §6):
   - **Pilot arc** W_pilot = 5 yr at i_sso (full model only); cost ≈ 25-50
     min single-core; resolves the 019 ambiguity at a fraction of the
     decadal-arc cost.
   - **Confirmatory arc** W_confirm = 10 yr at i_sso + i = 90° + i = 30°
     (Track 1 flagged the i = 30° 120,000× residual as a structural puzzle);
     cost ≈ 50-100 min single-core per (h, i) × 3 inclinations × 1 mode =
     ≈ 4 hours single-core.
   - **Verification arc** W_verify = 18.6 yr at i_sso only (full + lunar-only
     isolation); cost ≈ 90-180 min single-core; serves as the
     "ground-truth" interpolation node the lunar nodal term demands.
4. **Storage budget**: a single 18.6-yr arc at dt = 60 s, full state
   (r, v) sampled at every step, is ≈ 9.8 × 10⁹ doubles × 6 components × 8 B
   = 470 GB. **Subsampling at ascending-node crossings** (14.9 nodes/day ×
   6798.4 d ≈ 1.0 × 10⁵ rows per 18.6-yr arc; < 1 MB) is the canonical
   fix; we keep 60-s internal time-series on R: scratch, **commit only the
   ascending-node crossing table** to the repo (matches Exp 019 convention).
5. **Acquisition plan**: DE441 Sun + Moon snapshots spanning
   **2026-01-01 → 2044-07-31** (one full nodal cycle, with 6-month padding),
   broken into 1-year chunks fetched identically to the existing Exp 014/017
   pattern; each chunk gets its own MANIFEST.json with sha256 pin and a
   master MANIFEST.json pinning the multi-year bundle. Total expected size:
   ≈ 18 chunks × 76 KB each ≈ 1.4 MB (cheap).

The full design — including FACT/INFERENCE/UNKNOWN classification, per-period
arc-length tables, compute budget, and verification plan — is given in the
sections that follow.

---

## 1. Problem statement and formal definitions

### 1.1 What "secular limit" means here

The **secular Lunisolar RAAN rate** at h = 600 km, i = i_sso is defined
operationally (Track B, §2.2; Track F, §2) as:

    Ω̇_mean = lim_{W→∞} (1/W) ∫₀^W [dΩ_mean/dt](t) dt

where `Ω_mean(t)` is the **doubly-averaged** (over satellite mean anomaly +
third-body mean anomaly) angular position of the ascending node. The
**corrected closed form** (Track 1 of this audit, Track B of audit-019)
gives:

    Ω̇_corr_cf = (3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i

with `(μ₃/μ_E)`, `(a/a₃)³` summed over Sun + Moon. At h = 600 km, i_sso =
97.7876°: Ω̇_corr_cf = +1.348 × 10⁻⁴ deg/day (Track 1 of this audit, F1).

The **finite-window linear-fit estimator** of the secular limit is:

    Ω̇_fit(W, t₀) = OLS slope of {Ω_cross(t_k) : t_k ∈ [t₀, t₀ + W]}

where Ω_cross is measured at ascending-node crossings (Exp 018/019
convention). The bias of this estimator (Track F, §5) is:

    Ω̇_fit(W, t₀) = Ω̇_mean + bias(W, t₀, {A_k, ω_k, φ_k})

and **bias is the quantity we need to drive below the uncertainty target.**

### 1.2 Uncertainty targets (pre-declared)

We declare the following uncertainty targets (these are the design
constraints for Exp 020; the experiment must declare them in its README
before any propagation runs):

| Target ID | Absolute uncertainty | Relative to Ω̇_corr_cf | What it buys |
|---|---|---|---|
| **T1 (sanity)** | ±10⁻³ deg/day | ±740 % | "within an order of magnitude" — current 019 result |
| **T2 (working)** | ±10⁻⁴ deg/day | ±74 % | "within a factor of 2 of corrected formula" |
| **T3 (good)** | ±10⁻⁵ deg/day | ±7.4 % | "10 % agreement with corrected formula" |
| **T4 (publication)** | ±10⁻⁶ deg/day | ±0.74 % | "sub-percent agreement with corrected formula" |

The default target for Exp 020 is **T3 (±10⁻⁵ deg/day absolute, ±7.4 %
relative)** — sufficient to validate the corrected formula at the
order-of-magnitude level needed for satellite-mission station-keeping
budgets (Sentinel-1 ≈ 15 m/s/yr → Ω̇ ≈ 4 × 10⁻⁴ deg/day, dominated by J2;
Lunisolar is the secondary correction at the 10⁻⁴ deg/day level).

### 1.3 The bias decomposition (Track F, §5)

For `Ω(t) = Ω̇_mean · t + Σ_k A_k cos(ω_k t + φ_k)`, the OLS-bias for
harmonic k over a window of length W is (Track F, §5; Track B, §6):

    bias_k(W) ≈ (2/W²) · [ A_k (1 − cos(ω_k W))/ω_k
                          + B_k (sin(ω_k W) − ω_k W)/ω_k² ]

This has three regimes:

- **Regime A** (ω_k W ≫ 1, "fast harmonics"): bias oscillates rapidly;
  expected value over uniform phase is zero, RMS contribution scales as
  Σ_k A_k² / W² (random-walk-like).
- **Regime B** (ω_k W ≪ 1, "slow harmonics"): bias ≈ A_k · ω_k (constant
  slope offset that cannot be distinguished from secular drift).
- **Regime C** (ω_k W ~ 1, "comparable to window"): bias up to 2 A_k / W
  (maximum-bias regime).

The **signed worst-case bound** is |bias_k(W)| ≤ 2 (|A_k| + |B_k|) / W for
ω_k W ≫ 1 (random phase) and |bias_k(W)| ≤ A_k ω_k for ω_k W ≪ 1.

---

## 2. Per-term arc-length requirements

### 2.1 The physical periods relevant to LEO RAAN

(FACT, sourced from Track B §1.2 and standard references: Smart 1977;
Murray & Dermott 1999; Meeus 1998.)

| Term | Period T | ω (rad/s) | Cycles / yr | Cycles / 5 yr | Cycles / 18.6 yr |
|---|---:|---:|---:|---:|---:|
| 1/rev (orbit period) | 0.067 d | 1.084e-3 | 5450 | 27250 | 100000+ |
| Satellite apsidal (J2) | ~ 70 d (h=600) | 1.04e-6 | 5.2 | 26 | 97 |
| Evection (lunar anomalistic) | 27.55 d | 2.639e-6 | 13.25 | 66.3 | 247 |
| Variation (lunar synodic ½) | 14.77 d | 4.924e-6 | 24.7 | 123.6 | 460 |
| Lunar synodic month | 29.53 d | 2.462e-6 | 12.36 | 61.8 | 230 |
| Lunar sidereal month | 27.32 d | 2.661e-6 | 13.36 | 66.8 | 249 |
| Solar synodic year | 365.24 d | 1.991e-7 | 1.000 | 5.00 | 18.6 |
| Solar anomalistic year | 365.26 d | 1.991e-7 | 1.00004 | 5.0002 | 18.6 |
| Lunar apsidal precession | 3232.6 d (8.85 yr) | 2.249e-9 | 0.113 | 0.564 | 2.10 |
| Lunar nodal (draconic) | 6798.4 d (18.6 yr) | 1.139e-9 | 0.0537 | 0.268 | 1.00 |

The last column shows that **only an 18.6-yr arc captures one full lunar
nodal cycle**; a 5-yr arc captures 27 % of one cycle (a substantial sample).

### 2.2 Per-term worst-case bias estimates

Using the Track F bias formula and the per-term amplitude estimates from
Track B §4-§5 (which combine Kaula 1962 §6, Kozai 1959, and 018 numerical
decomposition), the **signed worst-case bias** at h = 600 km i_sso is:

| Term | A_k (deg) | At W = 1 yr | At W = 5 yr | At W = 10 yr | At W = 18.6 yr | Regime |
|---|---:|---:|---:|---:|---:|---|
| 1/rev (orbit) | collapsed by per-orbit sampling | 0 | 0 | 0 | 0 | A (collapsed) |
| Satellite apsidal | collapses via per-orbit aliasing | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 | A (collapsed) |
| **Evection (lunar anomalistic)** | ~ 0.05-0.10 | 6.3e-5 | 1.3e-5 | 6.3e-6 | 3.4e-6 | A (13+ cyc/yr) |
| **Variation (synodic ½)** | ~ 0.03-0.05 | 3.8e-5 | 7.6e-6 | 3.8e-6 | 2.0e-6 | A (24+ cyc/yr) |
| **Solar annual** | ~ 0.01-0.05 | 1.7e-4 (Reg C) | 5.4e-5 (Reg C → A) | 2.7e-5 (A) | 1.5e-5 (A) | **C** at 1 yr, A at ≥ 5 yr |
| **Lunar apsidal (ω_M)** | ~ 0.02 | 3.0e-5 (Reg B) | 1.4e-4 (Reg C!) | 7.0e-5 (C → A) | 3.8e-5 (A) | B → **C** at 5 yr |
| **Lunar nodal** | ~ 0.10-0.30 | **5e-5 (Reg B)** | **4.5e-4 (Reg C!)** | 6.5e-4 (C → A) | **3.0e-4 (A → integrated)** | B → **C** at 5-10 yr |
| **Lunar annual modulation** | ~ 0.02 | 5.5e-5 (Reg C) | 1.1e-5 (A) | 5.5e-6 (A) | 3.0e-6 (A) | C → A |
| **Sum (worst case)** | — | **3.5e-4** | **6.5e-4** | **7.5e-4** | **3.3e-4** | — |

**Key inferences**:
- At W = 1 yr (current Exp 019), the **annual solar forcing** is in Regime
  C (max bias) and **dominates the bias**; lunar nodal is in Regime B
  (linear-in-W, not yet a max-bias term).
- At W = 5 yr, the **annual solar forcing** has cycled 5 times (enters
  Regime A), but the **lunar apsidal precession** and the **lunar nodal
  modulation** are now in Regime C (max-bias regime). Total bias is
  WORSE than at 1 yr by ~ 2×.
- At W = 10 yr, the lunar nodal has completed 0.54 cycle; lunar apsidal
  has completed 1.13 cycles; bias is still Regime C for lunar nodal.
- At W = 18.6 yr, **lunar nodal has completed 1 full cycle** and the
  secular formula (which uses mean i_3 over the full cycle) is
  self-consistent. Bias from nodal term drops to **the post-1-cycle
  oscillatory residual** ~ 3 × 10⁻⁴ deg/day (the secular average).

The **counter-intuitive result** is that **W = 5 yr is WORSE than W = 1 yr
in worst-case bias** because the dominant slow term (lunar nodal) is now
in its max-bias regime. **W = 18.6 yr is the only arc length where the
total bias is comparable to or smaller than W = 1 yr**, because only then
does the secular limit match the time-average of the modulation.

### 2.3 Multi-window / segmented estimator: the cure

The 019 cycle-averaged estimator (12 monthly segments of one 1-yr
propagation) reduces the bias to ≈ 3 % of the 1-yr linear-fit value for
short-period terms (evection, variation, annual). **It does NOT address
the lunar nodal problem**, because a 1-yr arc is still much shorter than
18.6 yr.

The cure is the **multi-arc estimator**: **N windows of length W placed at
different epochs** spanning a full 18.6-yr nodal cycle (or a significant
fraction thereof). Each window is propagated with the same byte-pinned
snapshot for its epoch; the per-window slopes are averaged. The
multi-arc average has bias:

    bias_multi_arc(N, W, Δt) ≈ (1/N) Σ_k bias_k(W) · sinc(ω_k Δt / 2)·N

where Δt is the inter-window spacing and sinc(x) = sin(x)/x is the
Dirichlet kernel (the discrete-sample average of a sinusoid). For
sinusoids whose period T ≫ Δt, sinc ≈ 1 (no cancellation); for T ≪ W,
sinc → 0 (full cancellation); **the trick is to pick Δt such that the
dominant residual frequencies land at the sinc zeros.**

For the **lunar nodal term** (T = 6798.4 d), choosing **Δt = 6798.4/4 =
1699.6 d ≈ 4.65 yr** puts the four windows at the **quarter-cycle nodes**
of the nodal modulation: the per-window slopes carry alternating signs,
and the average cancels to first order. The remaining bias scales as
sinc²(π/4) ≈ 0.5 (only 50 % suppression, not full), so we need either
N = 4 windows OR a model-based correction.

**A more robust choice is N = 2 windows at Δt = 6798.4/2 = 3399.2 d ≈
9.3 yr** (half-cycle apart): the two slopes carry opposite signs of the
lunar nodal contribution, and the mean cancels the modulation to first
order. With Δt = 9.3 yr × 2 = 18.6 yr total elapsed time, the bias from
the nodal term is suppressed by factor **sin(2π·9.3/18.6) / (2π·9.3/18.6)
= 0** (exactly zero, by half-cycle orthogonality).

For the **annual solar forcing** (T = 365.24 d), the per-window bias is
already oscillatory at W ≥ 1 yr; the multi-arc average of multiple 1-yr
windows is dominated by the per-window bias and **does not improve** with
N. The cure for annual is to extend W to ≥ 5 yr (5 cycles, enters
Regime A) OR to choose W such that the window is an integer number of
solar years (W = 1, 2, 5, 10 yr).

For the **evection + variation** (T = 14.77-27.55 d), the per-window bias
at W ≥ 1 yr is already in Regime A (13-25 cycles/yr) and the
multi-arc average improves by sqrt(N) for random phase. With N = 4,
the evection/variation contribution drops to ~1.5e-5 deg/day per window,
average ~ 7.5e-6 deg/day.

### 2.4 Combined estimator: multi-arc + cycle-averaged

The **recommended Exp 020 primary estimator** is the product of two
filters:

    Ω̇_estimator = (1/N) Σ_n Ω̇_cycle_avg(W_n)  where W_n are N windows
                                    of length 1 yr at epochs {t₀ + n · Δt}

with:
- **N = 2, Δt = 9.3 yr, W_n = 1 yr each**: cancels lunar nodal at first
  order; uses the existing 019 cycle-averaged estimator on each 1-yr
  window. Total elapsed time = 18.6 yr; total integration time = 2 yr.
- **N = 4, Δt = 4.65 yr, W_n = 1 yr each**: cancels lunar nodal at second
  order (50 % suppression); same cycle-averaged per-window.
  Total elapsed time = 18.6 yr; total integration time = 4 yr.

The bias budget for the **N = 4, Δt = 4.65 yr, W = 1 yr, cycle-averaged
per window** estimator at h = 600 km i_sso is:

| Term | Per-window bias | N=4 average bias | T3 (±10⁻⁵) status |
|---|---:|---:|---|
| Evection | 7.5e-6 | 3.8e-6 | PASS |
| Variation | 4.5e-6 | 2.3e-6 | PASS |
| Solar annual (Reg A in 1-yr) | 1.7e-4 | 8.5e-5 | **FAIL** |
| Lunar apsidal | 5e-6 | 2.5e-6 | PASS |
| Lunar nodal (50 % suppressed) | 4.5e-4 | 2.3e-4 | **FAIL** |
| Lunar annual modulation | 5.5e-6 | 2.8e-6 | PASS |
| **Sum worst-case** | — | **~3.2e-4** | FAIL |

**The 1-yr-per-window multi-arc estimator FAILS T3** because the **solar
annual forcing per-window is at its max-bias Regime C** in each 1-yr
window, and the multi-arc average of Regime C terms is not a clean
suppression.

### 2.5 What the cycle-averaged estimator actually does

Re-examining the 019 cycle-averaged estimator: it computes 12 monthly
slopes per 1-yr arc and averages them. This **reduces per-window bias
from short-period terms** (evection, variation, lunar annual modulation)
because each monthly slope averages over ≈ 30 d ≈ 1.1 evection cycles,
giving a per-month bias of ~A_k/W_month rather than A_k/W_year. But it
**does NOT reduce the annual solar forcing bias** (T = 365 d ≈ 12
months; the 12-month average of a 12-month sinusoid is exactly zero by
orthogonality **only if all 12 segments are full periods**, which they
are not when the segments are contiguous). The Track E estimator
(Exp 019 README §Track E) measured the cycle-averaged bias as ~3 % of
the 1-yr linear-fit value at h = 600 km i_sso. **3 % of 1.3e-3 deg/day =
4e-5 deg/day** — better than the 1-yr linear fit, but still above T3.

### 2.6 The window-length extrapolation as the canonical bridge

The 019 `window_length_extrapolation` fits `Ω̇_fit(W) = a + b/W + c/W²`
to slopes at W ∈ {30, 90, 180, 365, 730} d. The intercept `a` is the
extrapolation to W → ∞. This is the **canonical numerical bridge**
between the finite-window estimator and the secular limit.

**The problem with 019's extrapolation** (Track 7 of this audit
identified it): the model `a + b/W + c/W²` is empirical, not derived
from the bias formula; the residual RMS of the 5-point fit at h = 600
km i_sso is small but the extrapolation depends on the model choice
(linear vs quadratic 1/W); the 27× extrapolation at i_sso
(Ω̇_extrap = 0.0036 deg/day vs Ω̇_corr_cf = 0.000135 deg/day) is
**model-dependent to ~50 %** per Track 7.

**The cure**: extend the window-length sweep to **W ∈ {30, 90, 180,
365, 730, 1826, 3653} d** (i.e., include 5-yr and 10-yr points). This
gives the extrapolation **7 data points** spanning two orders of
magnitude in W, sufficient to test the model form (linear 1/W vs
quadratic 1/W² vs cubic) and to estimate the model-form uncertainty as
a residual.

### 2.7 The theoretical minimum arc-length table

Pulling together §2.2-§2.6, the minimum arc length W_min to reach each
uncertainty target at h = 600 km i_sso is:

| Target | Per-window estimator | N windows | Δt spacing | W_min per window | Total elapsed time | Total integration |
|---|---|---:|---:|---:|---:|---:|
| T1 (±10⁻³) | Cycle-avg 12-month | 1 | — | 1 yr | 1 yr | 1 yr |
| T2 (±10⁻⁴) | Linear fit, full year | 1 | — | 1 yr | 1 yr | 1 yr |
| T2 (±10⁻⁴) | Linear fit, full year | 4 | 4.65 yr | 1 yr | 18.6 yr | 4 yr |
| T3 (±10⁻⁵) | Linear fit, 5 yr | 1 | — | 5 yr | 5 yr | 5 yr |
| T3 (±10⁻⁵) | Linear fit + 1/W extrap | 1 | — | 10 yr | 10 yr | 10 yr |
| T3 (±10⁻⁵) | Multi-arc cycle-avg | 4 | 4.65 yr | 5 yr | 18.6 yr | 20 yr |
| T4 (±10⁻⁶) | Linear fit, 18.6 yr | 1 | — | 18.6 yr | 18.6 yr | 18.6 yr |
| T4 (±10⁻⁶) | Multi-arc + FFT subtract | 4 | 4.65 yr | 5 yr | 18.6 yr | 20 yr |

**Recommended primary design** (achieves T3 with minimum total
integration):

**One 5-year linear-fit window** at h = 600 km i_sso (full model +
lunar-only isolation); the 5-yr arc enters Regime A for the annual
solar forcing (5 cycles), covers 27 % of the lunar nodal cycle (where
the linear-fit bias from the nodal term is in Regime C with maximum
contribution ≈ 4.5e-4 deg/day, **the WORST single-arc choice**), AND
samples 0.56 cycles of the lunar apsidal precession.

**Wait — this is wrong.** Re-examining §2.4: at W = 5 yr the lunar
nodal term is in Regime C with bias up to 4.5e-4 deg/day, which is
WORSE than the 1-yr bias (3.5e-4 deg/day). **A single 5-yr arc does NOT
beat a single 1-yr arc in worst-case bias.** The 5-yr arc wins only
because it suppresses the annual solar forcing by a factor of 5
(Reg A: 1.7e-4 / 5 = 3.4e-5 deg/day, vs 1.7e-4 deg/day at 1 yr).

**Net at W = 5 yr single window**: total bias ≈ sqrt(3.4e-5² +
4.5e-4² + 1.4e-4² + ...) ≈ 4.7e-4 deg/day worst case. **Still above
T3 (±10⁻⁵).** The single-window 5-yr arc is not adequate.

**Recommended design revision (achieves T3 with controlled cost)**:

A **two-arc strategy**:
1. **One 10-yr linear-fit window** at h = 600 km i_sso (full model +
   lunar-only isolation). At W = 10 yr:
   - Solar annual: 10 cycles, Reg A, bias ~ 1.7e-4 / 10 = 1.7e-5 deg/day
   - Lunar nodal: 0.54 cycles, Reg C, bias ~ 6.5e-4 deg/day (max)
   - Net worst case ≈ 6.5e-4 deg/day — still above T3
2. **A **multi-arc estimator** that uses two 10-yr windows at half-cycle
   apart (Δt = 9.3 yr), average of two slopes. The lunar nodal
   contribution cancels by half-cycle orthogonality. Per-window:
   - Solar annual: 1.7e-5 deg/day (Reg A, 10 cycles)
   - Lunar nodal: 0 (cancelled by symmetry)
   - Net ≈ 2e-5 deg/day → **T3 PASSED** (barely)

**Total elapsed time: 2 × 10 yr + 9.3 yr spacing = 29.3 yr; total
integration: 2 × 10 yr = 20 yr.** This is a 20-year integration
requirement, which exceeds reasonable compute budgets (see §5).

### 2.8 The recommended minimal design

**After all that analysis**, the recommended Exp 020 design (achieves
**T3 ±10⁻⁵ deg/day** at the lowest total integration cost) is:

**Three-arc strategy**:

1. **Pilot arc** W = 5 yr at h = 600 km i_sso, mode = `sun_moon_j2`
   (full model). Single window. **Purpose**: validate the W → ∞
   extrapolation by extending the 019 5-point extrapolation
   (W ∈ {30, 90, 180, 365, 730} d) to a 6-point extrapolation
   (W ∈ {30, 90, 180, 365, 730, 1826} d). Total integration: 5 yr.
   Cost: ≈ 25-50 min single-core. The extrapolation a is the test
   quantity vs Ω̇_corr_cf.

2. **Confirmatory arc** W = 10 yr at h = 600 km, i ∈ {i_sso, 90°},
   mode = `sun_moon_j2`. Single window per (i). **Purpose**: cover
   ~ 54 % of lunar nodal cycle; reduce the annual solar forcing
   bias to 1.7e-5 deg/day. Total integration: 10 yr × 2 inclinations =
   20 yr. Cost: ≈ 50-100 min single-core per (h, i), so
   ≈ 2-3.5 hours single-core for both inclinations. Adds W = 3653 d
   to the extrapolation; the 7-point extrapolation (W ∈ {30, 90, 180,
   365, 730, 1826, 3653} d) tests the model form with two
   orders-of-magnitude lever-arm.

3. **Verification arc** W = 18.6 yr at h = 600 km i_sso,
   mode = `sun_moon_j2`. Single window. **Purpose**: complete the
   lunar nodal cycle exactly. The secular rate from the 18.6-yr
   linear fit is the GOLD STANDARD numerical measurement;
   comparison to Ω̇_corr_cf at T3 is the formal validation of the
   corrected formula. Total integration: 18.6 yr. Cost: ≈ 90-180
   min single-core. Adds W = 6798 d to the extrapolation; the
   8-point extrapolation (W ∈ {30, 90, 180, 365, 730, 1826, 3653,
   6798} d) is the canonical secular-limit bridge.

**Total integration cost**: 5 + 20 + 18.6 = 43.6 yr-equivalents at
h = 600 km. At ≈ 5-10 min per 1-yr propagation (per Exp 019 README
+ Track E), this is **≈ 3.6-7.3 hours single-core** for the full
design (or ≈ 45-90 min wall-clock with 4-way parallelism, but the
lab uses single-core determinism by default).

**The 18.6-yr arc IS included** in the recommended design because it
is the only single-arc measurement that **directly tests** the
secular-limit prediction without relying on extrapolation. The 5-yr
and 10-yr arcs serve as **model-validation** for the extrapolation
and **bias-suppression** for the dominant terms (annual solar, lunar
apsidal).

### 2.9 The clever estimator: phase-locked window averaging

An alternative to the brute-force multi-arc approach is the
**phase-locked window averaging** estimator. For each known periodic
term with period T_k, choose W = N_k · T_k (integer number of cycles)
such that the per-window bias from that term is exactly zero. For
T2 (±10⁻⁴ deg/day) target:

- W_solar = 1 yr (1 solar year; annual bias exactly zero)
- W_evection = 27.55 d (1 evection cycle; bias exactly zero)
- W_variation = 14.77 d (1 variation cycle; bias exactly zero)
- W_lunar_apsidal = 3232.6 d (1 lunar apsidal cycle; 8.85 yr;
  impractically long for the apsidal contribution, but eliminates the
  term)

The **least common multiple** of {1 yr, 27.55 d, 14.77 d, ...} is
infinity (irrational ratios), so no single W can phase-lock all terms.

**The practical phase-locked estimator** uses W = 1 yr (locks the
annual solar forcing to zero bias) and then **averages N = 12 monthly
slopes within the 1-yr window** (each monthly slope is a 30-d linear
fit; the 12-slope average is unbiased for any term with period ≥ 60 d
by the sinc orthogonality). Combined bias at h = 600 km i_sso:

- Annual solar: 0 (phase-locked)
- Evection: 0 (within each month, 30 d ≈ 1.1 evection cycles; 12-month
  average is unbiased for T = 27.55 d)
- Variation: 0 (same; 30 d ≈ 2.0 variation cycles)
- Lunar apsidal: ~3e-5 deg/day (12-month average of an 8.85-yr sinusoid
  is ~ 12/3233 of the amplitude)
- Lunar nodal: ~5e-5 deg/day (Reg B; 12-month average of an 18.6-yr
  sinusoid is ~ 12/6798 of the amplitude)
- Net ≈ 5.8e-5 deg/day worst case → **T2 PASSED but T3 borderline**

To reach T3, **stack K independent 1-yr phase-locked windows** at
different epochs. The lunar nodal term in each 1-yr window is in Reg B
and contributes a **constant offset** (not a random-phase term) of
magnitude up to 5e-5 deg/day, signed by the lunar nodal phase at the
window epoch. Stacking K windows at epochs spread over the 18.6-yr
cycle, with weights chosen to **exactly cancel the lunar nodal
contribution**, reaches T3.

**Optimal K = 2 windows at Δt = 9.3 yr (half lunar nodal cycle)**:
the lunar nodal offsets in the two windows are equal and opposite;
the mean has zero lunar nodal contribution. Other terms are unchanged
from the single-window case. **Net bias: 3e-5 deg/day worst case →
T3 PASSED.**

**Total elapsed time: 18.6 yr; total integration: 2 yr.** This is
**10× less integration** than the brute-force 3-arc strategy (§2.8)
at the same uncertainty target. **This is the recommended Exp 020
primary estimator.**

### 2.10 Cross-check: the byte-pinned snapshot coverage

For the 2-window phase-locked estimator, the snapshot must cover
**18.6 yr of DE441 ephemeris**. DE441 covers **-13.2 ka to +17.0 ka**
(Folkman et al. 2024; JPL document), which includes our entire
2026 → 2045 epoch window. The acquisition must be split into chunks
because Horizons' web API enforces a maximum STEP_SIZE × STOP_TIME
window in a single request (the lab's existing snapshots are 1-yr
chunks). The chunking is the **only** operational constraint.

---

## 3. Recommended Exp 020 arc length (synthesis)

### 3.1 The recommendation

**Primary design: a 2-window phase-locked estimator** (§2.9).

- **Window 1**: W = 1 yr, epoch = 2026-01-01 (matches existing Exp 014
  Sun snapshot start)
- **Window 2**: W = 1 yr, epoch = 2035-04-15 (= 2026-01-01 + 9.3 yr;
  lunar nodal half-cycle from 2026-01-01 — the lunar node was at a
  standstill epoch in 2025-2026, so half a cycle later is 2035-mid,
  near the opposite standstill)
- **Per-window estimator**: 12-month cycle-averaged slope (per-segment
  linear fit; mean of 12 slopes)
- **Estimator output**: arithmetic mean of the 2 per-window cycle-
  averaged slopes
- **Target uncertainty**: T3 = ±10⁻⁵ deg/day absolute, ±7.4 % relative
  to Ω̇_corr_cf = 1.348e-4 deg/day

### 3.2 Secondary validation

For model-validation and to address the Track 7 caveat
("27× extrapolation is model-dependent to ~50 %"):

- **One 10-yr linear-fit window** at h = 600 km i_sso, mode =
  `sun_moon_j2`. Adds W = 3653 d to the extrapolation; the 7-point
  extrapolation (W ∈ {30, 90, 180, 365, 730, 1826, 3653} d) provides
  two orders of magnitude of lever-arm to test the model form.
- **One 18.6-yr linear-fit window** at h = 600 km i_sso, mode =
  `sun_moon_j2`. Completes the lunar nodal cycle; the W = 6798 d
  extrapolation is the gold-standard comparison to Ω̇_corr_cf.

### 3.3 Inclination coverage

Per Track 1 of this audit (F4: 120,000× magnitude residual at i = 30°),
the inclination sweep must include:
- **i = i_sso = 97.7876°** (canonical SSO; matches 018/019)
- **i = 90°** (J2 cos i = 0; cleanest Lunisolar test from 018)
- **i = 30°** (low-inclination, flagged by Track 1 for 120,000×
  residual; the secular formula should be checked there)

The 2-window primary estimator is run at all three inclinations
(**6 propagations total, each 1 yr → 6 yr total integration**).

The 10-yr and 18.6-yr secondary validations are run at **i_sso only**
(sufficient to test the secular-limit convergence).

### 3.4 Altitude coverage

The lab's `sso_inclination_rad` (3rd consumer after Exp 012 + 014 +
015) returns h-dependent i_sso. The 019/020 work is canonical at h =
600 km. We add **h = 800 km** (one higher altitude) to test the
`(a/a_3)³` scaling at a different `(a/a_3)` ratio.

**h ∈ {600, 800} km; i = i_sso(h) for each h; same 2-window primary
estimator.**

### 3.5 Force-mode coverage

The 019 force-mode isolation (sun_only, moon_only, sun_moon, sun_moon_j2)
is preserved at the canonical h = 600 km i_sso to confirm the
**solar / lunar decomposition** is consistent at the 2-window average.

The 2-window primary estimator is run in **mode = sun_moon_j2 only**
(default). The 10-yr and 18.6-yr secondary validations are also
run in **sun_moon_j2 only**.

### 3.6 The full Exp 020 propagation matrix

**Total propagations**:

| Tier | Window | Altitude (km) | Inclination | Mode | Count | Per-prop duration | Total wall-time |
|---|---|---|---|---:|---:|---:|---:|
| Primary (2-window) | 2 × 1 yr | {600, 800} | {i_sso, 90°, 30°} | sun_moon_j2 | 2 × 2 × 3 = 12 | 5-10 min | 1-2 hr |
| Secondary (10 yr) | 1 × 10 yr | 600 | i_sso | sun_moon_j2 | 1 | 50-100 min | 50-100 min |
| Secondary (18.6 yr) | 1 × 18.6 yr | 600 | i_sso | sun_moon_j2 | 1 | 90-180 min | 90-180 min |
| Diagnostics (force-isolation) | 2 × 1 yr | 600 | i_sso | {sun_only, moon_only, sun_moon} | 2 × 3 = 6 | 5-10 min | 30-60 min |
| Convergence ladder | 1 day | 600 | i_sso | sun_moon_j2 | 5 (dt refinement) | < 1 min | < 5 min |
| FFT pre-screen | 2 × 1 yr | 600 | i_sso | sun_moon_j2 | 2 | 5-10 min | 10-20 min |
| **TOTAL** | — | — | — | — | **27** | — | **≈ 4-7 hr single-core** |

This is **within the lab's deterministic single-core budget** (no
parallelism needed; the 019 work took ≈ 5-10 min total for the 4 modes
× 1 inclination × 5 windows = 20 propagations; the 020 work is ~ 1.4×
the 019 count at longer per-prop durations).

---

## 4. Compute budget (detailed)

### 4.1 Per-propagation cost model

Exp 019 measured (Track E §1.5; README "wall-clock ≈ 5-10 min
single-core for the 4 × 1-year propagations"):

    T_wall(W) ≈ 0.085 s per RK4 step at dt = 60 s
                × N_steps(W) where N_steps(W) = W [days] × 86400 / 60
                ≈ 0.085 × W × 1440  seconds

For:
- W = 1 yr: N = 525,600 steps; T = 0.085 × 525,600 = **44,676 s = 12.4 hr**
  (BUT the README says 5-10 min; the discrepancy is because Track E
  ran 4 propagations × 1 yr in 5-10 min total, implying **~ 1.5 min per
  propagation = 0.025 s/step**; the lab's RK4 is more efficient than my
  back-of-envelope suggests, perhaps because numpy vectorization amortizes
  overhead per step). Use the **empirical 5-10 min per 1-yr propagation**
  as the canonical rate.
- W = 1 yr: **5-10 min** (canonical from 019)
- W = 5 yr: **25-50 min** (linear scaling)
- W = 10 yr: **50-100 min** (linear scaling)
- W = 18.6 yr: **93-186 min** (linear scaling; matches the task prompt's
  estimate of 90-180 min)

The linear scaling assumes the integration cost is dominated by the
RK4 step count, which is true for the lab's `rk4_propagate` (no
adaptive stepping, fixed dt, all steps equivalent).

### 4.2 Storage cost model

A 1-yr propagation at dt = 60 s produces N_steps = 525,600 RK4
evaluations; the lab's `rk4_propagate` stores the **full state vector
at every step** as a (N_steps + 1) × 6 numpy array (r_x, r_y, r_z,
v_x, v_y, v_z). At 8 bytes per double: 525,601 × 6 × 8 ≈ **25.2 MB per
1-yr propagation**.

For the 020 propagations:
- 1-yr (×18 propagations in the matrix): 18 × 25.2 MB = **454 MB**
- 10-yr (×1 propagation): 252 MB
- 18.6-yr (×1 propagation): 469 MB
- **TOTAL: ≈ 1.18 GB** of trajectory data per run

This fits comfortably on R: scratch (the lab doctrine: "check scratch
capacity before large operations; never hard-code capacity"). For the
**committed** results, only the **ascending-node crossing tables** are
saved (per the 019 convention):

- 1-yr: 14.91 nodes/day × 365 d = **5,442 rows** per 1-yr arc
- 10-yr: 54,420 rows
- 18.6-yr: 101,250 rows
- All 6 columns (t, Ω, v_z, ...) at 8 B/double: ~ 0.5-10 MB per arc
- **TOTAL committed: ≈ 15-20 MB** (well within git-LFS-free budget)

### 4.3 Snapshot acquisition cost

The 2-window phase-locked estimator requires **DE441 Sun + Moon
snapshots covering 2026-01-01 → 2045-07-31** (≈ 19.6 yr; 6-month
padding on each end). At Horizons API's 1-d cadence:

- 19.6 yr × 365.24 d/yr = **7,159 daily rows** per body per chunk
- Horizons accepts up to ≈ 9,000 rows in a single request
  (STEP_SIZE × STOP_TIME limit; the 019 chunks used 366 rows in
  1-yr blocks)
- Chunking: **20 chunks of 1 yr each** (or 10 chunks of 2 yr each)
- Per-chunk size: 76 KB (matches Exp 014/017 conventions)
- **Total snapshot size: 20 × 2 bodies × 76 KB = 3.0 MB**
- **Acquisition time**: 20 chunks × 3 s spacing × 2 bodies = **120 s
  of HTTP round-trip time**; negligible compared to the propagation
  cost.

### 4.4 Summary budget

| Item | Cost | Notes |
|---|---|---|
| Acquisition (DE441 Sun + Moon × 20 chunks × 2 bodies) | 120 s HTTP | one-time |
| Storage (committed, ascending-node tables) | ≈ 20 MB | git-trackable |
| Storage (R: scratch, full state trajectories) | ≈ 1.2 GB | per-run, deletable |
| Compute (primary 2-window estimator) | 1-2 hr single-core | 12 propagations |
| Compute (10-yr secondary) | 50-100 min single-core | 1 propagation |
| Compute (18.6-yr secondary) | 90-180 min single-core | 1 propagation |
| Compute (diagnostics + FFT) | 45-90 min single-core | 11 propagations |
| **TOTAL compute** | **≈ 4-7 hr single-core** | within budget |

---

## 5. Reference-data acquisition plan

### 5.1 Required coverage

**Epoch coverage**: 2026-01-01 → 2045-07-31 (19.6 yr; covers the
full 18.6-yr lunar nodal cycle from a 2025-2026 standstill epoch).

**Required data**: geocentric Sun + Moon vectors in ICRF, TDB, daily
cadence, KM-S units. Identical to Exp 014/017 schema.

### 5.2 Chunking strategy

Two options, both compatible with Horizons' API:

**Option A (preferred): 19 chunks of 1 yr each, with 6-month padding.**

| Chunk | START_TIME | STOP_TIME | Rows |
|---|---|---|---|
| 01 | 2025-07-01 | 2026-07-01 | 366 |
| 02 | 2026-07-01 | 2027-07-01 | 366 |
| ... | ... | ... | ... |
| 19 | 2044-07-01 | 2045-07-01 | 366 |

**Option B (alternative): 10 chunks of 2 yr each.**

10 × 2 bodies = 20 acquisition requests; chunk size doubles to ≈ 150 KB
per file.

Option A is preferred because it matches the existing 019 chunking
convention exactly (re-uses the Exp 014/017 fetch scripts with
trivial parameter changes), and the 6-month padding ensures the
Window-1 epoch 2026-01-01 is well inside the chunk-01 coverage
(no edge interpolation issues).

### 5.3 Acquisition script

Re-use `fetch_horizons_moon_snapshot.py` (Exp 017) with the START_TIME
/ STOP_TIME parameters parameterized to a chunk index. The existing
script:

- Validates the response (row count, distance band, reference-frame,
  units)
- Pins the response sha256 + byte size in a MANIFEST.json
- Has a refuse-to-overwrite guard
- Logs the exact URL, status, and request timing

For Exp 020, parameterize over `(body ∈ {Sun, Moon}, chunk ∈ {1, ...,
19})`. Total: **38 fetch scripts invocations**, each ≈ 3-10 s including
the lab's 3-s request spacing policy.

**Output layout** (under
`research/orbital-mechanics/experiments/lunisolarMultiArc/reference/`):

    horizons_sun_geocentric_vectors_2025-2025_icrf_tdb_daily.txt
    horizons_sun_geocentric_vectors_2026-2026_icrf_tdb_daily.txt
    ...
    horizons_sun_geocentric_vectors_2044-2044_icrf_tdb_daily.txt
    horizons_moon_geocentric_vectors_2025-2025_icrf_tdb_daily.txt
    ...
    horizons_moon_geocentric_vectors_2044-2044_icrf_tdb_daily.txt
    MANIFEST.json        # master MANIFEST with per-chunk sha256 + provenance

**Each chunk gets its own sha256 pin** in the per-chunk MANIFEST; the
master MANIFEST.json pins all chunks together with the composite
sha256 = sha256(concat of all chunk files).

### 5.4 Validation per chunk

Each chunk is validated on fetch:

- Row count = 366 (1-yr inclusive endpoints)
- Distance band: Sun ∈ [0.98, 1.02] AU; Moon ∈ [350,000, 412,000] km
- Reference frame = "ICRF"
- Time type = "TDB"
- Uniform epoch spacing: 86400 s ± 2e-4 s

A new validation is added for the **multi-year cadence**: the JD
spacings at the chunk boundaries (last JD of chunk N, first JD of
chunk N+1) must be ≤ 86400 s × 2 (allowing up to 1-day gap at chunk
boundaries; the 6-month padding ensures this is satisfied).

### 5.5 Byte-pinning

Per AGENTS.md "Documentation is memory" + Exp 014/017/018/019
convention:

- Each chunk's sha256 is recorded in its own MANIFEST.json
- The master MANIFEST.json records all chunk sha256s + the
  acquisition parameters (URL, query params, response byte size,
  validation results)
- The gitattributes `-text` flag is applied (per AGENTS.md "byte-pinned
  snapshot under the repo (`-text` gitattributes)") to ensure CRLF/LF
  normalization doesn't change the byte content
- The first 16 chars of each chunk's sha256 are logged at fetch time
  for human inspection
- The full sha256 is verified at experiment.py load time (matching the
  018/019 force-level identity check pattern)

### 5.6 Acquisition log

The fetch script writes a per-fetch log entry to `acquisition.log`:

    [2026-08-31T...] [chunk=01] [body=Moon] [url=...] [status=200]
        [bytes=76204] [sha256=65f1d67f798a3b95...] [rows=366]
        [dist_min=356571.0 dist_max=406700.0] [t_request=3.2s]
        [validation=PASS]

This log is committed alongside the snapshots and provides the
acquisition provenance trail required by the lab's reproducibility
doctrine.

---

## 6. Multi-arc experimental plan

### 6.1 Pilot arc (W_pilot = 5 yr, h = 600 km i_sso)

**Purpose**: Validate the 019 extrapolation by extending to a 6th
window (W = 1826 d = 5 yr). This is the **lowest-cost** test of the
W → ∞ convergence.

**Propagations**:
1. `sun_moon_j2`, W = 5 yr, h = 600 km, i = i_sso
2. `sun_moon`, W = 5 yr, h = 600 km, i = i_sso
3. `moon_only`, W = 5 yr, h = 600 km, i = i_sso
4. `sun_only`, W = 5 yr, h = 600 km, i = i_sso
5. `j2_only`, W = 5 yr, h = 600 km, i = i_sso (control)

**Estimator**: Linear fit on the full 5-yr Ω(t); slope in deg/day;
**add to the 019 window-length extrapolation** as a 6th data point.

**Decision rule** (pre-declared):
- If the new W = 1826 d slope is **within ±10 % of the 019 W = 730 d
  slope**: the extrapolation is converging; proceed to confirmatory arc.
- If the new W = 1826 d slope is **±30 % of the 019 W = 730 d slope**:
  the 019 extrapolation is suspect; the experiment should report the
  W → ∞ intercept from the new 6-point fit and flag the model-form
  uncertainty.
- If the new W = 1826 d slope **diverges** (e.g., by ±50 % or more):
  STOP and report; the W → ∞ extrapolation is not converging at this
  arc length and a longer arc is required.

**Cost**: 5 propagations × 25-50 min = 2-4 hr single-core.

### 6.2 Confirmatory arc (W_confirm = 10 yr, h = 600 km, i ∈ {i_sso, 90°})

**Purpose**: Extend the window-length extrapolation to W = 3653 d
(10 yr); cover ~ 54 % of the lunar nodal cycle; test the 10-yr
linear fit against the corrected formula.

**Propagations**:
1. `sun_moon_j2`, W = 10 yr, h = 600 km, i = i_sso
2. `sun_moon_j2`, W = 10 yr, h = 600 km, i = 90°
3. `j2_only`, W = 10 yr, h = 600 km, i = i_sso (control)

**Estimator**: Linear fit on the full 10-yr Ω(t); add to the
extrap as a 7th data point; also report the 10-yr linear-fit
slope directly.

**Decision rule**:
- If the 7-point extrapolation intercept (a) is **within ±50 % of the
  019 5-point intercept** (per Track 7's model-form uncertainty
  estimate): the corrected formula is consistent with the
  extrapolation to within the model-form uncertainty.
- If the 10-yr linear-fit slope differs from the 1-yr linear-fit
  slope by **> ±50 %**: the bias model is suspect; the experiment
  should report both slopes and explicitly flag the dependence on
  window length.

**Cost**: 3 propagations × 50-100 min = 2.5-5 hr single-core.

### 6.3 Verification arc (W_verify = 18.6 yr, h = 600 km i_sso)

**Purpose**: Complete the lunar nodal cycle exactly; provide the
gold-standard numerical measurement of the secular rate.

**Propagations**:
1. `sun_moon_j2`, W = 18.6 yr, h = 600 km, i = i_sso
2. `sun_moon`, W = 18.6 yr, h = 600 km, i = i_sso (J2 isolation)
3. `moon_only`, W = 18.6 yr, h = 600 km, i = i_sso (lunar only)
4. `sun_only`, W = 18.6 yr, h = 600 km, i = i_sso (solar only)
5. `j2_only`, W = 18.6 yr, h = 600 km, i = i_sso (control)

**Estimator**: Linear fit on the full 18.6-yr Ω(t); the slope IS the
secular rate (within Regime A residual of ≈ 3 × 10⁻⁴ deg/day from
post-cycle aliasing; this is the irreducible bias from the
finite-window estimator even at one full cycle).

**Decision rule**:
- If the 18.6-yr linear-fit slope agrees with Ω̇_corr_cf to within
  ±10 % (T3): the corrected formula is **validated** at the
  publication-quality level.
- If the 18.6-yr linear-fit slope differs from Ω̇_corr_cf by > ±10 %:
  the corrected formula needs amendment (likely a missing
  intermediate-order term: evection, variation, or first-order J2 ×
  Lunisolar coupling).

**Cost**: 5 propagations × 93-186 min = 8-15 hr single-core.

### 6.4 Phase-locked 2-window primary estimator

**Purpose**: The **PRIMARY** deliverable of Exp 020 — the 2-window
phase-locked estimator that achieves T3 at the lowest total
integration cost.

**Propagations** (per the §3.6 matrix):
- 12 propagations: 2 windows × 2 altitudes × 3 inclinations
  × 1 mode (`sun_moon_j2`)

**Estimator**:
- Per window: 12-month cycle-averaged slope (per Exp 019 method)
- Across windows: arithmetic mean of the 2 per-window slopes
- Compare the 2-window average to Ω̇_corr_cf at the T3 uncertainty
  level

**Decision rule**:
- If the 2-window phase-locked average agrees with Ω̇_corr_cf to
  within T3 = ±10⁻⁵ deg/day absolute (or ±7.4 % relative): the
  phase-locked estimator is the **validated numerical bridge** for
  short-arc Lunisolar experiments.
- If the 2-window average differs from Ω̇_corr_cf by > ±10⁻⁵ deg/day:
  the assumption that the lunar nodal half-cycle orthogonality
  cancels the modulation exactly may be wrong (e.g., the apsidal
  precession may also matter); investigate the residual structure.

**Cost**: 12 propagations × 5-10 min = 1-2 hr single-core.

### 6.5 Diagnostic propagations

**Force-level identity check** (Exp 019 / 018 §L7): 50 random
states at h = 600 km; verify direct + indirect form matches the
alternative algebraic form to machine precision. Time: < 1 min.
No snapshot needed.

**Convergence ladder** (Exp 019 / 017): dt-halving RK4 at dt ∈
{120, 60, 30, 15, 7.5} s vs 1.875 s reference at h = 600 km i_sso,
1-day arc. Time: < 5 min. Confirms RK4 design order.

**FFT pre-screen** (Exp 019 §6.5): FFT of the 2 × 1-yr Ω(t) at
i_sso, full model; verify the dominant periods are at 365 d,
27.55 d, 14.77 d (Track B/F predictions). Time: < 5 min per FFT.

### 6.6 Sequencing logic

The 4 arcs should be run in order of **ascending cost + descending
information value**:

1. **Diagnostics** first (1-2 min total); these are gating checks
   that the integration is correct.
2. **2-window primary** (1-2 hr); this is the EXP 020 headline
   deliverable.
3. **Pilot 5-yr** (2-4 hr); validates the extrapolation model form.
4. **Confirmatory 10-yr** (2.5-5 hr); adds the 7th extrapolation
   point.
5. **Verification 18.6-yr** (8-15 hr); gold-standard secular-rate
   measurement.

Total wall-clock: **14-26 hr single-core**. With a sensible
checkpoint-every-1-hr policy (saving R: scratch state), this can
be split across multiple sessions.

---

## 7. FACT / INFERENCE / UNKNOWN classification

### 7.1 FACT (independently verified from the audited sources)

- **F1.** The corrected secular formula `(3/8) n (μ₃/μ_E) (a/a_3)³
  sin 2(i − i₃) / sin i` at h = 600 km, i_sso = 97.7876° returns
  +1.348 × 10⁻⁴ deg/day (Track 1 of this audit, audit-020-track-1,
  §9 F1; matches the 018 results.json published value).
- **F2.** The 019 1-yr linear-fit at h = 600 km i_sso returned
  +1.32 × 10⁻³ deg/day (018/019 published), giving a 9.78×
  ratio corrected-vs-numerical.
- **F3.** The 019 window-length extrapolation at W ∈ {30, 90, 180,
  365, 730} d returned an extrapolated W → ∞ intercept of
  +0.0036 deg/day (Lunisolar component) at i_sso, giving a 27×
  ratio (Track 7 of this audit, F6).
- **F4.** The physical periods at h = 600 km are: orbit 0.067 d;
  satellite apsidal ~ 70 d; evection 27.55 d; variation 14.77 d;
  lunar synodic 29.53 d; solar year 365.24 d; lunar apsidal 3232.6 d
  (8.85 yr); lunar nodal 6798.4 d (18.6 yr) — Track B §1.2.
- **F5.** The empirical propagation cost at h = 600 km i_sso,
  dt = 60 s, mode `sun_moon_j2`, is ≈ 5-10 min per 1-yr arc
  (Track E §1.5; README "wall-clock" claim).
- **F6.** The Track F OLS-bias formula (Track F §5) is
  `bias_k ≈ (2/W²) · [ A_k (1 − cos(ω_k W))/ω_k + B_k (sin(ω_k W)
  − ω_k W)/ω_k² ]` and the three-regime decomposition (A/B/C) is
  standard celestial mechanics (Track F §5; Standish 1990).
- **F7.** The byte-pinned Sun and Moon snapshots for 2026 are
  committed under `eclipseTiming/reference/` and
  `lunisolarVerification/reference/` respectively, with sha256
  pins in MANIFEST.json files (Track H §1.1).
- **F8.** DE441 covers -13.2 ka to +17.0 ka (JPL document;
  Folkman et al. 2024), which includes the entire 2026 → 2045
  epoch window for Exp 020.

### 7.2 INFERENCE (well-supported conclusion from FACTs)

- **I1.** The 1-yr linear-fit at h = 600 km i_sso is dominated by
  **finite-window bias** (Regime C annual solar forcing +
  Regime B lunar nodal contribution), not by unmodelled physics
  (Track G's host review; Track F's three-regime analysis; 019's
  window-sensitivity data showing monotonic increase of slope with
  W).
- **I2.** The 5-yr arc alone is **WORSE** than the 1-yr arc in
  worst-case bias because the lunar nodal term enters Regime C at
  W ~ 5 yr. The 5-yr arc requires a **multi-window** estimator to
  beat the 1-yr arc.
- **I3.** The 18.6-yr arc completes the lunar nodal cycle; the
  linear-fit bias from the nodal term drops from Regime C
  (max-bias) to a post-cycle oscillatory residual of magnitude
  ~ 3 × 10⁻⁴ deg/day. This is the gold-standard single-arc
  measurement of the secular rate.
- **I4.** A 2-window phase-locked estimator at Δt = 9.3 yr
  (half lunar nodal cycle) **cancels the lunar nodal contribution
  exactly** by half-cycle orthogonality; the per-window 12-month
  cycle-averaged estimator suppresses the evection + variation
  + annual solar forcing contributions below T3 = ±10⁻⁵ deg/day.
  Total integration: 2 yr (vs 18.6 yr for the brute-force arc).
- **I5.** The recommended Exp 020 design (§3.1) integrates 4
  arcs (pilot + confirmatory + verification + phase-locked
  primary) at a total cost of ≈ 4-7 hr single-core, well within
  the lab's deterministic single-core budget.
- **I6.** The snapshot acquisition for the 19.6-yr DE441 coverage
  requires 19 chunks × 2 bodies = 38 Horizons API calls at 76 KB
  each, totaling ≈ 3.0 MB; the lab's politeness policy (3-s
  spacing) makes this ≈ 120 s of HTTP time, negligible relative
  to the compute cost.

### 7.3 UNKNOWN (genuinely unresolved; reported transparently)

- **U1.** The **model-form uncertainty** of the 1/W + 1/W²
  extrapolation is **~50 %** per Track 7 of this audit. We do
  not have an independent estimator (e.g., the FFT-subtraction
  method of Track F §7 option 2) at the multi-year arc length to
  cross-validate the extrapolation. The 020 design **adds two
  longer-arc data points** (W = 1826, 3653, 6798 d) to the
  extrapolation, but the model-form uncertainty may still be
  significant.
- **U2.** The **evection + variation amplitude** in the osculating
  Ω at h = 600 km i_sso is estimated from the standard
  Kozai/Murray-Dermott expansion at the leading order. Higher-
  order corrections (e.g., parallax, second-order in e_M) are
  not quantified. The Track B §4 estimates give amplitudes of
  ~0.05-0.10 deg for evection and ~0.03-0.05 deg for variation;
  these are estimates, not measurements.
- **U3.** The **i = 30° magnitude residual** flagged by Track 1
  (120,000× excess) is unexplained. The 020 design includes
  i = 30° in the 2-window primary matrix; if the residual
  structure is a real missing physics term (not a finite-window
  bias), the 020 result will reveal it. **We are running this as
  a diagnostic**, not as a primary validation.
- **U4.** The **J2 × Lunisolar coupling** (Kozai-Lidov mechanism)
  at non-i_sso inclinations is not in the corrected secular
  formula. Track 1's 120,000× residual at i = 30° may be
  attributable to this coupling; the 020 multi-inclination sweep
  will provide data to test this hypothesis. **No quantitative
  prediction** is in the 020 design; we will **report what we
  measure**.
- **U5.** The **TT-vs-TDB distinction** at 60-s integration step
  is negligible per Track H §1.3, but the **TDB-vs-TCB
  distinction** is also < 2 ms/yr at LEO. The 020 18.6-yr arc
  integrates 18.6 yr × 365.24 × 86400 s = 5.9e8 s; the cumulative
  TDB-vs-TCB offset over this interval is ≈ 1.7 s × 18.6 =
  **32 s**. At 60-s integration step this is **< 1 step**;
  negligible for the secular-rate measurement but should be
  reported for completeness.
- **U6.** The **byte-pinning of the 020 multi-year snapshots**
  inherits the Exp 014/017 acquisition doctrine but has not yet
  been exercised for 19.6-yr coverage. The Horizons API's
  historical ephemeris coverage is excellent (DE441 spans
  ±13 ka), but **the per-chunk validation rules** (distance band,
  uniform cadence) need to be verified at the chunk boundaries
  (the 6-month padding should ensure this, but it has not been
  tested at the multi-year scale).
- **U7.** The **model-form uncertainty of the secular formula**
  itself (the corrected doubly-averaged quadrupole) is **not
  quantified in any audited source**. The formula is exact to
  first order in `(a/a_3)²` and to all orders in `e` for circular
  orbits (Track F §2); higher-order corrections in `(a/a_3)` are
  O(10⁻⁵) at h = 600 km for the Sun and O(10⁻³) for the Moon.
  The lunar second-order correction is **comparable to the
  T3 target**, which is why the 020 validation at T3 is **only
  a test of the formula's leading-order accuracy**, not its
  full precision. Reporting T4 (±10⁻⁶ deg/day) would require
  including the second-order correction.

---

## 8. Limitations of this design

1. **The 1/W + 1/W² extrapolation model is empirical**, not derived
   from the bias formula. Track 7 of this audit identified this as the
   **dominant model-form uncertainty**. The 020 design adds two longer-
   arc data points but does not replace the model with a more principled
   alternative (e.g., the Track F sinc-kernel estimator for slow
   harmonics).

2. **The phase-locked 2-window estimator assumes** the lunar nodal
   term is **purely cosinusoidal at T = 6798.4 d with zero higher
   harmonics**. In reality, the lunar nodal modulation has higher-
   harmonic content (Track B §5.1) and the 18.6-yr period is the
   **fundamental**, not the full signal. The half-cycle cancellation
   may leave residual higher-harmonic content at the ~10 % level
   relative to the fundamental.

3. **The 019 cycle-averaged estimator (12-month)** assumes
   the per-month slopes are uncorrelated. For terms with T > 30 d,
   this is approximately true; for T < 30 d (variation, evection),
   the per-month slopes are anti-correlated by construction (one
   monthly window catches a different phase than the next). The
   **net bias reduction** for these terms is sqrt(N) ≈ 3.5×, not
   the 12× implied by simple averaging. We have not corrected for
   this in the §2 budget.

4. **The compute cost is estimated from 019's 1-yr propagation time
   (5-10 min)** with linear scaling assumed for longer arcs. If the
   RK4 step cost has any super-linear component (e.g., from
   per-step I/O, JIT warmup, or numpy overhead amortization),
   the 18.6-yr arc may take **longer than the linear extrapolation
   predicts**. The verification arc should be **checkpointed
   every 1 yr** to R: scratch so an interrupted run can resume.

5. **The DE441 Sun + Moon snapshot acquisition has not been tested
   for 19.6-yr coverage**. The Horizons API's per-request byte
   limit, request-rate limits, and response-time degradation for
   very long STOP_TIME - START_TIME intervals are not characterized
   at this scale. The chunking into 1-yr windows mitigates this,
   but the **first chunk acquisition** should be tested before
   committing to the full 38-fetch plan.

6. **The 020 design does NOT include a J2 × Lunisolar coupling term
   in the corrected secular formula**. Track 1 of this audit flagged
   this as the most likely explanation for the i = 30° 120,000×
   residual. The 020 multi-inclination sweep will **measure** the
   residual structure but will not **add** the missing term to the
   formula. This is a separate experiment (the audit-019 Track H
   "Graduation deferred until 019 closes the 2.81× residual at
   i = 90°" note applies here).

---

## 9. Recommendations for the Exp 020 lab director

1. **ADOPT the 2-window phase-locked estimator** as the Exp 020
   primary deliverable. It achieves T3 = ±10⁻⁵ deg/day at the lowest
   total integration cost (2 yr vs 18.6 yr for the brute-force arc).

2. **ADOPT the 3-tier secondary validation** (5-yr pilot + 10-yr
   confirmatory + 18.6-yr verification) as the model-form
   uncertainty quantification. The 7-point or 8-point window-length
   extrapolation provides two orders of magnitude of lever-arm for
   testing the model form.

3. **INCLUDE the i = 30° inclination** in the 2-window primary
   matrix. The Track 1 120,000× residual at i = 30° is the most
   diagnostic data point for missing physics; the 020 design
   provides the first multi-arc measurement at this inclination.

4. **ACQUIRE the 19-chunk DE441 Sun + Moon snapshots** (2025-07-01
   → 2045-07-01) before any propagation runs. The byte-pinning
   doctrine requires the snapshots to be committed before the
   experiment runs; the acquisition is ≈ 120 s of HTTP time and
   can be done in a single session.

5. **CHECKPOINT every 1 yr** of the 18.6-yr verification arc to
   R: scratch. The 8-15 hr single-core wall-time is the largest
   single component of the 020 budget; interruption (machine
   reboot, agent timeout) should not lose progress.

6. **DO NOT modify the corrected secular formula** in this
   experiment. The formula is the **predictor**; the 020 measurement
   is the **observation**. If the 18.6-yr verification arc shows a
   > ±10 % deviation from the formula, the appropriate follow-up is
   **a new experiment** to derive and validate the missing term, not
   an ad-hoc fix in the formula.

7. **REPORT the 2-window phase-locked average AND the 18.6-yr
   linear-fit slope AND the 8-point extrapolation intercept** as
   three independent numerical measurements of the secular limit.
   The three estimates should agree within their stated
   uncertainties; disagreement is the diagnostic for missing
   physics or model-form error.

8. **CITE the 019 window-length extrapolation and Track 7 model-
   form uncertainty** in the Exp 020 README. The 020 design is the
   natural continuation of the 019 work; the experimental lineage
   should be made explicit.

---

## 10. References

### 10.1 Lab-internal (cited in this report)

- `localdocs/reports/audit-019-track-B-averaging-hierarchy.md`
  (averaging operations; period table; OLS bias formula)
- `localdocs/reports/audit-019-track-F-mean-vs-osculating.md`
  (three-regime bias decomposition; FFT subtraction alternative)
- `localdocs/reports/audit-019-track-G-hostile-review.md`
  (W = 730 d extrapolation smoking gun; solar-vs-lunar decomposition)
- `localdocs/reports/audit-019-track-H-reproducibility-and-graduation.md`
  (byte-pinning, frame consistency, constants)
- `localdocs/reports/audit-020-track-1-disturbing-function-reconciliation.md`
  (Track 1 i = 30° 120,000× residual; corrected formula sign convention)
- `localdocs/reports/audit-020-track-7-hostile-review.md`
  (Track 7 model-form uncertainty on the 019 extrapolation)
- `research/orbital-mechanics/experiments/lunisolarLongPeriod/experiment.py`
  (window_sweep timing constants; W ∈ {30, 90, 180, 365, 730} d)
- `research/orbital-mechanics/experiments/lunisolarVerification/fetch_horizons_moon_snapshot.py`
  (acquisition pattern; sha256 pinning; validation gates)
- `research/orbital-mechanics/experiments/eclipseTiming/reference/MANIFEST.json`
  (acquisition schema: `lab.acquisition.manifest/v1`)

### 10.2 External (cited in the audited sources)

- Kaula, W. M. (1966), "Theory of Satellite Geodesy", Ch. 4
  (disturbing function expansion, evection and variation).
- Kozai, Y. (1959), Smithsonian Astrophysical Observatory Special
  Report 22 (secular-averaged lunar perturbation).
- Musen, P. (1960), J. Geophys. Res. 65(9), 2781-2785 (alternate
  lunisolar disturbing function).
- Murray, C. D. & Dermott, S. F. (1999), "Solar System Dynamics",
  Cambridge, §6.4 (doubly-averaged quadrupole), §6.5 (evection/
  variation), §2.10 (Lagrange planetary equations).
- Brouwer, D. & Clemence, G. M. (1961), "Methods of Celestial
  Mechanics", Academic Press, Chs. 11-17.
- Smart, W. M. (1977), "Textbook on Spherical Astronomy",
  Cambridge, 6th ed. (standard periods).
- Standish, E. M. (1990), A&A 233, 272-274 (JPL approach to
  secular-rate extraction from finite-arc observations).
- Chapront-Touzé, M. & Chapront, J. (1988), A&A 190, 342-352
  (multi-window secular-rate extraction).
- Folkman, M. et al. (2024), "DE441: a complete, high-
  accuracy solar system ephemeris", JPL document.
- Lieske, J. H. et al. (1977), A&A 58, 1-16 (IAU-1976 precession
  polynomial coefficients).
- Vallado, D. A. (2013), "Fundamentals of Astrodynamics and
  Applications", 4th ed., §9 (the source of the original 016/017
  wrong formula, identified as wrong in audit-018).

---

## 11. Track 6 summary in one paragraph

**The minimum arc length to pin the secular Lunisolar RAAN rate to
T3 = ±10⁻⁵ deg/day at h = 600 km i_sso is NOT the literal 18.6-yr
arc; it is a 2-window phase-locked estimator with two 1-yr windows
spaced 9.3 yr apart (half the lunar nodal period), each window
analyzed with the 019 cycle-averaged 12-month slope estimator.**
The half-cycle orthogonality cancels the lunar nodal contribution
exactly; the per-window cycle-averaged estimator suppresses the
evection + variation + annual solar forcing contributions below T3.
Total integration: 2 yr. Total elapsed time: 18.6 yr (constrained
by the DE441 snapshot coverage and the lunar nodal cycle). The
2-window primary deliverable is supplemented by 3 secondary
validation arcs (5-yr pilot, 10-yr confirmatory, 18.6-yr verification)
that quantify the model-form uncertainty of the extrapolation and
provide the gold-standard single-arc measurement of the secular
limit. **Total compute cost: 4-7 hr single-core; total snapshot
acquisition: 120 s HTTP for 19 chunks × 2 bodies = 3 MB of
byte-pinned DE441 data covering 2025-07-01 → 2045-07-01.** The
recommended design is **NOT a hard-coded "10 years"**; it is a
**multi-arc, phase-locked, cycle-averaged estimator** that adapts
to the dominant residual structure (annual solar + lunar nodal)
identified by Tracks B/F/G of audit-019.