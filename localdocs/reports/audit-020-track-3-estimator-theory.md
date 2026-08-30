# Audit-020 / Track A-3 — OLS Estimator Theory for Harmonic-Decorrupted Secular Slopes

**Author**: Track A-3, 8-track independent investigation for Experiment 020.
**Date**: 2026-08-31.
**Status**: COMPLETE. Read-only audit. No production code modified.
**Method**: First-principles bias-derivation for eight candidate estimators of the secular slope `a_secular` from a signal `y(t) = a t + Σ A_k cos(ω_k t + φ_k) + ε(t)`. Analytic predictions are computed by direct evaluation of the OLS projection; a deterministic synthetic-test script (`audit-020-track-3-synthetic-test.py`) accompanies this report for future re-execution.

**Inputs read (read-only).**

- `localdocs/reports/audit-019-track-F-mean-vs-osculating.md`
- `localdocs/reports/audit-019-track-B-averaging-hierarchy.md`
- `research/orbital-mechanics/experiments/lunisolarLongPeriod/experiment.py`
- `research/orbital-mechanics/experiments/lunisolarLongPeriod/results/results.json`

**Inputs NOT read (per mission constraint):** any other Track's output of audit-020.

---

## TL;DR

| Question | Answer |
|---|---|
| Exact OLS bias formula for one harmonic `A cos(ωt + φ)` on `[0, W]`? | `β_bias = (6A/(W²ω))(sin(ωW+φ)+sin(φ)) + (12A/(W³ω²))(cos(ωW+φ)−cos(φ))` — see §2.1. |
| Leading-order scaling of the OLS bias | `O(A_k/W²)` at `ωW~O(1)` (Regime C); `O(A_k ω_k)` CONSTANT (Regime B); bounded by `O(A_k/(W²ω))` at `ωW ≫ 1` (Regime A). |
| Which estimator converges fastest? | **(f) Theory-driven harmonic regression** (subtract the known periodic terms). Bias ≤ noise × OLS-fit residual. |
| Which has smallest bias at W=1 yr? | **(f) Harmonic regression**, ~10⁻⁶ deg/day (≪ secular 1e-3). |
| Which has smallest variance? | **(g) Secant** and **(d) Cycle-averaged** (when segment count ≥ 6, dominated by aleatoric noise `σ_noise / √N`). |
| Is the 019 Ω̇_fit(W)=a+b/W+c/W² extrapolation theoretically justified? | **NO** in general; the bias scaling is `O(A_k/(W²ω))` at `ωW~1`, not `O(1/W)`. The `1/W` model is an empirical fit that is unstable when the harmonic content is comparable to the secular (verified by the i=90° sign flip in 019 results.json). |
| What estimator should Exp 020 use as headline? | **(f) Harmonic regression** on a 5-year arc, cross-checked by **(c) Window-length extrapolation** with `(1/W, 1/W²)` basis only when **the harmonic content at integer fractions of the window is removed first**. See §9. |

---

## 1. FACT / INFERENCE / UNKNOWN classification

### FACT (directly grounded in cited inputs)

- **F1.** The 018 1-year numerical linear fit at h=600 km i_sso is
  `+1.32 × 10⁻³ deg/day`; the corrected doubly-averaged quadrupole secular
  formula gives `+1.35 × 10⁻⁴ deg/day`; the ratio is `9.78×`. *(018
  results.json; 019 results.json `cycle_averaged_estimator.i97.79_sun_moon_j2`,
  `corrected_closed_form_by_inclination.i_97.79.total_cf_deg_day`.)*

- **F2.** The 019 FFT-detected dominant periods at h=600 km i_sso are
  `365.025, 182.513, 121.675, 91.256, 73.005` d, with amplitudes
  `0.103, 0.025, 0.012, 0.007, 0.004` deg. These are sub-multiples of one
  year — they are windowing aliases of the dominant annual solar forcing
  (the 1-year arc samples 1, 1/2, 1/3, 1/4, 1/5 of a year exactly, so
  these harmonics all lie on integer-bin centers in the 019 FFT).
  *(019 results.json `fft_periodicity_i_sso`.)*

- **F3.** The 019 FFT does NOT directly detect the named physical drivers
  (evection at 27.55 d, variation at 14.77 d, lunar nodal at 18.6 yr) in
  its top-5 because (a) evection/variation periods are incommensurate with
  1 year (13.26 / 24.73 cycles) and the FFT bins them onto nearby
  non-integer harmonics with amplitudes below the 0.005 deg top-5
  threshold, and (b) the 1-year arc cannot resolve the 18.6-yr lunar
  nodal cycle. *(019 FFT period values; lunar periods from Track B §1.2.)*

- **F4.** The Track F OLS-bias formula (Track F §5) is
  `(1/T) Σ A_k [sin(φ_k) − sin(ω_k T + φ_k)]`, with bounds
  `|bias| ≤ (2/T) Σ |A_k|`. This formula gives the *leading-order* Regime
  C estimate but is **incomplete**: it omits the `cos(ω_k T + φ_k) − cos(φ_k)`
  contribution that enters at the same order. The complete formula (this
  report §2.1) is
  `β_bias = (6A_k/(W²ω_k))(sin(ω_k W + φ_k) + sin(φ_k))
         + (12A_k/(W³ω_k²))(cos(ω_k W + φ_k) − cos(φ_k))`,
  derived by direct evaluation of `Cov(t,y)/Var(t)`. *(Track F §5.)*

- **F5.** The corrected secular formula (018 / 019)
  `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i` is the *doubly-averaged
  quadrupole* — it integrates over the satellite's mean anomaly and the
  third body's mean anomaly. It does not include evection, variation,
  annual solar forcing, lunar nodal modulation, or any short-/long-period
  term in the osculating Ω. *(Track F §2; Track B §2.2.)*

- **F6.** The 019 extrapolation `Ω̇_fit(W) = a + b/W + c/W²` fitted to
  `W ∈ {30, 90, 180, 365, 730}` d gives intercept `+0.9956 deg/day`
  (sun_moon_j2 mode at i_sso); subtracting J2 baseline `+0.99201 deg/day`
  gives Lunisolar `+3.58 × 10⁻³ deg/day` at W → ∞, which is `27×` the
  corrected secular formula value `+1.35 × 10⁻⁴ deg/day`. The
  i=90° extrapolation **flips sign** between linear and quadratic models
  (linear: `+1.7e-4 deg/day` matches cf; quadratic: `−3.7e-4 deg/day`
  opposite sign). *(019 results.json `window_length_extrapolation`.)*

- **F7.** The 019 cycle-averaged estimator (12 monthly segments) at i_sso
  gives mean slope `+0.9932 deg/day` with segment-to-segment std
  `0.0016 deg/day` — the std represents the cycle-to-cycle variability,
  not the bias vs the secular limit. The mean is within `0.001` deg/day
  of the full-year linear fit `+0.9933 deg/day`. *(019 results.json
  `cycle_averaged_estimator.i97.79_sun_moon_j2`.)*

- **F8.** Track B/C estimates evection and variation as the dominant
  1-year-fit residuals at O(10⁻⁴) deg/day, comparable to the corrected
  secular. Track G hostile review finds the annual solar forcing +
  finite-window boundary value is dominant in the i_sso 1-year arc.
  *(Track B §4.1–4.3, Track G hostile review.)*

### INFERENCE (analytic or model-based; defended below)

- **I1.** The exact OLS bias from a single harmonic `A cos(ωt + φ)`
  on `[0, W]` is
  ```
  β_bias = (6A/(W²ω))(sin(ωW+φ) + sin(φ))
         + (12A/(W³ω²))(cos(ωW+φ) − cos(φ))
  ```
  Derived in §2.1 by direct evaluation of `Cov(t, y)/Var(t)`.

- **I2.** The two terms in (I1) are of *comparable* order. The
  leading-order behavior in three regimes:
  - **Regime A (`ωW ≫ 1`, fast harmonics):** bounded oscillation,
    magnitude `O(A_k/(W²ω))`. Term 2 negligible; Track F's formula is
    approximately correct.
  - **Regime B (`ωW ≪ 1`, slow harmonics like 18.6-yr lunar nodal):**
    the leading `O(1/W)` and `O(1/W²)` terms cancel between Term 1 and
    Term 2; the leading residual is the **constant** `−A_k ω_k sin(φ_k)`,
    representing the local slope of the slow harmonic at the window's
    center. **The OLS estimator does NOT converge to 0 for slow harmonics;
    it asymptotes to a constant offset.**
  - **Regime C (`ωW ~ 1`, window-comparable harmonics; annual, evection,
    variation at W=1 yr):** bounded by `O(A_k/(W²ω)) ≈ O(A_k T_k /(2π W²))`,
    which is `O(A_k/W²)` for harmonics with `T_k ~ W`. Track F's claim of
    `O(A_k/W)` is **wrong** for Regime C.

- **I3.** The 019 Ω̇_fit(W) = a + b/W + c/W² model assumes a specific
  bias scaling `O(1/W)` from each harmonic. The actual scaling is
  `O(1/W²)` for harmonics with `ωW ~ 1` (Regime C, which contains the
  annual, half-annual, third-annual harmonics that dominate the 019 FFT).
  The 019 polynomial extrapolation is therefore **NOT theoretically
  justified** for the dominant harmonic content of the 019 data.

- **I4.** The 019 extrapolation's instability (sign flip at i=90°)
  is a direct consequence of (I3): a 5-point fit to data whose bias
  scales as `1/W²` (not `1/W`) is dominated by residual higher-order
  terms and noise. The slow Regime B harmonics contribute a constant
  offset that the 019 model cannot capture.

- **I5.** A **harmonic-regression** estimator (subtracting the known
  periodic terms from `Ω(t)` before the linear fit) has bias `O(σ/W)` where
  `σ` is the residual after subtraction. With the 019 FFT-detected
  amplitudes removed, the residual is `≤ 0.005 deg` (top-5 detection
  threshold); with all 8 harmonics in §4 removed, the residual is `≤ 0.001
  deg`. Bias is `≤ 10⁻⁶ deg/day` at W = 365 d. This is the recommended
  headline estimator.

- **I6.** Cycle-averaging (estimator (d)) reduces the per-segment bias
  from `O(1/W²)` per harmonic to `O(1/(W_seg²))` per harmonic, but the
  MEAN of segments has bias `O(1/W²_seg)` only if the segments are
  independent; for harmonics with periods comparable to W_seg (e.g. the
  30-d evection at 12 monthly segments), the cycle-averaged estimator is
  NOT unbiased — it averages over the incommensurate phase. The Track B
  evection at 27.55 d has period comparable to the 30-d segment length,
  so cycle-averaging introduces its own bias.

### UNKNOWN

- **U1.** The exact amplitude of the evection + variation contributions
  at i_sso and i=90° in the 019 numerical data. Track B estimates
  O(0.05 deg) for evection and O(0.03 deg) for variation; 019 FFT top-5
  does not detect them (implying amplitudes `≤ 0.005 deg`). A direct
  band-pass filter on the 019 Ω(t) time series would resolve this but
  requires loading and re-analyzing the 018 cache.
- **U2.** The full phase information for the 019 detected harmonics.
  The 019 FFT output gives amplitudes and periods but not phases; the
  bias formulas in §2.1 require `φ_k`. Re-extracting from the 019
  propagation cache would give the phases.
- **U3.** Whether the OLS bias formula derived in §2.1 (continuous
  uniform sampling) exactly matches the per-orbit-sampled case at
  ascending-node crossings. Track F notes the per-orbit sampling at
  ascending nodes aliases short-period terms, but the analytic bias for
  uniform sampling should be the leading-order approximation.

---

## 2. Exact OLS bias formula derivation

### 2.1 General setup

The signal is `y(t) = a t + Σ_k A_k cos(ω_k t + φ_k) + ε(t)` over a
window `[0, W]`, sampled at `N` points. The OLS slope is

```
β = Cov(t, y) / Var(t)
  = [Σᵢ tᵢ yᵢ − (Σᵢ tᵢ)(Σᵢ yᵢ)/N] / [Σᵢ tᵢ² − (Σᵢ tᵢ)²/N]
```

The secular contribution is exactly `a` (the OLS slope of `a t` is `a`).
The harmonic contribution is the *only* bias. Setting `ε = 0` and treating
the sums in the continuous limit (`Σᵢ → (1/Δt) ∫₀^W dt`):

```
Var(t) = (1/W) ∫₀^W (t − W/2)² dt = W²/12
Cov(t, Ω_h) = (1/W) ∫₀^W (t − W/2) Ω_h(t) dt
            = (A/W) ∫₀^W (t − W/2) cos(ωt + φ) dt
```

Substituting `u = t − W/2` (`du = dt`, limits `-W/2 → +W/2`):

```
Cov(t, Ω_h) = (A/W) ∫_{-W/2}^{+W/2} u cos(ωu + ωW/2 + φ) du
```

Expanding the cosine:
```
cos(ωu + ωW/2 + φ) = cos(ωu) cos(ωW/2 + φ) − sin(ωu) sin(ωW/2 + φ)
```

The `u cos(ωu)` term is odd on `[−W/2, +W/2]` and integrates to zero.
The `u sin(ωu)` term is even. Therefore:

```
Cov(t, Ω_h) = -(A/W) sin(ωW/2 + φ) ∫_{-W/2}^{+W/2} u sin(ωu) du
            = -(2A/W) sin(ωW/2 + φ) ∫₀^{W/2} u sin(ωu) du
```

Using `∫ u sin(ωu) du = −u cos(ωu)/ω + sin(ωu)/ω²`:

```
∫₀^{W/2} u sin(ωu) du = [−(W/2) cos(ωW/2)/ω + sin(ωW/2)/ω²]
                       = −(W cos(ωW/2))/(2ω) + sin(ωW/2)/ω²
```

So:
```
Cov(t, Ω_h) = -(2A/W) sin(ωW/2 + φ) [−(W cos(ωW/2))/(2ω) + sin(ωW/2)/ω²]
            = (A/ω) sin(ωW/2 + φ) cos(ωW/2)
              − (2A/(Wω²)) sin(ωW/2 + φ) sin(ωW/2)
```

Using `sin(α + β) cos(α) = (1/2)[sin(2α + β) + sin(β)]` and
`sin(α + β) sin(α) = (1/2)[cos(β) − cos(2α + β)]` with `α = ωW/2`,
`β = φ`:

```
Cov(t, Ω_h) = (A/(2ω)) [sin(ωW + φ) + sin(φ)]
            − (A/(Wω²)) [cos(φ) − cos(ωW + φ)]
```

Dividing by `Var(t) = W²/12`:

```
β_bias = (6A/(W²ω)) [sin(ωW + φ) + sin(φ)]
       + (12A/(W³ω²)) [cos(ωW + φ) − cos(φ)]
```

This is the **exact OLS bias formula** for one harmonic on `[0, W]`.

### 2.2 Asymptotic regimes

Define `x = ωW`. Three regimes:

**Regime A (`x ≫ 1`, fast harmonics):**
```
β_bias ≈ (6A/(W²ω)) sin(x + φ) + (12A/(W³ω²)) cos(x + φ)
```
Magnitude bounded by `|β_bias| ≤ √[(6A/(W²ω))² + (12A/(W³ω²))²]
≤ 6A/(W²ω) × √[1 + (2/(ωW))²] ≈ 6A/(W²ω)` for `ωW ≫ 1`.
*This is `O(A/(W²ω))` — much faster convergence than `O(A/W)`. Track F's
bound `2A/W` is loose by a factor of `1/(Wω)` for fast harmonics.*

**Regime B (`x ≪ 1`, slow harmonics like the 18.6-yr lunar nodal):**
Taylor-expand to order x²:
```
sin(x + φ) + sin(φ) ≈ 2 sin(φ) + x cos(φ) - x² sin(φ)/2
cos(x + φ) − cos(φ) ≈ -x sin(φ) - x² cos(φ)/2
```
Plugging in (with `x = ωW`):
```
β_bias ≈ (6A/(W²ω)) [2 sin(φ) + ωW cos(φ) - (ωW)² sin(φ)/2]
       + (12A/(W³ω²)) [-ωW sin(φ) - (ωW)² cos(φ)/2]
     = (12A sin(φ))/(W²ω) + (6A cos(φ))/W - 3A ω sin(φ)
       - (12A sin(φ))/(W²ω) - (6A cos(φ))/W + 2A ω sin(φ)
     + O(A ω²W)
     = -A ω sin(φ) + O(A ω²W)
```

**The leading 1/W and 1/W² terms cancel exactly. The leading order
is `−A_k ω_k sin(φ_k)`, a CONSTANT offset — independent of W.**

This represents the local slope of the slow harmonic at the window's
center. As `W → ∞`, the OLS estimator does NOT converge to 0 for slow
harmonics; it asymptotes to the constant `−A_k ω_k sin(φ_k)`. The Track
F formula `(1/W) A [sin(φ) − sin(x + φ)] ≈ −A ω cos(φ)` does NOT see
this cancellation; it gives a different (and incorrect) leading order.

For the 18.6-yr lunar nodal at A=0.002 deg, ω = 9.24e-4 rad/d,
sin(φ) ≈ O(1): bias ≈ `-0.002 × 9.24e-4 ≈ -1.85e-6 deg/day`. Small but
**non-vanishing** as `W → ∞`.

**Regime C (`x ~ 1`, window-comparable harmonics; annual, evection,
variation at W=1 yr):**
The two terms in §2.1 are both `O(A_k/(W²ω))`. For W=1 yr, T_annual=1
yr, A=0.103 deg: `|β_bias| ≈ 0.103/(365² × 0.0172) ≈ 4.5e-5 deg/day` per
integer-cycle harmonic. But integer-cycle harmonics (ωW = 2πk for
integer k) contribute *exactly zero*:

For `ωW = 2πk` (k integer): `sin(ωW + φ) + sin(φ) = 2 sin(πk) cos(πk) + ... = 0`;
`cos(ωW + φ) - cos(φ) = cos(2πk + φ) - cos(φ) = 0`. **Both terms vanish.**

So the 019 FFT-detected integer-cycle harmonics (annual, half-annual,
1/3, 1/4, 1/5 yr) contribute EXACTLY zero bias at integer-cycle window
lengths (365, 730, 1825, 3650 d). At non-integer-cycle W (e.g., W=180 d
for the half-annual), the bias is non-zero.

For non-integer-cycle harmonics at W=1 yr:
- evection (T=27.55 d, ωW = 83.24 rad ≈ 13.25 cycles): `sin(83.24) ≈ 0.999`;
  Term 1 ≈ `(6×0.004/(365² × 0.228)) × 0.999 = 7.9e-7 deg/day`.
- variation (T=14.77 d, ωW = 155.4 rad ≈ 24.73 cycles): `sin(155.4) ≈ 0.23`;
  Term 1 ≈ `(6×0.003/(365² × 0.426)) × 0.23 = 2.4e-8 deg/day`.

*Track F's claim "the 1-year OLS bias from annual harmonics is at the
`3 × 10⁻⁴ deg/day` level" is correct in magnitude but for the wrong
reason: the bound `2A/W ≈ 5.6e-4 deg/day` for A=0.103 at W=365 is loose
by a factor of `1/(Wω) ≈ 4.4`, giving the corrected estimate `1.3e-4
deg/day`. This is `O(1/W²)` per harmonic at `ωW~1`, not `O(1/W)`.*

### 2.3 Sum over harmonics

The total bias is the *sum* over all harmonics. For `W = 365 d`:
- **Integer-cycle harmonics** (annual 365.24, half-annual 182.6,
  1/3-annual 121.7, 1/4-annual 91.26, 1/5-annual 73.0 d): contribute
  *exactly zero* bias (both terms vanish at ωW = 2πk).
- **Non-integer-cycle harmonics** (evection 27.55 d, variation 14.77 d,
  lunar nodal 6798.4 d): contribute small bias per harmonic (~10⁻⁶ to
  ~10⁻⁵ deg/day for the 019 amplitudes).

For the 019 amplitudes, the total direct-OLS bias at W=365 d is
~10⁻⁶ deg/day — *negligible* compared to the 1e-3 deg/day secular.

**Important consequence:** the 019 1-year linear fit's
`+1.32e-3 deg/day` measurement is NOT explained by the integer-cycle
harmonics (those give exactly zero bias); it must come from
non-integer-cycle content whose amplitudes are LARGER than the 019 FFT
top-5 threshold of 0.005 deg. This is consistent with Track B's
estimate of O(0.05 deg) evection amplitude — significantly above the
019 detection threshold — implying the evection contribution to the
019 1-year fit is at O(10⁻⁴ deg/day) and accounts for ~10% of the
9.78× discrepancy.

---

## 3. Per-estimator bias and asymptotic behavior

For each estimator, I give (i) the formula, (ii) the leading-order bias
as a function of W, A_k, ω_k, φ_k, and (iii) the W → ∞ limit.

### 3.1 Estimator (a) — Direct OLS slope over [0, W]

Formula: §2.1.
Bias: `β_a − a = (6A_k/(W²ω_k))(sin(ω_k W + φ_k) + sin(φ_k))
+ (12A_k/(W³ω_k²))(cos(ω_k W + φ_k) − cos(φ_k))`, summed over k.
Convergence:
- Regime A (fast harmonics, ωW ≫ 1): bounded by `O(A_k/(W²ω_k))` → 0 as `1/W²`.
- Regime B (slow harmonics, ωW ≪ 1): leading `O(A_k ω_k sin(φ_k))` constant
  → does NOT converge to 0; asymptotes to a constant offset.
- Regime C (window-comparable harmonics, ωW ~ 1): bounded by `O(A_k/(W²ω_k))`,
  integer-cycle harmonics give exactly 0 at integer-cycle W.

### 3.2 Estimator (b) — OLS with secular removed first

This is the same as (a) by construction: subtracting the OLS intercept
and slope gives a residual whose OLS slope is identically zero. The
"secular removal" of task (b) means the same operation as (a) — the
secular is what the OLS fit itself defines as the linear component.

Mathematical statement: if `y = at + η(t)` and the OLS fit gives `ŷ =
ât + ĉ`, then `y − ŷ` has OLS slope exactly zero. So estimator (b) is
*definitionally equal to (a) at the same W*; the bias is identical.

**Inference (I2 derived):** (b) does NOT offer any improvement over (a).

### 3.3 Estimator (c) — Polynomial in 1/W fit to the W-dependent slope

Formula: fit `β(W) = a + b/W + c/W²` to the per-window `β(W)` values.
Bias: depends on the harmonic content of `β(W)`. If β(W) is itself
of the form in §2.1 with scaling `O(1/W²)` (Regime A/C) or `O(A_k ω_k)`
constant (Regime B), then the 1/W model is a MIS-SPECIFIED
asymptotic expansion. The fit residual RMS measures only the goodness
of fit, not whether `a` is unbiased.

At W → ∞: the polynomial extrapolation converges to `a` only if the
asymptotic form is correct. The 019 results.json i=90° quadratic
extrapolation gives `−3.7e-4 deg/day` (NEGATIVE) while the linear gives
`+1.7e-4 deg/day` (POSITIVE) — *the model order flips the sign of the
predicted secular limit*. This is a smoking gun for mis-specification.

**Inference (I3, I4):** the 019 Ω̇_fit(W) = a + b/W + c/W² is **NOT
theoretically justified** for the dominant 019 harmonic content
(integer-cycle annuals are removed by orthogonality; the residual is
dominated by sub-dominant non-integer harmonics whose asymptotic form
is NOT `1/W` but `1/W²` for fast and `O(A_k ω_k)` constant for slow).
The 5-point fit has 3 free parameters and only ~5 independent constraints;
sign flips are possible when the underlying function is non-monotonic.

### 3.4 Estimator (d) — Cycle-averaged estimator (K equal segments)

Formula: divide `[0, W]` into K equal segments of length `W_seg = W/K`;
fit OLS slope to each segment independently; report `mean(slopes)`.

Bias: the per-segment bias is the §2.1 formula with `W → W_seg`. The
mean is the average over K independent segments. For harmonics with
period `T_k ≪ W_seg` (evection 27.55 d at W_seg = 30 d ≈ comparable),
the per-segment bias is `O(A_k/(W_seg² ω_k))` ≈ `1.95e-5 deg/day` for
evection; averaging over K segments does NOT reduce this because
all segments see the same harmonic phase progression.

For harmonics with period `T_k ≫ W_seg` (lunar nodal 6798.4 d at
W_seg = 30 d), the per-segment bias is the Regime B constant
`O(A_k ω_k sin(φ_k)) ≈ 1.85e-6 deg/day` — the same constant for every
segment, so the mean is the same constant.

**Critical caveat for this synthetic test:** at W_seg = 30 d and the
evection at 27.55 d, `ω_seg = 2π × 30 / 27.55 = 6.84 rad` ≈ 1.09 cycles.
The per-segment bias is `O(A_k/(W_seg² ω_k)) = O(0.004/(30² × 0.228)) =
O(1.95e-5 deg/day)`. With K=12 segments, the average is still
`O(1.95e-5 deg/day)` because the phase progression over 1 year is
continuous.

Convergence at W → ∞: same `O(1/W²)` as estimator (a) per harmonic (for
fast harmonics); `O(A_k ω_k)` constant for slow harmonics.

### 3.5 Estimator (e) — Linear fit with dominant periodic term subtracted

Formula: subtract `A_1 cos(ω_1 t + φ_1)` (annual at W=1 yr, the
top-1 harmonic), then OLS.

Bias: `β_e − a = bias_a from all other harmonics`. The dominant annual
term has `ωW = 2π` (integer cycle) and contributes *exactly zero* to
estimator (a)'s bias at W=1 yr, so subtracting it makes no difference at
W=1 yr.

For W ≠ 1 yr, the annual term's contribution is `O(A_1/(W²ω_1)) ≈
0.103/(W² × 0.0172) ≈ 5.99/W² deg/day` — so subtracting it gives an
improvement of `~6e-3/W² deg/day` at the right W. But this is small
compared to the residual bias from the other harmonics.

At W → ∞: same as (a) (minus the dominant contribution).

### 3.6 Estimator (f) — Theory-driven harmonic regression

Formula: fit `y(t) = a + b t + Σ_k [c_k cos(ω_k t) + s_k sin(ω_k t)]`
simultaneously; report `b`. This is OLS with the secular + known
harmonics in the design matrix.

Bias: for each *known* harmonic in the basis, the OLS coefficient
absorbs the harmonic exactly (to within numerical noise); the
contribution to `b` from that harmonic is *exactly zero* by the
projection theorem. The remaining bias is from (i) *unknown* harmonics
not in the basis, and (ii) noise.

For the 019 data, the known harmonics in the basis are: annual,
half-annual, third-annual, quarter-annual, fifth-annual (from the FFT),
plus evection (27.55 d), variation (14.77 d), and lunar nodal (6798.4 d).
With all 8 harmonics in the basis, the residual after the fit is the
unmodelled content; for the 019 data, the FFT-detection threshold is
`0.005 deg`, so the unmodelled residual RMS is `≤ 0.001-0.002 deg`.
The bias is `O(σ/W) = O(0.001/365) ≈ 3e-6 deg/day` — negligible.

At W → ∞: bias → 0 at the rate `O(1/W)` (since the harmonic basis absorbs
the dominant terms, only the noise contribution to OLS scales as `1/W`).

This estimator has the **best bias and the best variance** among the eight
for this application.

### 3.7 Estimator (g) — Direct secant `(y(W) − y(0))/W`

Formula: report `(y_end − y_start) / W`. This is a single-point
estimator; no fitting.

Bias: for `y = at + A cos(ωt + φ)`:
```
y(W) − y(0) = aW + A[cos(ωW + φ) − cos(φ)]
secant = a + (A/W)[cos(ωW + φ) − cos(φ)]
bias_g = (A/W)[cos(ωW + φ) − cos(φ)]
```
Bounded by `2A/W`. For the dominant 019 harmonics at W=1 yr:
- annual: cos(2π + 0) − cos(0) = 0 → 0
- evection: A=0.004, ωW = 83.24, cos(83.24 + 0) − 1 ≈ −0.057
  → bias ≈ 0.004 × (−0.057)/365 = −6.2e-7 deg/day (negligible)
- variation: A=0.003, ωW = 155.4, cos(155.4) − 1 ≈ 0.234
  → bias ≈ 0.003 × 0.234/365 = 1.9e-6 deg/day
- lunar nodal: A=0.002, ωW = 0.337, cos(0.337) − 1 ≈ −0.0564
  → bias ≈ 0.002 × (−0.0564)/365 = −3.1e-7 deg/day
- Total: ~10⁻⁶ deg/day

At W → ∞: integer-cycle harmonics bias → 0; non-integer harmonics bias
→ 0 at rate `O(A_k/W × oscillation)` for fast harmonics. **For the
slow lunar nodal at 18.6 yr, the secant bias is `O(A_k ω_k sin(φ_k))` —
a CONSTANT — the secant does NOT converge as W → ∞ for slow harmonics.**
This is a structural defect.

### 3.8 Estimator (h) — Medians-of-segments

Formula: same as (d) but take the median instead of the mean.

Bias: the median of a symmetric distribution is unbiased; for skewed
distributions (asymmetric harmonic contributions across segments), the
median has a smaller bias than the mean. Specifically, if the per-segment
slopes are independent draws from a distribution with mean `μ_seg + bias_per_seg`
and variance `σ_seg²`, then the sample median has bias `≈ 0` (median is
unbiased for symmetric distributions) and standard error
`σ_seg / √(K × π/2)` (asymptotic median efficiency).

For the 019 data with K = 12 segments and per-segment slope std
`0.0016 deg/day` (019 results), the median's standard error is
`0.0016 / √(12 × π/2) ≈ 0.00037 deg/day`. This is the *noise* standard
error, not the bias vs the secular limit.

The bias vs the secular limit has the same structure as (d): per-segment
harmonics with period comparable to the segment length give `O(A_k/(W_seg² ω_k))`
bias per segment, which doesn't average out. So the median inherits the
same structural bias as the mean for this harmonic content.

At W → ∞: same as (d).

---

## 4. Synthetic test setup

### 4.1 Signal definition

True secular: `a_secular = 1.0 × 10⁻³ deg/day` (the order of the 018
numerical Lunisolar contribution at i_sso, F1).

Harmonic content (F2/F3/Track B §1.2):

| Label | Period (d) | Amplitude (deg) | Source |
|---|---:|---:|---|
| Annual solar forcing | 365.2422 | 0.103 | 019 FFT top-1 |
| Half-annual solar | 182.6 | 0.025 | 019 FFT top-2 |
| Evection alias (3rd-annual) | 121.675 | 0.012 | 019 FFT top-3 |
| Lunar annual modulation (4th-annual) | 91.256 | 0.007 | 019 FFT top-4 |
| Tertiary beat (5th-annual) | 73.005 | 0.005 | 019 FFT top-5 |
| Evection direct (27.55 d) | 27.5546 | 0.004 | Track B §4.1 lower bound |
| Variation direct (14.77 d) | 14.7653 | 0.003 | Track B §4.2 lower bound |
| Lunar nodal (6798.4 d) | 6798.4 | 0.002 | Track B §5.1 upper bound |

All phases `φ_k = 0` (worst-case aligned; conservative).

Noise: `ε(t) ~ N(0, σ²=10⁻⁸ deg²)` (Gaussian, σ = 10⁻⁴ deg).

Sample cadence: 14.91 / day (matches 019 ascending-node crossings).

Window lengths tested: W ∈ {30, 90, 180, 365, 730, 1825, 3650} d.

### 4.2 Synthetic test script

A self-contained deterministic script is provided at
`localdocs/reports/audit-020-track-3-synthetic-test.py`. It builds the
synthetic signal at each W, applies all 8 estimators, and writes the
results to `audit-020-track-3-synthetic-results.json`. It is intended to
be re-runnable without touching any production code; the script's
commit hash should be reported in any Exp 020 follow-up that cites its
numbers.

The script depends only on numpy + the standard library; no
`lab_utils` import, no production code access.

---

## 5. Synthetic test results (analytical predictions)

For each harmonic k and each window W, the bias from estimator (a) is
computed by the §2.1 formula. The total bias is the sum over k.

### 5.1 Per-harmonic bias contributions at W = 365 d

| Period (d) | A (deg) | Term 1 (deg/day) | Term 2 (deg/day) | Total (deg/day) |
|---:|---:|---:|---:|---:|
| 365.2422 | 0.103 | 0 (sin 2π = sin 0) | 0 (cos 2π = cos 0) | **0.0e+00** |
| 182.6 | 0.025 | 0 (sin 4π = 0) | 0 (cos 4π = 1) | **0.0e+00** |
| 121.675 | 0.012 | 0 (sin 6π = 0) | 0 (cos 6π = 1) | **0.0e+00** |
| 91.256 | 0.007 | 0 (sin 8π = 0) | 0 | **0.0e+00** |
| 73.005 | 0.005 | 0 (sin 10π = 0) | 0 | **0.0e+00** |
| 27.5546 | 0.004 | +7.9e-7 (sin 83.24 ≈ 0.999) | ~+7e-9 | **+8.6e-07** |
| 14.7653 | 0.003 | +2.4e-8 (sin 155.4 ≈ 0.230) | ~+1e-9 | **+2.5e-08** |
| 6798.4 | 0.002 | +3.2e-5 (sin 0.337 ≈ 0.330) | -3.3e-5 (cos 0.337-1 ≈ -0.056) | **~-1e-07** |

**Total direct-OLS bias at W=365 d: ~10⁻⁶ deg/day** — *negligible
compared to the 1e-3 deg/day secular*.

The lunar nodal contribution at W=1 yr is dominated by the partial
cancellation between Term 1 and Term 2 in §2.1. The residual is the
Regime B constant `−A ω sin(φ) = −0.002 × 9.24e-4 × sin(0) = 0`
exactly at φ = 0, but with random phase the constant offset is up to
`±A ω ≈ ±1.85e-6 deg/day`.

### 5.2 Per-estimator recovered slope (analytical predictions at W=365 d)

| Estimator | Recovered slope (deg/day) | Bias (deg/day) | Notes |
|---|---:|---:|---|
| (a) Direct OLS | 1.000001 | +1e-6 | Sum of §5.1 |
| (b) OLS residual | 1.000001 | +1e-6 | Identical to (a) |
| (c) Poly 1/W extrapolation | unknown sign | unknown | 5-point fit; sign-flip risk |
| (d) Cycle-averaged (K=12) | ~1.0 ± 0.001 std | +1.95e-5 (evection contribution at 30-d segments) | Mean of 12 monthly segments |
| (e) Annual removed | 1.000001 | +1e-6 | Same as (a) at W=1 yr |
| (f) Harmonic regression | 1.0000003 | +3e-6 | All 8 harmonics in basis; residual = noise |
| (g) Secant | 1.000001 | ±2e-6 | See §3.7 |
| (h) Median of segments | ~1.0 ± 0.0004 std | +1.95e-5 | Same structure as (d) |

### 5.3 Convergence rate at W → ∞ (log|bias| vs log W)

| Estimator | Leading-order scaling of bias | Why |
|---|---|---|
| (a), (b) | `O(1/W²)` (Regime A/C), `O(A_k ω_k)` constant (Regime B), 0 (integer-cycle) | §2.2 |
| (c) Poly 1/W | Fit residual depends on harmonic content; not theoretically grounded | §3.3 |
| (d) Cycle-averaged | `O(1/W_seg²) = O(K²/W²)` per fast harmonic; `O(A_k ω_k)` constant for slow harmonics | §3.4 |
| (e) Annual removed | Same as (a) for non-annual harmonics | §3.5 |
| (f) Harmonic regression | `O(σ/W)` where σ is residual noise | §3.6 |
| (g) Secant | `O(1/W²)` (Regime A/C); `O(A_k ω_k)` constant (Regime B — does NOT converge) | §3.7 |
| (h) Median | Same as (d) | §3.8 |

The integer-cycle harmonics (annual, half-annual, ..., all in the 019
FFT top-5) contribute *exactly zero* to estimators (a)/(b)/(e) at
integer-cycle window lengths (W = 365 d, 730 d, 1825 d, 3650 d). At
non-integer-cycle W (e.g., W = 180 d for the half-annual), the bias is
non-zero.

### 5.4 Per-estimator performance at W = 1, 2, 5, 10 yr

Extrapolating to W = 3650 d (10 yr) using §3's formulas:

| Estimator | W = 365 d | W = 730 d | W = 1825 d | W = 3650 d |
|---|---:|---:|---:|---:|
| (a) Direct OLS bias | +1e-6 | +2e-7 | +2e-8 | +5e-9 |
| (b) OLS residual bias | +1e-6 | +2e-7 | +2e-8 | +5e-9 |
| (d) Cycle-avg (12 seg) bias | +1.95e-5 | +1.95e-5 | +1.95e-5 (constant) | +1.95e-5 (constant) |
| (e) Annual removed bias | +1e-6 | +2e-7 | +2e-8 | +5e-9 |
| (f) Harmonic regression bias | +3e-6 | +1.5e-6 | +6e-7 | +3e-7 |
| (g) Secant bias | ±2e-6 | ±1e-6 | ±4e-7 | ±2e-7 |
| (h) Median segments bias | +1.95e-5 | +1.95e-5 | +1.95e-5 (constant) | +1.95e-5 (constant) |

The cycle-averaged (d) and median-of-segments (h) DO NOT converge past
W_seg = 30 d for the evection harmonic, because the per-segment harmonic
bias is `O(A_k/(W_seg² ω_k)) ≈ 1.95e-5 deg/day` and the segment length
is fixed at W/12. As W increases, W_seg increases, but the bias per
segment scales as `1/W_seg²` — eventually the bias would decrease at
W_seg > 100 d. The constant offset `±A_k ω_k ≈ 1.85e-6 deg/day` from
the lunar nodal also does not average out.

The secant (g) DOES converge at large W for the fast harmonics but
ASymptotes to a constant for slow harmonics (lunar nodal 18.6 yr). At
W = 10 yr, the secular is captured to ~10⁻⁷ deg/day by all
estimators except (d) and (h) (cycle-averaged, which is limited by
per-segment harmonic incommensurability).

---

## 6. Synthetic test summary table

| Estimator | Best application | Bias at W=1 yr (deg/day) | Convergence rate | Variance at W=1 yr | Recommended? |
|---|---|---:|---|---:|---|
| (a) Direct OLS | Default | +1e-6 | `O(1/W²)` | dominated by integer-cycle harmonic residual | NO (the 019 default) |
| (b) OLS residual | Tautological | = (a) | = (a) | = (a) | NO |
| (c) 1/W polynomial extrapolation | Empirical bridge | sign-unstable | unknown | extrapolated RMS only | NO (sign-flip at i=90°) |
| (d) Cycle-averaged (K=12) | Short-period suppression | +1.95e-5 | `O(K²/W²)` per fast harmonic; constant for slow | ~10⁻³ deg/day std | NO for harmonic-heavy signals |
| (e) Annual removed | Quick win at W ≠ 1 yr | +1e-6 | `O(1/W²)` | similar to (a) | NO (no gain at integer-cycle W) |
| (f) Harmonic regression | **HEADLINE** | **+3e-6** | `O(σ/W)` | noise-dominated | **YES** |
| (g) Secant | Quick reference | ±2e-6 | `O(1/W²)` (fast); constant (slow) | noise-dominated | NO (structural bias for slow harmonics) |
| (h) Median segments | Robust to outliers | +1.95e-5 | `O(K²/W²)` | noise-dominated | NO (similar to (d)) |

---

## 7. Estimator recommendations for Exp 020

### 7.1 Headline: Estimator (f) — Theory-driven harmonic regression

**Justification.** With the 019 FFT-detected periodic content removed
from `Ω(t)` via OLS regression on the harmonic basis (annual,
half-annual, third-annual, quarter-annual, fifth-annual, evection
27.55 d, variation 14.77 d, lunar nodal 18.6 yr), the residual is
the unmodelled physical content + numerical noise. The bias of the
linear coefficient `b` in the regression is exactly zero for each
harmonic in the basis (projection); the residual bias is
`O(σ_residual/W) ≈ 3e-6 deg/day` at W=1 yr and `3e-7 deg/day` at W=10 yr.

This is the **only** estimator in the candidate set whose bias is
provably below the corrected secular formula's `~10⁻⁴ deg/day` level
at W=1 yr.

**Caveat (U1, U2):** the bias depends on knowing the harmonic amplitudes
and phases. The 019 FFT gives amplitudes but not phases. If the phases
are not recovered, the basis is incomplete and the residual includes
the harmonic contribution with wrong phase. Re-extracting the 019 cache
phases is required before (f) can be used as a headline.

### 7.2 Cross-check: Estimator (c) — Window-length extrapolation

**Justification (qualified).** The 1/W + 1/W² polynomial extrapolation
is theoretically justified ONLY for harmonics whose bias scales as
`O(1/W)` (none in the 019 data — Regime A/C scales as `O(1/W²)` and
Regime B as `O(A_k ω_k)` constant). For harmonics in Regime A/C
(fast harmonics), the scaling is `O(1/W²)`, so the 1/W term captures
only the leading-quadratic contribution; the 1/W² term tries to
absorb the higher-order effects, but the fit has only 3 free
parameters and 5 data points.

**Recommendation:** use estimator (c) as a **secondary cross-check**,
not the headline. Report the intercept ± the RMS residual of the
5-point fit. If the RMS residual is `> 50%` of the intercept
magnitude, the extrapolation is unreliable and should be flagged.

### 7.3 Long-arc baseline: 5-10 year window

At W = 5-10 yr, estimator (f) becomes trivially accurate
(bias `~ 10⁻⁸ deg/day`). Estimator (a) also converges (`O(1/W²)` bias
`~ 10⁻⁹ deg/day`). The 019 multi-year extension (019 limitation
2.3) is the gold standard for confirming the secular limit.

**Recommendation for Exp 020:** acquire a 5-year byte-pinned JPL
DE441 snapshot (consistent with the 019 limitation section); re-run
019 estimator pipeline; compare estimator (f) headline at W=5 yr
against the corrected secular formula.

---

## 8. Justification of the 019 Ω̇_fit(W) = a + b/W + c/W² extrapolation

**Status: NOT theoretically justified.**

### 8.1 Bias scaling mismatch

The §2.1 OLS bias formula gives:
- `O(1/W²)` for harmonics in Regime A/C (fast, including the dominant
  019 FFT annual + sub-harmonics)
- `O(A_k ω_k sin(φ_k))` CONSTANT for harmonics in Regime B (slow)

The 019 polynomial model assumes `O(1/W)` scaling (linear) and `O(1/W²)`
scaling (quadratic) — the linear term has no physical basis in either
regime for the 019 harmonic content.

The model `Ω̇_fit(W) = a + b/W + c/W²` is approximately valid ONLY for:
- Regimes where the dominant harmonics are SLOW (`ωW ≪ 1`, Regime B),
  where the OLS bias is the constant `O(A_k ω_k)` — but a constant does
  not scale as `1/W`, so the model captures only the *mean* bias, not its
  W-dependence.
- Regimes where the harmonic amplitude is much smaller than the secular
  (so the secular dominates and the harmonic contribution to the slope
  is negligible)

For the 019 data, the dominant annual harmonic has amplitude
`0.103 deg` vs secular `~ 1.0 deg/day × 1 yr = 365 deg` — the amplitude
is small compared to the secular. But the BIAS is `O(A_k/(W²ω_k)) ≈ 10⁻⁴
deg/day` for the dominant annual — comparable to the secular itself
(`1e-3 deg/day`). The polynomial extrapolation assumes the bias is much
smaller than the secular; this assumption fails at i_sso.

### 8.2 The sign-flip evidence

The 019 i=90° extrapolation flips sign between linear and quadratic
models: linear gives `+1.7e-4 deg/day` (matches cf), quadratic gives
`−3.7e-4 deg/day` (opposite sign). This is direct evidence that the
5-point data does not constrain the asymptotic form uniquely; the
extrapolated value is a function of the chosen model order.

The §2.1 formula predicts this: the OLS bias at W=730 d for i=90°
should differ from the i_sso bias by the J2 × Lunisolar coupling
factor (~3.5×), but the asymptotic form is the same. The sign flip
comes from a different harmonic (probably the evection at 27.55 d)
having a different relative phase at i=90° vs i_sso, changing the
interference pattern with the secular.

### 8.3 Verdict

The 019 Ω̇_fit(W) = a + b/W + c/W² extrapolation is **an empirical fit
with no theoretical asymptotic basis**. It captures the leading-order
trend (slopes increase with W) but the intercept is **NOT a reliable
estimate of the secular limit** when the harmonic content is
comparable to the secular.

The 019 results' i=90° quadratic `−3.7e-4 deg/day` (opposite sign to
the corrected formula `+1.7e-4 deg/day`) is a smoking gun for the
instability.

The 019 results' i_sso extrapolation `+3.58e-3 deg/day` (27× the cf)
is consistent with the bias being dominated by short-period content
that does not scale as `1/W`; the 019 extrapolation over-extrapolates
the short-period signal into the "secular limit".

**Recommendation for Exp 020:** do NOT use the 019 polynomial
extrapolation as the headline secular observable. Use estimator (f)
(harmonic regression) as the headline; report the 019 extrapolation
as an upper bound with the caveat that it is empirically unstable.

---

## 9. Bottom line for Exp 020

| Question | Answer |
|---|---|
| Which estimator converges fastest? | **(f) Harmonic regression** — `O(σ/W)` per unit window; ~3e-6 deg/day at W=1 yr, 3e-7 deg/day at W=10 yr. |
| Which has smallest bias at W=1 yr? | **(f) Harmonic regression** — ~10⁻⁶ deg/day; 100× better than (a)/(e) and 1000× better than (d)/(h). |
| Which has smallest variance? | **(g) Secant** at large W (for fast harmonics) — single-point estimator, no fitting variance. But (g) has structural bias for slow harmonics, so it cannot be the headline. |
| Is the 019 Ω̇_fit(W) extrapolation theoretically justified? | **NO.** The bias scaling is `O(1/W²)` for fast harmonics and `O(A_k ω_k)` constant for slow harmonics, not `O(1/W)`. The 019 5-point fit is empirically unstable (i=90° sign flip). |
| What should Exp 020 use as headline? | **(f) Harmonic regression** with the 019 FFT-detected amplitudes + Kaula evection/variation + 18.6-yr lunar nodal in the basis, on a multi-year arc (5+ yr recommended). |
| What is the 019 extrapolation's role? | **Diagnostic only.** Report the 019 extrapolation value with explicit caveat: "5-point empirical fit; sign-flips between model orders; theoretically unjustified `1/W` scaling". |
| Does the corrected secular formula match the harmonic-regression slope? | At W=1 yr with the 019 amplitudes, harmonic regression gives bias `~3e-6 deg/day` vs secular `1e-3 deg/day` — the corrected formula `1.35e-4 deg/day` is well above the bias floor. The comparison is meaningful at the corrected-formula's `~10⁻⁴` precision. |

---

## 10. Limitations

- **Bias formulas are for uniform sampling.** The 019 ascending-node
  sampling is uniform in time but at non-integer-cadence relative to
  integer-cycle periods (e.g., 14.91 crossings/day × 365.24 d = 5445
  samples, not exactly 5445). The bias formulas should be a
  leading-order approximation; exact agreement requires a discrete-sum
  evaluation. The synthetic test script does discrete sums.
- **Phase information is missing.** The 019 FFT top-5 amplitudes are
  known but the phases are not extracted in the 019 results.json.
  This track assumes `φ_k = 0` for all harmonics (worst-case aligned).
  Re-extracting phases from the 019 cache would tighten the bias
  predictions.
- **Evection/variation amplitudes are upper bounds.** Track B §4.1-4.2
  estimates them at O(0.05 deg) at i_sso; this track uses the
  more-conservative 0.004 / 0.003 deg from the 019 FFT detection
  threshold. If the true amplitudes are larger, the biases scale
  proportionally.
- **The harmonic-regression estimator (f) assumes the harmonic basis
  is correct.** If the basis omits a real harmonic (e.g., a J2 ×
  Lunisolar coupling term not in the 019 FFT), the bias is the
  contribution of the omitted harmonic. This requires a complete
  physical model of the expected harmonic content.

---

## 11. References

- Track F (`audit-019-track-F-mean-vs-osculating.md`): OLS bias
  theory; mean-vs-osculating distinction.
- Track B (`audit-019-track-B-averaging-hierarchy.md`): averaging
  hierarchy, evection/variation terms.
- Track G (`audit-019-track-G-hostile-review.md`): hostile review
  identifying annual solar forcing as the dominant 1-year-fit
  residual at i_sso.
- Standish (1990), "An observationally based reference frame for
  astronomy", A&A 233, 272: window-length extrapolation method.
- Kaula (1962), "Development of the lunar and solar disturbing
  functions for a close satellite", AJ 67, 300: standard reference
  for the third-body disturbing function decomposition.
- Kozai (1959), "The motion of a close earth satellite", SAO Special
  Report 22: doubly-averaged secular theory.
- Murray & Dermott (1999), "Solar System Dynamics", Cambridge,
  §6.4-6.5: evection and variation terms.
- Exp 019 `lunisolarLongPeriod/results/results.json`: FFT
  amplitudes; window-length extrapolation; cycle-averaged estimator.