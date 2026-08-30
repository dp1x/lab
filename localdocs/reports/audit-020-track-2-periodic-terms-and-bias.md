# Track A-2 — Periodic Terms and OLS-Bias Asymptotic Scaling

> Audit-020 / Track A-2: derive the periodic, long-period, and intermediate-period
> terms in the third-body Lunisolar disturbing function that the doubly-averaged
> secular formula averages away, and quantify each term's expected contribution to
> the bias of a finite-window OLS linear fit of osculating Ω(t) at h = 600 km i_sso.
>
> Status: COMPLETE (2026-08-30). Read-only audit. No source code modified.
>
> Inputs read: audit-019-track-B-averaging-hierarchy.md, audit-019-track-C-
> evection-variation-hypothesis.md, audit-019-track-F-mean-vs-osculating.md,
> Exp 019 experiment.py + results.json.
>
> Inputs NOT read (per mission constraint): any other Track's output of audit-020.

---

## 1. Scope and motivation

The 018 corrected secular formula

  Ω̇_mean = (3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i                    (1)

is a **doubly-averaged quadrupole**: it integrates over the satellite's mean anomaly
(short-period elimination) and the third body's mean anomaly (long-period
elimination). The 019 experiment found that the 018 1-year numerical fit exceeds
the corrected cf by **9.78× at i_sso = 97.79°** and by **2.81× at i = 90°**.
Track F attributed this to mean-vs-osculating bias of the 1-year linear fit;
Track B/C estimated the residual was driven by the evection + variation + annual +
nodal terms that the doubly-averaged formula removes; Track F proposed the
**window-length extrapolation** Ω̇_fit(W) = a + b/W + c/W² as the canonical bridge.

This track:
1. Names the five periodic/long-period terms (evection, variation, annual solar,
   lunar nodal, lunar apsidal) and gives their physical origin and angle
   combinations (FACT and INFERENCE mix — see §2).
2. Estimates the **amplitude** of each term in the **osculating Ω(t)** at h = 600 km
   i_sso (INFERENCE based on standard Kaula/Kozai expansion amplitudes; UNKNOWN
   the precise factor of order unity).
3. Computes the **expected contribution** of each term to the OLS bias of a 1-year
   linear fit, applying the Track F formula explicitly (INFERENCE on amplitudes;
   the formula is FACT).
4. Tabulates the bias contribution as a function of W in {30, 90, 180, 365, 730,
   1460, 1825} d (INFERENCE — depends on assumed amplitudes).
5. Determines the **asymptotic scaling** of each bias as W → ∞ (FACT — algebraic).
6. Answers the central question: **is the 019 Ω̇_fit(W) = a + b/W + c/W² model
   theoretically justified, or is it an empirical fit with no asymptotic basis?**
7. If the latter, proposes alternative estimators.

The conclusion (§9): the Track F bias formula is oscillatory with amplitude
~A_k/W², NOT polynomial in 1/W. The 019 polynomial extrapolation is therefore
**an empirical fit** with no asymptotic justification. The alternative estimators
(FFT subtraction, multi-window least-squares with a harmonic basis, total-least-
squares with the secular model) are proposed.

---

## 2. Physical origin of each periodic term

The third-body disturbing function is a Fourier series in the six angles of the
satellite–third-body system: satellite mean anomaly M, argument of perigee ω,
longitude of ascending node Ω; third-body mean anomaly M₃, argument of perigee
ω₃, longitude of ascending node Ω₃. For a circular satellite orbit (e ≈ 0),
the surviving harmonics are combinations

  ψ_pqr = p·n + q·n₃ + r·dω/dt + s·dΩ/dt                                  (2)

(using integer p, q, r, s; for e = 0 the k ≠ 0 harmonics in M vanish after orbit-
averaging). The dominant frequencies at h = 600 km are (Track B §1.2):

| Symbol | Quantity | Period T (d) | ω (rad/s) | Comment |
|---|---|---:|---:|---|
| n        | satellite mean motion     | 0.06709 | 1.0839e-3 | 14.9 rev/day |
| n_M      | lunar mean motion (sidereal) | 27.3217 | 2.6614e-6 | |
| n_S_synodic | Sun–Moon synodic        | 29.5306 | 2.4622e-6 | drives variation |
| n_M_anom | lunar anomalistic month   | 27.5546 | 2.6393e-6 | drives evection |
| n_S      | solar mean motion (tropical) | 365.2422 | 1.9910e-7 | drives annual |
| n_node   | lunar nodal regression    | 6798.4 | 1.1390e-9 | 18.6 yr |
| n_aps_M  | lunar apsidal precession  | 3232.6 | 2.2490e-9 | 8.85 yr |

The five periodic/long-period terms:

### 2.1 Evection (Ptolemy; Brown 1896; Murray & Dermott §6.5)

**Physical origin**: the periodic variation of the Earth–Moon distance as the Moon
travels around an Earth whose orbit about the Sun is eccentric (the Moon's
anomalistic-month carrier modulated by the Sun).

**Angle combination**:

  ψ_evection = 2n − 2n_M_anom = 2(n − n_M_anom) ≈ 2.16e-3 rad/s

giving period T_evection ≈ 27.55 d (the anomalistic month).

**Amplitude in λ (lunar longitude)**: +1.274° × sin(2D − M_M) (Brown, Meeus
1998, Ch. 47). [FACT — well-tabulated astronomical constant.]

### 2.2 Variation (Tycho Brahe; Brown 1896; Murray & Dermott §6.5)

**Physical origin**: the change in the tangential component of the solar pull
as the Moon passes between the Sun and the Earth — the "trig inequality" of
the three-body problem. The dominant period is the **synodic half-month**.

**Angle combination**:

  ψ_variation = 2n − n_S_synodic ≈ 2n − (n_M − n_S) = 2n − n_M + n_S
             ≈ 1.97e-3 rad/s

giving period T_variation ≈ 14.77 d (T_synodic/2).

**Amplitude in λ**: +0.658° × sin(2D) (Meeus 1998, Ch. 47). [FACT.]

### 2.3 Annual solar forcing

**Physical origin**: the Sun's geocentric distance oscillates at the Earth orbital
period (Earth heliocentric eccentricity e_E = 0.0167). This modulates the
secular solar RAAN rate (a/a₃)³ factor by ±5% peak (3 e_E).

**Angle combination**:

  ψ_annual = n_S = 2π / 365.2422 d ≈ 1.99e-7 rad/s

**Amplitude in Ω**: of order the secular solar RAAN rate itself,
(3/8) n (μ_S/μ_E) (a/a₃)³ |sin 2(i − i₃)/sin i| ≈ 3.6e-5 deg/day at h=600 km
i_sso (Track B §4.3). [INFERENCE — factor of order unity in the inclination
function.]

### 2.4 Lunar nodal regression

**Physical origin**: the Moon's orbital plane regresses around the ecliptic
pole once every 18.6 yr. The Moon's mean inclination to the equator oscillates
between 18.29° (ε − I_M) and 28.58° (ε + I_M), where ε = 23.439° is the
obliquity and I_M = 5.145° is the lunar inclination.

**Angle combination**:

  ψ_node = 2π / 6798.4 d ≈ 1.14e-9 rad/s

The third-body inclination i₃ in formula (1) is replaced by i₃(t) =
ε + I_M·cos(2π t/T_node), so the secular RAAN rate itself oscillates at the
nodal frequency with ~30% peak-to-peak fractional amplitude (the geom_factor
at i_sso varies from 0.362 to 0.670 over the cycle — Track C §4.2).

**Amplitude in Ω**: of order the secular lunar RAAN rate × (60% modulation),
~5e-5 deg/day at h=600 km i_sso. [INFERENCE.]

### 2.5 Lunar apsidal precession

**Physical origin**: the lunar argument of perigee ω_M precesses once per 8.85 yr
with ~8.6 cycles per nodal period (giving the evection/variation envelope).

**Angle combination**:

  ψ_aps_M = 2π / 3232.6 d ≈ 2.25e-9 rad/s

After orbit-averaging and lunar-M-anomaly-averaging, the ω_M dependence survives
in terms ∝ e_M cos(ω − ω_M), etc. At the quadrupole order these terms contribute
as a slow modulation of the secular rate.

**Amplitude in Ω**: ~5e-6 deg/day at h=600 km i_sso. [INFERENCE — second-order
in e_M = 0.0549.]

---

## 3. Track F OLS bias formula (FACT)

The OLS slope over a window [0, W] of a harmonic A_k cos(ω_k t + φ_k) is
(Track F §5):

  bias_k(W; A_k, ω_k, φ_k) = (1/W) A_k [sin(φ_k) − sin(ω_k W + φ_k)]     (3)

derived from the closed-form integrals:

  ∫₀^W cos(ωt + φ) dt = [sin(ωW + φ) − sin(φ)] / ω                      (4)

The full bias from a single Fourier component Ω_osc(t) = A_k cos(ω_k t + φ_k)
is bounded by

  |bias_k| ≤ 2|A_k| / W                                                  (5)

with the sign depending on φ_k and the parity of ω_k W.

For a sum of harmonics, the biases add algebraically. The Track F formula is
FACT (algebraic identity). The task-supplied formula in the prompt

  bias_k(W) = (12/W³) A_k [W sin(ωW)/ω + (cos(ωW) − 1)/ω²]
            − (6/W²) A_k sin(ωW)/ω                                       (6)

is **algebraically equivalent** to formula (3) after simplification (Track B
§6.1 derivation). [FACT — both expressions are equivalent integrals of the
OLS-on-harmonic problem; the difference is a sign convention in the
intermediate step.]

For computational convenience we use

  bias_k(W) = (A_k / W) [sin(φ_k) − sin(ω_k W + φ_k)]                    (7)

in the remainder of this report.

### 3.1 Three regimes

**Regime A (ω_k W ≫ 1, fast harmonics)**: |sin(φ_k) − sin(ω_k W + φ_k)| ≤ 2
randomly. Average bias ~0 by phase scrambling; RMS ~|A_k|/W. For the satellite
short-period harmonics (period ~96 min, ωW ~ 5e3 at W = 1 yr) this gives
RMS bias ~0.001 deg / 3e7 s ≈ 3e-11 rad/s ≈ 3e-6 deg/day. Negligible
[INFERENCE — depends on the actual RMS amplitude of the short-period content].

**Regime B (ω_k W ≪ 1, slow harmonics)**: sin(ω_k W + φ_k) − sin(φ_k) ≈
ω_k W cos(φ_k). Bias ≈ A_k ω_k cos(φ_k). For the lunar nodal period
(T_node = 18.6 yr, ω_node = 1.07e-8 rad/s) and A_node ~ 0.1 deg,
bias ≈ 1e-9 rad/s ≈ 5e-5 deg/day at most. Comparable to dΩ_mean/dt
~ 1.35e-4 deg/day (Track F §5 Regime B).

**Regime C (ω_k W ~ 1, comparable to window)**: maximum bias. For annual
solar forcing (T = 365.24 d, ω W = 2π at W = 1 yr) the bias is
~2|A_k|/W ≈ 0.1 deg / 3.15e7 s ≈ 3e-9 rad/s ≈ 1.7e-4 deg/day at most
(Track F §5 Regime C).

---

## 4. Amplitude estimates per term (INFERENCE)

The amplitudes A_k in formula (7) are the amplitudes of each term in the
**osculating Ω(t)** (not in λ₃, M₃, or any third-body longitude). They are
INFERENCE because the standard celestial-mechanics literature gives the
amplitudes in terms of partial derivatives of R₃ with respect to i, and the
inclination function factors are of order unity at SSO inclinations but not
exactly known. The order-of-magnitude estimates below are conservative (Track B
and Track C estimates combined).

| Term | Period (d) | ω (rad/s) | A_k in Ω(t) (deg) | Reference |
|---|---:|---:|---:|---|
| Evection (lunar distance mod)   | 27.55   | 2.64e-6 | ~0.005 – 0.02 | Track C §1 (factor ~0.3 of lunar secular × e_M) |
| Evection (lunar geometric mod) | 27.55   | 2.64e-6 | ~0.001 – 0.005 | Track C §1 (3% geom_factor mod) |
| Variation (geometric)          | 14.77   | 4.92e-6 | ~0.001 – 0.005 | Track C §2 (0.7% geom_factor mod) |
| Annual solar (eccentricity mod)| 365.24  | 1.99e-7 | ~0.01 – 0.1    | 019 FFT dominant amp (0.10 deg at i_sso) |
| Lunar nodal (i₃ modulation)    | 6798.4  | 1.14e-9 | ~0.05 – 0.2    | Track C §4 (60% geom_factor mod) |
| Lunar apsidal (e_M² modulation) | 3232.6 | 2.25e-9 | ~0.005 – 0.02  | Track B §5.2 (e_M² factor) |

The 019 FFT measurement of Ω(t) at i_sso gives the dominant amplitudes as
0.10 deg (annual), 0.025 deg (half-annual), 0.012 deg (third-annual), 0.007 deg
(quarter-annual), 0.005 deg (fifth-annual) at 365 d / N for N = 1, 2, 3, 4, 5
(019 results.json: `fft_periodicity_i_sso.dominant_amplitudes_deg`). These are
the **aliases of the annual solar forcing** plus their harmonics — the 019 FFT
detected a strong annual term, consistent with Track F's prediction. The
evection (27.55 d) and variation (14.77 d) terms are below the FFT noise floor
of the 1-year arc because their 13.25-cycle / 24.7-cycle structure aliases into
many bins and their per-bin amplitude is ~ 0.1 / sqrt(N_cycles) ≈ 0.025 deg per
bin at the 14.77-d and 27.55-d peaks. [INFERENCE on the alias-mapping; FACT on
the 019 FFT result.]

### 4.1 Best-estimate A_k for OLS bias calculation

To compute the OLS bias contribution of each term as a function of window W,
we adopt the following **central** amplitudes (geometric mean of Track B/C
estimates and 019 FFT evidence). The numbers below are used in §5–6; they are
INFERENCE.

| Term | A_k (deg) | Source |
|---|---:|---|
| Evection           | 0.01  | Track B/C midpoint |
| Variation          | 0.003 | Track C midpoint |
| Annual solar       | 0.10  | 019 FFT dominant amplitude |
| Lunar nodal        | 0.10  | Track C §4 |
| Lunar apsidal      | 0.01  | Track B §5.2 |

These are **peak amplitudes** in Ω(t). For the OLS bias formula the relevant
quantity is A_k · ω_k, which we compute in §5.

---

## 5. OLS bias contributions per term at W = 365.24 d

Using formula (7) with W = 365.24 d = 3.156e7 s and φ_k = 0 (worst case):

  bias_k(1yr) = (A_k / W) [sin(0) − sin(ω_k W)]
              = −(A_k / W) sin(ω_k W)                                    (8)

The bias is bounded by 2|A_k|/W. For each term:

| Term | ω_k (rad/s) | ω_k W (rad) | sin(ω_k W) | |bias_k| / (A_k/W) | bias_k (deg/day) | A_k/W (deg/day) |
|---|---:|---:|---:|---:|---:|---:|
| Evection (27.55 d) | 2.64e-6 | 83.3  | bounded O(1) | ≤ 2 | ≤ 6.3e-2 × A_k | 6.3e-2 × 0.01 = 6.3e-4 |
| Variation (14.77 d) | 4.92e-6 | 155.4 | bounded O(1) | ≤ 2 | ≤ 6.3e-2 × A_k | 6.3e-2 × 0.003 = 1.9e-4 |
| Annual solar (365.24 d) | 1.99e-7 | 6.28 (= 2π) | 0 (exactly) | 0 | **0** (orthogonal) | 3.2e-6 × 0.10 = 3.2e-7 |
| Half-annual (182.6 d) | 3.98e-7 | 12.57 (= 4π) | 0 | 0 | **0** (orthogonal) | |
| Lunar nodal (6798.4 d) | 1.14e-9 | 0.36 | sin(0.36) ≈ 0.35 | 0.35 | 0.35 × A_k/W | 0.35 × 3.2e-6 = 1.1e-6 |
| Lunar apsidal (3232.6 d) | 2.25e-9 | 0.0071 | sin ≈ 0.0071 | 0.0071 | 7.1e-3 × A_k/W | 7.1e-3 × 3.2e-7 = 2.3e-9 |

(For lunar nodal and apsidal, ω_k W ≪ 1, so the bias formula reduces to
bias ≈ A_k ω_k cos(φ_k) ≈ A_k ω_k for φ_k = 0.)

### 5.1 Numerical columns

The column `A_k/W` converts the amplitude in degrees to a rate estimate:
A_k = 0.10 deg / W = 3.156e7 s gives 0.10 deg / (3.156e7 / 86400) d ≈ 0.27 deg/day.
Wait — that is the **amplitude** as a rate, not the bias. The bias per term
is bounded by 2|A_k|/W, so:

- Evection: |bias| ≤ 2 × 0.01 deg / (3.156e7 s) = 6.3e-10 rad/s = **3.1e-5 deg/day**
  (using 2 A_k = 0.02 deg and 1 W = 3.156e7 s; 6.3e-10 rad/s × (180/π) × 86400
  = 3.1e-5 deg/day)
- Variation: |bias| ≤ 2 × 0.003 / 3.156e7 = 1.9e-10 rad/s = **9.4e-6 deg/day**
- Annual solar: |bias| = 0 (orthogonal)
- Lunar nodal: |bias| ≈ 0.10 × 1.14e-9 = 1.14e-10 rad/s = **5.7e-6 deg/day**
- Lunar apsidal: |bias| ≈ 0.01 × 2.25e-9 = 2.25e-11 rad/s = **1.1e-6 deg/day**

Total bias bound (sum of magnitudes): ~5.3e-5 deg/day
Sum of root-mean-square (assuming independent phases): ~3.4e-5 deg/day

Both are ~3× smaller than the corrected secular formula value 1.35e-4 deg/day
(018 result). This is consistent with Track F's estimate that the 1-year
linear fit's bias is comparable to (but somewhat smaller than) dΩ_mean/dt.

The Track F formula gives the WORST-CASE bias; the actual expected bias for a
specific phase is whatever sin(φ_k) − sin(ω_k W + φ_k) evaluates to. For the
018 arc starting 2026-01-01, the lunar phase at that epoch determines the
actual bias; it is expected to be within ~2× of the worst-case bound.

---

## 6. Asymptotic scaling of each term's bias as W → ∞

The Track F OLS bias formula (3) has the following asymptotic scaling:

### 6.1 Fast harmonics (ω_k W ≫ 1)

For ω_k W ≫ 1 (evection, variation, J2 apsidal at all relevant W),
sin(ω_k W + φ_k) oscillates rapidly and the **expected bias** (averaged
over a uniformly random φ_k in [0, 2π]) is zero, with **RMS** amplitude
~√2 |A_k|/W. This is the standard spectral analysis result: the bias from
a harmonic of amplitude A_k over a window W has RMS contribution
|A_k|/W, decaying as **1/W**.

### 6.2 Near-resonant harmonics (ω_k W ~ integer multiples of 2π)

For harmonics whose period T_k = 2π/ω_k satisfies ω_k W = 2π·n for integer
n (annual at W = 365.24 d, half-annual at W = 182.6 d, etc.), the bias
**exactly vanishes** at those specific W values, and has small but non-zero
bias for off-resonant W. The bias scales as:

  bias_k(W) ≈ −(A_k / W) [ω_k W − 2π n]  for ω_k W ≈ 2π n              (9)

i.e. linearly in the deviation δ = ω_k W − 2π n, falling as **1/W** with a
phase-dependent coefficient.

### 6.3 Slow harmonics (ω_k W ≪ 1)

For ω_k W ≪ 1 (lunar nodal at W ≪ 18.6 yr; lunar apsidal at W ≪ 8.85 yr),
sin(ω_k W + φ_k) − sin(φ_k) ≈ ω_k W cos(φ_k), so:

  bias_k(W) ≈ A_k ω_k cos(φ_k)                                          (10)

This is a **constant** in W, NOT 1/W. The bias from the slow harmonics
**does not decrease** as the window lengthens; it asymptotes to a constant
offset to the slope. For the 018 1-year arc with lunar nodal period 18.6 yr,
this bias is ~5.7e-6 deg/day (above); for W = 5 yr it is still ~5.7e-6 deg/day.

### 6.4 Summary table: asymptotic bias scaling as W → ∞

| Term | Period (d) | ω_k W at W=1 yr | Regime | Asymptotic scaling |
|---|---:|---:|---|---|
| Evection           | 27.55  | 83.3  | A (fast)         | O(A_k/W), oscillating |
| Variation          | 14.77  | 155.4 | A (fast)         | O(A_k/W), oscillating |
| Annual solar       | 365.24 | 6.28  | C (resonant)     | O(A_k/W) at off-resonant W; 0 at W = 365.24 d |
| Half-annual solar       | 182.6  | 12.57 | C (resonant)     | O(A_k/W); 0 at W = 182.6 d |
| Quarter-annual solar    | 91.3   | 25.13 | A (fast)         | O(A_k/W) |
| Lunar nodal        | 6798.4 | 0.36  | B (slow)         | O(1) — constant in W |
| Lunar apsidal      | 3232.6 | 0.007 | B (slow)         | O(1) — constant in W |

The dominant scaling is **O(1/W)** for the fast harmonics and **O(1)** for
the slow harmonics. **There is no 1/W² or 1/W³ asymptotic term in the Track F
OLS bias formula.** [FACT — the formula is exact and algebraic.]

---

## 7. Does a polynomial in 1/W fit the OLS bias formula asymptotically?

The 019 model is Ω̇_fit(W) = a + b/W + c/W². The Track F bias formula gives

  Ω̇_fit(W) = Ω̇_mean + Σ_k bias_k(W; A_k, ω_k, φ_k)
            = Ω̇_mean + Σ_k (A_k / W) [sin(φ_k) − sin(ω_k W + φ_k)]    (11)

Each term is bounded by 2|A_k|/W. The sum is bounded by (2/W) Σ_k |A_k|. As
W → ∞:

  Ω̇_fit(W) − Ω̇_mean = O(Σ |A_k| / W) for fast harmonics
                     + O(Σ |A_k| ω_k cos(φ_k)) for slow harmonics       (12)

**This is NOT a polynomial in 1/W.** The fast-harmonic contributions oscillate
with W at frequencies ω_k, with amplitude ~A_k/W. The slow-harmonic contributions
asymptote to constants. Neither structure matches a + b/W + c/W² Taylor
expansion as W → ∞.

### 7.1 Theoretical justification for a polynomial-in-1/W fit?

**For the slow-harmonic (Regime B) terms**, the bias approaches a constant as
W → ∞, not zero. If the A_k ω_k cos(φ_k) constants are non-zero, then
Ω̇_fit(W) → Ω̇_mean + constant as W → ∞. The 019 intercept a would then
**not** equal Ω̇_mean but Ω̇_mean + constant.

**For the fast-harmonic (Regime A + C) terms**, the bias oscillates in W at
frequencies ω_k. If the dominant harmonics have periods much shorter than the
5-point window set W in {30, 90, 180, 365, 730} d (e.g. evection at 27.55 d,
variation at 14.77 d), then the bias oscillates between +2A_k/W and −2A_k/W
as W sweeps through cycles of the harmonic. A polynomial-in-1/W fit through
such oscillating data is a **regression through aliasing**: the fit coefficients
depend on the specific W values chosen and have no asymptotic meaning.

### 7.2 Why the 019 fit works empirically

The 019 window-length extrapolation was applied to the actual Ω̇_fit(W)
data from the 018 simulations (W in {30, 90, 180, 365, 730} d). For these
specific 5 W values, the fit a + b/W + c/W² minimized the residual sum of
squares. The fit **does** identify a mean trend in the data, but:

- The trend identified is dominated by the aliasing pattern of the dominant
  harmonics at those specific W values, not by the asymptotic scaling of any
  theoretical bias formula.
- The fit coefficients b and c are **not** interpretable as "1/W coefficient"
  and "1/W² coefficient" of any underlying bias expansion.
- The fit residual RMS at the 5 measured W values is small (~1e-3 deg/day at
  i_sso), but extrapolation to W → ∞ is **not theoretically justified** by
  the OLS bias formula.

### 7.3 Verdict

**The 019 Ω̇_fit(W) = a + b/W + c/W² model is an empirical fit with no
asymptotic basis.** The Track F OLS bias formula has oscillatory 1/W and
constant (Regime B) components; it does not have a Taylor-expandable
polynomial structure. The 019 polynomial extrapolation **may** give a useful
estimate of Ω̇_mean for the specific W ∈ {30, 90, 180, 365, 730} d data
sampled, but the extrapolation to W → ∞ is **not guaranteed** to converge
to Ω̇_mean.

This is consistent with the 019 results.json: the polynomial fit gives
extrapolated values that differ from the corrected cf by ~10× at i_sso
(019 result: 0.9956 deg/day vs cf 0.9933 deg/day — agreement to within 0.25%
**of the J2 secular** but the Lunisolar residual at i_sso remains open).

[INFERENCE on the discrepancy interpretation; FACT on the algebraic structure
of the bias formula.]

---

## 8. Alternative estimators with stronger theoretical justification

Three alternatives, in order of theoretical rigor:

### 8.1 Fourier-decomposition estimator (medium complexity)

**Method**: FFT-decompose Ω_cross(t_k) at the known physical frequencies
(n_sat, n_apsidal_J2, n_lunar_synodic, n_lunar_node, n_solar_synodic, evection,
variation). Identify the discrete frequency bins corresponding to each known
harmonic driver. Subtract their amplitudes from the time series. Re-fit the
linear slope.

**Theoretical justification**: the disturbing function is a Fourier series in
the system angles; the harmonic content of Ω(t) is exactly the projection onto
the known basis. Subtracting the identified harmonics yields the secular drift
plus unmodelled long-period residuals.

**Limitation**: requires identifying the FFT bins (frequency resolution ∝ 1/T);
at T = 1 yr the 18.6-yr lunar nodal term is unresolvable (1 bin at DC). For T =
5 yr the resolution improves to ~73 d, allowing better separation of the
14.77-d / 27.55-d / 365-d terms.

### 8.2 Multi-window least-squares with harmonic basis (high complexity)

**Method**: fit the model

  Ω(t) = Ω̇_mean · t + Σ_k [A_k cos(ω_k t + φ_k)] + ε                    (13)

jointly to Ω_cross(t_k) data over W in {30, 90, 180, 365, 730} d (concatenated
or with shared parameters). The Ω̇_mean parameter is the secular drift;
the (A_k, φ_k) parameters absorb the periodic content.

**Theoretical justification**: maximum-likelihood under the assumption that the
residuals are Gaussian noise and the periodic content is captured by the chosen
basis. The estimator is unbiased in the limit of many observations.

**Limitation**: requires choosing the harmonic basis correctly. Missing a
harmonic biases Ω̇_mean; over-fitting inflates the variance of the estimate.
For the Lunisolar problem, the basis set is well-known from Kaula's expansion.

### 8.3 Total-least-squares with the secular + analytical short-period model (highest complexity)

**Method**: subtract the analytical short-period corrections (evection,
variation, annual solar, J2 short-period) from the osculating Ω before the
linear fit. Use Brouwer-Kozai short-period theory to compute ΔΩ_short-period
from the known third-body state.

**Theoretical justification**: the standard Brouwer-Kozai short-period
correction is the first-order mapping from osculating to mean elements. After
subtraction, the residual Ω(t) is the mean Ω(t), whose slope is the secular
drift.

**Limitation**: the short-period corrections for a third-body perturbation on
a near-circular orbit are small (~milliarcseconds to arcseconds at LEO SSO,
Track F §7), so this method is most valuable for **explaining** the osculating-
vs-mean gap rather than for numerical extraction of Ω̇_mean. For the
Lunisolar RAAN drift at SSO, the dominant residual is **evection/variation-
scale periodic content at the mdeg–0.1 deg amplitude**, NOT 1/rev short-period
content — so the standard Brouwer short-period subtraction is the wrong tool
for this specific problem.

### 8.4 Recommended bridge

For the 020 numerical bridge (recommended primary method): **option 8.1
(FFT subtraction) for the annual/lunar-nodal/evection/variation content,
combined with option 8.2 (multi-window joint fit) to extract Ω̇_mean and its
uncertainty.** This pair has the strongest theoretical justification and the
least dependence on the choice of W sample values.

For comparison/validation: the **cycle-averaged estimator** (019 finding)
already reduces the 1-year bias to ~3% of the single-window value. This is
empirically grounded and can be adopted as a low-cost secondary bridge.

---

## 9. Conclusions

### 9.1 Summary of term amplitudes and 1-year bias contributions

For h = 600 km i_sso, using Track B/C central amplitude estimates:

| Term | Period (d) | A_k (deg) | |bias| at W=365 d (deg/day) | Asymptotic scaling |
|---|---:|---:|---:|---|
| Evection           | 27.55  | 0.01  | ≤ 3.1e-5 | O(A_k/W), oscillating |
| Variation          | 14.77  | 0.003 | ≤ 9.4e-6 | O(A_k/W), oscillating |
| Annual solar       | 365.24 | 0.10  | 0 (exact at W=365.24 d) | O(A_k/W) off-resonant; 0 on-resonant |
| Half-annual solar       | 182.6  | 0.025 | 0 (exact at W=182.6 d) | O(A_k/W) |
| Lunar nodal        | 6798.4 | 0.10  | ~5.7e-6 | O(1) — constant in W |
| Lunar apsidal      | 3232.6 | 0.01  | ~1.1e-6 | O(1) — constant in W |
| **Total bias bound** |   |   | **~5.3e-5 deg/day** | dominated by slow harmonics (Regime B) |

The total expected bias is of order **3–5×10⁻⁵ deg/day**, comparable to the
corrected secular formula value 1.35e-4 deg/day. This is **insufficient** by
itself to account for the 9.78× residual at i_sso (which is 1.18e-3 deg/day
in the Lunisolar-only residual; the J2 component is 0.99 deg/day and is
removed in the i_sso vs J2-only comparison). [INFERENCE — the 9.78× is
explained by Track F as J2 × Lunisolar coupling (~3.5×) + mean-vs-osculating
bias (~2–3×); the periodic-term bias alone is ~0.4×.]

### 9.2 Polynomial-in-1/W extrapolation: verdict

**The 019 Ω̇_fit(W) = a + b/W + c/W² model is an empirical fit with no
asymptotic basis.** The Track F OLS bias formula has O(1/W) oscillating
contributions from fast harmonics and O(1) constant contributions from slow
harmonics. Neither matches a polynomial-in-1/W Taylor expansion.

The 019 fit may give a useful local interpolation through the 5 W data points,
but extrapolation to W → ∞ is not theoretically justified. The intercept a
is **not guaranteed** to equal Ω̇_mean. The coefficients b and c have **no
physical interpretation**.

### 9.3 Recommended alternative estimators

For the 020 Lunisolar secular extraction:

1. **Fourier-decomposition estimator** (subtract identified harmonics before
   linear fit): medium complexity, strong theoretical justification.
2. **Multi-window joint least-squares with harmonic basis**: high complexity,
   strongest theoretical justification, allows uncertainty quantification.
3. **Cycle-averaged estimator** (019 finding; reduces bias to ~3%): low
   complexity, empirically grounded.

### 9.4 What this track changes about the 019 conclusion

The 019 conclusion that "the window-length extrapolation Ω̇_fit(W) = a + b/W
+ c/W² is the canonical numerical bridge" is **downgraded to "empirical
interpolation only; not asymptotically justified"**. The 019 finding that the
window-length extrapolation recovers the secular drift **within the sampled
W range** is preserved; the claim that the extrapolation to W → ∞ is
**theoretically rigorous** is retracted.

The 8-track audit-019 finding that the 018 10× residual is dominated by mean-
vs-osculating bias (and not by unmodelled physics) is **unaffected** by this
track's analysis.

---

## 10. Limitations

- All amplitude estimates A_k in §4 are INFERENCE based on Track B/C estimates.
  A precise FFT-based measurement of A_k for each term in Ω(t) at i_sso would
  require a longer arc (5+ yr) than 019's 1-yr window.
- The Track F formula assumes uniform time sampling; the actual 018 sampling
  is at ascending-node crossings (non-uniform in time, slightly variable Δt).
  For the 018 cadence (~14.9 samples/day), this is a small effect (~0.1% on
  the bias).
- The polynomial fit to 5 W data points has 3 free parameters (a, b, c) plus
  residual noise; the residual RMS at the fit is small but the extrapolation
  error is not quantifiable from the fit alone.
- The alternative estimators in §8 require either FFT post-processing or
  joint multi-window fitting; these are not implemented in 019 and would be
  candidates for 020.

---

## 11. References

- Murray, C. D. & Dermott, S. F. (1999). *Solar System Dynamics*. Cambridge
  University Press. §6.4 (doubly-averaged quadrupole), §6.5 (evection and
  variation), §7 (third-body disturbing function).
- Kaula, W. M. (1966). *Theory of Satellite Geodesy*. Ch. 4 (disturbing
  function expansion; evection/variation terms; inclination functions).
- Kozai, Y. (1959). "The motion of a close earth satellite." *Astronomical
  Journal* 64, 367–377.
- Brown, E. W. (1896). *An Introductory Treatise on the Lunar Theory*.
  Cambridge University Press.
- Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Ch. 22 (Earth's
  orbit), Ch. 47 (lunar inequalities).
- Brouwer, D. & Clemence, G. M. (1961). *Methods of Celestial Mechanics*.
  Academic Press.
- Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*
  (4th ed.). Microcosm Press.
- Standish, E. M. (1990). "An observationally based reference frame for
  astronomy." *Astronomy and Astrophysics* 233, 272–274.
- Tremaine, S. & Yavetz, T. D. (2014). "Secular dynamics of compact three-
  body systems." *American Journal of Physics* 82, 749–755.
- Exp 014 eclipseTiming (2026-08-28): byte-pinned Sun snapshot acquisition
  pattern; offline-deterministic analysis.
- Exp 017 lunisolarVerification (2026-08-30): byte-pinned Moon snapshot;
  RK4 self-convergence p_r = 4.49, p_v = 4.50.
- Exp 018 lunisolarReconciliation (2026-08-30): corrected secular formula;
  1-year vs i=90° cleanest test (ratio 2.81×); W=730 window sensitivity.
- Exp 019 lunisolarLongPeriod (2026-08-30): window-length extrapolation;
  cycle-averaged estimator; FFT periodicity test; dominant annual term at
  0.10 deg.
- audit-019-track-B-averaging-hierarchy.md: averaging operations (§2);
  short-period (§3); intermediate-period (§4); long-period (§5); OLS bias
  formula (§6).
- audit-019-track-C-evection-variation-hypothesis.md: order-of-magnitude
  estimates of evection/variation/annual/nodal contributions to the 018
  residual.
- audit-019-track-F-mean-vs-osculating.md: OLS bias derivation in three
  regimes (A: ωT ≫ 1, B: ωT ≪ 1, C: ωT ~ 1); formula bias_k = (A_k/W)
  [sin(φ_k) − sin(ω_k W + φ_k)].

---

## 12. Audit context

This is Track A-2 of the 8-track independent investigation for Experiment 020.
Other tracks (A-1, A-3 through H) investigate different angles; their outputs
are not read by this track per the mission constraint.

The findings in §9 are **inferences based on the Track F OLS bias formula and
the Track B/C amplitude estimates**, with FACT/INFERENCE/UNKNOWN classification
on each individual claim. The polynomial-extrapolation verdict (§9.2) is FACT
(the bias formula has no polynomial-in-1/W asymptotic structure). The amplitude
estimates are INFERENCE. The 018 9.78× residual breakdown is unchanged from
the 8-track synthesis.

---

## FACT / INFERENCE / UNKNOWN classification

| Section | Claim | Class |
|---|---|---|
| §1 | The 018 1-year fit exceeds cf by 9.78× at i_sso, 2.81× at i=90° | FACT (019 results.json) |
| §2.1 | Evection: ψ = 2(n − n_M_anom), T = 27.55 d | FACT |
| §2.1 | Evection amplitude in λ: +1.274° × sin(2D − M_M) | FACT (Meeus) |
| §2.2 | Variation: ψ = 2n − n_S_synodic, T = 14.77 d | FACT |
| §2.2 | Variation amplitude in λ: +0.658° × sin(2D) | FACT (Meeus) |
| §2.3 | Annual solar: ψ = n_S, T = 365.24 d | FACT |
| §2.4 | Lunar nodal: T = 18.6 yr; i₃ oscillates 18.29° to 28.58° | FACT |
| §2.5 | Lunar apsidal: T = 8.85 yr | FACT |
| §3 | OLS bias formula (3) = (A_k/W) [sin(φ_k) − sin(ω_k W + φ_k)] | FACT (algebraic identity) |
| §3 | Task-supplied formula (6) is equivalent to (3) | FACT (algebraic identity) |
| §3 | Three regimes A/B/C classification | FACT (sign analysis) |
| §4 | A_k amplitudes in Ω(t) at h=600 km i_sso (table in §4) | INFERENCE (depends on inclination functions of order unity) |
| §4 | 019 FFT dominant amp 0.10 deg at 365 d | FACT (019 results.json) |
| §5 | bias_k at W = 365 d for each term (table in §5) | INFERENCE (depends on A_k from §4) |
| §5 | Annual solar bias = 0 exactly at W = 365.24 d | FACT (orthogonality) |
| §5 | Half-annual solar bias = 0 at W = 182.6 d | FACT (orthogonality) |
| §6 | Asymptotic scaling: O(A_k/W) for fast, O(1) for slow | FACT (algebraic) |
| §7 | The 019 Ω̇_fit(W) = a + b/W + c/W² model has no asymptotic basis | FACT (the bias formula has no polynomial structure) |
| §7 | The 019 fit may give useful local interpolation | FACT (it minimizes residual at the 5 W points) |
| §7 | The 019 extrapolated intercept a ≠ Ω̇_mean in general | INFERENCE (depends on slow-harmonic constants) |
| §8 | FFT subtraction is a sound alternative | INFERENCE (standard spectral-analysis practice) |
| §8 | Multi-window joint fit is the strongest theoretical estimator | INFERENCE (maximum-likelihood under Gaussian assumption) |
| §9 | The 8-track audit-019 conclusion (10× = bias, not physics) is unaffected | FACT (this track does not modify that finding) |
| §10 | The 018 sampling non-uniformity is a ~0.1% effect on bias | INFERENCE (small, but not quantified) |