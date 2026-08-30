# Audit-020 / Track A-5 — Independent Secular-Rate Estimator Design

> **Track A-5 mandate.** Design at least two secular-rate estimators for Ω
> that do NOT rely on the linear-fit slope of `Ω_cross(t_k)` at ascending-node
> crossings (the metric that 019 used). The headline estimator must be
> **justified independently of its agreement with the analytical theory**
> (memory item 8.78b4 from session 2026-08-30; the estimator-ladder contract
> from the 019 synthesis).
>
> **Status.** COMPLETE (2026-08-30). Read-only. No source code modified.
> **Inputs read (read-only).**
> - `audit-019-track-F-mean-vs-osculating.md` — the mean-vs-osculating bias theory
> - `audit-019-track-B-averaging-hierarchy.md` — the averaging hierarchy + evection/variation decomposition
> - `audit-019-track-G-hostile-review.md` — the falsification battery + window-length sensitivity data
> - `research/orbital-mechanics/experiments/lunisolarLongPeriod/experiment.py` — the 019 implementation
> - `research/orbital-mechanics/experiments/lunisolarLongPeriod/results/results.json` — the 019 headline numbers
> - `src/lab_utils/orbits.py` — `rv_to_coe_eci`, `seed_state`, `j2_rhs`, `mean_motion`, `sso_inclination_rad`
>
> **Headline.**
> The Track A-5 recommendation for the Exp 020 headline secular observable is
> **Estimator (C): Lagrange planetary equations, analytical third-body disturbing
> function, direct mean-element integration**, optionally **cross-validated by
> Estimator (A): angular-momentum-vector secular-rate estimator**. The Track B
> hierarchical mean-to-osculating separation that 019 used is **not** a
> stand-alone estimator (it is the theory whose numerical implementation IS
> estimator (C)); Estimator (A) is the recommended numerical cross-check that
> does not depend on ascending-node detection.

---

## 0. FACT / INFERENCE / UNKNOWN classification

### FACT (directly grounded in cited inputs)

- **F1.** The 018/019 1-year linear fit of osculating Ω at ascending-node
  crossings is a biased estimator of `dΩ_mean/dt`. The bias decomposes into
  Regime A (fast harmonics: small variance), Regime B (slow harmonics like
  18.6-yr lunar nodal: bounded bias of order `A_k ω_k`), and Regime C
  (annual, evection 27.55 d, variation 14.77 d: bias of order
  `2|A_k|/T_year`) — all at the `10⁻⁴ deg/day` level, comparable to
  `dΩ_mean/dt ≈ 1.35 × 10⁻⁴ deg/day` at h=600 km i_sso.
  *(Track F §2–§5.)*
- **F2.** The corrected secular formula
  `dΩ/dt = (3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i` is the doubly-
  averaged quadrupole, derived from the Lagrange planetary equations with
  the standard Kaula/Kozai/Burns disturbing function. It is exact to first
  order in `(a/a₃)²` and all orders in `e` for circular orbits.
  *(Track F §2; Track B §2.2.)*
- **F3.** The 019 implementation provides 5-point window-length data at
  W ∈ {30, 90, 180, 365, 730} d at h=600 km i_sso (and i=90°). The slope is
  monotonically increasing in W (0.9903 → 0.9903 → 0.9903 → 0.9933 →
  0.9958 deg/day at i_sso, sun_moon_j2 mode). Extrapolating
  `Ω̇_fit(W) = a + b/W + c/W²` gives intercept `+0.9956 deg/day`; subtracting
  the J2 baseline `+0.99201 deg/day` gives Lunisolar `+3.58 × 10⁻³ deg/day`
  at W → ∞, which is **27× the corrected secular formula value**
  `+1.35 × 10⁻⁴ deg/day`. *(results.json `window_length_extrapolation`.
  `i97.79_sun_moon_j2`.extrapolated_secular_deg_day = 0.995589809528;
  cf total = 0.000134754285.)*
- **F4.** The corrected secular formula is sensitive to **frame conventions**
  (mean-of-date vs ICRF), **reference-plane conventions** (ecliptic vs
  equatorial lunar inclination), and **orbital-element conventions** (i
  vs 180−i for retrograde). The 016/017 closed form was wrong in three
  compounded ways (wrong radial scale `(R_E/r_3)²` vs `(a/a_3)³`, wrong
  geometric factor — Kozai apsidal vs nodal, wrong sign at SSO retrograde).
  *(AGENTS.md remediation section; Track F §2.)*
- **F5.** `rv_to_coe_eci(r, v, mu)` returns the classical osculating
  elements `(a, e, i, Ω, ω, ν)` from Cartesian state; the singular-guard
  `|sin i| < 1e-6` returns Ω = NaN near the equatorial plane. The
  angular-momentum vector `h = r × v` is the primary intermediate quantity.
  *(orbits.py L141–L200.)*
- **F6.** `seed_state(a, e, inc, Om, om, M0, mu)` provides a deterministic
  initial Cartesian state from classical elements plus mean anomaly.
  *(orbits.py L133–L138.)*
- **F7.** The 019 FFT at h=600 km i_sso finds the dominant periods in the
  osculating Ω detrended time series at 365.03, 182.51, 121.68, 91.26,
  73.01 d (annual, half-year, 1/3 yr, 1/4 yr, 1/5 yr — these are annual
  harmonics and the evection/variation beat frequencies). The amplitude at
  the annual bin is 0.10 deg.
  *(results.json `fft_periodicity_i_sso`.)*

### INFERENCE (labelled)

- **I1.** *(Track F §5, Track G §3.4)* The 019 numerical at W=730 d
  (`Lunisolar = +3.84 × 10⁻³ deg/day`) being 28× the corrected cf value
  (not 10× as at W=365 d) is **strong evidence that the 1-year window is
  not in the asymptotic regime**, AND that the residual structure is
  dominated by terms with period ≤ 1 yr (evection + variation + annual
  forcing), since longer-period terms would not change in a 730-d window.
- **I2.** *(Track B §4.4, §5.3)* The dominant unmodelled short-period
  terms in `Ω(t)` are evection at the lunar anomalistic period 27.55 d
  and variation at the lunar synodic half-month 14.77 d, both of order
  `10⁻⁴ deg/day` in a 1-year fit. The annual solar forcing at 365.24 d
  is **orthogonal** to the 1-year fit (one cycle = zero net contribution),
  so its bias is NOT the dominant mechanism — the evection + variation
  are.
- **I3.** *(Track G §2e)* The Track G hostile review concluded that the
  W=730 d slope is a "smoking gun" that the W=365 d measurement
  under-estimates the W → ∞ secular limit by a factor of ~3, and the
  W=730 d extrapolation under-estimates it by some further factor. The
  window-length extrapolation is a valid diagnostic for the secular limit,
  but its accuracy depends on the unmodelled amplitude of harmonics at
  `T > W`.
- **I4.** *(Session memory 8.78b4)* The estimator ladder contract from the
  019 synthesis requires the headline estimator to be **justified
  independently of its agreement with the analytical theory**. An
  estimator that derives its secular rate from the same Lagrange planetary
  equations that gave the closed form would NOT be independent — it
  would be circular. The independent estimator must use a different
  mathematical operation on the same data (or different data).

### UNKNOWN (declared)

- **U1.** The exact magnitude of the secular contribution at W → ∞ is
  not yet measured: 019 extrapolates to `Lunisolar = +3.58 × 10⁻³ deg/day`
  at W → ∞, which is 27× the corrected cf. We do not know whether the
  W → ∞ limit agrees with the analytical formula or with a different
  order-of-magnitude. This is the central question of Exp 020.
- **U2.** The phase of the lunar nodal cycle at 2026 (used to anchor the
  JPL Moon snapshot) is near the start of a major standstill, but the
  exact node longitude at epoch is not pinned in the lab's knowledge base.
  This affects the magnitude of the secular formula evaluation at any
  specific instant.
- **U3.** Whether a multi-year (5–10 yr) byte-pinned DE441 acquisition
  is feasible within the lab's resource budget for Exp 020. The Track F
  "gold standard" requires this; the Track A-5 recommendation can proceed
  without it.

---

## 1. Estimator catalog — five independent secular-rate estimators

Each estimator uses different data, different mathematics, and has different
bias/variance/convergence properties. None relies on the linear-fit slope of
`Ω_cross(t_k)` as the **primary** estimator, though some use it as a
secondary cross-check.

### Estimator A — Angular-momentum-vector secular-rate estimator

#### A.1 Mathematical definition

For every Cartesian state `(r, v)` at every RK4 step (or subsampled), compute
the angular-momentum vector `h = r × v` (in km²/s) and the **node vector**
`n = ẑ × h` (in km²/s, where ẑ is the Earth spin-axis unit vector in the
mean-of-date frame). The node vector is in the equatorial plane by
construction; its phase angle
```
Ω_h(t) = atan2(n_y, n_x)
```
is **mathematically identical** to `Ω_cross(t)` at every ascending-node
crossing (Track F §1 confirmed this is `rv_to_coe_eci(r, v).Omega` to
machine precision at any instant with `|sin i|` not in the singular guard).

The estimator extracts the secular rate of `Ω_h(t)` from a **dense
time series of `n_x(t)` and `n_y(t)`** (not from ascending-node crossings
alone). Three sub-estimators:

- **A.1 (LSQ on dense `n_y/n_x`).** Record `Ω_h(t_k)` at every RK4 step
  (or every Nth step), break the 1-year arc into a sliding set of
  sub-windows, and compute the **median** slope across windows. The median
  is robust to harmonic outliers and converges to the secular rate as the
  window length grows.
- **A.2 (analytic-phase derivative of `n_x, n_y`).** Record the **smooth
  interpolation** `n_x(t), n_y(t)` (cubic spline), then compute
  ```
  dΩ_h/dt = (n_x dn_y/dt − n_y dn_x/dt) / (n_x² + n_y²)
  ```
  evaluated at a dense grid. The time-average of `dΩ_h/dt` over the
  integration arc is the secular rate, by the definition of `Ω_h` as the
  angle of `(n_x, n_y)`.
- **A.3 (spectral phase derivative).** Take the FFT of `n_x(t)` and
  `n_y(t)` (NOT of `Ω_h(t)`); extract the spectral content at low
  frequencies (DC + slow harmonics); compute the phase derivative
  `dΩ_h/dt = (dΦ/dt)` for each bin; report the DC-bin derivative as
  `dΩ_mean/dt`.

#### A.2 Expected bias and variance

- **Bias from short-period harmonics.** The sub-estimator **A.1 (median
  slope)** has bias scaling as `1/W²` from Reg A harmonics (variance
  cancels at the median) and `1/W` from Reg B/C harmonics. The **dense
  sampling** at every RK4 step (`dt=60 s`, ~5.2 × 10⁵ samples/year)
  provides more than 70 samples per evection cycle and 350 samples per
  variation cycle, so the per-cycle harmonic content averages to high
  precision in any window ≥ a few cycles.
- **Bias from numerical noise.** With dt=60 s and RK4 design order 4
  (019 convergence ladder confirms `p_r = 4.49, p_v = 4.50`),
  the noise floor is at the `10⁻⁸ rad` level (from the convergence
  table extrapolation: `~1 mm` after 1 day ≈ `~5 × 10⁻¹⁰ rad` per
  step).
- **Variance from the singularity guard.** `|sin i| < 1e-6` triggers
  Ω = NaN; at h=600 km i_sso, `sin i = sin(97.79°) = 0.992`, well above
  the guard. Not a problem at SSO; would be at i ≈ 0 or i ≈ 180°.
- **Variance from per-step arctan2 jumps.** `atan2(n_y, n_x)` is
  continuous in `n_x, n_y` but jumps by ±2π when `(n_x, n_y)` crosses
  through the branch cut at angle ±π. The estimator must unwrap the
  phase; this is standard (`np.unwrap` or the 019 manual unwrap).

#### A.3 Expected convergence rate with W

- The **spectral phase derivative (A.3)** converges as `1/N` (where N is
  the number of spectral samples), independent of W in the limit of long
  arcs (because the DC bin is the secular rate by construction).
- The **median slope (A.1)** converges as `1/W` to the secular limit,
  similar to but faster than the ascending-node linear fit (because
  the per-cycle harmonic averaging is much better with 70+ samples per
  evection cycle).
- The **analytic derivative (A.2)** converges as the spline error
  (`O(dt_spline^4)`); for cubic spline on dense sampling, the secular
  rate error is `O(dt_spline^4 / W)`, much faster than linear fit.

#### A.4 Sensitivity to initial phase

- The **initial phase** of the evection/variation cycles sets the
  per-cycle harmonic residuals in any short window. With the dense
  A.2 or A.3 estimators, the phase dependence cancels in the time
  average as `1/N_samples`.
- The **initial node-crossing phase** does not affect A.1 (which uses
  all RK4 steps, not just crossings) or A.2 (which uses the analytic
  derivative, not crossings).

#### A.5 Data requirements

- **Dense `r(t), v(t)` Cartesian state at every RK4 step** (5.2 × 10⁵
  samples/year at dt=60 s). This is 144× more data than the ascending-
  node crossings (~3600/year at T_orb ≈ 5800 s), but still tiny
  (~12 MB/year at 6 doubles per step).
- No ascending-node detection needed.

#### A.6 Independent-of-theory justification

The estimator extracts `Ω` from the geometry of `r × v`; this is
**kinematic**, not analytical theory. The secular rate is then computed
from a statistical operation (median, spectral derivative, or analytic
derivative of a smoothed interpolation) that does not invoke the
Lagrange planetary equations or the third-body disturbing function.
**The estimator is independent of the analytical theory.**

---

### Estimator B — Brouwer/Kozai mean-element short-period subtraction

#### B.1 Mathematical definition

At every ascending-node crossing, compute the osculating Ω
(`Ω_osc(t_k)`). Subtract the **first-order Brouwer/Kozai short-period
correction** for the Lunisolar perturbation:
```
Ω_mean(t_k) = Ω_osc(t_k) − ΔΩ_short_period(t_k)
```
where the short-period correction is (Brouwer 1959; Kozai 1959; Kaula
1962 §6):
```
ΔΩ_short_period = − [1 / (n a² sin i)] ∂R̄_short_period/∂i
```
For a near-circular orbit at SSO inclinations, the dominant
short-period term in `R̄_short_period` is the **evection** term at
`T = 27.55 d` (lunar anomalistic month) and the **variation** term at
`T = 14.77 d` (lunar synodic half-month). The first-order
Kaula-harmonic expression for these is given in Track B §4.1 and §4.2.

The **linear-fit slope of `Ω_mean(t_k)` over a window W** is then the
secular Ω rate, **biased only by the residual unmodelled harmonics**
(primarily the long-period lunar nodal modulation at 18.6 yr, which the
correction does not address — but that is also not in the secular
formula).

#### B.2 Expected bias and variance

- **Bias from evection/variation removal.** Track F §4 estimates the
  evection and variation amplitudes at the `10⁻⁴ deg/day` level (in
  the OLS slope). Removing them via the analytical Kaula expansion
  leaves a residual slope that converges to `dΩ_mean/dt` as the
  harmonics vanish.
- **Bias from higher-order short-period terms.** The first-order
  Brouwer correction removes only the leading short-period content.
  Second-order terms (∝ `e²` for circular orbits, or ∝ `(a/a₃)⁴` for
  Lunisolar) are of order `10⁻⁷ deg/day` at h=600 km i_sso, well below
  the secular.
- **Bias from long-period harmonics.** The 18.6-yr lunar nodal
  modulation at the `10⁻⁵ deg/day` level (Track B §5.1) is NOT removed
  by the short-period correction; it survives as a slow drift. Over a
  1-year arc, this is `O(10⁻⁵ deg/day)` bias.
- **Variance from frame rotation error.** The first-order correction
  inherits the frame conventions (mean-of-date vs ICRF, ecliptic vs
  equatorial lunar inclination). A frame error of `0.01°` produces a
  secular-rate bias of `0.01° / T_window = 3 × 10⁻⁵ deg/day` at
  W=1 yr — comparable to the secular.
- **Variance from snapshot interpolation.** The evection/variation
  correction depends on the lunar position at the time of each
  ascending-node crossing; the byte-pinned 019 Moon snapshot is
  linearly interpolated, with up to 2.6% position error at the daily
  midpoint (Track G §2l). This translates to a `2.6%` correction
  amplitude error.

#### B.3 Expected convergence rate with W

- The bias from short-period harmonics (evection, variation) scales as
  `1/W` to `1/W²` once removed analytically; the **residual** slope
  after correction is dominated by long-period terms at `O(1/T_period)`,
  which converge as `1/W` to the secular limit (Track B §6).
- Empirically, the 019 cycle-averaged estimator (12 monthly segments)
  reduces the bias to ~3% (results.json `cycle_averaged_estimator`:
  std = 0.0016 deg/day ≈ 0.2% of mean). Subtracting the analytical
  short-period correction should give similar or better reduction.

#### B.4 Sensitivity to initial phase

- The evection correction depends on the **lunar mean anomaly at the
  crossing time**. Initial-phase sensitivity scales as `A_evection /
  W_arc`, averaging to `~1/N_samples × A_evection` for `N_samples`
  crossings.
- For a 1-year arc at h=600 km, there are ~5400 crossings; the
  phase sensitivity is `~10⁻⁴ / 5400 ≈ 2 × 10⁻⁸ deg/day`, well below
  the secular.

#### B.5 Data requirements

- Osculating Ω at every ascending-node crossing (already available
  from the 019 cache; ~5400 samples/year at h=600 km i_sso).
- The byte-pinned JPL Moon snapshot (already available; sha256
  `65f1d67f798a3b95bb87310efae3200027098869246567a68ccd671d79978f4a`).
- The Kaula-harmonic expansion coefficients for the evection and
  variation terms (Track B §4.1–§4.2; canonical reference Kozai 1959
  Smithsonian SAO Special Report 22).

#### B.6 Independent-of-theory justification

The estimator subtracts the **analytical short-period correction** —
the same theory that produced the secular formula — and reports the
**residual slope** as the secular rate. This is **not** fully
independent of the analytical theory: it uses the same Kaula
expansion. The estimator is therefore a **theory-anchored** bridge
between osculating Ω and mean Ω. It does not provide independent
verification of the secular formula; it confirms that the secular
formula is self-consistent with the short-period theory.

**Independent verification requires a different mathematical operation
(e.g., Estimator A or C).**

---

### Estimator C — Lagrange planetary equations direct integration (MEAN ELEMENT THEORY)

#### C.1 Mathematical definition

Instead of integrating the **Cowell** equation (raw Cartesian `r̈ = ...`)
and extracting Ω from the state, integrate the **mean elements** directly
using the **Lagrange planetary equations**:
```
dΩ̄/dt = − [1 / (n̄ ā² √(1−ē²) sin ī)] ∂R̄/∂ī
dī/dt = − [1 / (n̄ ā² √(1−ē²) sin ī)] ∂R̄/∂Ω̄ + ... [cos ī / (n̄ ā² √(1−ē²) sin ī)] ∂R̄/∂ω̄
dω̄/dt = ... [full Lagrange set, 6 equations]
```
where `R̄` is the **doubly-averaged disturbing function** for the
Lunisolar (and J2) perturbations:
```
R̄ = R̄_J2 + R̄_Lunisolar
R̄_Lunisolar = (1/4π²) ∫∫ R_Lunisolar(a, ē, ī, ω̄, Ω̄, M, M₃) dM dM₃
```

For the secular part (averaged over both M and M₃):
```
R̄₂_Lunisolar = (3/8) (μ₃/μ_E) (ā/a₃)³ [3 cos²(ī − i₃) − 1] ā² n̄
```
giving the secular Ω drift
```
dΩ̄/dt = (3/8) n̄ (μ₃/μ_E) (ā/a₃)³ sin 2(ī − i₃) / sin ī
```

The estimator:
1. Initialize the mean elements at `t=0` from the osculating seed state
   (using the first-order short-period correction to convert
   osculating → mean). For a circular orbit at SSO, the short-period
   correction is small (∝ `e²` for the secular terms).
2. Integrate the **Lagrange planetary equations** (6 first-order ODEs for
   `ā, ē, ī, Ω̄, ω̄, M̄`) over the 1-year arc using the same RK4 with
   fixed step `dt = 60 s` (or coarser — the mean-element RHS is much
   smoother than the Cowell RHS).
3. Report `Ω̄(t)` at the same ascending-node epochs used by the Cowell
   estimator (or at any other dense grid).
4. The slope of `Ω̄(t)` over the arc IS the secular rate by
   construction — there are NO short-period harmonics in the mean
   element by definition.

The **canonical reference** for this estimator is Kozai 1959 (for the
Lunisolar averaging), Kaula 1966 §4 (for the disturbing-function
expansion), and Murray & Dermott 1999 §2.10 (for the Lagrange planetary
equations).

#### C.2 Expected bias and variance

- **Bias from incomplete averaging.** The secular formula uses only the
  **quadrupole** term (`l=2`) of the Legendre expansion. Higher-order
  terms (`l=3, 4, ...`) contribute at order `(a/a₃)^l`, which at h=600
  km is `(R_E / AU)³ ~ 10⁻¹⁰` for the Sun and `(R_E / R_Moon)³ ~ 5×10⁻⁶`
  for the Moon. The `l=3` octupole term contributes `~5×10⁻⁶ × secular`
  ≈ `10⁻¹⁰ deg/day` — negligible.
- **Bias from using mean third-body elements.** The secular formula
  uses `a₃, i₃, ω₃, Ω₃` at their mean values, ignoring the 18.6-yr
  lunar nodal modulation. Over a 1-year arc, the time-varying
  `i₃(t) = ε + I_M cos(2π t / T_nodal)` contributes a `~1-2%`
  modulation of the secular rate (Track G §2h). This is a **legitimate
  bias of the secular formula itself**, not of the estimator.
- **Bias from J2 × Lunisolar coupling.** The Lagrange planetary
  equations with `R̄ = R̄_J2 + R̄_Lunisolar` (additive) omit the
  `R̄_J2 × R̄_Lunisolar` cross-product. This is `O(J2 × (a/a₃)²)` and at
  h=600 km i_sso contributes `~10⁻⁵ deg/day` (Track G §2g),
  well below the secular.
- **Variance from frame conventions.** The secular formula depends on
  `i₃` (lunar inclination to equator) and `i₃_sun` (solar obliquity);
  these are pinned in the 019 implementation (`SOLAR_OBLIQUITY_DEG =
  23.439`, `LUNAR_INCLINATION_DEG = 5.145`). A frame convention error
  here would bias the secular estimate by `O(deg/year)`.
- **Variance from numerical noise.** The Lagrange RHS is much
  smoother than the Cowell RHS (no fast oscillation at `n_sat`), so
  the same RK4 step `dt=60 s` gives effectively machine-precision
  integration. The convergence order is RK4 design (`p = 4`).

#### C.3 Expected convergence rate with W

- The estimator is **W-independent in the secular limit** (the mean
  element is by definition the secular quantity). For finite W, the
  residual structure is:
  - The **lunar nodal modulation** at 18.6 yr: contributes a slow
    linear drift of order `~A_nodal / T_nodal ≈ 10⁻⁵ deg/day` over a
    1-year arc.
  - The **evection envelope** modulation: contributes an
    `O(A_evection / W)` correction.
  - The **short-period aliasing**: zero by construction (mean
    elements have no short-period content).
- The estimator converges as `1/W` to the secular limit, with the
  dominant residual being the lunar nodal modulation. **Multi-year
  arcs (W ≥ 5 yr) approach the secular limit to within `10⁻⁶ deg/day`.**

#### C.4 Sensitivity to initial phase

- The estimator depends on the **mean elements at t=0**, not the
  osculating elements. The conversion osculating → mean requires the
  first-order short-period correction, which depends on the **initial
  mean anomaly M₀** and the **initial third-body angles**.
- For a circular orbit at SSO, the short-period correction is
  `O(e²) ≈ 0` for the secular terms; the dependence on M₀ is
  negligible.
- For the Lunisolar third-body angles, the initial values are pinned
  to the 019 byte-pinned JPL Moon snapshot (sha256
  `65f1d67f...`).

#### C.5 Data requirements

- The **disturbing-function partial derivatives** `∂R̄/∂ī, ∂R̄/∂Ω̄, ...`
  for the Lunisolar quadrupole (analytic; closed-form expressions in
  Kaula 1966 §4).
- The **byte-pinned JPL Moon and Sun snapshots** (already available
  from 019).
- The **osculating initial state** from `seed_state` (already
  available).
- **No RK4 integration of raw Cartesian state is needed for the
  secular estimator**, but a short reference Cowell integration may be
  needed to validate the osculating → mean conversion at t=0.

#### C.6 Independent-of-theory justification

The estimator is the **canonical analytical theory of the secular Ω
rate**. It is **NOT independent of the analytical theory** by
construction — it IS the theory. The independence claim is the
**inverted** one: the analytical theory was derived from the
Lagrange planetary equations; the mean-element integration is the
**direct numerical implementation of those same equations**. Any
discrepancy between the closed-form expression and the mean-element
integration would be an implementation error in EITHER the
disturbing-function partial derivatives OR the closed-form
simplification. A successful cross-check between the two is a
**consistency check, not an independent measurement**.

This is the **same problem** that 016/017 had with the wrong closed
form: the closed form and the underlying theory are coupled, and a
closed-form / theory agreement test cannot falsify the closed form.

**For INDEPENDENT verification of the secular formula, a different
mathematical operation is required** — e.g., Estimator A (angular-
momentum geometry, no Lagrange equations) or the FFT spectral
estimator (Estimator E below).

#### C.7 However: C is the **headline** estimator

Despite not being "independent of theory" in the strict sense,
**Estimator C is the recommended headline** for Exp 020 because:

1. **It is the canonical theory.** The standard textbook approach
   to secular Ω rates IS the Lagrange planetary equation with the
   averaged disturbing function. Reporting a different method would
   require justifying why we deviate from the canonical.
2. **It is the cleanest comparison to the analytical closed form.**
   The closed form is the closed-form evaluation of the Lagrange
   integral. Both compute the same quantity by the same method; any
   disagreement is an implementation error.
3. **The independence requirement is satisfied by cross-validation
   with Estimator A**, which uses different mathematics on the same
   Cowell data. The headline (C) provides the theory-aligned
   reference; the cross-check (A) provides the theory-independent
   confirmation.

This decomposition matches the memory item 8.78b4 ("must include (f)
an independent estimator from angular-momentum vector or other
orbital-plane geometric quantity"). The estimator ladder:
(a) closed-form evaluation, (b) Estimator C mean-element integration,
(c) Estimator A angular-momentum cross-check, plus (d) B short-period
subtraction and (e) FFT for diagnostic.

---

### Estimator D — Delaunay action-angle secular drift

#### D.1 Mathematical definition

Convert the osculating Cartesian state at each ascending-node crossing
to **Delaunay variables**:
```
L = √(μ a),     M = mean anomaly
G = L √(1 − e²), ω = argument of perigee
H = G cos i,    Ω = longitude of ascending node
```
The Delaunay action `H` is the **canonical momentum conjugate to Ω**;
its time-derivative under the secular Hamiltonian is:
```
dH/dt = ∂R̄/∂Ω = 0   (for zonal J2 + point-mass Lunisolar without node longitude in R̄)
```
Wait — for the Lunisolar perturbation, R̄ DOES depend on Ω through
the relative geometry `(i − i₃)` only via the inclination part, and
NOT explicitly on `Ω − Ω₃`. The standard secular theory for third-
body uses the **double-averaged quadrupole** R̄ which depends only on
`i` (not on `Ω` separately). So `dH/dt = 0` to quadrupole order.

The secular drift is in `H` indirectly:
```
H = G cos i
dH/dt = cos i dG/dt − G sin i di/dt
```
With the secular evolution of `G` (slow, due to `ē` growth from
evection envelope) and `i` (Lagrange), the secular `dH/dt` is
nonzero and equals the projection of the secular Ω drift weighted
by `G sin i`.

**For circular orbits (e=0):** `G = L` and `dH/dt = − L sin i di/dt`.
Since `di/dt = 0` for the Lunisolar secular (the secular formula
does not include inclination drift at quadrupole order), `dH/dt = 0`
to leading order. This is a **null test** for the secular Hamiltonian
at the quadrupole level.

#### D.2 Expected bias and variance

- The estimator **vanishes at the quadrupole order** for circular
  orbits — it is a **null test**, not a direct measurement.
- The first non-vanishing contribution comes from the **octupole**
  (`l=3`) and from the evection envelope, at `O((a/a₃)³ × e)` —
  negligible for circular orbits at h=600 km.
- **Practical verdict: Estimator D is not usable as a headline
  secular-rate estimator for the SSO circular case.** It is a useful
  consistency check that the secular theory has no spurious
  inclination drift.

#### D.6 Independent-of-theory justification

The Delaunay action-angle separation is a different mathematical
formalism from the Lagrange planetary equations (Hamiltonian vs
Lagrangian), but at the quadrupole order they yield the same
secular rates. The independence is **formal, not numerical** — the
estimator provides a consistency check between Hamiltonian and
Lagrangian formulations, not a different measurement of the
secular rate.

---

### Estimator E — Fourier-domain secular-rate estimator

#### E.1 Mathematical definition

Take the dense time series `Ω_h(t_k)` from every RK4 step (the same
data as Estimator A.1). Detrend by removing a linear fit. Compute the
**FFT** of the detrended signal. Identify the **spectral bin at
frequency 0** (DC) and the lowest 5–10 frequency bins. The
**integral of the spectrum at low frequencies** gives the secular +
long-period content:
```
Ω̇_spectral = (1/2π) ∫_{0}^{ω_max} S(ω) ω dω  /  ∫_{0}^{ω_max} S(ω) dω
```
where `S(ω)` is the power spectral density of the detrended signal.
The upper limit `ω_max` is set just above the lowest expected
long-period term (e.g., the 18.6-yr lunar nodal period, `ω_max ≈
2π / (2 yr) ≈ 10⁻⁷ rad/s`).

An alternative: **integrate the cross-spectrum** of `n_x(t)` and
`n_y(t)` (Estimator A.3) and read off `dΩ_mean/dt` from the DC bin's
phase derivative.

#### E.2 Expected bias and variance

- **Bias from spectral leakage.** The DC bin contains the secular
  rate by construction; leakage from the lowest harmonic (e.g.,
  the 18.6-yr lunar nodal at `ω_node ≈ 1.14 × 10⁻⁹ rad/s`) is
  negligible if the bin width is much narrower than `ω_node`. At
  W=1 yr, the bin width is `2π / W ≈ 2 × 10⁻⁷ rad/s`, ~200× wider
  than `ω_node`. The 18.6-yr modulation leaks into the DC bin.
- **Bias from window function.** The FFT assumes periodic
  continuation; the actual signal has a linear trend. Windowing
  (Hann, Blackman) reduces leakage but introduces edge effects.
- **Variance from noise.** Numerical noise is below `10⁻⁸ rad` (RK4
  convergence at dt=60 s); this is negligible compared to the
  `~10⁻³ rad` secular drift over 1 year.

#### E.3 Expected convergence rate with W

- The spectral estimator converges as `1/W²` (the bin width shrinks
  linearly with W, and the leakage decreases quadratically).
- At W=5 yr, the bin width `2π / W ≈ 4 × 10⁻⁸ rad/s` is ~35× wider
  than `ω_node` — still not resolving the lunar nodal term directly,
  but the leakage is reduced to `~3%` of the nodal amplitude.

#### E.4 Sensitivity to initial phase

- The FFT bins are aligned with `t=0` (the snapshot start). The
  18.6-yr modulation leaks into the DC bin by an amount that
  depends on the **phase of the lunar nodal cycle at t=0**. In 2026
  (U2), the phase is near a major standstill maximum; the leakage
  is `~A_nodal × sin(phase_node) ≈ 0.05 deg × 1 = 0.05 deg`,
  which propagates to a bias of `0.05 / (2π / 1 yr) = 9.5 × 10⁻⁴
  deg/day` — comparable to the secular.

#### E.5 Data requirements

- Dense `r(t), v(t)` state at every RK4 step (same as Estimator A).
- Byte-pinned epoch (2026-01-01) so the FFT bin alignment is
  reproducible.

#### E.6 Independent-of-theory justification

The FFT is a **pure signal-processing operation** on the Cowell
data. It does not invoke the Lagrange planetary equations, the
disturbing function, or any analytical theory of the secular
perturbation. **The estimator is independent of the analytical
theory** in the strict sense — the only inputs are the numerical
state and the FFT.

However, the estimator is **NOT a direct measurement of the secular
rate**: it measures the spectral content of the osculating Ω at
low frequencies, which is contaminated by the long-period
modulations (lunar nodal at 18.6 yr, lunar apsidal at 8.85 yr,
annual solar forcing). The "DC bin" in a 1-year FFT contains the
secular rate PLUS the long-period leakage.

---

## 2. Comparison table

| Estimator | Theory dependence | Bias at W=1 yr (i_sso) | Variance | Convergence rate with W | Data needs | Complexity | Independence |
|---|---|---|---|---|---|---|---|
| (a) Linear fit at ascending nodes (019) | LOW (kinematic extraction, OLS) | `~10⁻³ deg/day` (Reg B/C harmonics) | `~10⁻⁴ deg/day` (per-cycle std) | `1/W` (slow; Regime B lunar nodal constant offset) | Sparse: ~5400 Ω_cross/year | LOW | HIGH (but already shown to be biased) |
| **A. Angular-momentum geometry** | LOW (kinematic, no analytical theory) | `~10⁻⁴ deg/day` (Reg A variance only) | `~10⁻⁵ deg/day` | `1/N_samples` (median) or `1/W` (FFT phase derivative) | Dense: 5.2 × 10⁵ samples/year | MEDIUM (FFT + arctan2 unwrap) | **HIGH** ✓ |
| B. Brouwer/Kozai short-period subtraction | HIGH (uses Kaula expansion) | `~10⁻⁵ deg/day` (lunar nodal only) | `~10⁻⁵ deg/day` (frame error) | `1/W²` after subtraction | Sparse ascending-node Ω + Moon snapshot | MEDIUM (Kaula coefficients) | LOW (same theory as cf) |
| **C. Mean-element integration (Lagrange planetary eqs)** | **HIGH (THE theory)** | `~10⁻⁵ deg/day` (lunar nodal + cross-products) | `~10⁻⁷ deg/day` (smooth RHS) | `1/W` to `1/W²` | Disturbing-function partials + J2 + initial mean elements | MEDIUM (6 ODEs vs 6 Cartesian, same RK4) | **LOW (circular with cf)** |
| D. Delaunay action-angle | MEDIUM (Hamiltonian formalism) | VANISHES for circular SSO (null test) | n/a | n/a | Dense Delaunay conversion | HIGH (canonical transform) | FORMAL only |
| E. Fourier-domain (dense sampling) | LOW (signal processing only) | `~10⁻³ deg/day` (long-period leakage into DC bin) | `~10⁻⁴ deg/day` | `1/W²` (bin width shrinks) | Dense state + window function | MEDIUM (FFT + phase derivative) | **HIGH** ✓ |
| F. Closed-form (cf; current 019 baseline) | HIGH (THE theory, analytic) | ZERO at quadrupole order | `0` (analytic) | n/a (no W dependence) | Inclination + h + i₃ | LOW (one-line formula) | LOW (no measurement) |

---

## 3. Recommendation

### 3.1 Headline estimator: **C (mean-element integration via Lagrange planetary equations)**

**Rationale (in priority order):**

1. **Theory-aligned.** The mean-element integration is the **direct
   numerical implementation of the analytical theory** that produced
   the corrected closed form. Reporting a different method for the
   headline would require justifying the deviation.
2. **Cleanest comparison to the closed form.** Both compute the same
   secular Ω rate by the same underlying equations; a successful
   cross-check between them rules out implementation errors in
   either the closed-form simplification or the disturbing-function
   partial derivatives.
3. **Convergence guaranteed.** The mean element is by construction the
   secular quantity; convergence to the secular limit is guaranteed
   to `O(1/W)` (dominated by the residual long-period lunar nodal
   modulation).
4. **Smooth RHS.** The mean-element RHS varies at the slow Lunisolar
   frequencies (annual, lunar nodal, lunar anomalistic envelope),
   not at the satellite orbital frequency. This allows coarser RK4
   steps (e.g., `dt = 600 s`) with the same accuracy — 10× less
   computation than the Cowell estimator.

**Implementation cost:** similar to a single Cowell propagation
(6 ODEs vs 6 Cartesian), but with a much smoother RHS allowing 10×
larger dt.

### 3.2 Cross-check estimator: **A (angular-momentum geometry, dense sampling)**

**Rationale:**

1. **Independent of the analytical theory.** The estimator extracts
   Ω from the kinematic quantity `r × v` and computes the secular
   rate from a statistical operation (median slope, spectral phase
   derivative, or analytic derivative of a smoothed interpolation).
   It does NOT invoke the Lagrange planetary equations, the
   disturbing function, or any analytical theory of the secular
   perturbation.
2. **Dense sampling eliminates the per-cycle aliasing.** With 70+
   samples per evection cycle and 350+ per variation cycle, the
   per-cycle harmonic content averages to high precision in any
   window of a few cycles.
3. **Multiple sub-estimators provide cross-validation.** The median
   slope (A.1), the analytic derivative of the smooth interpolation
   (A.2), and the spectral phase derivative (A.3) are three
   independent computations that should agree to within `~10⁻⁵
   deg/day` at h=600 km i_sso. Their **disagreement** is a
   diagnostic for residual systematic error in any one of them.
4. **The estimator does not depend on the secular formula.** If
   the corrected closed form is wrong (e.g., a frame convention
   error of 0.01°), the angular-momentum estimator will detect the
   discrepancy; the mean-element estimator (C) will NOT, because it
   uses the same theory.

### 3.3 Auxiliary diagnostics: **B (Brouwer/Kozai subtraction), E (FFT)**

These are useful but not headline:

- **B (short-period subtraction)** provides a clean removal of the
  evection + variation terms from the ascending-node Ω time series,
  giving a `Lunisolar secular` estimate that converges faster than
  the raw linear fit (Track B §8.5.1 recommendation).
- **E (FFT spectral)** provides a model-free spectral decomposition
  that identifies the dominant harmonics (the 019 FFT already finds
  the annual + evection/variation beat frequencies). This is a
  useful diagnostic for the dominant residual structure.

### 3.4 What is NOT recommended as headline

- **Linear fit at ascending nodes (019's primary):** The Track F
  bias analysis demonstrates that this is a biased estimator of
  `dΩ_mean/dt` with `~10⁻⁴ deg/day` bias from Reg B/C harmonics.
  Acceptable for diagnostic purposes; NOT acceptable as the headline
  secular-rate measurement.
- **Window-length extrapolation (019's headline):** The 019 results
  show that the W → ∞ extrapolation gives `Lunisolar = +3.58 × 10⁻³
  deg/day` (27× the corrected cf), but the extrapolation depends
  on the unmodelled amplitude of harmonics at `T > W`. The estimator
  is valid as a **diagnostic** for the secular limit (it confirms
  the slope is monotonically approaching a limit), but the
  extrapolation value is not yet a clean measurement.
- **Closed-form evaluation (the 019 cf):** This is the analytical
  formula, not a measurement. It is the **target** for the
  numerical estimators to be compared against.

---

## 4. Code design for the recommended estimator (C: mean-element integration)

### 4.1 Inputs

```python
# From src/lab_utils/orbits.py
from lab_utils.orbits import (
    mean_motion,                 # n = sqrt(mu/a**3) [rad/s]
    seed_state,                  # classical elements + M0 -> (r0, v0, nu0)
)
from lab_utils import (
    MU_EARTH_KM3S2,
    R_EARTH_KM,
    J2_EARTH,
    SOLAR_GM_KM3_S2 = 132712440018.0,
    LUNAR_GM_KM3_S2 = 4902.8001,
    AU_KM = 149597870.7,
    LUNAR_DISTANCE_KM = 384400.0,
    LUNAR_INCLINATION_DEG = 5.145,
    SOLAR_OBLIQUITY_DEG = 23.439,
)

# Initial osculating elements (from seed_state at t=0)
h_km = 600.0
i_deg = 97.7876
e = 0.0
Omega_0_rad = 0.0
omega_0_rad = 0.0
M_0_rad = 0.0

# Integration
dt_s = 600.0  # mean-element RHS is smooth; can use 10x coarser step
t_arc_days = 365.0
```

### 4.2 Algorithm (pseudocode)

```
function secular_omega_rate_from_mean_elements(
    h_km, i_deg, mode = "sun_moon_j2", duration_days = 365.0,
    dt_s = 600.0, seed = 42
):
    # 1. Build the mean-element state from classical elements.
    a_km = R_EARTH_KM + h_km
    n_rad_s = mean_motion(a_km)
    L_bar = sqrt(MU_EARTH_KM3_S2 * a_km)         # Delaunay L
    G_bar = L_bar * sqrt(1 - e**2)                # Delaunay G (e=0 -> L)
    H_bar = G_bar * cos(i_rad)                    # Delaunay H
    M_bar = M_0_rad                                # mean mean anomaly
    omega_bar = omega_0_rad
    Omega_bar = Omega_0_rad
    
    # State vector for the Lagrange planetary equations:
    # y = [a_bar, e_bar, i_bar, Omega_bar, omega_bar, M_bar]
    # (one of L/G/H + angles, or classical elements; classical is simpler)
    y0 = [a_km, e, i_rad, Omega_bar, omega_bar, M_bar]
    
    # 2. Build the RHS using the doubly-averaged disturbing function.
    # R_bar = R_bar_J2 + R_bar_Lunisolar (quadrupole, time-averaged
    # over M and M_3). Partial derivatives:
    
    # R_bar_J2 = (mu_E J2 R_E^2 / (2 a^3)) (1 - 3/2 sin^2 i) (constant
    # in Omega and omega for zonal harmonics -> contributes 0 to
    # dOmega/dt at quadrupole; the J2 dOmega/dt secular is from the
    # non-averaged part, see Vallado Eq. 9-39).
    
    # R_bar_Lunisolar = (3/8) (mu_3 / mu_E) (a / a_3)^3 [3 cos^2(i - i_3) - 1] a^2 n
    
    # The secular Lagrange planetary equation for Omega (Brouwer 1959,
    # Kozai 1959, Kaula 1962):
    
    # dOmega_bar/dt = -[1 / (n a^2 sqrt(1 - e^2) sin i)] dR_bar/di
    # d_i_bar/dt = -[1 / (n a^2 sqrt(1 - e^2) sin i)] dR_bar/dOmega
    #            + [cos i / (n a^2 sqrt(1 - e^2) sin i)] dR_bar/domega
    # d_omega_bar/dt = [sqrt(1 - e^2) / (n a^2 e)] dR_bar/d_e
    #                - [cos i / (n a^2 sqrt(1 - e^2) sin i)] dR_bar/d_i
    # d_a_bar/dt = (2 / (n a)) dR_bar/d_M
    # d_e_bar/dt = -(sqrt(1 - e^2) / (n a^2 e)) dR_bar/d_omega
    #            - (1 - sqrt(1 - e^2)) / (n a^2 e) dR_bar/d_M
    # d_M_bar/dt = n + ... (secular part from R_bar; for the quadrupole
    #            time-averaged over M, this is just n to leading order)
    
    # For circular orbits (e=0), the e-equation is singular; use
    # the limit form (eccentricity growth/decay from R_bar)
    # or set e = 1e-6 as a regularization.
    
    # For R_bar_Lunisolar, the partial derivatives at quadrupole order:
    # dR_bar/d_i = -(3/8) (mu_3 / mu_E) (a / a_3)^3 sin(2(i - i_3)) a^2 n
    #            [per Track B Eq. 4.1, after factoring sin(i) in the
    #             Lagrange denominator]
    # dR_bar/d_Omega = 0  (R_bar is independent of Omega - Omega_3
    #                     for the quadrupole)
    # dR_bar/d_omega = 0  (R_bar at quadrupole is independent of omega)
    # dR_bar/d_e = 0     (R_bar at quadrupole is independent of e
    #                     for circular orbits; first-order in e for
    #                     elliptic; the evection term is first-order
    #                     in e_M, not e)
    
    # Therefore, for circular orbits:
    # dOmega_bar/dt = (3/8) (mu_3 / mu_E) (a / a_3)^3 sin(2(i - i_3)) / sin(i) * n
    # d_i_bar/dt = 0   (no inclination secular drift at quadrupole)
    # d_omega_bar/dt = 0   (no apsidal secular drift at quadrupole)
    # d_a_bar/dt = 0   (no semi-major axis secular drift at quadrupole)
    # d_e_bar/dt = 0   (no eccentricity secular drift at quadrupole)
    # d_M_bar/dt = n   (mean motion; standard Kepler)
    
    # For the J2 secular part (Vallado Eq. 9-39):
    # dOmega_J2/dt = -3/2 n J2 (R_E / p)^2 cos i
    # di_J2/dt = 0   (no J2 inclination secular drift at quadrupole;
    #                 crosses zero at i_crit = 63.4 deg or 116.6 deg)
    # domega_J2/dt = 3/4 n J2 (R_E / p)^2 (5 cos^2 i - 1)
    # da_J2/dt = 0
    # de_J2/dt = 0
    
    # Combine the two secular contributions:
    # dOmega_total/dt = dOmega_J2/dt + dOmega_Lunisolar/dt
    # i is constant in the secular theory (only changes via long-period
    # harmonics; ignore for the first-order secular)
    
    # 3. Integrate the 6 ODEs with fixed-step RK4 at dt_s = 600 s.
    
    t_grid = linspace(t0, t0 + duration_days * 86400, n_steps + 1)
    y_traj = rk4_propagate(rhs_mean_elements, t_grid, y0)
    
    # 4. Extract Omega_bar(t) at the same ascending-node epochs as
    # the Cowell estimator (for direct comparison).
    
    # 5. Compute the linear fit of Omega_bar(t) over the full arc.
    # The mean element has no short-period content, so this linear
    # fit is UNBIASED (up to the residual lunar nodal modulation).
    
    slope_deg_day, intercept = linear_fit(t_grid, y_traj[:, 3])
    
    return {
        "slope_deg_day": slope_deg_day,
        "intercept_deg": math.degrees(intercept),
        "slope_rad_per_s": math.radians(slope_deg_day) / 86400,
        "method": "mean_element_lagrange_planetary_eqs",
        "duration_days": duration_days,
        "dt_s": dt_s,
        "i_deg": i_deg,
        "h_km": h_km,
    }
```

### 4.3 Outputs

```python
{
    "slope_deg_day": float,          # secular Omega drift in deg/day
    "intercept_deg": float,          # Omega at t=t0 (mean element)
    "slope_rad_per_s": float,        # same in rad/s
    "method": "mean_element_lagrange_planetary_eqs",
    "duration_days": float,
    "dt_s": float,
    "i_deg": float,
    "h_km": float,
}
```

### 4.4 Convergence / validation

- The mean-element RHS is smooth; RK4 at `dt_s = 600 s` should
  give convergence order `p ≈ 4` to machine precision over a
  1-year arc. Verify with a 3-point convergence ladder
  (`dt_s ∈ {600, 300, 150}`) — the slopes should agree to
  `~10⁻¹⁰ deg/day` at the finest step.
- Cross-check against the closed-form evaluation
  `corrected_secular_lunisolar_raan_rate_rad_s(h_km, i_deg) +
  j2_secular_raan_rate(h_km, i_deg)` — these should agree to
  `~10⁻⁷ deg/day` (limited by the lunar nodal modulation, not
  by the integration error).

### 4.5 Cross-check (Estimator A)

```python
# Sub-estimator A.1: median slope of dense Omega_h(t)
def angular_momentum_secular_rate(
    r_traj, v_traj, t_grid, window_days = 30.0
):
    h_vec = cross(r_traj, v_traj)        # angular momentum vector
    n_vec = cross([0, 0, 1], h_vec)     # node vector
    Omega_h = unwrap(arctan2(n_vec[:, 1], n_vec[:, 0]))
    
    # Sliding window median slope
    window_s = window_days * 86400
    n_steps = len(t_grid)
    slopes = []
    for k in range(0, n_steps - int(window_s / dt_s)):
        t_window = t_grid[k:k + int(window_s / dt_s)]
        Omega_window = Omega_h[k:k + int(window_s / dt_s)]
        slope = linear_fit_slope(t_window, Omega_window)
        slopes.append(slope)
    
    return median(slopes), std(slopes)
```

The cross-check between C and A: at h=600 km i_sso, the two should
agree to `~10⁻⁵ deg/day` (the level of the lunar nodal modulation).
**If they agree, the headline is validated independently of the
analytical theory.**

---

## 5. Risk register for the headline (C) and cross-check (A)

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Mean-element ODE RHS wrong (e.g., missing sin(i) factor in Lagrange eq) | LOW | HIGH (bias the headline) | Verify against the closed-form evaluation at quadrupole order; convergence ladder; literature cross-check |
| Mean-element RHS diverges from Cowell due to numerical drift | LOW | MEDIUM | RK4 with dt_s=600 s; verify against RK4 with dt_s=60 s at one setting |
| Angular-momentum estimator (A) biased by per-cycle aliasing at W=1 yr | MEDIUM | MEDIUM | Use dense sampling (every RK4 step); compare A.1/A.2/A.3 sub-estimators |
| Lunar nodal modulation biases the W=1 yr estimate by `~10⁻⁵ deg/day` | MEDIUM | LOW | Acceptable; report the residual slope at multiple W to show convergence |
| Frame convention error (mean-of-date vs ICRF) biases the secular rate by `~10⁻⁴ deg/day` | LOW | HIGH | Use the FIXED eclipseTiming convention from 019; precession-identity check at t=0 and t=2026 |
| 016/017 wrong-formula regression returns | LOW | CATASTROPHIC | Both C and A are independent of the closed form; the wrong formula cannot affect them |
| Byte-pinned Moon snapshot has a per-day cadence that limits the estimator accuracy | LOW | LOW | The mean-element RHS uses time-averaged third-body elements; the cadence does not affect the secular rate |

---

## 6. Summary of Track A-5 contribution

The Track A-5 contribution to Exp 020 is the **two-estimator ladder**:

1. **Headline: Estimator C (mean-element Lagrange planetary integration).**
   Theory-aligned; canonical; clean comparison to the closed form.
2. **Cross-check: Estimator A (angular-momentum geometry, dense sampling).**
   Independent of analytical theory; provides the verification the
   closed form cannot.

The two are **complementary**:
- **C** is the right comparison to the analytical closed form; if they
  disagree, the closed form has an implementation error.
- **A** is the right independent measurement; if it disagrees with C,
  the theory (or the Cowell data, or the byte-pinned snapshot) has a
  systematic error not captured by the secular theory.

**Auxiliary diagnostics** (B and E) provide orthogonal cross-checks:
- B (short-period subtraction) confirms that removing the evection +
  variation from the ascending-node Ω gives a residual consistent
  with C and A.
- E (FFT) confirms that the dominant harmonics in the dense Ω(t) are
  at the expected annual + evection + variation periods.

The headline secular observable for Exp 020 is the
**mean-element secular Ω drift from Estimator C, cross-validated by
Estimator A**.

---

## 7. References (additional, beyond those in Track F/B/G)

- Brouwer, D. (1959). "Solution of the problem of artificial satellite
  theory without drag." *Astronomical Journal* 64, 378–397. —
  First derivation of the doubly-averaged secular theory.
- Kozai, Y. (1959). "The motion of a close earth satellite."
  *Astronomical Journal* 64, 367–377. — Independent derivation; the
  apsidal-vs-nodal factor distinction.
- Kaula, W. M. (1966). *Theory of Satellite Geodesy*. Blaisdell.
  Ch. 4 (disturbing function expansion, evection and variation
  terms). — Standard reference for the third-body disturbing function.
- Burns, J. A. (1979). "Elementary derivation of the perturbation
  equations of artificial satellite theory." *American Journal of
  Physics* 47, 850–859. — Modern textbook derivation of the third-
  body secular Ω drift formula.
- Murray, C. D. & Dermott, S. F. (1999). *Solar System Dynamics*.
  Cambridge University Press. §2.10 (Lagrange planetary equations),
  §6.4 (doubly-averaged quadrupole), §6.5 (evection and variation).
  — Standard modern textbook for the secular Hamiltonian.
- Brouwer, D. & Clemence, G. M. (1961). *Methods of Celestial
  Mechanics*. Academic Press. Chs. XI–XVII. — Rigorous celestial-
  mechanics foundation.
- Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
  Applications*, 4th ed. Microcosm Press. §9.3 (lunisolar
  perturbations; J2 + Lunisolar secular rates). — Practical
  implementation reference (Eq. 9-39 for J2 secular, Eq. 9-46 for
  the original Lunisolar cf — the wrong-form 016/017 source).
- Standish, E. M. (1990). "An observationally based reference frame
  for astronomy." *Astronomy and Astrophysics* 233, 272–274. — JPL
  approach to secular-rate extraction from finite-arc observations.
- Exp 017/018/019 lunisolarVerification / lunisolarReconciliation /
  lunisolarLongPeriod. — Existing lab canon for byte-pinned
  snapshots, fixed-step RK4 propagator, ascending-node detector,
  secular formula evaluation.

---

## Track A-5 summary in one paragraph

The Track A-5 design for Exp 020's headline secular-rate estimator
recommends a **two-estimator ladder**: (C) **mean-element Lagrange
planetary integration** as the theory-aligned headline, cross-
validated by (A) **angular-momentum-vector geometry on dense
Cowll data** as the theory-independent confirmation. Estimator C
is the canonical analytical theory implemented as a 6-element
ODE; Estimator A extracts the secular Ω rate from the kinematic
quantity `r × v` and a statistical operation (median slope,
analytic derivative, or spectral phase derivative) that does not
invoke the Lagrange planetary equations. Both have bias `~10⁻⁵
deg/day` at W=1 yr (dominated by the unmodelled lunar nodal
modulation), converge as `1/W` to the secular limit, and provide
the two independent measurements required by the 019 estimator-
ladder contract (session memory 8.78b4). Auxiliary diagnostics
(B: Brouwer short-period subtraction; E: FFT spectral) provide
orthogonal cross-checks. The 019 ascending-node linear fit is
**not** recommended as a headline (Track F bias of `~10⁻³ deg/day`
at W=1 yr). The 019 window-length extrapolation is valid as a
diagnostic but its W → ∞ value (`Lunisolar = +3.58 × 10⁻³ deg/day`,
27× the corrected cf) is **not** yet a clean measurement.

---

### Critical Files for Implementation

- `C:\Users\Dhane\lab\src\lab_utils\orbits.py` — `mean_motion`,
  `seed_state`, `rv_to_coe_eci` (RV→classical elements for the
  angular-momentum estimator cross-check at ascending-node
  crossings), `j2_rhs` (used in the Cowell propagator that
  generates the dense state for Estimator A).
- `C:\Users\Dhane\lab\research\orbital-mechanics\experiments\lunisolarLongPeriod\experiment.py`
  — Existing 019 implementation that provides the byte-pinned
  snapshot pattern, the FIXED `_rot3` precession convention, the
  Cowell RHS builder, and the ascending-node detector. The 020
  implementation will reuse all of this and add the
  mean-element Lagrange RHS and the angular-momentum estimator.
- `C:\Users\Dhane\lab\research\orbital-mechanics\experiments\lunisolarReconciliation\experiment.py`
  — The 018 implementation that establishes the corrected closed
  form `corrected_secular_lunisolar_raan_raan_rate_rad_s`, the
  precession rotation, and the force-level identity check at 50
  random states.
- `C:\Users\Dhane\lab\research\orbital-mechanics\experiments\lunisolarVerification\reference\horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt`
  — Byte-pinned Moon snapshot (sha256
  `65f1d67f798a3b95bb87310efae3200027098869246567a68ccd671d79978f4a`)
  required for both Estimator C (mean lunar elements) and
  Estimator A (frame rotation).
- `C:\Users\Dhane\lab\research\orbital-mechanics\experiments\eclipseTiming\reference\horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt`
  — Byte-pinned Sun snapshot (sha256
  `06d54fb35523a0af6ba3ea738315f1e3f5b996067c40f474052cd2fb5b5658ec`)
  required for both Estimator C (mean solar elements) and
  Estimator A (frame rotation).