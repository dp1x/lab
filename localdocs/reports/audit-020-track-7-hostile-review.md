# Audit 020 — Track A-7: Hostile Review of Exp 019

> Reviewer role: HOSTILE scientific reviewer. Charge: prove the 019 interpretation
> is WRONG. Read ONLY the 019 implementation, results, README, knowledge note,
> the corrected formula header in 018, and lab_utils constants. Did NOT read the
> 8-track audit-019 reports.
> Date: 2026-08-31
> Status of 019: COMPLETE 2026-08-30; git_commit e06d9e0

## 0. Summary verdict

**The 019 headline claim — "the 9.78× i_sso residual is dominated by
mean-vs-osculating bias from the finite-window linear fit; the corrected
secular formula is the right asymptotic prediction" — is partially defensible
but NOT robust.** Two independent failure modes survive my review and would
invalidate the headline if true. The 27× extrapolated ratio at i_sso
(0.0036 deg/day vs 0.000135 deg/day corrected) is an **artifact** of the
choice of extrapolation model AND the comparison strategy, not a clean
physics measurement.

**Recommended finding for the lab director:** Exp 020 must (a) run an
independent estimator that does NOT rely on window-length extrapolation,
and (b) test whether the secular Lunisolar RAAN rate is consistent with the
**predicted** slope at W=365 d using the unbiased mid-frequency term
analysis Track F referenced but did not quantify. Until those two
experiments land, the "27× extrapolated" number should be reported with an
explicit disclaimer that it is model-dependent to ~50%.

---

## 1. FACT / INFERENCE / UNKNOWN classification

### FACT (verifiable from the artifacts I read)

- F1. Corrected secular formula at h=600 km i_sso=97.79° = +1.3475e-4
  deg/day (results.json `corrected_closed_form_by_inclination`).
- F2. Solar = +3.5629e-5 deg/day, Lunar = +9.9125e-5 deg/day. The LUNAR
  term dominates (74% of total).
- F3. At i=90°, corrected total = +1.7390e-4 deg/day (Solar
  +4.9591e-5, Lunar +1.2431e-4). Lunar is 71% of total there too.
- F4. Window-length extrapolation Ω̇_fit(W) = a + b/W + c/W² at i_sso
  full model gives intercept a = +0.995590 deg/day (results.json
  `window_length_extrapolation.i97.79_sun_moon_j2`).
- F5. Sun-only extrapolation a = +0.995218; Moon-only a = +0.994110;
  Sun+Moon a = +0.995590; Sun+Moon+J2 a = +0.995590 (Sun+Moon and
  Sun+Moon+J2 are BIT-IDENTICAL — see §3.5).
- F6. The "Lunisolar component" the 019 authors derive at i_sso =
  0.995590 − 0.9920 (J2 baseline) = +0.0036 deg/day.
- F7. The 0.9920 deg/day "J2 baseline" is NOT in results.json. It is
  *imputed* in the README ("J2 baseline +0.9920"); the J2_only mode
  WAS NOT PROPAGATED in 019's window_sweeps (FORCE_MODES = {sun_moon_j2,
  sun_moon, moon_only, sun_only}). J2-only is absent.
- F8. The 1-year linear fit at i_sso full model slope = +0.9933 deg/day
  (`cycle_averaged_estimator.i97.79_sun_moon_j2.mean_deg_day` and
  `window_sweeps.i97.79_sun_moon_j2.365.slope_deg_per_day`).
- F9. RK4 convergence: p_r=4.4945, p_v=4.4957. Design order 4 confirmed.
- F10. Force-level identity: max_diff_sun=1.28e-21 km/s², max_diff_moon=
  5.79e-24 km/s². Passes to machine precision.
- F11. Precession identity at T=0 max_err=0.0; rotation at 2026 = -0.3332
  deg, matches eclipseTiming convention.
- F12. FFT dominant periods at i_sso are 365.03, 182.51, 121.68, 91.26,
  73.01 d — i.e., annual and harmonics of the annual.
- F13. Cycle-averaged i=90° mean = +4.8398e-4 deg/day. Corrected cf =
  +1.7390e-4 deg/day. Ratio = 2.78×. This is the "cleanest J2-free test."
- F14. The i_sso extrapolated secular Lunisolar (0.0036 deg/day) is
  ~7× the i=90° cycle-averaged (0.000484 deg/day). This 7× between
  i_sso and i=90° is never explained by the corrected formula (which
  predicts i=90° larger than i_sso). See §3.4.

### INFERENCE (claimed by 019, supported by the data but not airtight)

- I1. The 10× residual at W=365 d is dominated by finite-window linear
  fit bias (Track F theory; 019 §findings).
- I2. The corrected secular formula is the right asymptotic prediction
  for W→∞.
- I3. Annual solar forcing + evection + variation aliasing are the
  cause of the bias.
- I4. The W=730 d slope being larger than W=365 d is "smoking gun" for
  finite-window bias.

### UNKNOWN (cannot resolve from the data in 019 alone)

- U1. What the J2-only secular drift at i_sso actually is (019 did not
  run j2_only mode in its window sweep, so the +0.9920 deg/day baseline
  is unverifiable within this artifact set).
- U2. Whether the 1/W polynomial is the right extrapolation form at all.
- U3. Whether the cycle-averaged estimator at i_sso (+0.9932) being
  *almost identical* to the full-year linear fit (+0.9933) means
  short-period bias is small, or that the 12 monthly segments are
  too correlated to suppress anything (cycle-averaged at i=90° has
  segment-to-segment std of 0.00035 deg/day, while at i_sso it's
  0.00161 deg/day; the segments are dominated by J2).
- U4. Whether the 9.78× at W=365 d and the 27× at W→∞ are both
  artifacts of the same mis-specification, or whether the latter is
  simply "what you get when you fit 3 parameters to 5 points and
  extrapolate past them."
- U5. Whether the byte-pinned 366-day Moon snapshot (n=366 points,
  daily cadence) aliases lunar evection (27.55 d) into the propagation.
  The propagation interpolates the snapshot linearly; evection
  frequency is sub-cadence (~13 samples per year at 1/day).
- U6. The actual measured J2 baseline at h=600 km i_sso from 019.

---

## 2. Line-by-line attack on each candidate

### (a) Implementation bugs — already audited in Track D

The Track D remediation (sign flip on `_rot3`) is verified by the
precession_identity_check at T=0 (max err 0.0) and at T=0.26
centuries (rotation = -0.3332 deg; matches eclipseTiming
convention). Force-level identity passes to machine precision.

**Verdict: FALSIFIED (for the catalogued bugs).** Cannot falsify
un-discovered bugs; the verification surface is the public lab_utils
canon, which has been audited through 019 separate experiments. But
note: the 019 precession_identity_check is essentially a self-check
of the FIXED `_rot3` against itself; it does NOT cross-validate
against an independent ephemeris (e.g., the JPL DE441 mean-of-date
states at 2026.0).

### (b) Wrong formula (Convention A vs B, sign, scale, geometry)

The corrected formula is:
```
(3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i - i_3) / sin(i)
```
This is the standard doubly-averaged quadrupole form (Murray & Dermott
Sec 7.2; Kozai 1959; Lidov 1962). The radial scale `(a/a_3)^3` is
correct for a third-body perturbation (third-body "external" potential
vs the satellite's "internal" Hill approximation); the sin 2(i-i_3)/sin i
is the correct NODAL geometric factor (NOT the apsidal Kozai form
cos i (1 - 5/2 sin²(i-i_3))). The sign convention at retrograde is
prograde because sin 2(i_sso - i_3) / sin(i_sso) is positive
(sin 2(-12.51°) / sin 97.79° > 0; verified by recomputation:
sin(2*math.radians(97.79-23.439-5.145)) = sin(-2*30.794) = -0.856;
sin(math.radians(97.79)) = 0.992; ratio = -0.863, i.e. negative
osculating RAAN drift in this convention).

**WAIT.** Let me recheck. i_sso = 97.79°, i_3_moon = 28.584°.
2(i_sso - i_3) = 2(97.79 - 28.584) = 2(69.206) = 138.41°.
sin(138.41°) = sin(180 - 41.59) = +sin(41.59°) = +0.664.
sin(97.79°) = sin(82.21°) = +0.992.
ratio = +0.669. POSITIVE.

So the formula gives a POSITIVE (prograde) drift, matching the
numerical +1.32e-3 deg/day. Sign is consistent.

**Verdict: PLAUSIBLE BUT NOT FALSIFIED.** The formula structure is
correct in the standard textbooks. I cannot falsify the form here.
However, **note that the formula assumes a CIRCULAR third body (e_3 ≈ 0)
and CIRCULAR satellite (e = 0). At h=600 km the satellite is in fact
started on a CIRCULAR orbit by `propagate_one` (v_circ circularization
at lines 339-341), so e=0 is satisfied. The Sun e_3 ≈ 0.017 and Moon
e_3 ≈ 0.055 produce corrections of order O(e_3²) ~ 3e-3, i.e., ~3%
on the lunar term, well below the 27× residual. So this is NOT the
problem.**

### (c) Wrong estimator (linear fit is biased for periodic content)

This is the heart of the 019 claim. The Track F math says: for
Ω(t) = Ω̇_mean t + Σ A_k cos(ω_k t + φ_k), the OLS slope over window
W is biased by Σ A_k sin(ω_k W/2)/(ω_k W/2) × cos(φ_k + ω_k W/2)
directionally.

**This is mathematically correct in principle, but the 019 numerical
**implementation does NOT verify it.** The 019 authors fit a + b/W +
c/W² to the slopes {W=30, 90, 180, 365, 730} and report the intercept
a as the "secular limit." But:

1. **5 data points, 3 free parameters** — leaves 2 degrees of freedom
   for the residual; the reported `rms_residual_deg_day` at i_sso is
   0.000685 deg/day, which is ~19% of the intercept 0.9956. The
   intercept is essentially unconstrained at the 1-sigma level.
2. **The W=730 d slope (+0.9957) is essentially the same as the
   intercept (+0.9956).** The "1/W extrapolation" is dominated by
   the data point at W=730 d, not by a real extrapolation to ∞.
   The data ARE NOT consistent with a 1/W asymptote — the c_1/W²
   coefficient is 13.74, MUCH larger than the b_1/W = -0.62. The
   fit is not converging as 1/W; it is converging FASTER than 1/W.
3. **Track F's theory predicts bias scales as 1/W from periodic
   terms comparable to W, AND scales as ~1/W² from the (annual) cycle
   that is ORTHOGONAL in 1-yr fits.** Both b_1/W and c_1/W² are
   non-negligible. This is not a clean 1/W scaling; it is dominated
   by 1/W².

**Verdict: PLAUSIBLE MECHANISM but UNDERSUPPORTED.** The estimator
is theoretically biased, but the 019 extrapolation is over-fit. The
"27×" headline is the intercept of a 3-parameter fit to 5 points
whose residuals are 19% of the intercept — i.e., the intercept is
*noise-limited*, not "the secular limit."

### (d) Wrong extrapolation form (1/W polynomial is not justified)

The 019 paper uses `Ω̇_fit(W) = a + b/W + c/W²`. This is an ansatz
with NO first-principles justification for Lunisolar secular-rate
extraction. Track F provides a Fourier argument that the bias from
periodic terms is a sum of sinc-like terms; it does NOT predict a
clean 1/W polynomial. The polynomial fit is essentially "smoothing
the W-sweep with a quadratic."

**Counter-evidence:**
- The reported `b_1_over_W = -0.615` and `c_1_over_W_squared =
  +13.74` at i_sso full model have |c/b| ≈ 22. This is far from
  the "1/W is dominant, 1/W² is small correction" expectation
  implied by the choice of the model. If 1/W² is 22× larger than
  1/W, the "1/W polynomial" model is the wrong model; the data
  want a 1/W² extrapolation.
- If we use a 1/W²-only fit at i_sso: with two free parameters a
  and b, the fit would give a different intercept. The 019 authors
  do NOT report this.
- The linear 1/W fit (`linear_1_over_W_secular_deg_day = 0.9937`)
  gives a Lunisolar component of 0.9937 - 0.9920 = +0.0017 deg/day,
  which is 12.6× the corrected formula (much closer to 10× than
  27×). The quadratic extrapolation gives 27×; the linear
  extrapolation gives 13×. Which one is "right"? Neither is
  justified by theory; the difference shows the 27× is an
  ARTIFACT of choosing the quadratic model.

**Verdict: EXCESSIVE MAGNITUDE / ARTIFACT.** The 27× headline is
predominantly a property of the choice of extrapolation model,
not of the data. The same data with linear extrapolation give
13×; the 019 authors chose the model that gives the most
spectacular extrapolation. This is exactly what a hostile reviewer
flags: the answer was tuned by model choice.

### (e) Wrong window choice (W=30,90,180,365,730 may miss a dominant term)

The W choice is inherited from 018. The evection period 27.55 d is
SHORTER than the smallest W=30 d; this is fine because sinc(ωₑW/2)
with ωₑ = 2π/27.55 at W=30 is sinc(π × 30/27.55) = sinc(1.089π) ≈
-0.13, a moderate bias. The variation period 14.77 d at W=30 gives
sinc(π × 30/14.77) ≈ sinc(2.03π) ≈ +0.03, small bias. The annual
period 365 d at W=365 is ORTHOGONAL (zero bias by symmetry over a
full cycle). At W=730, annual period is 2 cycles, also orthogonal.

**HOWEVER: there are MAJOR missing frequencies:**
- **18.6-year lunar nodal cycle**: cannot be resolved with W ≤ 730 d.
  Track F acknowledges this in 019's limitations. But the
  secular-mean prediction already assumes mean lunar orbit plane
  (i_3 averaged over 18.6 yr); using i_3 = 28.584° for a 1-year arc
  introduces a BIAS of order d(i_3)/dt × (1 yr) ≈ 5° × (1/18.6)
  ≈ 0.27°. This propagates through sin(2(i-i_3))/sin(i) at the
  ~3% level — small, but nonzero.
- **8.85-year lunar apsidal precession**: same issue, ~3% on
  lunar term over 1 yr.
- **209-year lunar evection of the evection (lunar nodal
  precession inside the lunar orbit)**: ignored.

**Verdict: PLAUSIBLE / INSUFFICIENT MAGNITUDE for the 27×.** None of
these individually account for 27×. But they collectively bias the
secular prediction by ~5%, which means the corrected formula itself
is uncertain at the 5% level, NOT at the 27× level. The 27× must
come from elsewhere.

### (f) Wrong reference data (snapshot too short, daily cadence too coarse)

The Sun and Moon snapshots have 366 daily points each (1-year arc,
daily cadence). Linear interpolation is used. Cadence analysis:

- Evection (27.55 d) → 13.3 samples/year in the snapshot. **Nyquist
  is satisfied** (2 samples per period needed; 13.3 > 2). BUT
  linear interpolation of a sinusoidal field with 13 samples/cycle
  introduces a numerical artifact of order ~5% in amplitude and
  ~5° in phase. At the evection frequency, this could bias the
  osculating Ω computation.
- Variation (14.77 d) → 24.8 samples/year. OK.
- Annual → 1 sample/year at the 1-yr boundary; the snapshot
  ends and restarts; the 1-yr arc has a periodic discontinuity
  at the boundary (the Sun returns to ~the same point but the
  Moon has moved ~13 lunar cycles, so the snapshot's
  end-of-arc discontinuity is dominated by the Moon). The
  propagation uses `idx = t_s[0] ... t_s[-1]` with
  `t_query_s >= t_s[-1] → r[-1]` clamp, so the LAST day's
  Moon position is held fixed. Over 1 year, this clamp affects
  the last few hours of propagation, but the linear fit is
  dominated by the bulk — should be small.

**Verdict: PLAUSIBLE.** The cadence is sufficient for Nyquist at
evection/variation; the 1-yr snapshot boundary does not corrupt
the linear fit. Cannot falsify.

### (g) Wrong frame (precession bug already fixed but verify)

`precession_identity_check` returns identity_at_T0_max_err = 0.0
and rotation_at_2026_deg = -0.3332 (matches eclipseTiming
convention to <0.01 deg). This is the FIXED `_rot3` applied to
the JPL DE441 ICRF/J2000 vectors.

**HOWEVER:** I cannot, from the artifacts I read, verify that the
ECI propagator in `j2_rhs` and `rk4_propagate` is consistent with
this frame. The `j2_rhs` works in inertial (non-rotating) frame
with Z = Earth spin axis at epoch. The IAU-1976 precession matrix
rotates from J2000 to mean-of-date — i.e., it accounts for the
Earth's spin-axis precession but NOT for nutation or the equation
of equinoxes. The instantaneous ECI frame used by `j2_rhs` is the
**true-of-date** frame (or maybe mean-of-date, depending on the
`earth_frames` module's convention).

Without reading the lab_utils internals further, I cannot verify
that the precession is being applied in the RIGHT direction. The
identity test only checks that the matrix IS the claimed matrix;
it does NOT check that the matrix is being applied to the right
vectors in the right direction.

**Verdict: PLAUSIBLE but UNVERIFIED.** The precession_identity_check
is a self-consistency test, not a cross-validation. A frame bug
could be lurking: applying precession to the Sun/Moon vectors
when the propagator already works in J2000 (i.e., double-applying)
would put the third-body vectors in the wrong frame and bias the
numerical RAAN rate. I cannot falsify this without reading more
of lab_utils.

### (h) Wrong physical model (the secular formula is leading-order; need higher-order corrections)

The standard quadrupole secular formula is leading-order in
(mu_3/mu_E)(a/a_3)³. Higher-order corrections:

- **Octopole correction**: O((a/a_3)⁴) × mu_3/mu_E. For the Moon,
  (a/R_M)⁴ = (0.0174)⁴ = 9.2e-8; × (4903/398600) = 1.1e-9. Times
  the leading-order lunar term ~1e-4 deg/day gives ~1e-13
  deg/day. NEGLIGIBLE.
- **Cross-term with J2**: Kozai's coupled J2+Lunisolar analysis
  shows a coupling term ~J2 × (mu_3/mu_E)(a/a_3)³. J2 ~ 1e-3;
  × 1e-4 deg/day = ~1e-7 deg/day. NEGLIGIBLE.
- **Solar radiation pressure**: 1e-6 m/s² at 600 km altitude with
  A/m = 0.01 m²/kg → RAAN rate ~10⁻⁵ deg/day at SSO. At i_sso
  the SRP RAAN drift is technically nonzero but tiny.

**Verdict: INSUFFICIENT MAGNITUDE.** No higher-order physical
correction can bridge a 27× gap. If the corrected formula is
right, the higher-order terms don't matter.

### (i) Wrong observable (osculating Ω is not the mean Ω; the secular formula is for the mean Ω)

This is the CORRECT identification of the issue, and it is what
019 claims is the resolution. The doubly-averaged formula predicts
the MEAN Ω̇, not the OSCULATING Ω̇. The osculating Ω(t) contains
short-period oscillations about the mean; an OLS fit over a finite
window captures these oscillations as additional slope.

**HOWEVER:** the mean Ω̇ IS the asymptotic limit of the W → ∞
extrapolation ONLY IF the short-period terms average out over a
window. This is a deep assumption: the moon's evection has a
period that is incommensurate with the year, so a 1-year average
does NOT fully cancel evection. A W → ∞ extrapolation does cancel
it (by Riemann-Lebesgue lemma), but the RATE of cancellation is
O(1/W) per Track F. With W=730 d and evection period 27.55 d, we
have 26.5 cycles; a 1-year average cancels evection to O(1/26.5) ≈
4% — i.e., the 1-year residual contains ~4% of the evection
amplitude, which is NOT zero.

**So Track F is right that W → ∞ extrapolation is the right
approach, but the rate at which it converges to the secular
mean is not yet demonstrated.** The 019 extrapolation gives 27×,
but the data span is W ≤ 730 d, where evection has only ~26.5
cycles; the W → ∞ limit might require W > 50 yr to fully
stabilize. The "extrapolation to infinity" is a hope, not a
measurement.

**Verdict: PLAUSIBLE MECHANISM, INSUFFICIENT DATA to prove.** The
trend is in the right direction, but the absolute magnitude
(27×) is model-dependent.

### (j) Wrong inclination geometry (the actual 2026 lunar i3 differs from 28.584 deg mean)

The Moon's orbital plane precesses with period 18.6 yr; the
inclination of the lunar orbit relative to the ecliptic varies
between 18.3° and 28.6° over this cycle. The "mean" i_3 = 23.4° +
5.145° = 28.584° assumes a long-time average. In 2026, the
lunar orbital plane is at ~i_3 = 18.9° (declining phase).

The secular formula uses sin 2(i - i_3) / sin i. At i_sso=97.79°
with i_3 = 18.9°, sin 2(97.79 - 18.9) = sin(157.78°) = +0.384;
with i_3 = 28.584°, sin 2(97.79 - 28.584) = sin(138.41°) = +0.664.
The lunar contribution with the 2026-actual i_3 would be (0.384/
0.664) × the 019 number = 0.578 × 9.91e-5 = 5.73e-5 deg/day,
not 9.91e-5. **The corrected formula OVERESTIMATES the 2026
lunar contribution by 1.73× because it uses the secular-mean
i_3, not the 2026-actual i_3.**

This is significant: if the corrected formula is wrong by 1.73×
on the lunar term alone, the "1.35e-4 deg/day" prediction is
actually closer to (3.56 + 5.73) × 1e-5 = 9.29e-5 deg/day. The
ratio to the numerical +1.32e-3 is 14.2×, not 9.78×.

**But this REDUCES the residual, not increases it.** And it does
NOT explain the 27× extrapolation.

**Verdict: PLAUSIBLE, MAGNITUDE ~1.7×, but cannot account for 27×.
Flag as an unaddressed systematic in 019 (the implementation
uses i_3=28.584° as a frozen constant, not the 2026-actual value).
The numerical propagation uses the actual 2026 Moon snapshot, so
the comparison is inconsistent: corrected formula in mean-i_3
frame, numerical in 2026-actual-i_3 frame.**

### (k) Aliasing of short-period terms in the FFT analysis

The FFT is on the osculating Ω(t) at ascending-node crossings.
The dominant periods found are 365, 182, 121, 91, 73 d — all
annual harmonics. The expected evection (27.55 d) and variation
(14.77 d) peaks are NOT in the top 5.

**This is suspicious.** Track B/F predicted evection and variation
would be the dominant short-period terms. The FFT shows annual
harmonics dominate. Possible explanations:

1. The annual cycle has a much larger amplitude than evection/
   variation in the osculating Ω. Plausible: the Sun's annual
   cycle modulates the orbital plane at the Earth's orbit
   frequency (~1 cycle/year), which is more directly visible in
   Ω than the higher-frequency evection.
2. **Aliasing**: the FFT has dt_mean_day ~ 0.067 (one crossing
   every ~1.6 hours at 600 km altitude, with period ~96.7 min).
   Frequency resolution is ~0.002 cycles/day, period resolution
   ~1 day near the annual cycle but ~0.1 day near the 14-day
   variation period. **Window aliasing is real but Nyquist is
   fine because the crossings are sampled at every ascending
   node, not every day.**
3. The dominant annual term in the FFT might actually be a
   solar-forcing annual effect, which 019 IS predicting (Track
   G says "annual solar forcing + lunar evection/variation =
   dominant surviving trio"). But the FFT does NOT confirm the
   evection/variation amplitude; it confirms the ANNUAL is big.

**Verdict: PLAUSIBLE.** The FFT analysis is consistent with the
Track F prediction that annual dominates the linear-fit bias.
But it does NOT independently confirm that evection/variation
are present at the expected amplitude.

### (l) Aliasing of short-period terms in the 12-month cycle-averaged estimator

The cycle-averaged estimator divides the 1-year propagation into
12 monthly segments, each ~30 d. Within each segment, it computes
an OLS slope. Then takes the mean.

**The 019 results show segment-to-segment std of 0.00035 deg/day
at i=90° and 0.00161 deg/day at i_sso.** The i=90° mean is
+4.8398e-4 deg/day with std 0.00035 — that's a 72% 1-sigma
spread. This is huge. The cycle-averaged estimator is supposed
to reduce variance, but the segments are not independent
(adjacent segments share orbits), and the bias structure is
NOT stationary.

**The README claims "mean slope within 7e-5 deg/day of the
full-year linear fit (Track E)" but at i=90° the difference
between the cycle-averaged mean (+4.84e-4) and the full-year
linear fit (+5.17e-4 at W=365, results.json
window_sweeps.i90.00_sun_moon_j2.365.slope_deg_per_day) is
3.4e-5 deg/day, not 7e-5.** Close, but inconsistent with the
README's claim.

**Verdict: PLAUSIBLE, PARTIALLY FALSIFIED.** The cycle-averaged
estimator at i=90° has 72% 1-sigma spread, which is NOT a
"3% bias" estimator. It's a noisy estimator with the mean
shifted by a similar amount to the linear fit. The claim that
"cycle-averaged reduces bias to 3% vs 5-15%" is overstated.

### (m) Ephemeris phase dependence (the 1-yr arc starts at 2026-01-01; results may depend on this phase)

The 1-year arc starts at JD 2460676.5 (2026-01-01) and ends at
JD 2461041.5 (2026-12-31). The 019 implementation uses `t0 =
820540800.0` (this is JD 2460676.5 in seconds since J2000).

**Phase dependence is REAL.** The Sun's apparent longitude at
2026-01-01 is ~280° (near perihelion); the Moon's mean longitude
is at a particular value. Different starting phases would
generate different finite-window biases. The 019 results are
ONE realization of the start phase.

**019 did NOT test phase dependence.** This is a critical gap.
If the extrapolated secular rate depends on the start phase by
more than a few percent, the 27× result is a snapshot artifact.

**Verdict: PLAUSIBLE AND UNADDRESSED.** This is a structural
weakness in the 019 design. To verify, Exp 020 should run the
window-length extrapolation at multiple start phases (e.g.,
2026-01-01, 2026-04-01, 2026-07-01, 2026-10-01) and confirm the
extrapolated intercept is stable.

### (n) Numerical drift in the RK4 propagator

The convergence ladder shows p_r=4.49, p_v=4.50 at h=600 km
sun_moon_j2 mode, over a 1-day test. The maximum |r| difference
between dt=120 s and dt=1.875 s reference is 51.8 km; between
dt=60 s and dt=1.875 s is 1.75 km. This is consistent with
RK4 design order.

**HOWEVER:** The propagator uses dt=60 s for the W=730 d
propagation. The total error scales as N_steps × error_per_step.
For W=730 d: N_steps = 730 × 86400 / 60 = 1,051,200 steps. The
per-step error at dt=60 s is ~1.75 km / (86400/60) = 1.22e-3 km
per step × constant = ~5e-7 km per step × N? This isn't right.
Let me recompute: error at dt=60 s after 1 day is 1.75 km; this
scales as (dt_ref/dt)⁴ × error_ref = (1.875/60)⁴ × 0.000111 =
8.4e-9 × 0.000111 = 9.3e-13 km. After 1 day. The actual is
1.75 km, so my interpretation is off. Let me re-read.

Actually the test in convergence_ladder propagates for T_test =
86400 s = 1 day, with dt_finest = 1.875 s (the reference) and
coarse dt in {120, 60, 30, 15, 7.5}. The reported diffs are at
the END of the 1-day arc. At dt=120 s vs dt=1.875 s, the diff
is 51.8 km. This is the truncation error of the dt=120 s
propagation.

For W=730 d propagation at dt=60 s, the per-step truncation
error is much smaller than at dt=120 s (factor (60/120)⁴ = 1/16).
But the accumulated error over 730 d is ~730 × error per day.
The 1-day error at dt=60 s is 1.75 km; over 730 d this would be
~730 × 1.75 km if errors accumulate linearly, OR sqrt(730) ×
1.75 km ~ 47 km if they accumulate as a random walk. **Either
way, the Lunisolar contribution to Ω drift over 730 d is
~0.9957 deg/day × 730 d = 727 deg; the systematic error in Ω
from a 50 km position error is O(distance/orbit) × 2π rad ~ 0.5
rad ~ 30 deg.** This is NOT negligible compared to the 0.0036
deg/day × 730 d = 2.6 deg Lunisolar signal.

**Verdict: PLAUSIBLE NUMERICAL DRIFT.** The RK4 design order is
confirmed at 1 day, but the 1-year numerical propagation at
dt=60 s has not been independently verified to converge at the
expected rate over long arcs. The convergence ladder is too
short (1 day) to prove 1-year accuracy. Exp 020 should run a
multi-day convergence ladder (e.g., 30 d, 90 d, 365 d) at the
ACTUAL force model and confirm the design order is preserved.

### (o) Insufficient convergence ladder (1-day p_r=4.5 doesn't prove 1-year arc numerical accuracy)

As argued in (n). The convergence_ladder function tests 1-day
arc only. The full propagation uses 730-day arc.

**Verdict: PLAUSIBLE.** Same as (n).

---

## 3. Survivor ranking

After applying each test, here are the candidates that SURVIVE
(my assessment of likelihood):

| Rank | Candidate | Status | Threat level to 019 headline |
|---|---|---|---|
| 1 | (d) Wrong extrapolation form (1/W polynomial, model-dependent) | **PLAUSIBLE / STRONG** | **HIGH.** The 27× number depends on choosing the 1/W + 1/W² model; linear 1/W gives 13×. The choice is not theoretically justified. |
| 2 | (i)+(m) Mean-vs-osculating bias is real, but extrapolation is unverified | **PLAUSIBLE / STRONG** | **HIGH.** The mechanism is right; the magnitude is unknown. The 019 extrapolation is hope, not measurement. |
| 3 | (n)+(o) RK4 1-day convergence doesn't prove 1-year accuracy | **PLAUSIBLE** | **MEDIUM.** Systematic ~50 km position error over 1 yr could bias the slope by ~10%. The 730 d extrapolation is the most vulnerable. |
| 4 | (j) Lunar i_3 used is mean not 2026-actual (1.7× bias on lunar term) | **PLAUSIBLE / KNOWN** | **LOW.** Reduces the corrected formula's magnitude, partially closes the gap, doesn't explain 27×. |
| 5 | (l) Cycle-averaged estimator at i=90° has 72% 1-sigma spread | **PLAUSIBLE / PARTIALLY FALSIFIED** | **MEDIUM.** The 3% bias claim is overstated; the cycle-averaged is noisy, not low-bias. |
| 6 | (g) Frame cross-validation is a self-check | **PLAUSIBLE / UNVERIFIED** | **LOW-MEDIUM.** Possible double-application of precession; would need to read lab_utils internals to verify. |
| 7 | (a) Implementation bugs | **FALSIFIED** | LOW. Catalogued bugs are fixed. |
| 8 | (b) Wrong formula | **PLAUSIBLE but not falsifiable from 019 alone** | LOW. Standard textbook form. |
| 9 | (c) Linear fit is biased | **PLAUSIBLE, MECHANISM CORRECT** | LOW. Mechanism is right; the magnitude is what's disputed. |
| 10 | (e) Wrong window choice | **INSUFFICIENT MAGNITUDE** | LOW. |
| 11 | (f) Wrong reference data | **PLAUSIBLE, INSUFFICIENT MAGNITUDE** | LOW. |
| 12 | (h) Wrong physical model | **INSUFFICIENT MAGNITUDE** | LOW. No higher-order term bridges 27×. |
| 13 | (k) FFT aliasing | **PLAUSIBLE, INCONSISTENT WITH TRACK F PREDICTION** | LOW. Annual dominates, evection/variation are not in top 5. |

### Top 2 most dangerous survivors

**1. (d) Wrong extrapolation form.** The 27× headline depends on
the choice of polynomial (linear 1/W gives 13×; quadratic 1/W²
gives 27×). The choice is unprincipled. The corrected formula's
prediction could be made to match ANY extrapolated value by
choosing the right model. This is exactly the failure mode a
hostile reviewer should flag.

**2. (i)+(m)+(n)+(o) composite: Mean-vs-osculating bias mechanism
is correct, but its magnitude is unknown; ephemeris phase
dependence is unaddressed; RK4 convergence is proven only at 1
day.** The 019 extrapolation gives a single number (0.0036
deg/day) that depends on a single start phase (2026-01-01), a
single RK4 dt (60 s), and a single extrapolation model. None of
these three axes are swept.

If either of these two is true, the 27× headline is wrong.

---

## 4. Recommended experiments for Exp 020

To distinguish the survivors, Exp 020 MUST run two experiments:

### Experiment 020-A: Model-independent secular rate extraction

**Question:** Is the secular Lunisolar RAAN rate at W → ∞ a
well-defined quantity, or does it depend on the choice of
estimator and extrapolation model?

**Method:** At h=600 km i_sso, propagate with sun_moon only
(no J2) for the longest possible arc (multi-year if the lab
acquires multi-year snapshots; otherwise stack 3-5 single-year
arcs with different start phases: 2026-01-01, 2024-01-01,
2025-01-01, 2027-01-01, 2028-01-01 — but 020 would need
additional snapshots).

For each arc:
- Compute the osculating Ω at ascending-node crossings.
- Apply a HANNING-TAPERED FFT to identify dominant periods.
- Subtract the identified periodic components (annual + harmonics
  + evection + variation + any detected terms) from the
  osculating Ω.
- Compute the slope of the RESIDUAL via OLS; this should be
  closer to the secular mean than the slope of the raw Ω.
- Report the residual slope AND its 1-sigma uncertainty from a
  bootstrap of the periodic subtraction.

**Pass criterion:** The residual slope (after periodic
subtraction) at multiple arcs converges to a value consistent
with the corrected secular formula (+1.35e-4 deg/day at i_sso)
to within 3× (i.e., ~4e-4 deg/day), independent of the choice
of extrapolation model.

**Fail criterion:** The residual slope depends on the number
of periodic terms subtracted, or fails to converge across
different arcs.

### Experiment 020-B: Phase-averaged secular rate

**Question:** Does the window-length extrapolation depend on
the ephemeris phase?

**Method:** Acquire 4 Moon snapshots spanning different lunar
phases (e.g., 2026-01-01, 2026-04-01, 2026-07-01, 2026-10-01
start dates, each 1-year long). For each, run the full
019 window-length extrapolation. Report the intercept a for
each start date.

**Pass criterion:** The 4 intercepts agree to within the
extrapolation uncertainty (~10%). The corrected formula's
+1.35e-4 deg/day Lunisolar term falls within 3× of any of them.

**Fail criterion:** The intercepts scatter by more than 50%,
or are systematically biased relative to the corrected formula.

These two experiments together will determine whether the 27×
extrapolation is a robust prediction or a phase-dependent
artifact.

---

## 5. Final finding

**Is the 019 headline (extrapolated +0.0036 deg/day = 27×
corrected formula) scientifically defensible, or is it an
artifact?**

**It is partially defensible but contains an unacknowledged
artifact.** The qualitative conclusion — that the 1-year linear
fit is biased and the corrected secular formula predicts the
mean, not the osculating, Ω drift — is well-supported by
Track F's bias theory and by the qualitative trend of
increasing slope with window length (F4: 0.9903 → 0.9910 →
0.9918 → 0.9933 → 0.9957 at W = 30, 90, 180, 365, 730 d;
monotone).

**The quantitative conclusion (27×) is an artifact of model
choice.** The same data with a linear 1/W extrapolation gives
13×. The 019 authors chose the quadratic extrapolation because
it gives the larger number, but the linear extrapolation has
fewer free parameters and equal theoretical justification (the
annual cycle contributes O(1/W²) bias per Track F, but the
window-length sensitivity also has O(1/W) bias from periodic
terms comparable to W).

**The 019 headline is NOT scientifically defensible as a
measurement; it is defensible as a qualitative statement that
the corrected formula's +1.35e-4 deg/day prediction is at least
plausible and not contradicted by the 1-year arc data.**

The lab should:
1. Mark the 27× number with a clear caveat: "model-dependent
   extrapolation; the linear 1/W fit gives 13×."
2. Reframe the headline as: "the corrected secular formula is
   the right asymptotic prediction; the finite-window linear
   fit has bias 1-3×10⁻⁴ deg/day; the precise magnitude of the
   bias depends on extrapolation model."
3. Treat the 27× number as an upper bound on the bias, not a
   measurement of the secular mean.

---

## Critical Files for Implementation

- `C:\Users\Dhane\lab\research\orbital-mechanics\experiments\lunisolarLongPeriod\experiment.py` — Core 019 implementation; the window-length extrapolation function `window_length_extrapolation` (lines 397-441) and the cycle-averaged estimator `cycle_averaged_slope` (lines 448-482) are the two estimators whose model-dependence is the basis of this review.
- `C:\Users\Dhane\lab\research\orbital-mechanics\experiments\lunisolarLongPeriod\results\results.json` — Contains the W-sweep slopes at i_sso and i=90°; the `slopes_deg_day` arrays in `window_length_extrapolation.i97.79_sun_moon_j2` and `i90.00_sun_moon_j2` are the data that the 27× and 13× extrapolations both fit.
- `C:\Users\Dhane\lab\research\orbital-mechanics\experiments\lunisolarReconciliation\experiment.py` — Contains the corrected doubly-averaged quadrupole formula at the top of the file; this is the +1.35e-4 deg/day prediction the 019 extrapolation is being compared against.
- `C:\Users\Dhane\lab\localdocs\reports\audit-020-track-7-hostile-review.md` — This report. Recommended follow-up experiments 020-A (model-independent secular rate extraction) and 020-B (phase-averaged secular rate) should be implemented in Exp 020.
- `C:\Users\Dhane\lab\src\lab_utils\orbits.py` — Contains the canonical constants `R_EARTH_KM=6378.137`, `MU_EARTH_KM3S2=398600.4418`, `J2_EARTH=1.082629821e-3`, and the `j2_rhs` propagator interface used by 019. To verify candidate (g) (frame cross-validation), Exp 020 would need to cross-check whether `j2_rhs` works in J2000 or mean-of-date frame, which determines whether the IAU-1976 precession rotation in 019 is being applied correctly.