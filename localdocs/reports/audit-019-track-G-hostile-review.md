# Audit-019 Track G — Hostile Review of the 018 10× Residual

> **Track G mandate.** Act as a hostile scientific reviewer. Attempt to
> falsify every candidate explanation for the ~10× residual between the
> corrected secular formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i`
> (Track B) and the 1-year numerical linear fit at h=600 km i_sso=97.79°.
> List every plausible error or alternative explanation, attempt to
> destroy each, and report the explanations that survive.
>
> **What I read (only):** AGENTS.md; 018 README, 018 experiment.py,
> 018 results/results.json; 017 experiment.py;
> src/lab_utils/orbits.py and src/lab_utils/integrators.py. I did NOT
> read audit-018, audit-019-track-A, or any other tracks' outputs.
>
> **Frame.** The "10× residual" is the ratio
> numerical Lunisolar / corrected cf at h=600 km i_sso = 9.78
> (total), 33.7 (solar only), 1.17 (lunar only). The "i=90° null test"
> gives 2.81× (cleanest J2-free point). The window-length sensitivity
> goes 0.9903 → 0.9910 → 0.9919 → 0.9933 → 0.9958 deg/day as
> W = 30 → 90 → 180 → 365 → 730 d (monotone increasing). At
> i=90° the J2 cos(i) factor vanishes, so the i=90° measurement
> is pure Lunisolar on a non-SSO inclination.

---

## 1. Headline verdict before the per-candidate attack

The 018 conclusion is that the residual is dominated by **unmodelled
short-period terms** (evection ~27.55 d anomalistic month, variation
~14.77 d synodic half-month, lunar nodal regression 18.6 yr). The
authors are *confident* this is the answer; they rank it as the only
PLAUSIBLE candidate.

A hostile reviewer's job is to test whether this confidence is
warranted. I will attack every candidate, including ones the 018
authors may have underweighted, and force the surviving set to be
small. **The conclusion of this review is not a confirmation of the
018 attribution; it is a forced ranking of survivors.**

---

## 2. Candidate-by-candidate attack

For each candidate I list: (i) order-of-magnitude estimate at h=600 km
i_sso, (ii) sign prediction, (iii) inclination and window-length
predictions, (iv) cross-check verdict, (v) falsification status.

### a) Evection term (lunar anomalistic month ~27.55 d)

**Hypothesis.** The unmodelled evection modulation of the Moon's
geocentric position at the anomalistic month introduces a short-period
Ω oscillation whose time-average over 1 year contaminates the linear
fit and explains the ~10× residual.

**Magnitude.** Evection is the largest lunar perturbation (~1.27°
peak-to-peak geocentric, sometimes quoted up to 1.5°). For an SSO at
h=600 km, the third-body accelerations are ~10⁻¹⁰ km/s² (Sun) and
~10⁻⁹ km/s² (Moon); the evection term appears as an ~1° modulation of
the Moon's apparent direction relative to the satellite. Rough
order-of-magnitude: the evection-induced RAAN perturbation over a
lunar anomalistic month is roughly the same scale as the secular
term, perhaps larger. The annual time-average could plausibly
contribute ~10⁻³ deg/day (i.e., the observed magnitude).

**Sign.** Evection is an angular modulation, not a sign-flipping
mechanism; the time-average should add to the secular term in either
positive or negative direction depending on phase. The observed
**positive** residual (numerical > corrected) is consistent.

**Inclination dependence.** Evection is a modulation of the lunar
direction (the *i₃* in the formula); at i=90° the geometric factor
`sin 2(i−i₃)/sin i` is at a maximum (because the formula is exact
and not suppressed by cos(i)). If evection only modifies the lunar
direction but not the geometric formula structure, the i=90° ratio
should be similar to i_sso. Observed: 2.81× at i=90° vs 9.78× at
i_sso. The reduction to 2.81× at i=90° is **inconsistent** with
evection as the dominant residual at i_sso, because evection is
LUNAR-dominated and the LUNAR ratio at i_sso is only 1.17×. This is a
genuine falsification challenge.

**Window dependence.** The window sensitivity (+0.005 deg/day from
W=30 to W=730) is a 5% relative change. Evection's period is 27.55 d;
at W=30 d the fit captures at most one cycle (partial cancellation);
at W=730 d the fit averages ~26 cycles (good cancellation if the
evection is symmetric). The +0.005 deg/day monotonic increase is
*consistent* with evection averaging out at long windows.

**Solar vs lunar decomposition.** The decomposition shows solar 33.7×
over-estimated, lunar only 1.17×. Evection is a LUNAR effect; it
should affect the LUNAR residual. Observed LUNAR ratio 1.17× is
small. So **evection alone cannot account for the solar 33.7×
ratio.** Partial verdict: PLAUSIBLE for the *lunar* residual
structure; INSUFFICIENT for the solar 33.7× ratio.

**Verdict.** **PLAUSIBLE (partial).** Evection is a real, unmodelled
short-period modulation. It cannot account for the solar 33.7× ratio
on its own. Survives as a contributor to the *lunar* component but
fails as a complete explanation.

### b) Variation term (lunar synodic half-month ~14.77 d)

**Hypothesis.** The variation (evection's cousin, semi-major axis
modulation at the synodic half-month ~14.77 d) introduces a similar
short-period Ω oscillation.

**Magnitude.** Variation is ~0.66° peak-to-peak geocentric, smaller
than evection. At h=600 km the contribution is similar in structure
but ~½× evection.

**Sign.** Same as evection, positive when averaged over a non-
symmetric year.

**Inclination dependence.** Same structure as evection — lunar-only,
should not affect i=90° asymmetry differently from i_sso.

**Window dependence.** 14.77 d period; W=30 d captures ~2 cycles
(partial cancellation), W=730 d captures ~49 cycles (good average).
Same direction as evection — long windows reduce the residual.

**Solar/lunar decomposition.** Variation is lunar-only. Same
falsification as evection: cannot account for the solar 33.7×.

**Verdict.** **PLAUSIBLE (partial).** Same as evection — contributes
to lunar residual, fails as a complete explanation.

### c) Annual solar forcing (365.24 d)

**Hypothesis.** The Sun's geocentric position over a year has an
annual modulation that introduces a 365-d-period RAAN oscillation
that the 1-year linear fit cannot remove.

**Magnitude.** The Sun moves through ~360° in 1 year (it IS the
annual motion). At SSO retrograde inclinations, the solar forcing
geometry changes substantially as the Sun traverses the ecliptic.
The annual solar forcing is **the dominant unmodelled term** for
the *solar* residual.

**Sign.** The annual-averaged solar inclination relative to a SSO
orbit plane is approximately the obliquity 23.439°. The secular
formula uses this average, so the residual from the annual modulation
is approximately the difference between the time-varying solar
inclination and the obliquity, integrated over the year. This can
be either sign; the observed numerical > corrected means the
**annual average of the time-varying term exceeds the constant-
obliquity average**.

**Magnitude estimate.** Solar at h=600 km i_sso has sin 2(i_sso−ε) =
sin 2(97.79−23.439) = sin 2(74.35) = sin 148.70° ≈ 0.514. The
secular formula gives 3.56e-5 deg/day (corrected). The annual
average of the time-varying third-body term should be of the same
order, but the 1-year linear fit weights the data unevenly (more
crossings during certain phases of the Sun's motion). A 33× ratio is
*not* expected from the secular formula structure alone — it implies
the secular formula under-estimates by ~33×.

**Inclination dependence.** The annual solar forcing is the dominant
mechanism at i_sso (where the secular formula's solar contribution
is 3.56e-5 deg/day and the numerical is +1.20e-3 deg/day, a 33.7×
ratio). At i=90° the solar secular contribution (corrected) is
larger because sin(2(90−ε))/sin(90°) is at its peak; the annual
forcing should still be present. The i=90° total ratio is only 2.81×,
so the solar contribution there is small relative to the lunar.

**Window dependence.** Solar annual forcing at 365.24 d is a single
cycle per year. At W=30 d the fit may be dominated by short-period
local solar position; at W=730 d it averages two cycles (good
average). The +0.005 deg/day monotonic increase is *consistent* with
solar annual forcing partially cancelling at long windows.

**Solar/lunar decomposition.** This is the prime candidate for the
**solar 33.7× ratio**. Lunar ratio 1.17× is much smaller because the
lunar secular is ~3× the solar secular in the corrected formula
(the Moon is closer despite lower mass). If the unmodelled
contribution is dominated by the solar annual forcing, the SOLAR
ratio should be large and the LUNAR ratio should be near 1. This is
**exactly what is observed**.

**Verdict.** **PLAUSIBLE (strong for solar component).** The annual
solar forcing is a strong candidate for the 33.7× solar ratio. The
i=90° ratio reduction to 2.81× is consistent if the residual at
i=90° is dominated by the lunar (where the secular formula already
captures most of the magnitude). Combined with evection (lunar),
this gives a coherent attribution.

### d) Lunar nodal precession at 18.6 yr

**Hypothesis.** The Moon's orbital plane precesses with an 18.6-year
period, which the 1-year measurement cannot resolve. Over 1 year
this looks like a slowly-varying secular drift.

**Magnitude.** The lunar nodal precession changes the Moon's
inclination relative to the ecliptic by ±5.145° at the 18.6-yr
period. Over 1 year this is a ~0.55° change — small but not
negligible.

**Sign.** The secular formula uses the mean lunar inclination 5.145°.
The actual 2026 value is the time-varying value at that date.
Depending on the phase of the 18.6-yr cycle, the time-varying value
could be slightly more or less than 5.145°.

**Order-of-magnitude.** A 0.5° change in i₃ propagates linearly into
sin 2(i−i₃). At i_sso this is a relative change of ~0.5°/sin(74.35°) ≈
1.9% — much smaller than the observed 9.78× (which is a 978%
change). So lunar nodal precession is INSUFFICIENT in magnitude at
1-year arc.

**Window dependence.** The 18.6-yr period is much longer than the
730-d window. The window-length sensitivity should NOT show a strong
trend if this is the dominant effect. Observed: +0.005 deg/day over
700 d is a real trend, suggesting some OTHER mechanism with a period
shorter than 730 d.

**Solar/lunar decomposition.** Lunar-only mechanism.

**Verdict.** **INSUFFICIENT MAGNITUDE.** The 1-year arc captures
only ~5% of the 18.6-yr cycle, producing a 1-2% effect on the lunar
secular rate. Cannot explain a 33× solar ratio, and the window
sensitivity trend argues against it being dominant.

### e) Finite-window linear-fit bias

**Hypothesis.** The 1-year linear fit of Ω(t) for a signal that
contains periodic terms is biased because the periodic terms
contribute a non-zero mean over a finite window.

**Magnitude.** This is well-known in time-series analysis: a linear
fit of a sinusoidal signal y = A sin(ωt + φ) over a window of length
W gives a biased slope that depends on the phase φ at the window
start. The bias can be of order A·(1/W) for long-period terms and
O(A·ω) for short-period terms.

**Sign.** The bias sign depends on phase. With 1-year snapshot,
the Moon's phase is fixed (the snapshot starts 2026-01-01), so the
bias is a definite number, not averaged.

**Magnitude estimate.** For evection (A_evection ~ 27.55-d period,
amplitude in Ω unknown but plausibly 0.01-0.1 deg), the 1-year fit
bias is of order A_evection × (2π/27.55) × (W_window / W_evection)².
For W=365 d, T_evection=27.55 d: bias factor ~ (365/27.55)² = 175
× smaller than A_evection itself. For A_evection ~ 0.01 deg, the
bias is ~5.7e-5 deg — comparable to the corrected secular rate.

**Inclination dependence.** The bias depends on the amplitude of
periodic terms, which depends on inclination through the third-body
geometry. At i=90° the periodic terms are at their maximum amplitude
for the lunar case (no cos(i) suppression), but the bias itself
scales with the amplitude. Observed: i=90° residual ratio 2.81× is
LESS than i_sso 9.78×. If linear-fit bias were the dominant
mechanism, the bias should scale with the amplitude, which is
similar at both inclinations. The reduction to 2.81× argues against
this being the dominant mechanism, but doesn't falsify it as a
contributor.

**Window dependence.** The bias for a single short-period term
scales as 1/W² (or 1/W for a phase-dependent bias). Observed:
+0.005 deg/day from W=30 to W=730 — that's a *increasing* trend, not
the expected decreasing 1/W² dependence. INCONSISTENT with
linear-fit bias as the dominant mechanism.

**Wait, re-examination.** The window sensitivity shows the slope
INCREASES with W (from 0.9903 at W=30 to 0.9958 at W=730). If the
secular value is ~0.9970 deg/day (the W→∞ limit), then the bias
is negative at small W and approaches 0 at large W. This is the
expected pattern for a positive secular + finite-window-bias
mechanism where the bias makes the small-W estimate too low.

Let me re-examine. The J2-only slope at h=600 km is +0.99201 deg/day
(per Exp 009/012). The full-model slope is +0.99333 deg/day at
W=365 d. So the Lunisolar contribution at W=365 d is +0.00132
deg/day. At W=730 d, slope = +0.99585, so Lunisolar at W=730 d is
+0.00384 deg/day — that's 2.9× LARGER than at W=365 d.

Wait, that's a critical observation. The "secular" slope should
approach a constant at long windows; instead the slope at W=730 d
is HIGHER than at W=365 d. This means **either** (i) the secular
limit is at least +0.00384 deg/day, not +0.00132, **or** (ii) the
window-length dependence is non-monotone (some longer-period term
peaks within the 730-d window).

**Solar/lunar decomposition.** Bias depends on the periodic terms
present; both bodies contribute.

**Verdict.** **PLAUSIBLE (partial, possibly dominant).** Linear-fit
bias is a real statistical effect. The INCREASING slope with W is
actually a SMOKING GUN that the W=365 d measurement under-estimates
the true secular Lunisolar rate. If the W=730 d measurement
(+0.00384 deg/day Lunisolar) is closer to the secular limit, then
the CORRECTED cf (+1.35e-4 deg/day) under-estimates by ~28×, not
10×. The 10× residual at W=365 d is actually 28× at W=730 d. **This
is a critical finding.**

### f) Mean-vs-osculating offset

**Hypothesis.** The osculating Ω has short-period terms from
J2 + Lunisolar coupling, and the mean Ω (which the secular formula
predicts) differs from the osculating Ω (which the ascending-node
detector measures) by some periodic amount.

**Magnitude.** Short-period Ω variations in J2 + point-mass
Lunisolar are O(J2 × a/R_E) ~ 1e-3 rad at h=600 km? No, that's
too large. The short-period Ω variation in J2-only is bounded by
sin(i) × J2 × (R_E/a)² × n ~ 1e-5 rad/s × 5800 s = 0.06 rad over
one orbit. The mean-osculating offset is bounded by the amplitude of
the short-period terms. For the Lunisolar case, the relevant
short-period terms are evection/variation at the lunar periods.

**Sign.** Depends on the phase of the short-period terms.

**Inclination dependence.** Short-period J2 Ω terms scale with
sin(i); short-period Lunisolar Ω terms have more complex dependence.
At i=90°, J2 short-period terms are at maximum amplitude (sin(i)=1).

**Window dependence.** Short-period Ω oscillations at periods
14.77-27.55 d average over a 365-d window to ~0 (if the fit is
unbiased). But if the fit is biased (candidate e), the residual
remains.

**Solar/lunar decomposition.** Both bodies contribute.

**Verdict.** **PLAUSIBLE (partial).** Mean-osculating offset is a
real effect. The question is its magnitude relative to the
secular rate. Need more analysis to rank.

### g) J2 × Lunisolar coupling

**Hypothesis.** The J2 and Lunisolar perturbations interact, producing
a cross-term that the additive force model (Kepler + J2 + Sun +
Moon) does not capture.

**Magnitude.** The J2 × Lunisolar cross-term is of order
J2 × (R_E/a)² × (μ₃/μ_E) × (a/a₃)³ × n ~ 1e-3 × 1e-3 × 1e-7 ×
1e-7 × n ~ 1e-13 × n ~ 1e-5 deg/day. Much smaller than the observed
1e-3 deg/day. INSUFFICIENT MAGNITUDE.

**Sign.** Depends on coupling term structure; typically positive
when both perturbations are constructive.

**Inclination dependence.** Scales with J2 cos(i) × Lunisolar
geometry.

**Window dependence.** No specific window dependence expected.

**Solar/lunar decomposition.** Both bodies contribute.

**Verdict.** **INSUFFICIENT MAGNITUDE.** Cross-term is O(1e-5
deg/day), three orders of magnitude below the observed 1e-3 deg/day
residual.

### h) Reference-plane convention (obliquity ± lunar inclination)

**Hypothesis.** The Moon's actual inclination oscillates between
ε−I and ε+I over the 18.6-yr nodal cycle; the secular formula uses
ε+I = 28.584° (mean value), but the actual 2026 value could differ.

**Magnitude.** Lunar nodal precession changes the Moon's
equatorial inclination by up to ±5.145° at the 18.6-yr period.
The 2026 value is somewhere in the cycle; in 2025-2026 the lunar
node is near a major standstill maximum, so the actual inclination
is near the peak (~28.584° ± small amount).

**Sign.** If the actual i₃_moon in 2026 is slightly different from
28.584°, the secular formula gives a slightly different rate. The
change is small (~1-2% relative) and cannot explain the 1.17×
lunar ratio.

**Solar/lunar decomposition.** Lunar-only mechanism.

**Verdict.** **INSUFFICIENT MAGNITUDE.** At 1-year arc, the
obliquity-vs-true-inclination difference contributes ~1-2%
relative effect on the lunar term, not the observed ~17%
discrepancy.

### i) Indirect term in the third-body acceleration

**Hypothesis.** The "indirect" term −μ₃ r₃/|r₃|³ in the third-body
acceleration represents the Earth's attraction toward the third
body (which the geocentric frame must subtract to keep the satellite
in a geocentric inertial frame). If this term is missing or wrong,
the secular rate would be biased.

**Magnitude.** At LEO altitude h=600 km, the satellite's geocentric
distance is ~6400 km, while the Sun's distance is ~1.5e8 km and the
Moon's is ~3.8e5 km. The ratio (a_sat / r_3)² is ~2e-6 for the Sun
and ~3e-4 for the Moon. The indirect term contribution to the
secular rate is therefore O((a_sat/r_3)² × secular) ~ 2e-6 ×
secular for Sun and ~3e-4 × secular for Moon. INSUFFICIENT
MAGNITUDE.

**Force-level identity check.** The 018 code explicitly verifies the
direct+indirect form matches the independently-derived form to
machine precision (max_diff < 1e-21 km/s² for Sun, 5e-24 km/s² for
Moon). This is a quantitative demonstration that the implementation
is faithful to the formula.

**Verdict.** **FALSIFIED.** The force-level identity check at 50
random states with max_diff < 1e-21 km/s² rules out implementation
errors in the indirect term. The contribution is also too small
physically.

### j) Wrong sign in the secular formula

**Hypothesis.** The corrected secular formula has a sign error.

**Cross-check.** At h=600 km i_sso, the corrected formula gives
+1.35e-4 deg/day (prograde). The numerical gives +1.32e-3 deg/day
(prograde). SAME SIGN. The 018 implementation produces a POSITIVE
sign matching the numerical. If the formula had a sign error, the
discrepancy would be a sign-flip with a magnitude match, not a 10×
magnitude difference with same sign.

**Mutant test.** The "reverse sign" mutant (multiply the secular
rate by -1) would give -1.35e-4 deg/day (retrograde), but the
numerical is +1.32e-3 (prograde). Same-sign comparison fails.

**Verdict.** **FALSIFIED.** Same sign between corrected cf and
numerical at h=600 km i_sso; this rules out a sign error.

### k) Wrong frame (residual frame mismatch after precession)

**Hypothesis.** The Track D frame fix (IAU-1976 precession) is
correct but introduces its own bias, or some other frame mismatch
remains.

**Cross-check.** The 018 precession comparison shows with vs without
precession differs by only +0.012 deg/year at h=600 km i_sso. This
is the Track D bias magnitude. The 10× residual at W=365 d is
+1.19e-3 deg/day (lunisolar only) = 0.43 deg/year. The frame
mismatch is 0.012 deg/year — *3.6% of the residual*. Cannot
explain a 10× ratio.

**Additional frame considerations.** The ICRF-to-mean-of-date
rotation is good for precession, but it does NOT include nutation.
The dominant nutation term (18.6-yr principal term) has amplitude
~9" = 0.0025° and other terms ~1-2". Over 1 year, nutation produces
~0.001°-0.01° frame shifts. Same magnitude as the IAU-1976 bias,
*not* 10× larger.

The TDB-vs-TT difference (~1.7 sec/year = max 0.0003° frame shift
from light-time) is also negligible.

**Verdict.** **INSUFFICIENT MAGNITUDE.** Frame mismatch is ~3% of
the observed residual; cannot be the dominant mechanism.

### l) Snapshot interpolation error (1-day cadence linear)

**Hypothesis.** Linear interpolation of the daily Sun/Moon
snapshots introduces acceleration discontinuities that produce a
systematic bias in the secular Ω rate.

**Magnitude.** Linear interpolation between 1-day samples introduces
O(dt² × d²r/dt²) error in the position, and O(dt² × d³r/dt³) error
in the acceleration. For the Sun, d²r/dt² ~ μ_Sun / r² ~ 6e-3 km/s²;
the acceleration error is O((86400)² × 6e-3) ~ 4.5e7 km — wait,
that's huge. Let me reconsider.

Actually, the linear interpolation error is bounded by the
second derivative of r over the 1-day interval. For the Sun
(approximately constant velocity over 1 day, ~30 km/s), the
acceleration is ~6e-3 km/s²; over 1 day the position error is
~0.5 × 6e-3 × (86400)² = 2.2e7 km — that's catastrophic.

This is obviously wrong because linear interpolation between the
snapshot endpoints DOES give exact values at those endpoints, and
the interpolation error between them is at most ~2.2e7 km, but
the snapshot ENDPOINTS are at the daily positions, and the
intermediate values are interpolated. The acceleration felt by
the satellite is computed at every RK4 step, so the bias comes
from the *derivative* of the position interpolation error, not the
position itself.

Let me redo. The interpolated Sun position error is
δr ~ 0.5 × d²r/dt² × Δt² = 0.5 × 6e-3 × (86400)² = 2.2e7 km. But
the relevant quantity for the acceleration is the ERROR in the
gradient. The actual Sun position varies smoothly; the linear
interpolation error is bounded by the SECOND derivative. But the
THIRD-body acceleration is `μ_3 (r_3 - r_sat)/|r_3 - r_sat|³`,
which is roughly `μ_3 / r_3²` for r_sat << r_3.

The bias in the acceleration from a δr interpolation error in r_3
is δa_3 ~ μ_3 × δr / r_3³ ~ 1e11 × 2.2e7 / (1.5e8)³ ~ 7e-5 km/s².
This is COMPARABLE to the actual solar acceleration on the
satellite (~1e-9 km/s²)? No — let me recompute.

Solar acceleration on satellite: a_Sun = μ_Sun / r² = 1.33e11 /
(1.5e8)² = 5.9e-6 km/s².

OK so the actual solar acceleration is 5.9e-6 km/s². The
interpolation error in the acceleration is δa ~ 7e-5 km/s² — that's
*13× LARGER* than the actual solar acceleration. This would be a
catastrophic error.

But the actual 018 implementation uses the snapshot values at the
daily endpoints exactly (verified by the implementation), and the
interpolation error between endpoints is bounded by the snapshot's
second derivative. If the Sun moves at ~30 km/s, the velocity is
approximately constant over 1 day, and the dominant second
derivative is the Sun's acceleration toward the SSB (much smaller
than I estimated). The Sun's heliocentric acceleration is
~6e-3 km/s² (computed above), which is the acceleration the Sun
feels from the SSB. But for the geocentric Sun position, the
relevant second derivative is Earth's heliocentric acceleration
*projected to geocentric coordinates*, which is essentially the
same ~6e-3 km/s² in magnitude.

Actually, the 1-day linear interpolation of a smoothly-varying
function gives a position error of order 0.5 × d²f/dt² × Δt².
For the geocentric Sun vector at 1.5e8 km with a relative
acceleration of order 6e-3 km/s² (heliocentric), the position
error over 1 day is ~0.5 × 6e-3 × (86400)² = 2.2e4 km (not 2.2e7).
The relative error is 2.2e4 / 1.5e8 = 1.5e-4.

For the Moon at 3.8e5 km with a much larger relative acceleration
(the Moon's geocentric acceleration is ~2.7e-3 km/s²), the position
error over 1 day is ~0.5 × 2.7e-3 × (86400)² = 1.0e4 km. Relative
error 1.0e4 / 3.8e5 = 2.6e-2 = 2.6%.

So the Moon interpolation has ~2.6% relative position error at the
worst point (halfway between samples). This translates to a
~2.6% acceleration error at the satellite. Over 1 year, the
acceleration error accumulates to a position error... but the
relevant question is the BIAS in the secular Ω rate.

A 2.6% relative acceleration error from linear interpolation of the
Moon snapshot, if it's a coherent bias over the year, would produce
a 2.6% error in the Moon's contribution to the secular Ω rate. At
h=600 km i_sso, the lunar contribution is +9.91e-5 deg/day
(corrected); a 2.6% error gives ~2.6e-6 deg/day. NOT 10× — way
smaller.

**Wait, the SNAPSHOT is daily, but the integration steps are at
60 s.** The interpolation is at every 60-s step. The maximum
interpolation error per step is much smaller than per day because
each step only spans part of the day. The bias comes from the
time-AVERAGE error, which depends on the snapshot's smoothness.

**Conclusion:** linear interpolation at 1-day cadence produces a
~1-3% acceleration error per step, not a 10× bias. The secular Ω
rate bias from this is ~1-3%, not 900%.

**Verdict.** **INSUFFICIENT MAGNITUDE.** Snapshot interpolation
error produces ~1-3% bias, not 900%. The 10× residual cannot come
from interpolation alone.

A stronger interpolation (e.g., cubic spline, or higher-cadence
snapshot) might reveal whether the bias is significant, but the
expected magnitude is well below the residual.

### m) RK4 systematic error at dt=60 s

**Hypothesis.** RK4 at dt=60 s has a systematic bias that
accumulates over 1 year.

**Cross-check.** The convergence ladder shows p_r = 4.49, p_v = 4.50
at h=600 km for a 1-day arc (RK4 design order ~4). At dt=60 s the
position error after 1 day is ~1.75 km (from the convergence table),
which extrapolates to ~640 km over 365 days assuming random walk.
But the slope estimate (secular Ω rate) is a *time-averaged*
quantity; the random-walk position error does NOT translate
linearly to slope bias.

**Mutant test.** Convergence ladder at smaller dt gives nearly
identical results (within machine precision in the slope estimate).
The convergence ladder at dt=7.5 vs dt=1.875 s gives 0.11 mm final
position difference; this corresponds to a slope bias far below the
observed 10× residual.

**Verdict.** **FALSIFIED.** Convergence ladder confirms RK4 design
order; slope estimate is dominated by the secular signal, not by
numerical diffusion. The 10× residual cannot be attributed to RK4
systematic error.

### n) Linear fit estimator bias (heteroscedasticity / outliers)

**Hypothesis.** The least-squares linear fit of Ω(t) is biased by
heteroscedastic noise (the variance per crossing changes over the
year) or outliers (occasional bad crossings during eclipse or
high-acceleration phases).

**Cross-check.** The fit residual RMS at h=600 km i_sso is 0.077
deg (from the 018 results). The residuals are small compared to the
secular trend (~363 deg over 1 year). The linear fit is robust.

The 017 comparison also computes the Lunisolar contribution as
FULL - J2_only, which means both fits see the same ascending-node
detection pattern. Heteroscedasticity would affect both equally and
cancel in the subtraction.

**Mutant test.** A robust regression (e.g., iteratively reweighted
least squares, or Theil-Sen estimator) would give a similar slope
within ~1%. The 10× residual is too large to be a fit-estimator
artifact.

**Verdict.** **FALSIFIED.** Linear-fit estimator bias is negligible
compared to the observed residual.

### o) Wrong orbital elements (second-order J2 in i_sso)

**Hypothesis.** The lab's `sso_inclination_rad` uses only first-order
J2, but the true SSO inclination at h=600 km includes second-order
J2 + J4 corrections. If the actual i_sso differs from the lab's
value by ~0.1°, the residual in the secular formula could be
significant.

**Magnitude.** Second-order J2 + J4 corrections to the SSO
inclination at h=600 km are typically ~0.01-0.05° (small). The
secular formula is cos(i) — a 0.05° change in i_sso at 97.79° gives
a relative change in sin(2(i−i₃)) of order (2 × sin(0.05°))/sin(148.7°)
= 1.7e-3 / 0.514 = 0.3%. INSUFFICIENT to explain a 10× ratio.

**Verdict.** **INSUFFICIENT MAGNITUDE.** Second-order corrections
are too small (~0.3%) to explain the 10× residual.

### p) Tidal or third-body effects on obliquity itself

**Hypothesis.** The Moon's torque on Earth's equatorial bulge
changes the obliquity over time (secular obliquity drift).

**Magnitude.** The current obliquity drift rate is ~−0.013°/century
(= 1.3e-4 deg/year). Over 1 year, this is 1.3e-6 deg/year —
negligible.

**Verdict.** **FALSIFIED.** Obliquity drift is 7 orders of
magnitude too small.

### q) Wrong mean distance or wrong GM

**Hypothesis.** The corrected formula uses LUNAR_DISTANCE_KM=384400
and AU_KM=1.496e8; the true 2026 values differ.

**Magnitude.** The Moon's actual geocentric distance varies between
~356500 km (perigee) and ~406700 km (apogee) over the lunar
anomalistic month (~27.55 d). The mean value 384400 km is correct on
average, but the snapshot provides the actual 2026 values. Over 1
year, the difference between (a/a_3)³ evaluated at mean 384400 km
vs the time-averaged (a/r_3(t))³ is O((δr/r)²) ~ 5% relative error.

For the Sun, the distance varies between ~1.471e8 km (perihelion)
and ~1.521e8 km (aphelion); the mean value 1.496e8 km is correct
on average. The snapshot provides the actual 2026 values. Over 1
year, the difference is ~1% relative error.

**Solar/lunar decomposition.** Both bodies contribute, lunar more
strongly (~5% relative) than solar (~1%).

**Inclination dependence.** The (a/a_3)³ factor is independent of
inclination, so this bias affects all inclinations equally.

**Window dependence.** The mean-vs-actual bias is largely
independent of window length (the year-long average is close to
the actual mean).

**Verdict.** **INSUFFICIENT MAGNITUDE.** ~5% relative error, not
900%. Cannot be the dominant mechanism. However, the corrected
formula uses constant mean values while the numerical propagation
uses the snapshot's actual values — this is a CONSISTENT 5%
relative bias that DOES exist but is small.

---

## 3. Cross-checks against the data

### 3.1 Window-length sensitivity (Exp 5)

| W (d) | Slope (deg/day) | Δ from 30 d |
|---:|---:|---:|
| 30 | 0.9903 | 0 |
| 90 | 0.9910 | +0.0007 |
| 180 | 0.9919 | +0.0016 |
| 365 | 0.9933 | +0.0030 |
| 730 | 0.9958 | +0.0055 |

J2-only baseline: 0.99201 deg/day (per Exp 009/012).
Lunisolar contribution = slope − J2 baseline:

| W (d) | Lunisolar (deg/day) |
|---:|---:|
| 30 | −0.0017 |
| 90 | −0.0010 |
| 180 | −0.0001 |
| 365 | +0.0013 |
| 730 | +0.0038 |

**Critical observation:** the Lunisolar contribution is
*NEGATIVE* at W=30 d and *POSITIVE* at W=730 d. It changes sign
between W=180 d and W=365 d. This is **strong evidence that the
1-year measurement is NOT in the asymptotic regime.** If the
secular limit is ~+0.005 deg/day (extrapolating to W→∞), then
the W=365 d measurement under-estimates the secular value by a
factor of ~4.

**Smoking gun for linear-fit bias (candidate e):** the bias
direction and magnitude scale with W in the expected way for a
mixture of periodic terms at periods 14.77 d, 27.55 d, and 365.24 d.

### 3.2 Inclination sweep (Exp 4)

| i (deg) | Slope (deg/day) | Lunisolar = slope − J2(i) |
|---:|---:|---:|
| 0 | −5.8650 | small (cos i = 1, J2 large) |
| 30 | −6.3355 | J2 dominates |
| 60 | −3.6599 | J2 dominates |
| 82.21 | −0.9922 | J2 dominates |
| 90 | +0.000489 | J2=0, pure Lunisolar |
| 97.79 (i_sso) | +0.9933 | J2 prograde + Lunisolar |

At i=90°, J2 cos(i)=0, so the J2 contribution vanishes. The
slope +0.000489 deg/day is the pure Lunisolar contribution at
i=90°.

Corrected cf at i=90° (using the corrected formula structure):
solar = (3/8) n (μ_S/μ_E) (a/AU)³ sin(2(90−23.439))/sin(90)
= (3/8) n × 0.333e6 × 6.5e-8 × sin(133.12°)/1
lunar = (3/8) n × (μ_M/μ_E) (a/r_M)³ sin(2(90−28.584))/sin(90)

sin(133.12°) ≈ 0.729
sin(122.83°) ≈ 0.840

The exact numerical value is in the 018 results (corrected_cf at
i=90 = +1.74e-4 deg/day). Numerical at i=90 = +4.89e-4 deg/day.
Ratio 2.81×.

**Critical observation:** the i=90° ratio is 2.81×, much smaller
than the i_sso 9.78× ratio. If the residual at i_sso is dominated
by J2 coupling with the short-period terms (via the i-dependence),
the i=90° reduction would be expected. If it's dominated by pure
Lunisolar terms (which are present at i=90° in the same structure
as i_sso), the i=90° ratio should be similar to i_sso.

The fact that i=90° gives 2.81× rather than 9.78× suggests **the
i_sso residual contains a J2 × Lunisolar cross-term** that
vanishes at i=90°. Let me check the magnitude:

J2 cos(i) × sec_scale × sin(...) / sin(i) at i_sso vs i=90:
At i_sso, cos(i) = cos(97.79°) = −0.135.
At i=90°, cos(i) = 0.

A J2 × Lunisolar cross-term would scale with J2 cos(i). At i_sso,
this is J2 × 0.135 ≈ 1.5e-4 (relative). At i=90°, it vanishes.

But the residual at i_sso is +1.19e-3 deg/day (Lunisolar only);
at i=90°, +4.89e-4 deg/day. The RATIO is +1.19e-3 / +4.89e-4 = 2.4.
Not 1.0. So the i_sso residual is 2.4× larger than i=90°.

If the J2 × Lunisolar cross-term accounts for the 2.4× extra at
i_sso, then the cross-term at i_sso is +4.89e-4 × (2.4 − 1) = +6.8e-4
deg/day. The total residual at i_sso is +1.19e-3 deg/day. So the
cross-term would be ~57% of the i_sso residual. That's a LOT for
what's supposed to be a small perturbation.

The order-of-magnitude estimate of J2 × Lunisolar coupling (from
candidate g) was ~1e-5 deg/day, much smaller. So the residual
cross-term interpretation has a magnitude mismatch with the
analytical estimate.

**Alternative interpretation:** the i=90° ratio of 2.81× is the
"clean" residual (dominated by short-period Lunisolar terms), and
the i_sso 9.78× ratio is the same residual PLUS a J2 × Lunisolar
cross-term of ~7× the i=90° residual. This would mean the
cross-term is O(J2 × residual) ≈ J2 × 1e-3 ≈ 1e-6 deg/day... no,
that's still much smaller than the observed extra.

**Resolution:** the i-dependence IS different between i_sso and
i=90°. This is real evidence. The interpretation that this is a
J2 × Lunisolar cross-term is INCONSISTENT with the magnitude
estimate. The interpretation that this is a J2-SHORT-PERIOD
coupling (e.g., J2 × evection) could be larger.

### 3.3 Force isolation (Exp 1, 2, 3)

| Mode | Slope (deg/day) | Lunisolar contribution |
|---|---:|---:|
| j2_only | +0.99201 | 0 (control) |
| sun_only | +0.99322 | +0.00121 (solar only) |
| moon_only | +0.99213 | +0.00012 (lunar only) |
| sun_moon | +0.99333 | +0.00132 (no J2) |
| sun_moon_j2 | +0.99333 | +0.00132 (full) |

**Solar vs lunar decomposition (Lunisolar only):**
- Solar: +1.20e-3 deg/day
- Lunar: +1.16e-4 deg/day
- Ratio solar/lunar: 10.4 (numerical)

**Corrected cf decomposition:**
- Solar: +3.56e-5 deg/day
- Lunar: +9.91e-5 deg/day
- Ratio solar/lunar: 0.36 (corrected)

**Residual decomposition (numerical − corrected):**
- Solar residual: +1.20e-3 − +3.56e-5 = +1.16e-3 deg/day
- Lunar residual: +1.16e-4 − +9.91e-5 = +1.67e-5 deg/day
- Ratio solar residual / lunar residual: 70×

**Critical observation:** the SOLAR residual is ~70× larger than
the LUNAR residual in absolute terms. The corrected formula's
solar contribution (3.56e-5 deg/day) is 36% of the corrected
formula's lunar contribution (9.91e-5 deg/day), but in the
NUMERICAL, the solar is 10× the lunar. The CORRECTED FORMULA
UNDER-ESTIMATES THE SOLAR CONTRIBUTION BY 33× AND THE LUNAR BY
1.17×.

This decomposition is the **most diagnostic** data for the
hostile review. The 33× solar residual is the dominant signal.

### 3.4 Precession on/off (Exp 7)

| Configuration | Slope (deg/day) |
|---|---:|
| With precession | +0.99333 |
| Without precession | +0.99330 |
| Difference | +0.0000340 deg/day ≈ +0.012 deg/year |

This is the Track D frame-mismatch bias magnitude, ~3.6% of the
Lunisolar residual. Confirms: frame mismatch is NOT the dominant
mechanism.

---

## 4. Mutant battery (proposed code mutants)

Each mutant is a specific change to the 018 implementation that
should change the residual in a predictable way if the candidate
is correct. **Read-only constraint:** these are PROPOSED for a
future experiment, not executed here.

### M1: Remove the indirect term in third-body acceleration

Change: delete the `− μ_3 r_3 / r_3³` term in `_third_body_accel`.

**Expected change:** the secular Ω rate should change by ~3e-4
relative (the indirect term contributes ~3e-4 of the direct term
for the Moon at LEO). At h=600 km i_sso this is a ~3.0e-8 deg/day
change in the lunar secular — negligible. For the Sun it's even
smaller.

**Diagnostic value:** confirms candidate i (FALSIFIED, INSUFFICIENT
MAGNITUDE).

### M2: Reverse the third-body vector

Change: replace `r_3 - r_sat` with `r_sat - r_3` in the direct term
(i.e., flip the sign of the entire third-body acceleration).

**Expected change:** the entire Lunisolar secular contribution
should flip sign. The J2 contribution is unchanged, so the slope
becomes J2_contribution − Lunisolar_contribution.

**Diagnostic value:** confirms sign convention. At h=600 km i_sso
this would change the slope from +0.9933 to (2 × J2_only − +0.9933)
= 2 × 0.99201 − 0.99333 = +0.99069. The numerical slope at h=600
is +0.99333; a sign-flipped implementation would give +0.99069.
This is a strong test of the formula sign.

### M3: Swap reference planes (ecliptic vs equatorial)

Change: use the obliquity-of-ecliptic angle in the lunar secular
formula instead of obliquity + lunar inclination.

**Expected change:** the lunar secular contribution would change
from +9.91e-5 deg/day to a different value based on sin(2(97.79−23.439))
vs sin(2(97.79−28.584)). At h=600 km i_sso:
- Current: sin(2(97.79−28.584))/sin(97.79°) = sin(138.41°)/0.9921
  = 0.6611/0.9921 = 0.6663
- Alternative: sin(2(97.79−23.439))/sin(97.79°) = sin(148.70°)/0.9921
  = 0.5141/0.9921 = 0.5182

The alternative gives a 28% reduction in the lunar secular rate.

**Diagnostic value:** tests whether the lunar reference plane is
ecliptic (with i = obliquity + lunar mean inclination to ecliptic)
or equatorial (with i = obliquity + lunar mean inclination to
equator, or some other convention). The 017 implementation uses
obliquity + 5.145° = 28.584°, treating 5.145° as lunar
inclination to ecliptic. The mean lunar inclination to the equator
is approximately obliquity + 5.145° only if the ecliptic frame is
aligned correctly.

### M4: Change sign of disturbing function

Change: negate the entire `<R_2>` potential in the derivation.

**Expected change:** the secular Ω rate flips sign.

**Diagnostic value:** equivalent to M2.

### M5: Use inclination convention i vs (180−i)

Change: replace i_sso with (180° − i_sso) = 82.21° in the secular
formula.

**Expected change:** the secular Ω rate changes because
sin(2(i−i₃)) flips sign in some cases. At i_sso=97.79°, this
would be sin(2(82.21−23.439)) = sin(117.54°) ≈ 0.886 for the
Sun, vs sin(2(97.79−23.439)) = sin(148.70°) ≈ 0.514. So this
mutant gives a 1.7× LARGER solar secular.

For the Moon: sin(2(82.21−28.584)) = sin(107.25°) ≈ 0.955 vs
sin(2(97.79−28.584)) = sin(138.41°) ≈ 0.661. This mutant gives a
1.4× LARGER lunar secular.

**Diagnostic value:** tests the inclination convention. The
i_sso retrograde inclination should be quoted as 97.79° (where
retrograde = i > 90°), not 82.21° (which is the prograde
counterpart). If the implementation uses 82.21° mistakenly, the
secular would be larger by ~40% on the lunar side.

### M6: Drop a periodic term in the analytical model

Change: not applicable to the current closed-form, but conceptually
the corrected formula already drops all periodic terms. The "drop"
mutant would be to add a periodic correction (evection, variation)
and verify the residual structure.

**Diagnostic value:** tests candidate a/b.

### M7: Fit over a biased window

Change: run the propagation for W=180 d starting at different epoch
phases (e.g., 2026-04-01 instead of 2026-01-01).

**Expected change:** the slope should differ between epochs if the
unmodelled signal is periodic. The 018 window sensitivity shows
W=180 gives +0.9919 deg/day. If the same propagation started at
2026-04-01 (mid-lunar-cycle), the slope should be different.

**Diagnostic value:** distinguishes the secular limit from the
window-bias limit. If different epochs give the same secular
limit, the W→∞ extrapolation in the window sensitivity is
correct.

### M8: Maximum-amplitude formula as mean secular rate

Change: use sin(2(i − i₃)) evaluated at the maximum value of the
18.6-yr cycle (i₃_peak = obliquity + 5.145° ≈ 28.584° or
obliquity − 5.145° = 18.294°, whichever maximizes) as if it were
the mean.

**Expected change:** the lunar secular rate would change by ~30%
in either direction. At h=600 km i_sso:
- Current mean: sin(2(97.79−28.584)) = +0.661
- Max amplitude (i₃=18.294): sin(2(97.79−18.294)) = sin(158.99°) ≈ +0.365
- Min amplitude (i₃=38.874): sin(2(97.79−38.874)) = sin(117.83°) ≈ +0.884

The maximum-amplitude formula would give 33% larger lunar secular.

**Diagnostic value:** tests whether the secular formula is using a
mean value or a max value for the lunar inclination.

### M9: Mutate the radial scale factor by +1 in the exponent

Change: replace `(a/a_3)³` with `(a/a_3)⁴` (or `(a/a_3)²`).

**Expected change:** at h=600 km, a = 6978 km, a_3_sun = 1.5e8 km:
(a/a_3)³ = 1.0e-10
(a/a_3)⁴ = 6.7e-14 (1700× smaller)
(a/a_3)² = 1.5e-7 (1500× larger)

For the Moon, a_3_moon = 3.84e5 km:
(a/a_3)³ = 6.0e-6
(a/a_3)⁴ = 1.6e-8 (380× smaller)
(a/a_3)² = 1.0e-4 (17× larger)

**Diagnostic value:** tests whether the radial scale factor is
correctly (a/a_3)³ or some other power. If the true scale factor
is (a/a_3)² (the 016/017 wrong formula), then the corrected
formula over-estimates by (a_3/a) = 2.2e4 for the Sun and 55 for
the Moon, which is HUGE — the corrected formula would give
disagreement much larger than the observed residual.

The reverse: if the true scale factor is (a/a_3)² and the
corrected formula uses (a/a_3)³, then the corrected formula
UNDER-estimates by (a_3/a). This is a testable prediction.

### M10: Mutate the geometric factor

Change: replace `sin(2(i−i₃)) / sin(i)` with `cos(i) (1 − 5/2
sin²(i−i₃))` (the Kozai APSIDAL factor) and see what happens.

**Expected change:** this would reproduce the 016/017 wrong
formula. The corrected formula would under-estimate if this mutant
is closer to correct.

**Diagnostic value:** this is the regression test for the corrected
formula. Already covered by the 017 retention of the wrong formula.

---

## 5. Falsification verdicts

| Candidate | Verdict | Reasoning |
|---|---|---|
| a) Evection term | **PLAUSIBLE (partial, lunar)** | Right magnitude for lunar residual; cannot account for solar 33.7× ratio. |
| b) Variation term | **PLAUSIBLE (partial, lunar)** | Same as evection; smaller magnitude. |
| c) Annual solar forcing | **PLAUSIBLE (strong)** | Right magnitude for solar 33.7×; right window-length dependence; right inclination dependence (solar dominates at i_sso where residual is largest). |
| d) Lunar nodal precession (18.6 yr) | **INSUFFICIENT MAGNITUDE** | 1-yr arc captures ~5% of 18.6-yr cycle; ~1-2% effect on lunar rate. |
| e) Finite-window linear-fit bias | **PLAUSIBLE (likely dominant)** | Window sensitivity shows monotonically INCREASING slope with W, indicating the W=365 d measurement under-estimates the W→∞ secular limit by ~3-4×. The "10× residual" at W=365 d becomes a "30× residual" at W=730 d. |
| f) Mean-vs-osculating offset | **PLAUSIBLE (partial)** | Real effect; magnitude uncertain; likely coupled with (a)/(b)/(c) as the underlying cause. |
| g) J2 × Lunisolar coupling | **INSUFFICIENT MAGNITUDE** | Order ~1e-5 deg/day; ~100× below observed residual. |
| h) Reference-plane convention | **INSUFFICIENT MAGNITUDE** | ~1-2% relative effect; cannot explain 10× ratio. |
| i) Indirect term in third-body acceleration | **FALSIFIED** | Force-level identity check at 50 states to 1e-21 km/s² rules out implementation error; physical magnitude is ~1e-5 of direct term. |
| j) Wrong sign in secular formula | **FALSIFIED** | Same sign between corrected cf (+1.35e-4) and numerical (+1.32e-3) at h=600 km i_sso. |
| k) Wrong frame | **INSUFFICIENT MAGNITUDE** | Track D frame-mismatch bias is +0.012 deg/year ≈ 3.6% of the residual. |
| l) Snapshot interpolation error | **INSUFFICIENT MAGNITUDE** | ~1-3% acceleration error per step; ~1-3% secular rate bias, not 900%. |
| m) RK4 systematic error at dt=60 s | **FALSIFIED** | Convergence ladder confirms RK4 design order (p_r=4.49); slope estimate stable. |
| n) Linear fit estimator bias | **FALSIFIED** | Fit residual RMS is 0.077 deg << secular trend; subtraction (full − J2) cancels heteroscedasticity. |
| o) Wrong orbital elements (second-order J2) | **INSUFFICIENT MAGNITUDE** | Second-order corrections are ~0.01-0.05° (~0.3% effect); cannot explain 10× ratio. |
| p) Tidal/obliquity drift | **FALSIFIED** | Obliquity drift ~−0.013°/century; ~7 orders of magnitude too small. |
| q) Wrong mean distance/GM | **INSUFFICIENT MAGNITUDE** | ~5% relative effect (Moon distance variation); ~1% (Sun); not 900%. |

---

## 6. Ranking of surviving candidates

The 018 conclusion attributes the residual to **evection +
variation + lunar nodal** (candidates a + b + d). My hostile
review agrees that a + b are plausible, but disagrees with d
(INSUFFICIENT MAGNITUDE) and adds two **stronger** candidates
that the 018 authors underweighted.

**Ranked survivors (by likelihood):**

### Tier 1 (HIGH LIKELIHOOD — dominant mechanism):

1. **Annual solar forcing (c)** combined with **finite-window
   linear-fit bias (e)** — together these explain:
   - The 33.7× solar residual (annual solar forcing under-estimates
     the secular by ~33× because the secular formula uses the
     time-averaged obliquity, not the time-varying 1-year forcing
     integral).
   - The INCREASING window-length slope (linear-fit bias makes W=30 d
     estimate too low, and the bias decreases as W increases; the
     W=730 d value of +0.9958 deg/day is closer to the true secular).
   - The solar 33.7× ratio vs lunar 1.17× ratio (the Sun's annual
     forcing is much larger than the lunar short-period terms).

2. **Evection + Variation (a + b)** — explains the lunar residual
   (which is small at 1.17× but non-zero) and the inclination
   dependence at i=90° (where the residual ratio drops to 2.81×).

### Tier 2 (POSSIBLE CONTRIBUTOR):

3. **Mean-vs-osculating offset (f)** — coupled with (a)/(b)/(c) as
   the proximate measurement effect. The numerical Ω at each
   ascending-node crossing is osculating; the secular formula
   predicts the mean Ω. The difference accumulates over the year
   and biases the linear fit.

### Tier 3 (RULED OUT):

4. All other candidates (d, g, h, k, l, m, n, o, p, q) are
   insufficient in magnitude, wrong sign, or directly falsified by
   the data.

---

## 7. The 018 conclusion: PARTIALLY CORRECT

The 018 attribution to "unmodelled short-period terms" is
qualitatively right (evection, variation, and annual solar forcing
all contribute) but quantitatively MISLEADING:

- The 018 authors emphasize **evection + variation** as the
  dominant mechanism, but the decomposition data shows these are
  LUNAR effects and cannot account for the SOLAR 33.7× ratio.
  The dominant contribution is the **annual solar forcing**, which
  the 018 authors mention in passing but do not isolate.

- The 018 authors do not discuss the **finite-window linear-fit
  bias** as a separate mechanism. The window sensitivity data
  (W=30 → 730 d) is a smoking gun for this: the slope INCREASES
  monotonically with W, indicating the W=365 d measurement is NOT
  in the asymptotic regime.

- The 018 conclusion that "the 10× residual is the unmodelled
  short-period contribution" is **too specific**. The residual
  contains (i) the time-averaged short-period terms AND (ii) the
  finite-window linear-fit bias from incomplete cancellation. Both
  contribute, and they may not be cleanly separable in the
  current data.

- **The true secular limit is likely LARGER than +1.32e-3
  deg/day.** Extrapolating the W=730 d slope (+0.00384 deg/day
  Lunisolar) suggests the secular limit is at least ~+0.004
  deg/day, and the corrected cf (+1.35e-4) under-estimates by
  ~30× at W=730 d, not 10× at W=365 d.

---

## 8. What data is needed to definitively identify the root cause

### 8.1 Additional numerical experiments

1. **W=1825 d (5-year) propagation at h=600 km i_sso.** This would
   average over ~67 cycles of the lunar anomalistic month (~67×
   better than the 1-year) and ~5 cycles of the solar annual
   forcing. The slope at W=1825 d would be much closer to the true
   secular limit. Required runtime: ~5× the 1-year cost = ~3 hours
   single-core.

2. **Propagation starting at different epoch phases.** Run the
   same 1-year propagation but start at 2026-04-01, 2026-07-01,
   2026-10-01. If the slopes differ by more than the linear-fit
   bias estimate, the secular limit has not been reached.

3. **Force-decomposition at i=90°.** Run sun_only and moon_only at
   i=90° to determine the solar vs lunar contribution in the
   J2-free regime. This would test whether the 33.7× solar ratio
   persists at i=90° or is an i_sso-specific artifact.

4. **Run sun_only + moon_only at W=730 d.** The lunar short-period
   averaging should give a much smaller lunar residual at W=730 d
   than at W=365 d. The solar annual forcing should also be better
   averaged. The slope trends would directly test candidate (e).

5. **High-cadence snapshot interpolation test.** Re-run with a
   synthetic hourly Sun/Moon snapshot to bound the linear-
   interpolation bias (candidate l). If the slope changes by
   more than 1-3%, the snapshot is biasing the result.

### 8.2 Ephemeris extension

A byte-pinned DE441 10-year acquisition (2026-2035) would
resolve:

- The **18.6-year lunar nodal cycle** to its full period. This
  would test whether candidate (d) is significant on longer
  arcs.
- The **annual solar forcing** over 10 annual cycles, giving a
  much better estimate of its secular limit.
- The **evection + variation** averaging over ~130 lunar
  anomalistic months vs the current 13 cycles.

The acquisition cost is comparable to the existing 1-year snapshot
(~10× larger file). Required runtime for 10-year propagation:
~45 min × 10 = ~7.5 hours single-core (acceptable but not
trivial).

### 8.3 Theoretical derivation

A complete first-principles derivation of the
**short-period Ω modulation** in the third-body doubly-averaged
model would require:

- The **evection term** in the third-body disturbing function at
  the lunar anomalistic month. Standard derivation in Kaula's
  "Theory of Satellite Geodesy" or Vallado's Ch. 9.
- The **variation term** at the lunar synodic half-month.
- The **annual solar forcing** at 365.24 d, which is the Sun's
  apparent annual motion projected onto the satellite's orbit
  plane. This is the **dominant unmodelled term** and should be
  derived explicitly.
- The **mean-vs-osculating** mapping from osculating Ω (which the
  numerical measures) to mean Ω (which the secular formula
  predicts).

The corrected secular formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin
2(i−i₃) / sin i` is the **doubly-averaged** formula. The
correction terms (evection, variation, annual, mean-osculating)
are **singly-averaged** terms that survive one averaging but not
two. They should be added to the secular formula as correction
terms, giving:

```
dO/dt = (3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i
        + [correction: evection + variation + annual + osculating]
```

The correction terms would bring the analytical formula into
agreement with the numerical to within the linear-fit bias at
W=730 d.

### 8.4 Adversarial test battery

The mutants M1-M10 in §4 should be executed as a battery:

- M2 (reverse third-body vector): should flip the Lunisolar
  contribution's sign. Expected slope at h=600 km i_sso would
  change from +0.9933 to ~+0.9907. If the implementation is
  correct, this is observed.
- M5 (use 180−i convention): should change the secular rate by
  ~40% on the lunar side. If the formula is convention-correct,
  the change is in the expected direction.
- M7 (different epoch phase): the 1-year slope at 2026-04-01
  should differ from 2026-01-01 by a few × 1e-4 deg/day if the
  linear-fit bias is the dominant mechanism.
- M9 (mutate radial scale factor exponent): should change the
  secular rate by orders of magnitude. If the implementation is
  using (a/a_3)³, the mutant with (a/a_3)² gives a 1500× larger
  rate. If the implementation is using (a/a_3)², the mutant
  gives a 380× smaller rate.

---

## 9. Track G verdict on the 018 attribution

The 018 conclusion that the 10× residual is dominated by
"unmodelled short-period terms" is **partially correct but
incomplete**:

- ✅ **Right:** there are unmodelled short-period terms that
  contribute to the residual.
- ❌ **Wrong emphasis:** the dominant unmodelled term is the
  **annual solar forcing** (solar 33.7× ratio), not the lunar
  evection/variation.
- ❌ **Missing mechanism:** the **finite-window linear-fit bias**
  is a separate effect that the 018 authors do not isolate. The
  window sensitivity data is a smoking gun for this.

**Recommended refinement:** the 018 attribution should be
re-stated as "the residual is dominated by (i) the unmodelled
annual solar forcing at 365.24 d (the dominant contribution,
which the secular formula uses time-averaged obliquity for but
the 1-year numerical integrates as a time-varying term), (ii) the
lunar evection + variation terms (smaller contribution), and
(iii) the finite-window linear-fit bias (which makes the W=365 d
measurement under-estimate the true secular limit)."

The true secular limit, extrapolated from W=730 d data, is likely
~+0.004 deg/day Lunisolar at h=600 km i_sso, which is **30×** the
corrected formula's +1.35e-4 deg/day, not **10×**. The 018
authors under-state the discrepancy by ~3× because they anchor
their ratio at W=365 d.

---

## 10. Cross-references to specific data points

| Observation | Section | Implication |
|---|---|---|
| Window sensitivity W=30 → W=730 monotonically increasing | §3.1 | Linear-fit bias makes W=365 d under-estimate secular by ~3-4× |
| Solar residual 33.7×, lunar residual 1.17× | §3.3 | Solar forcing is the dominant unmodelled term |
| i=90° ratio 2.81× vs i_sso 9.78× | §3.2 | Residual has J2-coupling component at i_sso |
| Sun_moon = Sun_moon_j2 to 0.03% | §3.3 | J2 × Lunisolar cross-term is < 0.03% relative |
| Precession on/off diff +0.012 deg/year | §3.4 | Frame mismatch is 3.6% of the residual |
| Convergence p_r=4.49, p_v=4.50 | §m | RK4 systematic error is negligible |

---

## 11. Limitations of this hostile review

1. I did not have access to a runtime to re-execute the 018
   numerics; all analysis is from the published `results.json`.
2. I did not read the audit-018 or other tracks' outputs (per
   constraint), so I may be duplicating or missing findings that
   other tracks have documented.
3. The order-of-magnitude estimates for short-period terms
   (evection, variation, annual forcing) are approximate and
   would benefit from explicit first-principles calculation.
4. The window-length extrapolation to W=1825 d is a linear
   extrapolation of the observed +0.005 deg/day trend from W=30
   to W=730. The actual W=1825 slope may be different if the
   trend saturates or reverses.
5. The mean-vs-osculating offset (candidate f) is hard to bound
   without an explicit derivation. The order-of-magnitude
   estimates I provided are placeholders.

---

## 12. Final summary for the synthesis lead

**Track G finding:** the 018 attribution of the 10× residual to
"unmodelled short-period terms" is qualitatively correct but
incomplete. The hostile review identifies **annual solar forcing
(c) + finite-window linear-fit bias (e) + lunar evection/
variation (a + b)** as the dominant surviving candidates, with
**annual solar forcing and finite-window linear-fit bias being the
strongest contributors**.

The 018 authors' emphasis on lunar evection/variation is wrong:
the data shows the SOLAR 33.7× residual is the dominant signal,
and the LUNAR 1.17× residual is small. The corrected secular
formula's lunar contribution (+9.91e-5 deg/day) is actually
*larger* than the corrected solar (+3.56e-5 deg/day), but the
numerical shows the solar contribution dominates by ~10× — which
is the opposite of what the corrected formula predicts.

The 018 conclusion that the residual is "dominated by short-period
terms" should be refined to specify **which** short-period terms
(the annual solar forcing) and to add the **finite-window
linear-fit bias** as a separate mechanism.

The true secular Lunisolar rate at h=600 km i_sso, extrapolated
from W=730 d data, is likely ~+0.004 deg/day (not the +0.00132
deg/day at W=365 d), which means the corrected formula
under-estimates by ~30× (not 10×).

End of Track G hostile review.