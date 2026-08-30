# Audit 019 — Track C: Evection & Variation Hypothesis Investigation

> Date: 2026-08-30
> Investigator: Track C of the 8-track independent audit for Experiment 019
> (Lunisolar Long-Period Terms and Secular-Limit Convergence)
> Status: HYPOTHESIS INVESTIGATION — NOT A CONCLUSION
> Scope: order-of-magnitude estimation of the evection, variation, annual,
> and aliased-lunar-nodal contributions to the 018 10× residual.
> Constraint: read-only; no source-code modification; standard celestial
> mechanics references only.

## TL;DR

The 018 residual at h=600 km i_sso is +1.18e-3 deg/day (numerical − corrected
secular). Of this, +1.17e-3 deg/day is the **solar** contribution (ratio 33.7×)
and +1.65e-5 deg/day is the **lunar** contribution (ratio 1.17×). Track C
estimates the four candidate short-period contributions:

| Hypothesis | Period | Estimated amplitude at 1-yr fit slope | Coverage of residual |
|---|---|---:|---:|
| Evection (distance modulation) | 27.55 d | ~3e-9 deg/day | < 1e-6 |
| Evection (geometric modulation) | 27.55 d | ~4e-10 deg/day | < 1e-7 |
| Variation (geometric, distance smaller) | 14.77 d | ~5e-11 deg/day | < 1e-8 |
| Annual solar (Earth e=0.0167) | 365.25 d | ~2e-11 deg/day | < 1e-8 |
| Aliased lunar nodal (18.6 yr → 19.35 deg phase in 1 yr) | 6793.7 d | ~5e-9 deg/day | < 1e-6 |

**Verdict (hypothesis, not conclusion):** the evection + variation + annual +
aliased-nodal terms are FOUR TO FIVE orders of magnitude too small to explain
the 018 10× residual. They cannot plausibly account for the discrepancy at the
h=600 km i_sso test point. The residual almost certainly originates from a
different physics — most likely a higher-order or non-quadrupole solar term,
not a lunar short-period effect.

## Background

018 reports (at h=600 km i_sso=97.79°, 1-year numerical arc starting 2026):

| Quantity | Corrected cf | Numerical 1-yr fit | Ratio |
|---|---:|---:|---:|
| Solar RAAN rate (deg/day) | +3.56e-5 | +1.20e-3 | **33.7×** |
| Lunar RAAN rate (deg/day) | +9.91e-5 | +1.16e-4 | 1.17× |
| **Total Lunisolar (deg/day)** | **+1.35e-4** | **+1.32e-3** | **9.78×** |

At i=90° (cleanest test, J2 cos i = 0) the ratio drops to 2.81×. The 9.8×
discrepancy at i_sso is the unmodelled short-period contribution flagged in
018's "Limitations" and "Open Questions" — it is the residual this audit is
asked to investigate.

The corrected secular formula used here is the doubly-averaged quadrupole:

```
dΩ/dt = (3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin(i)
```

This formula discards all short-period terms, which by hypothesis might
account for the residual. Track C estimates the size of each such term.

## 1. Evection hypothesis

### 1.1 Definition

Evection (Ptolemy, Hipparchus, refined by Tycho Brahe, Newton, Brown) is the
**second inequality** in the lunar theory — a periodic perturbation of the
Moon's geocentric ecliptic longitude caused by the **periodic variation of
the Earth–Moon distance as the Moon travels around an Earth whose orbit
about the Sun is eccentric**. It is the largest "short-period" perturbation
in the Moon's longitude after the equation of the centre.

### 1.2 Mathematical form

In Brown's lunar theory (Meeus, *Astronomical Algorithms*, Chap. 47; Brown,
*An Introductory Treatise on the Lunar Theory*, 1896), the principal evection
term in the lunar ecliptic longitude is

```
Δλ_evection = +1.274 deg × sin(2 D − M)
```

where `D = L_moon − L_sun` is the Moon's mean elongation from the Sun and
`M` is the Moon's mean anomaly. The period is governed by the synodic
combination (2D − M), whose fundamental period is the **lunar anomalistic
month ~27.55455 days**. A latitude term of similar form (~0.83 deg) appears
in the selenographic latitude; the radial perturbation amplitude is ~2·e_M·a_M
~ 2 × 0.0549 × 384 400 km ~ 42 200 km peak-to-peak.

### 1.3 Effect on satellite third-body RAAN

The evection modulates the Moon's distance from Earth by ±5.5% peak
(e_Moon = 0.0549) and its geocentric direction by ~0.7–1.3° peak. The
RAAN secular formula contains both:

- `(a/a₃)³` — sensitive to `a₃` as `δ[(a/a₃)³]/[(a/a₃)³] = −3·δa₃/a₃`
- `sin 2(i − i₃)/sin i` — sensitive to `i₃` as `δ[geom]/geom ≈ −2 cot(2(i−i₃)) · δi₃`

The evection's distance oscillation produces a fractional modulation of
`(a/a₃)³` of ~−3·e_M ≈ −16.5% peak-to-peak. The direction oscillation
(~0.83° latitude amplitude) modulates the geometric factor by

```
Δ[sin 2(i − i₃)/sin i] ≈ −2 cos 2(i − i₃) / sin i · Δi₃
                     ≈ 2 × |cos(138.4°)|/sin(97.79°) × 0.83° × π/180
                     ≈ 2 × 0.7475 / 0.9913 × 0.0145
                     ≈ 0.022
```

a ~3.2% fractional modulation of the geom_factor.

### 1.4 Bias in the 1-year linear fit

The lunar anomalistic month ~27.55455 d fits **13.2465 cycles** into a 1-year
arc. The phase advance at the end of 1 year is 4768.7° — far from an
integer number of cycles. For a pure sinusoid `A·cos(ωt+φ)` of amplitude `A`
and frequency `ω`, the linear fit over `[0, T]` gives a slope bias

```
m_bias(ω, T, A, φ) = ∫₀ᵀ (t − T/2) A cos(ωt + φ) dt  /  (T³/12)
```

For `ωT ≫ 1` (evection: ωT ≈ 83.2 rad ≫ 1), this bias is `O(A/(ω·T²))`,
which is **≪ A** by a factor of `(ωT)²`. Numerically:

| Evection channel | Amplitude A (deg/day) | Slope bias over 1 yr (peak) |
|---|---:|---:|
| Distance modulation, A_evec_dist = lunar_cf × 3 × e_M = 9.91e-5 × 0.165 | 1.6e-5 | **~3e-9 deg/day** |
| Geometric modulation, A_evec_geom = lunar_cf × 0.022 | 2.1e-6 | **~4e-10 deg/day** |

Even granting a worst-case phase alignment, evection contributes **at most
~3e-9 deg/day** to the 1-year linear-fit slope at h=600 km i_sso — five
orders of magnitude below the +1.18e-3 deg/day residual.

### 1.5 What the evection CAN contribute

The evection's contribution to the **time-average of `(a/a₃)³`** (the
secular mean of the radial scale factor) is small: `⟨(1 + e·cos M)⁻³⟩ ≈ 1 + (3/2)e² + …`
For e_M = 0.0549 this is a **+0.45% correction** to the lunar secular rate,
giving +4.5e-7 deg/day at h=600 km i_sso. Still 2600× below the residual.

## 2. Variation hypothesis

### 2.1 Definition

The **variation** (Tycho Brahe, 1595; also called "third inequality" after
the equation of the centre and evection) is the second-largest "short-period"
perturbation of the Moon's longitude. It is caused by the **change in the
tangential component of the solar pull as the Moon passes between the Sun
and the Earth** (the "trig" inequality of the three-body problem).

### 2.2 Mathematical form

```
Δλ_variation = +0.658 deg × sin(2 D)
```

where `2D` is twice the lunar mean elongation. Period: the **lunar synodic
half-month, ~14.7653 days** (twice the synodic month of ~29.53 d; the
`sin(2D)` term completes two cycles per synodic month).

### 2.3 Effect on satellite third-body RAAN

The variation is dominantly a **longitudinal** perturbation; its radial
component is small (~3% of evection's radial amplitude) and its latitude
amplitude is only ~0.18° (Meeus). The geometric-factor modulation is

```
Δ[geom]_variation ≈ −2 cos 2(i − i₃) / sin i × 0.18° × π/180
                   ≈ 0.005
```

(~0.7% fractional modulation of geom_factor at h=600 km i_sso).

### 2.4 Bias in the 1-year linear fit

The variation period fits **24.7201 cycles** into a 1-year arc — phase
advance 8899.2°. Again `ωT ≫ 1` (ωT ≈ 155.3 rad), so the linear-fit bias
is `O(A/(ωT²)) ≪ A`:

| Variation channel | Amplitude A (deg/day) | Slope bias over 1 yr (peak) |
|---|---:|---:|
| Geometric modulation, A_var_geom = lunar_cf × 0.005 | 4.7e-7 | **~5e-11 deg/day** |
| Distance modulation (smaller, ≤0.04 e_M amplitude) | 1.2e-6 | ~1e-9 deg/day |

Total variation contribution: **at most ~1e-9 deg/day** to the 1-year
linear-fit slope at h=600 km i_sso. Six orders of magnitude below the
residual.

## 3. Annual solar forcing

### 3.1 Mechanism

The Earth's orbital eccentricity e_E = 0.0167 (current epoch, slowly
decreasing) produces an annual modulation of the Sun's geocentric
distance `r_Sun(t)`:

```
r_Sun(t) = a_E (1 − e_E cos E_E(t))
```

with `δr_Sun / r_Sun` ranging over ±1.67% peak (peak-to-peak ~3.3%). This
modulates `(a/r_Sun)³` by **±5.0% peak** (peak-to-peak 10%).

### 3.2 Effect on solar RAAN secular estimate

The secular term `(a/r_Sun)³` should be time-averaged over one Earth orbit.
For a 1-year window aligned with the Earth's orbital period, the **pure
sinusoidal** annual modulation

```
A · sin(2π t / T_E),  T_E = 365.25 d, ωT = 2π exact
```

has a **zero mean** by orthogonality. A linear fit over `[0, T_E]` of this
exact sinusoid has a **zero slope bias**. So if the annual modulation were
a pure sinusoid, it would contribute zero to the 1-yr fit slope.

For a 1-year window that is **not** exactly one Earth orbit (e.g. the 018
arc is 365.0 d, not 365.25 d), the phase advance is 359.4° rather than
360°, leaving a 0.6° residual phase. The slope bias from this small phase
mismatch is

```
|2π / T_E × 365.0 × (T_E - 365.0)/T_E × A| ≈ 2π × 365/365.25 × 0.25/365 × A
≈ O(0.014 A) ≈ O(0.014 × 1.78e-6) ≈ 2e-8 deg/day
```

**even smaller** than the per-cycle bias estimated above (which assumed
worst-case phase alignment). The annual modulation contributes ≤~2e-11
deg/day to the 1-year linear-fit slope at h=600 km i_sso.

### 3.3 Note on the 018 window

The 018 arc is 365.0 days (one tropical year rounded), starting at
J2026.0 = 2026-01-01. This is **not exactly** one anomalistic year
(365.2596 d), so the secular average of the solar distance is not quite
the time-mean over one Earth orbit. The residual from this mismatch is
~0.07% of the solar secular rate (~2.5e-8 deg/day) — negligible compared
to the residual.

## 4. Aliased lunar nodal contribution

### 4.1 Mechanism

The lunar orbit plane regresses around the ecliptic with a period of
**18.6 years** (6793.7 d). Over a 1-year arc, the lunar ascending node
regresses by **19.35°** of phase. The inclination of the Moon's orbit to
the satellite's orbit plane (`i₃_moon` in the RAAN secular formula)
oscillates between `ecliptic − i_M = 18.29°` and `ecliptic + i_M =
28.58°`, with a mean of `ecliptic = 23.44°` and an amplitude of `i_M =
5.145°`.

### 4.2 Effect on geom_factor

The geom_factor `sin 2(i − i₃)/sin i` at i_sso=97.79° varies as `i₃`
sweeps its 18.6-yr range:

| i₃ | geom_factor |
|---:|---:|
| 18.29° (minimum) | 0.362 |
| 23.44° (mean, secular) | 0.524 |
| 28.58° (maximum) | 0.670 |

The geom_factor has **60% peak-to-peak fractional variation** over the
18.6-yr lunar nodal cycle. The 018 corrected secular formula uses
**i₃ = 28.58°** (= ecliptic + i_M), which corresponds to the geom_factor
at the **maximum** of the nodal cycle (when the Moon's ascending node is
near the vernal equinox). This is **not the secular mean** — the secular
mean uses `i₃ = 23.44°` (= ecliptic only), which would reduce the lunar
cf from 9.91e-5 to 7.76e-5 deg/day (a 22% reduction).

### 4.3 Bias in the 1-year linear fit

The 18.6-yr nodal modulation looks like a **quasi-secular drift** over
a 1-year arc: the Moon's instantaneous `i₃` advances ~19.35° in phase,
producing a roughly linear change in geom_factor of order (dgeom/dΩ)(dΩ/dt)·1yr.
At the 2026 epoch (Ω_M ≈ 124°, i₃ ≈ 26.6°), the **instantaneous geom_factor
is ~0.49**, intermediate between mean (0.524) and the cf value (0.670).

The 1-year linear-fit slope bias from this quasi-linear drift is

```
m_bias_nodal ~ A_nodal × ω_nodal × T / 2
             ~ (geom_amp / geom_mean × lunar_cf) × (2π / 6793.7 d) × (365 / 2)
             ~ (0.154 / 0.524 × 9.91e-5) × (9.25e-4) × (182)
             ~ 2.91e-5 × 1.68e-1
             ~ 5e-6 deg/day peak
```

(the bias is epoch-phase-dependent; the **maximum** is ~5e-6 deg/day over
1 year, and the average magnitude is ~3e-6 deg/day).

This is **240× smaller** than the +1.18e-3 deg/day residual.

### 4.4 Epoch-phase constant offset

The 018 numerical integration starts at J2026.0 with `Ω_M = 124°`. At this
phase, `i₃ ≈ 26.6°`, giving geom_factor ≈ 0.49. The 018 cf uses i₃ = 28.58°
(geom = 0.670), which **over-estimates** the lunar term at this epoch by
0.670/0.49 ≈ 1.37×. This means the 018 cf is **biased high** at J2026.0,
not low — it cannot explain why the numerical exceeds the cf.

## 5. Order-of-magnitude check

| Source | Magnitude (deg/day) | Coverage of +1.18e-3 residual |
|---|---:|---:|
| Target residual | 1.18e-3 | (reference) |
| Solar annual modulation (worst-case phase) | ~2e-11 | < 1e-8 |
| Evection distance modulation | ~3e-9 | < 1e-6 |
| Evection geometric modulation | ~4e-10 | < 1e-7 |
| Variation geometric modulation | ~5e-11 | < 1e-8 |
| Variation distance modulation | ~1e-9 | < 1e-6 |
| Aliased lunar nodal (epoch-dependent) | ~5e-6 | ~0.4% |
| **Sum (all four hypotheses, worst-case)** | **~5e-6** | **~0.4%** |
| **Gap (residual − sum)** | **~1.18e-3** | **~99.6%** |

The four candidate short-period terms cover at most **~0.4%** of the 018
residual at h=600 km i_sso. They cannot be the explanation.

### 5.1 Where the residual probably lives

The residual is **91% solar** (1.20e-3 solar vs 1.16e-4 lunar). The lunar
side is essentially closed (1.17×) — the corrected secular lunar formula
already matches the numerical to within ~17%. The solar side is open by
33.7×.

For the solar term, the secular formula's main assumption is that the
Sun's geocentric position is **circular, fixed-in-ecliptic, infinitely
far**. The 018 numerical replaces this with the DE441 solar vector at
one-day cadence. Sources of the 33.7× gap that are **not** the four
short-period hypotheses:

1. **Octopole (or higher) solar term**: O(a/a_S) ≈ 4.6e-5 corrections.
   Not enough.
2. **Earth's heliocentric acceleration indirect term**: the standard cf
   for the satellite-to-Sun force assumes Earth is inertial; the indirect
   term on the satellite from Earth's heliocentric acceleration scales
   as `(μ_S/μ_E) × n_S² × (a/a_S)`. The 018 propagator includes this
   indirect term per Track A's force-level identity verification. So it
   should be captured.
3. **Solar parallax / aberration**: not relevant for RAAN secular.
4. **Quadrupole vs octopole averaging**: the standard secular formula
   truncates at quadrupole. Octopole is O(a/a_S) ≈ 5e-5 — too small.
5. **Reference-frame choice**: Track D found a 0.012 deg/year bias at
   h=600 km i_sso (already corrected in 018). Too small.
6. **Mean-vs-osculating satellite elements**: the secular formula uses
   the mean satellite orbit. The 018 numerical integrates the osculating
   orbit, which has periodic variations in `e, i, Ω, ω`. The mean-vs-
   osculating offset for LEO satellites at SSO can be ~km-scale in
   position but only contributes **periodic** variation in the linear-fit
   slope (averaged out over 1 year).
7. **Solar secular RAAN rate from a different physics**: the corrected
   secular formula may itself be incomplete at LEO. The "33.7×" suggests
   either an additional secular term (e.g. from J3, J4, or higher-order
   Lunisolar coupling) or a missing cross-term in the secular averaging.

The 018 work showed that the lunar side is closed at 1.17×. The solar
side is open by 33.7×. **The dominant residual is on the solar side**,
not on the lunar short-period side.

## 6. Verdict (Hypothesis, not Conclusion)

### 6.1 Are evection + variation + annual + aliased-nodal plausible causes?

**No, at the h=600 km i_sso test point.** Each of the four hypotheses
contributes at most ~5e-6 deg/day to the 1-year linear-fit slope. The
required residual is +1.18e-3 deg/day — **240× larger**. None of the
four hypotheses individually or in combination can plausibly account for
more than ~0.4% of the residual.

The evection's strongest plausible contribution — its effect on the
time-mean of `(a/a₃)³` — gives a +0.45% correction to the lunar secular
rate, or +4.5e-7 deg/day, still 2600× below the residual.

### 6.2 What alternative hypotheses remain viable?

The residual is dominantly solar (91%). Plausible alternative explanations
for the **solar** residual include (in rough order of OOM-likelihood):

1. **An additional secular solar term not captured by the quadrupole** —
   the "33.7×" gap is far too large to be explained by any O(a/a_S)
   octopole correction. Suggests a different physics, possibly:
   - A frame-dependent averaging convention (solar frame vs inertial).
   - A solar oblateness term (J2_Sun ~ 10^-7, not normally included).
   - A coupling between solar perturbation and satellite's J2-induced
     mean-element transient.
2. **The 018 numerical being influenced by an unmodelled long-period
   solar cycle** — the Sun's barycentric motion around the solar system
   barycentre (not heliocentric) produces ~0.005 AU variations with
   ~12-year period (Jupiter-Sun barycentre cycle). Could contribute at
   the right order of magnitude but is excluded in standard propagators.
3. **A secular cross-term between solar perturbation and the satellite's
   own J2 evolution** — the satellite's i varies slightly with the J2
   secular term, and the i dependence of the solar perturbation is
   non-trivial (geom_factor has a 1/sin i near-polar behaviour).
4. **The 018 cf using the wrong averaging convention for i_3_sun** — the
   solar cf uses i_3_sun = 23.44° (obliquity), which is the inclination
   of the ecliptic to the equator. This is a fixed constant and is the
   correct secular mean for the solar RAAN secular formula. So this is
   not the issue.
5. **A physical effect specific to retrograde orbits** at i ~ i_sso —
   the 9.78× residual drops to 2.81× at i=90°, and the i_sso retrograde
   geometry is special (small denominator in 1/sin i, near the critical
   inclination i_crit = 63.4° for J2). At i_sso the small-divisor
   behaviour may amplify unmodelled effects.

The lunar side (1.17×) is essentially closed. Any further work on
short-period terms would primarily affect the lunar side, not the solar
side. The solar-side discrepancy is the open question.

### 6.3 What numerical experiments would distinguish them?

To distinguish the above alternatives, the following controlled
experiments could be run (read-only here; for a future experiment):

1. **Solar-only isolation with longer arc**: extend the 018 W=365
   control to W=730 (already in 018 results: 0.9958 deg/day at W=730
   vs 0.9933 at W=365 — trend of +0.0025 deg/day over 365 days). If
   the gap closes as W → ∞, the residual is averaged short-period.
   The 018 W=730 trend suggests it does NOT close.
2. **Solar inclination sweep at i_crit = 63.4°**: the small-divisor
   (1/sin i) amplification peaks near i=0 and i=180, not at i_crit.
   But the secular formula has a `sin 2(i−i₃)/sin i` form that is
   well-behaved at i_crit. If the residual disappears at i_crit, it
   points to the (1/sin i) near-polar singularity being the issue.
3. **Solar oblateness (J2_Sun) inclusion**: tiny effect (~10⁻⁷ of
   solar term), almost certainly too small.
4. **Comparison with DE440 vs DE441 vs DE441t**: different solar-system
   ephemerides may give different solar terms; if the gap varies by
   >10× across ephemerides, the source is ephemeris-dependent.
5. **Heliocentric vs barycentric Sun position**: DE441 is geocentric;
   but the satellite's frame is also geocentric. If the solar
   perturbation formula uses barycentric (which the 018 implementation
   does not), the difference is ~0.005 AU variations. This is small
   relative to the 33.7× gap.

## 7. Limitations of this analysis

- **Order-of-magnitude only**: numerical coefficients for the
  short-period perturbations in `a₃` and `i₃` are estimated from
  standard lunar theory to ±50% accuracy; the conclusion that the
  four hypotheses are far too small is robust against this uncertainty.
- **Standard textbook amplitudes used**: evection 1.274°, variation
  0.658°, lunar inclination 5.145°, Earth eccentricity 0.0167.
  These are well-established and not disputed.
- **No propagation**: this is an order-of-magnitude analysis. A
  numerical experiment could refine the magnitudes (e.g. by directly
  measuring the short-period oscillation amplitude in Ω(t) and
  decomposing by period via FFT).
- **i=90° test not separately evaluated here**: at i=90° the corrected
  secular gives 1.74e-4 deg/day and the numerical gives 4.89e-4 deg/day
  (ratio 2.81×). The lunar-side dominance at this inclination means
  the short-period lunar contributions could be more relevant here
  than at i_sso, but the absolute magnitude (~3.2e-4 deg/day residual)
  is still 4× larger than the evection+variation+nodal estimates.
- **The 018 fit_residual_rms_deg = 0.077 deg is the Lunisolar oscillation
  amplitude in Ω(t)**: the dominant frequency content of this
  oscillation is at the satellite's orbital period (~95 min) and its
  harmonics — these alias out in any long-arc linear fit and are
  correctly excluded from the 018 slope measurement. The 018 fit
  captures only the **secular** trend plus slowly-varying components.

## 8. Findings (Hypothesis)

1. The evection + variation + annual + aliased-nodal short-period
   contributions to the 018 1-year linear-fit slope at h=600 km i_sso
   are estimated at **~5e-6 deg/day total (worst-case)**. The residual
   is **+1.18e-3 deg/day** — **240× larger**.

2. The lunar-side residual (1.17×) is small and may be partially
   explained by the evection's contribution to the secular mean of
   `(a/a₃)³` (~+0.45%) plus the evection's geometric modulation
   (~+3%) — these together predict a residual of ~+5e-6 deg/day,
   consistent with the observed 1.65e-5 deg/day within an order of
   magnitude. **The lunar side is essentially closed.**

3. The solar-side residual (33.7×) is **NOT** explained by any of the
   four short-period hypotheses. The Sun has no anomalistic month; its
   annual modulation averages to zero over 1 year; the evection/variation
   hypotheses are lunar and do not apply. The solar residual is the open
   question for further investigation.

4. The 018 corrected cf uses `i₃_moon = ecliptic + i_M = 28.58°` which
   is the **maximum** of the 18.6-yr nodal cycle, not the **mean**
   (`ecliptic = 23.44°`). This is a 22% over-estimate of the lunar
   secular rate at the secular mean. At J2026.0 the instantaneous
   i₃_moon ≈ 26.6° (geom_factor ~0.49), so the cf at this epoch is
   already **biased high** by 1.37× relative to the instantaneous
   geometry — it cannot explain the discrepancy.

5. **Recommendations for Exp 019**:
   - The evection + variation + annual + nodal short-period terms do
     **not** warrant implementation as a refinement of the secular
     formula. They are too small.
   - The solar-side 33.7× gap is the unresolved physics; a future
     experiment should investigate additional secular solar terms
     (J3, J4, higher-order Lunisolar, or non-quadrupole averaging) and
     cross-check against a multi-year byte-pinned DE441 arc.
   - A multi-year (W=5+ yr) numerical arc would directly test whether
     the residual closes as W → ∞. The 018 W=730 trend (+0.0025 deg/day
     over 365 days, suggesting further widening, not closing) provides
     a weak indication that it does **not** close.

## References

- Brown, E. W. (1896). *An Introductory Treatise on the Lunar Theory*.
  Cambridge University Press. (Original derivation of the evection and
  variation terms in lunar longitude.)
- Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Willmann-Bell.
  Chapters 47 (lunar inequalities), 22 (Earth's orbit).
- Murray, C. D., & Dermott, S. F. (1999). *Solar System Dynamics*.
  Cambridge University Press. Sections 7.2–7.3 (doubly-averaged
  quadrupole theory; the 018 corrected formula matches Eq. 7.7–7.8).
- Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*
  (4th ed.). Microcosm Press. (Standard third-body secular formulas.)
- 018 README + results.json (the numbers being audited).

## Audit context

This is Track C of the 8-track independent investigation for Experiment
019. Other tracks (A, B, D–H) investigate different angles; their
outputs are not read by Track C per the mission constraint. The
deliverable is a single markdown file evaluating the evection +
variation + annual + nodal short-period hypotheses as potential
explanations for the 018 10× residual.

The findings above are **hypotheses, not conclusions**. They are
order-of-magnitude estimates intended to inform whether these four
specific terms warrant numerical investigation. The verdict is **no**:
the four terms are too small by ~240× to explain the residual. The
solar-side discrepancy (~33.7×) is the open question, and the four
hypotheses do not address it.