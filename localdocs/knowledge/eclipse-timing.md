# Eclipse Timing & Launch Windows

Knowledge note for Exp 014, complementing the experiment card
(`research/orbital-mechanics/experiments/eclipseTiming/README.md`).

## What is the closed-form observation that makes event finding here trivial?

The lab's prior experiments standardized on **closed-form Kepler propagation**
via the canon (`solve_kepler`, `coe_to_rv_eci`, `true_anomaly_from_E`).
That means a state `r(t)` is available at *any* time, not just at integration
grid points. Therefore the shadow event function `g(t) = α_E(r(t), t) −
α_S(r(t), t) − θ(r(t), t)` can also be evaluated at any `t` — without
interpolation, without re-propagation, without a sub-step of any integrator.
This **decouples the event-time error from any integration step entirely**.

That changes the standard "ODE event finder" problem:
- Traditional ODE: bracket then refine with a *dense output* (Hermite or
  polynomial) approximation of the integrator, then re-bisect in the dense
  output. Event-time error is dominated by interpolation, then bisection,
  then integration.
- Lab case (Exp 014): bracket on a coarse grid, bisect to bracket width,
  period. Event-time error is dominated only by (a) bracket width (`Δw`
  in seconds, stops at 10⁻⁸ s) and (b) the absolute-epoch float-ULP floor
  (~1.5 s at J2000 + 26 yr × 365 d ≈ 1.2 × 10⁹ s).

The result is that the G4 convergence ladder can halve scan density 8×
and the entry event shifts by < 30 s — orders of magnitude tighter than
traditional integrators can claim for the same computational budget.

## Why dual formulations?

Route A (`g_A = α_E − α_S − θ`, apparent angles) and Route B (`g_B =
min(x, radius − ρ)`, shadow-axis algebra) are algebraically independent
presentations of the same spherical geometry:
- Different intermediate variables (angles vs distances)
- Different failure modes (degree/radian confusion in A; units/km ↔ AU
  in B's `tau = (R_S − R_E)/d_SUN`)
- Different cancellation patterns (A's dot-product in `cos θ` vs B's
  `ρ² = r² − x²`)

They must agree on event times because they describe the same surface.
When they disagree, the disagreement is either a real bug (caught) or a
shared failure mode (Sun model or state source — those need the
real-trajectory gate, not the routes themselves).

## Why both cylinder and cone?

The cylinder `R_cyl = R_E` is the d_SUN → ∞ limit of the cone `R_cyl
− x·tan(δ_u)`. The cone is more correct at high altitude; the cylinder
is more tractable as an analytic oracle and as a secondary check. They
are NOT the same model. The cylinder-vs-cone timing gap is a real
result, not an error bar:

| Altitude | Cylinder vs cone duration gap |
|---|---|
| LEO 420 km | ~8 s (umbra) |
| 1000 km | ~15 s |
| 2000 km | ~20 s |
| Molniya apogee (~46 200 km) | ~5 min + comparable partial-phase band |
| GEO 42 164 km | ~63 s per boundary; full umbra ~2.5 min shorter |

The conical correction grows with depth into the cone; a cylindrical
finder is acceptable as LEO primary and as universal negative control,
but at GEO it understates umbra by ~3.7 %.

## Why cylindrical is not "more correct" universally

The cylinder assumes parallel solar rays. The cone tracks the angular
extent of the Sun. At LEO, `tan(δ_u) = (R_S − R_E)/d_SUN ≈ 0.0046` rad,
so the cone radius deficit at LEO depth is 31 km — a fraction of a
percent. At GEO it's 194 km — a 3.7 % effect. At the cone tip itself
(r ~ 1.37 × 10⁶ km), the umbra ends; this is the "annular regime"
threshold beyond which the Earth disk can never cover the Sun disk
completely. None of the lab's cases reach this regime, but the model
must encode the threshold as a typed sentinel (e.g., a `NO_UMBRA`
return value) rather than fabricate an event.

## Why the Sun model lives in of-date, not ICRF

The Astronomical Almanac low-precision formulas return:
- mean longitude `L` and mean anomaly `g` (analytic functions of time)
- ecliptic longitude `λ` (after the equation of center)
- geocentric ecliptic unit vector `(cos λ, sin λ, 0)`
- rotated by mean obliquity of date `ε(t)` into equatorial of date

This is a *geometric* position in *mean equator/equinox of date*, NOT
in ICRF. The lab's pseudo-inertial ECI is J2000-fixed. The mismatch
causes a secular ~50.3″/yr precession bias (~0.14°/decade). For an
analysis window of a few days, the bias is < 0.5″, ~0.002 s event impact —
negligible. For a multi-decade mission, the bias grows linearly; the
contract documents this and the 2026 validation gate band (0.7 deg)
absorbs both the secular term and the dominant nutation term.

## Why the analytic Sun model and not a tabulated ephemeris

A truncated VSOP87 would give ~1″ accuracy at the cost of dozens-to-
hundreds of harmonic terms. The Almanac low-precision formulas give
~0.01 deg (= 36″) at zero tables, two trig calls per axis. The
quality factor is the **mapped** event-time error: at LEO, 36″
Sun-direction error maps through the mean-motion rate of 0.0658°/s
to ≤ 0.55 s worst case. That's inside the bracket-width convergence
floor (~1 s) and below the band target. The conclusion: at LEO event
tolerance, the analytic model is sufficient and the complex model is
wasted precision (AGENTS rule 7: complexity must justify itself).

## Why the launch-window constraint is "zero umbra entries in N revs"

A real mission wants eclipse-free operation in its early phase (post-
insertion through commissioning). Defining the window as a connected
component of the *negation* of "umbra entry in first N revs" gives:

1. A boolean indicator per launch time (`window_constraint(t_L) ∈ {0, 1}`)
2. A grid of indicator values → boolean transitions = window edges
3. Each transition is a one-dimensional root of the indicator; bisection
   refines to the resolution band

This avoids the "what does 'good launch time' mean" ambiguity that
plagues operational folklore. The "N revs" number is declared per case
(14 for the fine equinox day; 5 for the year sweep). The chosen N is
naturally mission-specific; the machinery is N-agnostic.

## Why the first 4 ISS events agree and the 3-day tail drifts

Exp 013 measured −7.08 km/day secular trend between the J2 model and the
TLE-provenance reference states. That's a boundary-rate of 2.65 km/s ×
2.7 s/day ≈ 7.1 km/day, consistent with the documented reference
envelope. Mapped to event times:
- Day 0: bias ~0 s (start)
- Day 1: bias ~3 s
- Day 2: bias ~6 s
- Day 3: bias ~9 s

The observed first-4-events agree to **5.5–13.5 s**, consistent with
the predicted day-0 to day-1 bias. The observed 308 s at the tail is
~30× larger than the secular prediction, suggesting a *real*
operational event in the ISS trajectory (reboost, dock/undock,
attitude change) that the radial second-difference gate catches
marginally (276 m over 240 s ≈ 1.1 m/s radial transient — likely a
reboost). Reported verbatim, not tuned away.

## Linkage to prior notes

- **Eclipse seasons** connect to **ground tracks** (Exp 008): the
  Sun's mean motion 0.9856°/day sets the cadence at which a given
  RAAN revisits the same β; that is the source of the two 7-week
  GPS eclipse seasons.
- **β-cutout threshold** `|β| ≤ arcsin(R_E/r)` is the same condition
  that appears in the J2 nodal-rate derivation (`-3/2 n J2 (R/p)² cos i`):
  the **critically inclined Molniya orbit** (i ≈ 63.4°) has
  `cos i = 1/5`, leading to zero apsidal precession and *maximum*
  eclipse dwell — opposite of "cutout" but in the same geometric family.
- **J2 mean-element secular rate** used in the launch-window
  constraint is the same first-order form whose rediscovery was
  Exp 009's headline. The N-rev-constraint evaluator folds the
  rate into a per-node RAAN (mean-element staircase), which is
  exact for a circular orbit under first-order J2.
- **Solar rate pinned at 0.985647332099 deg/day** (Exp 009/012 anchor)
  corresponds to `360 / 365.2422` (Gregorian mean-solar year) and is
  reused here without re-derivation.
- **Nodal regression rate** for the constraint is from the same first-
  order J2 formula as Exp 009 / orbitClasses.
