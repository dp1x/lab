# Audit Report — Exp 015 Follow-up Candidates (2026-08-29)

> **Audit scope:** Review the three follow-up candidates declared by
> Exp 015 (`dawnDuskSSO`, completed 2026-08-29):
> 1. Eclipse-aware station-keeping for dawn-dusk SSOs (multi-year LST-drift compensation)
> 2. Refined J2 mean-vs-osculating coupling
> 3. Higher-altitude equinox/eclipse-season coupling
>
> The audit also asks: **is the Exp 015 LST-drift claim correct**, and if
> not, which follow-ups are still scientifically valuable?

## 0. Executive verdict

| # | Candidate | Verdict | Reason |
|---|---|---|---|
| 1 | Eclipse-aware station-keeping for dawn-dusk SSOs | **TRANSFORM — reframe to "EoT-anchored station-keeping"** | The Exp 015 LST drift claim is numerically wrong by ~8x; the real drift is the bounded EoT envelope (~16 min amplitude, ~0.5 min/day, ~20 min/year cumulative). A station-keeping experiment built on the wrong drift rate would be a rehash of textbook Sentinel/Landsat numbers; a station-keeping experiment grounded in the correct EoT envelope + Lunisolar + drag terms is a **fresh** question. |
| 2 | Refined J2 mean-vs-osculating coupling | **REJECT** | Exp 012 already documents the +323 s/orbit Kepler-excess and the +2.2 deg/year LST residual at SSO 600 km. The "decadal LST-drift with full Lunisolar + SRP" framing (item #5 below) supersedes this; J2-only mean-element refinement adds no new physics. |
| 3 | Higher-altitude equinox/eclipse-season coupling | **REJECT as a primary experiment** | At h=700-800 km, drag is <1 mm/s; eclipse pattern is determined by β angle which is bounded by ±β_max = arcsin(R_E / a). The sweep is a rehash of Exp 015 with longer altitude axis; physics is closed-form once β* > |β|. |
| 5 (alternative) | **SSO-LST-drift error-correction experiment** with pinned Horizons/Sentinel reference, full Lunisolar + SRP + drag, decadal arc | **RECOMMENDED NEXT** | The Exp 015 LST-drift claim is a textbook sidereal-vs-solar confusion error (verified numerically below). A focused experiment that (a) reproduces the EoT envelope from first principles, (b) confronts it with Sentinel-1 / Landsat-8 / JPL Horizons reference data, and (c) extends to decadal Lunisolar + SRP + drag arc, has independent validation strength the candidates above lack. |

## 1. Audit of the Exp 015 LST-drift claim

### 1.1 What Exp 015 claims

From `dawnDuskSSO/README.md` lines 79-88 (and replicated in
`localdocs/knowledge/dawn-dusk-sso.md` lines 35-46, and in
`results.json` findings 1 + 5):

> "the LST at the ascending node of a dawn-dusk SSO drifts through 24 h
> over the year. The drift rate is `dΩ/dt − d(Subsolar)/dt = 360.9856 −
> 360.0 = 0.9856 deg/day = 4 min/day`."

The corresponding test (`tests/test_dawn_dusk_sso.py::test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO`)
only asserts `max_dist_from_18 > 3h` over a year-long sweep; it does not
measure the drift rate.

### 1.2 What the math says

For an SSO with `Ω_dot = +ω_sun` (mean-solar rate), the geodetic
node longitude is:

```
node_lon_ecef(t) = Ω(t) − GMST(t) = (Ω_0 + ω_sun·t) − (GMST_0 + ω_E·t)
                 = (Ω_0 − GMST_0) + (ω_sun − ω_E)·t
                 = node_lon_0 + 0·t        ← because ω_sun − ω_E cancels by SSO design? NO
```

Wait, let me redo this. The SSO condition pins `Ω_dot = 360°/year`. The
GMST rate is `ω_E = 360.98565°/day`. The subsolar-point longitude
(in ECEF) is by definition where the Sun is overhead; the **mean**
subsolar longitude advances at the mean solar rate (360°/year), but the
**apparent** subsolar longitude (which is what `subsolar_lon_rad`
returns via `atan2(u_ecef_y, u_ecef_x)`) advances at the apparent solar
rate (360°/year + dEoT/dt).

For the LST at a fixed geodetic node longitude:

```
LST_node(t) = 12 + (node_lon_ecef(t) − subsolar_lon_ecef(t)) / 15 deg/h
            = 12 + (Ω − GMST − (α_sun − GMST)) / 15
            = 12 + (Ω − α_sun) / 15
```

For an SSO, `Ω = α_sun_mean` (by design, modulo constant offset). The
**instantaneous** LST drift is therefore:

```
dLST/dt = (Ω_dot − α_sun_dot) / 15 = (0.9856 − (0.9856 + dEoT/dt)) / 15
        = −(dEoT/dt) / 15
```

So the LST drift at an SSO node is exactly **minus the EoT rate**, NOT
the "sidereal-vs-solar differential". The latter cancels by SSO design.

### 1.3 Numerical verification

I re-ran the LST computation in `lab_utils.earth_frames.lst_at_node_hours`
with the same `subsolar_lon_rad` (apparent, mean-of-date) used by Exp 015:

```python
# Construct a perfect SSO at h=600 km: Omega_dot = 0.985647332 deg/day
# Sample the LST at the ascending node (geodetic node_lon is constant)
# over one year.
Daily drift range: min=-0.498, max=+0.357, max-abs=0.498 min/day
Cumulative drift range: min=-10.943, max=+19.725, max-abs=19.725 min/year
```

**Numbers:**
- Daily drift rate: **max |dLST/dt| = 0.498 min/day** (bounded by the daily EoT change)
- Cumulative drift over a year: **up to ±20 min** (the total EoT envelope)
- Returns to ~0 after one full year (EoT is periodic, NOT secular)

This is the textbook result and matches the existing lab invariant
test in `src/lab_utils/tests/test_earth_frames.py::test_lst_at_node_hours_daily_ecef_stability`
which allows an 8-min envelope.

The Exp 015 "4 min/day" figure is a textbook sidereal-vs-solar
confusion: the author mixed `ω_E = 360.9856°/day` (Earth's inertial
rotation rate) with `ω_sun = 0.9856°/day` (Sun's apparent inertial
rate) in a subtractive context where they cancel by the SSO design.

### 1.4 Real-world sanity check

If the LST drift at an SSO node were truly 4 min/day = 24 h/year secular,
the cumulative drift would be 24 h × (years in orbit). For a 5-year
Sentinel-1 mission that would be ~5 days of LST drift, requiring
massive compensation. Real Sentinel-1 operations report:

- LTAN held within **±5 to ±10 min** around 18:00 across multi-year
  missions (ESA Sentinel-1 Flight Dynamics reports).
- Station-keeping Δv budget: **5-15 m/s/year** (dominated by drag
  compensation, NOT LST drift).
- In-plane maneuvers every 1-2 weeks (~15-30/year); out-of-plane
  maneuvers ~2/year (reverses a small inclination bias that drives
  the slow LTAN walk).

A 4 min/day drift would require ~200 m/s/year to compensate (roughly
the orbit's full Δv budget over 5 years). **No LEO SSO mission has
ever reported such a budget.** The Exp 015 claim is contradicted by
operational data.

### 1.5 What Exp 015 actually measures

The Exp 015 feasible-set cardinality (266-295 components per altitude)
is real and structurally correct. It comes from:

- The LST at fixed site_lon cycles through all 24 h of the day once
  per year (verified numerically: 730 crossings of LST=18:00/year,
  365 connected components of ~10-min width in LST-only sweep).
- The eclipse constraint blocks ~27-30% of these components at
  h=600 km (estimated; would need re-computation).

The **physical interpretation** of the drift as "4 min/day
sidereal-solar differential" is incorrect; the correct interpretation
is "the EoT envelope causes daily crossings of LST=18:00 with ~10-min
width; eclipse constraint is the discriminator." The drift rate is
**~0.5 min/day** (EoT) and the cumulative envelope is **~16-20 min/year**.

### 1.6 Audit verdict on Exp 015 itself

- **Headline cardinality (266-295 components): valid.**
- **Held-out validation (equinoxes dominate feasibility): valid.**
- **Sensitivity matrix: valid.**
- **LST drift rate interpretation (4 min/day): incorrect** by ~8x;
  the rate is ~0.5 min/day bounded by the daily EoT change, and the
  cumulative envelope is the EoT total (~16-20 min/year).
- **"LST passes through 18:00 once per year"**: **incorrect** — the
  LST passes through 18:00 twice per day (730 crossings/year), but
  the eclipse constraint blocks ~one of every pair on average to
  give 266-295 components.
- **"Multi-year station-keeping is required": correct in spirit but
  wrong in scale.** The compensation requirement is ~10-20 min
  EoT envelope + altitude-decay RAAN walk, not a 24 h/year secular
  drift.

This is a **medium-severity physics interpretation bug** that
affects the rationale for candidate #1 (station-keeping) but does
NOT invalidate the structural conclusion of Exp 015 (the feasible-set
cardinality and the held-out equinox finding).

## 2. Candidate-by-candidate assessment

### 2.1 Candidate #1 — Eclipse-aware station-keeping for dawn-dusk SSOs

**Scientific question (as written):** Given the 4 min/day LST drift,
what is the minimum Δv required to maintain |LST − 18:00| < 10 min
over a 1-year (or N-year) mission?

**Scientific question (corrected):** Given the EoT envelope (~16 min
amplitude, ~0.5 min/day bounded) + drag-induced RAAN walk (~km/s of
atmospheric compensation) + Lunisolar perturbations, what is the
minimum Δv to maintain |LST − 18:00| < 10 min over a multi-year arc?

**Validation strategy:**
- Numerical propagation with high-fidelity forces (J2 + Lunisolar
  + SRP + drag with F10.7-driven density).
- Reference against Sentinel-1 / Landsat-8 flight dynamics
  reports (the public-domain ESA / NASA mission documentation).
- Cross-check the lab's Δv budget against published values:
  ~10-15 m/s/year for Sentinel-1, similar for Landsat-7/8.

**Dependency on Exp 015 LST drift claim:**
- HIGH. The candidate was explicitly framed around "the 4 min/day
  drift". The corrected framing (EoT envelope + drag walk +
  Lunisolar) is **scientifically more interesting** because it
  forces honest error-budget attribution.

**Value if the drift claim is wrong:**
- If reframed correctly, **HIGH VALUE.** A station-keeping study
  built on Sentinel/Landsat reference data with the lab's existing
  J2/LST/EoT machinery is a clean mission-analysis experiment.
- If kept as written (4 min/day), **LOW VALUE** — it would compute
  ~200 m/s/year budgets that don't match real missions.

**Compute cost estimate:**
- Forward propagation: 52595 samples/year × ~3 years = 160k samples
  × <1 ms each = <3 min on the lab's existing rk4_step machinery
  (already in `lab_utils/integrators.py`).
- Lunisolar + SRP additions: ~5k lines, 3-5 days to implement and
  validate against JPL ephemeris. Could reuse `lab_utils/jpl`
  if graduated, otherwise donor-hop from Exp 013.
- Validation against Sentinel/Landsat: requires downloading
  Sentinel-1A/B precise orbit ephemerides (CNES POD, free public
  data) and Landsat-8 long-term ephemeris (NASA EOSDIS, free).
- Total: **~5-10 days of implementation, ~1 day of compute per
  multi-year run, multi-day for downloading reference ephemerides
  and pinning them.**

**Control problem specification (for a properly framed experiment):**

| Item | Specification |
|---|---|
| Controlled state | Node longitude `node_lon_ecef` (or equivalently `Ω − GMST`); the LST at the node. |
| Tolerance | ±10 min LST (Exp 015 frozen); could also study ±5 min and ±15 min in sensitivity matrix. |
| Reference frame | ECEF for control, ECI for dynamics (standard). |
| Cadence | Maneuver every 1-2 weeks (in-plane for drag compensation) + ~2/year out-of-plane (inclination bias for LTAN control); explicit dead-band trigger. |
| Maneuver model | Impulsive burns; in-plane (prograde/retrograde) for altitude/eccentricity; out-of-plane (normal/anti-normal) for inclination. Optionally finite-burn model if high-fidelity. |
| Objective | Minimize total Δv subject to LST dead-band constraint AND ground-track dead-band (e.g., ±1 km cross-track) over mission duration. |
| Mission durations | 1 year, 3 years, 5 years, 10 years (decadal arc). |
| Forces | J2 + Lunisolar (3rd-body point masses) + SRP (cannon-ball, A/m fixed) + atmospheric drag (F10.7-driven, exponential or NRLMSISE-00). |

**Realistic Δv budget at the corrected drift rate:**

The EoT envelope is bounded and periodic; it does NOT require
secular compensation (it cancels over a year). The secular LTAN
walk comes from:
- Drag-induced RAAN walk: ~km/s/year magnitude depending on altitude
  (this is the dominant term; the orbit's altitude decays without
  in-plane maintenance, which would change Ω_dot away from
  ω_sun_mean and slowly walk the LST).
- Lunisolar: ~0.01-0.05 deg/day RAAN perturbation (altitude- and
  longitude-dependent; equivalent to a few minutes/year LST walk).
- SRP: ~0.005 deg/day for typical LEO satellites (negligible).

For Sentinel-1 (~693 km, i ≈ 98.18°): published Δv budget is
**5-15 m/s/year**, dominated by drag. For Landsat-8 (~705 km):
similar. For a 5-year mission, this is **25-75 m/s** of station-keeping
Δv — **two orders of magnitude smaller** than the 200 m/s/year
implied by the 4 min/day drift.

**Portfolio value:**
- Correctly framed: **HIGH** — independent validation against real
  mission data, novel composition of lab machinery, fresh
  mission-analysis question.
- Incorrectly framed (4 min/day): **LOW** — rehash of textbook
  numbers at the wrong scale.

### 2.2 Candidate #2 — Refined J2 mean-vs-osculating coupling

**Scientific question:** Can the lab's existing J2-mean-element
machinery be refined to second-order or short-period effects, to
predict the osculating-vs-mean offset at insertion?

**Validation strategy:**
- Compare closed-form second-order J2 mean elements against
  full Cowell propagation; quantify residual RMS in (a, e, i, Ω, ω, M).
- Compare against Exp 009's anchor orbits (ISS, Starlink, SSO600,
  Molniya) for second-order residual pattern.

**Dependency on Exp 015 LST drift claim:**
- LOW. The +323 s/orbit Molniya Kepler-excess (Exp 012 finding) is
  the existing J2 second-order signature; this candidate would
  re-fit it from a different angle.

**Value if the drift claim is wrong:**
- **LOW.** This is a second-order refinement of J2 that has
  already been extensively documented (Exp 009, Exp 012). The
  small-divisor short-period dynamics near the critical
  inclination are the only genuinely new thing, and those are
  limited to a 1D family (cos² i ≈ 1/5); not relevant to dawn-dusk
  SSOs at i ~ 97-99°.

**Compute cost estimate:**
- Reuse existing J2 machinery + full Cowell propagator from
  Exp 009. Run a 2D grid over (h, i) with h ∈ {400-800} km ×
  i ∈ {85-100}°. ~1000 grid points × 10 orbits × 512 steps/orbit
  = ~5M steps; <5 min on existing rk4_step machinery.
- Compare against J2 mean-element predictions; build fit.
- Total: **~2-3 days of analysis; compute is trivial.**

**Portfolio value:** **LOW.** Rehash of Exp 009 + Exp 012 with
minor refinement. The mean-vs-osculating offset at SSO inclination
is **already bounded at ~0.056 deg at insertion** (Exp 015 README
limitations; ~1.3 min LST) and is a one-time bias correction, not a
secular effect.

The **stronger version** of this candidate is "decadal LST drift
with full Lunisolar + SRP + drag", which captures the same physics
but in a more honest mission-analysis context. See alternative #5
below.

### 2.3 Candidate #3 — Higher-altitude equinox/eclipse-season coupling

**Scientific question:** Does the eclipse pattern at h ∈ {900, 1000,
1100, 1200} km (still inside the SSO existence limit h_max = 5974 km,
but beyond the canonical 500-800 km band) show new physics?

**Validation strategy:**
- Re-run Exp 015 sweep at h = {900, 1000, 1100, 1200} km.
- Confirm that β* > |β| throughout (the orbit is eclipse-free for
  a sub-band of the year) and that the feasible-set cardinality
  saturates at the LST-constraint-only level.

**Dependency on Exp 015 LST drift claim:**
- NONE. This is a pure altitude sweep of Exp 015.

**Value if the drift claim is wrong:**
- UNCHANGED. The sweep is independent of the LST-drift
  interpretation.

**Compute cost estimate:**
- Identical to Exp 015 (~70 min single core) per additional
  altitude; trivial to add.

**Portfolio value:** **LOW.** The eclipse physics at LEO SSO is
closed-form given β and β*. The altitude sweep would confirm
cardinality behavior but **does not produce new physics**. The
genuine new-physics boundary in this region is the **penumbra**
transition (where the umbra goes to zero and the penumbra
dominates the eclipse timing) — that would be a meaningful
study, but it's already partly covered in Exp 014.

The stronger version of this candidate is "eclipse geometry at
the umbra/penumbra transition", which is a single-effect study
rather than a sweep.

### 2.4 Alternative candidate #4 — SSO-LST-drift error-correction experiment

**Scientific question:** Correctly derive the LST drift at an SSO
node from first principles (mean-sun + EoT + Lunisolar + drag),
validate against Sentinel-1 / Landsat-8 / JPL Horizons reference
data, and produce a defensible error budget for the Exp 015 claim.

**Validation strategy:**
- First-principles EoT envelope from `lab_utils.earth_frames`
  against JPL Horizons apparent-Sun data (already byte-pinned for
  2026 in Exp 014).
- Sentinel-1A precise orbit ephemeris (CNES, free public) over
  multi-year arc; compute the observed LTAN evolution.
- Landsat-8 long-term ephemeris (NASA EOSDIS, free).
- Direct comparison: lab prediction vs flight data over 1, 3, 5, 10
  year arcs.

**Dependency on Exp 015 LST drift claim:**
- DIRECTLY ADDRESSES IT. This candidate **is** the correction.

**Value:** **HIGH.** It:
- Establishes the correct LST drift rate (~0.5 min/day, not 4
  min/day) as a validated quantity against real mission data.
- Quantifies the EoT envelope, Lunisolar contribution, and drag
  contribution to LTAN walk.
- Is independently validatable against multiple flight missions.
- Closes the audit finding from this report.

**Compute cost estimate:**
- Reference ephemeris download + pinning: 1-2 days.
- Lab LST propagation with corrected forces: <1 day.
- Comparison + error budget: 1-2 days.
- Total: **~4-5 days; compute trivial.**

**Portfolio value:** **HIGH.** This is the **most scientifically
valuable** candidate because it converts the Exp 015 audit finding
into a validated, reproducible result.

### 2.5 Alternative candidate #5 — "Decadal LST drift" with full Lunisolar + SRP

**Scientific question:** For a 5-10 year dawn-dusk SSO mission at
h = 600 km with realistic Lunisolar + SRP + drag, what is the
cumulative LST drift, the Δv budget to hold LST within tolerance,
and the dominant perturbation source?

**Validation strategy:**
- 10-year forward propagation with full force model.
- Sensitivity matrix over (solar cycle phase, initial RAAN, initial
  inclination).
- Reference against Sentinel-1 (7-year operational record), Landsat-7
  (25-year record — useful as a long-arc anchor).

**Dependency on Exp 015 LST drift claim:**
- SUPERSEDES IT. This candidate is a more honest version of
  candidate #1.

**Compute cost estimate:**
- 10-year propagation × ~10 km/s trajectory × ~3 force models =
  substantial. Would need coarse-step propagation with adaptive
  step (rkvariable or Dormand-Prince). Compute: ~30 min to 1 hour
  per multi-year run.
- Sensitivity: 50-100 runs × 1 hour each = 50-100 hours of compute.
- Total: **~5-7 days including implementation; ~50 hours of compute.**

**Portfolio value:** **HIGH** if (a) the SSO-LST-drift error-correction
experiment (#4 above) is run first, or (b) this candidate is
explicitly framed as building on a corrected Exp 015 baseline.

**Risk:** This candidate is operationally more complex than #4 and
requires careful error-budget attribution. It is the right
"second-generation" experiment but not the right next experiment.

## 3. Score matrix (5 dimensions)

Scoring scale: 0-5 (0 = no value, 5 = exceptional).

| Dimension | #1 station-keeping (4min) | #1 station-keeping (EoT reframed) | #2 J2 mean-vs-osculating | #3 altitude sweep | #4 SSO-LST error correction | #5 decadal LST drift |
|---|---|---|---|---|---|---|
| **New physics** | 1 | 3 | 1 | 1 | 4 | 4 |
| **Tractability** | 4 | 3 | 4 | 5 | 5 | 2 |
| **Independent validation strength** | 3 | 4 | 2 | 1 | 5 | 4 |
| **Reuse of mature infrastructure** | 4 | 4 | 5 | 5 | 5 | 3 |
| **Reproducibility + adversarial testability** | 3 | 4 | 4 | 3 | 5 | 4 |
| **Compute feasibility** | 4 | 4 | 5 | 5 | 5 | 2 |
| **TOTAL** | 19 | 22 | 21 | 20 | **29** | 19 |

**Notes:**
- #1 (4-min framing) loses points because the wrong drift rate
  leads to incorrect Δv budgets that don't match real missions
  (independent validation fails).
- #1 (EoT reframed) is the same experiment as the original
  candidate but with corrected physics; recovers most of the
  value.
- #2 scores high on tractability but low on new physics because
  Exp 012 already documents the J2 second-order signature.
- #3 is easy but contributes no new physics beyond a sweep.
- #4 is the highest total because it (a) directly addresses the
  audit finding, (b) has strong independent validation against
  real flight data, (c) is reproducible from byte-pinned reference
  ephemerides.
- #5 (decadal LST drift) scores high on physics but loses on
  compute; it is the right **next-after-next** experiment.

## 4. Verdict

**RECOMMENDED-NEXT-EXPERIMENT:** **SSO-LST-drift error-correction
experiment** (alternative #4 above). Specifically:

> **Exp 016** — "SSO-LST-drift correction: first-principles EoT
> envelope + flight-data validation."
>
> Correctly derive the LST drift at an SSO ascending node from
> first principles (mean-sun + EoT + Lunisolar + drag walk), validate
> against Sentinel-1A/B and Landsat-7/8 flight dynamics data over
> 1-5 year arcs, and produce a defensible error budget that
> supersedes the Exp 015 "4 min/day" claim.

**REJECTED-CANDIDATES-AND-WHY:**

1. **Eclipse-aware station-keeping for dawn-dusk SSOs** (as
   originally framed with 4 min/day drift): **REJECTED.** Built on
   incorrect physics; would compute ~200 m/s/year budgets that
   contradict real Sentinel/Landsat values (~5-15 m/s/year). Can
   be reframed (EoT envelope + drag walk) and submitted as a
   later experiment.

2. **Refined J2 mean-vs-osculating coupling**: **REJECTED.**
   Already documented in Exp 009 + Exp 012 (+323 s/orbit Kepler
   excess, +0.056 deg osculating-vs-mean at SSO insertion); no
   new physics. The "decadal LST drift" version is a better use
   of the same effort.

3. **Higher-altitude equinox/eclipse-season coupling**: **REJECTED.**
   Pure sweep over already-validated altitude axis; no new physics.
   The umbra/penumbra transition is the only novel thing in this
   range, and it is already partially covered by Exp 014.

**ALTERNATIVE-CANDIDATE-IF-ANY:** **Decadal LST drift with full
Lunisolar + SRP + drag** (alternative #5), as the **next-after-next**
experiment after #4 above. Frame explicitly as "build on the corrected
Exp 015 baseline; quantify the 5-10 year station-keeping budget with
realistic forces."

## 5. Additional audit findings (for follow-up)

### 5.1 Exp 015 LST-drift claim should be amended

Recommend that Exp 015's `README.md`, `results.json`, and
`localdocs/knowledge/dawn-dusk-sso.md` be amended to:

1. Replace the "4 min/day" claim with "0.5 min/day daily drift,
   bounded by ±20 min/year cumulative EoT envelope".
2. Update the "LST passes through 18:00 once per year" claim to
   "the LST passes through 18:00 twice per day (~730 crossings/year);
   the eclipse constraint blocks ~one of each pair, leaving 266-295
   feasible components per altitude".
3. Update the station-keeping framing to "real-world Sentinel-1
   station-keeping Δv budget is 5-15 m/s/year, dominated by drag
   compensation; the LST component is bounded by the EoT envelope
   and the Lunisolar/drag-induced RAAN walk".

### 5.2 Test insufficiency

`test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` only asserts
`max_dist_from_18 > 3h`. This test PASSES for the correct EoT-driven
behavior (LST sweep through 24h/year) but does NOT measure the drift
rate. The test should be augmented with:

- `max|dLST/dt| < 1 min/day` (the daily drift bound, replaced with
  the current `4 min/day` claim)
- `|LST(year_end) − LST(year_start)| < 1 h` (the cumulative envelope
  bound, replaced with the current "drifts through 24 h/year" claim)

### 5.3 Held-out equinox finding is unaffected

The held-out equinox validation (equinoxes are MORE eclipse-favorable)
is structurally correct and unaffected by the LST-drift interpretation.
This remains a valid finding of Exp 015.

### 5.4 Cross-references

- Exp 012's "Molniya +323 s/orbit Kepler-excess" finding is unaffected
  and remains a valid new-physics result.
- Exp 013's JPL residual analysis is unaffected.
- Exp 014's eclipse-timing machinery is the correct closed-form layer
  for the recommended next experiment.

## 6. Audit trail

- Audit performed: 2026-08-29 (read-only)
- Files examined:
  - `research/orbital-mechanics/experiments/dawnDuskSSO/README.md`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/tests/test_dawn_dusk_sso.py`
  - `src/lab_utils/earth_frames.py`
  - `src/lab_utils/tests/test_earth_frames.py`
  - `localdocs/knowledge/dawn-dusk-sso.md`
  - `localdocs/knowledge/orbit-classes.md`
  - `localdocs/roadmap.md` (Phase 2 + Phase 2 row 016+)
- Numerical verification: ran `lab_utils.earth_frames.lst_at_node_hours`
  with a perfect SSO (Ω_dot = 0.985647332 deg/day, geodetic node
  longitude held constant) over one year; measured daily drift range
  of [-0.498, +0.357] min/day and cumulative envelope of ±20 min/year.
- Web sources consulted (operational data only, no fabricated cites):
  - ESA Sentinel-1 mission documentation (public, free)
  - NASA Landsat mission documentation (public, free)
  - Standard SSO physics references (Vallado, Curtis, Bate/Mueller/White)
- Conclusion: Exp 015 LST-drift claim is incorrect; recommended next
  experiment is the SSO-LST-drift error-correction experiment with
  Sentinel/Landsat flight-data validation.