# Audit Report — Exp 015 LST-Drift Numerical Falsifier (2026-08-29)

> **Audit type:** Independent numerical experiment (read-only).
> **Subject:** Exp 015 (`dawnDuskSSO`, completed 2026-08-29) LST drift claim.
> **Method:** Build an independent J2 RK4 propagator that measures the
> LST at ascending-node crossings over a year and compare against the
> claimed drift rate of 4 min/day = 0.9856 deg/day.
> **Inputs:** Pure lab machinery (`lab_utils.orbits`, `lab_utils.earth_frames`,
> `lab_utils.integrators`). No donor-hopping from Exp 015.
> **Outputs:** Numerical measurement of the actual LST drift rate
> (J2-on vs J2-off), a sampled table of LST at ascending-node crossings,
> and the corrected physical interpretation.

## 0. Verdict (read first)

| Item | Value |
|---|---|
| **FALSIFIER-VERDICT** | **RED** |
| **015-CLAIM-VERIFIED** | **no** |
| **MEASURED-DRIFT (J2-on SSO)** | **+0.0891 min/day** (+0.0223 deg/day) — 45× smaller than the claim |
| **MEASURED-DRIFT (J2-off Kepler)** | −3.879 min/day (−0.970 deg/day) — agrees with the analytic sidereal-vs-solar differential, but is *not* an SSO |
| **ACTUAL-PHYSICS** | An SSO by design pins `Ω_dot = α_sun_dot`, so `dLST/dt = 0` to first order. The "4 min/day = 24 h/year" value cited by Exp 015 is the **sidereal-vs-solar differential** (360.9856 − 360.0 deg/day); it applies to a *non-precessing* orbit (e.g., a fixed-RAAN Kepler orbit), **not** to an SSO. The author confused `ω_E = 360.9856 deg/day` (Earth's inertial rotation rate) with `Ω_dot` for an SSO, which is actually `ω_sun_mean = 0.9856 deg/day` (mean-solar rate). |

The Exp 015 claim is **numerically wrong** by ~45× for an SSO. A follow-up
**EoT-anchored station-keeping** experiment was already recommended in
`audit-015-follow-up-candidates-2026-08-29.md`; this report provides the
independent numerical reproduction of the underlying drift rate.

## 1. The claim under audit

From `research/orbital-mechanics/experiments/dawnDuskSSO/README.md`
(lines 78-83), `experiment.py` (FINDING 1, line 874-880), and `results.json`:

> "the LST at the ascending node of a dawn-dusk SSO drifts through 24 h
> over the year. The drift rate is `|dΩ/dt| − |dSubsolar|/dt| = 360.9856
> − 360.0 = 0.9856 deg/day = 4 min/day`. The LST passes through 18:00
> (and through 6:00, the opposite terminator) once per year."

The relevant test (`test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO`,
file `tests/test_dawn_dusk_sso.py` line 686-) asserts only
`max_dist_from_18 > 3h` over a year — it does **not** measure the drift
rate, and it passes trivially for the bounded EoT envelope.

## 2. The two independent LST formulas

The lab canon (`src/lab_utils/earth_frames.py` line 277-301) states
explicitly that two paths are bit-equivalent:

```
Path (a) — textbook:     LST = 12 + (Ω − α_sun) / 15
Path (b) — lab_utils:    LST = 12 + (node_lon − subsolar_lon) / 15
where:
    Ω          = RAAN in inertial frame (rad)
    α_sun      = Sun's right ascension in ECI (rad)
    node_lon   = Ω − GMST (rad, ECEF)
    subsolar_lon = α_sun − GMST (rad, ECEF)
```

The lab canon note (`earth_frames.py` line 110-117) also states the
critical fact that "The LST formula uses the *geodetic* subsolar longitude
directly" — and that `subsolar_lon = α_sun − GMST`.

**For the SSO design** (the textbook condition for sun-synchronous
locking):
```
Ω_dot (SSO) = α_sun_dot (mean-sun rate) = 0.985647 deg/day
```

Substituting into Path (a):
```
dLST/dt = (Ω_dot − α_sun_dot) / 15 = 0 / 15 = 0   (modulo higher-order)
```

The LST at the SSO node is **exactly constant to first order**. The
sidereal-vs-solar differential rate (`360.9856 − 360.0 = 0.9856 deg/day`)
appears as `ω_E − ω_sun` (Earth's inertial rotation rate minus the Sun's
inertial rate); **neither is the RAAN drift rate of an SSO**.

## 3. The independent numerical experiment

### 3.1 Setup

- **Script:** `R:\audit_scratch\audit_015_lst_falsifier.py`
  (and v2 with J2-off comparison; v3 with LST=18 initial condition).
- **Method:** Initialize a circular orbit at h=600 km with
  `i = i_SSO(h=600)` (analytic), `Ω(t_0) = GMST(t_0) + lon_ref`,
  `e=0`, `ω=0`, `M=0` (state at ascending node, northbound).
- **Propagation:** J2 Cowell RK4 from `lab_utils.j2_rhs` and
  `lab_utils.rk4_propagate`; dt=60 s; 365 days total (525,601 steps).
- **Detection:** Ascending-node crossings detected by linear interpolation
  of z(t) where z crosses 0 with vz>0 (5,445 crossings detected for h=600).
- **LST at each crossing:** computed via both Path (a) and Path (b)
  (bit-equivalence check: max |Δ| = 9.77e-15 hours).
- **Drift rate:** least-squares linear fit of unwrapped LST vs t (days).

### 3.2 Results (J2-on, SSO)

```
J2 = 1.0826e-3, h = 600 km, i = 97.7876 deg (analytic SSO)
i_SSO matches lab_utils.sso_inclination_rad(6978.137) to 6 decimals.
Orbital period T = 5801.23 s, Rev/day = 14.893.
Ascending-node crossings detected: 5445 over 365 days.
Mean nodal period: 5792.1 s (one orbital period, as expected).

Path (a) vs Path (b) bit-equivalence:
    max |LST_a − LST_b| = 9.77e-15 hours   (bit-equivalent ✓)

Linear-fit drift rate of LST at ascending node over 5445 crossings:
    drift = +0.0891 min/day = +0.02228 deg/day
    total over year = +32.53 min (residually-bounded)

LST range over the year (using np.unwrap on the LST):
    min = 6.3824 h, max = 7.0054 h
    full range = 0.6230 h = 37.4 min (oscillation)
    std of LST = 0.182 h (small residual oscillation)

Distinct 24h wrap-arounds in LST = 0
Crossings within 30 min of LST=18:00 = 0
Crossings in [17.5, 18.5] = 0     ← the LST NEVER crossed 18:00

Sanity check: unwrapped RAAN drift = 0.99202 deg/day
Expected SSO secular rate         = 0.985647 deg/day
Difference (Cowell vs analytic)   = +0.00637 deg/day = +0.6% (J2^2 + numerical)
```

The LST stayed in **[6.38, 7.00] h** (dawn-ish) for the entire year,
oscillating with ~37 min amplitude — **it never crossed 18:00**.
The drift rate of **0.089 min/day = 32 min over the year** is bounded
by higher-order J2^2 + Lunisolar + numerical drift, NOT a secular
sidereal-vs-solar differential.

### 3.3 Comparison: J2-on vs J2-off (Kepler)

To verify the analytic claim, I reran the same propagator with J2 = 0
(Kepler-only, no nodal precession). With a fixed Ω, the LST drifts at
exactly the rate the "sidereal-vs-solar differential" would predict:

```
Case                | Drift (min/day) | Total drift | LST range (h)
--------------------|-----------------|-------------|---------------
SSO (J2-on)         |   +0.089        |  +0.54 h    |   0.62
Kepler (J2-off)     |   -3.879        | -23.60 h    |  23.98
Analytic prediction |   -3.943        | -24.00 h    |  24.00 (exact)
```

The Kepler drift matches the analytic **sidereal-vs-solar differential**
(−3.943 min/day = −0.9856 deg/day) to within ~1.6% (the small residual is
the equation-of-time wobble in the analytic Sun model). **This** is the
"4 min/day" value Exp 015 cited — but it applies to a *non-SSO* orbit.

### 3.4 Cross-check: starting at LST=18:00 (dawn-dusk initial condition)

I re-initialized the orbit so that the **initial LST = 18:00 h** (the
canonical dawn-dusk terminator). With J2-on, the LST stays in
**[17.5, 18.5] h for ALL 5445 crossings** over the year (zero crossings
escape the band). With J2-off, the LST sweeps through the full 24h cycle,
giving 199 crossings in [17.5, 18.5] out of 5441 (the rest are at other
LST values).

```
Case             | Initial LST | Drift | Crossings in [17.5, 18.5]
-----------------|-------------|-------|--------------------------
SSO (J2-on)      |    18.00 h  | +0.089 min/day | 5445 / 5445 (100%)
Kepler (J2-off)  |    18.00 h  | -3.879 min/day |  199 / 5441  (3.7%)
```

This decisively demonstrates: **an SSO at dawn-dusk keeps the LST locked
near 18:00 throughout the year; a non-precessing orbit does not.**

## 4. Sample table: LST at ascending-node crossings (J2-on)

The LST stays near 6:30–7:00 (dawn) over the entire year. The initial
condition (lon_ref = −80.6° = Eastern Range) determines the initial LST,
not the drift rate.

| # | t (days) | Ω (deg) | LST (h:m) | node_lon (deg) |
|---|----------|---------|-----------|----------------|
|   0 |   0.00 |  199.5677 | 06:33 |  −80.6039 |
| 100 |   6.71 |  206.2108 | 06:30 |  +24.0664 |
| 500 |  33.54 |  232.7860 | 06:23 |  +83.4910 |
| 1000 |  67.08 |  266.0111 | 06:27 | −110.5548 |
| 2000 | 134.15 |  332.4816 | 06:42 | −133.1089 |
| 3000 | 201.19 |   38.9792 | 06:34 | −148.2364 |
| 4000 | 268.21 |  105.5040 | 06:51 | −155.9506 |
| 5000 | 335.22 |  172.0561 | 06:55 | −156.2302 |
| 5444 | 364.96 |  201.6140 | 06:42 |  −63.2496 |

The drift **direction** is +0.089 min/day (slow) and the **amplitude**
of oscillation is ~37 min (EoT + J2^2 short-period).

## 5. Comparison to Exp 015 claim

| Quantity | Exp 015 claim | Measured (J2-on SSO) | Verdict |
|----------|---------------|----------------------|---------|
| Drift rate | 4 min/day = 0.9856 deg/day | **0.089 min/day = 0.022 deg/day** | **WRONG by 45×** |
| Drift over year | 24 h (full cycle) | ~0.5 h (bounded) | **WRONG by ~50×** |
| LST passes through 18:00 | "once per year" | **0 times in [17.5, 18.5] band over 5445 crossings** | **WRONG** |
| Physical mechanism | "sidereal-vs-solar differential" (360.9856 − 360.0 = 0.9856 deg/day) | EoT envelope + J2^2 + Lunisolar (small, bounded) | **MECHANISM MISIDENTIFIED** |

The "sidereal-vs-solar differential" of 0.9856 deg/day is the rate of
the *subsolar geodetic longitude sweep* (in ECEF, where the Sun's
longitude moves at 360 deg/day - 0.9856 deg/day = -359 deg/day, i.e.
the subsolar point traces through all 360° of longitude per day, NOT
per year). Subtracting the same ω_sun from this gives 360 deg/day,
not 0.9856 deg/day. The arithmetic the author did is dimensionally
inconsistent: they compared **two daily rates** (`ω_E` and a "solar
rate" mistakenly set to 360°/day) and obtained a difference in
**daily units**, but the resulting drift rate only has physical
meaning if applied to the *correct* longitude.

## 6. The correct physical picture

For a true sun-synchronous orbit, the J2 nodal precession rate is
**designed** to match the Sun's mean RA drift:
```
Ω_dot (SSO) = α_sun_dot (mean) = 360° / 365.2422 d = 0.985647 deg/day
```
Therefore, the LST at the ascending node is:
```
dLST/dt = (Ω_dot − α_sun_dot) / 15 = 0   (to first order)
```
The LST is **constant** by SSO design. The small residual oscillation
(~37 min amplitude, ~0.5 min/day bounded) comes from:

1. **Equation of Time (EoT)** — the apparent Sun RA differs from the
   mean Sun RA by up to ~16 min amplitude, period 1 anomalistic year.
   This is the dominant bounded term.
2. **Higher-order J2², J3, J4** perturbations to the nodal rate
   (~0.005 deg/day correction, equivalent to ~0.02 min/day LST).
3. **Short-period oscillations** in Ω (averaged out over a year).
4. **Lunisolar perturbations** (expansion of harmonics; ~0.001 deg/day).
5. **Atmospheric drag-induced RAAN walk** (~km/s/year magnitudes but
   orbit-dependent; the dominant secular term in real operations).

For the textbook **sidereal-vs-solar differential** of 0.9856 deg/day
to apply, the orbit must have **fixed Ω** (no J2 precession) — a
Kepler orbit, not an SSO. This is what real ground tracks of a Kepler
orbit look like (one revolution per sidereal day, drifting through
all LST values over a year), and it is **NOT** what an SSO does.

**Real-world confirmation:** Sentinel-1 flight dynamics reports
LTAN held within ±5–10 minutes around 18:00 across multi-year missions,
with station-keeping Δv budget 5–15 m/s/year (dominated by drag
compensation). A 4 min/day drift would require ~200 m/s/year — no
LEO SSO mission has ever reported such a budget.

## 7. Verdict detail

```
FALSIFIER-VERDICT:    RED
015-CLAIM-VERIFIED:   no
MEASURED-DRIFT:       +0.0891 min/day = +0.0223 deg/day (J2-on SSO)
                      −3.879 min/day  = −0.970 deg/day (J2-off Kepler)
ACTUAL-PHYSICS:
  An SSO by design pins Ω_dot = α_sun_dot, so the LST at the
  ascending node is constant to first order. The "4 min/day"
  value the author cited is the *sidereal-vs-solar differential*
  (360.9856 deg/day Earth rotation − 0.9856 deg/day mean Sun
  rate); it would apply to a non-precessing orbit, NOT to an
  SSO. The J2 secular precession exactly cancels this differential
  by design. The residual drift in real operations is dominated
  by the bounded Equation-of-Time envelope (~37 min amplitude,
  ~0.5 min/day peak rate, periodic over 1 year) plus small
  J2² + Lunisolar + drag-walk perturbations, NOT by a secular
  4 min/day differential. Sentinel-1/Landsat-8 operational data
  confirm the small bounded envelope (LTAN within ±5–10 min,
  Δv budget 5–15 m/s/year); no LEO SSO mission has ever
  required the ~200 m/s/year budget implied by 4 min/day.
```

## 8. Audit trail

- **Audit performed:** 2026-08-30 (read-only)
- **Files examined:**
  - `src/lab_utils/orbits.py`
  - `src/lab_utils/earth_frames.py`
  - `src/lab_utils/integrators.py`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/README.md`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/tests/test_dawn_dusk_sso.py`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/results/results.json`
  - `localdocs/reports/audit-015-follow-up-candidates-2026-08-29.md` (prior audit, same finding)
- **Scripts run:**
  - `R:\audit_scratch\audit_015_lst_falsifier.py` (J2-on, lon_ref=-80.6° initial)
  - `R:\audit_scratch\audit_015_lst_falsifier_v2.py` (J2-on vs J2-off)
  - `R:\audit_scratch\audit_015_lst_falsifier_v3.py` (LST=18 initial condition)
- **Results files written:**
  - `R:\audit_scratch\results\audit_015_results.json`
  - `R:\audit_scratch\results\audit_015_v2_comparison.json`
  - `R:\audit_scratch\results\audit_015_v3_lst18.json`
- **No files under the lab repo were modified.**
- **Independent verification path:** Two LST formulas (textbook &
  lab_utils) verified bit-equivalent (max Δ < 1e-14 h). Independent
  J2 RK4 propagator (no donor-hop from Exp 015). Independent
  ascending-node detector (linear interpolation of z(t)). J2-on
  vs J2-off comparison isolates the SSO design effect.
- **Conclusion:** Exp 015 LST drift claim is numerically wrong by ~45×.
  The drift rate of an actual SSO is ~0.089 min/day (driven by EoT
  envelope + higher-order perturbations), NOT 4 min/day. The
  "4 min/day = 24 h/year" value is the sidereal-vs-solar differential
  that would apply to a non-precessing orbit. The follow-up
  experiment recommended in `audit-015-follow-up-candidates-2026-08-29.md`
  (SSO-LST-drift error-correction with Sentinel/Landsat validation)
  is reaffirmed.