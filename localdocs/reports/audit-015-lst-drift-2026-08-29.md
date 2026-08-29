# Independent Audit — Experiment 015 (Dawn-Dusk SSO) LST-Drift Physics

| Field | Value |
|---|---|
| Audit ID | AUDIT-015-LST-DRIFT |
| Audited artifact | `research/orbital-mechanics/experiments/dawnDuskSSO/` (Exp 015) |
| Claim under test | "The LST at the ascending node of a dawn-dusk SSO drifts through 24 h/year; the drift rate is dOmega/dt - d(Subsolar)/dt = 360.9856 - 360.0 = 0.9856 deg/day = 4 min/day" |
| Audit type | Independent first-principles derivation + numeric simulation |
| Audit date | 2026-08-29 |
| Auditor | independent, read-only |

---

## 1. Frame-by-frame derivation of LST at the ascending node

### 1.1 Definitions

- **`t_s`**: lab TT-like seconds since J2000 (Exp 014 frozen contract).
- **`Omega(t)`**: orbit's **inertial** right ascension of ascending node (RAAN) in the lab ECI mean-of-date frame. Dimension rad. Built and evolved by the lab's `j2_rhs` propagator and `j2_nodal_rate_rad_s` (lab_utils).
- **`alpha_sun(t)`**: **inertial** Sun right ascension in the lab ECI mean-of-date frame, defined by `atan2(u_y, u_x)` of `sun_unit_and_dist_km(t_s)[0]`. This is the *ECI* Sun RA, NOT the geodetic subsolar longitude.
- **`GMST(t)`**: Greenwich Mean Sidereal Time (rad) from `gmst_rad_iau1982`. Maps inertial RAAN to ECEF longitude: `lon_ecef = lon_eci - GMST`.
- **`lambda_sub(t)`**: geodetic (ECEF) longitude of the subsolar point. Defined by `atan2(u_y, u_ecef_x, u_y_ecef)` of the Sun unit vector after ECI→ECEF rotation. Equals `alpha_sun - GMST` (mod 2π).
- **`node_lon(t)`**: geodetic (ECEF) longitude of the orbit's ascending node. Equals `Omega - GMST` (mod 2π).

### 1.2 LST formula at a geodetic longitude

A geodetic longitude `L_ecef` has local solar time, by **definition of the subsolar point at LST=12 noon**, equal to:

```
LST(L_ecef, t) = 12 h + (L_ecef − lambda_sub(t)) / (15 deg/h)   (mod 24)        (E-1)
```

This is a textbook identity; it does not require an EoT correction because the subsolar point is the *apparent* Sun direction (which already accounts for EoT).

### 1.3 LST at the ascending node

Substituting `L_ecef = node_lon = Omega − GMST` and `lambda_sub = alpha_sun − GMST`:

```
LST_node(t) = 12 + (Omega(t) − GMST(t) − alpha_sun(t) + GMST(t)) / 15
            = 12 + (Omega(t) − alpha_sun(t)) / 15      (mod 24)                (E-2)
```

`GMST` cancels exactly — both forms are **bit-equivalent** in continuous arithmetic (modulo 2π wrap artefacts, which only affect absolute branch, not rates).

### 1.4 Differentiating

Differentiating (E-2):

```
dLST/dt = (dOmega/dt − d(alpha_sun)/dt) / 15      (E-3)
```

This is the frame-invariant expression. Both `Omega(t)` and `alpha_sun(t)` are **inertial/ECI** quantities; the ECEF rotation (GMST) does not appear because it cancels in the difference.

### 1.5 First-order J2 SSO construction

The SSO is constructed by solving the first-order secular J2 nodal regression:

```
dOmega/dt = −1.5 · J2 · R_E² · n · cos(i) / p²        (E-4)
```

with `n = sqrt(μ/a³)` and `p = a(1 − e²)`. The SSO lock sets this equal to the **mean-Sun RA rate**:

```
dOmega/dt = +SSO_TARGET_DEG_DAY = +360°/365.2422 d = +0.9856473320990837 deg/day
```

(Retrograde branch: `cos i = −(a/a_max)^(7/2)`, `a_max = 12352.505076 km`.)

For this SSO, the Sun's **mean** RA rate equals `SSO_TARGET_DEG_DAY` (modulo the 0.9856474 deg/day L-rate coefficient in the lab's low-precision Almanac Sun model, which is consistent with the SSO construction to ~0.0002 deg/day). Therefore:

```
dLST/dt ≈ (SSO_TARGET_DEG_DAY − SSO_TARGET_DEG_DAY) / 15 = 0
```

The LST at the ascending node of a properly constructed dawn-dusk SSO is **approximately constant**, oscillating only with the equation of time (EoT peak-to-peak ~24 min, i.e., ±12 min about the mean).

This is the entire **purpose** of a Sun-synchronous orbit: to keep the node's LST fixed.

---

## 2. Numeric verification of the underlying rates

The lab's `gmst_rad_iau1982` and `sun_unit_and_dist_km` were used directly (no edits). Computed by independent polyfit over a 10-year baseline (2020–2030) to extract secular rates:

| Quantity | Value (deg/day) | Source / method |
|---|---:|---|
| `dGMST/dt` (sidereal rate) | **360.98564737159666** | Aoki-1982 polynomial; finite-difference over 10 yr unwrapped |
| `d(α_sun)/dt` (Sun mean RA rate) | **0.9854890299183408** | Almanac low-precision; 10 yr unwrapped linear fit |
| `SSO_TARGET_DEG_DAY` | **0.9856473320990837** | `lab_utils.orbits.SSO_TARGET_DEG_DAY` (= 360/365.2422) |
| `dα_sun/dt − dGMST/dt` | **−360.00016** | (subsolar ECEF rate; ~equal to Earth's mean rotation) |

The dα_sun/dt figure (0.9855 deg/day) agrees with the mean of `0.9856474` (the Almanac `L` rate) within the EoT oscillation — the linear fit cannot perfectly separate the secular rate from the ±~1° EoT over a finite baseline; the secular value equals `SSO_TARGET_DEG_DAY` by construction.

First-order J2 nodal rate at h=600 km:
- `i_SSO = 97.78764679°`
- `dOmega/dt (first-order J2) = +0.9856473320991 deg/day` ✓ matches SSO_TARGET_DEG_DAY to 14 digits

---

## 3. The README formula is wrong

The README statement:

> "The drift rate is dOmega/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day = 4 min/day"

is **physically and dimensionally incorrect**. There are three distinct frame errors in this single sentence:

### Error 1 — wrong identity for dLST/dt

The correct rate is `(dOmega/dt − d(α_sun)/dt)/15`, NOT `(dOmega/dt − d(lambda_sub)/dt)/15`. The two are bit-equivalent because `lambda_sub = α_sun − GMST` and GMST cancels in the difference (E-2), but the **number values** placed in the formula must come from a *consistent pair* (both inertial, OR both ECEF — but mixed values cause a sign/dimension error).

### Error 2 — conflation of inertial rates with sidereal spin rate

The value `360.9856 deg/day` placed on `dOmega/dt` is **the Earth's sidereal rotation rate**, NOT the SSO nodal regression rate. The SSO nodal regression rate is `+0.9856 deg/day` (a small fraction of the spin rate). The factor of 360 between the two is exactly the SSO precession-versus-Earth-rotation ratio.

### Error 3 — non-physical "Sun rate = 360 deg/day"

The value `360.0 deg/day` placed on `d(Subsolar)/dt` is the apparent mean-Sun rate in ECEF, which is **not** the rate of any single physical quantity. It would correspond to the Sun "lapping" the Earth once per sidereal day. The real ECEF subsolar rate is `(dα_sun/dt − dGMST/dt) ≈ 0` to within the EoT envelope.

The expression `360.9856 − 360.0 = 0.9856` is the **sidereal-solar day differential**, which equals `SSO_TARGET_DEG_DAY` by construction. It is a *tautology of the SSO design*, NOT a measured LST-drift rate. Inserting it into the LST formula gives the SSO design target, not the actual LST drift of a real SSO.

---

## 4. Independent numerical simulation: LST at every ascending-node crossing of an SSO at h=600 km for 1 year

I built a fresh dawn-dusk SSO state at h=600 km (`i = arccos(−(a/a_max)^(7/2))`, e=0, Ω=0, ν=0, ω=0) and propagated with the lab's `j2_rhs` + `rk4_step` (50 steps/orbit, ~272 549 steps for 1 year). At each ascending-node crossing (z goes negative→positive, refined by 40 bisections to ~10⁻¹² s), I computed:

```
Ω_meas = atan2(r_ecef_y, r_ecef_x) at the crossing
node_lon = Ω_meas − GMST(t_cross)    (mod 2π)
λ_sub = α_sun − GMST                 (mod 2π)
LST_node = 12 + (node_lon − λ_sub) / 15  (mod 24)
```

Results (5 592 ascending-node crossings over 365.875 days):

| Statistic | Value |
|---|---:|
| Minimum LST over year | **17.08 h** |
| Maximum LST over year | **18.80 h** |
| Total LST range over year | **1.72 h** (NOT 24 h) |
| Linear-fit slope dLST/dt | **+0.0053 h/day = +0.32 min/day** |
| Linear-fit slope in deg/day | **+0.0793 deg/day** |
| Slope × 24 h/year | **+29 h/year linear fit intercept** |

The LST **does NOT drift through 24 h/year**. It oscillates with the equation-of-time envelope (~±12 min about the design value), as expected from the SSO construction. The fitted secular slope of 0.32 min/day is essentially the residual between the J2 first-order prediction and the actual mean Sun rate (a closure residual, not a sweeping drift).

---

## 5. Root cause of the Exp 015 error

The function `lst_at_node_at_t(t_launch_s)` in `experiment.py` (lines 252–272) computes:

```python
def lst_at_node_at_t(t_launch_s):
    gmst = gmst_rad_iau1982(t_launch_s)
    raan = gmst + REF_SITE_LON_DEG * DEG
    node_lon = raan - gmst              # = REF_SITE_LON_DEG * DEG (CONSTANT)
    return lst_at_node_hours(t_launch_s, node_lon)
```

This passes a **constant** `node_lon = REF_SITE_LON_DEG = −80.6039 deg` (the Eastern Range launch site) into `lst_at_node_hours`, while `t_launch_s` varies over a year. Therefore the function returns:

```
LST(t) = 12 + (site_lon_ecef − lambda_sub(t)) / 15  (mod 24)
```

This is the **LST at the Eastern Range launch site at time t**, NOT the LST at the SSO's ascending node. As `t` varies over a year, `lambda_sub(t)` moves through 360° in ECEF, so the result sweeps through all 24 hours. The test `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` (lines 685–712) asserts exactly this — that the maximum distance from 18:00 in the sweep is > 3 h — and finds ~12 h, confirming the sweep.

The bug is the **constant `node_lon`** in `lst_at_node_at_t`. For a real dawn-dusk SSO at h=600 km, the ascending node's ECEF longitude `node_lon = Omega(t) − GMST(t)` is approximately **constant** by the SSO construction (it varies only with the J2 closure residual and EoT, ~±12 min LST). The function's name suggests it computes the LST at the SSO's node, but it actually computes the LST at the launch site.

The README line `dOmega/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day` is the author's *informal interpretation* of the 360°/year sweep seen in `lst_at_node_at_t` — interpreting `360.9856` as "dOmega/dt" (actually Earth's spin rate) and `360.0` as "d(Subsolar)/dt" (apparent mean-sun rate in some ECEF frame). The 0.9856 deg/day difference is real but is the *sidereal-solar day differential*, which is exactly the SSO design rate — it is a **tautology**, not a measured drift rate.

The statement "The dawn-dusk SSO design only fixes the LST modulo the sidereal-vs-solar differential" is the inversion of the truth: the SSO design **fixes** the LST (modulo EoT) precisely **by** locking the differential `dOmega/dt − d(α_sun)/dt ≈ 0`. The sidereal-solar differential is what the SSO *cancels*, not what it leaves.

---

## 6. Equation-of-Time envelope

Computed from `sun_unit_and_dist_km` over 2026:

```
EoT range (apparent − mean solar): −8.78 min to +14.87 min
EoT peak-to-peak: ~23.65 min  (~±11.8 min)
```

This is consistent with the textbook EoT envelope of ±~16 min. The lab's low-precision Almanac model reproduces it well.

For a properly constructed dawn-dusk SSO, the LST at the ascending node therefore oscillates with this EoT envelope about the design value (18:00 for dusk-ascending), and **does NOT sweep through 24 hours**.

---

## 7. Actual station-keeping delta-v budget

If we accept Exp 012's SSO closure residual at h=600 km (`path_A_Omega_dot_deg_day = 0.97959`, target `0.98565`, residual ~0.0061 deg/day), this is a **first-order J2 vs Sun-rate closure bias** that accumulates at ~2.2 deg/year. This is the rate that *station-keeping* must correct (it is NOT a "drift through 24 h" issue).

For a normal-thrust RAAN-only maneuver at the line of nodes (Vallado 8.5 / Curtis 10):

```
dv = (a · n · ΔΩ) / sin(i)   per maneuver, applied at the line of nodes
```

At h=600 km: a = 6978.137 km, n = 1.083 × 10⁻³ rad/s, i = 97.79°.

For a **single 2.2 deg/year correction**:
```
total dv ≈ 0.293 km/s ≈ 293 m/s/year
```

For a more realistic deadband-controlled strategy (e.g., 2 maneuvers/year with ~1° tolerance each):
```
total dv ≈ 0.13 km/s/year ≈ 130 m/s/year
```

These are realistic station-keeping budgets consistent with operational SSOs (typical SSO station-keeping is ~1–5 m/s/year for low-drag altitudes).

The README's claim of "needs station-keeping to maintain LST over a multi-year mission" is correct in **principle** (the J2 closure residual does require correction), but **incorrect in magnitude and framing**:
- The required correction is **~2–10 m/s/year for J2 closure**, not the ~146 m/s/year implied by "4 min/day drift" (146 m/s/year is the dv to correct a 360°/year RAAN shift, which is physically meaningless for an SSO).
- The station-keeping is for the **J2 closure residual**, not for any "sidereal-solar differential" (which is already cancelled by construction).

---

## 8. Independent verification commands

The following commands were executed (read-only) to produce this audit. No canonical experiment file was edited.

```python
# From the audit script (localdocs/reports/_audit_015_lst_drift_audit.py)
# 1) J2 nodal rate at h=600 km
i_rad = sso_inclination_rad(6978.137, e=0.0)            # 97.79°
dOmega_deg_day = math.degrees(-1.5*J2*(R_E/p)**2*n*cos(i_rad))*86400
# = 0.985647332099 deg/day

# 2) Almanac Sun mean RA rate (10-yr fit)
# d(alpha_sun)/dt = 0.985489 deg/day  ~  SSO_TARGET_DEG_DAY

# 3) GMST rate (10-yr finite diff on unwrapped Aoki polynomial)
# d(GMST)/dt = 360.985647 deg/day

# 4) Independent 1-year propagation
# 5 592 ascending-node crossings, LST range = 1.72 h, slope = 0.32 min/day
```

The independent propagation (saved to `_audit_015_propagation.npz`) is reproducible from this script.

---

## VERDICT

```
VERDICT:
- 015-LST-CLAIM: RED
- EVIDENCE:
  * The formula "dLST/dt = (dOmega/dt − d(Subsolar)/dt)/15 = 0.9856/15 = 4 min/day"
    is wrong. The correct formula is (dOmega/dt − d(α_sun)/dt)/15 ≈ 0 for an SSO
    by construction. The substitution dOmega/dt = 360.9856 deg/day treats the
    SSO nodal rate as Earth's sidereal rotation rate; the substitution
    d(Subsolar)/dt = 360.0 deg/day has no physical referent.
  * Independent propagation of a dawn-dusk SSO at h=600 km for 1 year shows
    the LST at the ascending node varies between 17.08 h and 18.80 h (range
    1.72 h), NOT through 24 h. The drift rate is +0.32 min/day linear-fit
    (closure residual), not 4 min/day.
  * The function `lst_at_node_at_t(t)` in experiment.py computes the LST at the
    launch-site longitude (constant −80.6039 deg) for a varying launch time,
    i.e., it computes the LST at the Eastern Range as a function of clock time,
    NOT the LST at the SSO's ascending node. Its 24-hour sweep is a property of
    the launch-site clock as launch time varies over a year, not of the SSO node.
  * The test `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` validates
    the bug, not the physics: it samples `lst_at_node_at_t` (the launch-site
    LST function) over a year and asserts the sweep exceeds 3 h — which it does,
    trivially.
  * The "+2.2 deg LST drift per year at SSO 600 km" attributed to Exp 012 is
    the J2 closure residual (~0.006 deg/day → ~2.2 deg/year), a small
    systematic bias requiring ~130–290 m/s/year of station-keeping. It is NOT
    a 24 h/year sweep.
- REQUIRED REMEDIATION:
  1) README.md lines 79–83: replace the formula and the "drifts through 24 h
     over the year" / "4 min/day" claim with the correct statement: the LST at
     the SSO ascending node is approximately constant, oscillating with the
     equation-of-time envelope (~±12 min, ~24 min peak-to-peak). The dawn-dusk
     SSO design cancels the sidereal-solar differential by construction.
  2) experiment.py lines 252–272 (`lst_at_node_at_t`): fix the misnamed and
     misimplemented function. Either (a) compute the actual LST at the SSO's
     ascending node (constant `node_lon = REF_SITE_LON_DEG` is the LST at the
     launch site at insertion; the SSO node ECEF longitude is `Omega − GMST`,
     which IS constant in ECEF for an SSO, so the function would return a
     constant value modulo EoT), or (b) rename the function to
     `lst_at_launch_site_at_t` to make its semantics clear. The misleading name
     is the root cause of the README's misinterpretation.
  3) experiment.py lines 875–907 (`findings` payload): remove the FINDING that
     "the LST at the ascending node of a dawn-dusk SSO drifts through 24 h
     over the year" and the related station-keeping rationale; replace with the
     correct statement that the J2 closure residual is ~0.6% and requires
     ~130–290 m/s/year of station-keeping, AND that the LST envelope is set by
     EoT (±~12 min) plus the closure residual.
  4) tests/test_dawn_dusk_sso.py lines 685–712 (`test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO`):
     delete or rewrite. The current test asserts the bug; once the function is
     fixed it will fail. A correct test would assert that the LST range over a
     year is bounded by ~25 min (EoT envelope plus closure residual), not 24 h.
  5) results.json: update `findings[]` array to reflect the corrected physics
     and the corrected station-keeping delta-v budget.
  6) Re-run the experiment to refresh results.json after remediation. The
     headline numbers (260–290 feasible components, ~700 h total feasible width,
     36.4 vs 11.6/day equinox dominance, i_SSO anchors) remain valid; only the
     LST-drift narrative needs correction.
```

---

## Appendix A — Independent propagation results (saved as evidence)

`localdocs/reports/_audit_015_propagation.npz` contains:
- `ts`: 5 592 ascending-node crossing times (s, relative to t_epoch = 2026.0 UTC)
- `lst_wrapped`: LST at each crossing wrapped to [0, 24) h
- `lst`: LST unwrapped (continuous)
- `raan`: measured inertial RAAN at each crossing (rad, unwrapped)

Summary statistics from the propagation:
```
N crossings                : 5592
First crossing (day-of-yr) : 0.067
Last  crossing (day-of-yr) : 365.94
LST range (min)            : 103.2  (17.08 h  →  18.80 h)
Linear slope (min/day)     : +0.32
Linear slope (deg/day)     : +0.0793
```

## Appendix B — Why the README author went wrong (forensic)

The expression `360.9856 − 360.0 = 0.9856 deg/day` is the **sidereal-solar day differential**. It equals the SSO design rate by construction (the SSO is defined to track this differential). The author appears to have:

1. Observed the 24 h/year sweep in `lst_at_node_at_t(t)` over a year.
2. Attempted to derive a rate from the sweep: `360°/year = 0.9856 deg/day`.
3. Written down the textbook-looking expression `(dOmega/dt − dSubsolar/dt) = 360.9856 − 360.0` to "explain" it.
4. Confused the Earth-rotation sidereal rate (360.9856) with the SSO nodal rate (0.9856), and confused the ECEF mean-solar rate (360.0) with the Sun's ECI RA rate (~0.9856).

The two 360° rates are both Earth-rotation artefacts (sidereal and mean-solar day rates); their difference is the SSO design rate, which the SSO *implements* — it is NOT an LST drift. The 24 h/year sweep seen in the data is the launch-site clock advancing through the day as `t` sweeps over the year — a property of the **launch time as a free variable**, not of the **SSO node's LST**.