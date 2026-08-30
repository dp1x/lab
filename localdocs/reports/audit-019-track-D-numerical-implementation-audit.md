# Audit-019 Track D: Numerical-Implementation Audit

**Date:** 2026-08-30
**Scope:** Independently audit the numerical third-body implementation, center-of-mass / differential acceleration, Sun/Moon frame transformations, ephemeris interpolation, and epoch / time handling for Experiment 019 (Lunisolar Long-Period Terms and Secular-Limit Convergence), as built on Exp 017 (`lunisolarVerification`) and Exp 018 (`lunisolarReconciliation`).
**Method:** Read-only inspection of source code, manifests, snapshots, and lab utilities. No modification of source. All quantitative claims independently derived in this track; cited file:line locations are taken verbatim from the audited code.

---

## 1. Third-body acceleration verification

### 1.1 Derivation of the expected form (geocentric ECI)

For a satellite at geocentric position `r_sat`, attracted by a third body (Sun or Moon) at geocentric position `r_3` (both measured relative to Earth's center), the inertial-frame equation of motion in a geocentric, non-rotating frame is:

    a_sat = -mu_E * r_sat / |r_sat|^3                    (Kepler toward Earth)
           + mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3     (direct attraction of satellite by r_3)
           - mu_3 * r_3 / |r_3|^3                        (indirect: Earth is also pulled toward r_3)

The third term (the "indirect term") is the acceleration of Earth's center of mass by the third body, subtracted because we are working in a non-inertial geocentric frame. In barycentric coordinates this term would not appear; in geocentric ECI it appears with a minus sign because the geocentric frame is accelerated relative to the inertial frame by the third body's pull on Earth.

Equivalent form (the same physics, expressed differently):

    a_sat = -mu_E * r_sat / |r_sat|^3
           + grad_{r_sat} [mu_3 / |r_3 - r_sat|]         (direct gradient)
           + mu_3 * r_3 / |r_3|^3                        (indirect, +sign here because of how the disturbing potential U_3 = mu_3 (1/|r_3 - r_sat| - r_sat . r_3 / |r_3|^3) is written)

Both forms are equivalent and yield the same numerical result.

### 1.2 017 implementation audit

**Location:** `lunisolarVerification/experiment.py:418-441` (`_lunisolar_third_body_accel`).

The implementation reads (lines 428-439):

```python
r_sat_to_sun = r_sun_eci - r_eci_km
r3_sun = np.linalg.norm(r_sun_eci)
r3s_sun = np.linalg.norm(r_sat_to_sun)
a_sun = SOLAR_GM_KM3_S2 * (
    r_sat_to_sun / r3s_sun**3 - r_sun_eci / r3_sun**3
)
```

with the same pattern for the Moon on lines 433-439.

**Verdict:** This matches the expected form `mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3 - mu_3 * r_3 / |r_3|^3` exactly. Signs are correct: `r_sat_to_sun = r_sun - r_sat` (positive direct term, attraction toward the third body); the indirect term is subtracted (negative). The Moon term is identical in structure.

### 1.3 018 implementation audit

**Location:** `lunisolarReconciliation/experiment.py:269-289` (`_third_body_accel`).

The implementation reads (lines 279-289):

```python
r_sun = _interp_snapshot_precessed(t_s, sun_snap, apply_precession)
r_sat_to_sun = r_sun - r_eci_km
r3 = np.linalg.norm(r_sun)
r3s = np.linalg.norm(r_sat_to_sun)
a_total += SOLAR_GM_KM3_S2 * (r_sat_to_sun / r3s ** 3 - r_sun / r3 ** 3)
```

with the same pattern for the Moon.

**Verdict:** Identical algebraic structure to 017. Sign conventions correct. The only delta between 017 and 018 in this function is the `apply_precession` option, which rotates the Sun/Moon geocentric vectors before using them in the acceleration formula (see §2 for the precession audit; the bug there is in the rotation, not in this acceleration form).

### 1.4 Magnitude of the indirect term (why it is included even though small)

At LEO with `|r_sat| ~ 6978 km`, `|r_3| ~ 1.5e8 km` (Sun) or `~ 3.8e5 km` (Moon), the indirect/direct magnitude ratio is of order `|r_sat| / |r_3|`:

- Sun:  `|r_sat| / |r_3| ~ 6978 / 1.5e8 ~ 4.7e-5` (well below 1e-4; the 017 comment "order 1e-5 of the direct term" is correct in this bound)
- Moon: `|r_sat| / |r_3| ~ 6978 / 3.8e5 ~ 1.8e-2` (~2% — NOT 1e-5 as the 017 comment implies; this is a one-order-of-magnitude understatement in the comment, although for RAAN-rate purposes the indirect Moon term still contributes only at the few-percent level relative to the direct Moon term and is dwarfed by the short-period residuals that dominate the 1-year fit)

**Audit finding (cosmetic):** The 017 docstring at `lunisolarVerification/experiment.py:421` says the indirect term is "order 1e-5 of the direct term". This is the bound for the **Sun**. For the Moon the ratio is closer to 1.8e-2 (~2%). Including the indirect term is correct for both, and the bound quoted in the comment is the worst case (Sun). The Moon indirect term is genuinely larger than the comment suggests (~350x larger than the quoted bound) but is still small relative to the direct Moon term and is captured correctly. This is a documentation imprecision, not a code bug.

**Why the indirect term is still kept:** The differential geocentric formulation requires the indirect term to be self-consistent with the underlying physical model (the satellite and the Earth are both accelerated by the third body; in geocentric coordinates the Earth's acceleration must be subtracted to give the satellite's motion in the geocentric frame). Omitting it would bias the secular RAAN rate by ~2% of the Moon contribution, which is comparable to or larger than the corrected secular rate itself (see §8 for the quantified impact).

### 1.5 Force-level identity check (018 experiment.py:297-348)

The 018 implementation includes a `force_level_identity_check` (50 random states, fixed seed) that verifies the direct + indirect form equals an independently derived algebraic equivalent (lines 318-340) to machine precision. The reference implementation is at line 328-331:

```python
a_sun_b = -SOLAR_GM_KM3_S2 * (r_sat - r3_sun) / r3s_sun ** 3 - SOLAR_GM_KM3_S2 * r3_sun / r3_mag_sun ** 3
```

This is the alternative form `a = mu_3 * (r_sat - r_3) / |r_sat - r_3|^3 - mu_3 * r_3 / |r_3|^3` with the direct term sign flipped (since `r_sat - r_3 = -(r_3 - r_sat)`), which is the same physics written differently. The two forms must agree to machine precision; the test passes with `max_diff < 1e-15 km/s^2` per the published finding at `lunisolarReconciliation/experiment.py:639-642`.

**Verdict:** The acceleration implementation is correct. The force-level identity check is a strong self-consistency validator. **No bugs identified in the acceleration formula itself.**

---

## 2. Frame transformation audit (IAU-1976 precession)

### 2.1 017 vs 018 precession handling

**017 (lunisolarVerification):** Does NOT apply precession. The Sun and Moon vectors are used as-is (in ICRF/J2000), and the propagator treats the J2 axis as the lab's "ECI mean-of-date" frame. This is the frame mismatch identified by Track D in audit-018. At 2026, the ICRF-vs-MOD offset of the Sun and Moon geocentric directions is ~0.4 deg.

**018 (lunisolarReconciliation):** Applies the IAU-1976 precession rotation via `precession_j2000_to_mod(t_s)` at lines 147-156 to the Sun and Moon vectors before the third-body acceleration is computed (line 259 inside `_interp_snapshot_precessed`). This rotation is applied with `apply_precession=True` for all controlled experiments (`run_force_isolation`, `run_inclination_sweep`, `run_window_sensitivity`, `run_precession_comparison` with `apply=True`). Only `run_precession_comparison` with `apply=False` (line 537-538) bypasses precession to isolate the frame-mismatch bias.

### 2.2 Lieske polynomial coefficients verification

**Location:** `lunisolarReconciliation/experiment.py:147-156`.

```python
def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    T = t_s / (86400.0 * 36525.0)
    sec = math.radians(1.0 / 3600.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * sec
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * sec
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * sec
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)
```

The Lieske 1977 polynomial coefficients are:

| Coefficient | 018 uses | Standard (Lieske et al. 1977) |
|-------------|----------|------------------------------|
| zeta_A linear (arcsec/cy) | 2306.2181 | 2306.2181 ✓ |
| zeta_A quadratic | 0.30188 | 0.30188 ✓ |
| zeta_A cubic | 0.017998 | 0.017998 ✓ |
| z_A linear | 2306.2181 | 2306.2181 ✓ |
| z_A quadratic | 1.09468 | 1.09468 ✓ |
| z_A cubic | 0.018203 | 0.018203 ✓ |
| theta_A linear | 2004.3109 | 2004.3109 ✓ |
| theta_A quadratic | -0.42665 | -0.42665 ✓ |
| theta_A cubic | -0.041833 | -0.041833 ✓ |

**Verdict:** Polynomial coefficients match the IAU-1976 standard (Lieske, Lederle, Fricke 1977; see Wikipedia "Axial precession" / Bibliography, Lieske et al. 1977). Time convention `T = t_s / (86400 * 36525)` is the standard Julian-century conversion since J2000. **Coefficients are correct.**

### 2.3 Sign convention audit (the BUG)

**018 `_rot3` at line 132-135:**

```python
def _rot3(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
```

**eclipseTiming `_rot3` at `eclipseTiming/experiment.py:255-258` (the reference implementation):**

```python
def _rot3(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
```

The two matrices are **transposes** of each other (sign-flipped off-diagonals). Mathematically:
- `ecl._rot3(θ)` is the standard active rotation about +z by +θ, i.e., `R_z(+θ)` in active-vector convention.
- `018._rot3(θ)` is `R_z(-θ)` in active-vector convention.

Independent verification by computing both implementations at T=0.26 centuries (epoch 2026):
- 018 rotates X-axis by **+0.333 deg**
- eclipseTiming rotates X-axis by **-0.333 deg**
- Both are near-identity (det = 1, components within 0.005 of identity)
- 018's matrix is approximately the inverse-transpose of the standard Lieske J2000→MOD rotation

**Implication:** The precession rotation magnitudes in 018 are correct (Lieske polynomial values are applied correctly), but the **direction of rotation about the Z-axis is reversed**. Since the byte-pinned ICRF/J2000 Sun and Moon vectors are rotated by 018's matrix with the intent of "mapping J2000 → MOD", but the actual matrix maps J2000 → "MOD-inverse" (the time-reverse of MOD), the frame error in 018 is **NOT reduced** relative to 017; it is **reversed in direction**, leaving a residual ~0.66 deg (twice the ICRF/MOD difference) instead of 0 deg.

**Magnitudes (T = 0.26 centuries, January 2026):**
- zeta_A ≈ 0.1666 deg
- z_A ≈ 0.1666 deg
- theta_A ≈ 0.1448 deg
- Net rotation magnitude (dominant component about Z): ~0.33 deg (signed)

The eclipseTiming reference at `eclipseTiming/experiment.py:266-277` (the docstring claims "P = R3(-z) R2(theta) R3(-zeta)") is consistent with the standard Lieske J2000→MOD form: `R3(-zeta) R2(theta) R3(-z)` is the transpose of `R3(zeta) R2(-theta) R3(z)` (MOD→J2000), so taking the inverse gives the standard J2000→MOD rotation as `R3(-zeta) R2(theta) R3(-z)`. The 018 code computes the inverse of this.

**Defect classification:** Sign error in `_rot3` (or equivalently, sign error in the ordering of multiplications). The Lieske polynomial values are applied correctly to the angles; the matrix construction uses the wrong cos/sin pattern.

### 2.4 Impact on the secular RAAN rate

The 0.66 deg frame error (vs the 0.5 deg error in 017) introduces a fractional error in the third-body unit vector of ~0.66/57.3 = 1.15e-2.

For the Sun at h=600 km with `i_sso = 97.79 deg`, the secular RAAN contribution from the Sun is `(3/8) n (mu_S/mu_E) (a/AU)^3 sin(2(i - i_3)) / sin(i)`. The fractional change from a frame error δ in the Sun's apparent inclination is:

    d ln(rate) / d i_3 = -2 cos(2(i - i_3)) / sin(2(i - i_3))

At `i - i_3 = 97.79 - 23.44 = 74.35 deg`:
    `2(i - i_3) = 148.7 deg`
    `cos(148.7) = -0.854`
    `sin(148.7) = 0.520`
    `factor = -2 * (-0.854) / 0.520 = +3.28`

So a 0.66 deg shift in `i_3` produces `3.28 * 0.66/57.3 = 3.78e-2` fractional change in `sin(2(i-i_3))`. Multiplied by the ~6.7e-5 deg/day solar secular rate gives ~2.5e-6 deg/day solar bias. Moon is similar in magnitude (the Moon's secular rate is ~1.5e-4 deg/day, also prograde). Combined solar + lunar RAAN rate bias from the frame mismatch: **~7e-6 deg/day** = **~2.5e-3 deg/year** (PROGRADE).

This is **200x smaller** than the 0.5 deg/year frame-fix that the 018 precession rotation claims to deliver (line 27 of the 018 experiment.py docstring: "0.4 deg frame mismatch produces a ~0.5 deg/year bias"). The actual residual frame bias in 018 is **NOT zero** as the corrected implementation intends — it is ~0.0025 deg/year (prograde), which is below the 1-year linear-fit residual RMS in 017 (~0.4 deg; see 017 results.json `linear_fit_residual_rms_deg`) but above the 018 corrected-formula-vs-numerical agreement (1.35e-4 deg/day = 0.05 deg/year).

**Concrete impact on 018's headline numbers:**
- Corrected formula at h=600 km: +1.35e-4 deg/day (Track B derivation, no frame-error contribution)
- 1-year numerical fit: +1.28e-3 deg/day (prograde)
- Residual: 9.5x ratio, attributed to "unmodelled short-period terms" (evection + variation + lunar nodal)
- Bug contribution: +2.5e-3 deg/year ≈ +7e-6 deg/day, **3% of the corrected formula's magnitude**, much smaller than the 9.5x short-period residual
- The 9.78x residual at i_sso in the 018 headline is dominated by short-period physics (evection ~27.55 d, variation ~14.77 d, lunar nodal 18.6 yr); the sign-bug in the precession rotation adds a small ~3% bias on top of this and does **NOT** invalidate the 018 headline finding that the corrected formula agrees with the numerical in sign and within 10x in magnitude.

### 2.5 The 017 precession on/off comparison (Track D signal)

The 018 experiment runs both `apply_precession=True` (default) and `apply_precession=False` (017-equivalent) at h=600 km in `run_precession_comparison` (lines 530-540). The expected difference between the two configurations is:
- 017-equivalent (no precession): ICRF/MOD offset of ~0.5 deg → bias of ~7e-6 deg/day * 0.5/0.66 = ~5e-6 deg/day on the RAAN rate
- 018 (wrong-sign precession): MOD-reverse offset of ~0.66 deg → bias of ~7e-6 deg/day on the RAAN rate

The DIFFERENCE between these two configurations is `7e-6 - 5e-6 = 2e-6 deg/day = 7e-4 deg/year`, which is **below the 0.012 deg/year residual** that the 018 docstring estimates for the Track D finding. The experiment's `precession_comparison` will therefore measure a smaller bias than the docstring claim, but the bias is in the right direction (the wrong-sign precession makes things slightly worse than no precession, by ~0.4e-3 deg/year).

**Verdict:** The 018 implementation has a **sign bug in the precession rotation** (`_rot3` uses the transpose of the standard form). The magnitude of the resulting frame error is ~0.66 deg at 2026 (vs ~0.5 deg in 017 with no precession, and 0 deg with the correct-sign precession). The impact on the secular RAAN rate is ~2.5e-3 deg/year, which is small enough that the 018 headline finding (corrected formula vs 1-year numerical, 9.78x ratio dominated by short-period physics) is not invalidated by this bug, but the 018 contract claim that precession fixes the Track D frame mismatch is incorrect: the sign bug **preserves** a frame error similar in magnitude to the original 017 problem.

---

## 3. Sun/Moon snapshot integrity

### 3.1 SHA-256 verification

Verified independently via `Get-FileHash -Algorithm SHA256`:

| File | Computed | Manifest |
|------|----------|----------|
| `eclipseTiming/reference/horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt` | `06d54fb35523a0af6ba3ea738315f1e3f5b996067c40f474052cd2fb5b5658ec` | `06d54fb35523a0af6ba3ea738315f1e3f5b996067c40f474052cd2fb5b5658ec` ✓ |
| `lunisolarVerification/reference/horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt` | `65f1d67f798a3b95bb87310efae3200027098869246567a68ccd671d79978f4a` | `65f1d67f798a3b95bb87310efae3200027098869246567a68ccd671d79978f4a` ✓ |

Both snapshots match their manifests byte-for-byte.

### 3.2 SOE/EOE parsing

017 and 018 use identical parsing (017 at `lunisolarVerification/experiment.py:362-381`, 018 at `lunisolarReconciliation/experiment.py:209-228`):

```python
soe_idx = None
eoe_idx = None
for i, line in enumerate(lines):
    if "$$SOE" in line:
        soe_idx = i + 1
    if "$$EOE" in line:
        eoe_idx = i
        break
rows = []
for line in lines[soe_idx:eoe_idx]:
    s = line.strip()
    if not s:
        continue
    parts = [p.strip() for p in s.split(",")]
    jd_tt = float(parts[0])
    x = float(parts[2])
    y = float(parts[3])
    z = float(parts[4])
    rows.append((jd_tt, x, y, z))
```

Verified independently:
- **Sun snapshot:** 366 rows between $$SOE and $$EOE (matches MANIFEST.json `validation.rows: 366`)
- **Moon snapshot:** 366 rows between $$SOE and $$EOE (matches MANIFEST.json `validation.rows: 366`)
- All rows are non-empty (no blank-line skipping triggered)
- JD values are monotonic increasing
- Step sizes are uniform at exactly 86400.0 s (1-day cadence) — first 5 step sizes [86400, 86400, 86400, 86400, 86400], last 5 step sizes [86400, 86400, 86400, 86400, 86400]

**Verdict:** Parsing is correct, extracts all 366 daily rows, cadence is uniform.

### 3.3 Time base verification

**Sun snapshot first row:** `JD = 2461041.5`, last row: `JD = 2461406.5`. Per Horizons echo in the snapshot header, this is "A.D. 2026-Jan-01 00:00:00.0000 TDB" to "A.D. 2027-Jan-01 00:00:00.0000 TDB".

**Time conversion (017 line 384, 018 line 230):**

```python
t_s = (arr[:, 0] - JD_J2000) * 86400.0
```

with `JD_J2000 = 2451545.0` (from `lab_utils/earth_frames.py:67`).

Verified:
- Snapshot first row: `t_s = (2461041.5 - 2451545.0) * 86400.0 = 820497600.0 s`
- Lab convention `t0 = 820540800.0 s` (used in both 017 line 555 and 018 line 449)
- `t0 - 820497600.0 = 43200 s = 12 hours`

**Finding:** The snapshot's first row corresponds to 2026-01-01 00:00 TDB; the lab's `t0` corresponds to JD 2461042.0 = 2026-01-01 12:00 TT (per lab convention `lab_utils/earth_frames.py:67` and the 018 comment at line 449). The lab's `t0` sits at fractional index 0.5 between snapshot rows 0 and 1.

**Propagation window:**
- Lab propagator runs from `t0 = 820540800` to `t_end = t0 + 365 * 86400 = 852076800`
- Snapshot covers `t_s = [820497600, 852033600]`
- 018 propagator window: `[820540800, 852076800]` (lab)
- Snapshot window: `[820497600, 852033600]` (Horizons)
- Delta at start: lab is **+43200 s into the snapshot** (snapshot extends 12 hours before lab_t0)
- Delta at end: lab propagator runs **+43200 s past the snapshot's last row**

The 018 `_interp_snapshot_precessed` (line 247-249) clamps to endpoint values when the query time exceeds the snapshot range:

```python
if t_query_s <= t_s[0]:
    rv = r[0]
elif t_query_s >= t_s[-1]:
    rv = r[-1]
```

So the final 12 hours of the propagation use the clamped final-row values for Sun/Moon positions. This is a small but real error source — see §4 for the order of magnitude. The first 12 hours do NOT trigger the clamp because the snapshot starts 12 hours BEFORE t0 (the snapshot's first row at `t_s = 820497600` is earlier than `t0 = 820540800`).

**Subtle finding (sign reversal of the clamp):** The 018 code's clamp logic is:

```python
if t_query_s <= t_s[0]:
    rv = r[0]
```

This clamps to the FIRST snapshot row when `t_query_s` is BEFORE the snapshot. Since the snapshot starts at 2026-01-01 00:00 TDB and the propagator starts at 2026-01-01 12:00 TT (lab convention), the propagator's first query time is INSIDE the snapshot range — no problem here. But if the propagator were ever run with a `t0` BEFORE the snapshot's first row, the clamp would silently return the first row's value (NOT raise an error). The 017 clamp (line 396-399) is identical.

This is consistent with the 017 docstring (line 386-388): "Outside the snapshot range, clamp to the endpoint values (disclosed as a known limitation - the experiment runs strictly within [t0, t0 + 1 yr] where the snapshot fully covers)." For the actual experiment, this is correct. For the first 12 hours the snapshot extends BEFORE lab_t0; for the last 12 hours the propagator extends PAST the snapshot end. The "clamp" behavior is triggered ONLY for the last 12 hours, where the lab propagator queries `t > 852033600` and gets the final snapshot row clamped.

**Verdict:** Time base conversion is correct (`(JD - JD_J2000) * 86400` matches the standard TT/TDB convention at J2000 noon = 0). The 12-hour offset between the snapshot's first row and the lab's `t0` is internally consistent because the snapshot covers [lab_t0 - 12 h, lab_t0 + 1 yr - 12 h]. The last 12 hours of the 1-year arc use clamped endpoint values — see §4 for the impact.

### 3.4 Distance-band verification

Verified at three sample dates (k=0, k=183, k=365):

| Date | Sun \|r\| (km) | Moon \|r\| (km) |
|------|----------------|------------------|
| 2026-01-01 00:00 TDB | 147,103,578 | 361,026 |
| 2026-07-02 00:00 TDB | 152,083,483 | 399,054 |
| 2026-12-31 00:00 TDB | 147,106,071 | 391,495 |

Physical band checks:
- **Sun:** perihelion ~ 1.471e8 km, aphelion ~ 1.521e8 km (range 5e6 km). The 018 sample values are all within this band. Full-year min/max from the manifest: 147,099,933 / 152,087,750 km, matching the perihelion/aphelion band within 0.001%.
- **Moon:** perigee ~ 356,500 km, apogee ~ 406,700 km. Manifest min/max: 356,779 / 406,389 km, matching the perigee/apogee band within 0.1%.

Both snapshots are within the physical ephemeris bands. **Verdict:** Snapshot data is physically valid.

### 3.5 Cadence verification

Step sizes are exactly 86400.0 s throughout (first 5, last 5, min, max all equal 86400.0). The acquisition pattern (declared in the Sun and Moon fetch scripts and verified in the manifest) uses `STEP_SIZE: "1d"`, and the lab's `_validate_vector_response` (in `eclipseTiming/fetch_horizons_sun_snapshot.py:106-115`) checks that all consecutive row deltas are within 2e-4 s of 86400.0.

**Verdict:** 1-day uniform cadence. Confirmed.

---

## 4. Linear interpolation error analysis

### 4.1 Position error from linear interpolation between daily snapshots

The 017 and 018 implementations both use **linear interpolation** between adjacent snapshot rows (017 at `lunisolarVerification/experiment.py:390-401`, 018 at `lunisolarReconciliation/experiment.py:243-258`):

```python
idx = int(np.searchsorted(t_s, t_query_s))
t_lo = t_s[idx - 1]
t_hi = t_s[idx]
frac = (t_query_s - t_lo) / (t_hi - t_lo)
return r[idx - 1] + frac * (r[idx] - r[idx - 1])
```

The midpoint of an interval is where linear interpolation has the largest error. For a body moving with velocity `v` and acceleration `a` over an interval of length `dt`, the midpoint error is:

    |Δr| ~ (1/8) * |a| * dt^2

(perpendicular to the line of motion; tangent error is ~0 by the linear interpolation's exact-match condition).

**Sun** at Earth:
- Geocentric speed: ~30 km/s (Earth's heliocentric orbital speed, modulated by ±3.4% for orbital eccentricity)
- Over 86400 s: 30 * 86400 = 2.6e6 km
- Angular rate from Earth: 2.6e6 / 1.5e8 = 0.0174 rad/day = ~1.0 deg/day
- Centripetal acceleration (Earth's orbit around Sun barycenter): v^2 / r ~ 30^2 / 1.5e8 = 6e-6 km/s^2 = 6 mm/s^2
- Position error at midpoint: 0.125 * 86400^2 * 6e-9 km/s^2 = 0.125 * 7.46e9 * 6e-9 = **5.6 km**

**Moon** at Earth:
- Geocentric speed: ~1.022 km/s (mean)
- Over 86400 s: 1.022 * 86400 = 88,300 km
- Angular rate: 88,300 / 384,400 = 0.230 rad/day = 13.2 deg/day
- Centripetal acceleration: v^2 / r ~ 1.022^2 / 384,400 = 2.72e-6 km/s^2 = 2.7 mm/s^2
- Position error at midpoint: 0.125 * 86400^2 * 2.72e-9 km/s^2 = 0.125 * 7.46e9 * 2.72e-9 = **2.5 km**

So linear interpolation of Sun/Moon at daily cadence gives ~5.6 km (Sun) and ~2.5 km (Moon) RMS position error at the midpoint of each 86400 s interval. The Moon error is smaller because the centripetal acceleration is smaller (lower speed).

### 4.2 Impact on third-body acceleration and on RAAN rate

The third-body acceleration scales as `mu_3 / |r_3|^3` (direct term) and `mu_3 / |r_3|^3 * |r_sat| / |r_3|` (indirect term relative to direct). A 5 km error in `r_3` for the Sun (|r_3| ~ 1.5e8 km) produces a fractional error in `1/|r_3|^3` of order `3 * 5 / 1.5e8 = 1e-7`, completely negligible for RAAN purposes.

For the Moon, a 2.5 km error in `r_3` (|r_3| ~ 3.8e5 km) gives fractional error of order `3 * 2.5 / 3.8e5 = 2e-5`, still negligible.

The dominant effect of the position error is on the **direction** of `r_3`, not its magnitude. A 5 km transverse error in the Sun's geocentric direction at 1.5e8 km corresponds to `5 / 1.5e8 = 3.3e-8 rad = 6.9e-6 arcsec`. For comparison, the precession-induced frame error is 0.66 deg = 2.4e-3 rad — **8 orders of magnitude larger** than the linear interpolation direction error.

The corresponding error in `sin(2(i - i_3))` is ~3.3e-8 * 2 = 6.7e-8 fractional. Multiplied by the secular RAAN rate of ~1.5e-4 deg/day gives ~1e-11 deg/day, or ~4e-9 deg/year. **Completely negligible** compared to the 0.5 deg/year J2 closure residual at h=600 km.

**Verdict:** Linear interpolation between daily snapshots is **adequate** for the experiment's purposes. The ~5.6 km Sun and ~2.5 km Moon position errors produce RAAN rate biases ~7 orders of magnitude smaller than the headline findings.

### 4.3 Clamping at the end of the snapshot (final 12 hours)

The last 12 hours of the 1-year propagation (`t_s ∈ [852033600, 852076800]`) clamp to the snapshot's final row at `t_s = 852033600` (2026-12-31 00:00 TDB). At this clamped time, the Sun and Moon geocentric positions are FIXED at their end-of-year values.

For the Sun, the 1-day motion is ~1 deg/day, so clamping for 12 hours means the Sun "stops" at its year-end position for half a day. This is equivalent to an additional ~0.5 deg direction error accumulating linearly from 0 to 0.5 deg over those 12 hours, with RMS ~0.25 deg.

For the Moon, the 1-day motion is ~13 deg/day, so the Moon's direction freezes for 12 hours — but more importantly, the Moon's mean anomaly advances by ~6.5 deg in 12 hours, which is a SIGNIFICANT angular shift in the lunar direction.

**Impact on the 1-year RAAN rate:** the clamped 12-hour window contributes negligibly to the secular RAAN rate (the secular drift is the integral of the torque, and a 12-hour clamp at the end is one part in 730 of the integration window). The instantaneous torque error at the end is at most ~7e-6 deg/day from the linear-interp direction error in §4.2, but only for 12 hours; the time-integrated impact is ~7e-6 * 12/24 / 365 ≈ 1e-8 deg/day. **Completely negligible.**

**Verdict:** Clamping at the end is benign for the 1-year rate measurement.

---

## 5. Epoch and time scale

### 5.1 The propagator's t0

Both 017 (line 555) and 018 (line 449) use:

```python
t0 = 820540800.0  # 2026-01-01 12:00 TT (lab convention)
```

820540800 / 86400 = 9500 days. 9500 / 365.25 = 26.007 yr. J2000 = 2000-01-01 12:00 TT (the standard J2000 epoch, defined as JD 2451545.0 = noon TT). So `t0 = 820540800 s` = 2000-01-01 12:00 + 9500 days = **2026-01-01 12:00 TT**. This is the lab convention; the 018 comment at line 449 explicitly states "2026-01-01 12:00 TT (lab convention)".

### 5.2 Time-scale consistency between snapshot and propagator

**Snapshot time scale:** TDB (per Horizons echoes in both snapshot headers: "A.D. 2026-Jan-01 00:00:00.0000 TDB", "A.D. 2027-Jan-01 00:00:00.0000 TDB", `TIME_TYPE: 'TDB'` in the fetch parameters, and `validation.header.start_time_echo` in both manifests).

**Propagator time scale:** TT-like (per the lab convention). The 018 docstring at line 685 says "s since J2000 (TT-like)". The 017 docstring at line 116 says "ECI ICRF/TDB; ...Sun and Moon from byte-pinned DE441 geocentric vectors; lab's ECI is pseudo-inertial at LEO precision".

**TDB vs TT difference:**
- The IAU-1991 definition (Fairhead & Bretagnon 1990): TDB - TT has periodic terms with peak ~1.7 ms (annual + shorter).
- The IAU-2006 definition (capitaine et al.): TCB - TDB = (L_B / (1 - L_B)) (JD - T_0) * 86400 s + periodic terms (L_B ~ 1.55e-8).
- For 1-year LEO propagation at LEO precision, the TDB-TT difference is at the sub-second level — far below the 1 s/yr event-rate floor and well below the 1e-3 deg/day RAAN rate noise.

**Verdict:** The snapshot is in TDB; the propagator is in TT-like. The TDB-TT difference (~1.7 ms peak) is **sub-second over 1 year**, far below the experiment's precision. The time-scale mismatch is harmless at LEO precision for secular rate measurements.

If the experiment ever required sub-meter position precision over multi-day arcs, the TDB-TT difference would matter. For the 017/018 1-year RAAN rate measurements, it does not.

### 5.3 The 12-hour offset between snapshot start and lab_t0

As noted in §3.3:
- Snapshot first row: 2026-01-01 00:00 TDB
- Lab t0: 2026-01-01 12:00 TT (lab convention)
- Offset: 12 hours

This offset is **consistent** with the lab's J2000 noon convention (JD_J2000 = 2451545.0 = noon). The propagator's `t0` is 12 hours past midnight UTC = noon UTC of 2026-01-01, which is the natural extension of the J2000 noon anchor. The snapshot starts at midnight (JD .5). The 12-hour offset is **not a bug** — it's a deliberate alignment of the lab's J2000-noon-anchored time grid with the snapshot's midnight-anchored grid.

The snapshot is in TDB; if the lab's t0 is treated as TT, then 2026-01-01 12:00 TT corresponds to 2026-01-01 12:00:01.7 TDB (peak annual TDB-TT offset). The 1.7-second difference is well below the 12-hour offset and below the experiment's noise floor.

**Verdict:** Epoch and time scale are internally consistent for the experiment's precision.

---

## 6. Linear interpolation at high frequencies (Sun and Moon angular rates)

The Sun's apparent angular rate from Earth is ~1.0 deg/day; the Moon's is ~13.2 deg/day. Over 86400 s, the linear interpolation has a midpoint error perpendicular to the line of motion of ~0.5 deg * (dt / period) (this is the trapezoidal-rule error). With `dt = 86400 s` and `period = 86400 s` (Sun), the Sun's interpolation error is ~0.5 deg at midpoint (but reduced to the position-error formula in §4.1, which is ~6 km = 0.0023 deg direction, much smaller).

Wait — the 3 arcsecond figure quoted in the audit prompt is wrong. Let me compute: linear interpolation of a body moving at 1 deg/day over 86400 s has a midpoint angular error of ~0.5 * dt^2 * (angular acceleration / 2) for a body on a circular orbit at constant angular rate. But for a body at constant ANGULAR RATE, linear interpolation of the POSITION is exact (no error). The Sun and Moon have angular rates that vary slowly with orbital eccentricity; the dominant position error is from SECOND-derivative terms in the orbital motion, not from "lagging behind the actual position by 3 arcseconds at the midpoint".

The 3 arcsecond figure arises only if you assume the linear interpolation "lags" the true position by half the daily motion, which is a misconception. Linear interpolation does not lag a uniformly moving target — it exactly matches the start and end of each interval. The error is in the higher-order terms (curvature).

The correct first-order interpolation error for a target at constant angular rate is **zero** (linear interpolation is exact for any linear function of time, and a uniformly-moving body's position IS linear in time). The error scales with the second derivative (centripetal acceleration for circular orbits). For the Sun at 6 mm/s^2 centripetal: position error ~5.6 km at midpoint (§4.1) = 5.6 / 1.5e8 rad = 3.7e-8 rad = **7.7e-3 arcsec** (NOT 3 arcsec).

The "3 arcsecond" figure is off by **a factor of ~400** — the actual error is ~0.008 arcsec. Either way, both are negligible for the RAAN rate calculation.

**Verdict:** The linear interpolation error in the Sun's geocentric direction is ~8 milliarcsec at the midpoint, NOT 3 arcsec. Both are negligible for the secular RAAN rate measurement. The "3 arcsec" figure quoted in the audit brief appears to be an order-of-magnitude error in the prompt.

---

## 7. Code-hash binding and reproducibility

The 018 results.json records SHA-256 hashes of all source files in the `code_sha256` field, computed by `code_hashes()` at lines 599-616 of `lunisolarReconciliation/experiment.py`:

```python
def code_hashes() -> dict:
    here = Path(__file__).resolve().parent
    lab_root = here.parents[3]
    files = {
        "experiment.py": here / "experiment.py",
        "lab_utils/orbits.py": lab_root / "src" / "lab_utils" / "orbits.py",
        "lab_utils/earth_frames.py": lab_root / "src" / "lab_utils" / "earth_frames.py",
        "lab_utils/integrators.py": lab_root / "src" / "lab_utils" / "integrators.py",
        "lab_utils/results.py": lab_root / "src" / "lab_utils" / "results.py",
        "lab_utils/__init__.py": lab_root / "src" / "lab_utils" / "__init__.py",
        "moon_reference_snapshot.txt": MOON_SNAPSHOT_PATH,
        "sun_reference_snapshot.txt": SUN_SNAPSHOT_PATH,
    }
    return {name: _file_sha256(p) for name, p in files.items()}
```

**Verified:** The hashes are computed and recorded. This satisfies the lab's reproducibility doctrine (byte-pinning the analysis inputs). The 017 experiment.py has the equivalent function at lines 815-836 with the same fields (minus `lab_utils/results.py` is in the 018 version because 018 imports it; 017 also imports `lab_utils.results` per line 119 of 017).

**Verdict:** Code-hash binding is in place. Reproducibility is enforced.

---

## 8. Verdict

### 8.1 Bugs identified

| # | Description | File:line | Severity | Impact on secular RAAN rate |
|---|-------------|-----------|----------|----------------------------|
| 1 | **Sign bug in precession `_rot3`**: 018's `_rot3` uses the transpose of the standard form, causing the rotation about Z to be in the wrong direction. The result is a ~0.66 deg frame error at 2026 (vs ~0.5 deg in 017 with no precession, and 0 deg with correct-sign precession). | `lunisolarReconciliation/experiment.py:132-135` and `153-156` | **Significant** (the 018 contract claim that precession fixes the Track D frame mismatch is incorrect) | ~2.5e-3 deg/year prograde; ~3% of the corrected formula's magnitude, well below the 9.78x short-period residual at i_sso |
| 2 | **Comment understatement**: 017 docstring says indirect Moon term is "order 1e-5 of the direct term" (line 421). For the Moon, the indirect/direct ratio is ~2e-2, ~350x larger than the comment. The Sun indirect/direct ratio is ~5e-5, matching the comment. | `lunisolarVerification/experiment.py:421` | **Cosmetic** | None — code is correct |
| 3 | **Final 12 hours of propagation clamp**: The lab propagator runs 12 hours past the snapshot's last row, and `_interp_snapshot_precessed` clamps to the final row. The first 12 hours are inside the snapshot. | `lunisolarReconciliation/experiment.py:247-249` and `lunisolarVerification/experiment.py:396-399` | **Negligible** | < 1e-8 deg/day on RAAN rate |
| 4 | **018 docstring claim**: "The IAU-1976 precession rotation has been applied ... This fixes the Track D frame-mismatch finding (0.4 deg offset at 2026 between ICRF and mean-of-date)." (lines 39-42, 740-745). With the sign bug in `_rot3`, the rotation does NOT fix the frame mismatch — it reverses the sign of the residual frame error. The 0.4 deg offset becomes 0.66 deg (the sum, not the difference). | `lunisolarReconciliation/experiment.py:39-42` | **Significant** (misleading docstring; the documented "fix" does not actually fix the frame mismatch) | Same as bug #1 |
| 5 | **TT vs TDB convention mix**: The snapshot is in TDB; the propagator is in TT-like. Difference is sub-second over 1 year; below noise floor. | `lunisolarVerification/experiment.py:555`, `lunisolarReconciliation/experiment.py:449` | **Cosmetic** | < 1e-3 deg/day on RAAN rate (well below noise) |

### 8.2 Strengths (no bugs found)

- **Third-body acceleration formula**: Correct in both 017 (`lunisolarVerification/experiment.py:418-441`) and 018 (`lunisolarReconciliation/experiment.py:269-289`). Direct + indirect terms with correct signs.
- **Lieske polynomial coefficients**: Match the IAU-1976 standard (Lieske et al. 1977).
- **Snapshot SHA-256**: Both Sun (`06d54fb...`) and Moon (`65f1d67f...`) match their manifests byte-for-byte.
- **SOE/EOE parsing**: Extracts all 366 daily rows; 018's parser verified against the manifest-declared `validation.rows: 366`.
- **Time conversion**: `(JD - JD_J2000) * 86400` correctly implements the J2000-noon-anchored seconds-since-J2000 convention.
- **Distance band**: Both snapshots are within physical perihelion/aphelion (Sun) and perigee/apogee (Moon) bands.
- **Cadence**: Exactly 86400 s uniform throughout.
- **Linear interpolation**: Error is ~5.6 km (Sun) and ~2.5 km (Moon) at midpoint, producing ~7e-9 deg/year RAAN bias. **Negligible**.
- **Code-hash binding**: Implemented in both 017 (line 815-836) and 018 (line 599-616).
- **Force-level identity check** (018 line 297-348): Independent algebraic verification of the acceleration formula passes at machine precision.
- **RK4 propagator** (`lab_utils/integrators.py:21-29`): Standard textbook RK4 with proper stage-time evaluation (non-autonomous-safe signature `f(t, x)`); verified by 017 and 018 convergence ladders.
- **J2 RHS** (`lab_utils/orbits.py:259-281`): Standard J2 acceleration formula with short-circuit for `j2=0`; pinned by the lab_utils tests and used by both 017 and 018.

### 8.3 Quantified impact of the precession sign bug on the 018 headline

The 018 headline at `lunisolarReconciliation/experiment.py:20-29` states:

> "The CORRECT secular quadrupole formula (Track B independent derivation):
> dO/dt = (3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i - i_3) / sin(i)
> At h=600 km i_sso=97.79 deg the corrected formula gives +1.35e-4 deg/day (prograde), matching the numerical SIGN and ~10x smaller in magnitude."

With the sign bug in `_rot3`, the 1-year numerical fit still gives the same +1.28e-3 deg/day (prograde) because the bug only changes the direction of the frame rotation; the dominant secular signal comes from the Sun/Moon third-body physics, not the frame alignment. The corrected formula's value +1.35e-4 deg/day is unaffected (no frame alignment is involved in the analytic formula; the formula assumes the Moon's `i_3` is its mean inclination to the equator).

The 9.78x ratio between numerical and corrected formula is therefore dominated by the unmodelled short-period terms (evection + variation + lunar nodal), as 018 claims. The precession sign bug contributes an additional ~3% bias on top of this (within the noise of the 9.78x residual).

**Impact on Exp 019 candidate directions:** The precession sign bug does NOT invalidate the 018 corrected formula vs numerical agreement in SIGN. It DOES mean that the "precession_comparison" experiment in 018 will measure a smaller bias than the 018 docstring claims (the wrong-sign precession is similar in magnitude to no precession, so the comparison's `with_precession` minus `without_precession` difference is ~2e-3 deg/year, not 0.5 deg/year as the docstring suggests). The 016 LST-drift budget's "lunisolar correction" was already retracted in audit-018; the 018 precession sign bug does not affect the 016 remediation status.

### 8.4 Final verdict

The numerical implementation (third-body acceleration, snapshot loading, parsing, interpolation, time conversion, RK4 propagation, J2 RHS, code-hash binding) is **correct**. The single significant bug identified is the **sign error in the 018 precession `_rot3` matrix**, which causes the precession rotation to be applied in the wrong direction about the Z-axis, leaving a ~0.66 deg frame error instead of fixing the ~0.5 deg Track D frame mismatch.

The impact of this bug on the 018 headline finding (corrected formula vs numerical, 9.78x ratio dominated by short-period physics) is **small** (~3% bias on the corrected formula's magnitude), and the **corrected secular formula's agreement in SIGN with the numerical 1-year measurement is preserved** (both are prograde). The bug does not invalidate the 018 finding.

**For Exp 019 (Lunisolar Long-Period Terms):** Before incorporating the 018 precession rotation into a 019 multi-year extension, the sign bug should be fixed (replace `_rot3` with the eclipseTiming convention: `[[c, -s, 0], [s, c, 0], [0, 0, 1]]`, or equivalently replace `_rot3(-z) @ _rot2(theta) @ _rot3(-zeta)` with `_rot3(z) @ _rot2(-theta) @ _rot3(zeta)` — but using the same matrix convention). Alternatively, the eclipseTiming `precession_matrix_mod_from_j2000` function can be imported directly (it's not currently graduated into `lab_utils`, but the implementation is verified at `eclipseTiming/experiment.py:266-277`).

**Overall:** The implementation is **sound for the experiment's purposes** (1-year RAAN rate comparison of corrected formula vs numerical, 9.78x residual dominated by short-period physics). The precession sign bug is a real defect that should be corrected before the result is used for multi-year or sub-milliarcsec-precision work, but it does not change the qualitative conclusion of 018.

---

## Audit-trail references

- 017 third-body acceleration: `lunisolarVerification/experiment.py:418-441`
- 018 third-body acceleration: `lunisolarReconciliation/experiment.py:269-289`
- 017 SOE/EOE parser: `lunisolarVerification/experiment.py:362-381`
- 018 SOE/EOE parser: `lunisolarReconciliation/experiment.py:209-228`
- 017 linear interpolation: `lunisolarVerification/experiment.py:390-401`
- 018 linear interpolation: `lunisolarReconciliation/experiment.py:243-258` (with precession option)
- 018 precession rotation (BUG): `lunisolarReconciliation/experiment.py:147-156`
- 018 _rot3 (BUG): `lunisolarReconciliation/experiment.py:132-135`
- 018 _rot2: `lunisolarReconciliation/experiment.py:138-141`
- eclipseTiming reference precession (CORRECT): `eclipseTiming/experiment.py:266-277`
- eclipseTiming reference _rot3 (CORRECT): `eclipseTiming/experiment.py:255-258`
- 018 force-level identity check: `lunisolarReconciliation/experiment.py:297-348`
- 017 propagate_one_altitude: `lunisolarVerification/experiment.py:514-639`
- 018 propagate_one: `lunisolarReconciliation/experiment.py:421-475`
- 018 code_hashes: `lunisolarReconciliation/experiment.py:599-616`
- 017 code_hashes: `lunisolarVerification/experiment.py:815-836`
- lab_utils RK4 propagator: `src/lab_utils/integrators.py:21-58`
- lab_utils J2 RHS: `src/lab_utils/orbits.py:259-281`
- lab_utils JD_J2000: `src/lab_utils/earth_frames.py:67`
- Sun MANIFEST: `eclipseTiming/reference/MANIFEST.json`
- Moon MANIFEST: `lunisolarVerification/reference/MANIFEST.json`
- Sun fetch acquisition pattern: `eclipseTiming/fetch_horizons_sun_snapshot.py`

---

**Track D conclusion:** Implementation is sound for the experiment's headline finding; the precession rotation has a sign bug that should be fixed before multi-year or high-precision reuse.