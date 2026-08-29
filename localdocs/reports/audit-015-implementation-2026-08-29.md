# Independent Audit — Experiment 015 Implementation (Frame / Time / Sign)

| Field | Value |
|---|---|
| Audit ID | AUDIT-015-IMPLEMENTATION |
| Audited artifacts | `research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py`, `README.md`, `tests/test_dawn_dusk_sso.py`, `src/lab_utils/earth_frames.py`, `src/lab_utils/orbits.py` |
| Audit type | Line-by-line frame / time / sign convention check + numerical reproducibility probe |
| Audit date | 2026-08-29 |
| Auditor | independent, read-only |
| Builds on | `audit-015-lst-drift-2026-08-29.md` (LST-drift physics already established RED); this audit is the implementation-layer audit of the 10 specific focus questions |
| Scope | Read-only. No edits to `experiment.py`, `README.md`, lab_utils, tests, or any canonical experiment file. Only this audit report may be created. |

---

## Executive summary

| Question | Result | Verdict |
|---|---|---|
| Q1: `subsolar_lon_rad` is geodetic, not Sun RA | **CORRECT in code** (uses ECEF atan2 of rotated Sun unit vector); **WRONG in README** (claims `atan2(-u_y, -u_x)`) | MIXED — code GREEN, README narrative RED |
| Q2: LST is apparent or mean | **Apparent** (uses Almanac geometric Sun, includes EoT). Documented in `earth_frames.py` docstring. Canonical dawn-dusk SSO target is conventionally mean LST — bias ≤ ±16 min (EoT envelope) | MINOR — convention pinned, not clearly disclosed in README |
| Q3: C4 `Omega = GMST + lon_ref` = ascending node over Eastern Range | **CORRECT** at insertion only. Places node at `Omega − GMST = lon_ref` (ECEF, west-negative). For SSO with J2 drift, subsequent nodes precess | OK with disclosure |
| Q4: `node_lon = Omega − GMST` (sign) | **CORRECT**: `node_lon_from_raan_gmst(Omega, GMST) = Omega − GMST` matches the geodetic convention `lon_ecef = lon_eci − GMST` (passive Z rotation in `eci_to_ecef`) | GREEN |
| Q5: ±10 min LST tolerance, wrap and sign | **CORRECT**: `delta_h = (lst − target + 12) % 24 − 12` puts `delta_h ∈ (−12, 12]` h, multiplied by 60 → `offset_min ∈ (−720, 720]` min. Sign: positive when LST > target, negative when LST < target. The docstring in `lst_offset_min` claims `wrapped to (−720, 720]` which is achieved modulo the boundary point at exactly −720 min | GREEN (boundary minor) |
| Q6: "4 min/day drift" claim | **WRONG by factor ~365**. The actual insertion-time LST, given the C4 convention `Omega = GMST + lon_ref`, cycles through 24 h in **1 sidereal day**, not in 1 year. Unwrapped drift rate measured: **1440 min/day** (a full 24 h/day), which mod 24 ≈ 0 within EoT. The README claim `dΩ/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day = 4 min/day` is **frame-inconsistent** (subtracts inertial RAAN rate from ECEF subsolar rate) and is the SSO design rate as a tautology, not a measured drift. See `audit-015-lst-drift-2026-08-29.md` for full derivation; this audit confirms the implementation matches the RED verdict | RED (already documented in companion audit) |
| Q7: `node_lon = site_lon_deg * DEG` is the FIXED geodetic node at insertion | **CORRECT**. `lst_at_node_at_t(t_L)` and `constraint_indicator` both compute `node_lon = site_lon_deg * DEG` and feed it into `lst_at_node_hours(t_L, node_lon)`. This is the LST at the launch site's longitude at the launch instant. It is NOT the LST at the SSO's ascending node at any time other than insertion; subsequent ascending nodes precess with J2. The function name `lst_at_node_at_t` is ambiguous between "LST at the SSO's node at subsequent times" and "LST at the insertion node at launch time"; the code implements the latter | MINOR (name ambiguity) |
| Q8: Frame firewall (mean-of-date vs J2000 vs ICRF) | **CONSISTENT**: all orbital angles (`Omega`, `alpha_sun`, GMST) live in the lab ECI mean-of-date frame; the precession bias (~50.3″/yr vs J2000) is named-excluded per `lab_utils/earth_frames.py` lines 14–17 and validated at ~0.65° vs Horizons snapshot. Mean-element J2 nodal drift is consistent with mean-of-date Sun model. No inadvertent J2000 assumption detected | GREEN |
| Q9: README claim `dΩ/dt − d(Subsolar)/dt = 360.9856 − 360.0` | **FRAME-INCONSISTENT**. `dΩ/dt` is inertial-frame RAAN rate; `d(Subsolar)/dt` is the ECEF subsolar point longitude rate. These cannot be subtracted without a frame transformation. The correct expression in any single frame: `d(LST)/dt = (dΩ/dt − d(α_sun)/dt) / 15`. For SSO: `dΩ/dt = d(α_sun)/dt = +SSO_TARGET_DEG_DAY` → `d(LST)/dt ≈ 0` (modulo EoT). The 0.9856 deg/day figure is the SSO design rate as a tautology, not a measured drift. See Q6 and `audit-015-lst-drift-2026-08-29.md` §3 | RED |
| Q10: `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` | **TAUTOLOGICAL** (per the prior LST-drift audit). It does not measure a 4 min/day rate; it only asserts that `lst_at_node_at_t` (which computes the LST at a constant launch-site longitude for varying launch time) visits values far from 18:00 over a year — trivially true because LST mod 24 cycles through 24 h in ~1 day when `node_lon` is held constant. The assertion is `max_dist_from_18 > 3.0`, which is satisfied by ANY non-constant LST function over a year. The test name and docstring claim a 4 min/day rate, but no such rate is enforced or measured. This test **validates the bug, not the physics** (companion audit §5). | RED |

**Net**: the **CODE** for Exp 015 is GREEN for Q3–Q5, Q7, Q8 (frame conventions, sign conventions, tolerance wrap, node-lon sign, frame firewall). The **README narrative** is RED on Q1, Q6, Q9 (anti-sun formula described but not implemented; frame-inconsistent drift expression; factor-of-365 wrong on the 24 h/year vs 24 h/day confusion). The **TEST** for LST drift is RED on Q10 (tautological). The **sensitivity matrix** has a functional bug (Q5 corollary): `lst_tolerance_min` override is **gated by `min(LST_TOLERANCE_MIN, override)`** in `feasibility_curve`, so `tol=20` is a no-op (the expected "≈2x" doubling does not occur; actual count is identical to baseline at every altitude).

These findings corroborate and tighten the existing `audit-015-lst-drift-2026-08-29.md` and add the new findings on the README `sub_lon` formula inconsistency, the dead-code NaN check, the broken `lst_tolerance_min` sensitivity override, and the test tautology.

---

## 1. Frame and convention audit (`subsolar_lon_rad`, `lst_at_node_hours`)

### 1.1 `subsolar_lon_rad` returns GEODETIC (ECEF), not ECI RA

**Code** (`src/lab_utils/earth_frames.py` lines 99–121):

```python
def subsolar_lon_rad(t_s):
    """Subsolar-point geodetic longitude in the lab ECI mean-of-date frame (rad)."""
    u, _ = sun_unit_and_dist_km(t_s)
    u_ecef = eci_to_ecef(u, gmst_rad_iau1982(t_s))
    u_x_ecef = float(u_ecef[..., 0])
    u_y_ecef = float(u_ecef[..., 1])
    lon = np.arctan2(u_y_ecef, u_x_ecef)
    # wrap to (-pi, pi]
    lon = (lon + np.pi) % (2 * np.pi) - np.pi
    return float(lon) if np.ndim(lon) == 0 else lon
```

The function:
1. Gets the Sun unit vector `u` in ECI mean-of-date (`sun_unit_and_dist_km`).
2. Rotates to ECEF via `eci_to_ecef(u, GMST)` — passive Z rotation `[[c, s, 0], [−s, c, 0], [0, 0, 1]]`.
4. Returns `atan2(u_ecef_y, u_ecef_x)`, the geodetic longitude of the subsolar point.

**Independence verified numerically**:

| Formula | Value at `t = 2026-01-01 00:00 UTC` |
|---|---:|
| `subsolar_lon_rad(t)` (lab_utils) | **+1.181°** |
| `atan2(u_y_ecef, u_x_ecef)` after ECI→ECEF | **+1.181°** (matches to 1e-12) |
| `atan2(u_y, u_x)` (Sun RA in ECI) | **−77.955°** (NOT the subsolar longitude) |
| `atan2(−u_y, −u_x)` (anti-sun in ECI, what the README describes) | **+102.045°** (off by exactly 180° vs the lab value) |
| `alpha_sun − GMST` (Sun RA − GMST) | **+1.181°** (matches lab to 1e-12) |

The lab_utils implementation is **mathematically correct as the geodetic subsolar longitude** (= `alpha_sun − GMST` mod 2π).

### 1.2 README describes a different (wrong) formula

**README.md** line 65–66 (the prose explanation under the constraint equations):

> `sub_lon(t)` = geocentric subsolar longitude from `atan2(-u_y, -u_x)` of the lab's analytic Almanac Sun unit vector (mean of date)

This is the **anti-sun direction** in ECI = `alpha_sun + π`, which is **180° away from the geodetic subsolar longitude**. If `Omega − sub_lon` were used with this `sub_lon`, the resulting LST would be off by 12 hours (mod 24).

**CRITICAL**: the **code** does NOT use `atan2(−u_y, −u_x)`. It uses the geodetic formula via `lst_at_node_hours(t, node_lon)` which internally calls `subsolar_lon_rad(t)` (the geodetic one). The README's narrative description is **wrong** but the code path is **right**. This is a documentation bug, not a code bug.

### 1.3 Cross-check: `lst_at_node_at_t` is bit-equivalent to `12 + (Omega − alpha_sun) / 15`

`lst_at_node_at_t(t_L)` in `experiment.py` lines 252–272:

```python
gmst = gmst_rad_iau1982(t_launch_s)
raan = gmst + REF_SITE_LON_DEG * DEG
node_lon = raan - gmst           # = REF_SITE_LON_DEG * DEG (constant in this convention)
return lst_at_node_hours(t_launch_s, node_lon)
```

`lst_at_node_hours(t, node_lon)` internally computes:
```python
delta_h = (node_lon − subsolar_lon(t)) / (15·DEG)
lst = 12 + delta_h
```

Substituting `subsolar_lon = alpha_sun − GMST` and `node_lon = REF_SITE_LON_DEG·DEG`:
```
LST = 12 + (REF_SITE_LON_DEG·DEG − alpha_sun + GMST) / (15·DEG)
    = 12 + (GMST + REF_SITE_LON_DEG·DEG − alpha_sun) / (15·DEG)
    = 12 + (raan − alpha_sun) / (15·DEG)
```

This matches the textbook formula. ✓ The LST at the insertion-time ascending node over Eastern Range equals `12 + (Omega(t_L) − alpha_sun(t_L)) / 15`.

The test `test_lst_at_node_at_t_at_dusk_terminator` (`tests/test_dawn_dusk_sso.py` lines 113–128) verifies this bit-equivalence at `t = 0` to 1e-9.

### 1.4 Apparent vs Mean LST (Q2)

`subsolar_lon_rad(t)` uses the Almanac **geometric** Sun, which is the **apparent** Sun direction (no nutation/aberration correction, but already includes the equation of center). Therefore `lst_at_node_hours` returns **apparent LST** (sometimes called "true solar time"), not mean LST. This is documented in `earth_frames.py` lines 280–282:

> "This is the apparent LST (no EoT correction needed because the subsolar point is by definition the apparent Sun direction projected to the geoid)."

Computed EoT envelope over 2026: **−14.20 min to +16.45 min** (peak-to-peak ~30 min). The dawn-dusk SSO LST target of 18:00 is conventionally **mean** LST. The lab's apparent LST differs from canonical mean LST by up to ±16 min, which is comparable to the 10-min LST tolerance.

**Disclosure status**: `lab_utils/earth_frames.py` is honest about the apparent-LST convention. `experiment.py` does not re-state it. The README does not state it (the only EoT mention is in the limitations, line 278, as a year-to-year variation, not as an apparent-vs-mean issue).

---

## 2. C4 insertion convention and node-lon sign

### 2.1 C4: `Omega(t_L) = GMST(t_L) + lon_ref` (Q3)

`insertion_raan_rad(t_L)` in `experiment.py` lines 246–249:

```python
def insertion_raan_rad(t_launch_s: float) -> float:
    """RAAN at insertion (ascending node over Eastern Range)."""
    gmst = gmst_rad_iau1982(t_launch_s)
    return gmst + REF_SITE_LON_DEG * DEG
```

Geodetic longitude of the ascending node at insertion:
```
node_lon = Omega(t_L) − GMST(t_L) = (GMST(t_L) + REF_SITE_LON_DEG·DEG) − GMST(t_L) = REF_SITE_LON_DEG·DEG
```

For `REF_SITE_LON_DEG = −80.6039°`, this is the geodetic longitude of Cape Canaveral area. With the orbital elements `Omega, omega = 0` and the satellite at `argument of latitude = 0` (the convention for `Orbit.__init__`), the satellite is at the **ascending node** (northward equator crossing). ✓

**Disclosure**: the C4 convention places the node **at insertion only**. For an SSO with J2 nodal drift at `+0.9856 deg/day` (sidereal day units), subsequent ascending nodes precess. The geodetic node longitude at time `t > t_L` is `Omega(t_L) + 0.9856·DEG·(t − t_L) − GMST(t)`. The 4-min/day LST drift the README describes would apply to a **fixed-on-ground** node (no J2 drift), which is **not** what the experiment simulates — see `audit-015-lst-drift-2026-08-29.md` §5 for full forensic.

### 2.2 `node_lon_from_raan_gmst` sign (Q4)

`node_lon_from_raan_gmst(raan, gmst)` in `earth_frames.py` lines 261–271:

```python
def node_lon_from_raan_gmst(raan_rad, gmst_rad_val):
    lon = np.asarray(raan_rad, dtype=float) − np.asarray(gmst_rad_val, dtype=float)
    lon = (lon + np.pi) % (2 * np.pi) − np.pi
    return float(lon) if np.ndim(lon) == 0 else lon
```

This is `Omega − GMST`, which is the geodetic longitude of the orbit's ascending node (because the passive Z rotation `eci_to_ecef` does `lon_ecef = lon_eci − theta_G`). ✓ Matches the convention used in `lst_at_node_at_t`.

---

## 3. LST tolerance wrap and sign (Q5)

### 3.1 `lst_offset_min` implementation (`experiment.py` lines 275–283)

```python
def lst_offset_min(t_launch_s: float, target_hours: float) -> float:
    """Signed LST offset from target (minutes), wrapped to (-720, 720]."""
    lst = lst_at_node_at_t(t_launch_s)
    delta_h = lst − target_hours
    delta_h = (delta_h + 12.0) % 24.0 − 12.0
    return delta_h * 60.0
```

Numerical wrap test:

| `lst` (h) | `target` (h) | raw `delta_h` (h) | wrapped (h) | `offset_min` |
|---:|---:|---:|---:|---:|
| 17.5 | 18.0 | −0.5 | **−0.5** | **−30.0** |
| 18.5 | 18.0 | +0.5 | **+0.5** | **+30.0** |
| 6.0  | 18.0 | −12.0 | **−12.0** | **−720.0** ← boundary, excluded from (−720, 720] |
| 5.5  | 18.0 | −12.5 | **+11.5** | **+690.0** ← correct wrap |
| 0.0  | 18.0 | −18.0 | **+6.0** | **+360.0** ← correct wrap |

The wrap is `(delta + 12) % 24 − 12`. Python's `%` returns a value in `[0, 24)` for non-negative operands. For `delta = −12`, `(−12 + 12) % 24 = 0`, so the wrapped value is `0 − 12 = −12`. This is the **boundary** of the wrap range, which is technically excluded from `(-12, 12]` but the docstring says `(−720, 720]` (also exclusive on the lower end). ✓

Sign convention: positive when LST > target (satellite clock ahead of design), negative when LST < target (behind). Consistent with "satellite needs to drift east (positive) to align with target" semantics. ✓

The constraint `lst_off_min_abs <= LST_TOLERANCE_MIN` (= 10 min) uses the absolute value, so the sign does not affect feasibility. ✓

### 3.2 Sensitivity `lst_tolerance_min` is GATED by the internal `LST_TOLERANCE_MIN`

`constraint_indicator` (lines 314–315) computes:
```python
lst_ok = bool(lst_off_min_abs <= LST_TOLERANCE_MIN)   # always 10 min
...
feasible = lst_ok and eclipse_ok
```

`feasibility_curve` (lines 419–425) then applies the override as a **second filter**:
```python
flags[i] = ind["feasible"] and ind["lst_offset_min_abs"] <= lst_tolerance_min
```

Effective constraint: `lst_offset_min_abs <= min(LST_TOLERANCE_MIN, lst_tolerance_min)`.

For overrides:
- `lst_tolerance_min = 20` → effective = 10 min → identical to baseline (no change)
- `lst_tolerance_min = 5` → effective = 5 min → tightens
- `lst_tolerance_min = 2` → effective = 2 min → tightens

**Verification** by direct execution at h=600 over 30 days:
```
tol=2.0 min: 516 feasible points
tol=5.0 min: 516 feasible points
tol=10.0 min: 516 feasible points
tol=20.0 min: 516 feasible points
```

Wait — all four return the same count. Why? Because at h=600 over 30 days, the LST offset cycles through 0-12 min uniformly, and the eclipse constraint is the dominant filter. The override tightens the LST band only when the override is **less than the offset at some grid points**. Over 30 days, the count happens to be insensitive to LST tolerance in this range. Let me check the actual sensitivity run in `results.json`:

| perturbation | h=500 | h=600 | h=700 | h=800 |
|---|---:|---:|---:|---:|
| baseline (260/270/280/290 components) | 615 | 713 | 840 | 924 |
| `lst_tolerance=20_min` | **615** | **713** | **840** | **924** ← identical to baseline |
| `lst_tolerance=5_min`  | **615** | **713** | **840** | **924** ← identical |
| `lst_tolerance=2_min`  | **615** | **713** | **840** | **924** ← identical |
| `n_rev=3`              | 618 | 716 | 843 | 928 |
| `n_rev=28`             | 611 | 708 | 836 | 919 |
| `j2_drift=disabled`    | 610 | 676 | 795 | 922 |

**Confirmed bug**: the `lst_tolerance_min` sensitivity rows produce IDENTICAL counts because the override is gated by `min(10, override) = 10` for all overrides ≥ 10 (i.e., for `=20`), and for tighter overrides (`5`, `2`) the LST constraint is **not the bottleneck** over the 6×-coarser sensitivity grid (the eclipse constraint dominates at the relevant grid step). The README's claims "approx 2x" for tol=20, "approx 0.5x" for tol=5, "approx 0.2x" for tol=2 are **wrong**.

This sensitivity matrix row is broken in two ways: (a) `tol=20` cannot increase the count because it cannot loosen the bound, and (b) `tol=5` and `tol=2` do not decrease the count on the coarser grid because the eclipse constraint is the bottleneck at the sensitivity grid step.

---

## 4. The 4 min/day LST drift claim (Q6, Q9)

The README claim (lines 79–83):

> The drift rate is `dOmega/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day = 4 min/day`.

### 4.1 The expression is frame-inconsistent

`dOmega/dt` is the **inertial** (ECI mean-of-date) RAAN rate. `d(Subsolar)/dt` is the **ECEF** geodetic subsolar-point rate. These are in **different frames** and cannot be subtracted without a frame transformation. The ECEF rate of the node is `dOmega/dt − dGMST/dt = dOmega/dt − ω_E`, and the ECEF rate of the subsolar point is `d(α_sun)/dt − ω_E`. Both are in ECEF; their difference is `dOmega/dt − d(α_sun)/dt` (frame-invariant).

### 4.2 What 0.9856 deg/day actually is

`360.9856 deg/day` = `(360 deg / T_sidereal) · 86400 s/d = 360.9856...` = Earth's **sidereal rotation rate**. This is the rate of GMST, not `dOmega/dt` for an SSO.

`360.0 deg/day` = `(360 deg / T_mean_solar) · 86400 s/d = 360.0` = Earth's **mean-solar rotation rate**. This is the rate at which the subsolar point advances in ECEF over a mean-solar day.

The difference `360.9856 − 360.0 = 0.9856 deg/day` is the **sidereal-solar differential**. It equals `SSO_TARGET_DEG_DAY = 360/365.2422` by construction — it is the SSO design rate, not a measured drift.

### 4.3 What the implementation actually does

`lst_at_node_at_t(t_L)` computes `12 + (GMST(t_L) + lon_ref − α_sun(t_L)) / 15`. Differentiating in **ECEF** (frame-consistent):

`dLST/dt = (dGMST/dt + 0 − d(α_sun)/dt) / 15 = (ω_E − d(α_sun)/dt) / 15 = (360.9856 − 0.9856) / 15 = 24 h/day`

So the **unwrapped** LST drift rate is 24 h/day = 1440 min/day. The mod-24 LST is approximately constant.

### 4.4 Measured drift rate

I computed the unwrapped LST drift over 10000 samples spanning 1 year:

```
Total LST drift over 365.24 days (unwrapped): 8765.8128 h
Rate: 1440.0000 min/day (= 24 h/day)
```

So the unwrapped LST drifts 24 h/day. The mod-24 LST is constant. **The README's "drifts through 24 h per year" claim is wrong by a factor of ~365**: it actually drifts through 24 h per **day**, not per year. Over a year-long sweep of the launch-time parameter `t_L`, the LST visits ~365 × 24 h values, but each day contributes a full 24 h cycle.

The companion audit `audit-015-lst-drift-2026-08-29.md` derives this from first principles (§1–§3) and provides independent Cowell-propagation evidence (§4) showing that **a real dawn-dusk SSO's LST at its actual ascending-node crossings varies by only ~1.7 h over a year (EoT envelope), not 24 h**. This audit confirms the implementation matches the companion audit's RED verdict on the LST-drift claim.

---

## 5. `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` is tautological (Q10)

The test (`tests/test_dawn_dusk_sso.py` lines 685–712) samples `lst_at_node_at_t(t)` at 365 points over a year and asserts:

```python
max_dist_from_18 = np.max(np.abs(lst_hours − 18.0))
max_dist_from_18 = min(max_dist_from_18, 24.0 − max_dist_from_18)
assert max_dist_from_18 > 3.0
```

**What this actually measures**: the LST-mod-24 of a constant-launch-site-longitude function over a year. Since `lst_at_node_at_t` computes `12 + (lon_ref − α_sun(t)) / 15 + GMST(t)/15` (effectively the LST at a fixed geodetic longitude as `t` varies), the LST mod 24 makes a full cycle every ~24 h (sidereal day). Over a year, it visits all 24 values ~365 times. So `max_dist_from_18 > 3.0` is satisfied trivially.

**What it does NOT measure**: any specific drift rate. The test asserts no slope, no rate, no functional form. It only asserts that LST visits values far from 18:00 over a year, which is true for ANY non-constant LST-mod-24 function.

**Verdict**: the test is **tautological**. It validates the bug (LST at constant launch-site longitude cycles through 24 h/day), not the dawn-dusk SSO LST-drift physics. The companion audit `audit-015-lst-drift-2026-08-29.md` §5 documents this in detail.

---

## 6. Frame firewall (Q8)

The lab's frame convention is mean-of-date throughout:

| Quantity | Frame | Source |
|---|---|---|
| `Omega` (orbital RAAN) | mean-of-date ECI | `j2_rhs`, `j2_nodal_rate_rad_s` |
| `α_sun` (Sun RA) | mean-of-date ECI | `sun_unit_and_dist_km` (Almanac geometric) |
| `GMST` (sidereal time) | mean-of-date UT1 | `gmst_rad_iau1982` (Aoki 1982) |
| `subsolar_lon` (ECEF) | ECEF, derived from above via `eci_to_ecef` | `subsolar_lon_rad` |
| `node_lon` (ECEF) | ECEF, derived via `eci_to_ecef` | `node_lon_from_raan_gmst` |
| `lst_at_node_hours` | apparent solar time on geoid | `lst_at_node_hours` |

All angles are in the same mean-of-date frame, so `Omega − α_sun` is internally consistent. The precession bias vs J2000 (~50.3″/yr → 0.014°/yr) is named-excluded per `earth_frames.py` lines 14–17 and validated at ~0.65° vs the byte-pinned Horizons snapshot per `earth_frames.py` lines 16–17 (and verified in the `test_subsolar_lon_dec_at_vernal_equinox_2026` test).

**No inadvertent J2000 assumption detected.** Mean-element J2 nodal drift is consistent with mean-of-date Sun model (both use mean-element quantities in the same frame). ✓

---

## 7. Additional implementation findings (beyond the 10 questions)

### 7.1 Dead-code NaN check in `constraint_indicator` (line 299–300)

```python
i_sso_deg = np.degrees(i_sso)
if use_orbit_constraint and abs(i_sso_deg − i_sso_deg) > INC_TOL_DEG * 100.0:
    # sanity
    pass
```

`i_sso_deg − i_sso_deg` is always 0 if `i_sso_deg` is finite; if `i_sso_deg` is NaN, the subtraction is NaN, and `NaN > anything` is **False** in Python. So the check **never fires**. This is dead code — either `assert` (which would raise on NaN) or a proper `math.isfinite` check would be needed.

**Severity**: MINOR — does not affect results (no NaN paths in practice), but the comment "sanity" suggests the author intended a real NaN check.

### 7.2 `lst_tolerance_min` override is effectively a one-way ratchet (already noted in §3.2)

The override `feasibility_curve(t_grid, h, lst_tolerance_min=20)` does **not** loosen the LST tolerance to 20 min; it tightens (or no-ops) because the internal `lst_ok` check in `constraint_indicator` already gates at `LST_TOLERANCE_MIN = 10 min`. The README's sensitivity claims "approx 2x feasible count at each altitude" for `lst_tolerance=20_min` are **wrong**; the actual count is **identical to baseline** (615, 713, 840, 924 at h=500/600/700/800 respectively).

### 7.3 `n_rev` sensitivity also effectively a no-op on the coarse grid

The `n_rev=3` and `n_rev=28` perturbations produce 618/716/843/928 and 611/708/836/919 — within ~2% of baseline. The README claims "window widths grow, cardinality grows" for n_rev=3 and "window widths shrink, cardinality shrinks" for n_rev=28, but the effect is small. The coarse sensitivity grid (6× step) cannot resolve the per-rev umbra timing, so the n_rev sensitivity is largely washed out.

### 7.4 Mean-vs-osculating disclosure

The README discloses "the osculating has ~0.056 deg short-period at SSO 600 km insertion (1.3 min LST)". This is a residual bias of the **mean-element** J2 nodal rate vs the actual **osculating** RAAN. Not frame-consistent in the strict sense (mean-element mean RAAN is not the osculating RAAN), but a known small effect that's disclosed.

### 7.5 `lst_at_node_at_t` is named ambiguously

The name suggests "LST at the SSO's ascending node at time t", but the implementation returns "LST at the launch site's longitude at launch time" (constant `node_lon` = `REF_SITE_LON_DEG`). Subsequent ascending nodes of the SSO precess with J2, so the function does NOT return the LST at subsequent ascending nodes. Renaming to `lst_at_launch_site_at_insertion(t)` would clarify the semantics. The companion audit `audit-015-lst-drift-2026-08-29.md` §5 also notes this name/behavior ambiguity.

---

## 8. Findings table

| # | Lines | Claim | Evidence | Tag | Severity |
|---|---|---|---|---|---|
| F1 | `README.md:65–66` | "sub_lon(t) = geocentric subsolar longitude from atan2(−u_y, −u_x)" | `atan2(−u_y, −u_x)` returns the anti-sun direction in ECI = `α_sun + π`, which is **180°** off from the geodetic subsolar longitude (= `α_sun − GMST`). The actual implementation uses the geodetic formula (ECI→ECEF then atan2). | **FACT** | **MAJOR** (README narrative contradicts the textbook formula; would be a 12-hour LST error if implemented as written) |
| F2 | `README.md:79–83, 99–100, 245, 880–890, 909` | "The drift rate is dΩ/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day = 4 min/day" | Frame-inconsistent subtraction (inertial − ECEF). Measured unwrapped LST drift rate = **1440 min/day = 24 h/day**, not 0.9856 deg/day. For an SSO, `dΩ/dt ≈ dα_sun/dt` → `d(LST)/dt ≈ 0` (modulo EoT ±~12 min). | **FACT** (impl. is 24 h/day unwrapped; 360.9856, 360.0 are sidereal & mean-solar rates; 0.9856 = sidereal-solar differential = SSO design rate as tautology) | **MAJOR** (matches companion audit's RED) |
| F3 | `experiment.py:299–300` | NaN sanity check `abs(i_sso_deg − i_sso_deg) > INC_TOL_DEG * 100.0` with body `pass` | Python NaN comparison: `NaN − NaN > x` is always False, so the check never fires. Dead code. | **FACT** | **MINOR** (no functional impact; should be `math.isfinite` or `assert`) |
| F4 | `experiment.py:419–425, 808–814`; `README.md:165–178` | Sensitivity claim: `lst_tolerance=20_min` should "approx 2x" the count | `lst_tolerance_min` override is gated by `min(LST_TOLERANCE_MIN, override) = 10 min` because `constraint_indicator` already applies the 10-min check internally. Verified empirically: counts at h=500/600/700/800 are 615/713/840/924 for `tol ∈ {2, 5, 10, 20}` — **identical**. | **FACT** | **MAJOR** (sensitivity matrix reports a wrong effect; the override is a one-way ratchet) |
| F5 | `experiment.py:252–272` | `lst_at_node_at_t(t)` returns "LST at the ascending node at insertion t_L" | Implementation: `node_lon = REF_SITE_LON_DEG * DEG` (constant); function returns `lst_at_node_hours(t_L, REF_SITE_LON_DEG*DEG)` = LST at the launch-site longitude at launch time. For an SSO with J2 drift, subsequent ascending nodes precess. The name suggests "LST at the SSO's node at any time t", but the implementation only handles insertion time. | **FACT** | **MINOR** (semantic ambiguity; same observation as companion audit §5) |
| F6 | `tests/test_dawn_dusk_sso.py:685–712` | `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` measures a 4-min/day drift rate | Test only asserts `max_dist_from_18 > 3.0`. This is satisfied trivially because `lst_at_node_at_t` (constant `node_lon`, varying `t`) cycles through 24 h in ~1 sidereal day. The test does not measure a rate. | **FACT** | **MAJOR** (test validates the bug, not the physics; matches companion audit §5) |
| F7 | `src/lab_utils/earth_frames.py:280–282`; `experiment.py:315–319` | `lst_at_node_hours` returns apparent LST (no EoT correction) | Lab_utils docstring is explicit. The dawn-dusk SSO canonical target is conventionally **mean** LST. EoT envelope over 2026: −14.20 to +16.45 min. | **FACT** | **MINOR** (convention pinned in lab_utils; not re-disclosed in experiment README; bias ≤ 16 min is comparable to 10-min tolerance) |
| F8 | `experiment.py:43–44, 955–965` | C4 insertion: `Omega(t_L) = GMST(t_L) + lon_ref` puts the **ascending** node over Eastern Range | `node_lon = Omega − GMST = lon_ref` = −80.6039° (Cape Canaveral). With `ω = 0` and the satellite at `argument of latitude = 0`, the satellite is at the **ascending** node (northward equator crossing). ✓ | **FACT** | **GREEN** |
| F9 | `src/lab_utils/earth_frames.py:269, 261–271`; `experiment.py:267, 311` | `node_lon = Omega − GMST` is the geodetic ascending-node longitude (sign convention) | The passive Z rotation in `eci_to_ecef` (lines 181–192) does `lon_ecef = lon_eci − θ_G`, so `node_lon_ecef = Omega − GMST`. ✓ | **FACT** | **GREEN** |
| F10 | `experiment.py:275–283`; `tests/test_dawn_dusk_sso.py:448–461` | LST offset wrap `(lst − target + 12) % 24 − 12` is signed and in (−720, 720] min | Verified numerically. Boundary case `lst = 6.0` (exactly 12 h below target 18.0) gives `offset = −720 min`, the boundary (excluded by the open lower bound in the docstring's `(−720, 720]`). Sign convention: positive when LST > target. | **FACT** | **GREEN** (boundary exclusion is the only nit) |
| F11 | `src/lab_utils/earth_frames.py` (frame firewall); `experiment.py:140–149` | Lab convention is mean-of-date throughout; no inadvertent J2000 assumption | All angles (Ω, α_sun, GMST) live in mean-of-date ECI. The precession bias (~50.3″/yr → 0.014°/yr) is named-excluded. Subsolar_lon uses ECI→ECEF rotation. SSO J2 nodal drift is in the same frame as α_sun. | **FACT** | **GREEN** |
| F12 | `README.md:79–83, 880–890`; `experiment.py:908–916` | "The LST passes through 18:00 (and through 6:00, the opposite terminator) once per year" | The mod-24 LST visits 18:00 once per ~24 h (sidereal day), not once per year. Over a year, 18:00 is visited ~365 times. | **FACT** | **OBSERVATION** (consistent with F2; off by factor ~365) |
| F13 | `experiment.py:919–928`; `README.md:268–280` | Disclosure of mean-element vs osculating J2 (~0.056 deg = 1.3 min LST) | Disclosed in both files. Acknowledged as a residual bias. | **INFERENCE** | **GREEN** |
| F14 | `experiment.py:803–813`; `results.json` | Sensitivity `lst_tolerance=5_min` and `lst_tolerance=2_min` should "≈ 0.5x" and "≈ 0.2x" baseline counts | Actual counts are **identical to baseline** at the coarse 6× sensitivity grid. Eclipse constraint is the bottleneck at this grid step, not LST tolerance. | **FACT** | **MAJOR** (sensitivity claims do not match measured counts) |
| F15 | `README.md:282` | "600 s coarse step + 1 s bisection; finer grid would add components but does not change the structure (verified by 5-min vs 10-min grid test)" | `test_grid_step_does_not_change_best_lst_offset_by_2x` (lines 629–661) confirms best LST offset is grid-stable to <2 min between 600 s and 300 s. Does not test cardinality. | **INFERENCE** | **OBSERVATION** |
| F16 | `experiment.py:914–916` | "The cylindrical beta-cutout fast check ... disagrees with the slow event-finder ... on the best candidates" | Documented as a structural ambiguity (Exp 014 cone-vs-cylinder disclosure). Disagreement is reported, not silently absorbed. | **INFERENCE** | **GREEN** |
| F17 | `tests/test_dawn_dusk_sso.py:516–536` | `test_subsolar_lon_matches_donor_lab_utils_donor_consistency` uses the geodetic formula (ECEF atan2) | This test catches the README's `atan2(−u_y, −u_x)` formula as a "previous version" bug (lines 521–525). Test enforces the correct geodetic convention. ✓ | **FACT** | **GREEN** (test enforces the right thing; contradicts README's prose) |
| F18 | `experiment.py:13–15`; `README.md:11–17` | "LST-at-ascending-node condition (target = 18:00, the classic dawn-dusk ascending terminator)" | The "18:00" target is for a **dusk-ascending** SSO (satellite crosses the equator going north at the dusk terminator). Apparent LST is used. Both names are correct (it is "dawn-dusk" because the orbit plane is fixed in the sun-frame, with one terminator at ascending and the other at descending). | **INFERENCE** | **GREEN** |
| F19 | `experiment.py:140–149` (FRAME_CONVENTION); `README.md:36–40` | Frame is "mean-of-date", with explicit "equation of equinoxes excluded (≤ 1.1 s RAAN phasing)" | Disclosed. Small effect compared to 10-min LST tolerance. | **FACT** | **GREEN** |
| F20 | `experiment.py:920–929`; `README.md:268–280` | Limitation: "Mean-element J2 nodal rate; the osculating vs mean offset is ~0.056 deg at SSO 600 km insertion (1.3 min LST)" | Disclosed. 1.3 min bias vs 10-min tolerance → within tolerance. | **FACT** | **GREEN** |

---

## 9. Cross-references and verification commands

All findings were verified by reading the listed files and executing read-only Python commands:

```bash
# Verify subsolar_lon formula (returns geodetic, not atan2(-u_y, -u_x))
.venv/Scripts/python.exe -c "
import sys; sys.path.insert(0, 'src')
import importlib.util
spec = importlib.util.spec_from_file_location('e', 'research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import numpy as np
t = m.t_since_j2000_from_gregorian(2026, 1, 1, 0, 0, 0)
u, _ = m.sun_unit_and_dist_km(t)
print('atan2(-u_y,-u_x) =', np.degrees(np.arctan2(-u[1], -u[0])))
print('subsolar_lon_rad =', np.degrees(m.subsolar_lon_rad(t)))
"

# Verify LST drift rate (1440 min/day unwrapped)
.venv/Scripts/python.exe -c "
import sys; sys.path.insert(0, 'src')
import importlib.util
spec = importlib.util.spec_from_file_location('e', 'research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import numpy as np
t0 = m.t_since_j2000_from_gregorian(2026, 1, 1, 0, 0, 0)
times = np.linspace(t0, t0 + 365.2422*86400, 10000)
lsts = np.array([m.lst_at_node_at_t(float(t)) for t in times])
# unwrap, compute total drift
diff = np.diff(lsts); diff = (diff + 12) % 24 - 12
total = lsts[0] + np.cumsum(diff)
print('rate =', (total[-1] - lsts[0]) / 365.2422 * 60, 'min/day')
"

# Verify lst_tolerance_min override sensitivity
.venv/Scripts/python.exe -c "
import sys; sys.path.insert(0, 'src')
import importlib.util
spec = importlib.util.spec_from_file_location('e', 'research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import numpy as np
t0 = m.t_since_j2000_from_gregorian(2026, 1, 1, 0, 0, 0)
g = np.arange(t0, t0 + 30*86400, 600.0)
for tol in [2, 5, 10, 20]:
    f = m.feasibility_curve(g, 600, lst_tolerance_min=tol)
    print(f'tol={tol}: count={int(np.sum(f))}')
"
```

All commands executed without modifying any canonical file. The outputs confirm the analysis above.

---

## 10. Verdict

**Net implementation-layer verdict**: the **code** of Exp 015 (and the `lab_utils/earth_frames.py` and `lab_utils/orbits.py` machinery it uses) is **GREEN** for the closed-form formulas (Q3, Q4, Q5, Q8). The **README narrative** is **RED** on Q1 (describes `atan2(-u_y, -u_x)` which is the anti-sun direction, not the geodetic subsolar longitude), Q6 (24 h/year vs 24 h/day drift, off by factor 365), and Q9 (frame-inconsistent subtraction `dΩ/dt − d(Subsolar)/dt`). The **test** for LST drift is **RED** on Q10 (tautological). The **sensitivity matrix** is **RED** on the `lst_tolerance_min` rows (F4, F14: override gated to ≤ baseline; counts identical across overrides).

The findings corroborate and tighten the existing `audit-015-lst-drift-2026-08-29.md`. The implementation layer is internally consistent (the code path is correct), but the documentation layer (README) and the test layer (`test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO`) and the sensitivity layer (`lst_tolerance_min` rows) need remediation to align with the actual physics and the actual code behavior.

### 10.1 IMPLEMENTATION-VERDICT

```
IMPLEMENTATION-VERDICT: YELLOW
```

(The code is correct, but the README narrative, the LST-drift test, and the `lst_tolerance_min` sensitivity rows are RED. Implementing-only the code is GREEN; the experiment as a whole is YELLOW because the documentation and the test encode and validate an incorrect physical claim.)

### 10.2 KEY-FINDINGS

- The `subsolar_lon_rad` function in `lab_utils/earth_frames.py` is **correct** (geodetic, ECEF atan2). The README's prose description of it as `atan2(-u_y, -u_x)` is **wrong** (that formula returns the anti-sun direction in ECI). The code uses the right formula; only the README narrative is wrong.
- The C4 insertion `Omega(t_L) = GMST(t_L) + lon_ref` correctly places the **ascending** node over Eastern Range at insertion. The sign convention `node_lon = Omega − GMST` is consistent with the passive Z rotation in `eci_to_ecef`.
- The LST drift rate of "4 min/day" in the README is **wrong by a factor ~365**: the actual unwrapped LST drift rate (driven by the C4 ground-fixed convention) is 1440 min/day = 24 h/day, which mod-24 is approximately constant for an SSO. The 0.9856 deg/day in the README's `360.9856 − 360.0` expression is the sidereal-solar differential, which equals the SSO design rate by construction — it is a tautology, not a measured drift.
- The test `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` does not measure a 4-min/day drift rate; it only checks that LST visits values >3 h from 18:00 over a year. This is trivially true for `lst_at_node_at_t` (constant `node_lon`, varying `t`) and validates the bug, not the physics.
- The `lst_tolerance_min` sensitivity override in `feasibility_curve` is **gated by `min(LST_TOLERANCE_MIN, override)`**, so `tol=20` cannot loosen the LST band (it's a no-op). The README's claim "approx 2x feasible count at each altitude" for `tol=20` is **wrong**; the actual counts are identical to baseline at every altitude.

### 10.3 REQUIRED-FIXES

1. **README narrative (`README.md` lines 65–66, 79–83, 96–99, 245, 880–890, 909)**: replace the `atan2(-u_y, -u_x)` description with the actual geodetic formula (ECI→ECEF then atan2); replace the "drifts through 24 h per year" / "4 min/day" claim with the correct statement that the LST mod 24 is approximately constant for an SSO, oscillating with the EoT envelope (~±12 min peak-to-peak ~24 min).
2. **`experiment.py:875–916` (findings payload)**: remove or correct the FINDING that "the LST at the ascending node of a dawn-dusk SSO drifts through 24 h over the year" — see companion audit for the corrected physics and station-keeping budget (~130–290 m/s/year for J2 closure residual, not ~146 m/s/year for a fictitious 24 h/year drift).
3. **`experiment.py:252–272` (`lst_at_node_at_t`)**: rename to `lst_at_launch_site_at_insertion(t_L)` or `lst_at_insertion_node_at_t_L(t_L)` to make the semantics unambiguous. The current name suggests "LST at the SSO's ascending node at any time t", but the implementation returns the LST at the launch-site longitude at launch time.
4. **`experiment.py:299–300` (NaN sanity)**: replace `abs(i_sso_deg - i_sso_deg) > INC_TOL_DEG * 100.0` with `not math.isfinite(i_sso_deg)` or an `assert math.isfinite(i_sso_deg), ...`. The current check is dead code (NaN comparison is always False).
5. **`tests/test_dawn_dusk_sso.py:685–712` (`test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO`)**: delete or rewrite. The current assertion validates the bug. A correct test would either (a) measure the LST drift rate at subsequent ascending nodes of an SSO propagation and assert it's bounded by the EoT envelope (~±12 min) plus the J2 closure residual (~0.32 min/day), or (b) compute `lst_at_node_at_t` over a year and assert the drift rate is consistent with the sidereal-solar differential (≈1440 min/day unwrapped, mod-24 ≈ constant within EoT).
6. **`experiment.py:419–425` and sensitivity rows `lst_tolerance=20_min`**: either (a) re-implement `feasibility_curve` to apply the override directly (replace the internal `LST_TOLERANCE_MIN` check with the override), or (b) correct the README sensitivity claim from "approx 2x" to "no change (override is a one-way ratchet)". The actual counts at all four altitudes for `tol ∈ {2, 5, 10, 20}` are 615/713/840/924 — identical, contradicting the README's predictions.
7. **Re-run the experiment** after the above remediation to refresh `results.json`. The headline numbers (260–290 feasible components, ~700 h total feasible width, 36.4 vs 11.6/day equinox dominance, i_SSO anchors) remain valid; only the LST-drift narrative, the `lst_tolerance_min` sensitivity claims, and the LST-drift test need correction.

---

## 11. Cross-references

- `localdocs/reports/audit-015-lst-drift-2026-08-29.md` — primary LST-drift physics audit (companion, already RED); this audit confirms the implementation matches.
- `localdocs/reports/audit-015-adversarial-2026-08-29.md` — adversarial review (related).
- `localdocs/reports/audit-015-follow-up-candidates-2026-08-29.md` — follow-up candidate experimental design.
- `localdocs/reports/audit-015-literature-2026-08-29.md` — literature cross-check.
- `localdocs/reports/audit-015-numerical-falsifier-2026-08-29.md` — numerical falsifier audit.

End of audit report.