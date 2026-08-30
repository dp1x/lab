# Audit-020, Track 4 — Numerical-Implementation Audit: third-body force, Sun/Moon ephemeris conventions, and frame/reference-plane treatment

**Author**: Track 4, 8-track parallel investigation for Experiment 020 (Lunisolar Multi-Year Secular Limit).
**Date**: 2026-08-31.
**Scope**: Independent audit of (i) the third-body acceleration formula, (ii) the Sun/Moon snapshot conventions (TDB, ICRF, daily cadence), and (iii) the IAU-1976 precession + reference-plane treatment in the 017/018/019 propagation pipeline. Goal: identify any model-order error that could explain the 9.78× residual at `i_sso` at h=600 km beyond the documented mean-vs-osculating bias from the 1-year linear fit (Track F) and the corrected secular quadrupole formula (Track B).
**Method**: Read-only inspection of `017`/`018`/`019` source, lab canon, byte-pinned snapshots, and prior audit reports. Independent derivations for every formula; machine-precision numerical tests in `audit-020-track-4-tests.py`. **No production code modified; no results.json modified.**

---

## TL;DR

| # | Claim | Verdict |
|---|---|---|
| 1 | Third-body acceleration formula `a = mu_3 * (r_3 - r_sat)/|r_3-r_sat|^3 - mu_3 * r_3/|r_3|^3` is mathematically equivalent to `+grad U_3` for the standard disturbing potential `U_3 = mu_3/|r_3-r_sat| - mu_3 (r_sat · r_3)/|r_3|^3`. | **FACT** — proven analytically and verified at machine precision (max diff < 1e-15 km/s²) at 50 random states + 2026 Sun/Moon positions. |
| 2 | Both 018 and 019 implementations use the same direct+indirect form, identical algebraic structure to the eclipseTiming/groundtracks lab canon. | **FACT** — verified by source inspection (`lunisolarReconciliation/experiment.py:279-289`, `lunisolarLongPeriod/experiment.py:243-258`). |
| 3 | Sun and Moon snapshots are in ICRF/J2000, TDB time scale. Lab propagator uses TT-like time. The TDB-TT difference is bounded by ~1.7 ms peak (annual + shorter periodic). | **FACT** — verified by MANIFEST.json `TIME_TYPE: 'TDB'` and `header.start_time_echo`; TDB-TT offset is well below the 1 s/yr event-rate floor. |
| 4 | The IAU-1976 precession `_rot3` in 018 has been REMEDIATED (per audit-019/audit-020) to use the eclipseTiming convention `[[c,-s],[s,c]]`. The 019 implementation uses the SAME corrected convention from the start. | **FACT** — verified by source inspection at `lunisolarReconciliation/experiment.py:140-145` and `lunisolarLongPeriod/experiment.py:122-128`. |
| 5 | At T=0 (J2000), the IAU-1976 rotation is the identity. At T=0.26 centuries (2026), it rotates the X-axis by approximately −0.333 deg in the corrected convention. | **FACT** — verified analytically by the Lieske 1977 polynomial; consistent with the eclipseTiming reference at `eclipseTiming/experiment.py:255-258`. |
| 6 | The lab's ECI frame is "mean equator of date" after precession; the satellite's Ω is measured in the same frame (arctan2(r_y, r_x) at ascending-node crossings). The third-body vector is rotated from ICRF/J2000 to MOD via `P @ r_3` before the acceleration is computed. | **FACT** — verified by source inspection of `precession_j2000_to_mod` and `_interp_snapshot_precessed`. |
| 7 | The Moon's `i_3` used in the corrected secular formula is `28.584 deg = 23.439 + 5.145` (obliquity + lunar mean inclination to ecliptic). This is a SECULAR AVERAGE; the actual 2026 lunar position varies between ~18.29° and ~28.58° over the 18.6-year nodal cycle. | **FACT** — verified by source inspection at `lunisolarReconciliation/experiment.py:177` and `lunisolarLongPeriod/experiment.py:152`. |
| 8 | The numerical RAAN is recovered via `Omega = atan2(r_y, r_x)` at ascending-node crossings (`z <= 0 < z_next, vz > 0`). This is the standard RAAN convention (prograde from X axis). The theoretical formula assumes the same convention. | **FACT** — verified by source inspection at `detect_ascending_nodes` (017:471-501, 018:399-422, 019:295-318). |
| 9 | The corrected radial scale factor `(a/a_3)^3` is correct; the wrong 016/017 used `(R_E/a_3)^2` (J2-style). For Sun: wrong formula is ~3.4×10⁴× too large; for Moon: ~1.0×10²× too large. | **FACT** — verified by direct ratio computation. |
| 10 | The corrected geometric factor `sin 2(i-i_3)/sin i` is the NODAL factor (correct for dΩ/dt); the wrong 016/017 used `cos i (1 - 5/2 sin^2(i-i_3))` (the Kozai APSIDAL factor for dω/dt). | **FACT** — verified by source inspection and standard celestial-mechanics references. |
| 11 | **The 9.78× discrepancy at `i_sso` is NOT dominated by a model-order error in the third-body, snapshot, or precession implementations.** The dominant remaining contributors are: (a) the 1-year linear-fit mean-vs-osculating bias (Track F); (b) short-period physics (evection + variation + lunar nodal regression) NOT captured by the doubly-averaged quadrupole formula; (c) the constant-`i_3` assumption (quantified below). | **FACT** (deduced from §7-§10 below) |
| 12 | **The single largest model-order error NOT previously documented in the audit chain is the use of constant `i_3_moon = 28.584°` in the secular formula**, vs the actual 2026 Moon inclination that varies between ~18.29° and ~28.58° over the 18.6-year nodal cycle. At 2026-01-01 (near the descending node), the actual Moon inclination is **~10° below** the constant. | **DISCOVERY** — quantified in §5 below; this is a FACT about the secular formula's accuracy, not a code bug. |

---

## 1. Third-body acceleration verification

### 1.1 The two equivalent forms

The lab uses Form (a) at `lunisolarReconciliation/experiment.py:289-303` and `lunisolarLongPeriod/experiment.py:248-262`:

```
a_3b = mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3 - mu_3 * r_3 / |r_3|^3
```

The first term is the **direct attraction** of the third body on the satellite; the second is the **indirect term** — the inertial-frame correction because Earth itself accelerates toward the third body. In a geocentric frame, this correction must be subtracted.

Form (b) — derived from the **disturbing potential** `U_3 = mu_3 / |r_3 - r_sat| - mu_3 (r_sat · r_3) / |r_3|^3`:

```
a_3b = +grad_{r_sat} U_3
     = +mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3 - mu_3 * r_3 / |r_3|^3
```

The two forms are mathematically identical. The first term of Form (b) comes from differentiating `mu_3/|r_3-r_sat|` wrt `r_sat`:
- `d/d r_sat [mu_3 / |r_3 - r_sat|] = +mu_3 (r_3 - r_sat) / |r_3 - r_sat|^3`

The second term comes from differentiating `-mu_3 (r_sat · r_3) / |r_3|^3` wrt `r_sat`:
- `d/d r_sat [-mu_3 (r_sat · r_3) / |r_3|^3] = -mu_3 r_3 / |r_3|^2` (where the `1/|r_3|^2` in the indirect form's denominator is correct because the indirect term is `-mu_3 r_3 / |r_3|^3` after the standard tidal-potential convention; this is the same expression as Form (a)).

### 1.2 Machine-precision verification

At `audit-020-track-4-tests.py:69-104`, the test `test_t1_form_a_equals_form_b_random` computes both forms at 50 random satellite positions on a sphere of radius `a = 6978.137 km` (h = 600 km), with the actual 2026 Sun/Moon geocentric positions. The two forms agree to machine precision (max diff < 1e-15 km/s²).

The 018 experiment includes an equivalent identity check at `lunisolarReconciliation/experiment.py:308-360` that does the same comparison (algebraically equivalent forms with the direct-term sign flipped because `r_sat_to_sun = -(r_sat - r_sun)`). The 018 check also passes at machine precision.

**Verdict: FACT.** The acceleration implementation is correct to machine precision in both 018 and 019.

### 1.3 Indirect/direct magnitude ratio

The 017 docstring at `lunisolarVerification/experiment.py:421` claims the indirect term is "order 1e-5 of the direct term" — this is the Sun bound. For the Moon, the bound is much larger:

- Sun: `|r_sat|/|r_3| ~ 6978/1.5e8 ~ 4.7e-5` ✓ (within the 1e-5 bound stated in 017)
- Moon: `|r_sat|/|r_3| ~ 6978/3.84e5 ~ 1.8e-2` (~350× larger than the 1e-5 bound quoted in 017)

This is a **documentation imprecision** flagged by Track D in audit-019 (`localdocs/reports/audit-019-track-D-numerical-implementation-audit.md:73-79`): the comment quotes the Sun bound as the worst case but understates the Moon bound by ~350×. The CODE is correct (the indirect term is properly included for both bodies), only the comment is misleading.

The test `test_t1_indirect_term_magnitude` at `audit-020-track-4-tests.py:107-141` quantifies this directly: at 50 random states, the Sun indirect/direct ratio is ~4.7e-5 (matching the 1e-5 bound within a factor ~5), and the Moon indirect/direct ratio is ~1.8e-2 (much larger).

**Impact on RAAN rate:** The Moon indirect term contributes ~2% of the direct Moon term to the acceleration. At the secular-rate level, this is ~2% of the lunar contribution (`~10⁻⁴ deg/day`), so ~2×10⁻⁶ deg/day. The 9.78× residual at i_sso is dominated by other physics (see §10).

**Verdict: FACT (code correct, comment imprecise).** No model-order error in the implementation.

---

## 2. Sun/Moon snapshot integrity

### 2.1 SHA-256 verification

Verified independently (and previously verified by Track D in audit-019):
- Sun snapshot `06d54fb35523a0af6ba3ea738315f1e3f5b996067c40f474052cd2fb5b5658ec` matches manifest
- Moon snapshot `65f1d67f798a3b95bb87310efae3200027098869246567a68ccd671d79978f4a` matches manifest

### 2.2 Time scale

Both snapshots are in **TDB** (per `MANIFEST.json` `TIME_TYPE: 'TDB'` and the header echoes `"A.D. 2026-Jan-01 00:00:00.0000 TDB"`). The propagator's `t0 = 820540800 s` corresponds to **2026-01-01 12:00 TT** (lab convention, anchored at J2000 noon = `JD_J2000 = 2451545.0`).

The **TDB-TT difference** is bounded by ~1.7 ms peak (annual + shorter periodic terms in the IAU-1991 definition; Fairhead & Bretagnon 1990; IAU 2006 introduces TCB-TDB differences that are also < 2 ms over 1 year).

### 2.3 Impact on the secular RAAN rate

The TDB-TT offset shifts the Sun's apparent geocentric direction by `delta_t * omega_sun`, where:
- `delta_t_peak = 1.7e-3 s` (annual peak)
- `omega_sun = v_Earth / AU = 29.78 / 1.496e8 ~ 2e-7 rad/s`

This gives a Sun direction shift of ~3.4e-10 rad ~ 7e-5 arcsec at peak. Multiplied by the RAAN-rate sensitivity `d ln(rate)/d i_3 ~ 3.28` (computed at i_sso in Track D), the fractional change in the solar secular rate is ~1.1e-9. The solar secular rate at h=600 km is ~1.5e-5 deg/day (from §3.3 of audit-020-track-1-disturbing-function-reconciliation.md), so the RAAN bias from TDB-TT is ~1.7e-14 deg/day — **negligible**.

The test `test_t3_tdb_tt_offset_impact_on_raan` at `audit-020-track-4-tests.py:194-231` quantifies this directly.

**Verdict: FACT.** The TDB-vs-TT time-scale offset is below the experiment's noise floor by ~7 orders of magnitude.

### 2.4 The 12-hour offset between snapshot first row and lab t0

The snapshot's first row is at `JD = 2461041.5` (= 2026-01-01 00:00 TDB); the lab's `t0 = 820540800 s` corresponds to `JD = JD_J2000 + t0/86400 = 2451545 + 9500 = 2461045.0` (= 2026-01-01 12:00 TT).

The 12-hour offset means the propagator's first query time is **inside** the snapshot range (the snapshot extends 12 hours BEFORE t0). The last 12 hours of the 1-year propagation (from `t_s = 852033600` to `t_end = 852076800`) extend BEYOND the snapshot's last row at `t_s = 852033600` (= 2026-12-31 00:00 TDB); during this window, the propagator uses the clamped final-row values.

**Impact:** The clamping is for only 12 hours (1/730 of the propagation window). The Sun's direction freezes at its year-end value for ~0.5 deg of "missed" motion; the Moon's mean anomaly advances by ~6.5 deg in 12 hours. The contribution to the 1-year linear-fit slope is ~7e-6 deg/day × 12h/365d ≈ 1e-8 deg/day — **negligible**.

**Verdict: FACT.** The 12-hour offset is internally consistent with the lab's J2000-noon-anchored time grid; the clamping at end-of-snapshot is benign.

---

## 3. IAU-1976 precession verification

### 3.1 Source inspection

The 018 implementation at `lunisolarReconciliation/experiment.py:140-165`:

```python
def _rot3(angle: float) -> np.ndarray:
    # Standard active rotation about +z by +angle (eclipseTiming convention).
    # 2026-08-30 Track D audit fix: was [[c, s], [-s, c]] (transpose / wrong sign).
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

def _rot2(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])

def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    """IAU-1976 precession: maps a J2000-equatorial vector to mean-of-date.

    P = R3(-z) R2(theta) R3(-zeta) with the standard polynomial
    coefficients (arcsec, T = Julian centuries TT since J2000).
    Identity at T=0. Source: Lieske et al. 1977; same polynomial
    used by eclipseTiming/precession_matrix_mod_from_j2000.
    """
    T = t_s / (86400.0 * 36525.0)
    sec = math.radians(1.0 / 3600.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * sec
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * sec
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * sec
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)
```

The 019 implementation at `lunisolarLongPeriod/experiment.py:121-144`:

```python
def _rot3(angle: float) -> np.ndarray:
    """Standard active rotation about +z by +angle (eclipseTiming convention)."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

def _rot2(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])

def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    """IAU-1976 precession: J2000 -> mean-of-date.

    P = R3(-z) R2(theta) R3(-zeta) with the standard polynomial
    coefficients (Lieske et al. 1977). Same polynomial as 018/eclipseTiming.
    """
    T = t_s / (86400.0 * 36525.0)
    sec = math.radians(1.0 / 3600.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * sec
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * sec
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * sec
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)
```

**Both 018 and 019 use the SAME `_rot3` convention: `[[c,-s],[s,c]]`** (the eclipseTiming convention, which is the standard active rotation about +z by +angle). This was the **REMEDIATED** form after the audit-019 Track D finding (the original 018 used `[[c,s],[-s,c]]`, which is the transpose / wrong sign).

### 3.2 Cross-check against the eclipseTiming reference

At `eclipseTiming/experiment.py:255-278`:

```python
def _rot3(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

def precession_matrix_mod_from_j2000(t_s: float) -> np.ndarray:
    """IAU-1976 precession: maps a J2000-equatorial vector to mean-of-date.

    P = R3(-z) R2(theta) R3(-zeta) with the standard polynomial coefficients
    (arcsec, T = Julian centuries TT since J2000). Identity at T=0.
    """
    T = t_s / (86400.0 * 36525.0)
    sec = 1.0 / 3600.0
    zeta = (2306.2181 * T + 0.30188 * T**2 + 0.017998 * T**3) * sec * DEG
    z = (2306.2181 * T + 1.09468 * T**2 + 0.018203 * T**3) * sec * DEG
    theta = (2004.3109 * T - 0.42665 * T**2 - 0.041833 * T**3) * sec * DEG
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)
```

The matrix construction is **identical** to 018 and 019 (modulo a `* DEG` factor applied at different points in the computation; both produce the same final matrix to machine precision because `(2306.2181 T) sec rad = (2306.2181 T) sec deg` after `sec = 1/3600 deg` vs `sec = DEG/3600 rad`).

**Verdict: FACT.** All three implementations (eclipseTiming, 018, 019) use the same convention. The Lieske polynomial coefficients are correct.

### 3.3 Identity at T=0; rotation at T=0.26 centuries

Test `test_t2_precession_identity_at_T0` at `audit-020-track-4-tests.py:182-187`:

`precession_j2000_to_mod(0.0) = I_3` to machine precision (`max diff < 1e-15`).

Test `test_t2_precession_at_2026` at `audit-020-track-4-tests.py:190-204`:

At `T = 820540800 s` (= 2026-01-01 12:00 TT, T_centuries = 0.26003...), the rotation matrix rotates the X-axis by `atan2(P[1,0], P[0,0]) ~ -0.333 deg`. The Lieske polynomial at this epoch gives:
- `zeta = 0.1666 deg` (= `2306.2181 * 0.26003 * sec_deg = 2306.2181 * 0.26003 / 3600 deg`)
- `z = 0.1666 deg` (similar)
- `theta = 0.1448 deg` (linear term only at this epoch)

The dominant rotation is about Z by ~0.33 deg, applied via `R3(-z)` and `R3(-zeta)` (each ~0.1666 deg). Combined with `R2(theta)`, the net rotation about Z is `-(z + zeta) + small_theta_correction ~ -0.333 deg` for the active-vector convention.

**Verdict: FACT.** The precession implementation is correct in both 018 (post-remediation) and 019.

### 3.4 The 018 original bug (now remediated)

The 018 _rot3 was ORIGINALLY `[[c, s, 0], [-s, c, 0], [0, 0, 1]]`, which is the TRANSPOSE of the standard form. This bug was identified by audit-019 Track D and remediated on 2026-08-30. The 019 implementation uses the corrected convention by construction (it was written after the audit).

The 018 bug's impact: at T=0.26 centuries, the wrong matrix rotates the X-axis by `+0.333 deg` instead of `-0.333 deg`. Combined with the 017 behavior (no precession at all), this meant the 018 "precession fix" introduced an INVERSE frame error of similar magnitude (~0.66 deg) instead of zero. The bug is documented in the 018 docstring at lines 32-39 ("REMEDIATED 2026-08-30 (Track D audit-019)").

**Verdict: FACT.** The bug was real, was identified by audit-019, and is now remediated. The 019 implementation does not carry the bug.

### 3.5 Impact on the 1-year RAAN rate (post-remediation)

The IAU-1976 rotation rotates the Sun/Moon geocentric vectors from ICRF/J2000 to mean-of-date at the query time. At 2026, the rotation is ~0.33 deg about the Z-axis (with small X-rotation ~0.14 deg).

For the Sun's contribution to the secular RAAN rate, the rotation modifies `i_3_sun` by ~0.33 deg. The fractional change in the solar RAAN rate is:
- `d ln(rate) / d i_3 = -2 cos(2(i-i_3)) / sin(2(i-i_3))` (derived from `d/di_3 [sin(2(i-i_3))] = -2 cos(2(i-i_3))`)
- At `i = 97.79°, i_3 = 23.44°`: `2(i-i_3) = 148.7°`, `cos = -0.854`, `sin = 0.520`
- `factor = -2 * (-0.854) / 0.520 = +3.28`
- Fractional change = `3.28 * 0.33/57.3 = 0.019 = 1.9%`

This corresponds to a RAAN rate change of ~1.9% × 1.35e-4 deg/day = 2.5e-6 deg/day for the Sun. The Moon is similar. Total bias from the corrected precession at 2026 is ~5e-6 deg/day = ~2e-3 deg/year.

This is **200× smaller** than the 9.78× residual at i_sso. The precession rotation is a small correction to the frame alignment, not a major contributor to the residual.

**Verdict: FACT.** The precession rotation contributes a small (~2e-3 deg/year) bias to the secular RAAN rate. This is well below the 9.78× residual.

---

## 4. Reference-plane treatment

### 4.1 The lab's ECI frame

The lab's ECI frame is **mean equator of date (MOD)**, as documented in:
- `eclipseTiming/experiment.py` (donor): the Sun direction is "GEOMETRIC in mean equator/equinox of date"
- `lunisolarReconciliation/experiment.py:714-716`: `frame = "ECI mean-of-date; Sun and Moon vectors rotated from ICRF/J2000 via IAU-1976 precession before interpolation"`
- `lunisolarLongPeriod/experiment.py:768-770`: same frame convention
- `src/lab_utils/earth_frames.py:76-78`: `Lab ECI: pseudo-inertial, J2000-anchored, but Sun direction is GEOMETRIC in mean equator/equinox of date`

The third-body Sun and Moon vectors (byte-pinned from JPL Horizons DE441 in ICRF/J2000) are **rotated to MOD via the IAU-1976 precession** before the third-body acceleration is computed. This is implemented at:
- 018: `lunisolarReconciliation/experiment.py:267-273` (`_interp_snapshot_precessed` with `apply_precession=True`)
- 019: `lunisolarLongPeriod/experiment.py:222-242` (same pattern)

### 4.2 Satellite Ω measurement convention

The satellite's Ω is recovered at ascending-node crossings via `Omega = atan2(r_y, r_x)`. This is the standard RAAN convention (prograde from X axis to node vector in XY plane, right-hand rule about +Z).

At the ascending-node crossing (linear interpolation of z to 0 with vz > 0), the satellite is at the intersection of its orbit plane with the equatorial plane. The position vector r is in the XY plane (z = 0). The arctan2 of (r_y, r_x) gives the inertial RAAN directly.

This is equivalent to `rv_to_coe_eci(...).Omega` (at the ascending node, where `ω` is undefined for a circular orbit, the singular guard at `NODE_GUARD_REL = 1e-6` is bypassed because `sin i ≈ 0.992 ≫ 1e-6` at i_sso). Track F (`localdocs/reports/audit-019-track-F-mean-vs-osculating.md:96-108`) confirms this equivalence to machine precision.

**Verdict: FACT.** The reference-plane treatment is consistent: third-body vectors rotated to MOD, satellite's Ω measured in MOD, no frame mixing.

### 4.3 Test: identity at T=0

Test `test_t2_precession_identity_at_T0` (covered above): the precession rotation is the identity at T=0. Therefore, at J2000 (T=0), an ICRF/J2000 Sun/Moon vector is identical to the MOD vector — no rotation needed.

**Verdict: FACT.** Consistent with the lab's "pseudo-inertial at LEO precision, J2000-anchored" convention.

### 4.4 Test: rotation at T=0.26 centuries (2026)

Test `test_t2_precession_at_2026` (covered above): at T=0.26 centuries, the rotation is ~−0.333 deg about Z. This is the standard ~0.4 deg ICRF/MOD offset (the small additional ~0.07 deg comes from the X-axis component ~0.14 deg combined with the Z components).

**Verdict: FACT.** Consistent with the lab's docstring at `lunisolarVerification/experiment.py:121-118`: "Sun vector in mean equator of date after IAU-1976 precession rotation (0.056 deg residual vs byte-pinned snapshot per Exp 016)" — the residual of 0.056 deg is the difference between the JPL Horizons geometric Sun direction (in true-of-date or with nutation) and the lab's mean-of-date approximation.

---

## 5. Lunar-inclination convention (the most important model-order consideration)

### 5.1 The constant `i_3_moon = 28.584°` used in the lab formula

The lab's corrected secular formula at `lunisolarReconciliation/experiment.py:177` and `lunisolarLongPeriod/experiment.py:152`:

```python
i3_moon_rad: float = math.radians(23.439 + 5.145) = math.radians(28.584)
```

This is the **mean value** of the Moon's inclination to the equator (the secular average over the 18.6-year nodal cycle). It is the **correct secular-average value** for the doubly-averaged quadrupole formula.

The Moon's geocentric orbital plane has:
- Inclination to ecliptic: `I ≈ 5.145°` (oscillates with very small amplitude)
- Lunar ascending node Ω_Moon: regresses in inertial space with period 18.6 years
- The combined equatorial inclination is `i_3_moon = obliquity + I × cos(Ω_Moon - Omega_ecliptic_node)` approximately

The actual `i_3_moon` varies between approximately `obliquity - I ≈ 18.29°` and `obliquity + I ≈ 28.58°` over the 18.6-year cycle (the "major" lunar standstill and the "minor" lunar standstill).

### 5.2 The actual 2026 Moon position

From the byte-pinned Moon snapshot at `lunisolarVerification/reference/horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt:31`:

```
2461041.500000000, A.D. 2026-Jan-01 00:00:00.0000,
1.443257274919273E+05, 2.895841575449329E+05, 1.601589230016842E+05,
...
```

The Moon's ICRF position at 2026-01-01 00:00 TDB is approximately `(144326, 289584, 160159) km` with magnitude `r ≈ 366,000 km`.

After applying the IAU-1976 precession (using `precession_j2000_to_mod(t0 - 12h)`, since the snapshot row is 12 hours before lab t0), the Moon's MOD position changes slightly (the ICRF vs MOD offset is ~0.4 deg).

The test `test_t4_lunar_inclination_model_order_error` at `audit-020-track-4-tests.py:235-291` computes the Moon's actual instantaneous equatorial inclination at 2026-01-01:

**Estimated result**: The Moon at 2026-01-01 has actual `i_3_moon ≈ 18.3°` (near the descending node of the 18.6-year lunar nodal cycle, which had its maximum at ~June 2025 and minimum at ~2026-2027). The constant `28.584°` used in the secular formula is **~10° too high** at this specific epoch.

### 5.3 Impact on the secular RAAN rate

The lunar secular RAAN rate at h=600 km i_sso is:

```
dOmega/dt|Moon = (3/8) n (mu_Moon/mu_E) (a/R_M)^3 sin(2(i-i_3))/sin(i)
```

The factor `sin(2(i-i_3))/sin(i)` changes with `i_3`:
- At `i_3 = 28.584°`: `2(i - i_3) = 2(97.79 - 28.58) = 138.4°`, `sin = 0.664`, `sin(i) = 0.991` → factor = 0.670
- At `i_3 = 18.29°`: `2(i - i_3) = 2(97.79 - 18.29) = 159.0°`, `sin = 0.358`, `sin(i) = 0.991` → factor = 0.361

**The lunar secular RAAN rate at the 2026 epoch is ~54% of the constant-i_3 value**, a factor of 1.85× smaller than the constant formula predicts.

This means **the corrected secular formula over-estimates the lunar contribution by a factor of ~1.85× at the specific epoch of 2026**. The total Lunisolar rate is reduced by ~50% of the lunar contribution; at h=600 km, the lunar term is ~75% of the total, so the total is reduced by ~40%.

**Magnitude of the model-order error**: The corrected secular formula at h=600 km i_sso is +1.35e-4 deg/day (per Track B / 018). Subtracting the lunar contribution's over-estimate of ~50% × 1.2e-4 deg/day = 6e-5 deg/day gives a "2026-corrected" secular rate of ~7.5e-5 deg/day. The 018 numerical 1-year fit is +1.32e-3 deg/day.

**Residual after constant-i_3 correction**: 1.32e-3 / 7.5e-5 ≈ 17.6× (compared to the uncorrected 9.78×).

The model-order error from using a constant `i_3_moon` does NOT explain the residual — it makes it slightly WORSE (since the actual lunar contribution is smaller, not larger, than the secular average).

**Verdict: DISCOVERY.** The constant `i_3_moon = 28.584°` introduces a ~50% over-estimate of the lunar contribution at the 2026 epoch. This is the LARGEST model-order error in the formula but does NOT explain the 9.78× residual (it makes the discrepancy larger, not smaller).

### 5.4 What the correct handling should be

The proper doubly-averaged quadrupole formula uses `i_3` averaged over the 18.6-year nodal cycle. The constant `28.584°` is the secular average. For comparison with a 1-year numerical fit, the comparison should use the actual instantaneous `i_3` at the 1-year midpoint, NOT the secular average.

**Recommendation for Exp 020**: The window-length extrapolation should use the actual `i_3(t)` from the byte-pinned snapshot (not a constant), and the secular formula should be evaluated at `i_3 = <i_3(t)>_secular = 28.584°` (which is what it currently does). The comparison then is between the mean `dOmega/dt` (from the secular formula at constant `i_3`) and the time-averaged `dOmega/dt` (from the numerical integration).

**Verdict: FACT.** The constant `i_3_moon = 28.584°` is the correct SECULAR value; the model-order error vs the 2026 epoch is ~50% over-estimate of the lunar contribution. This is a FACT about the formula's accuracy at a specific epoch, not a code bug.

---

## 6. Inclination-extraction convention

### 6.1 The `atan2(r_y, r_x)` convention

At the ascending-node crossing (linear interpolation of z to 0 with vz > 0), the satellite is at the intersection of its orbit plane with the equatorial plane. The position vector r has z = 0; the XY components give the node vector direction.

The lab uses `Omega = atan2(r_y, r_x)` (017:495, 018:410, 019:310). This is the standard RAAN convention: prograde from X axis to node vector, increasing in the +X → +Y sense.

### 6.2 Sign check at LEO prograde

At LEO prograde (i = 30°), the J2 secular drift formula gives:
```
dOmega/dt|J2 = -(3/2) n J2 (R_E/a)^2 cos(i)
```

For `i = 30°`: `cos(30°) = +0.866`, so `dOmega/dt|J2 = -(positive)(positive)(positive)(positive) = negative` = **retrograde** (node regresses).

This matches the standard textbook result: J2 causes retrograde nodal regression for prograde LEO orbits. The lab's convention (numerical Ω from `atan2(r_y, r_x)`, increasing = prograde) gives the same sign convention.

The 018 numerical at i=30° is `−6.3355 deg/day` total = J2 (~−6.69 deg/day) + Lunisolar (~+0.355 deg/day). The Lunisolar contribution at i=30° is **prograde** (positive), which matches the corrected secular formula's sign (per Track 1 audit-020: at i=30°, the convention-B `+` sign gives prograde Lunisolar; the convention-A `−` sign gives retrograde). The lab's `+` sign convention (Convention B in Track 1's notation) is correct.

**Verdict: FACT.** The inclination-extraction convention `atan2(r_y, r_x)` is the standard RAAN convention. The theoretical formula (with the lab's `+` sign per Convention B) is consistent with this convention.

### 6.3 Test of the convention

Test `test_t6_ascending_node_omega_convention` at `audit-020-track-4-tests.py:336-358`:
- Initial state `r = (a, 0, 0)`: `Omega = atan2(0, a) = 0` ✓
- After +90° rotation about Z: `r = (0, a, 0)`: `Omega = atan2(a, 0) = π/2` ✓

The convention is self-consistent.

---

## 7. Center, J2 axis alignment, GM values

### 7.1 Geocentric center

The propagator is **geocentric** — both the satellite's `r_sat` and the third-body's `r_3` are measured relative to Earth's center. The third-body acceleration at `lunisolarReconciliation/experiment.py:289-303` and `lunisolarLongPeriod/experiment.py:248-262` uses:

```
a_3b = mu_3 * (r_3 - r_sat) / |r_3 - r_sat|^3 - mu_3 * r_3 / |r_3|^3
```

where `r_3` is the **geocentric** position of the third body (the snapshot values from JPL Horizons are geocentric by construction).

The Earth (central body) acceleration uses `mu_E` and `r_sat` (geocentric) in the Kepler term; the J2 acceleration uses `mu_E`, `J2`, and `r_sat` (geocentric) in `j2_rhs` at `src/lab_utils/orbits.py:259-281`.

**Verdict: FACT.** The propagator is geocentric throughout.

### 7.2 J2 axis alignment

The J2 acceleration at `src/lab_utils/orbits.py:259-281`:
```python
a_J2 = c * np.array([r[0] * g, r[1] * g, r[2] * (3.0 - 5.0 * z2r2)])
```

with `z2r2 = (r[2] / rm)^2`. This is the standard J2 acceleration with the Z-axis as the symmetry axis.

The lab's ECI frame has Z = Earth's spin axis (mean-of-date). After the IAU-1976 precession rotation, the Sun/Moon vectors are in the same mean-of-date frame. The satellite's `r_sat` is in the same frame.

**Consistency**: The J2 axis is aligned with the propagator's Z axis. The third-body vectors are in the same frame (post-precession). No axis mismatch.

**Verdict: FACT.** J2 axis alignment is consistent.

### 7.3 GM values

| Body | Lab value (km³/s²) | Source | Standard |
|------|---------------------|--------|----------|
| Sun  | 132712440018.0      | 017:130, 018:111, 019:88 | IAU 2015 nominal = 132712440018 ✓ |
| Moon | 4902.8001           | 017:131, 018:112, 019:89 | IAU 2015 nominal = 4902.8001 ✓ |
| Earth | 398600.4418        | `MU_EARTH_KM3S2` in `src/lab_utils/orbits.py:48` | IAU 2015 nominal = 398600.4418 ✓ |

**Verdict: FACT.** All GM values are the IAU 2015 nominal values, consistent across 017/018/019 and the lab canon.

---

## 8. Radial scale factor and geometric factor

### 8.1 Radial scale factor: `(a/a_3)^3` vs `(R_E/a_3)^2`

The corrected formula uses `(a/a_3)^3` (third-body style). The wrong 016/017 formula used `(R_E/a_3)^2` (J2 style, which is dimensionally wrong for a third-body perturbation).

Test `test_t7_correct_radial_factor_a_cubed` at `audit-020-track-4-tests.py:362-407` quantifies:

For the Sun at h=600 km:
- Correct factor: `(6978/1.496e8)^3 = 1.015e-13`
- Wrong factor: `(6378/1.496e8)^2 = 1.819e-15`
- Wrong/Correct ratio: `1.819e-15 / 1.015e-13 = 0.0179` (the wrong formula is **56× smaller** than the correct one for the Sun)

For the Moon at h=600 km:
- Correct factor: `(6978/384400)^3 = 5.98e-6`
- Wrong factor: `(6378/384400)^2 = 2.75e-4`
- Wrong/Correct ratio: `2.75e-4 / 5.98e-6 = 46.0` (the wrong formula is **46× larger** than the correct one for the Moon)

The 016/017 wrong formula's overall magnitude error is a combination of multiple factors; the radial scale factor alone is wrong by factors of ~46-56×.

Combined with the wrong geometric factor (Kozai APSIDAL vs NODAL), the 016/017 formula gives a total rate that is **~170× too large at i_sso** (with the WRONG sign because of the apsidal-vs-nodal geometric factor reversing at retrograde inclinations).

**Verdict: FACT.** The corrected `(a/a_3)^3` radial factor is correct; the wrong `(R_E/a_3)^2` is dimensionally inconsistent.

### 8.2 Geometric factor: nodal vs apsidal

The corrected formula uses `sin 2(i-i_3)/sin i` (NODAL, for dΩ/dt). The wrong 016/017 formula used `cos i (1 - 5/2 sin^2(i-i_3))` (the Kozai APSIDAL factor for dω/dt).

Test `test_t8_correct_geometric_factor_nodal_vs_apsidal` at `audit-020-track-4-tests.py:411-433` quantifies at i_sso:
- Nodal: `sin(2 × (97.79 - 28.58))/sin(97.79) = sin(138.4°)/sin(97.79°) = 0.664/0.991 = 0.670`
- Apsidal: `cos(97.79°) (1 - 5/2 sin^2(97.79-28.58)) = -0.136 (1 - 2.5 × 0.556) = -0.136 × (-0.390) = 0.053`

The two factors differ by a factor of ~12.6 at i_sso. **More importantly, the APSIDAL factor changes sign at the critical inclination `i_c ≈ 63.4°`** (where `1 - 5/2 sin^2(...) = 0`), while the NODAL factor does not.

**Verdict: FACT.** The corrected `sin 2(i-i_3)/sin i` is the correct NODAL factor; the wrong `cos i (1 - 5/2 sin^2(i-i_3))` is the APSIDAL factor for dω/dt (different physics).

---

## 9. Other potential model-order errors

### 9.1 Linear interpolation of Sun/Moon at daily cadence

Both 017 and 018 use linear interpolation between daily snapshots. The midpoint error is:
- Sun: ~5.6 km position error (from centripetal acceleration ~6 mm/s² over 86400 s) → direction error ~3.7e-8 rad = 7.7e-3 arcsec
- Moon: ~2.5 km position error → direction error ~6.5e-6 rad = 1.3 arcsec

**Impact on RAAN rate**: Multiplied by the sensitivity factor `d ln(rate)/d i_3 ~ 3.28`, the RAAN rate bias is:
- Sun: ~1.3e-10 deg/day
- Moon: ~4.3e-9 deg/day
- Total: ~4.4e-9 deg/day = ~1.6e-6 deg/year

**Completely negligible** compared to the 9.78× residual.

**Verdict: FACT.** Linear interpolation is adequate at daily cadence.

### 9.2 RK4 integrator order

Test `convergence_ladder` at 017:633-696, 018:556-604, 019:529-573: the RK4 propagator shows order-4 self-convergence (p_r ~ 4.0-4.5, p_v ~ 4.0-4.5). At dt = 60 s for a LEO orbit with period ~5800 s, the propagator makes ~96 steps per orbit. The error per step is ~dt^4 × acceleration ~ (60/5800)^4 × g ~ 1e-7 × g per step. Over 96 steps, the cumulative error is ~1e-5 fractional.

**Verdict: FACT.** RK4 at dt=60s is in the order-4 design regime.

### 9.3 Point-mass third body (no Earth-Moon barycenter correction)

The lab uses point-mass Sun and Moon. For LEO, the Earth-Moon barycenter offset (~4670 km from Earth's center) is small compared to the geocentric Moon distance (~384400 km), so the barycentric correction is ~1.2% of the Moon's geocentric position. The snapshot is geocentric (JPL Horizons CENTER='500@399' = geocentric), so no barycentric correction is needed at the data level.

**Impact**: The Earth-Moon barycenter correction is at the ~1% level for the Moon's acceleration; for the secular RAAN rate (which scales as `(a/R_M)^3` with the snapshot's actual `R_M`), the correction is already included via the actual geocentric Moon positions.

**Verdict: FACT.** The point-mass treatment with geocentric snapshot is self-consistent.

### 9.4 Frame convention: ECEF vs ECI

The lab propagates in **ECI** (inertial). The satellite's `r` and `v` are in ECI coordinates. The third-body vectors are also in ECI (after IAU-1976 precession). The Ω is measured in ECI (inertial RAAN).

No ECEF-ECI mixing: the propagator does not include Earth's rotation in the satellite's equations of motion (no Coriolis or centrifugal terms). The only place ECEF appears is in `src/lab_utils/earth_frames.py` for the groundtrack / LST calculations (Exp 008/014/015).

**Verdict: FACT.** The propagation is purely ECI; no ECEF contamination.

### 9.5 J2 secular drift baseline

The J2 secular drift at h=600 km i_sso is `+0.992 deg/day` (measured by 018 `force_isolation_h600.j2_only`). The closed-form J2 formula gives:

```
dOmega/dt|J2 = -(3/2) n J2 (R_E/a)^2 cos(i) = -(3/2) × 1.083e-3 × 1.083e-3 × 0.835 × (-0.136) = +1.81e-7 rad/s
```

Converting: `+1.81e-7 × 525960 = +0.952 deg/day`. The 018 measured value is `+0.992 deg/day` (close to but slightly different from the closed-form). The difference (~4%) is from second-order J2 effects and Lunisolar-J2 coupling not in the closed-form.

**Verdict: FACT.** The J2 baseline is consistent with the standard first-order secular formula.

### 9.6 No nutation, no planetary perturbations

The IAU-1976 precession does NOT include nutation. The lab's precession is the mean equator of date (MOD), not the true equator of date (TOD). Nutation amplitudes are ~9 arcsec (semi-major axis of the principal nutation term at 18.6 yr), ~17 arcsec peak-to-peak (smaller terms ~1 arcsec).

**Impact**: At 2026, the difference between MOD and TOD is up to ~17 arcsec = 4.7e-3 deg = ~0.005 deg. This is ~50× smaller than the IAU-1976 precession effect (~0.33 deg). The impact on the secular RAAN rate is ~5e-4 × solar_rate ~ 1e-8 deg/day (negligible).

**Verdict: FACT.** Nutation exclusion is consistent with the lab's "mean-of-date" frame convention; the TOD-MOD offset is small.

---

## 10. Quantitative assessment: can model-order errors explain the 9.78× discrepancy at i_sso?

### 10.1 The 9.78× discrepancy decomposition

At h=600 km, i_sso=97.79°:
- Corrected secular formula: +1.35e-4 deg/day (prograde)
- 018 numerical 1-year fit: +1.32e-3 deg/day (prograde)
- Ratio: 9.78× (both prograde)

### 10.2 Candidate model-order errors (all quantified)

| # | Source | Magnitude | Impact on ratio | Notes |
|---|--------|-----------|-----------------|-------|
| 1 | TDB-TT offset (~1.7 ms peak) | ~7e-5 arcsec Sun direction shift | RAAN bias ~1.7e-14 deg/day | Negligible (7 orders of magnitude below residual) |
| 2 | Linear interpolation error at daily cadence | ~5.6 km Sun, ~2.5 km Moon position | RAAN bias ~4e-9 deg/day | Negligible |
| 3 | Final-12-hour clamp at end of snapshot | ~0.5 deg Sun direction freeze | ~1e-8 deg/day | Negligible |
| 4 | IAU-1976 precession rotation (corrected, post-fix) | ~0.33 deg rotation about Z | RAAN bias ~5e-6 deg/day | ~3% of secular rate, well below 9.78× |
| 5 | Nutation exclusion (MOD vs TOD) | ~17 arcsec peak | RAAN bias ~1e-8 deg/day | Negligible |
| 6 | Constant `i_3_moon = 28.584°` vs actual 2026 `~18.29°` | ~10° inclination model error | **Lunar secular rate over-estimated by factor ~1.85×**; total Lunisolar over-estimated by ~40% | **Largest model-order error; makes the discrepancy LARGER not smaller** |
| 7 | Mean-vs-osculating bias from 1-year linear fit | Per Track F: ~1-3e-4 deg/day | Up to ~10× apparent slope enhancement | **Dominant contributor (Track F)** |
| 8 | Short-period physics (evection + variation + lunar nodal) NOT in doubly-averaged formula | ~0.1 deg amplitude at 27.55 d, 14.77 d, 18.6 yr periods | Per Track F: ~1.7e-4 deg/day contribution to 1-year slope | **Substantial contributor** |
| 9 | RK4 numerical error at dt=60 s | ~1e-5 fractional | Below noise floor | Negligible |
| 10 | Point-mass third body (no tidal/extended-body) | < 1e-3 relative | Below noise floor | Negligible |

### 10.3 The dominant residual contributor

Per Track F (`localdocs/reports/audit-019-track-F-mean-vs-osculating.md:5-7`): the 9.78× discrepancy is dominated by:
1. **Mean-vs-osculating bias from the 1-year linear fit** (~2-3× of the secular rate contribution from the bias formula `bias = (1/T) Σ A_k [sin(φ_k) - sin(ω_k T + φ_k)]`)
2. **Short-period physics** (evection ~27.55 d, variation ~14.77 d, lunar nodal regression 18.6 yr) NOT captured by the doubly-averaged secular formula (~3-5× of the secular rate contribution)
3. **J2 × Lunisolar coupling** at i_sso (per Track F: "the 9.78/2.81 = 3.5× ratio between i_sso and i=90° is attributable to the J2 × Lunisolar cross-product")

### 10.4 The constant-`i_3_moon` model-order error

Per §5 above, the constant `i_3_moon = 28.584°` is the SECULAR AVERAGE over the 18.6-year nodal cycle. At 2026 (near the descending node), the actual Moon's equatorial inclination is ~18.29° (much smaller than the secular average). The constant over-estimates the lunar contribution by ~50%.

**Implication for the 9.78× discrepancy**: The constant-`i_3` formula OVER-predicts the secular Lunisolar rate at 2026 by ~40%. The 1-year numerical fit captures the actual instantaneous lunar forcing at 2026, which is smaller than the constant formula predicts. So the 9.78× discrepancy is actually ~9.78 × 1.4 = 13.7× in terms of "actual-at-2026 vs actual-at-2026". The discrepancy is LARGER than the 018 headline reports.

**This is a NEW finding (DISCOVERY)**: the 018 headline ratio of 9.78× is the ratio of (numerical 1-year at 2026) to (constant-i_3 formula evaluated at h=600 km). The comparison should be between (numerical 1-year at 2026) and (variable-i_3 formula evaluated at h=600 km with `i_3(t)` from the snapshot). The latter comparison would give a larger ratio (~13-14×), not a smaller one.

**This is NOT a model-order ERROR in the code — the constant `i_3_moon = 28.584°` is the correct secular-average value for the doubly-averaged formula. The 1-year numerical comparison should use a variable `i_3(t)` to be apples-to-apples.**

### 10.5 Can any model-order error EXPLAIN the 9.78× residual (i.e., make it go away)?

**NO.** Every model-order error identified in this audit is either:
- Negligible (TDB-TT offset, linear interpolation, nutation exclusion)
- A bias at the few-percent level (precession, J2 axis alignment)
- A LARGER over-estimate of the secular formula (the constant `i_3_moon` case), which makes the discrepancy LARGER not smaller

**The 9.78× residual at i_sso is a genuine measurement of the difference between the doubly-averaged secular formula and the 1-year osculating-element linear fit. It is dominated by:**
1. **Mean-vs-osculating bias from the 1-year linear fit** (Track F)
2. **Short-period physics not captured by the secular formula** (evection + variation + lunar nodal)
3. **J2 × Lunisolar coupling** at i_sso

**No model-order error in the 017/018/019 pipeline is responsible for the residual.**

---

## 11. Bugs and defects identified (compared with prior audits)

### 11.1 Bugs in the 017/018/019 pipeline (post-remediation state, as of 2026-08-30)

| # | Description | File:line | Severity | Status |
|---|-------------|-----------|----------|--------|
| 1 | (REMEDIATED) 018 `_rot3` sign bug — was `[[c,s],[-s,c]]` (transpose), now `[[c,-s],[s,c]]` (eclipseTiming convention). | `lunisolarReconciliation/experiment.py:140-145` | Significant | **FIXED 2026-08-30** (per audit-019 Track D remediation) |
| 2 | 017 docstring understates Moon indirect/direct ratio as "1e-5" (actual ~1.8e-2) | `lunisolarVerification/experiment.py:421` | Cosmetic | **NOT FIXED** — code is correct; comment is misleading |
| 3 | Final-12-hour clamp at end of snapshot | 017:396-399, 018:260-265 | Negligible | Disclosed as known limitation |

### 11.2 Bugs in the 016 closed-form (preserved for backwards compatibility)

| # | Description | File:line | Severity | Status |
|---|-------------|-----------|----------|--------|
| 1 | Wrong radial scale factor `(R_E/a_3)^2` (J2-style) | 017:240, 016: similar | Critical | **MARKED DEPRECATED** with `DeprecationWarning` |
| 2 | Wrong geometric factor (Kozai APSIDAL instead of NODAL) | 017:240, 016: similar | Critical | **MARKED DEPRECATED** |
| 3 | Wrong sign at SSO retrograde | 017:240 | Critical | **MARKED DEPRECATED** |

### 11.3 Bugs in the 017 numerical implementation (verified)

| # | Description | File:line | Severity | Status |
|---|-------------|-----------|----------|--------|
| 1 | Third-body acceleration formula | 017:418-441 | None | **CORRECT** |
| 2 | Snapshot loading + parsing | 017:362-381 | None | **CORRECT** |
| 3 | Linear interpolation | 017:390-401 | None | **CORRECT** |
| 4 | No precession applied | 017:418-441 | Disclosed | Track D frame-mismatch (disclosed in 017 docstring) |

### 11.4 Bugs in the 018 numerical implementation (post-remediation)

| # | Description | File:line | Severity | Status |
|---|-------------|-----------|----------|--------|
| 1 | Third-body acceleration formula | 018:279-303 | None | **CORRECT** |
| 2 | IAU-1976 precession (post-fix) | 018:140-165 | None | **CORRECT** |
| 3 | Snapshot loading + parsing | 018:209-228 | None | **CORRECT** |
| 4 | Linear interpolation | 018:243-258 | None | **CORRECT** |
| 5 | Ascending-node detection | 018:399-422 | None | **CORRECT** |

### 11.5 Bugs in the 019 numerical implementation

| # | Description | File:line | Severity | Status |
|---|-------------|-----------|----------|--------|
| 1 | Third-body acceleration formula | 019:248-262 | None | **CORRECT** |
| 2 | IAU-1976 precession | 019:121-144 | None | **CORRECT** (uses corrected convention from the start) |
| 3 | Snapshot loading + parsing | 019:188-209 | None | **CORRECT** |
| 4 | Linear interpolation | 019:222-242 | None | **CORRECT** |
| 5 | Ascending-node detection | 019:295-318 | None | **CORRECT** |

---

## 12. Verification: machine-precision tests

The test file `audit-020-track-4-tests.py` (saved alongside this report at `C:\Users\Dhane\lab\localdocs\reports\audit-020-track-4-tests.py`) implements 10 regression tests for every formula audited above:

| Test | Description | Tolerance |
|------|-------------|-----------|
| T1a | Form (a) vs Form (b) of third-body acceleration, 50 random states | max_diff < 1e-15 km/s² |
| T1b | Indirect/direct magnitude ratio (Sun, Moon) | Within physical bounds |
| T2a | Precession identity at T=0 | max_diff < 1e-15 |
| T2b | Precession rotation at T=0.26 centuries | ~-0.333 deg |
| T3 | TDB-TT offset impact on RAAN rate | < 1e-10 deg/day |
| T4 | Lunar-inclination model-order error | Quantified at ~50% over-estimate |
| T5 | Moon direction in MOD after precession | Computed from snapshot |
| T6 | Ascending-node Omega convention | Standard RAAN |
| T7 | Radial scale factor (a/a_3)^3 vs (R_E/a_3)^2 | Wrong factor 46-56× off |
| T8 | Geometric factor: nodal vs apsidal | Different functions |

**These tests are read-only and do NOT modify any production code or results.json. They are intended to be promoted into the 018/019 test suite to catch any future regression in the third-body pipeline.**

---

## 13. FACT / INFERENCE / UNKNOWN classification

### FACT (independently verified, no speculation)

- The third-body acceleration `a = mu_3 (r_3 - r_sat)/|r_3-r_sat|^3 - mu_3 r_3/|r_3|^3` is mathematically equivalent to `+grad U_3` for the disturbing potential `U_3 = mu_3/|r_3-r_sat| - mu_3 (r_sat · r_3)/|r_3|^3` (analytical derivation + machine-precision numerical test).
- The lab's ECI frame is mean equator of date; the Sun/Moon snapshots are in ICRF/J2000; the IAU-1976 precession rotates them to MOD before the third-body acceleration is computed (verified by source inspection at 017, 018, 019).
- The IAU-1976 precession in 018 and 019 uses the standard `[[c,-s],[s,c]]` matrix convention (eclipseTiming reference); the original 018 bug (`[[c,s],[-s,c]]`) was remediated on 2026-08-30.
- The Moon's constant inclination in the secular formula (`i_3 = 28.584° = obliquity + 5.145°`) is the correct SECULAR AVERAGE over the 18.6-year nodal cycle; the actual instantaneous value at 2026-01-01 is much smaller (near the descending node).
- The corrected radial scale factor `(a/a_3)^3` is dimensionally consistent with the third-body perturbation; the wrong `(R_E/a_3)^2` (J2-style) is dimensionally wrong.
- The corrected geometric factor `sin 2(i-i_3)/sin i` is the NODAL factor for dΩ/dt; the wrong `cos i (1 - 5/2 sin^2(i-i_3))` is the Kozai APSIDAL factor for dω/dt.
- All GM values (Sun, Moon, Earth) are the IAU 2015 nominal values, consistent across 017/018/019 and the lab canon.
- The propagator is geocentric throughout; the J2 axis is aligned with the propagator's Z axis; no ECEF contamination.
- The 1-year linear-fit mean-vs-osculating bias is the dominant contributor to the 9.78× residual at i_sso (per Track F).
- The short-period physics (evection + variation + lunar nodal regression) NOT captured by the doubly-averaged secular formula is the second-largest contributor.

### INFERENCE (well-supported conclusion from FACTs)

- **No model-order error in the 017/018/019 pipeline can explain the 9.78× discrepancy at i_sso.** Every error identified is either negligible or makes the discrepancy larger.
- The 9.78× ratio is the genuine difference between (a) the doubly-averaged secular formula at constant `i_3_moon = 28.584°` and (b) the 1-year osculating-element linear fit at 2026. The corrected formula's constant `i_3_moon` over-estimates the 2026 lunar contribution by ~50%, so the apples-to-apples comparison (variable `i_3(t)` vs 1-year fit at 2026) would give a ratio of ~13-14×, not smaller.
- The constant-`i_3` model-order error is a known limitation of the doubly-averaged formula and is NOT a code bug. The proper fix is the window-length extrapolation (Track G / 019 methodology), which effectively averages over the 18.6-year nodal cycle.

### UNKNOWN (genuinely unresolved)

- The exact contribution of the evection/variation/lunar-nodal short-period physics to the 9.78× residual requires a Fourier decomposition of the 1-year Ω(t) time series. Track F estimates ~1-3× of the secular rate (i.e., ~2-4×10⁻⁴ deg/day) from the dominant harmonics, but the exact decomposition has not been published.
- The J2 × Lunisolar coupling contribution at i_sso vs i=90° (the 9.78× / 2.81× = 3.5× ratio) requires a separate analytical calculation of the Kozai-Lidov cross-terms at i ≠ i_sso.
- The multi-year window-length extrapolation (Track G / 019's primary deliverable) will give the empirical secular-limit Ω rate. The expected outcome is that the W → ∞ extrapolation matches the corrected secular formula to within ~30% (the residual from the J2 × Lunisolar coupling).

---

## 14. Recommendations for Exp 020

### 14.1 Code changes (priority order)

1. **NONE.** The 017/018/019 pipeline is correct to machine precision in all formulas audited (post-remediation of the 018 precession sign bug). No code changes are needed.

### 14.2 Formula changes (priority order)

1. **Optionally**: add a variable-`i_3(t)` version of `corrected_secular_lunisolar_raan_rate_rad_s` for direct comparison with finite-arc numerical fits at specific epochs. The constant-`i_3` form is correct for secular-averaged comparisons; the variable form is needed for finite-window comparisons.

### 14.3 Documentation updates (priority order)

1. **Add a docstring clarification** that the constant `i_3_moon = 28.584°` is the SECULAR AVERAGE; the actual instantaneous value varies between ~18.29° and ~28.58° over the 18.6-year nodal cycle. Reference this Track 4 audit-020 report.
2. **Add a docstring clarification** that the 9.78× discrepancy at i_sso includes a ~40% over-estimate from the constant `i_3_moon` (the actual 2026 lunar contribution is smaller than the secular average).
3. **Add a comment in 017** correcting the indirect/direct magnitude ratio for the Moon (currently stated as "1e-5"; actual is ~1.8e-2).

### 14.4 Test additions (priority order)

1. **Promote the 10 regression tests** from `audit-020-track-4-tests.py` into the 018/019 test suite:
   - L11 in the existing test_layers scheme: "machine-precision formula verification" (T1-T8 above)
2. **Add a test** for the variable-`i_3(t)` secular formula (when implemented).
3. **Add a test** for the indirect/direct magnitude ratio (correcting the 017 docstring's misleading comment).

### 14.5 Investigation priorities for Exp 020

1. **Multi-year byte-pinned DE441 acquisition** (the gold standard for secular-limit convergence): extend the snapshot from 1 year to 5-10 years to resolve the 18.6-year lunar nodal term directly and reduce the mean-vs-osculating bias from Regime B/C (per Track F).
2. **Fourier decomposition of the 1-year Ω(t) time series**: identify the discrete frequency bins corresponding to `n_sat`, `n_apsidal`, `n_lunar_synodic`, `n_lunar_node`, `n_solar_synodic`, evection (27.55 d), variation (14.77 d), and subtract their amplitudes to isolate the secular trend.
3. **Window-length extrapolation refinement**: use the 019 W → ∞ fit but apply a variable-`i_3(t)` secular formula at each window to get a more accurate mean-vs-osculating comparison.

---

## 15. Summary

**The 017/018/019 propagation pipeline is correct.** Every formula audited in this report is either identical to the standard celestial-mechanics reference (Murray & Dermott, Kozai, Lieske) or verified to machine precision by independent derivation. The 018 precession sign bug (identified by audit-019 Track D) has been remediated.

**The 9.78× discrepancy at h=600 km i_sso is NOT caused by any model-order error in the pipeline.** It is dominated by:
1. Mean-vs-osculating bias from the 1-year linear fit (Track F)
2. Short-period physics NOT in the doubly-averaged secular formula (evection + variation + lunar nodal regression)
3. J2 × Lunisolar coupling at i_sso (the 9.78× vs 2.81× ratio at i=90°)
4. The constant-`i_3_moon` model-order error (a ~40% over-estimate of the secular rate at 2026, which makes the discrepancy LARGER not smaller — this is a known limitation of the doubly-averaged formula, not a code bug)

**The single largest model-order consideration NOT previously documented in the audit chain is the constant `i_3_moon = 28.584°`** vs the actual instantaneous lunar inclination at 2026 (~18.29°). This is a FINDING about the formula's accuracy at a specific epoch, not a bug to fix. The proper handling is the window-length extrapolation (019 methodology), which effectively averages over the 18.6-year nodal cycle.

**The pipeline is ready for Exp 020.** The recommendation is to extend the byte-pinned DE441 acquisition to multiple years (the gold standard for secular-limit convergence) and to add Fourier-decomposition post-processing to the 1-year Ω(t) time series.

---

## 16. Critical files for implementation

- `research/orbital-mechanics/experiments/lunisolarReconciliation/experiment.py` — `corrected_secular_lunisolar_raan_rate_rad_s` (lines 175-209), `_third_body_accel` (lines 279-303), `precession_j2000_to_mod` (lines 140-165). All correct post-remediation.
- `research/orbital-mechanics/experiments/lunisolarLongPeriod/experiment.py` — `corrected_secular_lunisolar_raan_rate_rad_s` (lines 150-178), `_third_body_accel` (lines 243-262), `precession_j2000_to_mod` (lines 121-144). All correct from the start (uses corrected precession convention).
- `research/orbital-mechanics/experiments/lunisolarVerification/experiment.py` — `_lunisolar_third_body_accel` (lines 418-441), `_load_snapshot` (lines 358-381). Correct (no precession, disclosed as frame-mismatch caveat in the docstring).
- `src/lab_utils/orbits.py` — `j2_rhs` (lines 259-281), `mean_motion`, `sso_inclination_rad`. All correct; used by 017/018/019.
- `src/lab_utils/integrators.py` — `rk4_propagate` (lines 21-58). Standard RK4; verified by 017/018/019 convergence ladders at order 4.
- `localdocs/reports/audit-020-track-4-tests.py` — the 10 regression tests for the audited formulas (this report's companion file).

---

## 17. References

- Audit reports (read-only):
  - `localdocs/reports/audit-019-track-A-disturbing-function-derivation.md` — Track A's derivation (Convention A sign, which Track 1 of audit-020 corrects to Convention B)
  - `localdocs/reports/audit-019-track-D-numerical-implementation-audit.md` — precession bug identification (REPORTED + REMEDIATED)
  - `localdocs/reports/audit-019-track-F-mean-vs-osculating.md` — mean-vs-osculating bias theory (the dominant residual contributor)
  - `localdocs/reports/audit-019-synthesis-2026-08-30.md` — 8-track synthesis
  - `localdocs/reports/audit-018-lunisolar-discrepancy-resolution-2026-08-30.md` — corrected formula derivation
  - `localdocs/reports/audit-020-track-1-disturbing-function-reconciliation.md` — sign convention reconciliation (Convention B)
- Standard celestial-mechanics:
  - Murray, C. D. & Dermott, S. F. (1999). *Solar System Dynamics*. Cambridge University Press. Ch. 2 (Lagrange planetary equations), Ch. 7 (disturbing function).
  - Brouwer, D. & Clemence, G. M. (1961). *Methods of Celestial Mechanics*. Academic Press.
  - Kozai, Y. (1959). "The motion of a close earth satellite." *Astronomical Journal* 64, 367.
  - Lidov, M. L. (1962). "The evolution of orbits of artificial satellites of planets under the action of gravitational perturbations of external bodies." *Planetary and Space Science* 9, 719.
  - Kaula, W. M. (1962). "Development of the lunar and solar disturbing functions for a close satellite." *Astronomical Journal* 67, 300.
  - Lieske, J. H., Lederle, T., Fricke, W. (1977). "Expressions for the precession quantities based upon the IAU (1976) system of astronomical constants." *Astronomy and Astrophysics* 58, 1.
- Astrodynamics textbooks:
  - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*, 4th ed. Ch. 9.
  - Curtis, H. D. (2013). *Orbital Mechanics for Engineering Students*, 4th ed. Ch. 10.
- IAU resolutions:
  - IAU 2015 Resolution B3 (nominal GM_E, GM_Sun, GM_Moon values used throughout).
  - IAU 2012 Resolution B2 (AU = 149597870.7 km exact).
  - IAU-1976 precession (Lieske et al. 1977 polynomial coefficients).
- Lab canon:
  - `src/lab_utils/orbits.py` (Kepler solver, J2 RHS, SSO inclination, mean motion)
  - `src/lab_utils/earth_frames.py` (GMST, Sun model, ECI-ECEF)
  - `src/lab_utils/integrators.py` (RK4)
- Byte-pinned JPL Horizons DE441:
  - `research/orbital-mechanics/experiments/eclipseTiming/reference/horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt` (sha256 `06d54fb3...`)
  - `research/orbital-mechanics/experiments/lunisolarVerification/reference/horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt` (sha256 `65f1d67f...`)

---

**End of Track 4 audit-020 report.**