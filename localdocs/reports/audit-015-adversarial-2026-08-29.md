# Adversarial Audit Report — Exp 015 "4 min/day LST drift" Claim

> **Scope:** Hostile adversarial review of the LST-drift claim in
> `research/orbital-mechanics/experiments/dawnDuskSSO/` (Exp 015, 2026-08-29).
>
> **Claim under test:**
> "The drift rate is `dOmega/dt − d(Subsolar)/dt = 360.9856 − 360.0 =
> 0.9856 deg/day = 4 min/day`. The LST at the ascending node of a dawn-dusk
> SSO drifts through 24 h over the year."
> (README.md L79-88; replicated in `results.json` findings 1+5 and in
> `localdocs/knowledge/dawn-dusk-sso.md` L35-46.)
>
> **Audit method:** Try every plausible interpretation of the formula; decide
> which (if any) gives 4 min/day; check whether the surviving interpretation
> is the one Exp 015 actually uses; run mutants to break the claim; verify
> against lab_utils directly.

## 0. Executive verdict

| Item | Value |
|---|---|
| **ADVERSARIAL-VERDICT** | **RED** |
| **015-CLAIM-SURVIVES-MUTANTS** | **no** |
| **EXACT-SOURCE-OF-ERROR** | The 015 formula substitutes `dOmega/dt = 360.9856 deg/day` (the **Earth's sidereal rotation rate** = `dGMST/dt`) and `d(Subsolar)/dt = 360.0 deg/day` (NOT the Sun's apparent inertial RA rate, which is 0.9856 deg/day) into a frame-mixed subtraction. The difference `0.9856 deg/day` is the real *sidereal-vs-solar-day* differential (1 sidereal day − 1 mean solar day = 1/365.2422 d ≈ 3.93 min), which is a property of **Earth rotation about its own axis vs. Earth orbit around the Sun** — it has nothing to do with the satellite's orbit. The LST at a fixed ground point on Earth drifts at this rate; the LST at the **orbit-plane node of a true SSO is constant** (0 min/day by SSO design). The 015 code measures the fixed-ground-point LST because its C4 insertion convention `Omega = GMST + site_lon` forces the geodetic node longitude `node_lon = Omega − GMST = site_lon` to be a CONSTANT (not the moving orbit-plane node), so the 24 h/year wraparound is observed at the launch-site geodetic point, not at the satellite's orbit-plane node. The README then ascribes this Earth-rotation artifact to "dawn-dusk SSO physics" and falsely attributes a "+2.2 deg LST drift per year at SSO 600 km" finding to Exp 012, which contains no such quantity. |
| **CORRECT-MAGNITUDE** | **0 min/day drift at the orbit-plane node of a true dawn-dusk SSO**, with a bounded ±0.5 min/day daily swing from the equation-of-time (EoT) and a cumulative EoT envelope of **±16-20 min/year** that returns to zero after one mean-solar year (i.e. it is periodic, NOT secular). Verified directly with `lab_utils.earth_frames.lst_at_node_hours` using a perfect SSO (Ω_dot = 0.985647332 deg/day, geodetic node longitude free) over 2026. The prior audit report `localdocs/reports/audit-015-follow-up-candidates-2026-08-29.md` reaches the same numerical conclusion (max daily drift 0.498 min/day, cumulative envelope ±20 min/year). |

---

## 1. The claim, restated precisely

`dawnDuskSSO/README.md` L80-88:

> A *crucial* physics finding emerged during the analysis (and was caught by
> the hostile adversarial review): the LST at the ascending node of a
> dawn-dusk SSO drifts through 24 h over the year. The drift rate is
> `dOmega/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day = 4 min/day`.

The same paragraph claims the drift comes from "the LST at the ascending
node of a dawn-dusk SSO" and asserts a "multi-year dawn-dusk mission needs
station-keeping" as a direct consequence.

`dawnDuskSSO/experiment.py` L881-883 and `results.json` finding #1 add:

> "consistent with the Exp 012 finding of +2.2 deg LST drift per year at
> SSO 600 km."

`localdocs/knowledge/dawn-dusk-sso.md` L177-180 repeats the attribution.

A targeted grep across the entire lab found **no occurrence of "2.2 deg
LST drift per year at SSO 600 km" anywhere except in Exp 015's own
artifacts**. `orbitClasses/experiment.py`, `orbitClasses/README.md`,
`orbitClasses/results.json`, and `localdocs/knowledge/orbit-classes.md`
contain no such figure. **The Exp 012 attribution is fabricated.**

---

## 2. The four interpretations

A priori, `dOmega/dt` and `d(Subsolar)/dt` can each be evaluated in three
natural frames (inertial ECI, geodetic ECEF, or the unwrapped inertial RA
of the Sun). All four natural interpretations and their resulting drifts:

| Interpretation | `dOmega/dt` (deg/day) | `d(Subsolar)/dt` (deg/day) | Resulting drift | Frame | Verdict |
|---|---|---|---|---|---|
| **A** Both inertial: dΩ/dt = SSO nodal rate; d(sun RA inertial)/dt = sun's mean motion | +0.985647332 | +0.985647332 | **0** | inertial LST formula `LST = 12 + (Ω − α_sun) / 15` | The 015 README's formula text suggests this frame, but then puts 360.9856 and 360.0 in, which contradict the labels. **Inconsistent within the README itself.** |
| **B** dΩ/dt in inertial, d(subsolar lon_ecef)/dt in ECEF | +0.985647332 | −360.0 | **+361.97 deg/day** (or +7.88 min/day) | geodetic LST formula `LST = 12 + (node_lon − subsolar_lon_ecef) / 15` with `node_lon = Ω` (NOT Ω − GMST) | Doesn't match 0.9856. |
| **C** dΩ/dt mislabeled as Earth's sidereal rotation (360.9856), d(Subsolar)/dt mislabeled as 360.0 | +360.9856 | +360.0 | **+0.9856 deg/day = 4 min/day** ✓ | mixed / non-physical | This is the formula the README literally writes. **The result matches the claimed 4 min/day only because the formula uses the wrong labels.** |
| **D** d(node_lon_ecef)/dt for a non-SSO orbit, in ECEF | d(Ω−GMST)/dt = 0 − 360.9856 = −360.9856 for inertial-Ω orbit; d(subsolar_lon_ecef)/dt = α_sun_dot − GMST_dot = 0.9856 − 360.9856 = −360 | −360 | **−0.9856 deg/day = 4 min/day** ✓ | geodetic LST formula applied to a *non-SSO* orbit's node | This is the *real* physics of the 4 min/day figure — but it is the LST drift of the **geodetic node longitude of an inertial orbit** (e.g. a Molniya or non-SSO LEO), NOT a dawn-dusk SSO. |

**Interpretation C is what gives 4 min/day, and it corresponds to a
NON-SUN-SYNCHRONOUS orbit's geodetic node longitude drift, not a true SSO.**

For a **true dawn-dusk SSO** (where `dΩ/dt = dα_sun/dt = +0.985647332
deg/day` by design), the LST at the orbit-plane node — using *any* of the
LST formulas consistently — is **constant** (modulo the EoT envelope):

```
LST = 12 + (Ω − α_sun) / 15      →  dLST/dt = (dΩ/dt − dα_sun/dt)/15 = 0
LST = 12 + (node_lon_ecef − subsolar_lon_ecef)/15  with node_lon_ecef = Ω − GMST
                                     →  dLST/dt = (dΩ − dGMST − dα_sun + dGMST)/15
                                     →  dLST/dt = (dΩ − dα_sun)/15 = 0
```

The 015 README's own formula, with correct labels plugged in (dΩ/dt =
0.9856, dα_sun/dt = 0.9856), gives **0 min/day**, not 4 min/day.

---

## 3. What the 015 code actually computes

`dawnDuskSSO/experiment.py` L262-272:

```python
def lst_at_node_at_t(t_launch_s):
    gmst = gmst_rad_iau1982(t_launch_s)
    raan = gmst + REF_SITE_LON_DEG * DEG
    node_lon = raan - gmst                    # = REF_SITE_LON_DEG * DEG (CONSTANT)
    return lst_at_node_hours(t_launch_s, node_lon)
```

The C4 insertion convention `Ω(t_L) = GMST(t_L) + lon_ref` forces the
**geodetic node longitude at insertion** to be the launch-site longitude,
a CONSTANT. The 015 code therefore measures `LST(t, site_lon)` — the LST
at a fixed geodetic point on Earth, which is a property of Earth rotation
and Earth orbit, NOT of the satellite's orbit.

A numerical sweep with `lab_utils.earth_frames.lst_at_node_hours` at fixed
geodetic longitude = -80.6039° over 2026:

```
LST at fixed -80.6° lon: range = [0.018, 23.981] h
Wraparounds: 365 (one per day)
Continuous drift: 24.0000 h/year = 4 min/day (average) = 0.9856 deg/day
```

The 015 README's headline number "24 h/year" and "4 min/day" both come
from this Earth-rotation property at a fixed ground point.

**The 4 min/day drift is a property of the LAUNCH SITE, not the ORBIT.
It applies to ANY fixed geodetic longitude on Earth** (Greenwich, Paris,
the equator, anywhere). Verified directly:

```
LST at Greenwich (lon = 0) over 2026: range = [0.032, 23.991] h  (also 24 h/year)
LST at -80.6° lon over 2026:        range = [0.018, 23.981] h  (also 24 h/year)
```

A real **dawn-dusk SSO at the orbit-plane node** has the LST range bounded
by the equation-of-time, NOT 24 h/year. Verified directly with `lab_utils`
using a perfect SSO at h=600 (Ω_dot = 0.985647332 deg/day):

```
Omega0 =   0°  →  LST at orbit-plane node range = [17.020, 17.531] h  (EoT envelope ~31 min)
Omega0 =  90°  →  LST at orbit-plane node range = [23.020, 23.531] h
Omega0 = 180°  →  LST at orbit-plane node range = [ 5.020,  5.531] h
Omega0 = 270°  →  LST at orbit-plane node range = [11.020, 11.531] h
Total drift in every case: 0.0000 h/day  (NOT 4 min/day)
```

The orbit-plane-node LST of a true dawn-dusk SSO is **constant** (with a
±~16 min EoT envelope and no secular drift). The 015 claim is **wrong by
construction**.

---

## 4. The Mutant Battery

The 015 README asserts "the dawn-dusk SSO design only fixes the LST
modulo the sidereal-vs-solar differential." A genuinely SSO-locked LST
should be **invariant under changes in dΩ/dt, sign conventions, year
type, and the SSO target itself** (the SSO is a *definition*, not a
thing that "tracks" the Sun by approximation). Mutants:

### M1: Set `dΩ/dt = 0` (no J2, inertially fixed orbit plane)

For an orbit plane inertially fixed, the LST at the orbit-plane node
should drift at the Sun's mean motion: `dLST/dt = −dα_sun/dt = −0.9856
deg/day = 4 min/day WEST`. **This IS the 015 headline number.**

But the 015 code's `lst_at_node_at_t` returns the LST at the **fixed
ground point** `site_lon`, which is independent of `dΩ/dt`. The mutant
is therefore **invisible to the 015 measurement**. Confirmed numerically:
`lst_at_node_hours(t, site_lon*DEG)` returns the same 24h/year wraparound
regardless of whether J2 is enabled in the eclipse constraint.

The 015 measurement **cannot detect the M1 mutant**, so it cannot
distinguish a true SSO from an inertially fixed orbit. The claim of
"4 min/day" is therefore a property of the measurement, not the orbit.

### M2: Set `dΩ/dt = 2 × SSO_TARGET` (orbit plane precesses twice as
fast as the Sun)

For a true non-SSO orbit with `dΩ/dt = 1.9712 deg/day`, the LST at the
orbit-plane node drifts at `(1.9712 − 0.9856)/15 = +0.0657 h/day =
+3.94 min/day EAST`.

Again, 015's `lst_at_node_at_t` returns the fixed-ground-point LST, which
is unchanged. **015 cannot detect M2 either.**

### M3: Use sidereal-year rate (360/365.25636) instead of mean-solar-year
rate (360/365.2422) as the SSO target

Sidereal year gives `dΩ/dt = 0.98560912 deg/day`; mean-solar gives
`0.98564733 deg/day`. For a true SSO at sidereal-year rate, the LST
drift at the orbit-plane node would be `(0.98560912 − 0.98564733)/15 =
−2.55e-6 deg/day = −0.15 s/day` — a slow westward drift, present but tiny.

The 015 measurement returns the same 24 h/year wraparound at the fixed
ground point. **015 cannot detect M3.**

### M4: Wrong sign conventions (negate `dΩ/dt` or `dα_sun/dt`)

For a retrograde precession, `dLST/dt = (−0.9856 − 0.9856)/15 = −0.131
h/day = −7.88 min/day WEST`. Twice the rate in the wrong direction.

015's measurement is again unchanged. **015 cannot detect M4.**

### M5: Remove C4 insertion convention (let the geodetic node longitude
vary with the orbit)

If we drop C4 and let the geodetic node longitude vary with the orbit,
the LST at the orbit-plane node becomes 0 min/day for a true SSO and
−0.9856 deg/day for an inertially fixed orbit. The 015 code does not
have this measurement, but the lab_utils machinery does — verified in
Section 3.

**Summary:** The 015 measurement is *conventionally* 4 min/day for ANY
orbit at ANY altitude from ANY launch site. It is not an orbit property
and cannot validate (or invalidate) the SSO design. The mutants M1-M4
all pass undetected, demonstrating that the test
`test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` (which only
checks `max_dist_from_18 > 3 h`) is **blind to the physics it claims to
verify**.

---

## 5. The station-keeping conclusion (revisited)

The 015 README and `localdocs/knowledge/dawn-dusk-sso.md` conclude that
"a multi-year dawn-dusk mission needs station-keeping" based on the
4 min/day drift. If the 4 min/day were real, the implied Δv budget would
be ~200 m/s/year to compensate (24 h/year × ~30 m/s per minute of
equatorial plane change at SSO altitude is on the order of 100s of m/s).
**No real LEO SSO mission has ever reported a 200 m/s/year station-keeping
budget** — Sentinel-1's published budget is 5-15 m/s/year (dominated by
drag compensation, NOT LST drift), and Landsat-8 is similar.

The 015 claim is contradicted by operational mission data.

The *correct* station-keeping framing is:
- The LST at the orbit-plane node is bounded by the **EoT envelope**
  (peak-to-peak ~16-20 min/year, max daily drift ~0.5 min/day, periodic).
- The secular LTAN walk comes from **drag-induced RAAN walk** (dominant
  at SSO altitudes; ~10-50 m/s/year) and **Lunisolar + SRP** (~0.01-0.05
  deg/day, ~few min/year).
- A 4 min/day figure is the **sidereal-vs-solar-year differential** in
  Earth-rotation physics, not an orbit property.

---

## 6. What "Exp 012 +2.2 deg LST drift per year" actually represents

Grep across the entire lab for "+2.2 deg LST drift" finds the phrase ONLY
in Exp 015 artifacts. `orbitClasses/experiment.py` (Exp 012) does have a
related quantity — `Molniya apogee_longitude_drift` of `−1.133368081634
deg/day` westward — but that is for a **Molniya orbit** (i=63.4°, e=0.74,
a=26561 km), not an SSO, and it is `deg/day` not `deg/year`. The
`0.02683 deg/day` GEO "stationarity residual" is also irrelevant. There is
no `2.2 deg/year` figure in Exp 012.

**The "Exp 012 +2.2 deg LST drift per year at SSO 600 km" citation in
`dawnDuskSSO/README.md`, `experiment.py`, `results.json`, and
`localdocs/knowledge/dawn-dusk-sso.md` is fabricated.** It is not in the
Exp 012 source code, README, results.json, or knowledge note. It does
not exist anywhere in the lab except as a self-reference inside Exp 015.

(For completeness: 0.9856 deg/day × 365.2422 d / 160 = 2.25 deg/year,
which is where "2.2" might have come from if someone divided the
sidereal-solar differential by 160 by mistake. Or it could be a
re-statement of "4 min/day" in different units. In any case, the
attribution is not in Exp 012.)

---

## 7. Summary of the adversarial case

| Question | Answer |
|---|---|
| Which interpretation of the 015 formula gives 4 min/day? | **Interpretation C** (dΩ/dt mislabeled as Earth's sidereal rotation, d(Subsolar)/dt mislabeled as 360.0). |
| Is Interpretation C the physics of a true dawn-dusk SSO? | **No.** It is the geodetic-node-longitude drift of a non-SSO orbit, with the LST taken at a fixed ground point. |
| What does the 015 code actually measure? | LST at the launch-site geodetic longitude (fixed by C4). This wraps 365 times/year because the site is fixed and the Earth rotates. |
| What is the LST drift at the orbit-plane node of a true dawn-dusk SSO? | **0 min/day** (modulo a bounded ±0.5 min/day daily EoT swing and a cumulative ±16-20 min/year EoT envelope, verified directly with `lab_utils`). |
| Is the 4 min/day "station-keeping" conclusion correct? | **No.** A 4 min/day drift would require ~200 m/s/year to compensate; real Sentinel-1 / Landsat-8 budgets are 5-15 m/s/year. The 015 number is contradicted by operational data. |
| Does the "Exp 012 +2.2 deg LST drift per year" citation exist? | **No.** Fabricated. Grep across the lab finds the phrase only in Exp 015 artifacts. |
| Does the test `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` catch this? | **No.** The test only asserts `max_dist_from_18 > 3h`, which is true for the correct EoT-driven behavior but says nothing about the drift rate. The test is silent on the magnitude, sign, and physical interpretation. |

---

## 8. Recommended corrections (in order of severity)

1. **Red.** The 015 README, `experiment.py`, `results.json` findings, and
   `localdocs/knowledge/dawn-dusk-sso.md` all assert "4 min/day LST
   drift at the ascending node of a dawn-dusk SSO" as a finding. This
   is wrong by ~8× (actual daily drift at the orbit-plane node of a
   true SSO is 0 min/day, bounded by ±0.5 min/day EoT) and ascribes
   an Earth-rotation property to the orbit. **Correct the finding and
   the headline numbers.**
2. **Red.** The 015 artifacts cite an "Exp 012 +2.2 deg LST drift per
   year at SSO 600 km" finding that does not exist. **Remove the
   citation or replace with the correct Exp 012 disclosure** (which
   would be the Molniya apogee drift of −1.133 deg/day, not a 2.2
   deg/year SSO quantity).
3. **Yellow.** The test
   `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` should be
   augmented with:
   - `max|dLST_orbit_node/dt| < 1 min/day` (correct EoT-bounded drift)
   - `|LST(year_end) − LST(year_start)| < 1 h` (cumulative EoT envelope)
   Both are required to make the test mean what it claims.
4. **Yellow.** The station-keeping framing in the README "Next
   Question" should be reframed to "EoT envelope + drag-induced RAAN
   walk + Lunisolar/SRP contribution", as recommended in
   `localdocs/reports/audit-015-follow-up-candidates-2026-08-29.md`
   (the prior audit). A station-keeping experiment built on the
   4 min/day drift would compute ~200 m/s/year budgets that don't
   match real missions.
5. **Green (no action).** The feasible-set cardinality (266-295
   components per altitude, monotone in h) and the held-out equinox
   finding (36.7 vs 11.9 feasible/day) are unaffected by the LST-drift
   interpretation bug. They remain valid findings of Exp 015.
6. **Green (no action).** The 015 code correctly identifies that
   the LST at a fixed geodetic longitude wraps 365 times/year and
   uses this to set up the 12-hour bug fix and the L1 closed-form
   identity test. The bug it fixes (12-hour `subsolar_lon_rad`
   sign error) is real, and the fix is correct.

---

## 9. Audit trail

- **Audit performed:** 2026-08-29 (read-only, no files modified except this report)
- **Files examined:**
  - `research/orbital-mechanics/experiments/dawnDuskSSO/README.md`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/tests/test_dawn_dusk_sso.py`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/results/results.json`
  - `localdocs/knowledge/dawn-dusk-sso.md`
  - `localdocs/reports/audit-015-follow-up-candidates-2026-08-29.md`
  - `src/lab_utils/earth_frames.py` (`lst_at_node_hours`, `subsolar_lon_rad`,
    `gmst_rad_iau1982`, `eci_to_ecef`, `node_lon_from_raan_gmst`)
  - `src/lab_utils/orbits.py` (`SSO_TARGET_DEG_DAY`, `j2_nodal_rate_rad_s`,
    `sso_inclination_rad`)
  - `research/orbital-mechanics/experiments/orbitClasses/experiment.py`,
    `README.md`, `results/results.json`, `localdocs/knowledge/orbit-classes.md`
    (confirmed no "+2.2 deg LST drift per year" finding)
- **Numerical verification:** ran `lab_utils.earth_frames.lst_at_node_hours`
  with (a) the 015 setup (fixed geodetic longitude = site_lon) and (b) a
  perfect SSO (Ω_dot = 0.985647332 deg/day, geodetic node longitude free)
  over 2026. Case (a) gives 24 h/year wraparound with 4 min/day average
  drift; case (b) gives 0 min/day drift with ±16-20 min/year EoT envelope.
  Both are reproducible directly from the lab_utils interface.
- **Cross-check against prior audit:** this report reaches the same
  numerical conclusion as
  `localdocs/reports/audit-015-follow-up-candidates-2026-08-29.md` §1.2-1.3
  (max daily drift 0.498 min/day, cumulative ±20 min/year). The two
  reports differ in focus: the prior one audits the follow-up candidates
  in light of the bug; this one attacks the headline claim itself and
  dissects the four-interpretation matrix.

---

## 10. Verdict (restated for the host)

**ADVERSARIAL-VERDICT: RED**

**015-CLAIM-SURVIVES-MUTANTS: no**

**EXACT-SOURCE-OF-ERROR:** The 015 README writes
`dOmega/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day` and
identifies this as "the LST at the ascending node of a dawn-dusk SSO".
The two subtracted terms are NOT the SSO nodal precession rate and NOT
the Sun's apparent RA rate; they are the Earth's sidereal rotation rate
(`dGMST/dt = 360.9856 deg/day`) and a synthetic 360.0 deg/day that
appears nowhere in the lab's Sun or Earth-rate canon. The difference
`0.9856 deg/day` is the textbook sidereal-vs-solar-day differential,
which is a property of Earth rotation about its own axis vs. Earth orbit
around the Sun — the same 4 min/day that makes solar time and sidereal
time diverge by 4 minutes per day at every ground point on Earth. The
015 code's C4 insertion convention (`Ω = GMST + site_lon`) forces the
geodetic node longitude at insertion to be the launch-site longitude
(constant), so the code measures the LST at a fixed ground point, where
the 4 min/day drift is real but irrelevant. The correct physical
quantity — the LST at the **orbit-plane node** of a true dawn-dusk SSO —
is constant to first order, with a bounded ±0.5 min/day daily EoT swing
and a ±16-20 min/year cumulative EoT envelope that returns to zero
after one mean-solar year. The "Exp 012 +2.2 deg LST drift per year at
SSO 600 km" cited as corroboration in `README.md` L87, `experiment.py`
L883, `results.json` finding #1, and `localdocs/knowledge/dawn-dusk-sso.md`
L48-49 + L179 does not exist anywhere in the lab except as a
self-reference inside Exp 015 itself.

**CORRECT-MAGNITUDE:** **0 min/day drift at the orbit-plane node of a
true dawn-dusk SSO**, with a daily swing bounded by the EoT rate
(peak |dLST/dt| ≈ 0.5 min/day) and a cumulative envelope of ±16-20
min/year that is periodic, not secular. This is the textbook result
for a sun-synchronous orbit and matches the operational behavior of
real Sentinel-1 / Landsat-8 missions (~5-15 m/s/year station-keeping
budget, dominated by drag, NOT by LST drift). The 015 number
(0.9856 deg/day ≈ 4 min/day) is correct for a fixed ground point on
Earth but irrelevant to the orbit-plane-node LST of a dawn-dusk SSO;
the ascription of that number to "dawn-dusk SSO physics" is the bug.
