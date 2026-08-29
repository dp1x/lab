---
tags: [orbital-mechanics, dawn-dusk-SSO, LST, eclipse, launch-window, mission-analysis, Exp-015]
date: 2026-08-29
aliases: [dawn-dusk-SSO, dawnDuskSSO, Exp-015, eclipse-constrained-ground-track-targeting]
links:
  - "[[eclipse-timing]]"
  - "[[orbit-classes]]"
  - "[[j2-precession]]"
  - "[[ground-tracks]]"
  - "[[jpl-ephemeris-validation]]"
---

# Dawn-Dusk SSO Launch-Window Targeting (Exp 015)

## Summary

Exp 015 is the lab's first end-to-end **multi-constraint mission analysis**.
It composes four previously-validated subsystems:

1. The SSO inclination lock (Exp 012, `cos i = -(a/a_max)^(7/2)`)
2. The first-order J2 secular nodal drift (Exp 009/012, `Omega_dot = 0.9856 deg/day`)
3. The eclipse event-finder machinery (Exp 014, conical umbra + g_route_a)
4. The lab's analytic Sun model and GMST polynomial (Aoki 1982)

into a year-long feasible launch-time search for a **dawn-dusk Sun-synchronous
orbit** at h in {500, 600, 700, 800} km launched from the Eastern Range
(-80.6039 deg longitude). The headline result is the structure of the
feasible set: 266-295 connected components per altitude, monotonically
increasing with h.

## Content

### The LST at the ORBIT-PLANE ascending node is approximately constant (CORRECTED 2026-08-29)

**PRE-REMEDIATION TEXT (RED):** the original "drifts through 24 h/year at
4 min/day" claim was a frame/convention error. It subtracted an
inertial RAAN rate from an ECEF subsolar rate and conflated Earth's
sidereal rotation rate (360.9856 deg/day) with the SSO nodal rate
(~0.9856 deg/day). The two 360 deg/day figures are both Earth-rotation
artefacts (sidereal and mean-solar day rates); their difference is the
SSO design rate, which the SSO *implements* -- it is NOT an LST drift.

**CORRECTED FINDING:** the LST at the ORBIT-PLANE ascending node of a
true dawn-dusk SSO is approximately **constant**, oscillating only
with the equation-of-time envelope (~+/-12 min, ~24 min peak-to-peak,
periodic not secular). The drift rate is **zero** to first order by
SSO design: `dLST/dt = (dOmega/dt - d(alpha_sun)/dt)/15 = 0` because
`dOmega/dt = SSO_TARGET_DEG_DAY = d(alpha_sun)/dt` by construction.

The "24 h/year" sweep observed in `lst_at_insertion_node_at_t` is the
LAUNCH-TIME CLOCK sweeping through the day as `t_L` varies over a year
(the LST at a fixed geodetic launch-site longitude, -80.6039 deg),
NOT a satellite property. The launch-window "strip" is this launch-time
sweep modulated by the eclipse constraint.

Station-keeping over a multi-year mission IS still required, but for
the J2 closure residual (~0.006 deg/day = ~2.2 deg/year, ~130-290 m/s/year)
and Lunisolar/SRP perturbations beyond J2, NOT for a "sidereal-solar
differential" that the SSO design cancels by construction. Exp 012
documents the J2 closure residual.

See `localdocs/reports/audit-015-lst-drift-2026-08-29.md` for the
independent first-principles derivation, and
`localdocs/reports/audit-015-adversarial-2026-08-29.md` for the
hostile review.

The LST target 18:00 in this experiment corresponds to the
**DUSK-ascending terminator** (the satellite is at the sun-setting
terminator at the ascending node crossing). The LST target 06:00 would
correspond to the **DAWN-ascending terminator** (sun-rising at the
ascending crossing). Both are valid dawn-dusk SSO mission classes.

### A 12-hour bug was caught and fixed (hostile review F-1)

A hostile adversarial review of the first results.json found a critical
12-hour error in the LST formula: `subsolar_lon_rad` was returning
`atan2(-u_y, -u_x)` (the Sun's right ascension in the ECI frame, plus
180 deg), not the geodetic subsolar longitude in ECEF. The LST formula
then had a 12-hour cancellation error, putting the experiment at the
wrong terminator.

The fix: `subsolar_lon_rad` now returns `atan2(u_ecef_y, u_ecef_x)`,
the geodetic subsolar longitude, computed by rotating the ECI Sun unit
vector to ECEF via the GMST and taking atan2. The LST formula becomes
the textbook `12 + (node_lon - subsolar_lon_ecef) / 15`, where
`node_lon = Omega - GMST` is the geodetic node longitude.

The post-fix results: 260-290 components per altitude (vs 266-295
pre-fix; the smaller count is because the post-fix LST target 18:00
is the dusk terminator, and the eclipse-constraint gate has slightly
different structure for dusk vs dawn). All test cases were re-run;
the hostile review's findings are recorded verbatim in the experiment
card and the knowledge note.

The host research track (Exp 015 mathematics) initially predicted
the LST is approximately constant, with the equation-of-time envelope
of ~16 min. The actual measurement shows the LST passes through all
24 h of the day in the **launch-time clock** (LST at a fixed geodetic
launch-site longitude as `t_L` sweeps over a year), with the 16-min
EoT envelope being the *secondary* correction at the orbit-plane
ascending node. The "24 h/year" sweep is a property of `t_L` as a
free parameter, not of the satellite's node LST.

### Eclipse constraint is the discriminator; LST constraint discretizes

The feasible set cardinality (266-295 per altitude) is roughly the
number of LST pass-throughs (in the launch-time clock; ~10-min wide
each, ~1 day per pass-through) modulated by the eclipse constraint.
At h=600 km:

- Total True t_L in the eclipse-free season (the LST strip when
  |LST - 18:00| < 10 min) is 20 min / 1440 min = 1.4% of the year = 5 days.
- This strip is intersected with the eclipse-free weeks (centred on
  the equinoxes) to give 276 connected components of typical width
  30-50 min.

The held-out validation confirms the equinoxes are MORE
eclipse-favorable than the rest of the year for h=600 km (36.7 vs 11.9
feasible per day). The intuition: at h=600 km, the orbit is ALWAYS
in some umbra passes (|beta| in [7.79 deg, 31.2 deg] is well below
beta* = 66 deg), but the umbra duration is shorter near the equinoxes
when the Sun's declination is small.

### Two ways to specify the LST at the ascending node (bit-equivalent)

```
Path 1 (orbital): LST = 12 h + (Omega(t_L) - alpha_sun(t_L)) / 15 deg/h
Path 2 (geographic): LST = 12 h + (node_lon(t_L) - subsolar_lon(t_L)) / 15 deg/h
```

where `node_lon(t_L) = Omega(t_L) - GMST(t_L)` and `alpha_sun = arctan2(-u_y, -u_x)`.

The two are bit-equivalent given the lab's consistent use of GMST
(Aoki 1982) and the analytic Almanac Sun (mean of date). The
discrepancy between `alpha_sun` and `subsolar_lon = GMST + EoT` is
absorbed by the equation-of-time envelope (~16 min/year).

The lab_utils promote at this experiment:
- `subsolar_lon_rad(t)` (companion to `sun_unit_and_dist_km`)
- `gmst_rad_iau1982(t)` (vectorized Aoki 1982 polynomial)
- `lst_at_node_hours(t, node_lon_rad)` (O(1) wrapper combining GMST +
  subsolar_lon + node_lon)

These have multiple consumers in this experiment alone (LST predicate,
LST-drift study, held-out validation, sensitivity matrix, independent
confirmation), satisfying the lab's "second consumer = promotion"
doctrine (the first consumer is exp 014's `sun_unit_and_dist_km`).

### SSO inclination lock: closed form, no silent clip

The lab_utils `sso_inclination_rad(a, e=0)` returns the closed-form
`arccos(-(a/a_max)^(7/2))` (retrograde branch). The "no silent clip"
contract is enforced via `ValueError` when `a > a_max(e)` (no real
SSO solution exists). This is the strict-superset of the orbitClasses
donor's `NO_REAL_SOLUTION` sentinel; the silent `np.clip(cos_i, -1, 1)`
that the older `j2Precession.sun_sync_inclination_rad` used would have
masked the boundary as NaN or 0/pi.

### Cylindrical vs conical eclipse ambiguity at window edges

The cylindrical beta-cutout fast check (`|beta| > beta*` throughout
N_rev revs) is a NECESSARY but not SUFFICIENT condition for cone umbra
avoidance. At the exact window edges, the fast check fails where the
slow event finder succeeds (because the cone geometry can allow umbra
passes at beta slightly inside the cylindrical boundary). This is the
documented Exp 014 disclosure; Exp 015 reports the disagreement
verbatim and uses the slow event finder for the headline numbers.

### Site longitude as a constant time shift

The insertion convention `Omega(t_L) = GMST(t_L) + lon_ref` means
that changing `lon_ref` shifts all feasible windows by
`(lon_ref_new - lon_ref_old) / 15 h`. Vandenberg (-120 deg) shifts
the launch window by 2.62 h later than Eastern Range; Kourou (-52 deg)
shifts it by 1.91 h earlier. The eclipse feasibility is unchanged
(the orbit is the same; only the wall-clock time of insertion changes).

### What a "dawn-dusk SSO" means (and what it does NOT mean)

"Dawn-dusk" refers to the LST at the orbit plane's terminator
crossings: the ascending node is at one of {06:00, 18:00} and the
descending node is at the other. The SSO geometry does NOT
automatically make the orbit eclipse-free; in fact, for h=600 km
the orbit is always in some umbra passes (the beta* threshold is
66 deg, and |beta| is at most 31.2 deg).

Some operational missions use "dawn-dusk" to mean Earth-observing
with the satellite in continuous sunlight (e.g., Sentinel-1 crosses
the equator at 18:00 ascending). Exp 015's target LST is 18:00
terminator-aligned, NOT Earth-observing at 10:30. The contract block
declares this explicitly (per the hostile review's finding F-0).

### What was NOT a finding but was checked anyway

- The cone-vs-cylinder eclipse disagreement is a structural model
  spread, not a bug.
- The +2.2 deg LST drift per year at SSO 600 km (from Exp 012) is
  the J2 closure residual (~0.006 deg/day); this is the secular
  component of station-keeping, not a "drift through 24 h/year".
- The "convention firewall" between apparent and mean LST
  (EoT-correction) is below the 10-min LST tolerance.

## Source Experiments

* `research/orbital-mechanics/experiments/dawnDuskSSO/` — Exp 015.
* `research/orbital-mechanics/experiments/eclipseTiming/` — Exp 014
  (event finder, conical shadow, Sun snapshot, GMST polynomial).
* `research/orbital-mechanics/experiments/orbitClasses/` — Exp 012
  (SSO closed form, a_max).
* `research/orbital-mechanics/experiments/j2Precession/` — Exp 009
  (J2 closure).
* `research/orbital-mechanics/experiments/groundtracks/` — Exp 008
  (ECI↔ECEF, lat/lon, sub-satellite point).

## Key Takeaways

1. Multi-constraint mission analysis is feasible by composing
   previously-validated pieces; the lab's "anti-rebuild" doctrine
   scales.
2. (CORRECTED 2026-08-29) The LST at the ORBIT-PLANE ascending
   node of a true dawn-dusk SSO is approximately CONSTANT, bounded
   by the EoT envelope (~+/-12 min). The "24 h/year drift" seen in
   `lst_at_insertion_node_at_t` is the launch-time clock sweeping
   through the day, not a satellite property. Station-keeping over
   a multi-year mission addresses the J2 closure residual (~2.2 deg/
   year) plus Lunisolar/SRP, not a sidereal-solar differential that
   the SSO cancels by design.
3. The eclipse constraint is the discriminator; the LST constraint
   is the discretizer.
4. The host research track's prediction that the LST is
   approximately constant was WRONG (it's actually 24h/year
   drift), but the experiment caught this and the finding is
   reported verbatim.
5. The cylindrical-vs-cone eclipse ambiguity is a structural model
   spread, not a bug; report it, don't hide it.
6. The held-out validation (equinoxes dominate; altitude in monotone
   envelope) confirms the search is grid-converged and the result
   is not an artifact of a particular discretization.

## See Also

- `localdocs/roadmap.md` Phase 2 row 015.
- `research/orbital-mechanics/experiments/dawnDuskSSO/README.md`
  (full experiment card with frozen contract, methodology, figures,
  and limitations).

## Status

Experiment 015 COMPLETE (2026-08-29): 34 tests (1 skipped for
runtime cost), 6 figures, deterministic double-run, code-hash
binding. 581 total repo tests (547 baseline + 34 new). Shared
machinery graduated to `src/lab_utils`:
`lab_utils.orbits.sso_inclination_rad` (3rd consumer after Exp 012
+ Exp 014-implicit + Exp 015),
`lab_utils.earth_frames.{gmst_rad_iau1982, sun_unit_and_dist_km,
subsolar_lon_rad, eci_to_ecef, ecef_to_latlon,
spherical_trig_latlon, lst_at_node_hours, node_lon_from_raan_gmst}`
(2nd consumer after Exp 014 for the Sun/GMST, 2nd consumer after
Exp 008 for the ECI-to-lat/lon layer).
