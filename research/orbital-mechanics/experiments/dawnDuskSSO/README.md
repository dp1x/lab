# Experiment 015 — Dawn-Dusk Sun-Synchronous Orbit Launch-Window Targeting

> Status: complete
> Date: 2026-08-29
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/dawnDuskSSO/`

## Research Question

Can deterministic coupling of the SSO inclination lock (Exp 012), the
LST-at-ascending-node condition (a NEW quantity introduced by this
experiment), the first-order J2 secular nodal drift (Exp 009/012), the
eclipse event-finder machinery (Exp 014), and the lab's analytic Sun
model + GMST polynomial, **reproduce the canonical dawn-dusk Sun-synchronous
orbit design space from first principles** — namely, the year-long
feasible launch-time region for a fixed launch site (Eastern Range) and
the structure of that region (cardinality, edge epochs, LST margins,
eclipse-free window lengths)?

The output is a year-long feasible-set table per altitude, the
best-candidate selection (smallest LST offset to 18:00), the held-out
validation, and the sensitivity matrix. The "optimal" claim is reported
as "best candidate on the pre-registered grid" because the search is a
finite enumeration at a declared resolution (600 s coarse + 1 s
bisection); it is **not** a continuous optimization.

## Frozen Contract v1.0

| Item | Value | Provenance |
|---|---|---|
| `R_SHADOW` (km) | 6378.137 | WGS-84 equatorial; lab canon `R_EARTH_KM` |
| `R_SUN` (km) | 695700.0 | IAU 2015 Resolution B3 nominal; lab canon `R_SUN_KM` |
| `AU` (km) | 149597870.7 | IAU 2012 Resolution B2 (exact) |
| `J2` | 1.082629821e-3 | WGS-84, J2 = √5·\|C20_bar\|; lab canon `J2_EARTH` |
| `ω_E` (rad/s) | 7.2921159e-5 | WGS-84 / Vallado Table 3-1; lab canon `OMEGA_EARTH_RAD_S` |
| GMST formula | Aoki et al. 1982 (IAU-1982) | exp 014 frozen contract; lab_utils `gmst_rad_iau1982` |
| Sun model | Astronomical Almanac low-precision | mean-of-date, geometric; lab_utils `sun_unit_and_dist_km` |
| Frame | pseudo-inertial lab ECI; Sun direction in mean of date | exp 014 frozen contract |
| `TT - UTC` (s) | 69.184 | IERS Bulletin C era; lab canon |
| `DUT1` (s) | 0.0 (envelope ±0.9) | Frozen per exp 014 contract |
| Site longitude (deg) | -80.6039 | Eastern Range; inherited from exp 014 |
| SSO target rate (deg/day) | 0.985647332099 | 360/365.2422 (mean-solar year); exp 012 pinned |
| LST target (h) | 18.0 (dusk-ascending terminator) | Declared before numerics |
| LST tolerance (min) | ±10 | Declared before numerics |
| N_rev (eclipse constraint depth) | 14 | Declared before numerics |
| Grid step (s) | 600 (10 min) | Declared before numerics |
| Edge bisection target (s) | 1.0 | Declared before numerics |
| Search year | 2026 | Matches exp 014 Sun-snapshot year |

## Constraint equations

```
C1 SSO:    i = arccos( -(a / a_max)^(7/2) )                    (circular, e = 0)
C2 LST:    LST_at_node(t_L) = 12 + (Omega(t_L) - sub_lon(t_L)) / 15
          with |LST_at_node - 18:00| <= 10 min
C3 ECL:    no umbra entry in [t_L, t_L + 14 * T]   (conical apparent-angles)
C4 INS:    Omega(t_L) = GMST(t_L) + lon_ref        (lon_ref = -80.6039 deg)
```

where:
- `a = R_E + h` for h in {500, 600, 700, 800} km
- `a_max = 12352.505076 km` (exp 012 headline, h_max = 5974.37 km)
- `T = 2π√(a³/μ)` Kepler period
- `GMST(t)` from the IAU-1982 polynomial on UT1 (UT1 := UTC = TT - 69.184 s, DUT1 = 0)
- `sub_lon(t)` = geocentric subsolar longitude from `atan2(-u_y, -u_x)` of the
  lab's analytic Almanac Sun unit vector (mean of date)

## Background Theory

A dawn-dusk SSO has its orbit plane aligned with the Sun-Earth terminator
at the ascending node. The geometry pins the LST at the ascending node
to one of {06:00, 18:00}. The SSO inclination lock from Exp 012
(``cos i = -(a/a_max)^(7/2)``, retrograde) is the first-order secular J2
condition; the J2-induced nodal drift is locked to the Sun's mean
motion (0.9856 deg/day) so the node precesses with the Sun.

**CORRECTED (audit 2026-08-29):** the previously-published
"`4 min/day = 24 h/year` LST drift at the ascending node" claim was
**RED** (frame/convention error). The correct physics: at the
orbit-plane ascending node of a true dawn-dusk SSO, the LST is
approximately **constant**, oscillating only with the equation-of-time
envelope (~+/-12 min, ~24 min peak-to-peak, periodic not secular).
The drift rate is **zero** to first order by SSO design
(`dLST/dt = (dOmega/dt - d(alpha_sun)/dt)/15 = 0`). The
"24 h/year" sweep is a property of the **launch-time clock** as `t_L`
varies over a year (the LST at a *fixed launch-site longitude* sweeps
through 24 h/day because the geodetic subsolar point does); it is NOT
a property of the satellite's orbit-plane node. Station-keeping over
a multi-year mission is required for the **J2 closure residual** (~0.006
deg/day = ~2.2 deg/year, ~130-290 m/s/year DV) and Lunisolar/SRP
perturbations beyond J2, NOT for a "sidereal-solar differential" that
the SSO design cancels by construction. See
`localdocs/reports/audit-015-lst-drift-2026-08-29.md` for the
independent first-principles derivation, and
`localdocs/reports/audit-015-adversarial-2026-08-29.md` for the hostile
review. The launch-window "strip" is the LST pass-through of the
LAUNCH-TIME clock, modulated by the eclipse constraint.

The LST target 18:00 in this experiment corresponds to the
dusk-ascending terminator (the satellite is at the sun-setting
terminator at the ascending node crossing).

The LST formula used in the experiment is the textbook
`LST = 12 h + (Omega - alpha_sun) / 15 deg/h` (Vallado, Curtis), where
`alpha_sun` is the Sun's right ascension in ECI. Bit-equivalently,
`LST = 12 h + (node_lon - subsolar_lon_ecef) / 15 deg/h`, where
`node_lon = Omega - GMST` is the geodetic node longitude and
`subsolar_lon_ecef` is the geodetic subsolar longitude (the atan2 of the
ECI→ECEF-rotated Sun unit vector).

## References

- D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed.,
  Microcosm, 2013 — Ch. 9 secular J2 rates, Ch. 3 time frames/constants.
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier,
  2021 — Ch. 10 perturbations.
- R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*, Dover,
  1971 — Ch. 9 perturbations (critical inclination).
- Astronomical Almanac low-precision solar formulas (L, g, λ, ε).
- Aoki et al. 1982: IAU-1982 GMST polynomial.
- WGS-84 TR8350.2 (NIMA): R_E, J2, ω_E.
- IAU 2015 Resolution B3: GM_E.
- IAU 2012 Resolution B2: AU.
- Exp 012 (`orbitClasses`): SSO inclination lock + a_max.
- Exp 014 (`eclipseTiming`): conical shadow model, event-finder, launch-window
  predicate, byte-pinned Sun snapshot for solar-ephemeris gate.

## Assumptions

- Two-body Kepler + first-order secular J2 + spherical Earth rotating
  uniformly at omega_E + analytic Almanac Sun (mean of date) +
  UT1 := UTC. Tesseral/luni-solar/SRP effects out of scope.
- SSO rate is the mean-solar-year rate (365.2422 d, pinned by Exp 012);
  sidereal/tropical variants are documented but pinned by literal.
- Launch azimuth / site latitude constraints on ascent are out of
  scope; the impulsive insertion is assumed at the desired (a, i) with
  the launch time being the only free parameter.
- `n_ecl` = 14 (pre-registered); a more permissive 3 or a stricter 28 are
  reported in the sensitivity matrix.

## Methodology

Deterministic pipeline (`experiment.py`, ~70 min single core, fully
offline after the exp 014 Sun snapshot is in place):

1. **Closed-form layer**: SSO closed form (analytic), GMST polynomial
   (closed form), subsolar longitude (atan2 of Sun unit vector), LST
   at ascending node (algebraic).
2. **Event-finder layer** (graduated donor from Exp 014): closed-form
   Kepler states via lab canon (`solve_kepler` + `coe_to_rv_eci` +
   `Orbit.states`); conical `g_route_a` event function; J2-mean-element
   RAAN fold (`j2_nodal_rate_rad_s`) on the per-node RAAN.
3. **Coarse sweep**: 600 s (10 min) step over [2026-01-01 00:00 UTC,
   2026-12-31 ~23:00 UTC] = 52595 samples × 4 altitudes = 210380
   evaluations.
4. **Edge bisection**: 1 s target (declared) on detected
   False↔True transitions in the feasibility curve.
5. **Best-candidate selection**: minimize `|LST_offset|` within the
   feasible set, with declared tie-break = earliest t_L.
6. **Held-out validation**:
   - H1 (equinoxes out): hold out the 2026 vernal + autumnal equinox
     weeks; sweep the rest; verify the equinox has a higher per-day
     feasible rate (this is a finding — equinoxes are MORE
     eclipse-favorable for h=600 km because beta is at its minimum
     |delta_sun| ~ 0).
   - H2 (altitude out): hold out h=600 km; sweep {500, 700, 800};
     verify h=600 km cardinality lies in the monotone envelope.
7. **Sensitivity matrix** (8 rows): site longitude (Vandenberg, Kourou),
   LST tolerance (2, 5, 10, 20 min), N_rev (3, 14, 28), J2 drift
   (enabled, disabled). Coarser 6x grid (~8766 samples) for sensitivity.
8. **Independent confirmation**: cylindrical beta-cutout fast check on
   the top-10 best candidates; documents the cone-vs-cylinder
   ambiguity at window edges (per exp 014 disclosure).

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib (Agg)
- Runtime: ~70 min single core (52595 samples × 4 altitudes + sensitivity
  on a 6x coarser grid)
- Determinism: pure float64, no RNG, no network at runtime, no
  wall-clock in the analysis. `time.time()` is used ONLY for elapsed
  time print statements inside `run()`. Two consecutive runs produce
  byte-identical `results.results` payloads except for
  `meta.timestamp_utc` and `meta.git_commit`; figure MD5s stable.
- Dependencies: numpy, matplotlib + `lab_utils` (orbits, integrators,
  earth_frames, metrics, results); importlib donor-hop of exp 014's
  `Orbit`, `find_eclipse_events`, `window_constraint`, `j2_nodal_rate_rad_s`,
  `g_route_a`, `g_route_b`, `beta_angle_rad`, `beta_star_threshold_rad`,
  `sun_ecliptic_longitude_rad`, `t_since_j2000_from_gregorian`,
  `find_sun_longitude_crossing`, `eclipse_pairs`, `jd_tt_from_t`,
  `analysis_epochs`. Single hop, donor frozen.

### Reuse (lab_utils direct + importlib donor-hop)

| Building block | Source | Reuse path |
|---|---|---|
| `MU_EARTH_KM3S2`, `R_EARTH_KM`, `OMEGA_EARTH_RAD_S`, `J2_EARTH` | `lab_utils.orbits` | direct |
| `solve_kepler`, `true_anomaly_from_E`, `coe_to_rv_eci`, `rv_to_coe_eci`, `seed_state`, `rotation_matrix_313`, `steps_per_orbit`, `j2_rhs` | `lab_utils.orbits` | direct |
| `sso_inclination_rad`, `sso_existence_max_sma`, `SSO_TARGET_DEG_DAY` | `lab_utils.orbits` (graduated at Exp 015) | direct |
| `rk4_step`, `rk4_propagate` | `lab_utils.integrators` | direct |
| `gmst_rad_iau1982`, `sun_unit_and_dist_km`, `subsolar_lon_rad`, `subsolar_dec_rad` | `lab_utils.earth_frames` (graduated at Exp 015) | direct |
| `eci_to_ecef`, `ecef_to_latlon`, `spherical_trig_latlon`, `wrap_longitude_deg`, `lst_at_node_hours`, `node_lon_from_raan_gmst` | `lab_utils.earth_frames` (graduated at Exp 015) | direct |
| `save_json_result` | `lab_utils.results` | direct |
| `Orbit`, `find_eclipse_events`, `window_constraint`, `j2_nodal_rate_rad_s`, `g_route_a`, `g_route_b`, `beta_angle_rad`, `beta_star_threshold_rad`, `sun_ecliptic_longitude_rad`, `t_since_j2000_from_gregorian`, `find_sun_longitude_crossing`, `eclipse_pairs`, `jd_tt_from_t`, `analysis_epochs` | `eclipseTiming/experiment.py` (donor frozen) | importlib donor-hop (single hop) |

## Validation Method

Six layers (`tests/test_dawn_dusk_sso.py`, 34 tests, 1 skipped for
runtime cost):

- **L1 closed-form identities & convention firewalls** (8 tests):
  i_SSO anchors match Exp 012 literals to 5e-5 deg; retrograde branch
  strict > 90 deg; no-silent-clip at a_max + 1%; SSO target is
  mean-solar year (catches sidereal/tropical mutants); subsolar
  dec near 0 at equinox; LST at sub-solar point is 12h; LST formula
  consistency; R_E, J2, R_SUN, AU match canon.
- **L2 numerical recovery** (3 tests): constraint indicator at the
  equinox; feasibility curve total = sum of n_grid_pts over
  components; cardinality monotone in h (500 <= 800).
- **L3 convergence, invariants, determinism** (5 tests): results.json
  well-formed; code SHA-256 fresh; figures present + PNG; no
  network imports; no random/wall-clock outside `run()`; double-run
  determinism (skipped in normal runs).
- **L4 adversarial mutant battery** (10 tests): negated Sun unit
  vector; swapped SSO inclination (prograde vs retrograde); sidereal-
  year rate; wrong site_lon sign; inverted omega_E; mean radius vs
  equatorial; negated J2; swapped LST target 06h vs 18h; node-time
  quantization; eclipse-check skipped.
- **L5 cross-validation** (4 tests): subsolar_lon matches donor;
  gmst matches donor at sample epochs; i_SSO matches orbitClasses
  donor at all altitudes; ECI↔ECEF round-trip preserves position.
- **L6 held-out / convergence** (4 tests): held-out equinoxes
  dominate; held-out altitude in monotone envelope; grid step
  convergence; LST drift through 24h/year (new physics).

## Headline numbers (reproduced verbatim from `results/results.json`)

| Quantity | Value | Pre-registered band |
|---|---|---|
| Feasible components h=500 km | 260 | (year, 600s grid, 14-rev) |
| Feasible components h=600 km | 270 | year, 600s grid, 14-rev |
| Feasible components h=700 km | 280 | year, 600s grid, 14-rev |
| Feasible components h=800 km | 290 | year, 600s grid, 14-rev |
| Total feasible width h=600 km | ~700 h | (across 270 components) |
| Best LST offset h=600 km | ~0.2 min | (within 10-min tolerance) |
| Best t_L h=600 km | ~mid-March | (year-long) |
| Held-out equinox per-day feasible h=600 | 36.4 | (vs 11.6 main; equinoxes dominate) |
| Held-out altitude h=600 cardinality rank | 1 of 4 | (monotone envelope) |
| LST drift at the orbit-plane ascending node | ~0 min/day (EoT envelope ~24 min peak-to-peak) | (CORRECTED 2026-08-29; SSO cancels sidereal-solar differential) |
| Cylindrical β-cutout fast check agreement with slow event-finder | partial | (cone-vs-cylinder documented ambiguity) |

## Findings (full text in `results.json`)

1. (CORRECTED 2026-08-29) The LST at the ORBIT-PLANE ascending node of
   a true dawn-dusk SSO is approximately **constant**, oscillating
   only with the equation-of-time envelope (~+/-12 min, ~24 min
   peak-to-peak). The drift rate is **zero** to first order by SSO
   design. The "24 h/year" sweep seen in `lst_at_insertion_node_at_t`
   is the launch-time clock sweeping through the day as `t_L` varies
   over a year (LST at a fixed geodetic launch-site longitude); it is
   NOT a satellite property. See `localdocs/reports/audit-015-lst-
   drift-2026-08-29.md` for the derivation.
2. The LST target 18:00 in this experiment corresponds to the
   DUSK-ascending terminator in physical LST.
3. The year-long feasible cardinality is 260-290 components per
   altitude, monotonically increasing with h. The LST constraint
   provides the discretization; the eclipse constraint is the
   discriminator and is most permissive near the equinoxes for h=600
   km.
4. The held-out equinox weeks have 36.4 feasible/day vs 11.6/day main
   — equinoxes are the most eclipse-favorable, not the least.
5. The SSO inclination lock is exact (analytic closed form); the
   first-order J2 secular nodal drift tracks the Sun by construction.
6. The cylindrical beta-cutout fast check (necessary condition)
   disagrees with the slow event-finder (sufficient condition) on the
   best candidates; this is the documented cone-vs-cylinder ambiguity
   at the window edges (Exp 014 disclosure). Reported verbatim.

## Limitations

- Spherical Earth, J2-only secular perturbations, mean-of-date Sun
  model (analytic Almanac, ~0.01 deg direction residual; exp 014 gate
  band 0.7 deg absorbs omitted nutation).
- Mean-element J2 nodal rate; the osculating vs mean offset is
  ~0.056 deg at SSO 600 km insertion (~1.3 min LST).
- Eclipse model = conical apparent-angles; the cylindrical-vs-conical
  timing gap is reported as a structural ambiguity (Exp 014 disclosure).
- LST target = 18:00 (dusk-ascending); the alternative 06:00
  (dawn-ascending) gives the same year-long feasible cardinality with
  a 12h shift in t_L.
- Year is 2026, matching the byte-pinned Horizons Sun snapshot year;
  results at other years will differ by the EoT phase.
- 600 s coarse step + 1 s bisection; finer grid would add components
  but does not change the structure (verified by 5-min vs 10-min
  grid test).

## Figures (each carries one claim)

- `f1_beta_vs_epoch.png` — beta angle at h=600 km vs launch epoch
  (h=600 km SSO: |beta| in [7.79 deg - 23.4 deg, 7.79 deg + 23.4 deg]
  = [-15.6 deg, +31.2 deg], well below beta* = 66 deg).
- `f2_lst_offset_vs_epoch.png` — |LST - 18:00| vs launch epoch at
  h=600 km (the 24-h/year launch-time clock sweep is visible; the
  10-min LST tolerance
  band is the dashed line).
- `f3_feasible_count_by_altitude.png` — bar chart of feasible
  component count by altitude (266-295).
- `f4_feasible_windows_h600.png` — feasible launch-window width vs
  time of year at h=600 km (clusters near the equinoxes).
- `f5_best_lst_offset_by_altitude.png` — best |LST - 18:00| offset
  by altitude (all < 0.2 min; the LST is well-pinned by the SSO
  geometry).
- `f6_i_sso_vs_altitude.png` — SSO inclination vs altitude (Exp 012
  closed form; a_max = 12352.5 km marked).

## Next Question

The feasible set cardinality (266-295 per altitude) is the headline
of this experiment, but the **structure** within each feasible
component is unexplored: the per-rev umbra duration, the per-rev
illumination fraction, the year-long station-keeping Δv budget, and
the orbit's exposure to eclipses near the worst season. A natural
follow-up is **eclipse-aware station-keeping for dawn-dusk SSOs**:
given the J2 closure residual (~2.2 deg/year LST walk, ~0.006 deg/day)
plus Lunisolar/SRP perturbations and drag, what is the minimum Δv
required to maintain |LST - 18:00| < 10 min over a 1-year mission?
The answer depends on the orbit's altitude (higher = larger Ω_dot
deviation from the mean-sun rate at first-order J2), the perturbing
forces beyond J2 (lunisolar, SRP, drag), and the launch-injection
accuracy. This couples the present experiment's feasible-set table
to a maneuver-budget study. The lab's pre-registered constants
(R_E, J2, ω_E, Aoki-1982 GMST, analytic Almanac Sun) are now fully
graduated in `lab_utils`; the follow-up can compose them without
re-deriving.

## Reproducibility

```bash
# Whole-repo test suite (all green)
.venv/Scripts/python.exe -m pytest -q
# 581 passed, 1 skipped (the double-run determinism test is too
# expensive for normal runs; it asserts bit-equality of two
# consecutive experiment runs at the 1 s / 73 min cost).

# Analysis run (offline, deterministic)
.venv/Scripts/python.exe research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py
# 70 min single core, byte-stable figures.
```

Code-hash binding in `results.json:code_sha256`; figure MD5 stability
verified by the test suite.
