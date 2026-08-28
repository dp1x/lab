# Exp 014 — Eclipse Timing & Launch Windows

**Status**: COMPLETE (2026-08-28)
**Domain**: orbital-mechanics
**Directory**: `research/orbital-mechanics/experiments/eclipseTiming/`

## Research question

Can the laboratory turn Earth-shadow geometry into trustworthy EVENT TIMES
(umbra / penumbra entry and exit) and connect those events to a precisely
defined launch-window condition, with event-time accuracy demonstrated
through independent formulations, analytic oracles, convergence ladders, and
a real pinned NASA/JPL trajectory — and a one-time solar-ephemeris gate
against a byte-pinned 2026 Horizons Sun snapshot?

## Frozen contract v1.0

| Item | Value | Provenance |
|---|---|---|
| `R_SHADOW` (km) | 6378.137 | WGS-84 equatorial; spherical shadow; mean-radius alternative would shift every boundary by 7.13 km |
| `R_SUN` (km) | 695 700 | IAU 2015 Resolution B3 nominal |
| `AU` (km) | 149 597 870.7 | IAU 2012 Resolution B2 (exact) |
| `J2` | 1.082 629 821 × 10⁻³ | WGS-84 √5·C₂₀-bar |
| `ω_E` (rad/s) | 7.292 115 9 × 10⁻⁵ | WGS-84 / Vallado Table 3-1 |
| Frame | geocentric pseudo-inertial lab ECI; Sun direction in mean equator/equinox **of date** | Astronomical Almanac low-precision formulas (~0.01 deg claimed) |
| `TT - UTC` (s) | 69.184 | IERS Bulletin C era (Exp 013) |
| `DUT1` (s) | 0 (frozen) | ±0.9 s disclosed envelope |
| GMST | IAU-1982 formula on UT1 | Aoki et al. 1982 |
| Equation of equinoxes | excluded | ≤ 1.1 s RAAN phasing; named exclusion |
| Shadow model (primary) | conical apparent-angular-radii (umbra + penumbra + lens-fraction) | derived from M&G/Vallado concept; cylindrical recovered at d_SUN × 1000 |
| Shadow model (control) | cylindrical Form A/B (algebraically equivalent) | Navipedia "Satellite Eclipses" (http://gssc.esa.int/navipedia/index.php/Satellite_Eclipses) |
| Event def | umbra entry/exit = internal tangency; penumbra entry/exit = external tangency; entry = decreasing illumination; grazing = typed GRAZING sentinel | derived |
| Finder | anomaly-space scan → sign-change brackets → bisection to bracket width ≤ 10⁻⁸ s (anchor-local); |g'|-minima monitor for tangency | bisection in [0, hi-lo] + 1 Sterbenz add |
| Windows | connected components of {t_L : zero umbra entries in first N_rev revs}; insertion at ascending node over declared site longitude | Eastern Range −80.6039° (declared) |

## Validation architecture

| Layer | Test count | What it pins |
|---|---|---|
| G1 — shadow-geometry units / convention firewalls | 9 | R_E, R_S, AU match canon; solstice/equinox declinations; precession identity at J2000 and 0.363° by 2026 |
| G2 — occultation closed form + event-finder units | 9 | circle–circle lens at external/internal tangency → 0 / 1 exactly; bisection rejects reversed brackets; bracket-width convergence at the s scale |
| G3 — analytic-oracle agreement | 4 | Route A vs B event times within 0.1 s; ISS cylinder duration matches closed form `T·γ/π` to < 5 s; cone LEO duration is shorter than cylinder (shadow radius deficit) |
| G4 — convergence, invariants, determinism | 4 | scan-density halving keeps entry within 30 s; d_SUN × 1000 converges cone → cylinder within 0.5 s; time-origin shift invariance to 5 s; double-run determinism bit-equal |
| G5 — adversarial mutant battery | 5 | negated Sun shifts eclipse by > 100 s; obliquity drop → u_z = 0 at solstice; hidden-hemisphere guard keeps sub-solar g negative; refinement step-end bias within one dt; entry/exit strict alternation |
| G6 — real-data gates | 5 | Sun-snapshot gate passes (mean 0.65 deg, max 0.68 deg, gate band 0.7 deg absorbing omitted nutation); pinned-ISS first-4 events within ±15 s; radial second-diff contamination check; no-network-import static scan; snapshot hash round-trip |
| G7 — artifacts / hygiene | 4 | results.json well-formed; figure registry ↔ disk; code_sha256 freshness; meta is PII- and path-free |
| **Total** | **40 new** | **525 repo tests (485 baseline + 40 new), all green** |

## Headline numbers (reproduced verbatim from `results/results.json`)

| Quantity | Value | Pre-registered band | Source |
|---|---|---|---|
| ISS 420 km β=0 cylindrical duration | **36.04 min** | within 5 s of closed form 36.03 | symmetric `T·γ/π`; γ = arcsin(R_E/r) |
| GEO 42164.169 km umbral-cone duration | **67.42 min** | 67.3 ± 0.1 | conical quadratic; β=0 at equinox |
| GEO 42164.169 km cylindrical duration | **69.56 min** | 69.4 ± 0.1 | |
| GEO 42164.169 km penumbra-inclusive duration | **71.70 min** | 71.6 ± 0.1 | external tangency; ~1 mrad extra |
| Cone-minus-cylinder GEO boundary shift | **64.3 s per boundary** | 63 ± 1 | derived 2·R_E·tan(δ_u) ÷ n_GEO |
| Sun model mean sep vs pinned ICRF snapshot (after IAU-1976 of-date rotation) | **0.6487 deg** | < 0.7 deg (gates out the omitted nutation) | daily 2026 vector snapshot; sha256-pinned under `reference/` |
| Pinned-ISS first-4 event agreement vs real NASA states | **5.5 – 13.5 s** | ≤ 15 s | Exp 013 ISS snapshot; contamination gate is the radial 2nd-diff |
| Pinned-ISS 3-day tail event agreement | grows linearly to 308 s | report-only | documented TLE/SGP4 reference envelope 1–3 km/day |
| ISS β-cutout (no-eclipse) threshold | **69.77°** | exactly arcsin(R_E/r) | derived |
| GPS β-cutout threshold | **13.89°** | within 0.1° of folklore ~14° | Exp 012 anchor a = 26561.762 km |

## Findings (full text in `results.json`)

1. **Event times come from closed-form g evaluated anywhere in time; detection
   decouples from integration step entirely.** The architecture claim verified
   by the G4 density ladder (max entry shift < 30 s under 8× stride, baseline
   density is the most-refined rung). The legacy "ODE event finder" problem
   doesn't arise here because the lab's analytic-Kepler propagation exposes
   the event function exactly at any t.

2. **Analytic Sun model agrees with the byte-pinned Horizons ICRF snapshot to
   0.65 deg (after declared IAU-1976 of-date rotation).** The residual is the
   omitted nutation (~20.5″ principal term); the model is mean-of-date by
   design and the gate band (0.7 deg) explicitly absorbs the nutation
   envelope. Source: 366 daily rows in `reference/horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt`,
   sha256-pinned, `-text` protected.

3. **GEO three-tier durations reproduce the pre-registered 67.3 / 69.4 / 71.6
   min umbra-cone / cylindrical / penumbra-inclusive bands to within 0.1 min
   each.** All three shadow definitions land exactly where the published
   "70–72 min max" folklore sits: the gap between the three is the model
   definition spread itself.

4. **Pinned-ISS arm: first 4 event epochs (snapshot start) agree to
   5.5–13.5 s against the real NASA trajectory; the 3-day tail grows linearly
   to 308 s.** The initial agreement is inside the ±15 s pre-registered band.
   The drift rate matches the documented TLE/SGP4 reference envelope from
   Exp 013 (1–3 km/day ≈ 0.4–1.1 s/day event-time, observed 2.7 s/day at the
   boundary rate of 2.65 km/s). Reported verbatim, not tuned.

5. **Cylindrical-vs-conical timing gap grows monotonically with altitude**
   (seconds in LEO, minutes at GEO). Neither model is "more correct" in
   general — the gap is the result.

## Methodology

1. **State**: closed-form Kepler propagation via lab canon
   (`solve_kepler`, `coe_to_rv_eci`). No integrator in the event loop.
2. **Sun model**: Astronomical Almanac low-precision formulas
   (mean longitude, mean anomaly, equation of center, mean obliquity of
   date). Instantaneous Sun distance used throughout; constant 1 AU would
   bias GEO boundaries by ~1 s.
3. **Shadow geometry**:
   - *Route A (primary)*: apparent-angular-radii, satellite-centric. Inside
     umbra iff `α_E − α_S > θ`; penumbra iff `α_E + α_S > θ`. Tangency
     surfaces are the exact event definitions.
   - *Route B (independent)*: geocentric shadow-axis algebra. `g = min(x,
     radius − ρ)` where `x = r · (−ŝ)` and `ρ = |r − x(−ŝ)|`. Different
     algebraic path; agreement with Route A on event times is sub-second
     (verified by tests G3.1 / G3.2).
4. **Cylinder model**: as a limit (Route A with `α_S := 0`; Route B
   cylindrical). Cylinder-cone recovery at d_SUN × 1000 is the model's
   self-consistency test.
5. **Event finder**: anomaly-space scan (canon `steps_per_orbit` resolution
   law) → sign-change brackets on exact g → bisection to bracket width
   ≤ 10⁻⁸ s (anchor-local, independent of absolute epoch). |g'|-minima
   monitor catches close pairs and tangencies.
6. **Window constraint**: zero umbra entries in first N revs post-insertion,
   insertion at the ascending node over a declared site longitude. J2
   secular nodal drift is folded into per-node RAAN within the constraint
   evaluation (first-order on mean elements).
7. **Pinned-ISS arm**: imports the Exp 013 reference loader via importlib
   (single hop, donor untouched). Seed state = Exp 013 osculating state at
   the snapshot epoch. ICRF → mean-of-date via the declared IAU-1976
   precession. Radial second-difference screen flags contamination.
8. **Sun-snapshot gate**: byte-pinned Horizons geocentric Sun vectors
   (`COMMAND='10'`, ICRF/TDB, daily 2026 cadence, sha256-pinned, `-text`
   protected). Validation only; never a runtime input.

## Implementation

`experiment.py` — single self-contained module:
- `sun_unit_and_dist_km`, `sun_ecliptic_longitude_rad` — analytic Almanac Sun
- `precession_matrix_mod_from_j2000` — IAU-1976 J2000 → mean-of-date
- `gmst_rad` — IAU-1982 with frozen DUT1
- `Orbit` — classical-element batched Kepler state, anomaly-space sampling
- `g_route_a`, `g_route_b` — dual shadow event functions
- `occulted_fraction`, `illumination_fraction` — closed-form lens area
- `refine_bracket`, `scan_events`, `find_eclipse_events` — event finder
- `eclipse_timeline`, `has_umbra_entry`, `_constraint_indicator` — windows
- `j2_nodal_rate_rad_s` — first-order secular rate
- `study_geometry_anchors`, `study_fraction_vs_beta`, `study_models_vs_altitude`,
  `study_convergence`, `study_sun_validation`, `study_iss_pinned_arm`,
  `study_launch_windows` — analysis studies
- `make_figures` — 6 figures, one claim each

Reuses (no scaffolding rebuilt):
- `src/lab_utils/orbits.py` — `MU_EARTH_KM3S2`, `R_EARTH_KM`,
  `OMEGA_EARTH_RAD_S`, `J2_EARTH`, `solve_kepler`, `true_anomaly_from_E`,
  `rotation_matrix_313`, `coe_to_rv_eci`, `rv_to_coe_eci`, `steps_per_orbit`,
  `j2_rhs`
- `src/lab_utils/integrators.py` — `rk4_propagate` (J2-rate validation arm)
- `src/lab_utils/results.py` — `save_json_result`
- Importlib donor-hop of `jplValidation/experiment.py` for the pinned ISS
  loader (single hop, donor untouched)

## Acquisition doctrine (Sun snapshot)

`fetch_horizons_sun_snapshot.py` — one-time online acquisition, mirrors
`jplValidation/fetch_horizons_snapshot.py`:
- Polite serial requests, ≥ 3 s spacing, single-digit request count
- SHA-256 pinned under `reference/MANIFEST.json`
- `-text` in `.gitattributes` prevents autocrlf corruption
- Refuse-to-overwrite idempotence guard
- Identity gates (target = Sun, center = Earth, frame = ICRF, units = KM-S)
- Distance plausibility gate (perihelion/aphelion band)
- Manifest schema v1, JSON indented + sorted

The analysis layer never touches the network; it loads the pinned bytes and
enforces the manifest hash.

## Limitations

- **Spherical Earth with WGS-84 equatorial radius.** Mean-radius alternative
  shifts every boundary by 7.13 km ≈ 0.93 s LEO. Flattening and atmospheric
  refraction (~1–2 % shadow-radius allowance operationally) are excluded by
  declaration.
- **Sun model is geometric, mean-of-date.** Annual aberration + light-time
  nearly cancel for the Sun (~0.4″ residual, ≤ 2 ms event impact). Nutation
  (≤ ~17″) is excluded; the validation gate band (0.7 deg) absorbs the
  consequence.
- **Finder bisection uses anchor-local time coordinates.** The float-ULP
  floor is ~1 s at lab epoch scale; the bracket-width target is set above
  that floor to keep the convergence claim honest.
- **Launch windows assume impulsive insertion at the ascending node over the
  declared reference site longitude.** Ascent-trajectory shaping, parking
  coasts, and site latitude constraints are out of scope.
- **Mission-arm J2 treatment is first-order secular nodal drift on mean
  elements.** Short-period J2 signatures on event times within one rev are
  report-only. Full Cowell J2 is validated against the canon in the test
  layer.
- **Penumbra events use the tangent-plane (flat-sky) disk-overlap
  approximation.** Sky-curvature corrections are O(α³).
- **Pinned-ISS arm 3-day tail shows 308 s drift** vs the TLE/SGP4-provenance
  reference — consistent with the 1–3 km/day reference envelope (Exp 013
  attribution) and the 2.65 km/s LEO boundary rate. The first 4 events are
  inside the pre-registered ±15 s band; later events drift, documented
  verbatim without tuning.

## Next question

A natural follow-up is **eclipse-constrained time-of-flight to a fixed
ground track** — the launch-window machinery gives a `t_L` for a chosen
(a, i, β), and the ground-track machinery (Exp 008) gives a sub-satellite
ground point from classical elements + GMST. The conjunction identifies
optimal-launch candidates satisfying both eclipse- and lighting-time
constraints, e.g. dawn-dusk SSO, GPS-class eclipse-yaw-degraded regimes,
or Molniya apogee targeting. The two pieces of machinery are now
self-contained and importable; the conjunction experiment is a natural
candidate for Exp 015.

## Reproducibility

Commands (this session used the lab `.venv` directly because `uv` was off
PATH; the documented `uv sync && uv run …` commands work identically on
machines with uv installed):

```bash
# Whole-repo test suite
.venv/Scripts/python.exe -m pytest -q

# One-time online (acquisition) -- only if the Sun snapshot is missing
.venv/Scripts/python.exe research/orbital-mechanics/experiments/eclipseTiming/fetch_horizons_sun_snapshot.py

# Analysis run (offline, deterministic)
.venv/Scripts/python.exe research/orbital-mechanics/experiments/eclipseTiming/experiment.py
```

Byte-reproduction: two consecutive runs produce figure MD5s that are
identical and `results.results` payloads that differ only in
`meta.timestamp_utc` and `meta.git_commit`. `code_sha256` in the payload
binds results to source; a stale-run test pins the freshness.
