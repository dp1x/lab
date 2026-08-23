# Experiment Card: Ground Tracks — Spherical Earth, Uniform Rotation, Keplerian Orbits

> Status: complete
> Date: 2026-08-21
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/groundtracks/`

## Research Question

What is the deterministic ground track (geocentric latitude φ(t) and longitude λ(t)) of an Earth-orbiting satellite on a **spherical Earth** with **uniform sidereal rotation**, given Keplerian orbital elements (a, e, i, Ω, ω, M0) propagated either by the closed-form Kepler solution or by the verified 3-D Cowell integrator — and do two independent spherical-geometric derivations, geometric invariants, and numerical propagation agree to machine-grade tolerances for generic, polar, equatorial, retrograde, and highly elliptical (Molniya) real orbit parameters?

Sharpened from the roadmap prompt “Path a satellite traces over Earth | Spherical geometry, real orbit params” (`localdocs/roadmap.md:42`). The question is deliberately framed to be falsifiable: every invariant has a closed-form value to beat.

## Background Theory

### Coordinate frames

* **ECI** (Earth-centered inertial, J2000 pseudo-inertial): origin at geocenter, +Z = Earth rotation axis (CIP), +X = vernal equinox γ at epoch, +Y = Z × X. Inertial: Newton’s `r¨ = −μ r/|r|³` is exact in this frame (Curtis §4.2, Vallado §3.2).
* **ECEF** (Earth-centered, Earth-fixed, WGS-84): same origin and +Z, +X′ = Greenwich meridian, +Y′ = Z × X′. Rotates eastward about Z at `ω_E` (Vallado §3.5, Curtis §5.2). For the spherical calculation the frames share origin and Z; only the Z-rotation `R_z(θ_G)` separates them.

This experiment uses no precession, nutation, or polar motion — handled as idealization and contrasted with realism in the limitations.

### Inertial orbit position `r_ECI(t)` from Keplerian elements

Given `(a, e, i, Ω, ω, M0)` at epoch `t0`:

```
n = √(μ/a³)                                          mean motion               (Vallado Eq.2-23)
M(t) = M0 + n·(t−t0)                                  mean anomaly              (Curtis Eq.2.34)
M = E − e·sin E  → Newton solve for E                  Kepler equation           (Murray & Dermott Eq.2.38)
cos ν = (cos E − e)/(1 − e·cos E)
sin ν = √(1−e²)·sin E /(1 − e·cos E), ν = atan2(sinν,cosν)  true anomaly   (Curtis Eq.3.13)
r = a(1 − e·cos E) = p/(1+ e·cos ν), p = a(1−e²)      radius
r_pf = [ r·cos ν,  r·sin ν,  0 ]ᵀ                      perifocal
Q = R_z(Ω)·R_x(i)·R_z(ω)                               3-1-3 rotation            (Curtis Eq.4.47, Vallado Alg.9)
r_ECI(t) = Q·r_pf(t)                                   (1)
v_pf = (μ/h)[−sin ν, e+cos ν, 0], h = √(μp)           perifocal velocity
v_ECI = Q·v_pf                                           for propagation seed
```

The closed form (1) is the truth model against which the Cowell propagation is tested. The Kepler solver is Newton from `E0 = M + e·sin M`, tol 1e-14, cap 100 iterations — identical to the verified solver at `research/orbital-mechanics/experiments/keplerOrbitValidation/experiment.py:133`.

### Earth rotation

```
θ_G(t) = θ_G0 + ω_E·(t−t0)                              (2)   GMST
ω_E = 2π/T_sid ,  T_sid = 86164.0905 s (WGS-84 sidereal day, NOT 86400 solar)
θ_G0 = 0 at t0                                            (Greenwich ≡ γ at epoch)
```

`T_sid = 86164.090530… s`, `ω_E = 7.2921159×10⁻⁵ rad/s` (WGS-84, Vallado Table 3-1). The solar day 86400 s would give a 0.986°/day systematic longitude error (360°/year) and fail the GEO stationary test — an adversarial check in the suite. `θ_G0 = 0` is an idealized epoch alignment (real J2000 GMST ≈ 280.46°) documented as idealization; it makes GEO exactly stationary and removes a −280° lane bias that would otherwise mask the sign/direction test.

### ECI → ECEF

Passive Z rotation (west-positive longitude drift):

```
R(θ) = [[ cosθ  sinθ  0],
        [−sinθ  cosθ  0],
        [ 0     0     1]]   = R_zᵀ(θ) in active convention       (3)
r_ECEF(t) = R(θ_G(t))·r_ECI(t)                           magnitude preserved: |r_ECEF|=|r_ECI|
```

Then `lon_ECEF = atan2(y_ECEF, x_ECEF) = lon_ECI − θ_G` (mod 360°). Inverting the sign (`R_z(+θ)` active) flips the ground-track drift eastward and is caught by the GEO and equatorial west-drift tests.

### Latitude / longitude on the sphere

Spherical Earth radius `R_E = 6378.137 km` (WGS-84; lab nominal 6378.1 km differs by 37 m, 5.8 ppm). Geocentric ≡ geodetic for a sphere:

```
r = |r_ECEF| = √(x²+y²+z²)
φ_gc = arcsin(z / r) = atan2(z, √(x²+y²))   ∈ [−90°,+90°]          (4a)
λ    = atan2(y_ECEF, x_ECEF)                 ∈ (−180°,+180°] wrapped (4b)
h  = r − R_E  (height, not used for sub-satellite point)
```

Use `arcsin`/`atan2`; `arccos(z/r)` loses sign at the equator. At |φ| > 89.999° (pole) `λ = atan2(0,0)` is undefined — returned as 0.0 by convention and masked in longitude-error metrics; plotting inserts a NaN gap at antimeridian.

### Independent spherical-trig form

Let `u = ω + ν` (argument of latitude). Then from `x_ECI+iy_ECI = r·e^{iΩ}(cos u + i·cos i·sin u)`:

```
sin φ = sin i·sin u          → φ = arcsin(sin i·sin u)                     (5a)
lon_ECI = Ω + atan2(cos i·sin u, cos u)                                     (5b)
λ = lon_ECI − θ_G   (wrapped)                                                (5c)
```

Equations (4) via the rotation matrix and (5) via spherical trig are algebraically identical but share no code — the tightest L1 cross-check (machine precision). For `i = 0`, `φ ≡ 0`, `λ = Ω+ω+ν−θ_G` linear in `t` for `e = 0`; for `i = 90°`, `λ` is piecewise constant between poles.

### Periodicity and repeat ground tracks

Kepler period (nodal ≡ Kepler for pure Kepler):

```
T = 2π√(a³/μ)                                                                (6)
Δλ = −ω_E·T   (mod 360°, west-negative) per orbit at fixed latitude        (7)
Repeat after m orbits and n sidereal days: m·T = n·T_sid  (coprime)       (8)
GEO: a_geo = (μ·T_sid²/4π²)^{1/3} = 42164.169 km → T = T_sid → Δλ = −360° → 0 (stationary)
12-hour orbit (26561.76 km) → T = T_sid/2 → 2 orbits = 1 day repeat
```

For LEO 400 km, `T ≈ 5553 s → Δλ ≈ −23.20°`; ISS 420 km `T ≈ 5578 s → Δλ ≈ −23.31°`. Inclined vs equatorial sharing the same `a` share the same `Δλ` (latitude-independent).

### Invariants and symmetries

* `max|φ| = min(i, π−i)` for retrograde `i > 90°` (e.g., `i = 98° → max 82°`). Polar `i = 90°` reaches ±90°.
* Equatorial `i = 0° or 180°` → `φ ≡ 0`, `Ω` degenerate; handled by setting `Ω = 0` at the equator.
* GEO equatorial circular at `a_geo, i = 0` → `φ ≡ 0`, `λ ≡ Ω−θ_G0` constant; inclined GEO (`i = 5°`) traces a figure-8 lat ±5°.
* `|r_ECEF| = |r_ECI|` exactly (rotation preserves magnitude).
* `λ(t+T) = λ(t) − ω_E·T` wrapped; the unwrapped `λ(t)` satisfies `λ(t+T) = λ(t) + 360° − ω_E·T` for prograde.

## References

* H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier, 2021 — Ch.2 Kepler equation, Ch.3 ν↔E, Ch.4 elements & Q, §5.2 GMST; the 300 km LEO→GEO worked example is reused for constants context.
* D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed., Microcosm, 2013 — §2.2 frames, Eq.3-10 geodetic, Eq.3-34 ECI↔ECEF, §3.5 GMST/ω_E, Table 3-1 ω_E, §9.8 repeat ground tracks.
* R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*, Dover, 1971 — Ch.1-3 two-body, Kepler solution.
* C. D. Murray & S. F. Dermott, *Solar System Dynamics*, Cambridge UP, 1999 — Ch.2 Kepler solution.
* NIMA, *Department of Defense World Geodetic System 1984 (WGS-84)*, TR8350.2 — `R_E = 6378.137 km`, `ω_E = 7.2921159×10⁻⁵ rad/s`.
* IAU 2015 Resolution B3 (Mamajek et al., arXiv:1510.07674) — nominal `GM_E = 398600.4418 km³/s²`; JPL DE440 planet-only `GM_E = 398600.435507 km³/s²` differs by 1.5×10⁻⁸ relative (documented).
* JPL Publication 82-43, *Spacecraft Trajectory Design and B-Plane targeting* (B-plane definition referenced for rotation conventions, not used directly).
* Exp 002 Kepler validation (`research/orbital-mechanics/experiments/keplerOrbitValidation/`), Exp 004 Hohmann, Exp 005 bi-elliptic, Exp 006 plane-change (3-D Cowell RK4), Exp 007 gravity assist — reused machinery and constant provenance.

## Assumptions

* Spherical Earth with `R_E = 6378.137 km` (WGS-84) — **idealization** (oblate `f = 1/298.257` would require geodetic latitude; error ≤0.2° at mid-latitudes, handled as limitation; Exp 009 adds J2).
* Uniform Earth rotation at `ω_E = 7.2921159×10⁻⁵ rad/s`, `θ_G0 = 0` at epoch, no precession/nutation/polar motion — **idealization** (documented; real GMST at J2000 is 280.46°).
* Two-body Kepler propagation (`μ = GM_E` only, no J2, drag, lunisolar, SRP) — **idealization** (standard; Exp 009 is J2 precession, Exp 010 drag).
* Impulsive Kepler elements fixed at epoch; no maneuver, no element osculation — **idealization**.
* Scale-free `μ = 1` normalization not used; physical units `km, km³/s², s, rad` throughout — **verified** (km vs m bug proactively tested).
* Longitude wrapping `(-180°,180°]` with NaN-gap insertion at antimeridian — **convention, documented** (alternative 0–360° equally valid).
* Real orbit parameters are screening conventions (mean TLE-like `a, e, i`), not physical laws — **convention** (cf. `r_p ≥ 1.02 R_eq` in Exp 007).

## Methodology

1. **Closed-form core:** implement `solve_kepler`, `true_anomaly_from_E`, `rotation_matrix_313`, `eci_to_ecef`, `ecef_to_latlon`, `spherical_trig_latlon` as pure functions. Provide both matrix (4) and trig (5) paths for the same physics.
2. **Orbit anchors:** 8 representative Kepler orbits with real parameters — ISS (420 km, 51.6°, e=0.0003), Equatorial LEO (400 km, 0°), Polar LEO (500 km, 90°), SSO-like (600 km, 98°), GEO (42164.17 km, 0°, exactly tuned so `T = T_sid`), inclined GEO (5°), Molniya (a=26560 km, e=0.74, 63.4°, ω=270°), Retrograde LEO (180°). Each carries precomputed `T = 2π√(a³/μ)` and `Δλ = −ω_E·T` (table in Results).
3. **Validation layers:**
   * L1 analytic cross-check: trig (5) vs matrix (4) on identical `r_ECI` over 2 orbits (720/orbit) — expect `<1e-10 deg`.
   * L2 invariants: `max|φ| = min(i,180−i)` to sampling limit, `|r_ECEF|=|r_ECI|` to `1e-14`, equatorial `φ≡0`, polar symmetry, GEO stationary `Δλ=0`, 12-hour 2-orbit repeat, antimeridian step `<5°`, retrograde vs prograde sign, `ω_E` sidereal vs solar distinction.
   * L3 numerical propagation: seed `r0,v0` from `coe_to_rv_eci` at `M0`, integrate raw `r¨ = −μr/r³` via verified Exp 006 3-D Cowell RK4 (`propagate_3d_rk4`, fixed-step, deterministic `dt = T/512` circular, `T/2048` for Molniya), project both analytic and propagated `r_ECI(t)` through the same `eci_to_ecef→lat/lon` chain; measure `max|Δφ|,|Δλ|`.
   * L4 convergence: halve step (128→256→512→1024) on ISS-circular 1-orbit; expect RK4 order 4 (`error ∝ h⁴`, ratio ~16).
   * L5 real anchors: publish `T` and `Δλ` for each anchor; GEO `T` vs `T_sid` to `8.4e-16` relative.
   * L6 pathological: sweep `i ∈ {0,0.01,30,51.6,63.4,89.9,90,90.1,98,120,179.9,180}°` × `e ∈ {0,1e-12,0.3,0.6,0.74,0.8}` — all finite, no NaN/Inf, antimeridian gap handling, pole guard.
4. **Determinism:** pure float64, no RNG, `matplotlib.use("Agg")`, fixed `t` grids, `importlib` single-hop reuse of `propagate_3d_rk4`; two runs byte-identical apart from timestamp.
5. **Figures:** equirectangular ground-track map (3 orbits each), lat & lon vs time (2 orbits), `Δλ` vs altitude, RK4 convergence log-log — all generated deterministically from `results.json` via `make_figures`.
6. **Sweep vs summary:** time series would be large at high cadence; this experiment stores the analytic/invariant summary in `results.json` (compact, committed) and samples figures at 720/orbit (per Exp 002 `STEPS_PER_ORBIT = 512` convention extended for eccentric sampling `(1−e)^{−3/2}`).

## Implementation

* Script: `experiment.py`
* Language/runtime: Python 3.12, numpy 2.5.1, matplotlib 3.11.1, mpmath 1.4.1 (mpmath not required at runtime; used only in prior exps)
* Runtime: `$REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/groundtracks/experiment.py`  (~15 s single core)
* Determinism: pure float64, no RNG, fixed grids, `Agg`, single-hop `importlib` reuse of `planeChangeManeuvers/experiment.py:516` `propagate_3d_rk4` (identical vetting as Exp 007 reuse of same function)
* Reuse: `src/lab_utils/results.py:53` `save_json_result` (provenance envelope), `src/lab_utils/metrics.py:30` `convergence_rate`, Kepler solver logic from `keplerOrbitValidation/experiment.py:133`, 3-D Cowell RK4 from `planeChangeManeuvers` — no scaffolding rebuilt.
* Dependencies: numpy, matplotlib (both already in `uv.lock` from Exp 001–007); no `cartopy`/`poliastro`/`astropy` (would add >200 MB, non-determinism, scope creep — see Exp 009/013).

## Validation Method

* Unit tests: `tests/test_groundtracks.py` (pytest, 31 tests). L1/L2 parametrized over 8 orbits; L3 propagation vs analytic on 5 cases; L4 convergence order; L6 pathological sweep; determinism and unit tests for `ω_E`, `Q` orthonormality, independent bisection Kepler oracle, degree/radian firewall, km/m guard, sign/wrap checks. Tests derive expected from theory or independent math, never by calling the same helper to produce expected (adversarial).
* Analytic identities: trig vs matrix `<1e-10 deg`, `|r_ECEF|=|r_ECI| <1e-14` relative, `max|φ|−i` to sampling limit, GEO stationary `<1e-09 deg`, 12-hour repeat wrapped `<1e-09 deg`.
* Independent integration: `propagate_3d_rk4` (Exp 006, ≤1e-11 h-drift there) vs analytic Kepler; `max|Δφ|,|Δλ|` reported per case (ISS 5 orbits @512 → `9.9e-07 deg` lat, `1.97e-06 deg` lon; Molniya 3 orbits @2048 → `1.87e-04 deg` lat, `6.49e-04 deg` lon — larger due to periapsis resolution law).
* Published anchors: `T` and `Δλ` table compared to closed-form `2π√(a³/μ)` and `−ω_E·T`.

## Results

See `results/results.json` and figures. Headline numbers (machine-readable in `results.json:headline`):

* **Dual-algebra agreement (L1):** `max|Δφ| = 2.27e-13 deg`, `max|Δλ| = 1.14e-13 deg`, great-circle `2.11e-08 rad` (≈4 arcsec — floating `sin/cos` noise floor, not physics) over 8 anchors × 1440 points.
* **Invariants (L2):** `max|φ|` vs `i` error ≤ `1.30e-05 deg` (ISS, `e=0.0003` sampling limit; circular polar/equatorial exactly 0); `Δλ` vs `−ω_E·T` wrapped error exactly `0.0` for all anchors; `|r_ECEF|−|r_ECI|` relative `2.6e-16`; antimeridian max step `0.47°` at 720/orbit (<5° threshold); no NaN/Inf over `12×6` pathological grid.
* **Propagation vs analytic (L3):**
  | case | pts/orbit | orbits | `max|Δφ|` [deg] | `max|Δλ|` [deg] | max `|Δr|/r` |
  |------|-----------|--------|------------------|------------------|---------------|
  | ISS (51.6°, 420 km) | 512 | 5 | 9.90e-07 | 1.97e-06 | 9.28e-10 |
  | Polar (90°, 500 km) | 512 | 5 | 1.26e-06 | 5.68e-14 | 9.28e-10 |
  | Equatorial (0°,400km)|512|5|0.0|1.26e-06|9.28e-10|
  | GEO (0°,42164km)     |512|5|0.0|1.26e-06|9.28e-10|
  | Molniya (63.4°,0.74)|2048|3|1.87e-04|6.49e-04|1.46e-06|
  Molniya larger error is the `(1−e)^{−3/2}` periapsis-resolution law (Exp 002) — not a bug; steps/orbit must scale as `720/(1−e)^{1.5}` (≈5430 at `e=0.74`) for equal accuracy.
* **Convergence (L4):** ISS-like circular 1-orbit, `max|Δφ|,|Δλ|` at 128→256→512→1024 steps: `8.63e-05 → 5.04e-06 → 3.04e-07 → 1.87e-08 deg`; measured order per interval `4.10, 4.05, 4.03`, mean `4.06` (theory 4), confirming RK4 order-4.
* **Repeat (L2):** GEO `T` vs `T_sid` relative `8.44e-16`; GEO 5-orbit closure `1.28e-12 deg` (wrapped); 12-hour 2-orbit closure `4.55e-13 deg` wrapped; GEO wrapped longitude variation `1.48e-12 deg` (stationary to machine).
* **Anchors (real params, idealized Kepler):**
  | orbit | `a` [km] | `e` | `i` [°] | `T` [min] | `Δλ = −ω_E·T` [°] | `max|φ|` [°] |
  |-------|----------|-----|---------|-----------|-------------------|--------------|
  | ISS | 6798.14 | 0.0003 | 51.6 | 92.97 | −23.31 | 51.6 |
  | Equatorial LEO | 6778.14 | 0 | 0 | 92.56 | −23.20 | 0 |
  | Polar LEO | 6878.14 | 0 | 90 | 94.62 | −23.72 | 90 |
  | SSO | 6978.14 | 0 | 98 | 96.69 | −24.24 | 82 (180−98) |
  | GEO | 42164.17 | 0 | 0 | 1436.07 | −360.00 → 0 (stationary) | 0 |
  | Molniya | 26560.0 | 0.74 | 63.4 | 717.96 | −179.98 | 63.4 |
  ISS `T = 5578.22 s` → `15.54 rev/sidereal day` (solar `86400/T = 15.49 rev/day`); canonical ISS TLE mean `92.68 min` is within `0.3%` of the spherical-Kepler value (J2 not modeled — see Limitations).

Figures (`results/figures/`):

* `ground_tracks_map.png` — equirectangular map with 3-orbit tracks for Equatorial, ISS, Polar, SSO, GEO (point at lon 0), and Molniya (1 orbit) with NaN-gap antimeridian splitting and `max|φ|` dashed bounds.
* `lat_lon_vs_time.png` — latitude (with `±i` bounds) and longitude wrapped/unwrapped vs time for ISS, Polar, Equatorial over 2 orbits.
* `delta_lon_vs_altitude.png` — analytic `Δλ = −ω_E·2π√(a³/μ)` vs altitude (200–40000 km), GEO zero crossing at 35786 km highlighted, ISS/Polar/GEO anchors overplotted.
* `rk4_convergence.png` — log-log `max|Δφ|,|Δλ|` vs steps/orbit with order-4 reference.

## Limitations

* Spherical Earth only: `f = 0`, `φ_gc = φ_gd`; real geodetic latitude differs by ≤0.2° at mid-latitudes (WGS-84 `f = 1/298.257`). Exp 009 adds J2 and oblate corrections.
* Uniform rotation: no precession, nutation, UT1−UTC, or Chandler wobble; `θ_G0 = 0` idealized (real J2000 GMST 280.46° would bias absolute `λ` by that constant if not documented).
* Two-body Kepler only: no J2 nodal regression, drag, lunisolar, SRP. GEO and SSO repeat conditions are Kepler-only; real J2 causes SSO to precess at `−2.06°/day` and shifts LEO `T` by ~0.1%.
* Cowell RK4 is not symplectic; energy drift `~1e-09` over 5 orbits at 512/orbit (central-force, not round-off exact) — adequate for ground-track geometry but not for year-long ephemeris (Exp 013).
* Propagation dt is uniform in time, not adaptive on anomaly; high `e` (Molniya 0.74) needs `(1−e)^{−3/2}` step scaling or dense sampling to keep periapsis error bounded — documented and tested, but the default 512/orbit is insufficient for `e > 0.6`.
* `r < R_E` (re-entry/retro) still defines lat/lon but height negative; no atmosphere model (Exp 010).
* Frame is J2000 pseudo-inertial with GMST0=0; linking to a real TLE epoch requires adding the true GMST offset and verifying against `skyfield`/`poliastro` — out of scope for this spherical demonstration.

## Future Improvements

* J2 secular rates (Exp 009): add nodal precession `Ω_dot` and adjust `T` vs `T_nodal`, repeat condition becomes `m(T − ΔT_J2) = n·T_sid`; validate SSO `i ≈ 98°` repeat 233/16.
* Oblate latitude (Exp 009): add `φ_gd = atan2(z + e_e²·N·sinφ_gd, p)` iteration (Vallado Eq.3-10) and publish `φ_gd − φ_gc`.
* Ephemeris-level reconstruction (Exp 013): ingest a historical TLE epoch, propagate with SDP4/SGP4 via `skyfield` and compare lat/lon RMS vs this Kepler sphere.
* Drag analysis (Exp 010): add `B*` and integrate `a_dot` for decay/re-entry timeline.
* Time-accurate GMST: replace uniform `θ_G0 + ω_E·t` with IAU 2006/2000A GMST/GST plus `UT1`.
* Cartographic polish: add `cartopy` coastlines only if/when the 200 MB dependency is justified (Exp 013), kept out here by design.

---

### Reproducibility Notes

* `uv.lock` pins exact dependency versions.
* Command to reproduce: `$REPO_ROOT\.venv\Scripts\python.exe -m pytest && $REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/groundtracks/experiment.py`
  (generic `uv sync && uv run pytest && uv run python research/orbital-mechanics/experiments/groundtracks/experiment.py` when `uv` is on PATH)
* Figures regenerate deterministically via `make_figures` from `results.json` (Agg, fixed `t` grids, no RNG).
* Provenance: `results/results.json:meta` stores `name`, `description`, `timestamp_utc`, `git_commit`, `python_version` via `src/lab_utils/results.py:53` `save_json_result`.

