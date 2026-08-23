---
tags: [orbital-mechanics, ground-track, spherical-geometry, kepler, earth-rotation, validation]
date: 2026-08-21
aliases: [ground-tracks, ground-tracks-spherical, sub-satellite-point]
links:
  - "[[kepler-orbit-validation]]"
  - "[[combined-transfer-plane-change]]"
  - "[[ode-integration-basics]]"
---

# Ground Tracks on a Spherical Earth: Dual-Algebra Validation and Invariants

## Summary

The ground track (geocentric latitude φ(t) and longitude λ(t)) for a Keplerian orbit about a spherical Earth with uniform sidereal rotation is exactly

```
sin φ = sin i·sin u,    u = ω + ν
lon_ECI = Ω + atan2(cos i·sin u, cos u),   λ = lon_ECI − θ_G(t),   θ_G = θ_G0 + ω_E·t
φ = arcsin(z_ECEF/|r_ECEF|),   λ = atan2(y_ECEF, x_ECEF)
```

with `r_ECEF = R(θ_G)·r_ECI`, `R = [[cosθ, sinθ, 0],[−sinθ, cosθ,0],[0,0,1]]` and `T = 2π√(a³/μ)`, `Δλ = −ω_E·T` (west-negative) per orbit. Two independent implementations — (i) rotation-matrix `R·r_ECI → φ,λ` and (ii) spherical trig (5a–5c) — agree to `2.3×10⁻¹³ deg` in φ and `1.1×10⁻¹³ deg` in λ (great-circle `2.1×10⁻⁸ rad` noise floor) over 8 real anchors × 1440 points. Geometric invariants prove correctness: `max|φ| = min(i,180−i)` to `1.3×10⁻⁵ deg` (sampling limit; 0 exactly for circular), `Δλ` vs `−ω_E·T` to `0.0` wrapped, `|r_ECEF|=|r_ECI|` to `2.7×10⁻¹⁶`, GEO stationary `1.5×10⁻¹² deg`, and pathological `12×6` grid all finite with antimeridian step `0.47°`. The verified Exp 006 3-D Cowell RK4 reproduces the Kepler ground track to `≤1.97×10⁻⁶ deg` (ISS 5 orbits @512) and `≤6.49×10⁻⁴ deg` (Molniya 3 orbits @2048), with RK4 order `4.06` (128→256→512→1024). The pattern is reusable: any ECI orbit + `R(θ_G)` + spherical projection yields a verifiable ground track without rebuilding propagators.

## Content

### The two algebras

For `p = a(1−e²)`, `r = p/(1+e cos ν)`, `h = √(μp)`:

* **Matrix path (Curtis Eq.4.47, Vallado Eq.3-34):** build perifocal `r_pf = [r cos ν, r sin ν, 0]`, rotate by `Q = R_z(Ω)R_x(i)R_z(ω)` to `r_ECI = Q·r_pf`, then `r_ECEF = R(θ_G)·r_ECI` (passive-Z, west drift). Latitude/longitude via `φ = arcsin(z/r)`, `λ = atan2(y,x)` wrapped `(-180,180]`. Magnitude preserved by construction; `φ` sign from `arcsin` (not `arccos`), `λ` quadrant from `atan2`.

* **Trig path (independent derivation):** `x_ECI + i y_ECI = r·e^{iΩ}(cos u + i cos i sin u)` gives `sin φ = sin i sin u` and `lon_ECI = Ω + atan2(cos i sin u, cos u)` with `u = ω+ν`. Then `λ = lon_ECI − θ_G` wrapped. Shares no matrix code with the first path except `sin/cos` primitives — a non-trivial cross-check that catches sign, order `R_z(Ω)R_x(i)R_z(ω)` vs permuted, degree/radian, `μ` units, and `R_z(+θ)` vs `R_z(−θ)` bugs simultaneously.

Agreement is at the `10⁻¹³ deg` level because both are `O(1)` trig with the same rounding budget; great-circle `arccos(sin φ₁ sin φ₂ + cos φ₁ cos φ₂ cos Δλ)` inflates to `2×10⁻⁸ rad` due to `arccos(1−ε)` sensitivity near 1, a known noise floor — tighter thresholds would be dishonest.

### Earth-rotation convention (the dominant prior-bug family)

* `ω_E = 7.2921159×10⁻⁵ rad/s` (WGS-84, `T_sid = 86164.0905 s`); solar `2π/86400 = 7.2722×10⁻⁵` is `0.986°/day` systematic error, fails GEO stationary by `0.25°/day` and is caught by the `sidereal_vs_solar` test.
* `R(θ) = [[cosθ, sinθ],[−sinθ, cosθ]]` yields `λ_ECEF = lon_ECI − θ_G` (west). The transpose gives east drift and is caught by `gmst_sign_west_drift` (equatorial prograde must increase unwrapped `336°` per LEO orbit, not decrease) and by GEO requiring `0` not `+720°` over 5 orbits.
* `θ_G0 = 0` (Greenwich ≡ γ at epoch) is idealization; real J2000 `280.46°` would bias absolute `λ` by that constant. GEO stationary test fixes the sign; a TLE-anchored run would need the true GMST offset (Exp 013).

### Periodicity and repeat

`T = 2π√(a³/μ)`; `Δλ = −ω_E·T` (wrapped `(-360,0]`). Equatorial LEO 400 km: `5553.62 s → −23.20°`; ISS 420 km: `5578.22 s → −23.31°`; Polar 500 km: `5676.98 s → −23.72°`; Molniya `43077.76 s → −179.98°`; GEO `86164.09 s → −360° → 0` wrapped. Inclination cancels: `Δλ` is `a`-only for Kepler sphere. Repeat `m·T = n·T_sid` gives `m=1,n=1` GEO and `m=2,n=1` 12-hour (26561.76 km) closures to `4.5×10⁻¹³ deg` wrapped.

### Invariants that are really tests

* `max|φ| = min(i,π−i)` for retrograde (98° → 82°); equatorial `φ≡0` exactly; polar reaches `±90°` exactly at `u=90°` (handled by `asin` clipping `z/r∈[−1,1]` and pole guard `hypot(x,y) < 1e-12·r → λ:=0`).
* `λ(t+T) = λ(t) − ω_E·T` wrapped; the unwrapped `λ(t+T) = λ(t) + 360° − ω_E·T` for prograde is the trap that catches degree/radian and `μ` unit bugs.
* Wrapped longitude step `<5°` at 720/orbit prevents the `359°` jump bug (`|Δλ|>180°` split inserts NaN for plotting).
* Rotation matrix `Q` orthonormal `Q·Qᵀ = I`, `det Q = 1` to `1e-14`.

### Numerical propagation (L3)

Seed `r0,v0 = coe_to_rv_eci(a,e,i,Ω,ω,ν0)` with `ν0` from `M0` via `solve_kepler`, integrate raw `r¨ = −μr/r³` with Exp 006 `propagate_3d_rk4` (fixed step, deterministic). Circular LEO at 512/orbit → `≤1.97×10⁻⁶ deg` over 5 orbits; high `e=0.74` Molniya needs `(1−e)^{−3/2}` scaling — 2048/orbit → `≤6.49×10⁻⁴ deg` over 3 orbits (residual scales linearly with orbit count and with periapsis resolution; 512 would be `100%` at `e=0.85` per Exp 002). Convergence `128→256→512→1024` gives orders `4.10,4.05,4.03`, mean `4.06` (theory 4), identical to Exp 002 `4.07` and Exp 007 `16×` halving.

### Real-orbit anchor arithmetic

`μ = 398600.4418 km³/s²`, `R_E = 6378.137 km` → ISS `a = 6798.137 km → T = 92.97 min → 15.54 rev/sidereal day` (solar `15.49`); canonical TLE `92.68 min` is `0.3%` low — J2 not modeled (Exp 009). GEO `a = (μ T_sid²/4π²)^{1/3} = 42164.169 km` matches `R_E+35786 km` by construction; period match `8.4×10⁻¹⁶` relative is the `a_geo` tuning check. These are screening conventions, not physical laws, as in Exp 007 `r_p ≥1.02 R_eq`.

### Reuse and scope

Uses Exp 002 Kepler solver logic (`solve_kepler`, `true_anomaly_from_E`), Exp 006 `propagate_3d_rk4` via single-hop `importlib`, `lab_utils.results.save_json_result` and `lab_utils.metrics.convergence_rate`. No new integrator, no `cartopy`/`poliastro` (200 MB, non-determinism, scope creep — deferred to Exp 013). Figures are `Agg` deterministic from `results.json`.

## Source Experiments

* `research/orbital-mechanics/experiments/groundtracks/` — spherical-Earth ground tracks: Kepler ECI to ECEF lat/lon, dual-algebra `2.3e-13 deg`, invariants, 3-D RK4 `≤6.5e-04 deg`, order-4, pathological grid. Runnable: `$REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/groundtracks/experiment.py` (or `uv run python ...`)
* `research/orbital-mechanics/experiments/keplerOrbitValidation/` (Exp 002) — RK4 foundation, `solve_kepler`, `orbital_elements`, `STEPS_PER_ORBIT ∝ (1−e)^{−3/2}` law, IAU constants.
* `research/orbital-mechanics/experiments/planeChangeManeuvers/` (Exp 006) — 3-D Cowell `propagate_3d_rk4` reused here.
* `research/orbital-mechanics/experiments/gravityAssist/` (Exp 007) — B-plane rotation and `importlib` reuse pattern, cancellation-safe trig.

## Key Takeaways

* A spherical ground track has an exact, compact form: `sin φ = sin i sin u`, `λ = Ω + atan2(cos i sin u, cos u) − θ_G(t)`, with `u = ω+ν`, `ν(E(M(t)))`; the matrix `R·r_ECI` path is identical and the two agree to machine precision — use this as the primary automated regression.
* The hard bugs are conventions, not physics: sidereal vs solar day (`0.986°/day` drift), `R_z(+θ)` vs `R_z(−θ)` sign (east vs west), `Ω+ω+ν−θ_G` vs `Ω+ω+ν+θ_G`, degree/radian, km/m, `μ` value, `θ_G0` epoch, and longitude wrapping. Each has a one-line test that fails if the sign/rate is wrong (GEO stationary, equatorial west-drift, 90° RAAN shift, 12-hour repeat).
* Periodicity `Δλ = −ω_E·T` is `a`-only for a Kepler sphere; `max|φ|` is `i`-only; both are powerful zero-cost oracles that do not require propagation.
* RK4 ground-track error follows the `h⁴` law and the `(1−e)^{−3/2}` periapsis-resolution law — steep `e` demands adaptive steps or dense sampling; quoting a single `512/orbit` tolerance for all `e` is unsound.
* Reusing the verified 3-D Cowell propagator via explicit `importlib` is the correct laboratory practice: it preserves the `≤1e-11` invariant drift proven in Exp 006 and avoids rebuilding scaffolding.
* Document every constant (`μ`, `R_E`, `ω_E`, `T_sid`, `θ_G0`, frame `R`) with source and note any idealization; the “why” of `θ_G0=0` (Greenwich≡γ) vs J2000 `280.46°` belongs in the card, not hidden in code.

## See Also

* [[kepler-orbit-validation]] — Kepler laws, `solve_kepler`, IAU constants, periapsis-resolution law, sector-area `dA/dt = h/2`.
* [[combined-transfer-plane-change]] — `R_z(Ω)R_x(i)R_z(ω)` rotation, 3-D Cowell validation to `1e-11`, non-unimodality lesson.
* [[ode-integration-basics]] — integrator order vs geometric fidelity, why RK4 is not symplectic.
* [Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier, 2021]
* [Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed., Microcosm, 2013]
* [Bate, Mueller, White, *Fundamentals of Astrodynamics*, Dover, 1971]
* [Murray & Dermott, *Solar System Dynamics*, Cambridge UP, 1999]
* [WGS-84, NIMA TR8350.2] and [IAU 2015 Resolution B3, arXiv:1510.07674]
