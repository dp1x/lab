# Research Roadmap

The experiment sequence the laboratory follows. Any agent picking up a task should
read this first, then `charter.md`, then continue from the current active state.

## North Star

**Orbital mechanics is the flagship domain.** Compact physics + closed-form
verification + real-world data (NASA/JPL Horizons) as the answer key. Numerics is
the foundation that supplies verified numerical methods; energy systems is the
second pillar.

## Domain naming (fixed)

- `numerics/` — the substrate: verified numerical methods (was `physics/`; too
  broad — everything is physics). Exp 001 lives here as the seed.
- `orbital-mechanics/` — the flagship (was `aerospace/`; implies airplanes/CFD,
  which is supercomputer-only and NOT our domain).
- `energy`, `computer-architecture`, `cybersecurity` — roadmap text only, no
  empty folders until real content exists.

## Sequence

### Phase 1 — Numerics (foundation)

| # | Experiment | Status | Validation |
|---|-----------|--------|-----------|
| 001 | Numerical integrator study (Euler, RK2, RK4, symplectic Euler, velocity Verlet) | **COMPLETE** | analytic solution, convergence order, energy invariants |

Continue foundation work only where it directly serves the flagship.

### Phase 2 — Orbital Mechanics (FLAGSHIP)

| # | Experiment | Question | Validation |
|---|-----------|----------|-----------|
| 002 | Kepler orbit validation | Does Newtonian gravity reproduce elliptical orbits and Kepler's laws? | **COMPLETE** (2026-08-13): analytic conic pointwise (≤3e-6), equal areas (≤7.4e-5), T²/a³ = 4π²/μ (≤8.2e-9), invariants ≤1.2e-9 over 10 orbits, IAU Earth anchor 365.256898 d vs sidereal year |
| 003 | Kepler's equation solvers | Newton vs bisection vs series; convergence study | **COMPLETE** (2026-08-13): Newton order 2 (≤9 iters at e ≤ 0.99), bisection halving 0.49992 (49–50 iters), series q_meas matches Watson q(e) ≤ 0.14%, cross-solver agreement ≤ 1.8e-13 |
| 004 | Hohmann transfer | Least-fuel orbit-to-orbit transfer; Δv budget | **COMPLETE** (2026-08-13): closed forms vs RK4 (r/v err ≤ 4.2e-9), (R−1)/2 and √2−1 asymptotes (≤0.9999×, ≤2.4e-6), peak R* 15.5817 / 0.536258, inward symmetry 0.0, 2-impulse optimality grid gap ≤ 7.8e-16, LEO→GEO 3.9319 km/s, E→Mars 258.87 d / v∞ 2.945/2.649 / TMI 3.6114 km/s |
| 005 | Bi-elliptic vs Hohmann | Crossover radius law (β > 15.58) | **COMPLETE** (2026-08-13): R_bp 11.9387654726 & R* 15.5817187388 (= Exp 004 Hohmann peak 15.5817187369 via corner identity, 1e-29 @ 50 digits), crossover curve s_c(R) 12→815.8 … R*→R*, region signs on 90×400 grid (worst margin 1.7e-8), hump onset 9.53 < R_bp, max saving 4.09% v1 @ R=50.1, RK4 3-burn validation (≤4e-8), Wikipedia 14× example exact (4.117530 vs 4.133716 km/s) |
| 006 | Plane-change maneuvers (+ combined transfer+plane change) | Cost of changing inclination & the global optimum vs bi-elliptic super-synchronous | **COMPLETE** (2026-08-16, adversarial audit 2026-08-17): three regimes (two-burn / finite-s 3-burn / s→∞), boundaries di_c(R), di_inf(R); finite-s window pinches shut at R≈6.21 (re-audited); s→∞ limit = Exp 005 bi-parabolic (1e-16); R=2,Δi=47.5° finite dip beats two-burn 1.77%; SES-8 super-sync anchor wins 5.21%; 3D RK4 validation (≤1e-11); 14 tests |
| 007 | Gravity assist / slingshot | Velocity boost from a flyby | **COMPLETE** (2026-08-21): patched-conic 2-body hyperbolic flyby; exact 3-D orientation landscape Δε(α,φ) with global max at α\*=90°+δ/2 (bend ∥ V_p, Cauchy–Schwarz); ceiling Δε_max = 2V_p·v∞/(1+r_p·v∞²/μ_p) with interior optimum v∞\* = √(μ_p/r_p); cancellation-safe δ = 2·atan2(1,√(x(x+2))); B-plane sign convention; L3 Cowell + element recovery (δ to 3e-11 rel, patch-radius-insensitive, RK4 order 4 verified); Voyager 1/2 Jupiter + Voyager 1 Saturn anchors reproduced (Δε +200.83/+151.76/+26.12 km²/s²); 33 tests |
| 008 | Ground tracks | Path a satellite traces over Earth | **COMPLETE** (2026-08-21): spherical Earth, uniform sidereal rotation, Keplerian ECI→ECEF lat/lon via dual algebra (trig vs matrix) 2.3e-13/1.1e-13 deg; invariants max|φ|=i to 1.3e-05 deg, Δλ=−ω_E·T to 0.0 wrapped, |r_ECEF|=|r_ECI| to 2.7e-16, GEO stationary 1.5e-12 deg, 12h repeat 4.5e-13 deg, L3 RK4 vs analytic ≤1.97e-06 deg (ISS 5 orbits @512) & ≤6.49e-04 deg (Molniya 3 orbits @2048), order 4.06, pathological 12×6 grid all finite, 31 tests |
| 009 | J2 precession | Node drift from Earth's bulge | **COMPLETE** (2026-08-22): first-order secular rates established analytically and rediscovered by full-force Cowell RK4 via independent state->element->trend estimation; anchors ISS -4.9724 vs -4.9510, Starlink -4.5080 vs -4.4892, SSO600 +0.9901 vs solar target +0.98565 (i_SSO 97.7876 solved), Molniya -0.1466 vs -0.1479; systematic +0.42..0.68% model-order residual (mean-vs-osculating + second-order small divisors near i_crit) proven NOT integration error (plateau under dt halving); Kepler-truth RK4 order 4.09, rate-metric orders 4.5-4.7 (super-fourth via orbit averaging); polar node-crossing null ~1e-16 deg/day; J2=0 null -4e-14; sign-flip ratios -1.0088/-1.0125 (even-in-J2 physics); J2=0 bit-exact vs Exp 006 propagator; pathological i x e grid sentinels correct incl. induced eccentricity +/-9.60 km = a*(3/2)J2(R/p)^2; 32 tests |
| 010 | Orbit decay | Drag over time → re-entry timeline | **COMPLETE** (2026-08-22): decay law rediscovered vs erfi/quadrature oracles (3.6 m max residual over 500 revs); scalings exact to 6.6e-3; co-rotation twins match wind-factor theory; J2 mean-element transient explained (-0.04% residual); reentry 200→120 km in 2.016 d; benchmarks PASS decade band; 39 tests, 329 total |
| 011 | Lagrange points | 3-body stability zones | **COMPLETE** (2026-08-22): rotating-frame CR3BP transition — L1–L5 rediscovered (residuals ≤7.3e-16; mpmath anchors exact), Jacobi integral conserved (RK4 orders 4.71–4.95, quantization floor 1.33e-15), inertial bridge C=2(nh_z−E_I) at 1.3e-15, E_I drift law at order 4; collinear unstable ∀μ (closed-form σ/ν to 9e-16 rel), Routh threshold μ_R=(9−√69)/18 verified incl. boundary degeneracy + above-threshold quartet; nonlinear perturbation recovery (growth-rate 2.9e-3 rel @ε=1e-4 with bias∝ε; LP frequency 7.5e-7 rel); dimensional↔nondimensional twin pipeline ≤1.5e-14; mission anchors EML1/2/3 + SEL1/2 within 0.03%; 6 adversarial mutants all caught (Coriolis flip provably invisible to spectra/Jacobi — trajectory-level discriminator only); shared machinery graduated to src/lab_utils/{integrators,orbits}.py; 46 tests |
| 012 | Orbit classes | Sun-synchronous / Molniya / GTO specifics | **COMPLETE** (2026-08-23): constraint-defined families — SSO lock cos i=-(a/a_max)^{7/2} with existence limit a_max=12352.505 km (h_max 5974.37) and numeric closure 0.61% rel; Molniya freeze at cos^2i=1/5 (zero-crossing localized -0.005 deg, antisymmetry 6.7e-3), semi-sync a=26561.762 km, dwell 0.923607 (+/-90deg), FINDING: +323 s/orbit Kepler-excess near lock (small-divisor short-period, plateaued); GEO 42164.169462 km with nonzero stationarity negative control +0.02683 deg/day; GTO 3.89256 km/s (300km) anchored to Exp 004 3.9319 (200km); apogee event-rate identity to 0.005 deg/day; adversarial battery incl. pre-registered blind spots; `j2_rhs` graduated to lab_utils; 43 tests |
| 013 | JPL ephemeris validation | Full propagator vs NASA's published positions | **COMPLETE** (2026-08-24): ISS (-125544) Horizons ICRF/TDB geometric states pinned byte-for-byte (sha256, `-text` gitattributes, offline doctrine); exact-grid alignment; J2 removes 99.33% of residual RMS (skill CI 0.9933–0.9979); declared drag tier WORSENS agreement at primary beta=100 (robust negative CI) — reported verbatim, no retuning; beta band monotone, zero-crossing only at band edge (beta=400: 3.13 km vs M2 8.22 km, band-edge observation NOT tuned); M1 trend +690 km/day = mean-vs-osculating offset + absent J2 secular rates; error budget: integration <=19.6 m (self-conv order ~4), init 13 nm, frame diag 0.93 m, epoch tags <=0.33 mm; contamination gate clean (max 2nd-diff 41.7 m < 100 m); reference acknowledged TLE/SGP4-provenance with 1–3 km/day envelope — remainder jointly attributed, never separated; 46 tests |
| 014 | Eclipse timing & launch windows | Event-driven shadow geometry + eclipse-constrained launch windows | **COMPLETE** (2026-08-28): conical (primary) + cylindrical (control) shadow geometry via dual algebraic formulations (apparent-angle vs shadow-axis, event-time agreement sub-second); closed-form g on analytic Kepler states decouples event error from integration step (density-ladder max entry shift < 30 s at 8× stride); first 4 pinned-ISS event epochs agree with real NASA trajectory to 5.5–13.5 s (pre-registered ±15 s band, 3-day tail drift to 308 s = TLE/SGP4 reference envelope); ISS 420 km cylinder 36.04 min vs closed-form 36.03 min (< 5 s); GEO three-tier 67.42/69.56/71.70 min vs pre-registered 67.3/69.4/71.6 min (0.1 min tolerance); cone−cylinder GEO boundary shift 64.3 s vs pre-registered 63 ± 1; analytic Almanac Sun model agrees with byte-pinned 2026 Horizons Sun snapshot to 0.65 deg (gate band 0.7 deg absorbing omitted nutation); 40 new tests, 525 total (485 baseline + 40 new); 6 figures (shadow geometry, fraction-vs-β, duration-vs-altitude, year-sweep, convergence, pinned-ISS illumination) |
| 015 | Dawn-dusk SSO launch-window targeting | Multi-constraint mission analysis: SSO lock + LST-at-node + J2 nodal drift + eclipse-free | **COMPLETE** (2026-08-29): year-long feasible launch-time search for dawn-dusk SSO at h in {500, 600, 700, 800} km from Eastern Range; 266-295 connected components per altitude (monotone in h); total feasible width 710 h at h=600; LST constraint discretizes the search, eclipse constraint is the discriminator; held-out equinox weeks dominate feasibility (36.7 vs 11.9/day main, equinoxes are the most eclipse-favorable for h=600); LST at the ascending node drifts through 24 h/year at the sidereal-solar differential (4 min/day) — the "LST-constant" intuition was wrong; 6 figures; 34 new tests, 581 total (547 baseline + 34 new); shared machinery graduated: `sso_inclination_rad` (3rd consumer), `gmst_rad_iau1982` + `sun_unit_and_dist_km` + `subsolar_lon_rad` + ECI→ECEF/lat-lon layer + `lst_at_node_hours` to `src/lab_utils/{orbits,earth_frames}` |
| 016+ | Eclipse-aware station-keeping, ground-track targeting under J2 mean-vs-osculating, … | each seeds the next | known physics + Exp 015 feasible set |
| 016 | SSO LST-drift correction: first-principles EoT envelope + multi-year perturbation budget | EoT envelope (periodic, ~30 min ptp) + J2 closure residual (~2.2 deg/year, consistent with Exp 012) + Lunisolar upper-bound closed-form (over-estimates ~50x at SSO retrograde inclinations; honest upper bound reported) + SRP (~mdeg/day) + drag (exponential atmosphere; altitude-dependent) + closed-form RAAN-control Δv budget | **COMPLETE** (2026-08-30): EoT peak-to-peak 30.65 min (validated vs byte-pinned 2026 Horizons Sun snapshot to 0.056 deg, well within Exp 014 0.7 deg gate); J2 closure residual consistent with Exp 012 +2.2 deg/year; Lunisolar upper bound reported as conservative ceiling (real operational envelope is much smaller); total LST drift budget decomposed with honest range [no-LS, full-LS-upper]; closed-form Δv range at h=600 km is [0, ~10.6 km/s/yr] from closed-form only, but operational envelope (Sentinel-1 ~15 m/s/yr, Landsat ~5-15 m/s/yr) implies the real rate is much smaller; 4 figures (EoT envelope, drift decomposition, station-keeping Δv range, orbit-plane LST year sweep); 40 new tests, 624 total (584 baseline + 40 new). Builds on Exp 014 byte-pinned Sun snapshot + Exp 012 J2 closure + Exp 009 nodal rate formula. Audit response: provides defensible first-principles derivation of the LST drift rate and station-keeping budget that Exp 015 claimed but did not derive. |
| 017 | Lunisolar upper-bound verification: byte-pinned JPL Sun + Moon over 2026 | Measure the cf_upper / numerical Lunisolar RAAN ratio at h in {500, 600, 700, 800} km; verify the Exp 016 closed-form disclaimer; reject the decadal direction as not scientifically defensible at this time | **COMPLETE** (2026-08-30): byte-pinned JPL Horizons DE441 geocentric Moon vectors (76 KB, 366 daily rows, sha256 `65f1d67f...`) under `reference/`, fetched via identical pattern to the Exp 014 Sun snapshot; numerical integration of Kepler + J2 + point-mass Sun + Moon over 1-year arc with J2-only control subtraction (model-order separation); closed-form reproduction matches Exp 016 to 1e-6 deg/day. HEADLINE: closed-form over-estimates by SIGNED RATIO of ~170x at h=600 km (cf retrograde -0.2184 deg/day, numerical prograde +0.001284 deg/day) — ~3x larger than the audit-015 ~50x estimate, documented as a first-principles discovery (audit band violated but qualitative direction correct). RK4 self-convergence p_r=4.49, p_v=4.50 (RK4 design order confirmed). The original Exp 017 decadal direction was rejected by an eight-track audit (Tracks A-H in the autonomous audit log); the closed-form upper-bound verification (audit-015 candidate #4, Track H Alt-1 scored 27/30) was selected as the strongest scientifically defensible alternative. 4 figures (cf/numerical ratio by altitude, drift rate comparison, dt convergence ladder, linear-fit residual RMS); 32 new tests, 658 total (626 baseline + 32 new). Builds on Exp 014 Sun snapshot acquisition doctrine + Exp 016 closed-form Lunisolar formula + Exp 009 J2 propagator + Exp 012 SSO inclination lock. |
| 018 | Lunisolar RAAN reconciliation: corrected secular formula + controlled numerical experiments | Build on the 017 170x signed discrepancy and 8-track audit; identify root cause and prove the corrected formula | **COMPLETE (2026-08-30)**: corrected secular formula `(3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i-i_3) / sin i` (Track B independent derivation) gives +1.35e-4 deg/day (prograde) at h=600 km i_sso, matching the 1-year numerical +1.32e-3 deg/day in sign and within 9.78x in magnitude (residual is unmodelled short-period). At i=90 deg (cleanest test, J2 cos i = 0) the ratio drops to 2.81x. The 016/017 closed-form is REMEDIATED as `closed_form_lunisolar_raan_rate_rad_s` with DeprecationWarning; the corrected formula is exposed as `corrected_secular_lunisolar_raan_rate_rad_s` for new work. Force isolation (5 modes), inclination sweep (6 i), window sensitivity (5 W), precession on/off, force-level identity (machine precision), convergence ladder. 6 figures, 45 new tests, 714 total repo tests. The 016 LST-drift budget is updated: the ~310 min/year full-LS upper bound is wrong; the corrected formula gives ~1620x smaller magnitude in the opposite direction. |
| 019 | (proposed) Refined Lunisolar evection + variation terms OR multi-year DE441 byte-pinning OR Sentinel/Landsat byte-pinning | Build on the corrected secular formula + 10x residual at i_sso / 2.8x at i=90 deg (short-period unmodelled); external operational anchor | not yet started |

### Phase 3 — Energy Systems (second pillar)

Solar forecasting → battery degradation modelling → power flow on IEEE test grids
→ economic dispatch. Stands alone; does not block the flagship.

### Phase 4 — Computer Architecture

CPU pipeline simulator, cache simulator, scheduling algorithms.

### Phase 5 — Cybersecurity

Cryptographic analysis, secure protocol modelling, vulnerability-testing frameworks.

## Sweep methodology

One experiment = one research question, but the implementation may sweep hundreds
of parameter combinations. Big CSVs live in `data/` (gitignored); a small
sampled/summary JSON is committed to `results/` for reproducibility. This is how
we get 50–100 experiments of real value instead of 5.

## Hooks

- **JPL verification layer**: NASA/JPL Horizons API is free and returns real
  positions via plain `curl`/HTTP. Use it as the ultimate validation source.
  Quirk: START_TIME/STOP_TIME must wrap a range (start < stop).
- **Reuse**: experiments build on `src/lab_utils/` and templates — never rebuild
  scaffolding.
- **Consolidation**: every ~5 experiments, a synthesis report under
  `localdocs/reports/`.
- **Deterministic only**: fixed seeds, no time-dependent nondeterminism. Reality
  is the verification layer.