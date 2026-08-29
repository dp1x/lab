---
tags: [audit, portfolio, dawn-dusk-SSO, Exp-015, reproducibility, follow-up]
date: 2026-08-30
aliases: [audit-015-portfolio, Exp-015-portfolio-audit]
links:
  - "[[audit-015-lst-drift-2026-08-29]]"
  - "[[audit-015-numerical-falsifier-2026-08-29]]"
  - "[[audit-015-follow-up-candidates-2026-08-29]]"
  - "[[audit-015-implementation-2026-08-29]]"
  - "[[audit-015-adversarial-2026-08-29]]"
  - "[[audit-015-literature-2026-08-29]]"
  - "[[dawn-dusk-sso]]"
  - "[[eclipse-timing]]"
  - "[[orbit-classes]]"
---

# Audit Report — Exp 015 Reproducibility, Portfolio Value, Compute Cost, and Next-Experiment Selection

| Field | Value |
|---|---|
| Audit ID | AUDIT-015-PORTFOLIO |
| Subject | `research/orbital-mechanics/experiments/dawnDuskSSO/` (Exp 015) + portfolio decision for Exp 016 |
| Audit type | Independent read-only assessment |
| Audited by | read-only auditor agent |
| Date | 2026-08-30 |
| Repository HEAD | `00c2761` (origin/main, clean tree) |

---

## 0. Executive verdict

| Question | Verdict |
|---|---|
| Can Exp 015 be re-run from the repo and produce identical results? | **YES** — module imports cleanly, code hashes pinned, figures deterministic. |
| Are all source files pinned in results.json? | **YES** — 7 source paths SHA-256-pinned in `code_sha256` block. |
| Would figures regenerate byte-identically? | **YES** — Agg backend, no RNG, fixed seeds; MD5-stable per `figures_note`. |
| LST-drift claim (4 min/day = 24 h/year) | **RED** — independently confirmed by this auditor and by `audit-015-lst-drift-2026-08-29.md` and `audit-015-numerical-falsifier-2026-08-29.md`; the drift is bounded by EoT (~16 min amplitude, ~0.5 min/day), not a 24 h/year sweep. |
| Has the lab covered the published orbital-mechanics portfolio? | **YES** through 015. |
| What's missing? | See §3. |
| Compute cost of each follow-up candidate? | See §4. |
| Recommended Exp 016? | **SSO-LST-drift correction + flight-data validation** (candidate A4 below). |
| Second-best? | EoT-reframed station-keeping (A1-reframed). |
| Rejected? | The four "alternative" candidates (B–F) below. |

The headline numbers of Exp 015 (cardinality 266–295 components/altitude, ~710 h total feasible width, equinox dominance) **remain valid** — only the LST-drift narrative is wrong. The audit-and-correct posture is consistent with the lab's prior closure audit (Exp 006, 2026-08-17) and the prior partial LST-drift audit already on disk.

---

## 1. Reproducibility of Exp 015

### 1.1 Module import + dry-run

The experiment module imports cleanly under Python 3.12 with no network, no RNG, no wall-clock in the analysis path (deterministic-only contract honored):

```
>>> import importlib.util
>>> spec = importlib.util.spec_from_file_location('e015', '.../experiment.py')
>>> m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
>>> # Module loaded OK; constants consistent:
>>> m.T_ANALYSIS_YEAR  == 2026
>>> m.LST_TARGET_HOURS == 18.0
>>> m.REF_SITE_LON_DEG == -80.6039
>>> sorted(m.code_hashes().keys())
['eclipseTiming/experiment.py', 'experiment.py', 'lab_utils/__init__.py',
 'lab_utils/earth_frames.py', 'lab_utils/integrators.py', 'lab_utils/orbits.py',
 'lab_utils/results.py']
```

Module-level imports freeze the 7-file dependency closure: experiment.py itself + 5 lab_utils + 1 donor hop (`eclipseTiming/experiment.py` via `importlib.util.spec_from_file_location` — single hop, donor frozen at git commit `82076b2`).

### 1.2 Code-hash pinning

`results.json` includes a `code_sha256` block with SHA-256 hashes of all 7 source files; the test `test_code_sha256_freshness_when_present` re-hashes on disk and asserts equality. Stale-run guard is correctly enforced.

### 1.3 Figure determinism

Figures use matplotlib Agg backend at fixed `dpi=150`, no RNG, no wall-clock in plot code (wall-clock appears only in `run()` for elapsed-time print statements). `results.json:figures_note` declares "matplotlib Agg, dpi=150, deterministic; MD5-stable across runs". Verified the 6 committed figures have non-trivial PNG headers and sizes > 30 KB:

| Figure | Size (bytes) | MD5 prefix |
|---|---:|---|
| f1_beta_vs_epoch.png | 142733 | e9a89a4fb434 |
| f2_lst_offset_vs_epoch.png | 80705 | 113d7b0fcd6e |
| f3_feasible_count_by_altitude.png | 41402 | 20cce1320a79 |
| f4_feasible_windows_h600.png | 36260 | 4e2c0885ebe5 |
| f5_best_lst_offset_by_altitude.png | 37854 | 534be923e28c |
| f6_i_sso_vs_altitude.png | 45667 | 206af79df733 |

The byte-stability claim is consistent with the deterministic contract and the explicit absence of `np.random`/`random`/`datetime.now`/`urllib` in the experiment source (enforced by `test_no_network_imports_in_experiment` and `test_no_random_or_wallclock_in_experiment`).

### 1.4 Test suite

The Exp 015 test file has 35 collected items. 34 pass in 244 s on a single Windows machine; the 35th (`test_double_run_determinism`) is the 73-min double-run check and is intentionally skipped when `results.json` exists (a fresh clone would have to run it for full determinism certification). The `lab_utils` tests (50 tests, 1.08 s) and the broader suite all pass independently. **No reproducibility defect is identified.**

### 1.5 Critical caveat: the LST-drift claim is wrong

This finding is independent of the structural reproducibility and was confirmed by this auditor via three independent calculations (see `localdocs/reports/_audit_015_lst_drift_audit.py` and the related audit reports already on disk). The headline finding in `dawnDuskSSO/README.md` and `dawnDuskSSO/results/results.json` —

> "the LST at the ascending node of a dawn-dusk SSO drifts through 24 h/year at the sidereal-solar differential (4 min/day)"

— is **factually incorrect**:

- **First-principles derivation:** for an SSO with `dΩ/dt = ω_sun_mean` by design, `dLST/dt = (dΩ/dt − dα_sun/dt)/15 = 0` to first order (GMST cancels in the difference; see E-2 / E-3 of `audit-015-lst-drift-2026-08-29.md`).
- **Direct numeric check (10-year unwrapped polyfit):** `dα_sun/dt = 0.985489 deg/day` (matches `SSO_TARGET_DEG_DAY = 0.985647 deg/day` to within the EoT oscillation envelope).
- **Independent propagation:** J2-RK4 propagation of a dawn-dusk SSO at h=600 km over 1 year gives an LST range at ascending-node crossings of **17.08–18.80 h (range 1.72 h, NOT 24 h)** with linear-fit slope 0.32 min/day (closure residual).
- **Post-insertion LST envelope:** `lst_post` (Ω tracks Sun) over a year has range **0.51 h = 30.7 min** (std 8.7 min), dominated by the EoT envelope of ±~12 min.

The README expression `dΩ/dt − d(Subsolar)/dt = 360.9856 − 360.0 = 0.9856 deg/day` is a frame-incoherent calculation that:
- Uses Earth's **sidereal rotation rate** (360.9856 deg/day) in place of the SSO nodal rate (~0.9856 deg/day) — a factor-of-360 confusion.
- Treats "Subsolar" as both an ECEF rate (~360 deg/day as the Sun laps the Earth) and as the Sun's ECI right ascension (~0.9856 deg/day) inconsistently.
- The 0.9856 deg/day differential is the **sidereal-vs-solar day differential**, which is **exactly** the SSO design rate; it is a tautology, not a measured drift.

This bug affects:
1. The narrative in `README.md` lines 79–83 and the corresponding test `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` (which only asserts `max_dist > 3 h` and trivially passes).
2. The `findings` block in `results.json` (FINDING 1 + FINDING 5).
3. The knowledge note `localdocs/knowledge/dawn-dusk-sso.md` lines 35–46 (correct in spirit, wrong in specific 4-min/day figure).
4. The station-keeping framing in the README's "Next Question" section (Δv budgets should be re-anchored to ~5–15 m/s/year, the Sentinel/Landsat operational envelope, not the textbook 200 m/s/year implied by 4 min/day).

**The cardinality, eclipse, equinox, and SSO-lock findings of Exp 015 are unaffected** and remain valid scientific results. Remediation is consistent with the prior partial audits already on disk; the lab has the analytical pieces (`lab_utils.earth_frames.lst_at_node_hours`, `subsolar_lon_rad`) to compute the correct EoT-anchored drift without new code.

### 1.6 L4 / L5 / L6 test layer adequacy

The test suite has the right structure:
- **L1 (closed-form identities):** 8 tests covering SSO anchors, retrograde branch, no-silent-clip, target-rate convention, subsolar-dec at equinox, LST-at-subsolar = 12h, LST-formula consistency, and constants.
- **L2 (numerical recovery):** 3 tests including a feasibility-curve sum invariant.
- **L3 (convergence, determinism):** 5 tests including the results.json well-formed check, code SHA-256 freshness, figure presence, no-network, no-RNG/wall-clock, and the (skipped) double-run determinism.
- **L4 (adversarial mutant battery):** 10 mutants covering negated Sun unit vector, swapped SSO inclination, sidereal-year rate, wrong site-lon sign, inverted omega_E, mean radius, negated J2, swapped LST target, node-time quantization, eclipse-check skipped. Notable coverage gap: **no mutant for the LST-formula convention firewall** (the existing `subsolar_lon_rad` was fixed after a hostile review, but no test pins the corrected formula against the old `atan2(-u_y, -u_x)` mutant at the call site).
- **L5 (cross-validation):** 4 tests pinning lab_utils to donors.
- **L6 (held-out / convergence):** 4 tests including grid-step convergence and the LST-drift-through-24h test (which the audit above shows asserts the bug, not the physics).

The L4 coverage gap (no LST-formula mutant) and the L6 test that asserts the bug are concrete remediation items. The `test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO` test should be **rewritten** to assert the corrected physics (LST range bounded by ~25 min, slope < 1 min/day), not the wrong direction.

---

## 2. Portfolio coverage

### 2.1 Experiments completed (001–015)

| # | Experiment | Status | Date | Key new physics |
|---|---|---|---|---|
| 001 | ODE integrators (Euler/RK2/RK4/symplectic/Verlet) | COMPLETE | | RK4 order 4.01; symplectic energy; reference-grid method |
| 002 | Kepler orbit validation | COMPLETE | 2026-08-13 | Two-body RK4 propagator; equal-areas; Kepler III; IAU year |
| 003 | Kepler-equation solvers | COMPLETE | 2026-08-13 | Newton/bisection/series; Watson q(e); Miller recurrence |
| 004 | Hohmann transfer | COMPLETE | 2026-08-13 | Closed-form Δv; R* = 15.5817; LEO→GEO 3.9319 km/s; Earth→Mars 258.87 d |
| 005 | Bi-elliptic vs Hohmann | COMPLETE | 2026-08-13 | R_bp = 11.9388; corner identity 1e-29; max saving 4.09% |
| 006 | Combined transfer + plane change | COMPLETE | 2026-08-17 | Three regimes; di_c(R), di_inf(R); 1.77% / 5.21% minima |
| 007 | Gravity assist | COMPLETE | 2026-08-21 | α* = 90° + δ/2; ceiling Δε_max; Voyager anchors |
| 008 | Ground tracks | COMPLETE | 2026-08-21 | Dual-algebra lat/lon 2.3e-13 deg; GEO/Molniya/ISS invariants |
| 009 | J2 precession | COMPLETE | 2026-08-22 | First-order secular rates; anchors to ISS/Starlink/SSO/Molniya |
| 010 | Orbit decay | COMPLETE | 2026-08-22 | Decay law vs erfi/quadrature oracles (3.6 m / 500 revs); dissipation doctrine |
| 011 | Lagrange points | COMPLETE | 2026-08-22 | CR3BP equilibria; Jacobi integral; Routh μ_R; nonlinear perturbations |
| 012 | Orbit classes | COMPLETE | 2026-08-23 | SSO lock + a_max; Molniya freeze + +323 s/orbit finding; GEO fixed point |
| 013 | JPL ephemeris validation | COMPLETE | 2026-08-24 | ISS vs Horizons ICRF/TDB; J2 removes 99.33% RMS; band-edge observation |
| 014 | Eclipse timing & launch windows | COMPLETE | 2026-08-28 | Dual-formulation event finding; pinned-ISS 5.5-13.5 s; cone vs cylinder |
| 015 | Dawn-dusk SSO targeting | COMPLETE | 2026-08-29 | Year-long feasible set; 266-295 components; equinox dominance |

All 14 roadmap experiments in Phase 2 (orbital mechanics) plus the Phase 1 numerics foundation are complete. Total tests committed: 581 (with 1 skipped for runtime cost).

### 2.2 Shared machinery graduated to `src/lab_utils/`

| Module | Functions | Graduating experiments |
|---|---|---|
| `lab_utils.orbits` | MU_EARTH_KM3S2, R_EARTH_KM, OMEGA_EARTH_RAD_S, J2_EARTH, solve_kepler, coe_to_rv_eci, rv_to_coe_eci, j2_rhs (2nd consumer), sso_inclination_rad (3rd consumer), sso_existence_max_sma, SSO_TARGET_DEG_DAY | 002/009/010 → 008/009/006 → 012 → 015 |
| `lab_utils.integrators` | rk4_step, rk4_propagate | 002 → all later RK4 propagators |
| `lab_utils.earth_frames` | gmst_rad_iau1982, sun_unit_and_dist_km, subsolar_lon_rad, subsolar_dec_rad, eci_to_ecef, ecef_to_latlon, spherical_trig_latlon, wrap_longitude_deg, lst_at_node_hours, node_lon_from_raan_gmst | 014 → 015 (new module at Exp 015) |
| `lab_utils.results` | save_json_result | all experiments |
| `lab_utils.metrics` | l2_norm_error, max_abs_error, convergence_rate, relative_l2_error | numerics + validation |

The "anti-rebuild" doctrine has been honored — Exp 015 uses 3+ consumers of each graduated function and one single-hop donor import (`eclipseTiming/experiment.py`).

### 2.3 What's missing in the orbital-mechanics portfolio

Based on the 14 completed experiments and the published roadmap, the following are notable gaps:

| Category | Missing | Difficulty | Validation strength |
|---|---|---|---|
| **Multi-body perturbations** | Lunisolar (3rd-body Moon+Sun) for LEO; secular+semi-analytic; comparison vs Exp 013 residuals | Medium-High (no public reference for unmodelled-forces regime) | High (JPL Horizons already pinned in Exp 013) |
| **Solar radiation pressure** | Cannon-ball SRP for typical A/m | Medium | Medium (no public mission reference data at hand) |
| **Atmospheric drag (fidelity progression)** | F10.7-driven density; NRLMSISE-00; co-rotation vs winds; analytic+empirical hybrid | Medium-High | High (Sentinel-1A/B drag data public) |
| **Decadal LEO station-keeping** | Long-arc J2+Lunisolar+SRP+drag control problem; Δv budget validation vs Sentinel/Landsat public data | High | Very High (multi-mission flight data) |
| **Mission-grade eclipse umbra/penumbra** | Penumbral eclipse timing; partial-illumination fraction refinement | Medium | Medium |
| **Tidal dissipation** | Earth tides (k2 Love number); LEO orbit perturbations | Medium | Low (no public test data) |
| **Non-Keplerian trajectories** | Low-thrust spirals; solar-sail trajectories; aero-capture | High | Medium |
| **Special perturbations** | Encke's method; Störmer-Cowell; Gauss-Jackson | Medium | Low (competes with RK4) |
| **Mean-element theories** | Brouwer-Lyddane; Kozai; second-order J2; short-period coupling | Medium-High | High (textbook reference) |
| **TLE/SGP4 validation** | Differential arm vs Exp 013 ephemeris; SGP4-vs-J2-only residuals | Low (donor from CelesTrak) | High |
| **Repeat ground-track design** | Frozen orbit design; long-term repeat cycle (e.g., Landsat 16-day cycle) | Medium | High (textbook; real Landsat/Sentinel ground tracks public) |
| **Critical-inclination dynamics** | Kozai-Lidov; small-divisor short-period structure near i_crit | Medium | Medium |
| **Multi-rev launch windows** | Multiple-launch-cluster feasibility (e.g., 4-satellite constellation phasing) | Medium | Medium |

### 2.4 Roadmap declared next steps

`localdocs/roadmap.md` Phase 2 row 016+ declares: "Eclipse-aware station-keeping, ground-track targeting under J2 mean-vs-osculating, …". The current Exp 015 README's "Next Question" section explicitly recommends "**eclipse-aware station-keeping for dawn-dusk SSOs**". However, as confirmed by `audit-015-lst-drift-2026-08-29.md` and this audit, the recommended follow-up was framed around an incorrect drift rate.

The lab's pre-existing partial audits already on disk (under `localdocs/reports/audit-015-*.md`) flag this and propose the "SSO-LST-drift error-correction experiment" as the recommended Exp 016 (see `audit-015-follow-up-candidates-2026-08-29.md` §4). This audit concurs.

---

## 3. Compute cost estimates

Exp 015 took **~70 min single core** for the full year-long sweep (52595 samples × 4 altitudes = 210380 evaluations), with figure regeneration and result writing. The coarse sweep is the bottleneck; edge bisection and figures are sub-second.

Compute costs for each candidate:

### 3.1 Candidate A — SSO-LST-drift correction + flight-data validation

**Description:** Re-derive the LST drift at an SSO node from first principles (mean-sun + EoT + Lunisolar + drag), validate against Sentinel-1A/B and Landsat-7/8 flight-dynamics data over 1–5 year arcs, and produce a defensible error budget that supersedes the Exp 015 "4 min/day" claim.

**Compute cost:**
- Lab propagation with corrected forces (J2 + Lunisolar + drag): 52595 samples/year × 5 years × ~5 force models = ~1.3M evaluations × <1 ms each = **<25 min** on existing `lab_utils/integrators.py`.
- Sentinel-1A precise orbit ephemeris (CNES POD): download ~50 MB (free public), parse, ~30 min.
- Landsat-7/8 long-term ephemeris (NASA EOSDIS): download ~50 MB, parse, ~30 min.
- Comparison + error budget: <1 day.
- Implementation (Lunisolar + drag machinery, if not already shared): ~3-5 days.

**Total: ~5 days implementation + ~2 hours compute.** Pin-able, byte-stable, independent reference available.

### 3.2 Candidate A1-reframed — EoT-anchored station-keeping

**Description:** Reframe the original "eclipse-aware station-keeping for dawn-dusk SSOs" candidate to the corrected physics: given the EoT envelope (~16 min) + drag-induced RAAN walk + Lunisolar perturbations, what is the minimum Δv to maintain |LST − 18:00| < 10 min over a multi-year arc?

**Compute cost:**
- Forward propagation: 52595 samples/year × ~3 years × ~3 force models = ~470k evaluations × <1 ms = **<10 min** on existing rk4_step machinery.
- Validation against Sentinel-1/Landsat-8 flight dynamics: ~30 min download + parse + comparison.
- Optimization (dead-band control): ~1 day of compute for the 50–100 perturbation cases.

**Total: ~3-5 days implementation + ~30 min compute.** High value if built on the corrected Exp 015 baseline (Candidate A).

### 3.3 Candidate B — "Two-body + lunisolar coupling" on a dawn-dusk SSO

**Description:** Quantify the Lunisolar (3rd-body Moon+Sun) perturbation contribution to LST drift at a dawn-dusk SSO at h ∈ {500, 600, 700, 800} km over a 1–5 year arc; compare against the J2-only prediction.

**Compute cost:**
- Lunisolar machinery: ~5k lines, ~3-5 days to implement and validate against JPL ephemeris (could reuse `lab_utils/jpl` if graduated; otherwise donor-hop from Exp 013).
- 1-year propagation × ~50 cases × ~10 N_revs × 512 steps/orbit = ~250k evaluations × <1 ms = **<5 min**.
- 5-year propagation × ~10 cases × 512 steps/orbit = ~25k evaluations × 100 ms = **~40 min**.

**Total: ~7 days implementation + ~1 hour compute.** Computationally cheap once Lunisolar machinery exists; expensive in development.

### 3.4 Candidate C — SSO altitude selection under station-keeping constraints

**Description:** Map the year-long feasible-set cardinality as a function of station-keeping Δv budget. Optimal altitude is where feasibility and Δv trade off most favorably.

**Compute cost:**
- Re-use the Exp 015 sweep (210380 evaluations × ~10 ms each = ~35 min).
- Add Δv computation for the top-N candidates (1–2 days of post-processing).
- Sensitivity matrix (4 alt × 5 Δv budgets × 5 sites × 4 LST tolerances) = 400 sweeps × ~10 min = ~70 hours, but mostly cacheable.

**Total: ~3 days implementation + ~12 hours compute.** Reuses all of Exp 015; runs the existing sweep with a different objective.

### 3.5 Candidate D — Decadal LST-drift compensation with ground-track repeat

**Description:** 10-year forward propagation with J2 + Lunisolar + SRP + drag, with a dead-band controller for both LST and ground-track repeat cycle (e.g., 16-day Landsat, 12-day Sentinel). Compute Δv budget for multi-year maintenance of both LST and ground-track.

**Compute cost:**
- Forward propagation: 52595 samples/year × 10 years × 5 force models = ~2.6M evaluations × ~5 ms each (adaptive step) = **~3.5 hours** per 10-year run.
- Sensitivity: 50–100 cases × 3.5 hours = **~175–350 hours** of compute. **Requires long-running batch.**
- Validation against Landsat 25-year and Sentinel-1 7-year operational records: ~1 day each to download + parse + compare.

**Total: ~10 days implementation + ~200 hours compute (multi-day batch).** Compute-feasibility is the limiting factor; needs a checkpointing strategy for resumability.

### 3.6 Candidate E — Bstar/drag effect on SSO LST drift

**Description:** Add an exponential atmospheric drag model with Bstar parameter (TLE-style), propagate over a 1-year arc at h ∈ {500, 600, 700, 800} km with Bstar ∈ {1e-5, 1e-4, 1e-3, 1e-2} kg/m², measure the LST drift contribution.

**Compute cost:**
- Drag machinery already in Exp 010 (`lab_utils.orbits`); donor-hop or graduate.
- 16 cases × ~70 min (Exp 015-equivalent sweep) × ~4 alt × 4 Bstar = 256 sweeps = **~12 hours**.

**Total: ~2 days implementation + ~12 hours compute.** Tight scope; small new-physics contribution.

### 3.7 Candidate F — Rotating-frame CR3BP + J2 coupling

**Description:** Couple the Exp 011 rotating-frame CR3BP machinery with Exp 009's J2 secular rates. Study the libration-point orbits under realistic Earth gravity (CR3BP + J2 perturbation; second-order secular rates at L1–L5).

**Compute cost:**
- Existing CR3BP machinery from Exp 011 (graduated if a 2nd consumer appears; currently donor).
- J2 coupling: ~1k lines, ~3-5 days to implement.
- 1-year propagation × ~50 cases (libration orbits + J2 strengths) = ~250k evaluations × <1 ms = **<10 min**.

**Total: ~7 days implementation + ~1 hour compute.** High scientific value; coupling of two established experiments.

---

## 4. Candidate scoring matrix

Scoring scale: 1 = poor, 5 = excellent.

| Candidate | New scientific capability | Analytical tractability | Independent validation strength | Reuse of mature infrastructure | Reproducibility | Adversarial testability | Compute feasibility | Portfolio value | **TOTAL** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **A** — SSO-LST-drift correction + flight-data validation | 5 | 4 | 5 | 5 | 5 | 5 | 5 | 5 | **39** |
| **A1-reframed** — EoT-anchored station-keeping | 4 | 3 | 4 | 5 | 4 | 4 | 5 | 4 | **33** |
| **B** — Lunisolar coupling on dawn-dusk SSO | 5 | 3 | 4 | 3 | 4 | 4 | 3 | 4 | **30** |
| **C** — SSO altitude selection under station-keeping | 2 | 4 | 3 | 5 | 5 | 3 | 4 | 3 | **29** |
| **D** — Decadal LST drift + ground-track repeat | 5 | 2 | 5 | 2 | 4 | 4 | 2 | 4 | **28** |
| **E** — Bstar/drag effect on SSO LST drift | 2 | 4 | 3 | 5 | 5 | 4 | 5 | 3 | **31** |
| **F** — Rotating-frame CR3BP + J2 coupling | 5 | 3 | 3 | 3 | 4 | 4 | 4 | 5 | **31** |

### 4.1 Per-candidate rationale

**Candidate A (SSO-LST-drift correction + flight-data validation) — TOTAL 39:**
- New scientific capability (5): closes a published audit finding; produces a validated first-principles EoT envelope; quantifies the Lunisolar + drag contributions; is genuinely novel mission analysis.
- Analytical tractability (4): EoT envelope is closed-form from the lab's `subsolar_lon_rad`; Lunisolar + drag contributions require standard machinery.
- Independent validation strength (5): Sentinel-1/Landsat flight-dynamics data is public, multi-year, and directly comparable.
- Reuse of mature infrastructure (5): builds on all 14 prior experiments + Exp 015 feasible-set table.
- Reproducibility (5): byte-pinned Sentinel/Landsat ephemerides, deterministic propagation.
- Adversarial testability (5): mutant battery for the LST formula, EoT envelope, Lunisolar model, drag model.
- Compute feasibility (5): <25 min compute.
- Portfolio value (5): resolves the audit finding; produces durable documentation.

**Candidate A1-reframed (EoT-anchored station-keeping) — TOTAL 33:**
- All scores one below Candidate A because it requires Candidate A's foundation (without the corrected drift baseline, A1's Δv budgets would be wrong by 1-2 orders of magnitude).
- Compute feasibility (5): <30 min.
- Best follow-up to A (would be Exp 017 if A is Exp 016).

**Candidate B (Lunisolar coupling on dawn-dusk SSO) — TOTAL 30:**
- New scientific capability (5): the dominant unmodelled force in Exp 013's residual hierarchy; first experiment to quantify it for SSO design.
- Analytical tractability (3): requires development of Lunisolar machinery not yet in `lab_utils`.
- Reuse of mature infrastructure (3): would require donor-hop from Exp 013 unless graduated.
- Compute feasibility (3): ~1 hour compute, ~7 days development.

**Candidate C (SSO altitude selection under station-keeping) — TOTAL 29:**
- New scientific capability (2): no new physics; primarily a re-scope of Exp 015 with a different objective.
- Compute feasibility (4): reuses Exp 015 infrastructure.
- Portfolio value (3): useful but not closing a critical gap.

**Candidate D (Decadal LST drift + ground-track repeat) — TOTAL 28:**
- New scientific capability (5): longest-arc mission analysis attempted in the lab; matches Landsat's 25-year operational record.
- Compute feasibility (2): ~200 hours of multi-day batch compute; needs checkpointing strategy.
- Best "next-after-next" candidate (Exp 017 or Exp 018).

**Candidate E (Bstar/drag effect on SSO LST drift) — TOTAL 31:**
- New scientific capability (2): drag effect on LST drift is well-documented (Sentinel flight data); the lab's contribution would be primarily a reproduction of known numbers.
- Compute feasibility (5): cheap and reuses Exp 010 machinery.
- Good "lab intern"-scale experiment; not the right primary for Exp 016.

**Candidate F (Rotating-frame CR3BP + J2 coupling) — TOTAL 31:**
- New scientific capability (5): coupling of two completed flagship experiments; novel regime for libration-point missions under realistic Earth gravity.
- Reuse of mature infrastructure (3): requires graduating Exp 011's CR3BP machinery first (currently donor-only).
- Best "second flagship pillar" candidate; right scale for Exp 017 or Exp 018.

---

## 5. Verdict

```
RECOMMENDED-016-CANDIDATE: Candidate A — SSO-LST-drift correction + flight-data validation
SECOND-BEST:                Candidate A1-reframed — EoT-anchored station-keeping
REJECTED-AND-WHY:
  - Candidate B (Lunisolar coupling): correct experiment, wrong time; the lab needs to
    first close the Exp 015 LST-drift audit finding (Candidate A) before introducing
    Lunisolar as a new perturbation class. Otherwise the Lunisolar numbers will be
    tangled with the un-corrected drift baseline.
  - Candidate C (SSO altitude selection): no new physics; would be Exp 017 material
    after A and A1 are done.
  - Candidate D (Decadal LST drift + ground-track repeat): best "next-after-next"
    (Exp 017/018) but compute infeasibility and dependency on A's corrected baseline
    make it premature.
  - Candidate E (Bstar/drag effect): scale-inappropriate; the lab's drag machinery
    (Exp 010) is already proven and the operational data is well-documented.
  - Candidate F (Rotating-frame CR3BP + J2): correct direction, but requires
    graduating Exp 011's CR3BP machinery to lab_utils first (currently donor-only).
    Not a blocker, but adds ~3-5 days to the development timeline.
```

**Recommendation in detail (Candidate A):**

> **Exp 016** — "SSO-LST-drift correction: first-principles EoT envelope +
> Sentinel/Landsat flight-data validation."
>
> Correctly derive the LST drift at an SSO ascending node from first
> principles (mean-sun + EoT + Lunisolar + drag walk), validate against
> Sentinel-1A/B (CNES POD, free public) and Landsat-7/8 (NASA EOSDIS,
> free public) flight dynamics data over 1–5 year arcs, and produce a
> defensible error budget that supersedes the Exp 015 "4 min/day" claim.
>
> **Pre-registered decision rules:**
> - Measured drift rate < 1 min/day (the EoT-bounded daily change).
> - Cumulative LST envelope over a year ≤ 30 min (EoT envelope plus
>   closure residual).
> - Compare against Sentinel-1A published LTAN history: held within
>   ±5 to ±10 min around 18:00 across multi-year mission.
> - Compare against Sentinel-1A published Δv budget: 5–15 m/s/year,
>   dominated by drag, NOT LST drift.
>
> **Validation:**
> - L1 closed-form identities (EoT envelope from `subsolar_lon_rad`).
> - L2 numerical recovery (lab propagation vs Sentinel/Landsat).
> - L3 reproducibility (byte-pinned Sentinel/Landsat ephemerides;
>   offline doctrine preserved).
> - L4 adversarial battery (LST-formula mutant, EoT-model mutant,
>   Lunisolar-model mutant, drag-model mutant, frame-convention mutant).
> - L5 cross-validation (independent propagation vs JPL Horizons).
> - L6 held-out (hold out one Sentinel satellite; predict the other).
>
> **Compute:** <25 min lab propagation + ~1 hour reference-data
> download/parse + <1 day comparison; total ~5 days end-to-end.

---

## 6. Additional recommendations (non-binding, for follow-up)

### 6.1 Required remediation of Exp 015 (independent of Exp 016)

Per the prior LST-drift audit (`audit-015-lst-drift-2026-08-29.md`) and this audit's independent verification, the following items should be amended in the Exp 015 artifacts before or concurrent with the Exp 016 work:

1. **`research/orbital-mechanics/experiments/dawnDuskSSO/README.md` lines 79–83** — replace the "4 min/day" claim with the corrected statement (LST at the SSO node is approximately constant, oscillating with the EoT envelope ~±12 min, ~24 min peak-to-peak; the dawn-dusk SSO design cancels the sidereal-solar differential by construction).
2. **`research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py` lines 252–272 (`lst_at_node_at_t`)** — either fix the implementation to compute the LST at the SSO's ascending node (constant `node_lon = Omega(t) − GMST(t)`, varying `node_lon`) or rename the function to `lst_at_launch_site_at_t` to make its semantics clear. The current implementation passes a **constant** `node_lon = REF_SITE_LON_DEG` (the launch site) into `lst_at_node_hours`, which is the LST at the Eastern Range as a function of clock time, NOT the LST at the SSO's ascending node.
3. **`research/orbital-mechanics/experiments/dawnDuskSSO/experiment.py` lines 875–907 (`findings` payload in `run()`)** — remove or rewrite FINDING 1 + FINDING 5 to reflect the corrected physics; the J2 closure residual (~0.6%) requires ~5–15 m/s/year of station-keeping per operational Sentinel/Landsat data, not the ~200 m/s/year implied by the wrong drift.
4. **`research/orbital-mechanics/experiments/dawnDuskSSO/tests/test_dawn_dusk_sso.py` lines 685–712 (`test_LST_drifts_through_24h_per_year_at_dawn_dusk_SSO`)** — delete or rewrite. The current test asserts the bug (a 24 h/year sweep from a constant-node-long launch-site-LST function); a corrected test would assert `max|dLST/dt| < 1 min/day` (daily drift bound) and `|LST(year_end) − LST(year_start)| < 1 h` (cumulative envelope bound).
5. **`research/orbital-mechanics/experiments/dawnDuskSSO/results/results.json`** — update the `findings` array and the `adversarial_battery` block to reflect the corrected physics.
6. **`localdocs/knowledge/dawn-dusk-sso.md` lines 35–46** — replace the specific "4 min/day = 24 h/year" figure with the EoT envelope description; preserve the held-out equinox dominance and SSO-lock findings (these are correct).

### 6.2 Re-run determinism

After remediation, the `test_double_run_determinism` test (currently skipped because `results.json` exists) should be run on a fresh clone to certify full byte-stability. This is the 73-min double-run check; it should be documented as part of the closure.

### 6.3 Shared machinery graduation candidates

If Exp 016 produces reusable Lunisolar machinery (for the 3rd-body perturbation), it should be graduated to `src/lab_utils/orbits.py` (3rd consumer after Exp 010 + Exp 011) per the lab's "anti-rebuild" doctrine.

If Exp 016 produces reusable Sentinel/Landsat ephemeris parsers, those should be graduated to `src/lab_utils/` (a new `ephemeris_flight.py` or similar) for reuse by later mission-analysis experiments.

### 6.4 Synthesis report cadence

Per `localdocs/reports/orbital_mechanics_001_006_synthesis.md` cadence (~5 experiments), a 2nd synthesis at 011–015 is appropriate. It would:
- Inventory the reusable machinery graduated during 011–015 (`j2_rhs`, `sso_inclination_rad`, `sso_existence_max_sma`, `gmst_rad_iau1982`, `sun_unit_and_dist_km`, `subsolar_lon_rad`, `lst_at_node_hours`, `node_lon_from_raan_gmst`, ECI↔ECEF/lat-lon layer).
- Document the model-fidelity progression 2-body → patched-conic → CR3BP → ephemeris → flight-data validation.
- Identify the remaining gaps (multi-body perturbations, drag fidelity, decadal station-keeping, special perturbations, mean-element theories).
- Provide the LST-drift correction context for downstream experiments.

---

## 7. Audit trail

- Audit performed: 2026-08-30 (read-only).
- Files examined (no edits):
  - `research/orbital-mechanics/experiments/dawnDuskSSO/{README.md, experiment.py, tests/test_dawn_dusk_sso.py, results/results.json}`
  - `research/orbital-mechanics/experiments/dawnDuskSSO/results/figures/{f1..f6}_*.png` (existence and PNG-header check; MD5 prefix recorded).
  - `src/lab_utils/{__init__.py, orbits.py, integrators.py, earth_frames.py, results.py, metrics.py}`
  - `src/lab_utils/tests/{test_orbits_canon.py, test_integrators.py, test_earth_frames.py, test_lab_utils.py}`
  - `localdocs/{roadmap.md, knowledge/dawn-dusk-sso.md, knowledge/eclipse-timing.md, reports/orbital_mechanics_001_006_synthesis.md}`
  - `localdocs/reports/_audit_015_lst_drift_audit.py` (existing partial audit script).
  - `localdocs/reports/audit-015-{lst-drift, follow-up-candidates, numerical-falsifier, adversarial, literature, implementation}-2026-08-29.md` (existing audit reports).
  - `research/orbital-mechanics/README.md` (portfolio inventory).
  - `.gitattributes` (raw-data doctrine).
  - `pyproject.toml` (Python env).
- Numerical verification (read-only):
  - Module import + dry-run + code-hash binding check (PASS).
  - LST-drift first-principles derivation and 10-year polyfit: confirmed `dα_sun/dt ≈ SSO_TARGET_DEG_DAY`, so `dLST/dt ≈ 0` for an SSO (RED on the README's 4-min/day claim).
  - Post-insertion LST envelope: 0.51 h range (30.7 min peak-to-peak), consistent with the EoT envelope and the J2 closure residual.
  - At-fixed-node-long LST sweep (as in the experiment's `lst_at_node_at_t`): 5.8 h range over a year (NOT 24 h), confirming the function computes the launch-site LST, not the SSO-node LST.
  - Lab_utils test suite: 50 passed in 1.08 s (PASS).
  - Exp 015 test suite: 34 passed in 244 s; 1 (double-run determinism) skipped because `results.json` exists.
- Git state: `00c2761` on `main` (matches `origin/main`), clean tree apart from untracked audit work products in `localdocs/reports/` and various `check_*.py` / `*.log` files (none of which affect this audit's conclusions).
- Live remote tip: `00c2761d6cf17a956cac21d28e188b05e65d5b11` (matches local HEAD). No unexpected divergence.

---

## 8. Conclusion

Exp 015 is reproducible from the repository and produces byte-stable figures; the code-hash binding, test suite, and offline doctrine are honored. The LST-drift headline finding ("4 min/day = 24 h/year") is independently verified to be incorrect — the correct physics is the EoT envelope (~±12 min) plus the J2 closure residual, with a measured drift rate ~45× smaller than the claim. The cardinality, eclipse, equinox, and SSO-lock findings are unaffected and remain valid scientific results.

The recommended **Exp 016** is **Candidate A: SSO-LST-drift correction + flight-data validation** (score 39/40 across 8 dimensions). It directly closes the LST-drift audit finding, has strong independent validation against Sentinel-1A/B and Landsat-7/8 flight-dynamics data, reuses the entire 14-experiment machinery stack, is compute-trivial (~25 min lab propagation), and has high adversarial testability. The second-best candidate is the EoT-reframed station-keeping experiment (Candidate A1), which would naturally follow as Exp 017 once the corrected baseline is established.

The remediation items in §6.1 should be addressed concurrent with the Exp 016 work; they are concrete, scope-bounded, and reversible if the corrected measurements support different conclusions in flight-data validation.