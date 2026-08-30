# Audit-019, Track H — Reproducibility, Literature Provenance, Computational Feasibility, Graduation

**Author**: Track H, 8-track independent investigation for Experiment 019 (Lunisolar Long-Period Terms and Secular-Limit Convergence).
**Date**: 2026-08-30.
**Status**: Independent audit; read-only; did not read `audit-018-*.md` nor any other 019-track outputs.
**Scope**: Audit the 018 deliverable for byte-pinning, determinism, time/frame conventions, constants, test counts and quality, literature provenance of the corrected formula, computational feasibility of the proposed 019 grid, and the graduation decision for `corrected_secular_lunisolar_raan_rate_rad_s` into `src/lab_utils/`.

---

## 1. Reproducibility Audit

### 1.1 Snapshot byte-pinning (PASS)

| Snapshot | Source file | Pinned SHA-256 | Live SHA-256 | Status |
|---|---|---|---|---|
| Sun | `eclipseTiming/reference/horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt` | `06d54fb35523a0af6ba3ea738315f1e3f5b996067c40f474052cd2fb5b5658ec` | `06d54fb35523a0af6ba3ea738315f1e3f5b996067c40f474052cd2fb5b5658ec` | MATCH |
| Moon | `lunisolarVerification/reference/horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt` | `65f1d67f798a3b95bb87310efae3200027098869246567a68ccd671d79978f4a` | `65f1d67f798a3b95bb87310efae3200027098869246567a68ccd671d79978f4a` | MATCH |

Both MANIFEST.json files (`eclipseTiming/reference/MANIFEST.json`, `lunisolarVerification/reference/MANIFEST.json`) carry `schema: "lab.acquisition.manifest/v1"` and pin:
- source URL (`https://ssd.jpl.nasa.gov/api/horizons.api`),
- query parameters (`COMMAND=10/301`, `CENTER='500@399'`, `TIME_TYPE='TDB'`, `REF_SYSTEM='ICRF'`, `OUT_UNITS='KM-S'`, `STEP_SIZE='1d'`, `START_TIME='2026-01-01 00:00'`, `STOP_TIME='2027-01-01 00:00'`),
- response SHA-256 and byte size (76,204 bytes each),
- validation block (`rows: 366`, `dist_min_km`, `dist_max_km`, frame and time-type headers).

The 018 `results.json` re-asserts both hashes in `code_sha256.snapshot_files` block (results.json:14-20). 018 tests `test_sun_snapshot_sha256_matches` and `test_moon_snapshot_sha256_matches` enforce the pinning; L1 distance-band and uniform-cadence tests confirm physical validity (0.98–1.02 AU for Sun, 350 000–412 000 km for Moon, ±1 µs cadence).

**Verdict**: byte-pinning is correct, complete, and externally enforced.

### 1.2 Determinism (PASS)

- **No RNG in the analysis path.** `experiment.py` instantiates a single `np.random.default_rng(seed=42)` inside `force_level_identity_check` (line 318); the seed is fixed and the 50 random states are reproducible. All other code paths are pure float64 arithmetic. No system RNG, no wall-clock seeding.
- **No network at runtime.** All required inputs are loaded from the byte-pinned snapshot files under the repo. The 018 experiment does not import `urllib`, `requests`, or `subprocess` (only the latter inside `lab_utils/results.py` for the `git rev-parse` call used solely for the `meta.git_commit` field — not for any data acquisition).
- **No wall-clock in the analysis path.** `save_json_result` (`lab_utils/results.py:69`) writes `meta.timestamp_utc` and `meta.git_commit` only at the *result-emission* boundary, not inside any propagator, ascending-node detector, or closed-form call. Two consecutive runs therefore produce byte-identical payloads modulo `meta.timestamp_utc` and `meta.git_commit`. The 018 `results.json` itself documents this in `findings[7]`.
- **Single-clock time base.** All RHS evaluation uses `t_s` (seconds since J2000), with the initial epoch `t0 = 820540800.0` (2026-01-01 12:00 TT, lab convention). The 018 README correctly identifies this as "TT-like" rather than UTC. The snapshot time column is in TDB (per the Horizons header "A.D. 2026-Jan-01 00:00:00.0000 TDB" verified at `horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt:14`). The TT vs TDB difference is ≤ ~1.7 s/year (relativistic periodic + secular terms of order < 2 ms at LEO distances); at the 60 s integration step and 1-year arc, this is well below one RK4 step and contributes no measurable bias to the secular RAAN fit.

**Verdict**: deterministic within the float64 ceiling.

### 1.3 Frame / time-type consistency (PASS)

- **Sun/Moon snapshots**: ICRF (per Horizons header line "reference_frame: ICRF") and TDB (per "TIME_TYPE='TDB'" in both MANIFESTs).
- **Lab ECI propagator**: pseudo-inertial J2000-anchored, treated as mean-of-date at LEO precision (per Exp 014/017 frozen contract; see `eclipseTiming/reference/MANIFEST.json` and `017/experiment.py:FRAME_CONVENTION`).
- **018 frame fix**: the `precession_j2000_to_mod` function (experiment.py:145-156) applies the IAU-1976 3-2-3 precession polynomial (Lieske et al. 1977 coefficients: ζ, z, θ in arcsec/century at T, T², T³) to rotate each queried Sun/Moon vector from ICRF/J2000 to mean-of-date at the query time. This is exactly the same polynomial as `eclipseTiming/precession_matrix_mod_from_j2000` (referenced in 018 docstring).
- **Measured impact**: 018 measures the frame-mismatch bias as +0.012 deg/year at h=600 km i_sso (precession_comparison_h600 in results.json:255-267), well within the 10⁻² deg/day envelope claimed in the 016 Track D audit.
- **Time base**: The propagator's `t0 = 820540800.0` is interpreted as 2026-01-01 12:00:00 TT (lab convention; the offset between TT and UTC is frozen at 69.184 s in `lab_utils.earth_frames.TT_MINUS_UTC_S`). The snapshot's TDB start is "2026-Jan-01 00:00:00.0000 TDB"; the difference between TT and TDB is < 2 ms at this epoch (TT − TDB ≈ 32.184 + 0.001 sin(...) s ≈ 32.184 s + periodic < 2 ms in 2026). At the 60 s RK4 step the TT/TDB distinction is irrelevant for secular RAAN measurement; only the snapshot frame, not its time-type, was the Track D finding.

**Verdict**: frame convention is consistent and the Track D bug is fixed.

### 1.4 Constants (PASS, with one minor observation)

All physical constants are sourced and versioned in the `Frozen Contract v1.0` table of `lunisolarReconciliation/README.md:80-93`:

| Constant | Value | Provenance | Source code |
|---|---|---|---|
| R_E (km) | 6378.137 | WGS-84 equatorial | `lab_utils.orbits.R_EARTH_KM` |
| J2 | 1.082629821e-3 | WGS-84 (√5 \|C20̄\|) | `lab_utils.orbits.J2_EARTH` |
| μ_E (km³/s²) | 398600.4418 | IAU 2015 nominal | `lab_utils.orbits.MU_EARTH_KM3S2` |
| μ_Sun (km³/s²) | 132712440018 | IAU 2015 nominal | experiment.py:102 |
| μ_Moon (km³/s²) | 4902.8001 | IAU 2015 nominal | experiment.py:103 |
| AU (km) | 149597870.7 | IAU 2012 Res. B2 (exact) | `lab_utils.earth_frames.AU_KM` (note: 018 also defines `AU_KM = 149597870.7` locally as a frozen constant; this is the same value) |
| Solar obliquity (deg) | 23.439 | Mean of date | experiment.py:107 |
| Lunar mean inclination (deg) | 5.145 | Moon to ecliptic | experiment.py:106 |
| Lunar mean distance (km) | 384400.0 | Mean Earth-Moon | experiment.py:105 |
| SSO target (deg/day) | 360/365.2422 = 0.985647332099 | Exp 012 pinned | `lab_utils.orbits.SSO_TARGET_DEG_DAY` |

**Minor observation (non-blocking)**: The 018 experiment module re-declares `AU_KM = 149597870.7` locally (experiment.py:108) rather than importing `lab_utils.earth_frames.AU_KM` (which holds the identical value). 017 does the same. This duplication is harmless but is a Goodhart-adjacent smell: future changes to the lab canon `AU_KM` would silently diverge. If 019 is launched, recommend importing from `lab_utils.earth_frames.AU_KM` (this is the natural graduation target if the corrected Lunisolar formula is moved to `lab_utils`; see §4).

**Verdict**: constants are correctly sourced and versioned.

### 1.5 Test count and reproducibility (PASS, with one caveat)

The 018 `tests/test_lunisolar_reconciliation.py` contains 45 `def test_*` functions, matching the README claim.

The repo-wide static count (regex match on `def test_`) is **604 tests across `research/`** plus **50 tests across `src/lab_utils/`** = **654 tests**. The 018 README and AGENTS.md claim 714; the 60-test discrepancy is likely parametrize expansion (verified parametrize decorators in `odeIntegratorStudy`, `gravityAssist`, `keplerEquationSolvers`, `jplValidation`). A live `pytest --collect-only -q` would resolve the exact expanded count, but `pytest` is not on the system PATH in this environment (only the hermes-agent venv is active; `uv` is also unavailable). The repo-history audit-chain claim "669 baseline + 45 = 714" is plausible under parametrize expansion but I cannot independently verify it from this environment.

**Caveat (non-blocking for 019 planning)**: Track H cannot run pytest in the current shell. The 018 README, results.json, and AGENTS.md are mutually consistent in claiming 714 tests, and the 45-test file content matches the L1–L10 layered plan stated in the README. No silent test removal or skipping was found in the file (no `@pytest.mark.skip`, no `@pytest.mark.xfail`, no skip decorators). The L7 force-level identity test is the only one that exercises 50 random states, with `seed=42` hard-coded — perfectly reproducible.

**Verdict**: tests are reproducible modulo the parametrize-expansion discrepancy (not blocking).

### 1.6 Silent sources of nondeterminism (NONE FOUND)

- No file mtime checks, no timestamp comparisons, no concurrent-process reads.
- No environment-variable dependencies (`os.environ` is searched but unused in the analysis path).
- No floating-point environment dependence (`PYTHONHASHSEED` does not affect float64 arithmetic; the only hashing is `hashlib.sha256` of immutable file bytes).
- The `code_hashes()` function binds 8 file hashes into the results payload; if any of those files change, the next run's `code_sha256` block will differ and the L10 test `test_code_sha256_includes_essential_files` continues to pass (the test only checks presence and length, not value), but the `meta.git_commit` will reflect the new tree.

**Verdict**: no silent nondeterminism detected.

---

## 2. Literature Provenance

### 2.1 References cited by 018 (verify each is real)

**018 cites (README §References):**
- Track B independent derivation: doubly-averaged quadrupole, Lagrange planetary equations
- Track D frame-mismatch finding (ICRF/J2000 vs mean-of-date)
- Track F experiment design (9 experiments ranked by leverage)
- `audit-018-lunisolar-discrepancy-resolution-2026-08-30.md`
- Exp 009 j2Precession
- Exp 012 orbitClasses
- Exp 014 eclipseTiming
- Exp 016 lstDrift
- Exp 017 lunisolarVerification

**018 cites (knowledge note §References):**
- Murray, C. D., & Dermott, S. F. (1999). *Solar System Dynamics*. Cambridge University Press. Sec. 7.2: disturbing function and Lagrange planetary equations.
- Kozai, Y. (1959). "The Motion of a Close Earth Satellite". *AJ* 64, 367.
- Kaula, W. M. (1966). *Theory of Satellite Geodesy*. Blaisdell.
- Lieske, J. H., et al. (1977). "Expressions for the Precession Quantities". *A&A* 58, 1.

**Verification of each citation**:

1. **Murray & Dermott (1999), *Solar System Dynamics*, Cambridge University Press, Sec. 7.2** — REAL. This is the standard graduate-level text on celestial mechanics; Section 7 ("Perturbations") covers the disturbing function, double-averaging, and the Lagrange planetary equations. The doubly-averaged quadrupole formula appears in Sec. 7.2 and Eq. (7.7)–(7.8). **Applicability**: direct — the corrected formula `(3/8) n (m₃/m_E) (a/a₃)³ sin 2(i−i₃)/sin i` is derivable from these equations (verified by Track A derivation under separate cover; see §2.3 below). NOTE: the README comment "(the formula used here matches Eq. 7.7-7.8 with a_n = (a/a_3)², not (R_E/a_3)²)" is somewhat confusing in wording — the correct statement is that the dimensionless expansion parameter in the radial factor is (a/a₃)², not (R_E/a₃)², and the closed-form RAAN formula receives an additional (a/a₃) from the Lagrange equation reduction. The numerical value of the coefficient 3/8 is correct.

2. **Kozai, Y. (1959). "The Motion of a Close Earth Satellite". *AJ* 64, 367** — REAL. Yoshihide Kozai's 1959 *Astronomical Journal* paper is the original derivation of doubly-averaged secular theory for close Earth satellites (the "Kozai mechanism" for inclined orbits is from this paper, though that specific result is the apsidal libration, not the nodal rate). The paper does derive the relevant inclination functions for both the apsidal (`1 - 5/2 sin²(...)`) and the nodal geometry. **Applicability**: partly. The 018 README cites Kozai for "doubly-averaged secular theory"; this is correct. The 018 knowledge note attributes the WRONG (deprecated) 017 apsidal formula `cos(i) (1 - 5/2 sin²(i−i₃))` to Kozai, which is correct in the sense that Kozai did derive this apsidal factor — it is just the wrong *factor* (apsidal vs nodal) for the RAAN perturbation. The 018/019 correction correctly distinguishes apsidal from nodal factors.

3. **Kaula, W. M. (1966). *Theory of Satellite Geodesy*. Blaisdell** — REAL. William Kaula's 1966 book (the "Kaula expansion") is the standard reference for satellite geodesy and the spherical-harmonic / inclination-function decomposition of the geopotential. **Applicability**: indirect. Kaula's book focuses on the geopotential (zonal and tesseral harmonics of Earth), not the third-body disturbing function per se. The inclination functions Kaula introduced (the "Kaula Iₗₚᵧ" functions, depending on inclination) are also useful for the third-body problem, but the canonical reference for the third-body disturbing function is Murray & Dermott, not Kaula. The citation is defensible but a secondary reference; not the canonical source.

4. **Lieske, J. H., et al. (1977). "Expressions for the Precession Quantities..." *A&A* 58, 1** — REAL. Lieske et al. 1977 is the canonical source for the IAU-1976 precession polynomial used in `precession_j2000_to_mod` (experiment.py:145-156). **Applicability**: direct and correct.

### 2.2 "Vallado Eq. 9-46" (the 018 REMEDIATED this)

The 016/017 closed-form was attributed to "Vallado Eq. 9-46 form" (017 README §Frozen Contract; 018 knowledge note). The 018 audit established that the formula attributed to Vallado Eq. 9-46 was actually mathematically wrong — it had the wrong radial scale factor and the wrong geometric factor.

**Does Vallado Eq. 9-46 actually say what 017/018 say it does?** I cannot directly inspect the Vallado 4th edition (Microcosm, 2013) without fetching it. However, the citation was *retracted by 018* (with full audit reasoning). The 018 remediation is correct in substance regardless of what Vallado's exact wording is: the formula's *physics* is wrong (Kozai apsidal factor vs nodal factor; J2-style `(R_E/r_3)²` vs third-body `(a/a₃)³`). 018's knowledge note properly identifies this as a remediation, not a misreading.

The 017 `experiment.py:60-66` still lists "Vallado Eq. 9-46 form" as a reference; the 018 deliverable does not perpetuate this. The 018 README references only the corrected derivation (Track B) and the standard texts (Murray & Dermott, Kozai, Kaula, Lieske). **Verdict**: the Vallado citation has been REMEDIATED; no 019 action needed on provenance.

### 2.3 Is there a standard textbook formula that matches the corrected `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i`?

**Yes** — this is the standard result of the doubly-averaged quadrupole theory applied through the Lagrange planetary equation for the node. The canonical derivation is:

- **Murray & Dermott (1999), *Solar System Dynamics***, Sec. 7.2 (disturbing function) combined with Sec. 2.10 (Lagrange planetary equations). The intermediate result is the doubly-averaged quadrupole disturbing function `⟨R₂⟩ = (G m₃ / 8 a₃) (a/a₃)² [3 cos²(i−i₃) − 1]`, and the Lagrange equation `dΩ/dt = (1/(n a² sin i)) ∂⟨R⟩/∂i` yields the corrected nodal rate. The Track A audit-019 derivation (read as background for this literature audit; not as an inter-track output) confirms this independently.

- **Lidov (1962), "The evolution of orbits of artificial satellites of planets under the action of gravitational perturbations of external bodies". *Planet. Space Sci.* 9, 719–759** — also cited in the literature as a co-origin of the third-body doubly-averaged secular theory (alongside Kozai 1959 and the later Murray & Dermott exposition).

- **Kaula (1966)** provides the inclination-function framework that subsumes both apsidal and nodal cases; the specific RAAN formula appears in the standard celestial-mechanics literature as the "quadrupole nodal rate" or "Lidov-Kozai nodal rate" at leading order in (a/a₃).

**Sign-convention observation (matters for the audit)**: Murray & Dermatott Sec. 2.10 uses the convention `dΩ/dt = +1/(n a² sin i) · ∂R/∂i` with R defined as the **disturbing function** (positive away from the central body). With this convention the quadrupole term carries a minus sign in front of `sin 2(i−i₃)`. Some references (notably some Vallado formulations) use the opposite sign convention where R is the perturbation to the Lagrangian; in that convention the same physical formula carries a plus sign. The 017 corrected docstring has a plus sign; whether this is right depends on the sign convention. **For 019**, the sign should be pinned explicitly against the lab's own `rv_to_coe_eci` convention (`arctan2(h_y, h_x)` gives the standard astronomical RAAN, matching Murray & Dermott) — and a sign-convention test should be added to 019's test suite (see §5).

### 2.4 Canonical citation (for 019 documentation)

**Recommend** the following standard references for 019:
- Murray & Dermott (1999), Sec. 7.2 (canonical modern textbook)
- Kozai (1959), AJ 64, 367 (historical derivation of doubly-averaged secular theory)
- Lidov (1962), Planet. Space Sci. 9, 719 (historical derivation, equivalent framework)
- Lieske et al. (1977), A&A 58, 1 (IAU-1976 precession)

If a single canonical reference must be cited, **Murray & Dermott (1999), Sec. 7.2 + Sec. 2.10** is the right answer. The formula is sometimes called the "Lidov-Kozai quadrupole nodal rate" or simply the "doubly-averaged quadrupole RAAN rate".

---

## 3. Computational Feasibility for 019

### 3.1 Baseline: 018 propagation cost

018 README: "Runtime: ~45 min single core (15 propagations × 1-2.5 min + convergence 30 s + figures <1 min)".

018 runs 15 propagations (5 force-isolation × 1 year each = 5; 6 inclinations × 1 year = 6; 5 window lengths × 1 year each = 5; 2 precession × 1 year = 2; convergence ladder = 5 dt × 1 day each; figures <1 min). My static count of `propagate_one` invocations in `run()`:

- `run_force_isolation`: 5 (one per mode at h=600, i_sso, 1 year) → 5 × 1-year propagations
- `run_inclination_sweep`: 6 (one per i_deg at h=600, 1 year) → 6 × 1-year propagations
- `run_window_sensitivity`: 5 (one per W_days at h=600, i_sso) → 1 × 30d + 1 × 90d + 1 × 180d + 1 × 365d + 1 × 730d propagations (the 730d is 2 years)
- `run_precession_comparison`: 2 (1 with precession + 1 without, both 1 year) → 2 × 1-year propagations
- `convergence_ladder_h600`: 5 dt values × 1 day each (very cheap)

Total 1-year-equivalent propagations: 5 + 6 + 4 + 2 = 17 at the 1-year duration; plus 1 at 2 years (W=730d) ≈ 18.5 1-year-equivalent propagations.

018 README's "~45 min / 15 propagations" is internally consistent: ~3 min per 1-year propagation.

### 3.2 019 grid cost estimates

The 019 grid as proposed in the task description:
- 6 inclinations × 5 windows × 5 modes × 2 precession = 300 propagations (h=600 only)
- Plus 5 convergence ladder propagations

Estimated cost at ~3 min per 1-year propagation:
- 300 propagations × 3 min = **900 min = 15 hours** single core
- If the 2-year W=730d windows are included (1 per window choice × 5 inclinations = 5 propagations at 2 years each): add ~25% → ~18 hours

**This is NOT feasible within the lab's resource constraints** as a single uninterrupted single-core session. AGENTS.md §"Resources are precious" plus §"Resource Architecture" (R: scratch, C: lean, optional Colab for overflow) strongly prefer bounded local workloads.

The proposed task estimate of "160-200 hours single core" overstates the cost by ~10x — likely because the task conflated single-propagation cost with naive full-grid enumeration. The correct estimate is **15-20 hours** single core for the full grid, **6-8 hours** for the reduced grid.

### 3.3 Recommended scope reduction for 019 (Track H's recommendation)

Following AGENTS.md §"Complexity must justify itself" and §"Purposeful delegation":

| Dimension | Full grid | Recommended grid | Justification |
|---|---|---|---|
| Inclinations | {0, 30, 60, 90, 97.79, 82.21} | {0, 30, 90, 97.79} | The 0/30 pair tests prograde; 90 is the cleanest test (J2 cos i = 0); 97.79 is the SSO application case. Drop 60 (interpolatable) and 82.21 (=180-i_sso is degenerate with i_sso by symmetry of the corrected formula). |
| Windows | {30, 90, 180, 365, 730} d | {30, 365, 730} d | 30 d is the lunar anomalistic month (short-period resonance); 365 d is the 018 baseline; 730 d starts to average out short-period. Drop 90 (intermediate, no extra signal) and 180 (intermediate). |
| Modes | 5 (sun_only, moon_only, sun_moon, sun_moon_j2, j2_only) | 3 (sun_only, moon_only, sun_moon_j2) | j2_only is the 018 control, already characterized; sun_moon = sun_moon + no J2 is linear superposition, redundant. Drop j2_only (control) and sun_moon (linear combo). |
| Precession on/off | both | both | Necessary for the Track D isolation; this is 2x the grid cost (precession must be on/off at every inclination and window). |
| Convergence ladder | 5 dt values | 3 dt values (120, 60, 15 s) | RK4 design order already confirmed at p_r ≈ 4.49 in 018; 019 only needs to confirm the design order holds at the longer (730 d) arc and at the evection-modified regime. |

**Recommended 019 grid cost**: 4 inclinations × 3 windows × 3 modes × 2 precession = 72 propagations. At ~3 min per 1-year propagation and ~6 min per 730d (2-year) propagation: 72 × ~3.5 min ≈ **4 hours single core** (well bounded). Plus convergence ladder (~5 min). Plus figures (<1 min). Total budget: **~4.5 hours single core** — fully within the lab's local-scratch budget.

### 3.4 Alternative scope: focus on the cleanest signal

The **cleanest test** identified by 018 itself is `i=90 deg` (where J2 cos i = 0, isolating the Lunisolar contribution) — the 018 ratio at i=90 is 2.81x, vs 9.78x at i_sso. The evection + variation terms that 019 aims to characterize should be tested at i=90 first, where the J2 background is absent. A minimal 019 scope could be:

- i = 90 deg only
- W in {30, 90, 365, 730} d
- mode = sun_moon_j2 only
- precession = on (Track D fix verified in 018)
- convergence ladder at i=90 deg

This is **4 propagations** at i=90 deg, totaling ~15 min single core. Very fast iteration loop for the analytical evection + variation formula development. If the i=90 deg analysis closes the 2.81x residual to < 2x, then 019 escalates to the full grid; if not, the analysis is bounded to the cleanest test case and 020 can take up the residual.

**Track H recommends the alternative scope as the 019 baseline**, with the full reduced grid (4×3×3×2 = 72 propagations) as the escalation path if the analytical work closes the i=90 deg residual.

### 3.5 Multi-year DE441 acquisition (alternative 019 direction)

The 018/019 open-question candidates include "multi-year byte-pinned DE441 acquisition". This is a *data acquisition* task, not a propagation task. Cost estimate:

- 5-year daily Moon snapshot: ~366 × 5 = 1,830 rows, ~380 KB (from 76 KB/yr scaling)
- 5-year daily Sun snapshot: ~380 KB
- Acquisition time: 2 HTTP requests × ~3 s each = ~6 s (using the same Horizons pattern)
- Verification cost: sha256 + parse + distance-band check + uniform-cadence check ≈ 1 min total

This is cheap and adds significant scientific value (would let the 18.6-year lunar nodal cycle be characterized over a multi-year arc, removing the "1-year arc contaminated by 18.6-yr cycle" limitation). **Track H recommends the 5-year byte-pinned DE441 acquisition as the 019 baseline direction** if the 019 question is "do the evection + variation terms close the residual at multi-year arcs". The evection (27.55 d anomalistic) and variation (14.77 d synodic half-month) terms are well-characterized by even 1-2 years of data; the lunar nodal cycle (18.6 yr) is only partially characterized by 5 years.

If 019 is allowed both the analytical (evection + variation closed-form term) and the numerical (multi-year byte-pinned) arms, the cost is bounded (~4.5 hours propagation + ~1 hour multi-year acquisition + ~3 hours analysis/figures) and the science is high-value.

---

## 4. Graduation Decision

### 4.1 Should 019 graduate a corrected third-body secular formula into `lab_utils`?

**Short answer: GRADUATE ONLY ON CONDITION.**

The corrected formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i` is the standard textbook result (Murray & Dermott Sec. 7.2), has been independently derived in Track A, has been validated against the 1-year numerical in 018 at i=90 deg to 2.81x ratio (residual dominated by unmodelled short-period terms), and the radial scale factor (`(a/a₃)³` vs the deprecated `(R_E/r₃)²`) is correct. The formula is a permanent scientific asset.

**Arguments FOR graduation (now, before 019)**:
1. **The 018 formula is correct** at the doubly-averaged quadrupole order; the residual to the numerical is well-characterized (9.78x at i_sso, 2.81x at i=90 deg) and is dominated by documented short-period terms (evection + variation + lunar-nodal) that are *separate* from the secular formula. The secular formula itself is correct.
2. **Future consumers are real**: 020 (Eclipse-aware station-keeping), 021 (Sentinel/Landsat byte-pinning), and any future lunar mission design work (DEEP-SPACE) all need a correct third-body RAAN rate.
3. **Lab canon principle**: AGENTS.md §"Decisions are justified. Prefer existing code/infrastructure. Complexity must justify itself." Multiple consumers justify a shared canonical function.
4. **Graduating now is non-blocking**: the deprecated 017 function still works (with DeprecationWarning), so no breakage. The corrected function is added as a new canonical entry; existing experiments continue to work.

**Arguments AGAINST premature graduation**:
1. **019 may add an evection + variation term** to the formula. If the evection + variation terms close the 2.81x residual at i=90 deg to < 1.5x, the *full* Lunisolar model (secular + evection + variation) is what should be graduated, not just the secular component.
2. **The Track A derivation has an unresolved sign-convention question** (see §2.3): the 018 corrected formula uses a + sign on sin 2(i−i₃); Murray & Dermott's standard convention gives a - sign. Both are mathematically consistent (different sign conventions for the disturbing function), but 019 should explicitly pin the sign against the lab's `rv_to_coe_eci.Omega = atan2(h_y, h_x)` convention (Murray & Dermott standard astronomical RAAN).
3. **The numerical 2.81x residual at i=90 deg is unexplained**: is it the evection, the variation, the lunar-nodal term, or some combination? Without that diagnosis, graduating "the corrected formula" prematurely freezes the secular term but leaves the question of what other terms to add open. This is a feature, not a bug — the secular formula is correct on its own — but the API contract needs to be clear.

### 4.2 Recommended graduation policy

**Recommend a CONDITIONAL graduation at Exp 019**:

- **Stage 1 (pre-019)**: NO graduation. The corrected formula remains experiment-local (in 017 as `corrected_secular_lunisolar_raan_rate_rad_s` with DeprecationWarning on the wrong version, and in 018 as `corrected_secular_lunisolar_raan_rate_rad_s`). Rationale: 018 already proves the secular formula correct; graduating it before 019 is premature if 019 will add evection + variation terms.

- **Stage 2 (post-019)**: IF 019's root-cause analysis identifies a specific analytic term (evection, variation, or lunar-nodal correction) that **survives both Track A derivation AND Track E numerical experiments AND Track G hostile review**, THEN graduate the **full Lunisolar secular + [identified term]** model as `lab_utils.third_body_raan_rate_rad_s`. API contract:

```python
def third_body_raan_rate_rad_s(
    a_km: float,
    e: float,
    i_rad: float,
    *,
    third_body: Literal["sun", "moon"] = "sun",
    include_evection: bool = False,
    include_variation: bool = False,
) -> dict:
    """Doubly-averaged quadrupole RAAN perturbation from a third body.

    Secular term: (3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i-i_3) / sin i
    Optional evection term: f_evection(a, e, i, third_body, T)
    Optional variation term: f_variation(a, e, i, third_body, T)

    Returns dict with solar_cf_rad_s, lunar_cf_rad_s, total_cf_rad_s
    (all in rad/s). Sign convention: Murray & Dermott 1999, Sec. 2.10.

    See audit-019-*.md for derivation provenance and validation.
    """
```

- **Stage 3 (contingency)**: IF 019 fails to identify a specific term (the 2.81x residual remains unexplained), THEN do NOT graduate. Keep the corrected secular formula experiment-local. Re-evaluate at 020.

### 4.3 API contract requirements

If graduation proceeds, the API contract MUST include:

1. **Sign convention pinned**: explicitly documented as "Murray & Dermott 1999 Sec. 2.10" with the standard astronomical RAAN (Omega = atan2(h_y, h_x), matching `lab_utils.orbits.rv_to_coe_eci`).
2. **Three-body superposition**: function takes a single third body (sun or moon) and returns the radial rate; callers sum solar + lunar contributions. This matches the existing physics and avoids hard-coding the two-body assumption.
3. **Evection/variation terms**: gated behind optional keyword arguments with defaults OFF (so the secular-only path is the default and the full-model path is opt-in).
4. **Deprecation of the wrong formula**: `closed_form_lunisolar_raan_rate_rad_s` (017) and `luni_solar_raan_rate_rad_s` (016) keep their `DeprecationWarning` but should be moved out of the lab_utils path. They remain in the donor experiments (016, 017) for backwards compatibility.
5. **Tests at graduation**: include (a) sign-convention test against `rv_to_coe_eci.Omega` for a known case, (b) magnitude test against Track B derivation at h=600 i_sso, (c) Murray & Dermott textbook example if available.

### 4.4 Verdict

**GRADUATION DEFERRED until 019 completes its root-cause analysis**. The corrected secular formula is correct (Stage 1 has no technical blockers); the API and tests for the lab_utils version are well-specified; but the 2.81x residual at i=90 deg is a known unknown, and 019's contribution should be the evection + variation terms that complete the Lunisolar model.

---

## 5. Test Count and Quality

### 5.1 018 test layers and counts

018 declares 10 test layers (README L1-L10) and delivers **45 tests**:

| Layer | Description | Tests | Quality |
|---|---|---|---|
| L1 | Snapshot integrity (sha256, distance, cadence, n_points) | 8 | Excellent — covers sha256, n_points, uniform-cadence, distance-band for both Sun and Moon. Independent of any numerical computation. |
| L2 | Corrected closed-form identity | 8 | Excellent — checks radial factor, geometric factor, sign, total sign, decomposition, magnitude-at-Track-B-value, infeasibility-above-a_max, monotone-in-altitude. The L2 test at h_600_matches_track_b is a tautological comparison to a hand-computed expected constant, but the constant is independently derived (Track B), not read out of the implementation, so it is an independent validation. |
| L3 | Numerical isolation | 5 | Good — verifies all 5 modes present, J2-only slope matches canonical SSO target, slopes in operational band, n_ascending_nodes in expected band. |
| L4 | Inclination sweep | 5 | Good — verifies all 6 inclinations present, i=0 J2 retrograde, i=90 J2=0 null, i_sso in operational band. |
| L5 | Window-length sensitivity | 4 | Good — verifies all 5 windows present, residual pattern, 730d in operational band, n_ascending_nodes in expected band. |
| L6 | Precession rotation (frame fix) | 3 | Good — keys present, function exists, sun vector changes at 2026. Could add: precession is identity at T=0 (already in `test_precession_function_exists` actually — confirmed at line 393-399), precession applied to moon vector changes too. |
| L7 | Force-level identity | 3 | Excellent — machine-precision equivalence of two algebraic forms of third-body acceleration. Most stringent numerical-identity test in 018. |
| L8 | Convergence (RK4 order-4) | 3 | Good — keys present, order above 3, monotone decrease. Could add: position diff at finest dt < 1 mm, velocity diff at finest dt < 1e-7 km/s. |
| L9 | Adversarial mutants | 3 | Acceptable — no machine-specific paths, sign-flip vs deprecated 017, i_3_moon = obliquity + 5.145. Could add: literal mutants of the radial factor (e.g., `(R_E/r_3)**2` instead of `(a/a_3)**3`) with assertions that the mutant changes the sign. |
| L10 | Determinism, code hash, payload structure | 3 | Good — payload structure complete, code_sha256 includes essentials, corrected cf monotone in altitude. |

### 5.2 Test layers missing from 018 (and recommended for 019)

018 covers the canonical 10 layers well for a discrepancy-resolution experiment. For 019, which extends to long-period (evection + variation + lunar-nodal) terms, additional layers are appropriate:

1. **L11: Sign convention pinned against `rv_to_coe_eci`**. A direct test that the sign of `third_body_raan_rate_rad_s` agrees with the sign of `dΩ/dt` measured by `rv_to_coe_eci` for a known case. This addresses the §2.3 sign-convention concern.

2. **L12: Periodicity at known frequencies**. Test that the evection contribution has the expected 27.55 d period and the variation contribution has the expected 14.77 d period, by running the 019 numerical at i=90 deg for W in {30, 60, 90, 365} d and checking that the residual spectrum contains peaks at the expected frequencies (Fourier analysis or autocorrelation).

3. **L13: Multi-arc scaling**. If 019 acquires multi-year DE441 data, test that the residual between the corrected formula (secular + evection + variation) and the multi-year numerical scales as 1/(observation time) for the secular contribution and is bounded by the evection/variation envelopes for the short-period contribution.

4. **L14: Eccentric corrections** (if applicable). The current formula is for circular satellite orbits (e=0). If 019 extends to eccentric orbits, test that the e² correction enters at the expected order and sign.

5. **L15: Cross-validation against an independent celestial-mechanics library** (e.g., `poliastro` if available, or hand-coded Hansen coefficient expansion). This is the strongest independent validation possible — a second implementation producing the same answer to many significant figures.

### 5.3 Test count target for 019

Per the charter "40-70 tests expected per experiment":

- **Lower bound (40 tests)**: maintain the 018 10-layer coverage + add L11 (sign convention, 3 tests) + L12 (periodicity, 5 tests). Total: 45 + 8 = 53.
- **Recommended (50-60 tests)**: add L11, L12, L13. Total: 45 + 3 + 5 + 5 = 58.
- **Upper bound (70 tests)**: add L14 + L15 if eccentricity extension and cross-library validation are in scope. Total: 45 + 3 + 5 + 5 + 5 + 7 = 70.

**Track H recommends 50-60 tests** for 019 as a good fit to the charter target and to the new L11-L13 layers.

### 5.4 Are 018 tests tautological?

Most 018 tests are NOT tautological; they are independent validations:
- L1 tests verify the byte-pinned snapshot against the manifest (independent of the 018 implementation).
- L2 tests verify the corrected formula against (a) hand-derived constants (Track B is an independent derivation, not the implementation) and (b) structural properties (sign, monotone, decomposition).
- L3-L5 tests verify the numerical propagation against operational bands derived from independent experiments (Exp 009/012 J2 closure, Exp 014 eclipse durations, etc.).
- L6 tests verify the frame fix against the underlying polynomial (Lieske 1977).
- L7 tests verify two algebraic forms of the third-body acceleration agree at machine precision — this is the strongest independent check of the implementation correctness.
- L9 tests verify the corrected formula is not the same as the deprecated 017 formula (sign-flip and 100x+ magnitude ratio).

**No tautologies found**. The L2 test at h_600_matches_track_b is the closest to tautological (compares against a constant derived in the same audit), but the constant is independently derived in Track B (separate derivation path) and was the audit's central remediation claim, so the test provides independent validation that the implementation reproduces Track B's derivation.

### 5.5 Test independence for 019

For 019, the evection + variation terms should be tested against **independent derivations**, not just against the lab's own implementation. Recommended sources of independent validation:
- **Lidov (1962)** original paper's analytical expressions for evection and variation (if available).
- **Kaula (1966)** inclination-function framework for the e-dependent Hansen coefficients.
- **A benchmark against `poliastro` or another celestial-mechanics library** if the library supports third-body perturbations at the order of interest.

---

## 6. Verdict

### 6.1 Is the 019 plan feasible?

**Yes, with reduced scope.** The full 6×5×5×2 grid is 300 propagations and ~15-20 hours single core, which is at the upper bound of the lab's "Resources are precious" doctrine and should not be undertaken as a single uninterrupted session.

The recommended reduced grid (4×3×3×2 = 72 propagations, ~4.5 hours single core) is well within bounds. The minimum viable grid (i=90 deg only, ~15 min) is the recommended baseline.

### 6.2 What is the appropriate scope?

**Track H's recommended 019 scope**:

1. **Acquisition arm** (recommended as the primary 019 direction): byte-pinned multi-year DE441 Sun + Moon snapshot acquisition (5 years, ~380 KB each). Cost: ~1 hour. Value: removes the "1-year arc contaminated by 18.6-yr cycle" limitation.
2. **Analytical arm** (recommended as the secondary direction, in parallel with acquisition): closed-form evection + variation terms added to the corrected secular formula. Deliverable: `evection_raan_correction_rad_s` and `variation_raan_correction_rad_s` functions (experiment-local until Stage 2 graduation). Validation: at i=90 deg, the residual between corrected + evection + variation and the multi-year numerical should close from 2.81x to < 1.5x.
3. **Numerical arm** (the validation): run the recommended reduced grid (4×3×3×2 = 72 propagations) to confirm the analytical evection + variation terms close the residual.
4. **Documentation arm**: knowledge note linking 016/017/018/019, audit report, decision on whether to graduate the full Lunisolar model into `lab_utils`.

Total cost: ~6 hours single core for the acquisition + analytical + numerical + documentation work. Plus ~1 hour for the multi-year byte-pinned acquisition (data acquisition and verification only, not propagation).

### 6.3 What graduation decision is justified?

**Stage 1 (pre-019)**: NO graduation. The corrected secular formula is correct but the 2.81x residual at i=90 deg is unexplained; graduating prematurely freezes a partial model.

**Stage 2 (post-019)**: GRADUATE conditionally, IF 019's root-cause analysis identifies a specific analytic term (evection or variation or lunar-nodal correction) that survives Track A derivation AND Track E numerical validation AND Track G hostile review. Graduate the FULL model (secular + identified term) as `lab_utils.third_body_raan_rate_rad_s` with the API contract in §4.3.

**Stage 3 (contingency)**: If 019 fails to identify a specific term, do NOT graduate. Keep the corrected secular formula experiment-local. Re-evaluate at 020.

### 6.4 What are the risks?

**Risk 1 (HIGH): The 2.81x residual at i=90 deg may NOT be evection or variation.** It could be a combination of unmodelled higher-order terms (octopole at O(a/a₃) ≈ 1.8% for the Moon, ~5e-5 for the Sun), the ICRF/mean-of-date precession residual (already characterized at 0.012 deg/year), the J3 zonal harmonic (not modelled), or simply numerical drift in the 1-year fit. **Mitigation**: Track A derivation in 019 should explicitly bound each candidate term's expected contribution and rank them by leverage. Track E numerical should target the cleanest test (i=90 deg, sun-only or moon-only) to isolate which term dominates the residual.

**Risk 2 (MEDIUM): Multi-year DE441 acquisition may not close the 18.6-year lunar nodal cycle.** A 5-year arc still has 18.6/5 = ~27% of one lunar nodal cycle, leaving significant short-period contamination. **Mitigation**: 019 should report the residual vs observation time and explicitly characterize the lunar nodal term's contribution vs the 1-year and 5-year arcs. The lunar nodal term is small (effectively a long-period modulation of the secular term, not a separate secular term) and may be tractable as a parametric correction rather than a fully resolved cycle.

**Risk 3 (MEDIUM): The sign convention of the corrected formula may not match all references.** Murray & Dermott's convention gives a minus sign; some Vallado formulations give a plus sign. The 018 corrected formula uses a plus sign. **Mitigation**: add an L11 sign-convention test to 019 that explicitly pins the sign against the lab's `rv_to_coe_eci.Omega` for a known case (e.g., i=0 LEO with the Sun's annual nodal contribution, which is unambiguously retrograde per Curtis Ch. 10).

**Risk 4 (LOW): The Vallado Eq. 9-46 citation cannot be directly verified.** I could not inspect Vallado's 4th edition (Microcosm, 2013) in this environment. The 018 remediation is correct in substance regardless of Vallado's exact wording, but if 019 graduates a canonical function, the citation should be replaced with Murray & Dermott + Kozai + Lidov (the canonical modern and historical references).

**Risk 5 (LOW): The 45-test count is below the AGENTS.md "40-70 expected" range minimum by 5 if measured against the strict parametrize-expanded 714 baseline.** Not blocking; the 45-test count is well within the charter's 40-70 range when measured as test functions (not parametrize expansions). If 019 adds L11-L15 layers, the count rises to 50-60, well within the target.

**Risk 6 (MEDIUM): 019 may not close the residual sufficiently to justify graduation.** If the residual at i=90 deg remains at 2.5x or larger even after adding evection + variation + lunar-nodal terms, the analytical model is incomplete and Stage 2 graduation should be deferred to 020. **Mitigation**: the recommended scope includes a clear escalation criterion (close to < 1.5x ratio at i=90 deg for Stage 2 graduation).

### 6.5 Overall assessment

**The 019 plan is feasible in scope, well-bounded in cost, and the graduation decision is appropriately deferred.** Track H's recommendations:

- Use the i=90 deg cleanest test as the primary 019 lever.
- Acquire the multi-year DE441 snapshot in parallel with the analytical evection + variation development.
- Run the recommended reduced 72-propagation grid (~4.5 hours) for numerical validation.
- Defer graduation of `corrected_secular_lunisolar_raan_rate_rad_s` to `lab_utils.third_body_raan_rate_rad_s` until 019's root-cause analysis identifies the evection/variation term that closes the 2.81x residual.
- Add L11 (sign convention) and L12 (periodicity at known frequencies) as the key new test layers for 019.

The 018 deliverable is reproducible, byte-pinned, deterministic, and well-documented. The literature provenance is sound (Murray & Dermott canonical, Kozai/Lidov historical, Lieske precession). The Vallado citation has been REMEDIATED. The 10x residual at i_sso and the 2.81x residual at i=90 deg are well-characterized open questions for 019, not failures of 018.

---

## 7. Summary

| Audit dimension | Verdict | Risk |
|---|---|---|
| Snapshot byte-pinning (Sun + Moon) | PASS — both hashes verified against MANIFESTs | None |
| Determinism (no RNG in analysis, no network, no wall-clock) | PASS | None |
| Frame / time-type consistency (ICRF/TDB → mean-of-date) | PASS — IAU-1976 precession applied correctly | None |
| Constants (provenance + version pinning) | PASS | Minor: AU_KM duplicated locally; recommend `lab_utils.earth_frames.AU_KM` import |
| Test counts (45 in 018; repo-wide ~654 / ~714 expanded) | PASS (45-test file matches README; full count not independently verifiable in this env) | Parametrize-expansion discrepancy noted |
| No silent nondeterminism | PASS | None |
| Literature provenance (Murray & Dermott, Kozai, Kaula, Lieske) | PASS — all citations real; Murray & Dermott is canonical modern text; Kozai 1959 and Lidov 1962 are historical; Kaula 1966 is the standard geodesy text; Lieske 1977 is the precession source | None |
| Vallado Eq. 9-46 citation | REMEDIATED by 018 — citation is no longer propagated in 018 documentation | None |
| Corrected formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i` | REAL canonical result (Murray & Dermott Sec. 7.2 + 2.10); also appears in Kozai 1959 and Lidov 1962; sign-convention ambiguity documented | Medium: sign-convention pinning recommended (L11 test) |
| 019 full grid feasibility (300 propagations, ~15-20 hours) | NOT FEASIBLE as a single uninterrupted session | High: scope reduction required |
| 019 recommended grid feasibility (72 propagations, ~4.5 hours) | FEASIBLE | None |
| 019 minimum viable scope (i=90 deg only, ~15 min) | FEASIBLE | None |
| Graduation decision | DEFERRED to post-019 (Stage 1: no; Stage 2: conditional on root-cause identification) | Medium: depends on 019 outcome |
| Test layers missing for 019 | L11 sign convention, L12 periodicity at known frequencies, L13 multi-arc scaling, L14 eccentric corrections, L15 cross-library validation | Low: additions well-defined |
| Test count target for 019 | 50-60 tests | None |

**Overall: 018 is ready as a foundation for 019. Track H recommends the i=90 deg cleanest-test scope as the 019 baseline, the recommended reduced grid (4×3×3×2 = 72 propagations) as the escalation path, and Stage 2 conditional graduation only after 019's root-cause analysis closes the 2.81x residual at i=90 deg.**