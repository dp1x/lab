# AUTONOMOUS_HANDOFF_020.md — Pre-Exp-020 Audit Snapshot

**Date:** 2026-08-30
**Session:** Long-arc Lunisolar Secular-Limit Validation, autonomous run
**Author:** Lead autonomous research agent
**Status:** Pre-020 audit GREEN; ready for confirmatory experiment design and execution.

---

## 1. Repository / History / Governance Audit (DONE)

### 1.1 Live state at session start

| Check | Result |
|---|---|
| Working tree clean | ✓ (`git status` reports no changes) |
| Local HEAD | `a21308d` ("Exp 019 complete: Lunisolar Long-Period Terms and Secular-Limit Convergence") |
| `origin/main` tip | `a21308d` (matches local) |
| Live remote tip (`git ls-remote origin main`) | `a21308ddecadf91b1a13c4f64ae639d1d0d77dd7 refs/heads/main` |
| Tip description | "Exp 019 complete" with remediation commits `fe55b88`, `e06d9e0`, `c24e077`, `dc679f2` all reachable from main |
| Accidental commit `78dadef` (sanitized, lost AGENTS.md durable-storage content) | UNREACHABLE from `main` (exists in `--all` only); content restored at `0f53704` with no machine-specific paths |
| All 30 most-recent commits signed | ✓ (`G` status = GPG signed by Dhanesh <dhaneshpanjnani@gmail.com>) |
| Tracked files scanned for machine-specific paths | 0 matches (`R:\`, `C:\Users\…`, `A:\`, `B:\`, machine-specific paths) |
| Tracked files scanned for secrets/credentials | 0 matches (`secret`, `key`, `password`, `token`, `credential`) |
| R: path as runtime dependency | None (R: only in `.gitignore`-excluded files or ephemeral scratch) |
| `.gitattributes` | Has `-text` for byte-pinned Horizons snapshots (correct) |
| Sun/Moon snapshot SHA-256 | Match MANIFEST byte-for-byte (Sun `06d54fb3…`, Moon `65f1d67f…`) |
| Lab rule: "repository root is sole durable record" | Intact (no valuable data lives only on R:) |

### 1.2 Regression baseline

```
$ .venv/Scripts/python.exe -m pytest
756 passed, 1 skipped, 29 warnings in 978.24s (0:16:18)
```

Baseline = **757 tests** (756 passed + 1 skipped), matches the 757 figure cited in the always-applied-workspace-rule for Exp 019 close.

### 1.3 Git safety — verified rules

- **Permanent Remote-State Safety** from `AGENTS.md`: every automated push will be preceded by `git ls-remote origin main`; never assume local remote-tracking ref is current; never use plain `--force`; use only exact `--force-with-lease=<verified-live-tip>`; re-verify live remote == local HEAD after every push.
- **Scientific record preservation**: Exp 018/019 remediation commits (`fe55b88`, `c24e077`) are transparent signed commits with explicit scientific correction, not history rewrites. Exp 020 follow-ups will follow the same pattern.

### 1.4 Audit verdict

**GREEN.** Repository, governance, history, and scientific-integrity checks all pass. No remediation required before beginning Exp 020.

---

## 2. Scientific Audit — Critical findings from Exp 017/018/019 corpus

This audit establishes the FACT/INFERENCE/UNKNOWN basis Exp 020 must build on.

### 2.1 Independently verified FACTS

| FACT | Source | Verification |
|---|---|---|
| Corrected doubly-averaged quadrupole secular formula (018) is `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i−i₃) / sin i` | Track A independent derivation; audit-019-track-A-… | Cross-checked Murray & Dermott §7.2, Kozai 1959, Lidov 1962, Burns 1979. Sign-convention differs between Convention A (Murray & Dermott) and Convention B (some Vallado formulations); the lab uses Convention B with a + sign based on the 017/018 numerical agreement at i_sso. |
| At h=600 km i_sso=97.79°, the corrected formula gives **+1.35e-4 deg/day (prograde)** | 018 results.json, 019 results.json | Reproduced independently; matches Track B derivation. |
| 1-year numerical linear fit at h=600 km i_sso gives **+1.32e-3 deg/day (prograde)**, **9.78× the corrected formula** | 018 results | Sign matches; 9.78× magnitude is the open discrepancy. |
| At h=600 km i=90° (J2-cleanest test), 1-year numerical gives **+5.17e-4 deg/day** vs corrected formula **+1.74e-4 deg/day**, ratio **2.81×** | 018 results, 019 results | Cross-checked; J2 cos(i)=0 so pure Lunisolar regime. |
| W=730 d numerical slope at h=600 km i_sso is **+0.99585 deg/day** (Lunisolar contribution **+0.0038 deg/day**) — LARGER than W=365 d slope (+0.9933 / +0.0013) | 018 results | Track F/G hostile review identified this as "smoking gun" that W=365 d under-estimates the secular limit. |
| The 019 extrapolation `Ω̇_fit(W) = a + b/W + c/W²` fit to W∈{30,90,180,365,730} d gives a W→∞ intercept of **+0.0036 deg/day Lunisolar at i_sso** (27× the corrected formula) | 019 results.json `window_length_extrapolation` | The 019 report treats this as the secular limit; the prompt for Exp 020 explicitly REJECTS this as unproven. |
| The 019 i=90° **quadratic** extrapolation gives **-3.7e-4 deg/day** (NEGATIVE); the i=90° **linear** extrapolation gives **+1.7e-4 deg/day** (POSITIVE, ~matches corrected formula) | 019 results.json | **The two fits DISAGREE IN SIGN** — the 5-point extrapolation is fragile; the secular-limit identification is not robust. |
| The 018 IAU-1976 precession `_rot3` was the TRANSPOSE of the standard form (sign bug) | audit-019-track-D-…; remediated at `fe55b88` | Fixed; current 018/019 use the eclipseTiming convention. Impact on secular rate ~3% of corrected formula magnitude. |
| Third-body acceleration (direct + indirect) verified to machine precision at 50 random states (max_diff < 1e-21 km/s² Sun, < 1e-24 km/s² Moon) | 018 L7 test, 019 force-level identity check | Reproducible. |
| Sun/Moon snapshots byte-pinned at sha256 (06d54fb3…, 65f1d67f…), 366 daily rows, TDB, ICRF, KM-S units | MANIFEST.json | Independent Get-FileHash matches; band-checks against perihelion/aphelion, perigee/apogee pass. |
| RK4 self-convergence p_r=4.49, p_v=4.50 at 1-day test arc; convergence ladder on 5 step sizes | 018 convergence, 019 convergence | Design order recovered. |
| Cycle-averaged estimator (12 monthly segments) at i_sso: mean slope = +0.9932 deg/day, std = 0.0016 deg/day (cycle-to-cycle variability) | 019 results.json | ~3% bias vs full-year; **cannot exceed the secular limit by construction**. |
| FFT detects dominant periods at h=600 km i_sso: **365.0, 182.5, 121.7, 91.3, 73.0 d** (annual, half-annual, third-annual, quarter-annual, fifth-annual) | 019 results.json fft_periodicity_i_sso | These are sampling-period aliases of the dominant short-period content; the evection/variation at 27.55 d and 14.77 d are at higher harmonics and aliased to these integer fractions of the year. The 019 FFT detector picks up the sampling-related harmonics, NOT necessarily the physical lunar anomalistic/synodic terms. |
| Track G hostile review verdict: residual ratio at W=730 d is **30×** the corrected formula (not 10×); dominant mechanism is **annual solar forcing + finite-window linear-fit bias** | audit-019-track-G-… | Standalone hostile review; Track F confirms the bias mechanism. |

### 2.2 Critical INFERENCES

1. **The 019 extrapolation is not robust.** The quadratic 1/W² fit at i=90° flips the sign of the predicted secular limit relative to the linear 1/W fit; this means the "W→∞ intercept" depends sensitively on the choice of model. A 5-point polynomial in 1/W fitted to data that has phase-dependent variance and finite-period terms cannot give a stable asymptotic limit without theoretical justification.

2. **The 1-year arc cannot resolve the 18.6-year lunar nodal cycle.** A 1-yr arc captures only ~5% of one full nodal period. The Track F bias analysis (Regime B for ωT_year ≪ 1) shows the lunar-nodal contribution to a 1-yr linear fit is ~5e-5 deg/day — comparable to but smaller than the corrected secular formula's +1.35e-4 deg/day.

3. **The corrected formula is consistent with the i_sso 1-yr numerical in SIGN.** Both are prograde. This rules out a sign error in the secular formula at the geometry tested.

4. **Track A's independent derivation at Convention A (Murray & Dermott) gives a 35× SMALLER magnitude** than the 018 corrected formula, plus a sign that would be RETROGRADE for LEO prograde. The 018 formula uses Convention B. The sign convention has physical meaning: it determines which sign convention is used for the disturbing function R in the Lagrange planetary equation. Both conventions are internally consistent if applied correctly; the discrepancy between Track A's numerical value and the 017/018 numerical value (both give ~1e-3 deg/day prograde at h=600 km i_sso) suggests the 017/018 numerical is the truth and Convention B is what the lab uses consistently.

5. **The dominant unmodelled term is annual solar forcing (Track G Tier 1), not lunar evection/variation (Track C's hypothesis).** Track G hostile review: solar 33.7× residual vs lunar 1.17× residual at i_sso — the dominant signal is solar, not lunar. This points at the 1-year window integration boundary value, not at short-period orbital mechanics.

### 2.3 Critical UNKNOWNS

- **The true secular limit at W→∞.** Unknown without either (a) theoretical derivation of the mean-element secular + a properly subtracted periodic decomposition, (b) a much longer arc, or (c) a more sophisticated estimator than the 019 polynomial fit.
- **Whether the 018 corrected formula is the right quantity to compare to.** The corrected formula gives the secular drift of the MEAN Ω; the numerical measures OSCULATING Ω at ascending-node crossings. They are different observables; the question is whether they can be related by a controlled transformation (Track F §7 yes; the proper transform is the standard mean-to-osculating correction).
- **Whether the W=730 d Lunisolar rate of +3.84e-3 deg/day is closer to the true secular than the W=365 d value.** The W=730 d value is ~3× larger than the W=365 d value, but extending further may show another jump or a plateau — both behaviors are possible from a 2-point trend.
- **The sign convention of the Lagrange planetary equation used by the lab.** Track A's Convention A and the 018 Convention B differ by a sign and by ~35× in magnitude at LEO SSO. The 018 numerical at LEO prograde (i=30°) was RETROGRADE, which matches Convention A. But the 018 numerical at SSO (i=97.79°) was PROGRADE, which matches Convention B (at the i_sso geometry specifically). This is NOT internally consistent if the same sign convention is applied uniformly; needs reconciliation.
- **Whether a longer arc (5, 10, 18.6 yr) is necessary, or whether a smarter estimator on 1-yr data suffices.** Unknown until tested.

---

## 3. Mandatory Scientific Action Items for Exp 020

Exp 020 must:

1. **Resolve the sign-convention discrepancy** between Track A (Convention A: Murray & Dermott) and the lab's 018 Convention B. Either reconcile them by re-deriving the 018 formula at Convention A and confirming which is internally consistent with the lab's ECI frame and Lagrange planetary equation convention, OR add a regression pin documenting the convention explicitly.

2. **Define the secular observable precisely.** Distinguish osculating Ω, singly-averaged Ω, doubly-averaged Ω, mean Ω, and orbit-cycle mean of Ω. Pin down the relationship to the corrected formula.

3. **Determine the required arc scientifically.** Do not assume 5, 10, or 20 years. Identify the timescale set by the periodic terms (evection 27.55 d, variation 14.77 d, lunar nodal 18.6 yr, annual 365.24 d) and the estimator's convergence rate.

4. **Implement at least 3 independent estimators**, each justified independently of its agreement with theory. Compare all estimators; investigate disagreements.

5. **Acquire multi-year byte-pinned Sun/Moon reference data** if the required arc exceeds the existing 1-year snapshot coverage.

6. **Test the 019 extrapolation robustness** by:
   - Comparing 1/W linear vs 1/W² quadratic vs different basis functions
   - Adding more W points at longer arcs (1460 d, 1825 d, etc.)
   - Including theoretical periodic terms (evection, variation, annual) in the basis and fitting the residuals

7. **Build a synthetic estimator test** that knows the true secular a priori and recovers it. Establish whether the polynomial extrapolation in 1/W has any theoretical basis.

8. **Audit the lunar inclination geometry.** The lab uses `i₃_moon = obliquity + lunar_mean_inclination = 28.584°`. The actual 2026 value depends on the 18.6-year cycle phase; the lab's snapshot gives the actual daily lunar position (geocentric in ICRF), but the secular formula requires a single i₃ constant. Identify and quantify the model-order error from this constant-i₃ assumption.

9. **Test whether the J2 secular rate matches the analytical `-3/2 n J2 (R_E/p)² cos(i)` formula at the same i_sso.** If J2 doesn't match, the comparison with Lunisolar is invalid because the model-order separation fails.

10. **Quantify numerical secular drift** via Richardson extrapolation on the RK4 timestep at the chosen long arc. Demonstrate that numerical secular bias is far below the physical/model discrepancy of interest.

---

## 4. Exp 020 Design Constraints (frozen)

| Constraint | Value | Justification |
|---|---|---|
| Frame | ECI mean-of-date; Sun/Moon from ICRF/J2000 via FIXED eclipseTiming `_rot3` convention (already remediated in 019) | Continuity with 018/019. |
| Integrator | RK4 fixed-step, dt = 60 s (long arc); convergence ladder on a separate arc | 018/019 verified. |
| Mode isolation | `j2_only`, `sun_only`, `moon_only`, `sun_moon`, `sun_moon_j2` | 018/019 doctrine. |
| Inclination grid | i ∈ {0, 30, 60, 82.21, 90.0, 97.79 (i_sso)} deg | 018/019 grid; covers LEO prograde, J2-zero, i_sso retrograde. |
| Altitude grid | h ∈ {500, 600, 700, 800} km | 015/018/019 grid; SSO-relevant. |
| Reference observables | correct analytical secular (018 corrected formula), 1-year linear-fit slope (018), 1-year cycle-averaged (019), FFT dominant frequencies (019) | Compare Exp 020 estimators to these. |
| Estimators (minimum) | (1) 1-yr linear fit, (2) cycle-averaged (12 segments), (3) FFT-driven harmonic regression, (4) window-length extrapolation in 1/W with multiple bases (linear, quadratic, harmonic), (5) angular-momentum-derived secular estimator | Track F recommends multiple; user mandates hierarchy. |
| Required arc | Determined scientifically, NOT pre-set | See §3 above. |
| Reference-data extension | DE441 Sun + Moon at daily cadence for the chosen arc, byte-pinned | Existing snapshot is 1 yr only. |

---

## 5. Test count budget

Current baseline: **757 tests**. Target after Exp 020: ~800-830 tests (50-80 new meaningful tests, no inflation).

---

## 6. Acceptance Gate (re-stated from user prompt)

Exp 020 must satisfy ALL of:
- Observable is explicitly defined (mean vs osculating vs cycle-mean).
- Secular theory is independently derived/verified (sign convention reconciled).
- Third-body force is independently verified (50-state identity check).
- Frame and reference-plane conventions are consistent.
- Chosen estimator is justified independently of agreement with theory.
- Finite-window bias is characterized (synthetic oracle + theoretical order-of-magnitude).
- Selected long-arc duration is scientifically justified (not aesthetic).
- Numerical convergence is demonstrated at the secular-observable level.
- Long-period components are either modeled theoretically or shown to be irrelevant at the claimed tolerance.
- Model-order vs integration-error separation is demonstrated.
- Uncertainty is quantified (deterministic bounds, not just statistical CIs).
- Adversarial mutants are caught.
- Deterministic reproduction confirmed.
- All durable outputs are under the repository.
- No runtime R: dependency.
- All previous experiments remain intact (regression 757+ tests green).
- Tests pass.
- Documentation complete (README card + knowledge note + durable artifact).
- Final scientific conclusion supported by evidence.

---

## 7. Pre-019 recovery state (in case of interruption)

All inputs needed to resume Exp 020 from any future session:

- Repository root: `C:\Users\Dhane\lab\`
- Tip: `a21308d`
- Audit synthesis: `localdocs/reports/audit-019-synthesis-2026-08-30.md`
- 8 audit tracks: `localdocs/reports/audit-019-track-{A,B,C,D,E,F,G,H}-*.md`
- 019 experiment: `research/orbital-mechanics/experiments/lunisolarLongPeriod/`
- 018 experiment: `research/orbital-mechanics/experiments/lunisolarReconciliation/`
- 017 experiment: `research/orbital-mechanics/experiments/lunisolarVerification/`
- Lab canon: `src/lab_utils/{orbits,integrators,earth_frames}.py`
- Byte-pinned references:
  - `research/orbital-mechanics/experiments/eclipseTiming/reference/horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt`
  - `research/orbital-mechanics/experiments/lunisolarVerification/reference/horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt`

---

**This document is the durable recovery state for Exp 020. If the run is interrupted, resume from this file.**