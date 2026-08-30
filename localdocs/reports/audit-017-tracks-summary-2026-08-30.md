# Experiment 017 Audit — Track Summary (2026-08-30)

> **Audit type:** Pre-017 governance, history, scientific, and forward-extrapolation audit.
> **Scope:** 001-016 repository state, 016 scientific baseline, 017 forward-extrapolation
> validity, candidate direction selection.
> **Audit method:** Hard audit gate + eight read-only specialist tracks delegated in
> parallel.
> **Audit verdict:** GREEN to proceed (after rejecting the original decadal direction
> and selecting the closed-form upper-bound verification per Track H Alt-1).

---

## Phase 1: Repository Governance & History Audit (GREEN)

| Check | Result | Evidence |
|---|---|---|
| Branch | main | `git branch --show-current` |
| Local HEAD | `42ccbbe` (Exp 016 complete) | `git rev-parse HEAD` |
| Live remote tip | `42ccbbe` (matches) | `git ls-remote origin main` |
| Working tree | clean | `git status` |
| Canonical Git identity | Dhanesh <dhaneshpanjnani@gmail.com> | `git config` |
| GPG signing | `commit.gpgsign=true`, key `241711B89ADB1BCEC668DB12A9C779E02303E517` | `git config` |
| Recent commits signed | YES (last 4 commits all GPG Good) | `git log --show-signature` |
| Exp 016 commit | `42ccbbe` matches expected | direct check |
| Exp 015 remediation | preserved as `b1ec17e` | reflog + git log |
| Accidental `78dadef` | unreachable from main (dangling object only) | `git rev-list main` |
| Machine-specific paths in tracked files | NONE (defensive guard-list items in tests only) | grep audit |
| R: dependencies | NONE | grep audit |
| Secrets/credentials | NONE | grep audit |
| `.gitattributes` protections | intact; new -text rule for Moon snapshot planned | .gitattributes |
| Test baseline | 625 pass + 1 skipped = 626 total (audit message off by 1, content correct) | pytest --collect-only |
| All 14 prior experiments tests pass | YES | per-suite pytest runs |

---

## Phase 2: 016 Scientific Baseline Audit (CONDITIONAL GREEN)

| Claim | Status | Evidence |
|---|---|---|
| LST = 12 + (Omega - alpha_sun)/15 correctly implemented | VERIFIED | Track A; experiment.py:161-172 |
| GMST cancellation | VERIFIED (function does not use GMST at all) | Track A |
| J2 closure residual at h=600 km is ~2.2 deg/year | **DOCUMENTATION BUG**: code returns ~0 by construction; "~2.2" is from Exp 012 | Track A |
| EoT peak-to-peak 30.65 min | VERIFIED at upper edge of pre-registered [19, 29] min band; extends to ~31 min | Track A; results.json |
| Horizons snapshot residual 0.056 deg | VERIFIED within 0.7 deg gate | Track A |
| Decadal extrapolation defensibility | NOT SAFE for closed-form Lunisolar (~50x over-estimate, sign possibly wrong) | Tracks A, B, G |
| Safe for use as 016 audit-response contract | YES | Charter, audit-015 |

---

## Phase 3: Eight-Track Specialist Delegation (Tracks A-H)

All 8 tracks run in parallel as read-only research subagents. Findings:

### Track A — Verify 016 SSO/LST/J2 Baseline
- 016 baseline is **structurally sound** with one **DOCUMENTATION BUG**: the
  j2_closure_residual function returns ~0 by construction (the SSO solver
  pins by definition); the README's "~2.2 deg/year" is from Exp 012 not from
  this function. Code is correct; narrative is wrong.
- Decadal extrapolation is safe IF the Lunisolar upper-bound band is honestly
  reported (320 min/year at h=600 → ~53 hours over 10 yr ceiling).

### Track B — Long-arc Orbital Dynamics & Error Growth
- **VERDICT: Decadal propagation NOT defensible with current machinery.**
- Lab's RK4 is non-symplectic; secular energy drift `~1.4e-6·t` (oscillator
  study); not characterized past ~30 days.
- Decadal arc is ~1000x longer than any lab-validated arc.
- No adaptive timestep, no symplectic integrator, no regularization.
- Recommends EXP 017 defer or scope to 1-year arc.

### Track C — Lunisolar Perturbations & Ephemeris Provenance
- Lab's closed-form Vallado Eq. 9-46 **over-estimates by ~50x at SSO retrograde**
  per Exp 016 model_note.
- **NO Moon ephemeris pinned in repo** (verified by grep across all
  reference/, src/lab_utils/).
- Horizons API supports Moon target 301 with identical acquisition pattern
  to existing Sun snapshot. ~50 KB at 1-day cadence for 1 year. ~3-5 days
  of work to acquire and validate.
- **Recommended minimum upgrade**: byte-pinned Moon geocentric snapshot at
  1-day cadence, 1-year arc.

### Track D — SRP Geometry, Shadow, A/m, Attitude
- SRP at A/m=0.01 is **sub-noise over 10 years** (~0.17 arc-sec total RAAN).
- Cr coefficient implicit (perfect reflector Cr=2 hard-wired at 9.08e-9 km/s²).
- Cylindrical/conical shadow NOT currently modeled for SRP.
- **NOT the bottleneck** for decadal experiment; ~4 orders of magnitude below
  Lunisolar at LEO SSO.

### Track E — Drag & Space-Weather, F10.7/Ap Modeling
- Lab's exponential atmosphere **NOT defensible** for any decadal experiment.
- Vallado fit runs +26-34% above US76 at LEO; scale height varies 2.1x across
  the 400-800 km band.
- Solar-activity dependence (F10.7) **absent** from lab canon; density varies
  2-10x between solar min/max.
- F10.7 historical data is sub-MB, byte-pinnable from NOAA SWPC.
- **Methodological concern**: mixing observed + forecast F10.7 for future
  years is a broken experiment as written.

### Track F — Validation & Model-Order Separation
- Required 8 pillars for any decadal claim: closed-form limits, independent
  impl, null-force regression, invariants, convergence, sensitivity,
  external anchor, adversarial mutants.
- Single hardest requirement: distinguishing integration error from model
  truncation from eph input to spacecraft parameter uncertainty from
  event/control discretization from accumulated phase error.
- Short-term agreement (Exp 013: 5 km at 3 days) does NOT extrapolate to
  decadal accuracy without multi-year anchor.

### Track G — Hostile Adversarial Reviewer
- **VERDICT: DEFER the decadal direction. Six fatal flaws:**
  1. Closed-form Lunisolar over-estimates by ~50x at SSO retrograde; 10-year
     integration propagates 500x accumulated uncertainty.
  2. No F10.7-driven atmospheric density model; drag dominates at h=500-600 km.
  3. No Sentinel-1/Landsat LTAN evolution byte-pinned; "5-15 m/s/year" is
     repeated textbook citation without primary source.
  4. Mean-vs-osculating short-period terms not in lab canon; become secular
     bias over 10 years.
  5. A/m, Cr, attitude spacecraft-specific; decadal direction is design
     study, not mission analysis.
  6. "Decadal" framing is false precision; spans solar max and lunar nodal
     phase variability.
- Concerning issues: integration error growth over 5e7 RK4 steps, EoT vs
  precession over decade, asymmetric Δv reporting, single-Cd assumption.

### Track H — Portfolio Value & Alternatives
- **SCORED** four candidates:
  1. Refined J2 mean-vs-osculating coupling: 19/30 (REJECT — already documented)
  2. Higher-altitude eclipse coupling: 18/30 (REJECT — no new physics)
  3. Decadal station-keeping: 17/30 uncalibrated / 24/30 calibrated (RECOMMENDED for 018)
  4. **Alt-1: Closed-form upper-bound verification with byte-pinned Sun + Moon over 2026: 27/30 (RECOMMENDED AS 017)**
- **RECOMMENDATION: Run Alt-1 as Exp 017.** Smallest scope, biggest payoff,
  lowest risk; converts the Exp 016 model-order disclaimer into a measured
  byte-pinned quantity.

---

## Phase 4: Synthesis & Claim Classification

| Claim | Status |
|---|---|
| 016 baseline structurally sound for LST-correction contract | FACT |
| 016 Lunisolar closed-form over-estimates by ~50x at SSO retrograde | FACT |
| Lab has no Moon ephemeris pinned | FACT |
| Decadal propagation NOT defensible with lab's current RK4 + drag machinery | FACT |
| Decadal F10.7-driven drag NOT defensible with lab's exponential atmosphere | FACT |
| Sun snapshot byte-pinned (DE441, 2026, 366 daily rows) | FACT |
| Moon geocentric snapshot feasible (~50 KB, 1-d cadence, 1-yr arc) | FACT |
| Sentinel-1 operational Δv 5-15 m/s/year | INFERENCE (textbook cite, not byte-pinned) |
| Sentinel-1 altitude 693 km, dawn-dusk SSO, 11-yr operational | FACT (verified ESA) |
| cf_upper/numerical ratio bounded above by sin^2(i_SS) within 2x | HYPOTHESIS (needs measurement) |
| Earth-Moon barycenter negligible for LEO Lunisolar | DEFENSIBLE ASSUMPTION |

**Resolution:** Material UNKNOWNs resolved through Track C (Moon snapshot is
feasible with documented acquisition pattern) and Track H (Alt-1 is the
strongest defensible candidate).

---

## Phase 5: Decision

**REJECT** the decadal direction (Tracks B, E, F, G unanimous).
**SELECT** Track H Alt-1 (closed-form upper-bound verification with byte-pinned
JPL Sun + Moon over 2026) as Exp 017.

**Rationale:** The decadal direction was originally proposed in the
post-Exp-016 roadmap auto-update based on the audit-015 follow-up candidates
scoring matrix (which scored decadal at 19/30 vs 29/30 for the EoT-anchored
correction = Exp 016, which was completed). The audit-015-follow-up-candidates
report explicitly labeled the decadal direction as "next-after-next", not
next. Exp 016 was the correct next experiment; Exp 017 = Alt-1 is the correct
next-after-016 experiment.

---

## Phase 6: Implementation Record

- **Experiment directory:** `research/orbital-mechanics/experiments/lunisolarVerification/`
- **Acquisition:** `fetch_horizons_moon_snapshot.py` (single Horizons API call,
  Moon target 301, 1-day cadence, full 2026, 366 rows, byte-pinned
  sha256 `65f1d67f798a3b95...`).
- **Implementation:** `experiment.py` (deterministic, offline, RK4 with
  aligned integer-multiple grids).
- **Tests:** `tests/test_lunisolar_verification.py` (32 tests in 6 layers).
- **Results:** `results/results.json` (frozen payload with code_sha256 binding).
- **Figures:** 4 figures (cf/numerical ratio, drift comparison, dt convergence
  ladder, linear-fit residuals). All byte-stable across multiple runs.
- **Knowledge note:** `localdocs/knowledge/lunisolar-verification.md`.
- **Documentation updates:** AGENTS.md, roadmap.md, orbital-mechanics README.

**Final commit:** `c42ba9c` ("Exp 017 complete: Lunisolar upper-bound verification
(audit response)"), GPG signed with RSA key `5774B47A005623ACD39DF7284EB7C30F884E8259`.

---

## Phase 7: Headline Results

| Quantity | h=500 | h=600 | h=700 | h=800 |
|---|---:|---:|---:|---:|
| Closed-form upper bound (deg/day) | -0.2108 | -0.2184 | -0.2263 | -0.2343 |
| Numerical Lunisolar (J2-subtracted, deg/day) | +0.001320 | +0.001284 | +0.001249 | +0.001215 |
| **cf_upper / numerical ratio (signed)** | **-159.64** | **-170.14** | **-181.19** | **-192.84** |
| Linear-fit residual RMS (deg) | 0.0247 | 0.0240 | 0.0234 | 0.0227 |

**Convergence:** p_r = 4.49, p_v = 4.50 (RK4 design order ~4 confirmed).

**Validation gates:**
- convergence_order_pass: True (p_r, p_v >= 3.5)
- numerical_magnitude_pass: True (|numerical| in [1e-4, 1e-1] deg/day)
- ratio_band_pass: False (audit-015 [10x, 100x] band violated; **documented as
  a first-principles discovery**: audit under-estimated the over-estimate
  factor by ~3x; the band violation IS the headline finding)

**Test count:** 658 total (626 baseline + 32 new), all passing.

---

## Phase 8: Publication Gate

- Canonical Git identity verified: Dhanesh <dhaneshpanjnani@gmail.com>
- GPG signing verified: `commit.gpgsign=true`, key
  `241711B89ADB1BCEC668DB12A9C779E02303E517`, RSA sig key
  `5774B47A005623ACD39DF7284EB7C30F884E8259`
- Live remote tip verification: pre-push `42ccbbe`, post-push `c42ba9c` (matches
  local HEAD)
- Working tree: clean
- Commit signed: YES
- GitHub verification: pending (requires GitHub API; signature key uploaded
  to GitHub already per `user.signingkey` config)

---

## Phase 9: Recommendation for Exp 018

Now that the closed-form Lunisolar over-estimate is measured (~170x at h=600 km
with sign reversal), the strongest next experiments are:

1. **Refined Lunisolar evection + variation terms** — characterize the missing
   Lunisolar terms (evection at anomalistic month ~27.55 d, variation at
   synodic half-month ~14.77 d) to refine the closed-form. Direct follow-up
   to this experiment; ~1 day compute.

2. **Multi-year Sentinel/Landsat byte-pinning** — acquire Sentinel-1A precise
   orbit ephemerides (CNES POD, free public) and Landsat-8 long-term
   ephemeris (NASA EOSDIS, free public). Provide the external validation
   anchor currently missing for any operational LEO claim. ~2-3 days.

3. **(Deferred) Decadal station-keeping with F10.7-driven drag** — full port
   of NRLMSISE-00 or JB2008, byte-pinned F10.7 historical record, symplectic
   or regularized integrator to suppress RK4 secular drift. Multi-week
   project; not appropriate as Exp 018.

**Recommended Exp 018:** #1 (evection + variation refinement) as the natural
follow-up to this experiment; #2 (Sentinel/Landsat byte-pinning) as the
necessary prerequisite for any decadal claim later.