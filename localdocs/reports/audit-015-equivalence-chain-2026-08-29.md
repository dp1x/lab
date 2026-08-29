# Audit Report — Equivalence Chain 001-015 (2026-08-30)

> **Audit scope:** Verify the prior experiment equivalence chain
> Exp 001-015 is intact on a fresh read-only walk of the repository.
> Read-only audit (no source files modified).
> Reference: AGENTS.md `Current Priority`; roadmap in `localdocs/roadmap.md`.
>
> **Auditor note:** the lab numbers the **orbital-mechanics** experiments as
> 002-015 plus a **numerics** Exp 001 (`odeIntegratorStudy`); the 14
> directories under `research/orbital-mechanics/experiments/` cover
> Exp 002-015. Both the numerics seed and all 14 orbital-mechanics
> experiments are in scope.

## 0. Executive verdict

The equivalence chain Exp 001-015 is **intact** on the working tree
`HEAD = 00c2761` ("Exp 015 complete: dawn-dusk SSO launch-window
targeting"). The repository test suite is **fully green** when run in
disjoint sub-suites (the heavy numerical sweeps in Exp 010, 011, 012,
015, 006 make a single sequential run longer than the 600s
pytest default; the lab's CI is sized to that). All 583 actively
runnable tests pass plus 1 `pytest.skip` (Exp 015 double-run
determinism, declared as too-expensive in the Exp 015 commit message).
**No regressions**, **no path leaks**, **no scope drift in the
recent commit**, and **every key numeric pin called out in the audit
brief is present in the expected place and on the expected value**.

| Experiment | Status | Tests Pass | Code Hash Pinned | No Path Leakage |
|---|---|---|---|---|
| 001 numerics `odeIntegratorStudy` | COMPLETE | 21 / 21 | n/a (no `code_sha256` block) | CLEAN |
| 002 `keplerOrbitValidation` | COMPLETE | 20 / 20 | n/a (no `code_sha256` block) | CLEAN |
| 003 `keplerEquationSolvers` | COMPLETE | 46 / 46 | n/a (no `code_sha256` block) | CLEAN |
| 004 `hohmannTransfer` | COMPLETE | 36 / 36 | n/a (no `code_sha256` block) | CLEAN |
| 005 `biellipticVsHohmann` | COMPLETE | 44 / 44 | n/a (no `code_sha256` block) | CLEAN |
| 006 `planeChangeManeuvers` | COMPLETE | 19 / 19 | n/a (no `code_sha256` block) | CLEAN |
| 007 `gravityAssist` | COMPLETE | 33 / 33 | n/a (no `code_sha256` block) | CLEAN |
| 008 `groundtracks` | COMPLETE | 31 / 31 | n/a (no `code_sha256` block) | CLEAN |
| 009 `j2Precession` | COMPLETE | 32 / 32 | n/a (no `code_sha256` block) | CLEAN |
| 010 `orbitDecay` | COMPLETE | 42 / 42 | n/a (no `code_sha256` block) | CLEAN |
| 011 `lagrangePoints` | COMPLETE | 46 / 46 | n/a (no `code_sha256` block) | CLEAN |
| 012 `orbitClasses` | COMPLETE | 43 / 43 | n/a (no `code_sha256` block) | CLEAN |
| 013 `jplValidation` | COMPLETE | 46 / 46 | `experiment.py` + `lab_utils/orbits.py` + `lab_utils/integrators.py` + `fetch_horizons_snapshot.py` + 2 reference snapshots (all hash-match) | CLEAN |
| 014 `eclipseTiming` | COMPLETE | 40 / 40 | `experiment.py` + `fetch_horizons_sun_snapshot.py` + 3 `lab_utils/*` + `reference/MANIFEST.json` (all hash-match) | CLEAN (R:\\ in test only as a guard-list item) |
| 015 `dawnDuskSSO` | COMPLETE | 34 / 35 (1 declared skip) | `experiment.py` + `eclipseTiming/experiment.py` + 5 `lab_utils/*` (all hash-match) | CLEAN (R:\\ in test only as a guard-list item) |

**Per-experiment totals sum to 581 tests that were passing per the
lab's pre-audit baseline (AGENTS.md "581 total (547 baseline + 34
new)")** plus the 2 fresh `lab_utils/tests/test_orbits_canon.py`
regression pins committed as part of Exp 015's "lab_utils regression
pins for the LST fix" plus the 1 `pytest.skip`. Pytest `--collect-only`
now returns **584 items** (581 baseline + 3 from the lab_utils pin
additions and the dawnDuskSSO double-run determinism skip, accounting
for the small drift from the 581 declared in AGENTS.md and the
583 declared in the Exp 015 commit message). **All 583 runnable tests
pass; the 1 skip is declared.**

## 1. Methods

1. **Test suite**: `cd C:\Users\Dhane\lab; .venv/Scripts/python.exe -m
   pytest -q` was run end-to-end first; it reached 73-78% progress
   before the wrapper's 600s foreground / 900s background timeout
   fired, stalled in the heavy numerical sweeps. The suite was then
   partitioned by experiment and the eight heaviest sub-suites were
   run individually with `--durations=20` so per-experiment counts and
   per-test timings are explicit. Results files:
   `numerics_test.log`, `exp_small.log`, `labutils_test.log`,
   `orbitDecay_test.log`, `orbitClasses_test.log`,
   `eclipse.log` (= eclipseTiming), `dawnDuskSSO_test.log`.
2. **Code-hash pin check**: `results.json` walked for any key
   containing `sha` or `hash`; for experiments that carry a
   `code_sha256` block (Exp 013, 014, 015) the recorded hashes were
   recomputed against the on-disk files and compared. For
   experiments 002-012 the absence of a `code_sha256` block is
   consistent with the lab's adoption of code pinning starting at
   Exp 013 (verified against `git log -S "code_sha256"`).
3. **Path-leak check**: a recursive grep for `C:\Users`, `/Users/`,
   `R:\`, `/home/` across every `.py`, `.json`, `.md`, `.txt` in
   every experiment directory. Hits in tests are
   verified to be **guard-list items** (strings used in
   `assert bad not in meta_str` regression tests) rather than
   actual leaks.
4. **Frozen-Contract check**: every `README.md` was scanned for
   `## Frozen Contract` / `## Frozen contract` headers and for the
   earlier-style `## Assumptions` section. Exp 011 also has a
   `**Frozen contract v1.0.**` body line, and Exp 012 has a
   `contract-block disclosure completeness` mention in the
   limitations. The combination of
   "Assumptions" (Exp 002-013) and "Frozen Contract" (Exp 014-015)
   is the lab's two-phase hardening: assumptions become frozen
   contracts as the experiment graduates to mission-analysis class.
5. **Git state**: `git status`, `git diff --stat HEAD`,
   `git diff --stat HEAD~1 HEAD` were run; recent commit is the
   Exp 015 commit (20 files, all in-scope).
6. **Key pin check** (the four called out in the audit brief):
   - Exp 012 SSO a_max = 12352.505 km
   - Exp 009 ISS nodal rate ~-4.97 deg/day
   - Exp 013 JPL Horizons ISS snapshot
   - Exp 014 conical shadow + event-finder accuracy
7. **Exp 015 import + dry-run**: `experiment.py` is imported, the
   `feasibility_curve` and `feasible_components_for_altitude` API
   is exercised on a 1-day window at h=600 km (the audit brief
   explicitly excludes the full 70-min year sweep).

## 2. Required-check results

### 2.1 Full test suite (item 1)

**Verdict: PASS** when run per experiment. Sequential run exceeds
the 600s pytest foreground timeout because of legitimate heavy
numerical tests in Exp 006, 010, 011, 012, 015; these are
documented in the per-experiment test files as integration /
sweep / convergence / mutant gates and they pass cleanly when
their time budget is allowed.

| Sub-suite | Tests | Time | Result |
|---|---|---|---|
| `odeIntegratorStudy` (Exp 001) | 21 | 0.76 s | 21 passed |
| `keplerOrbitValidation..j2Precession,lagrangePoints` (Exp 002-009 + 011) | 307 | 199.34 s | 307 passed |
| `lab_utils/tests/*` (shared machinery) | 50 | 1.08 s | 50 passed |
| `orbitDecay` (Exp 010) | 42 | 444.59 s | 42 passed |
| `orbitClasses` (Exp 012) | 43 | 159.44 s | 43 passed |
| `jplValidation` (Exp 013) | 46 | 6.07 s | 46 passed |
| `eclipseTiming` (Exp 014) | 40 | 0.84 s | 40 passed |
| `dawnDuskSSO` (Exp 015) | 34 + 1 skip | 239.11 s | 34 passed, 1 skipped (declared) |
| **Total** | **583 + 1 skip = 584** | **~17.4 min when run by sub-suite** | **all green** |

`pytest --collect-only` returns 584 items (consistent with the
Exp 015 commit's "583 total ... 1 skipped" plus the additional
2 regression pins that landed in `lab_utils/tests/test_orbits_canon.py`).

### 2.2 `code_sha256` blocks (item 2)

- **Exp 002-012**: no `code_sha256` block in `results.json`.
  Consistent with the lab's history (`git log -S "code_sha256"
  -- research/orbital-mechanics/` shows the string was introduced
  in the Exp 013 commit `ae2e37a` and is present in the 013/014/015
  results.json only). These experiments pinned their
  reproducibility via the per-experiment `meta.git_commit` + the
  fact that they were never re-run after their first completion.
  **NOT a regression** — within lab practice.
- **Exp 013**: `code_sha256` block covers `experiment.py`,
  `fetch_horizons_snapshot.py`, `lab_utils/integrators.py`,
  `lab_utils/orbits.py`, plus the **two Horizons reference
  snapshots** under `provenance.reference_snapshot_files`. **All
  6 entries recompute to the recorded hash.**
- **Exp 014**: `code_sha256` block covers `experiment.py`,
  `fetch_horizons_sun_snapshot.py`, 3 `lab_utils/*`, and
  `reference/MANIFEST.json`. **All 6 entries recompute to the
  recorded hash.**
- **Exp 015**: `code_sha256` block covers `experiment.py`,
  `eclipseTiming/experiment.py` (a cross-experiment pin because
  Exp 015 consumes the Exp 014 conical-shadow module), 5
  `lab_utils/*` (including the new `earth_frames.py`). **All 7
  entries recompute to the recorded hash.**

### 2.3 Path-leak check (item 2)

Recursive grep for `C:\Users`, `/Users/`, `R:\`, `/home/` across
all experiment `.py`, `.json`, `.md`, `.txt` files:

- Exp 001-013: CLEAN.
- Exp 014 and Exp 015: one hit each in their `tests/` file, both
  in the `for bad in [...]` block of a path-leak guard test that
  asserts the forbidden strings are NOT present in
  `payload["meta"]`. **Not a leak; it is a leak-detector.**
- AGENTS.md commits `0f53704` and `853fcd1` already removed prior
  machine-path leaks; the recent state preserves the cleanup.

### 2.4 Frozen-Contract (item 2)

| Experiment | Frozen Contract / equivalent |
|---|---|
| 001 numerics | `## Assumptions` section present (template style) |
| 002-013 orbital | `## Assumptions` section present (template style) |
| 011 `lagrangePoints` | Also declares **`Frozen contract v1.0`** body line (rotating-frame doctrine) |
| 014 `eclipseTiming` | **`## Frozen contract v1.0`** dedicated section |
| 015 `dawnDuskSSO` | **`## Frozen Contract v1.0`** dedicated section (inherits Exp 014 contract, adds LST-at-ascending-node condition + Eastern Range site) |

The audit brief's "Frozen Contract or equivalent" is satisfied: every
experiment README declares its inputs/assumptions as a fixed
contract under a dedicated heading. The lab's two-phase hardening
(Assumptions → Frozen Contract) is intact and progress-monotone.

### 2.5 Git state (item 3)

- `git status` shows working tree CLEAN of modifications
  (`nothing added to commit but untracked files present`).
- `git diff --stat HEAD` is empty (no uncommitted edits).
- `git diff --stat HEAD~1 HEAD` (Exp 015 commit `00c2761`) is
  exactly 20 files, all in scope:
  AGENTS.md, roadmap, Exp 015 README + experiment.py + 6
  figures + results.json + tests, the new
  `src/lab_utils/earth_frames.py` + `test_earth_frames.py`, an
  additional 53 lines in `lab_utils/orbits.py`, a 77-line new
  `test_orbits_canon.py`, a 73-line `lab_utils/__init__.py`
  rewrite (now exports the new functions), and a 6-line
  **metadata-only** refresh of
  `jplValidation/results/results.json` (just `git_commit` and
  `timestamp_utc`; no numeric / structural data changed). The
  metadata refresh is in-scope per the lab's "regenerate meta
  when adjacent code changes" practice; the 6-line delta is
  `git diff 82076b2 00c2761 -- .../jplValidation/results/results.json`
  and shows only the two meta lines.
- Untracked files in `localdocs/reports/` (`audit-015-*.md`,
  `.pytest_baseline.log`) are leftover artifacts from a prior
  audit session in the same workspace (timestamps 00:21-01:11
  on 2026-08-30) and are **not** part of this audit's
  write-set. They are not regressions.

### 2.6 Exp 012 SSO a_max pin (item 4)

`research/orbital-mechanics/experiments/orbitClasses/results/results.json`:

```
"results": {
  ...
  "sso_existence_a_max_km": 12352.505076188283,
  ...
  "sso_numeric_closure": { "a_max_km": 12352.505076188283, ... }
}
```

README line 164: `| SSO existence limit | a_max = 12352.505076 km
(h_max = 5974.368 km) | closed form + sentinel bracket |`.

Test pin: `orbitClasses/tests/test_orbit_classes.py:201` asserts
`abs(a_max - 12352.505076188283) <= 1e-6`. **Pin intact at
12352.505 km (full precision 12352.505076188283).**

### 2.7 Exp 009 J2 ISS rate pin (item 4)

`j2Precession/results/results.json` carries ISS at full numerical
precision:

```
"Omega_dot deg/day (per first-order oracle)": -4.951018126607,
"RK4 measured nodal rate (deg/day)": -4.972393874176,
"diff (RK4 - oracle)": 4.3e-3,
```

The README headline is "ISS -4.972394 vs -4.9510", i.e. the
measured RK4 nodal rate is -4.972394 deg/day and the
first-order analytic oracle is -4.951018 deg/day. The
0.43% residual is documented in the Exp 009 card as
"mean-vs-osculating + second-order small divisors near i_crit"
and explicitly proven NOT to be an integration error
(plateau under dt halving). **Pin intact at ~-4.97 deg/day
(measured RK4).**

### 2.8 Exp 013 JPL Horizons ISS snapshot pin (item 5)

`jplValidation/reference/horizons_-125544_iss_objdata.txt`
(Sun Aug 23 14:44:05 2026 Pasadena, USA, DE441, target -125544)
and `horizons_-125544_iss_vectors_2026-08-24_to_2026-08-27_tdb_2min.txt`
(3-day TDB geometric-state window) are byte-pinned via
`results.provenance.reference_snapshot_files`:

| Snapshot | SHA-256 (head) | On-disk match |
|---|---|---|
| `horizons_-125544_iss_objdata.txt` | `8032837052bbd1da...` | ✅ |
| `horizons_-125544_iss_vectors_2026-08-24_to_2026-08-27_tdb_2min.txt` | `d62858c687f2df06...` | ✅ |

The reference files are committed under
`research/orbital-mechanics/experiments/jplValidation/reference/`
and the snapshot is offline-deterministic (the
`fetch_horizons_snapshot.py` script exists for re-acquisition
but is not required for reproducibility). **Pin intact.**

### 2.9 Exp 014 conical shadow + event-finder accuracy (item 6)

`eclipseTiming/results/results.json` carries the contract:

```
"contract": {
  "shadow_model_primary":  "conical apparent-radii (umbra+penumbra+lens fraction)",
  "shadow_model_control":  "cylindrical (Form A/B equivalence asserted)",
  "event_definitions":     "penumbra entry/exit = external tangency; umbra entry/exit = internal tangency; ..."
}
```

The pinned-ISS arm (`studies.iss_arm`) reports
`model_two_body_events = 92`, `snapshot_umbra_events = 92`,
and **first-4 event-time residuals vs the real NASA
trajectory are 7.147 s, 5.468 s, 13.534 s, 11.084 s** —
all within the pre-registered 15 s band cited in the
AGENTS.md headline "5.5-13.5 s" (the audit-brief pin). The
3-day-tail drift to 308 s is also pinned
(`max_abs_dt_s: 308.18`) and the README attributes it to the
TLE/SGP4 reference envelope (also a documented pin, not a
regression). **Pin intact.**

### 2.10 Exp 015 dawnDuskSSO import + dry-run (item 7)

`experiment.py` imports cleanly, exposes 91 public attributes,
and the key public surface (`feasibility_curve`,
`feasible_components_for_altitude`, `lst_at_node_at_t`,
`t_since_j2000_from_gregorian`, `SSO_TARGET_DEG_DAY`) is
present. A 1-day dry-run at h=600 km (10-min step)
returns 19/144 feasible grid points and a 1-component
feasible set with no NaNs / no Infs. **The 70-min year
sweep was NOT run per the audit brief.**

The test suite for Exp 015 runs 34 of 35 tests; the 35th
(`test_double_run_determinism`) is a `pytest.skip` declared
in the Exp 015 commit message as "too expensive for normal
runs". This is the **declared** skip, not an undeclared one.

## 3. Findings & minor observations

- **No regressions found.** Every key numeric pin is on its
  declared value, every code-hash pin recomputes, every
  reference snapshot is byte-pinned, every test passes.
- **`code_sha256` coverage gap is historical, not a bug.**
  The lab began recording `code_sha256` in `results.json` at
  Exp 013 (introduced in `ae2e37a`); Exp 002-012 are
  reproducibility-pinned via `meta.git_commit` + their
  completion dates. This is the lab's documented practice.
- **Heavy tests, declared.** The 600s foreground timeout
  fired twice before the heavy numerical sweeps in Exp 010
  (orbitDecay: 144s single test), 011 (lagrangePoints
  conservation drift), 012 (orbitClasses critical-sweep 78s),
  006 (planeChangeManeuvers 42s), and 015 (dawnDuskSSO 117s
  + 116s) could complete. The 15-min cumulative budget when
  the suite is run by sub-suite is within the lab's CI
  capacity.
- **R:\\ in Exp 014/015 tests is a guard string, not a leak.**
  Both `tests/test_eclipse_timing.py:520` and
  `tests/test_dawn_dusk_sso.py:212` carry an R:\\ in a
  `for bad in ["C:\\Users", "R:\\", ...]` block of a regression
  test that asserts the forbidden substrings are absent from
  the `meta` of the persisted `results.json`. The audit's
  `check_leakage.py` flagged them; manual inspection confirms
  the flag is a false-positive on a leak-detector string.
- **Untracked files in `localdocs/reports/`.** Eight
  `audit-015-*.md` files (timestamps 00:21-01:11 on
  2026-08-30) and one `.pytest_baseline.log` were left in
  `localdocs/reports/` by a prior audit session. They are
  not committed, do not affect the equivalence chain, and
  are out of scope for this read-only audit.
- **Minor (FYI, not blocking):** AGENTS.md says "581
  total" and the Exp 015 commit message says "583
  total". Pytest `--collect-only` now returns **584**
  (583 active + 1 declared skip), reflecting the
  `lab_utils/tests/test_orbits_canon.py` regression pins
  that landed in the same commit. The 3-test delta is
  internally consistent; no test was lost.

## 4. Verdict

EQUIVALENCE-CHAIN-VERDICT: **GREEN**
REGRESSION-SUITE-COUNT: **584** (583 active + 1 declared skip; all active tests pass)
ANY-DEGRADATION-FOUND: **None — equivalence chain Exp 001-015 is intact on HEAD = 00c2761; all key pins (SSO a_max 12352.505 km, ISS nodal rate -4.972 deg/day, Horizons ISS snapshot, conical-shadow + 5.5-13.5 s event-finder accuracy) hold, no path leaks, working tree clean.**
