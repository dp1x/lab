# Scientific stack audit — what the Lab hand-rolls, and what should be replaced

**Date:** 2026-10-02
**Type:** Capability audit (read-only; no dependency was installed, no source file was edited).
**Scope:** the whole repository — 96 tracked `.py` files, 41,110 lines, 29 test files, 828 collected tests.
**Authority:** `LAB_CONSTITUTION.md` §7 (knowledge model), §9 (FET gate 6, capability-advancing), AGENTS.md rule 7 ("decisions are justified").

---

## 1. FACT — the dependency surface (as originally found)

At the time of this audit the Lab depended on **four packages**, declared in `pyproject.toml`:

```toml
dependencies = ["numpy>=2.0", "matplotlib>=3.9", "mpmath>=1.4.1"]
[dependency-groups] dev = ["pytest>=8.0"]
```

> **SUPERSEDED 2026-10-02 (later the same day).** `scipy>=1.18.1` has since been
> added via `uv add scipy`, so `pyproject.toml` and `uv.lock` now both list five
> runtime packages. `uv` 0.12.22 was installed to the user profile to do it
> properly, keeping the lockfile in sync. See §9 for what was actually built and
> §5 for why the addition is hygiene rather than risk reduction. The evidence
> below is preserved as originally recorded.

Resolved in `.venv` (Python 3.12): `numpy 2.5.1`, `matplotlib 3.11.1`, `mpmath 1.4.1`,
`pytest 9.1.1` (+ transitive `pillow`, `fonttools`, `cycler`, `pyparsing`, `python-dateutil`,
`contourpy`, `kiwisolver`, `six`, `packaging`, `iniconfig`, `pluggy`, `Pygments`, `colorama`).

**`scipy` has never been a dependency.** Verified two ways: `pyproject.toml` has no entry,
and `uv.lock` contains **zero** `name = "scipy"` records. Same negative result for `astropy`,
`pandas`, `sympy`, `poliastro`, `scikit-learn`, `numba`, `jax`, `uncertainties`.

Import census across `research/` + `src/` (top-level imports, by frequency):
`lab_utils` 98, `numpy` 82, `pathlib` 78, `math` 62, `__future__` 61, `json` 55,
`matplotlib` 51, `importlib` 41, `sys` 32, `hashlib` 23, `time` 20, `pytest` 20,
`multiprocessing` 9, `re` 7, `mpmath` 7, `subprocess` 6, `datetime` 6.

**Interpretation:** every numeric routine in the Lab is either numpy, mpmath, or
hand-rolled. This is the finding — not a defect, but it is a *deliberate, undocumented*
position that has never been stated in a card or report.

## 2. FACT — what that position has cost

Twenty numbered experiments and four post-roadmap missions completed on **8 CPU cores,
no GPU, no cluster**, ≤4 hr wall per mission. That satisfies FET gate 1
(reasoning-constrained, not compute-constrained) and gate 7.

The cost is concentrated in a small number of hand-rolled numeric primitives, all of
which flow into headline numbers.

### 2.1 The estimator (`np.linalg.lstsq`, 22 sites)

The canonical implementation is
`mission_j2_precession_modulated_lunisolar_term/mission_experiment.py:265`:

```python
def ols_slope(t_s: np.ndarray, y_rad: np.ndarray) -> float:
    """Slope in rad/s."""
    A = np.column_stack([np.ones_like(t_s), t_s])
    coef, *_ = np.linalg.lstsq(A, y_rad, rcond=None)
    return float(coef[1])
```

It returns **a point estimate with no uncertainty**. Every "±" and every SNR in the
lunisolar chain is produced by a *different*, local re-derivation of the covariance —
e.g. `analyze_phase_b.py:74-80`:

```python
sigma2 = np.sum((luni_arr - A @ coeffs) ** 2) / (n_data - n_coeffs)
cov = sigma2 * np.linalg.inv(A.T @ A)
```

This is the textbook OLS covariance, which is **only valid under an i.i.d.-noise
assumption**. The inputs are not noisy: they are deterministic outputs of a
deterministic RK4 integrator. The residuals are *model misfit*, not measurement error.
The resulting ratio `|a11| / sigma` was reported as "SNR 6.89" and was treated as a
statistical significance claim.

This is recorded here as a **methodology finding**, not as a Mission 2 criticism —
the mis-specification was found and corrected by Mission 3, which is the evidence
doctrine working. But `scipy.stats.linregress` returns slope, standard error, and a
p-value from one tested call, and its stderr uses the same i.i.d. assumption — so
**adopting scipy would not have prevented this error.** The error was a *modelling*
error. See §5.

Site distribution: `lunisolarLongPeriod/experiment.py` (3), `mission_lunisolar_closure/experiment.py` (2),
`mission_j2_lunisolar_coupling/{mission_experiment,analyze_phase_b}.py` (2 each),
`lunisolarSecularLimit/experiment.py` (2), plus test files.

### 2.2 Bootstrap confidence intervals (12 sites)

Concentrated in `jplValidation/experiment.py` (8) and its tests (3). Used for the
skill-score CI quoted in Exp 013 ("J2 removes 99.33% of residual RMS, CI 0.9933–0.9979").

Hand-rolled resampling is a well-understood, low-risk routine. It is nonetheless
~50 lines of Lab-authored code whose only validation is self-consistency.

### 2.3 `np.fft.rfft` on unevenly-sampled data (1 site, 2 calls)

`lunisolarLongPeriod/experiment.py:505-506`:

```python
spec = np.abs(np.fft.rfft(om_detrended))
freqs = np.fft.rfftfreq(n, d=dt_mean)  # cycles/day
```

This assumes **uniform sampling**. The data are ascending-node crossings, whose spacing
is *not* uniform: across the 1-yr Phase 0 run, `n_cross` was 5445 / 5449 / 5462 for
different force modes — a ~0.3% spacing variation over ~5,400 samples.

**Sizing this honestly:** a 0.3% spacing jitter is small. The uniform-grid FFT is
*approximately* correct here and this is **a refinement opportunity, not a bug**. It was
previously flagged as a potential smoking gun; that flag was itself an overclaim and is
withdrawn. The correct tool is `scipy.signal.lombscargle`, which is designed for exactly
this sampling.

### 2.4 `np.polyfit` (9 sites)

`orbitClasses` (2), `keplerOrbitValidation` (2), `jplValidation` (2), `eclipseTiming` (1),
plus tests. All are simple linear/quadratic trend fits where `np.polyfit` is
already the right, documented numpy API. **No replacement needed.**

### 2.5 `erfi` / quadrature (16 sites)

`orbitDecay/experiment.py` (12), `jplValidation` (1), plus tests. These use `mpmath`
special functions for the drag integral. **`mpmath` is already the correct choice** —
arbitrary-precision, well-tested for exactly this. **Keep.**

### 2.6 The integrator (hand-rolled RK4, `lab_utils.integrators`)

`rk4_step` / `rk4_propagate` — used by all 20 experiments. Fixed step (dt = 60 s
throughout the lunisolar chain).

`scipy.integrate.solve_ivp` would add adaptive stepping, dense output, and native
event detection. **However** — see §4. Changing the integrator would invalidate the
`equivalence pinning` that `AGENTS.md` and the constitution require for graduated
machinery, and would change every committed result number in the repository. This is
the **single most dangerous** proposed substitution and is explicitly **rejected**.

---

## 3. The three-bucket classification

### Bucket 1 — REPLACE with a tested library (justified)

| Item | Sites | Replacement | Justification |
|---|---|---|---|
| Bootstrap CI machinery | 12 | `scipy.stats.bootstrap` | Removes ~50 Lab-authored lines of resampling code; gains documented method dispatch and returned CI object. Moderate value. |
| Node-crossing event detection | 4 missions | `scipy.integrate.solve_ivp(events=...)` | *Partially* — see §4. Event detection via root-finding on a dense interpolant is strictly better than linear interpolation between RK4 steps. Can be adopted **without** changing the integrator by keeping `rk4_step` and solving the crossing on a local cubic. Medium value, medium risk. |
| Periodogram on node samples | 1 | `scipy.signal.lombscargle` | Correct tool for uneven sampling (§2.3). Low risk, contained to one function. |

### Bucket 2 — KEEP Lab-owned (justified)

| Item | Reason |
|---|---|
| `erfi` / drag quadrature (`mpmath`) | Already the right library. High-precision special functions, exactly the use case. |
| `np.polyfit` trend fits | Already a documented, tested numpy API. Replacing it adds a dependency for no gain. |
| Fixed-step RK4 (`lab_utils.integrators`) | Deterministic, order-verified (Exp 001/008/009 report measured order 4.06–4.95), equivalence-pinned across experiments. Swapping it would change every published number. |
| `ols_slope` point estimate | A two-parameter least-squares slope is not a gap. The *covariance* around it is the gap — and that is Bucket 1 / §5. |
| DE441 snapshot loader + interpolation | Regression-pinned after a real bug; byte-pinned inputs with sha256 manifests. Domain-specific, validated, load-bearing. |

### Bucket 3 — AMBIGUOUS (decision required, not taken)

| Item | Question |
|---|---|
| Centralized OLS-with-uncertainty helper | 22 `lstsq` sites each hand-roll their own covariance. Centralizing into `lab_utils.estimation.ols_slope_ci()` would give one tested implementation and one documented caveat about i.i.d. assumptions. This is also **constitution backlog candidate #1** ("Estimation Doctrine Graduation"). **Recommended** — but it must expose the assumption explicitly, not hide it. |
| `np.fft` → `lombscargle` | Whether to *re-run* Exp 019's periodicity arm, or record it as a known refinement. Re-running changes a committed result. |

---

## 4. What must NOT be replaced

`scipy.integrate.solve_ivp` as a **replacement integrator** is rejected. Reasons:

1. `AGENTS.md` requires graduated machinery to be "equivalence-pinned vs donors".
   Replacing the integrator makes that pinning impossible — there is no donor to pin to.
2. Every committed result JSON in the repository was produced by `rk4_step`.
   Changing it silently invalidates 34 result artifacts.
3. `mission_mean_element_leakage` is an **active mission with a pre-registered decision
   rule and a computed, pinned Phase 0 result**. Changing the numeric stack underneath a
   pre-registered mission breaks the provenance chain the mission rests on.

Adaptive stepping is nonetheless a *legitimate future capability mission* on its own
merits — it would matter for cislunar/NRHO arcs where a fixed dt is wasteful. It is
**not** justified as a retrofit to the LEO chain.

---

## 5. The uncomfortable finding

**Adding scipy would not have prevented the errors that actually occurred.**

| Error | Would scipy have caught it? |
|---|---|
| `use_j2 = mode != 'kepler_only'` force-mode contamination | **No** — logic bug |
| DE441 interpolation sign/weight error (a11 moved 5×) | **No** — data-pipeline bug |
| IAU-1976 precession `_rot3` transpose | **No** — convention bug |
| `a11` SNR 6.89 from a mis-specified 2nd-order polynomial | **No** — `linregress` has the same i.i.d. covariance assumption |
| Exp 015 sidereal-vs-SSO rotation-rate confusion | **No** — conceptual bug |
| Exp 017 three-compound closed-form error | **No** — derivation error |

Every serious error in this repository's history was a **physics, convention, or
modelling** error. Not one was caused by a weak numeric primitive.

**Therefore:** the library upgrade is *hygiene*, not *risk reduction*. Its real value is
(a) removing Lab-authored code that has no independent validation, and (b) giving one
audited place where "how do we estimate a rate, and what does its error bar mean?" is
written down. It does **not** raise the Lab's ability to catch the errors that have
actually bitten it. That capability came from the 8-track audit pattern, the
pre-registration, and the sanity-check gates — all of which are already present.

This finding is recorded here so that the dependency addition is justified on its true
merits and not oversold.

---

## 6. Stage D — `astropy` decision

**DECISION: DO NOT ADD `astropy`.** Written decision, per plan.

**Arguments for:** `astropy.time` is the professional answer for TDB/TAI handling;
`astropy.coordinates` handles frame transforms the Lab hand-rolls (IAU-76 precession
`_rot3`, GMST polynomial). Both classes of bug have already occurred in this repository
(a time/epoch question in Exp 013, a frame-convention question in Exp 019/018).

**Arguments against (decisive):**
1. The Lab works in **ICRF/TDB with a byte-pinned DE441 text snapshot**. There is no
   leap-second, UT1, or frame ambiguity in that setup — Exp 013 already pinned the
   epoch tags to ≤0.33 mm. The pain astropy would remove is pain the Lab does not have.
2. `astropy` is a heavy dependency (~30 transitive packages) for a single-domain Lab.
3. The hand-rolled frame code is **tested against a byte-pinned Horizons Sun snapshot
   to 0.056 deg** (Exp 016) and against pinned ISS states (Exp 013). It is validated.
4. AGENTS.md rule 7: "Complexity must justify itself." It does not, today.

**Revisit condition:** if a mission enters a regime where time-scale or frame
convention becomes a live scientific question — **cislunar/NRHO (mission backlog #3)
is the likely trigger**, since CR3BP + libration-point work is where frame definitions
actually start to bite — re-run this decision then.

---

## 7. Summary

| Finding | Status |
|---|---|
| Lab runs on numpy + matplotlib + mpmath + pytest only | FACT, verified against `pyproject.toml` + `uv.lock` |
| scipy absent historically, not merely uninstalled | FACT |
| 3 low-risk replacements justified (bootstrap, lombscargle, node event detection) | RECOMMENDATION |
| Integrator replacement rejected | RECOMMENDATION (hard) |
| astropy not adopted; revisit at cislunar | DECISION |
| Library upgrade is hygiene, **not** risk reduction | FINDING — see §5 |

**Stage A is complete. Stage C (dependency addition) is BLOCKED on toolchain — see §8.**

## 8. Stage C — RESOLVED 2026-10-02

Originally this stage was **blocked**:

- **`uv` was not installed** on this machine.
- **`pip` was not present inside `.venv`** (`No module named pip`).
- PyPI was reachable (HTTP 200).

`AGENTS.md` mandates "Python managed by uv — `uv sync` first" and "add deps with
`uv add`". Installing into the venv with a bootstrapped `pip` would have left
`uv.lock` out of sync with `pyproject.toml` — a reproducibility regression in a
repository whose reproducibility status is otherwise clean.

**Resolved:** `uv` 0.12.22 was installed to the user profile (official installer,
no admin), then `uv add scipy`. `uv sync` reconciles clean and `uv lock --check`
reports no drift. `scipy 1.18.1` is installed alongside `numpy 2.5.1`,
`matplotlib 3.11.1`, `mpmath 1.4.1`, `pytest 9.1.1`.

---

## 9. What was actually built from this audit

`scipy` arrived; the substantive deliverable was not the dependency. It was
**two graduated canonical modules**, each equivalence-pinned against the donor it
was moved from, plus a recorded finding about the cost of graduation.

### 9.1 `src/lab_utils/ephemeris.py` (new, ~230 lines)

Moved faithfully from Mission 3's `mission_experiment.py`: `load_snapshot`,
`interp_snapshot`, `interp_snapshot_buggy`, `precession_j2000_to_mod`,
`third_body_accel`, `circular_ic`, plus `de441_snapshot_paths`.

**Why:** every mission in the chain was importing these from a *completed
mission's folder*. That is exactly the shared-ancestry arrangement that let one
DE441 interpolation defect corrupt four missions simultaneously. The fix for that
class is one canonical home plus pinning — not four copies.

**Pinned by** `tests/test_ephemeris.py`: bitwise equality with the donor for all
six functions on the real byte-pinned snapshots; sha256 of both snapshots;
orthonormality and determinant of the precession rotation; **and a regression
guard that the buggy interpolator still differs from the fixed one** (if they ever
agree, someone has silently reverted the fix that moved `a11` by 5x).

### 9.2 `src/lab_utils/estimation.py` (new, ~200 lines)

Collapses the six divergent `ols_slope` implementations found in §2.1.
`ols_slope` is bitwise-identical to the donors; `ols_slope_ci`, `bootstrap_ci`,
`cadence_uniformity`, and `dominant_period_days` (Lomb-Scargle) are new, with the
i.i.d. assumption stated once in one docstring rather than re-derived per mission.

Also constitution `LAB_CONSTITUTION.md` §12.3 backlog candidate #1.

### 9.3 Three real defects found while testing the new code

1. **`dominant_period_days` returned sidelobes.** Taking the top-N peaks by power
   reported one tone as 40.00 / 39.98 / 40.01 d — useless for the job the
   function exists for. Fixed by selecting distinct local maxima with a minimum
   log-period separation. Now recovers injected 365.25 / 182.6 / 27.3 d harmonics
   to <2% — the Exp-019 shape.
2. **`lombscargle` raised `MemoryError`.** scipy broadcasts to
   `(n_samples, n_freqs)`; 8572 x 20000 = 1.28 GiB, and `plan.md` records this host
   has **no pagefile**. Fixed by explicit chunking with a documented 32 MB cap.
3. **Graduating into `src/lab_utils/__init__.py` broke a committed provenance
   hash.** Exp 015's `results.json` pins `code_sha256` for `lab_utils/__init__.py`;
   the re-export invalidated it and `test_code_sha256_freshness_when_present`
   correctly failed. Rather than hand-edit the hash — which would forge
   provenance — the re-export was reverted and the new modules are consumed as
   submodules (`from lab_utils import ephemeris, estimation`).

**Finding 3 is the general lesson and is recorded for future graduations:** adding
to `lab_utils/__init__.py` is not a free action. It invalidates downstream
provenance records and must be paired with regenerating every dependent result
artifact, or it must be avoided. Prefer new submodules; graduate the re-exports in
a dedicated provenance-aware change.

### 9.4 Deliberately NOT done

- **The integrator was not replaced** (§4). Rejected on equivalence-pinning and
  provenance grounds.
- **Missions 1–3 were not rewired** to the new modules. The equivalence tests make
  that migration safe, but it is per-mission work belonging in a mission's own
  change, not in a cleanup session.
- **Exp 019's uniform-grid FFT was not migrated** to Lomb-Scargle. That is a
  *scientific* change altering a committed result and needs a pre-registered rule.

### 9.5 A fourth defect, found at commit time — **NOT fixed**, recorded for the provenance task

This one was not found while testing the new code; it was found while verifying
that the artifacts about to be committed could actually be checked.

`phase0_baseline_breathing.py` fingerprints itself with
`hashlib.sha256(p.read_bytes())` — over the **worktree** bytes. This repository sets
`core.autocrlf=true`, and `.gitattributes` marks only the raw `*.txt` snapshots as
`-text`, so `*.py` is stored LF in the blob and checked out CRLF.

| Bytes hashed | sha256[:16] |
|---|---|
| `71b2448` blob, LF as stored | `ac673b31cd9fd5be` |
| `71b2448` blob, CRLF — i.e. the worktree at run time | **`49d20cf9e62de55f`** |
| `phase0_1yr.json` → `provenance.code` | **`49d20cf9e62de55f`** |
| current worktree file | `646d68f7ecce4097` |
| `phase0_18.6yr.json` → `provenance.code` | **`646d68f7ecce4097`** |

Two conclusions, one reassuring and one not:

- **Reassuring:** the 1-yr artifact really was produced by the committed `71b2448`
  code. Its provenance is sound in substance.
- **Not reassuring:** the fingerprint is computed over bytes that **no longer exist in
  the repository**. Checking out the blob and hashing it yields `ac673b31…`, not the
  recorded `49d20cf9…`. Every result artifact carrying a `provenance.code` field is
  therefore *unverifiable from the repo alone*, and the failure mode is a contributor's
  `core.autocrlf` setting rather than anything about the science.

This is precisely the gap that audit §4's "reproducible or it did not happen"
attitude should catch, and it is a clean example of a claim that looks rigorous until
you try to reproduce it. **Not fixed here** — repairing it means re-running campaigns
to regenerate fingerprints, which would touch pre-registered evidence. Minimum fix:
normalise line endings before hashing (or hash the git blob), and add `*.py text eol=lf`
to `.gitattributes` so worktree and blob agree.

---

**Note on BLAS:** `plan.md` requires "BLAS threads pinned to 1". `numpy 2.5.1` links
its own OpenBLAS. The 18.6-yr campaign was launched with `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` exported.