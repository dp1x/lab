# Autonomous Session Report — mission_lunisolar_closure

**Date**: 2026-09-03
**Session**: First clean scientific execution after LAB_CONSTITUTION adoption
**Branch**: main
**Final commit**: 62288cb (mission_lunisolar_closure complete: 18.6-yr DE441 arc refutes 018/020 at i_sso and i=30)
**Signed**: ✓ (RSA key 5774B47A005623ACD39DF7284EB7C30F884E8259; "Good signature from Dhanesh")

---

## 1. What was accomplished

### Mission executed end-to-end
- **Mission**: `mission_lunisolar_closure` (POST_ROADMAP_PROBE §13.1)
- **FET verdict**: PASS on all 7 gates
- **Status**: COMPLETE — `PARTIALLY-VERIFIED-WITH-OPEN-QUESTION`

### Scientific outcome (refutes prior claims)
At h = 600 km, 18.6-yr direct RK4 arc (one full lunar nodal cycle), J2+Sun+Moon minus J2-only:

| Inclination | Numerical (deg/day) | Corrected cf | Ratio | Sign match |
|---|---|---|---|---|
| i_sso (97.79°) | -2.29e-2 (retrograde) | +1.35e-4 (prograde) | **-170×** | **NO** |
| i=90 | +4.70e-3 (prograde) | +1.74e-4 (prograde) | +27× | YES |
| i=30 | -3.47e-4 (retrograde) | +4.55e-5 (prograde) | -7.6× | **NO** |

**The Exp 018/020 conclusion that the corrected doubly-averaged quadrupole formula gives the correct SIGN is REFUTED at i_sso and i=30 at the 18.6-yr arc.**

The Exp 019 polynomial-in-1/W extrapolation to W → ∞ is also REFUTED (the 18.6-yr numerical is opposite sign from the extrapolation).

### Implementation completed
- Mission moved from repo root to constitutional path `research/orbital-mechanics/missions/mission_lunisolar_closure/` per LAB_CONSTITUTION §2.3
- MANIFEST.json updated to drop R: path dependency (durable-storage violation fixed)
- `experiment.py`: streaming RK4 + 4 estimators + IAU-1976 precession; `code_hashes()` rewritten to walk up to repo root robustly; buggy `propagate_streaming` stub removed; `phase_locked_two_window` added with synthetic recovery test
- `run_parallel_campaign.py`: multiprocessing.Pool with 7 workers on 8 cores; 6 propagations in parallel; ~67 min wall-clock vs ~6 hr sequential. CPU utilization 80-100%.
- `run_focused_campaign.py`: legacy single-threaded orchestrator (still works, slower)
- `run_smoke.py`: pipeline validation script
- `make_figures.py`: 5 publication-quality figures from results.json (schema fixed)
- 19-yr DE441 Sun + Moon snapshots (sha256-pinned) preserved and verified

### Tests
- **13 new tests** (snapshot integrity, formula pin, synthetic oracle, force-level identity, phase-locked estimator with synthetic recovery + degenerate input, idealized bridge, headline decision rule)
- **784 total repo tests** (771 baseline + 13 new), all green
- Decision-rule test asserts the actual sign-disagreement finding (sign contract: numerical retrograde × cf prograde → product < 0 at i_sso and i=30; sign agreement at i=90)

### Documentation
- Mission card (`README.md`) per LAB_CONSTITUTION §2.3: question, FET verdict, hypothesis, frozen protocol, decision rule, force model, estimator theory, limitations, status, recommended next action
- Knowledge note (`localdocs/knowledge/lunisolar-closure-021.md`): full findings + supersession record
- Final scientific report (`localdocs/reports/mission-lunisolar-closure-2026-09-03.md`)
- AGENTS.md current-priority + roadmap.md updated

---

## 2. Evidence obtained

### Reproducibility
- Campaign re-run produces byte-identical Lunisolar contributions (verified twice)
- All 4 estimators (direct OLS, secant, harmonic regression, node-vector OLS) agree within ~4% at each inclination
- Synthetic oracle test: harmonic regression estimator recovers known secular to machine precision (bias 7×10⁻²⁰ deg/day)
- Force-level identity check: machine precision (max_diff = 0.0 km/s²) at 50 random states

### Adversarial-survival
- Pre-registered decision rule (50% magnitude tolerance, sign match required)
- Byte-pinned DE441 snapshots (sha256: `f2c4f048` for sun, `aee85099` for moon)
- Independent estimators (4 separate code paths) agree
- Phase-locked estimator tested against synthetic slow-harmonic oracle

---

## 3. What remains unresolved

The secular limit at W → ∞ remains UNRESOLVED, with the additional finding
that the leading-order corrected formula is NOT a valid asymptotic predictor
of the OSCULATING-element secular rate at LEO SSO under real DE441
ephemerides.

### Candidate explanations (open)
1. **J2 × Lunisolar coupling**: the leading-order secular formula treats J2
   and Lunisolar as independent. The cross-product term is unmodelled.
2. **Higher-order Lunisolar secular terms**: octupole `(a/a₃)⁴` omitted.
3. **Real-ephemeris vs doubly-averaged theory**: lunar orbit has e ≈ 0.05
   and i varies ±5° on the 18.6-yr nodal cycle; forced-mode secular
   contributions from evection (~31.8 d) and variation (~14.8 d) do not
   average to zero over one nodal cycle.
4. **Mean-vs-osculating bias**: corrected formula predicts MEAN element
   secular rate; numerical measures OSCULATING Ω at ascending-node
   crossings.

---

## 4. Repository changes

### Git state
- **HEAD**: 62288cb (signed)
- **Remote**: pushed to `origin/main`
- **Working tree**: clean
- **Commits ahead of origin**: 0 (after push)

### Files changed (35 staged)
- Renamed: 19 tracked files from `mission_lunisolar_closure/` (repo root) → `research/orbital-mechanics/missions/mission_lunisolar_closure/` (constitutional path)
- New: README.md (mission card), results/results.json, results/figures/fig{1..5}_*.png, run_parallel_campaign.py, run_smoke.py, reconstruct_results.py
- Modified: experiment.py (code_hashes fix, phase_locked_two_window addition), make_figures.py (schema fix), tests/test_mission_lunisolar_closure.py (decision rule), MANIFEST.json (no R: paths), .gitattributes (new paths)
- Updated: AGENTS.md (current priority), localdocs/roadmap.md (mission_lunisolar_closure entry)
- New knowledge: localdocs/knowledge/lunisolar-closure-021.md (full report)
- New report: localdocs/reports/mission-lunisolar-closure-2026-09-03.md

### Resource usage
- CPU utilization: 80-100% across 8 cores during campaign (user-requested utilization level achieved)
- Wall-clock for campaign: ~67 min (vs ~6 hr sequential without parallelization)
- Memory: streaming propagator uses ~50-60 MB per worker (no full-trajectory storage)
- Disk: no R: drive dependency; all artifacts in repo (under 5 MB total for snapshots + figures + JSON)

---

## 5. Test counts and status

- **Before session**: 771 tests
- **After session**: 784 tests (+13 from mission_lunisolar_closure)
- **All green**: confirmed via `pytest -q` (exit code 0, ~16 min wall-clock on commodity hardware)
- **Mission-specific**: 13 tests in
  `research/orbital-mechanics/missions/mission_lunisolar_closure/tests/test_mission_lunisolar_closure.py`

---

## 6. Final commit

```
62288cb (HEAD -> main, origin/main)
mission_lunisolar_closure complete: 18.6-yr DE441 arc refutes 018/020 at i_sso and i=30

First clean scientific execution after LAB_CONSTITUTION adoption
(commit 9d9a495). Mission at constitutional path
research/orbital-mechanics/missions/mission_lunisolar_closure/
...
[full commit body in git log]
```

---

## 7. Recommended next mission or follow-up question

Per `POST_ROADMAP_PROBE.md §13.2-13.5`:

1. **J2 × Lunisolar coupling derivation** (highest priority follow-on):
   - Hypothesis: the leading-order formula treats J2 + Lunisolar as
     independent, but the cross-product term at i_sso (where J2 cos(i) ≈ 0
     but its derivative is non-zero) may dominate
   - Goal: derive the J2 × Lunisolar secular coupling term; test whether
     it accounts for the i_sso sign disagreement
   - FET: gates 1, 2, 3, 4, 6, 7 pass; gate 5 (adversarial-survivable)
     requires 8-track audit
   - Estimated compute: 1-2 days; ~30 new tests

2. **Repeat-Ground-Track Targeting** (composition mission, POST_ROADMAP §13.2):
   - Composes Exp 008 ground tracks + Exp 012 SSO lock + Exp 015 launch
     windows into Landsat-7 16-day / 705 km reference orbit
   - Goal: end-to-end repeat-cycle targeting methodology
   - Estimated compute: 1-2 days; ~30 new tests

3. **Estimation Doctrine Graduation** (capability mission, §13.3):
   - Extract audit-020 lessons (OLS bias formula, harmonic regression,
     phase-locked windowing) into reusable `src/lab_utils/estimation.py`
   - Goal: codify the lab's estimation doctrine for future missions

The lab should NOT continue the lunisolar thread beyond this mission
without first investigating the J2 × Lunisolar coupling as the most likely
explanation for the sign disagreement. Otherwise the chain would continue
indefinitely without producing durable improvement.