# plan.md — plan-of-record for mission_mean_element_leakage
**Date:** 2026-09-24. **Status:** ACTIVE. **Authority:** README.md (frozen protocol + pre-registered decision rule §4).

## Phases
0. **Baseline-breathing pre-check (1-yr, 3 inclinations × 5 modes).** Streaming RK4 recording osculating `a, e, i, Ω` at every ascending-node crossing; compute `D = ⟨Ω̇_J2(a_f,i_f) − Ω̇_J2(a_j,i_j)⟩` (the spurious rate the constant-baseline subtraction injects) and the measured residual `R` from the comparable campaign estimator; report `D/R` per inclination and the implied mean-element drifts `Δi, Δa`. Gate: propagator must reproduce the predecessor's `j2_only` rate; otherwise stop and debug.
   - 1-yr first (~10 min wall on 4 workers). 18.6-yr only if the 1-yr arm is informative (~40 min wall).
1. **Estimator ladder.** (A) angular-momentum-vector secular rate (no equator-crossing detection); (C) mean-element integration (upgrade `averaged_propagator.py` from the predecessor); bridge estimators FFT subtraction + multi-window joint fit (audit-020 Track 2 §8.4). Synthetic oracle gate: recover a known injected drift <1%.
2. **Mean-element re-test.** 6 modes × 3 inclinations, 1-yr and 18.6-yr, all estimators, with/without breathing correction; altitude ladder if H-baseline is confirmed.
3. **Real-geometry arm (conditional).** Re-derive lunisolar `A_k` with time-varying lunar `i3(t), e3(t)` from the pinned DE441 Moon — only if the joint artifact correction leaves >50% of `R`.
4. **Docs + staged commits.** Card → Phase 0 (+tests+results) → ladder (+tests) → re-test → report + knowledge note + roadmap row correction + AGENTS entry.

## Compute budget (4 workers, BLAS threads pinned to 1)
- Phase 0 1-yr: 15 propagations × ~1.5–3 min single-core ≈ 10 min wall.
- Phase 0 18.6-yr: 6 propagations × ~35 min single-core ≈ 40 min wall (2 modes × 3 inclinations).
- Estimator ladder + oracles: minutes (1D integrations, synthetic data).
- Phase 2 re-test: reuses Phase 0/1 outputs; extra Cowell runs only if the altitude ladder is triggered (~15 min).
- **Host constraint:** this machine currently has NO pagefile (20 GB commit ceiling). Keep `Pool(4)`, never 7; if MemoryError appears on small arrays, reduce to `Pool(2)` rather than editing test tolerance.

## Recovery
Every phase writes `results/phase*.json` with provenance (snapshot sha256 + code hashes). A fresh session resumes from the last committed JSON + the README decision rule. `R:` is scratch only; nothing durable lives there. No timing telemetry in any result file or document.

## Staged-commit contract
One commit per completed stage: (1) card + plan + Phase 0 script + tests, (2) Phase 0 results + verdict, (3) estimator ladder + oracle gate, (4) re-test + verdict, (5) report/knowledge/roadmap/AGENTS. Each commit is signed with the canonical identity and pushed only after the working tree is verified clean and the local suite is green.
