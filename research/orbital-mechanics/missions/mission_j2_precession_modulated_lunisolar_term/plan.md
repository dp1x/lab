# plan.md — plan-of-record for mission_j2_precession_modulated_lunisolar_term
**Date:** 2026-09-05. **Status:** COMPLETE (2026-09-22) — verdict FALSIFIED (README §10). **Authority:** README.md (frozen protocol + decision rule).
**Phase outcomes:** 1-3 DONE (theory, reduced models, short-arc campaigns); 4 DONE for dt ladder (4-estimator/4-phase gates not executed — documented limitation); 5 partially DONE (octupole bound + interp-bug quantification; mutant battery not executed); 6 DONE (18.6-yr reproduction); 7 adversarial review NOT completed (subagent provider failures — recorded); 8 DONE (card/report/knowledge/roadmap); 9 DONE (suite, signed commit, push, post-push verify). Control-semantics remediation: passive rotating-frame attempt 2026-09-12 reverted to active open-loop controls with bias-invariance guard; unfrozen rows bit-identical.

## Phases
1. **Theory (lead DERIVATION.md + theory_modulation.py + independent subagent tracks).** DONE (lead); subagents async; reconciliation before final report.
2. **Reduced analytic/numeric (seconds–minutes).** Synthetic OLS-bias oracles, quadrature amplitudes, averaged-element propagator vs Cowell spot-checks, detection-limit quantification. Gate: synthetic bias formula <1% vs oracle; averaged path reproduces H-mod sign.
3. **Short-arc Cowell (90-d scaling + prescribed-precession + frozen-plane, 365-d headline at 3–5 inclinations × 2–3 altitudes).** Reuses corrected 6-mode isolation; parallelized; streaming (no full-trajectory storage). Gate: reproduce prior a11 sign + R order before interpreting; verify bug-class regression green.
4. **Robustness (dt ladder, 4 estimators, 4-phase, cadence decimation, small-scale symmetry).** Gate: §6(e) tolerances or declare failure.
5. **Alternatives (octupole bound, synthetic-vs-real Moon, solar-lunar cross, frame/sign mutants, mean/osculating synthetic).** Gate: §6(f).
6. **Final 18.6-yr validation (6 propagations: 2 modes × 3 i, dt=60 s, daily DE441).** Only after 2–5 pass/fail is recorded. Uses theory `w(omega·W)` prediction (no refit). ~30–70 min parallel.
7. **Adversarial review (dedicated hostile subagent + mutant battery).** Must attempt falsification; blind spots recorded.
8. **Report + knowledge + lifecycle (roadmap update, supersession explicit, no timing telemetry in docs).**
9. **Full suite green → pre-commit audit → remote-tip recheck → signed commit → push → post-push verification.**

## Compute budget (single-core equivalents)
- Reduced/synthetic: <5 min. Amplitudes quadrature cached per (h,i).
- 90-d scaling (28+symmetry ≈ 40 propagations × ~40 s) ≈ 25 min single-core, ~5 min on 7 workers.
- 365-d headline (6 modes × 5 i × 1 alt ≈ 30 × ~3 min) + altitudes (6 modes × 3 i × 2 extra alt ≈ 36 × ~3 min) ≈ 3 hr single-core, ~25 min parallel. Trimmed to headline-first; altitude second.
- Robustness/phase/cadence deltas reuse headline states where possible; extra ~1 hr single-core.
- Final 18.6-yr (6 × ~35 min single-core) ≈ 3.5 hr single-core, ~35 min parallel.
- Total ≈ 8 hr single-core, ~1.2 hr wall parallel — within 10-hr envelope. Checkpoint after each phase (results JSON).

## Recovery
Each run script writes `results/<phase>.json` with provenance (code sha, snapshot sha, constants, estimator config). Fresh session resumes from last committed JSON + README decision rule; R: never required.
