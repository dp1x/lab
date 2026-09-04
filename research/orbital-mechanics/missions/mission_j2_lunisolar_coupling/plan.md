# mission_j2_lunisolar_coupling — Plan of Record

**Date**: 2026-09-03
**Status**: ACTIVE — execution in progress
**Mission card**: README.md
**Constitutional authority**: LAB_CONSTITUTION.md §13.1

---

## 0. Mission status (as of session open)

| Phase | Status | Output |
|-------|--------|--------|
| §1. FET evaluation | DONE | README.md §2 |
| §2. Literature research | DONE | See §3 below |
| §3. Reduced-model numerical experiments | IN PROGRESS | experiment.py |
| §4. Perturbative scaling experiment | PENDING | |
| §5. Full 18.6-yr discriminating test | PENDING | |
| §6. Adversarial test battery | PENDING | |
| §7. Final report + knowledge note | PENDING | |
| §8. Commit + push + post-push audit | PENDING | |

---

## 1. Mission question (verbatim from mission card)

The `mission_lunisolar_closure` mission found that the leading-order doubly-averaged quadrupole secular formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i` predicts the OSCULATING-element RAAN rate WRONG in sign at i_sso and i = 30° and off by 27× in magnitude at i = 90° (at h = 600 km, 18.6-yr arc). The mission investigates **J2 × Lunisolar coupling as a candidate mechanism**, with a pre-registered willingness to accept a NEGATIVE result.

---

## 2. FET verdict

PASSED on all 7 gates. See README.md §2.

---

## 3. Authoritative literature research

### 3.1 Sources consulted (via web_search; 4 distinct authoritative queries)

1. Murray & Dermott (1999), *Solar System Dynamics* — referenced in all sources
2. Brouwer & Clemence (1961), *Methods of Celestial Mechanics*
3. Kaula (1966), *Theory of Satellite Geodesy*
4. Kozai (1959), "The motion of a close earth satellite"
5. Cook (1962) — classic lunisolar third-body treatments
6. Vallado (2013), *Fundamentals of Astrodynamics and Applications*
7. Vagners, R. (no canonical reference found — "Vagners" was not a standard name; the relevant terms are "evection/variation forced mode" + inclination functions)
8. Wnuk, Kudryavtsev — second-order Lie transform treatments
9. Lyddane (1963) — Poincaré-equinoctial reformulation

### 3.2 Physics summary (FACT / INFERENCE / UNKNOWN)

**FACT (verified by ≥1 independent source):**

1. The J2 secular rates are: `Ω̇_J2 = -(3/2) n J2 (R_E/p)² cos i`, `ω̇_J2 = (3/4) n J2 (R_E/p)² (5 cos²i − 1)`. Source: Murray & Dermott §2, Kaula Ch. 4.

2. The first-order doubly-averaged quadrupole Lunisolar secular formula is `Ω̇_3b = (3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i` (Convention B / Murray & Dermott Eq. 7). Source: lab 018 (verified by audit-020 Track 1).

3. The DIRECT Lie-transform J2 × 3b cross term scales as `J2 × (n₃/n)²`. At LEO SSO: J2 ≈ 1e-3, (n₃/n)² for Moon ≈ (2.66e-6/1.08e-3)² ≈ 6×10⁻⁹. **Result: relative scaling 6×10⁻⁹ — too small to explain a 170× discrepancy by 8 orders of magnitude.** Source: standard Lie-transform perturbation theory.

4. The J2² second-order secular term scales as `J2² (R_E/a)⁴`. At h = 600 km: J2² ≈ 1.2e-6, (R_E/a)⁴ ≈ 0.7. **Result: relative scaling 8×10⁻⁷ — also too small.**

5. The KL (Kozai-Lidov) mechanism conserves H_z = √(1−e²) cos i and creates coupled e−i oscillations when the mutual inclination is in (39°, 141°). At LEO SSO (i = 97.79° ≈ 98°, in the KL "danger zone"), J2 detunes ω̇ so strongly that KL amplitudes remain modest for circular orbits. Source: KL classic theory.

6. The "forced-secular lunar nodal mode" appears in the Kaula expansion as the term `cos i × cos(Ω − Ω₃)` with magnitude ~ 3/4 n (n₃/n)² × (inclination function) — does NOT average to zero over one anomalistic cycle; does average over one nodal cycle IF the average is correctly taken over the slow oscillation. **Source: standard Kaula expansion.**

**INFERENCE (analytic; supported by ≥1 source):**

7. The dominant candidate explanation for the 18.6-yr discrepancy at i_sso is the **forced-secular lunar nodal mode** — the doubly-averaged formula assumes BOTH satellite AND third-body mean anomalies are averaged; in reality, the slow lunar nodal regression (18.6 yr) creates a residual secular term that survives the short-period averaging but is NOT fully captured by the doubly-averaged formula. The 18.6-yr arc DOES average over the nodal cycle, so this term should average to zero OVER 18.6 yr IF the only remaining contribution is the cos(Ω − Ω₃) modulation. But it does NOT — the i₃(t) itself varies during the nodal cycle, so the "secular" formula at constant i₃ ≠ "secular" formula at time-varying i₃(t).

8. The "J2 × Lunisolar" cross coupling (H1) can also occur INDIRECTLY via KL: the J2 secular rate is proportional to cos i / (1−e²)²; if the third-body secular Hamiltonian modulates e and i on long timescales, the J2 rate itself is modulated. But at LEO SSO with J2 detuning, the modulation amplitude is bounded.

**UNKNOWN (cannot resolve without numerical experiment):**

9. The exact magnitude of the J2 × Lunisolar cross term at h = 600 km, i = 97.79° (the prompt's central question).

10. Whether the forced-secular lunar nodal mode, the higher-order Lunisolar terms, or the mean-osculating bias dominates the discrepancy.

11. The magnitude of the i₃(t) secular-modulation contribution (the formula's i₃ is held at constant 28.584°, but DE441 has the actual time-varying i₃).

### 3.3 The mission's pre-registered prediction based on literature

**Before running any numerics**, the literature says:

- Direct J2 × 3b cross term scales as 10⁻⁹ — **predict H1 will be falsified** in its strongest form (coupling is not large enough to explain the 170× discrepancy)
- Forced-secular lunar nodal mode contributes a cos i × cos(Ω − Ω₃) term that, averaged over 18.6 yr, vanishes IF the average is taken correctly — but the residual after a finite-window fit may be substantial
- Higher-order Lunisolar terms (octupole) scale as (a/a₃)⁴ ≈ 10⁻⁷ of the quadrupole — too small alone but additive

The mission is therefore most likely to **falsify H1** and **identify the forced-secular lunar nodal mode as the dominant residual mechanism**. But the mission does NOT take this prediction as truth; it executes the investigation and lets the numerics decide.

---

## 4. Phase plan

### Phase A — Reduced-model experiments (compute: ~30 min)

#### A.1 Idealized circular-Moon geometry
- Build a synthetic Moon at fixed i₃, e₃ = 0, RAAN₃ = 0
- Propagate J2 + Sun + Moon for 1 yr at h = 600 km, i = 97.79°
- Compare Lunisolar residual against the case with real DE441
- Isolates: lunar eccentricity + inclination variation effects (H0b, H0c)

#### A.2 Synthetic vs real Moon scaling
- Compare full + J2-only differences for (a) synthetic circular Moon at 384400 km, (b) real DE441 Moon
- The difference, if any, is the lunar eccentricity / inclination-variation contribution
- This is a "perturbative difference" — small if the lunar e/i effects are small; large if they dominate

#### A.3 Mode isolation at 1 yr (continuity with Exp 018)
- Compute: 2-body, J2-only, Sun-only, Moon-only, Sun+Moon, J2+Sun+Moon
- Each at h = 600 km, i = 97.79°
- Compute R_J2x3b = (full) − J2-only − Sun-only − Moon-only + 2-body
- This is the cross-term residual, computed numerically

### Phase B — Perturbative scaling experiment (compute: ~1-2 hr)

#### B.1 λ_J2 × λ_3body sweep
- For λ_J2 ∈ {0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0} (7 values)
- For λ_3body ∈ {0, 0.5, 1.0, 2.0} (4 values)
- Total: 7 × 4 = 28 propagations at 1-yr arc, h = 600 km, i = 97.79°
- Each propagation streams ascending-node crossings

#### B.2 Two-dimensional polynomial fit
- The RAAN rate `f(λ_J2, λ_3body)` is fit to the model:
  `f = a₁₀ λ_J2 + a₀₁ λ_3body + a₁₁ λ_J2 λ_3body + a₂₀ λ_J2² + a₀₂ λ_3body²`
- The cross coefficient `a₁₁` is the perturbation-scaling discriminator
- If |a₁₁| > 3σ and sign is consistent with the H1 prediction, H1 is supported
- If |a₁₁| ≤ 3σ but |a₂₀| or |a₀₂| is significant, the residual is higher-order single-perturbation (H0)

### Phase C — Full 18.6-yr discriminating test (compute: ~3-4 hr)

Only executed if Phase A and Phase B give **consistent** results that warrant the 18.6-yr arc.

#### C.1 Force-mode decomposition at 18.6 yr
- All 6 modes at 3 inclinations × 18.6 yr = 18 propagations
- Streaming RK4 with the byte-pinned DE441

#### C.2 Cross-term residual extraction
- Compute R_J2x3b at all three inclinations
- Compare against analytical expectations from the Lie-transform theory

#### C.3 Forced-secular lunar nodal mode test
- Compute the analytical cos i × cos(Ω − Ω₃) term
- Compare against R_J2x3b at 18.6 yr

### Phase D — Adversarial battery (compute: ~30 min)

Tests the historical 015–020 bug classes:
- Sign-flip mutants (deliberate sign error in J2 acceleration; deliberate sign error in third-body formula)
- Frame transpose mutants (`_rot3` transpose bug, like audit-019 Track D)
- Scaling mutants (wrong power of a/a₃; wrong J2 scaling)
- Estimator mutants (raw OLS vs phase-locked)
- Subtraction cancellation mutants (force-mode isolation order swapped)

### Phase E — Final report + knowledge note

Per `mission_lunisolar_closure` template.

---

## 5. Compute budget

| Phase | Wall-clock (single-core) | Wall-clock (8 cores parallel) |
|-------|--------------------------|--------------------------------|
| A | ~30 min | ~10 min |
| B | ~1.5 hr (28 propagations × ~3 min each) | ~12 min |
| C | ~3-4 hr (18 propagations × ~12 min each at 18.6 yr) | ~30 min |
| D | ~30 min | ~10 min |
| **Total** | ~5-6 hr | **~1 hr** |

The 8-core parallel budget is well within the local-envelope constraint (`LAB_CONSTITUTION.md §6.3`: < 10 hr single-core, < 1 GB RAM). Streaming propagator at ~50-60 MB per worker.

---

## 6. Decision-rule execution

After Phase A + B + C, the mission:
1. Computes the cross-term residual R_J2x3b at 18.6 yr, 3 inclinations
2. Fits the perturbative scaling law and extracts `a₁₁` with uncertainty
3. Tests sign compatibility with H1 (cross term must contribute RETROGRADE at i_sso)
4. Tests the forced-secular lunar nodal mode as an alternative
5. Declares the final state:
   - **H1-SUPPORTED** if all five conditions in README §4.1 are met
   - **H1-FALSIFIED** if README §4.2 conditions hold
   - **H1-PARTIALLY-SUPPORTED** if some conditions hold
   - **UNRESOLVED** if the discrimination fails

---

## 7. Acceptance gate

The mission is COMPLETE when:
- All 7 mission-specific tests pass
- The full repo test suite passes (baseline 784 + new mission tests)
- Results.json is written with provenance (code hashes, snapshot hashes, decision-rule output)
- Figures are generated
- The scientific report is written in `localdocs/reports/mission-j2-lunisolar-coupling-2026-09-03.md`
- The knowledge note is written in `localdocs/knowledge/j2-lunisolar-coupling.md`
- A signed commit is pushed to origin/main with post-push verification

---

## 8. Risk and recovery

- If Phase A reveals an unexpected numerical artifact (e.g., the cross-term residual at 1 yr is dominated by subtraction cancellation rather than physics), the mission pauses and investigates the cancellation. This is consistent with the prompt's caution about subtraction cancellation.
- If Phase B's perturbative scaling fit gives ambiguous results (e.g., `a₁₁` is not 3σ significant), the mission falls back to Phase C with the alternative explanations explicitly tested.
- If compute budget exceeds 8 hr single-core or 1 hr parallel, the mission is paused and the work is checkpointed as an `AUTONOMOUS_HANDOFF_<date>.md` per the constitution's recovery-state rule.
- A 12-hr wall-clock pause is a hard stop per LAB_CONSTITUTION §4.3; if encountered, the mission writes an `AUTONOMOUS_HANDOFF_<date>.md` and stops.
