# mission_j2_lunisolar_coupling — J2 × Lunisolar Coupling as Candidate Explanation for the 18.6-yr DE441 Residual

**Mission type:** Discrepancy mission (per LAB_CONSTITUTION.md §2.5)
**Status (2026-09-03):** active; mission card written; literature review complete; numerical investigation in progress
**Constitutional authority:** `LAB_CONSTITUTION.md §13.1`
**Roadmap position:** direct follow-on to `mission_lunisolar_closure` (recommended next action in the 2026-09-03 report)

---

## 1. Mission question

The `mission_lunisolar_closure` mission (2026-09-03) executed an 18.6-yr direct RK4 arc with byte-pinned DE441 Sun + Moon snapshots at h = 600 km, i ∈ {97.79° (i_sso), 90°, 30°}, and found that the leading-order doubly-averaged quadrupole secular formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i` (the "corrected formula") predicts the OSCULATING-element RAAN rate WRONG in sign at i_sso and i = 30° and off by 27× in magnitude at i = 90°.

The mission's open question — left by `mission_lunisolar_closure` itself — is: **what physical mechanism is responsible for the discrepancy, and is "J2 × Lunisolar coupling" a competitive explanation?**

This mission investigates J2 × Lunisolar coupling as a candidate mechanism, **with a pre-registered willingness to accept a NEGATIVE result** if the coupling turns out to be too small to explain the discrepancy. The prompt explicitly invited this: "A negative result is completely acceptable and should be preferred over a forced explanation."

---

## 2. FET verdict (passed before mission execution)

Per `LAB_CONSTITUTION.md §9`, this mission passes all 7 gates:

| Gate | Verdict | Notes |
|------|---------|-------|
| 1. Reasoning/compute ratio | PASS | High reasoning content: literature review (Brouwer, Cook, Kaula, Lie-transform theory) + controlled force-model decomposition + perturbative scaling experiment + 18.6-yr validation. Compute: streaming RK4 with mode isolation; < 5 hr single-core. |
| 2. Independent validation | PASS | Two independent routes: (a) analytical derivation of expected J2 × 3b scaling from first-principles Hamiltonian theory; (b) computational residual isolation via mode subtraction + perturbative scaling law fit. Byte-pinned DE441 + synthetic oracle. |
| 3. Durable scientific knowledge | PASS | Bounds the actual magnitude of J2 × Lunisolar coupling at LEO SSO; falsifies or quantifies the leading hypothesis from `mission_lunisolar_closure`; identifies the most likely alternative explanation (forced-secular lunar nodal mode) if J2 coupling is too small. |
| 4. Hypothesis-distinguishing | PASS | Multiple competing hypotheses: (H1) J2 × Lunisolar coupling is large enough to explain the discrepancy; (H0) coupling is negligible and the discrepancy comes from something else (forced-secular lunar nodal mode, octupole, evection/variation, mean/osculating bias). The mission declares a pre-registered decision rule to discriminate. |
| 5. Adversarial-survivable | PASS | Pre-registered decision rule with falsifiable quantitative criteria; the prompt explicitly invited the test design that catches 015–020 bug classes (sign errors, scaling errors, frame mismatches, estimator artifacts). |
| 6. Capability-advancing | PASS | Reuses streaming RK4 + IAU-1976 precession + ascending-node detection from `mission_lunisolar_closure`; introduces perturbative-scaling experiment design as a reusable methodology. |
| 7. Deterministic on modest resources | PASS | Streaming propagator ~50 MB per worker; ~2 MB repo-resident data; < 5 hr single-core for the full mission (smaller for the discriminating tests). |

ROI tuple (private): (information_gain: HIGH — directly addresses the residual; capability_advance: MEDIUM — perturbative scaling is reusable; validation_strength: E4/E5 — byte-pinned DE441 + analytical theory; adversarial_survival: pre-registered decision rules; compute_cost: ~5 hr single-core; attention_cost: ~one human interaction at final report).

**Mission selected**: highest-ROI candidate per `mission_lunisolar_closure`'s "recommended next action" section. No FET-passing competitor with higher information_gain currently exists in the queue.

---

## 3. Hypotheses to discriminate

The mission pre-registers the following hypotheses (all evaluated against the byte-pinned 18.6-yr DE441 + J2 numerical Lunisolar contribution to RAAN at h = 600 km, i_sso = 97.79°):

### H1 (J2 × Lunisolar coupling IS the dominant explanation)
A secular or long-period cross-coupling term proportional to J2 × (μ₃/μ_E) (a/a₃)³ × (inclination function) is large enough, and of the correct sign and inclination dependence, to materially explain the discrepancy observed at the 18.6-yr arc.

**Pre-registered prediction**: If H1 holds, the residual after subtracting (J2_only + Sun_only + Moon_only − 2-body) from the full propagation should be:
- Order of magnitude ~ J2 × 3b_signal (the cross term should be a sizable fraction of the J2 contribution, not orders of magnitude smaller)
- Survives numerical convergence (dt ladder)
- Scales as λ_J2 × λ_3body under perturbative scaling
- Sign compatible with the observed retrograde at i_sso

### H0 (J2 × Lunisolar coupling is NOT the dominant explanation)
The cross-coupling term is present but too small (or wrong sign) to explain the discrepancy; the actual mechanism is something else. Strongest candidates per the literature and prior audits:
- (a) **Forced-secular lunar nodal mode**: the "cos i × cos(Ω − Ω₃)" term that survives averaging over one anomalistic cycle but not the full 18.6-yr nodal cycle (the doubly-averaged formula assumes BOTH averages)
- (b) **Higher-order lunar terms**: octupole, evection/variation-induced secular contributions
- (c) **Real lunar eccentricity / inclination**: the doubly-averaged theory assumes e₃ = 0, i₃ constant
- (d) **Mean-vs-osculating bias**: the corrected formula predicts MEAN-element secular; the numerical measures OSCULATING Ω at ascending-node crossings

### H-uncertainty (discrimination fails)
The mission fails to discriminate between H1 and H0 within the budget; the residual structure is ambiguous.

---

## 4. Pre-registered decision rule

The mission defines its quantitative decision rule BEFORE executing any final-scale numerical work:

### 4.1 Discriminator: perturbative scaling fit

For the mission to **support H1**:
- (a) The combined Lunisolar residual (full − J2-only − Sun-only − Moon-only + 2-body) at h = 600 km, i_sso, 18.6-yr arc must have magnitude > 10⁻⁴ deg/day (comparable to the J2-only contribution, NOT negligible compared to it)
- (b) The residual must scale with the product λ_J2 × λ_3body in the perturbative scaling experiment, with the fit coefficient on the cross term exceeding the residuals on the linear terms by > 3σ
- (c) The residual must have sign compatible with the observed retrograde at i_sso (i.e., the J2 × coupling at i_sso must contribute RETROGRADE to the total RAAN drift)
- (d) The residual at i = 90° (where J2 cos i = 0 — the J2-clean test) must be MUCH smaller than at i_sso, demonstrating the J2-modulated structure
- (e) The residual must survive the numerical convergence ladder at dt ∈ {60, 30, 15, 7.5} s with stable sign and magnitude (NOT an integration artifact)

**All five conditions must be met for H1 support.** Failure of any single condition falsifies H1 and the mission declares the result.

### 4.2 Discriminator: magnitude bound

For the mission to **falsify H1**:
- (a) The residual is below numerical noise OR has the wrong sign OR
- (b) The residual scales as λ_J2² alone (or λ_3body² alone), NOT as the cross product λ_J2 × λ_3body — indicating it is NOT a coupling term but a second-order single-perturbation term

### 4.3 Discriminator: alternative identification

If H1 is falsified, the mission identifies which alternative explanation (H0a–H0d) is competitive by:
- (a) Computing the analytical forced-secular lunar nodal mode term and comparing to the residual (test H0a)
- (b) Computing the octupole Lunisolar secular formula and comparing (test H0b)
- (c) Testing whether the residual structure survives in an idealized circular-Moon model (test H0c)
- (d) Testing whether the residual correlates with the mean/osculating bias estimators (test H0d)

The mission's final scientific statement is one of: **H1-SUPPORTED**, **H1-FALSIFIED**, **H1-PARTIALLY-SUPPORTED**, or **UNRESOLVED** (with reason).

---

## 5. Scientific protocol (frozen)

### 5.1 Numerical setup
| Item | Value | Justification |
|---|---|---|
| Frame | ECI mean-of-date; Sun/Moon rotated from ICRF/J2000 via FIXED IAU-1976 precession | Continuity with `mission_lunisolar_closure` |
| Integrator | RK4 fixed-step, dt = 60 s | 018/019/020 verified; design order p ≈ 4.5 |
| Ephemeris | Byte-pinned DE441 Sun + Moon, 2026-01-01 → 2045-01-01, daily, ICRF/TDB | Reuse from `mission_lunisolar_closure` |
| Inclinations | i ∈ {97.79° (i_sso), 90°, 30°} | Continuity with prior; i=90° is the J2-clean control |
| Altitude | h = 600 km | Lab SSO reference |
| Modes (force isolation) | 2-body only, J2-only, Sun-only, Moon-only, Sun+Moon, J2+Sun+Moon | 6 modes × 3 inclinations × (headline + perturbative scaling variants) |
| Output cadence | Ascending-node crossings + subsampled node-vector samples | Streaming, no full-trajectory storage |

### 5.2 Perturbative scaling experiment
A RESEARCH-ONLY force model with multipliers `λ_J2 ∈ {0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0}` and `λ_3body ∈ {0, 0.5, 1.0, 2.0}` is added to the propagation; the Lunisolar RAAN rate is measured at each (λ_J2, λ_3body) pair. A two-dimensional polynomial response `f(λ_J2, λ_3body) = a₁₀ λ_J2 + a₀₁ λ_3body + a₁₁ λ_J2 λ_3body + a₂₀ λ_J2² + a₀₂ λ_3body² + ...` is fit by least squares. The fitted coefficient `a₁₁` on the cross term is the perturbation-scaling discriminator.

This is the **critical diagnostic** the prompt specified. If `a₁₁` is statistically significant (exceeds its uncertainty by > 3σ), there is a genuine cross term. If `a₁₁` is consistent with zero AND `a₂₀` and `a₀₂` are significant, the residual is a higher-order single-perturbation term (e.g., J2²).

### 5.3 Reduced-model experiments (preliminary, before full 18.6-yr arc)

Per the prompt's guidance ("develop smaller idealized/reduced-model experiments to answer derivation questions efficiently, then run the full DE441 case only for the final discriminating tests"), the mission executes:

1. **Idealized circular-Moon geometry** (fixed i₃, e₃ = 0, fixed RAAN3 = 0) — clean test of J2 × Lunisolar coupling in the absence of all periodic effects; if H1 holds, the residual here should be CLEANLY a J2 × Lunisolar signal.
2. **Real-DE441 / synthetic DE441 comparison** — same mode set, but with synthetic circular orbit replacing the Moon; isolates the lunar-eccentricity contribution.
3. **Multiplicative scaling** at 90-d arc (NOT 18.6-yr) — fast iteration on the perturbative-scaling fit before committing to 18.6-yr compute.

Only after the reduced-model experiments show consistent discrimination does the mission commit to the full 18.6-yr arc.

---

## 6. Force-model decomposition (the non-additive residual)

The central numerical observable is:

    R_J2x3b(Ω) = Ω̇_full - Ω̇_J2 - Ω̇_Sun - Ω̇_Moon + Ω̇_2body

where each `Ω̇_mode` is the secular RAAN rate from a propagation in that mode only.

If H1 holds, R_J2x3b is a genuine coupling signal; if H0 holds, R_J2x3b is a combination of:
- Subtraction cancellation (each mode contributes ~10⁻³ deg/day; differencing them may amplify relative error)
- Mean-osculating bias differences between single-mode and combined-mode propagations
- Higher-order single-perturbation terms (J2², third-body octupole)

The mission quantifies each potential contribution separately.

---

## 7. Limitations (declared upfront)

1. The corrected doubly-averaged quadrupole formula is the first-order secular; it omits J2² and the cross term. We do NOT have an independent closed-form for J2 × Lunisolar coupling to compare to the numerical — we must derive it ourselves and check that it matches the observed residual structure.
2. The 18.6-yr arc measures OSCULATING Ω at ascending-node crossings; the corrected formula is for MEAN elements. The mission explicitly tests for mean-osculating bias.
3. The lunar ephemeris is byte-pinned DE441 daily; interpolation is linear. Evection (27.55 d) and variation (14.77 d) at the daily-cadence resolution may be aliased — though the mission's purpose is the secular rate, not the periodic terms.
4. SRP, drag, and higher-order geopotential are omitted as in `mission_lunisolar_closure`.

---

## 8. Implementation summary

- **`experiment.py`** — streaming RK4 propagator with mode isolation + perturbative scaling multipliers + 4 estimators + IAU-1976 precession + ascending-node detection (largely reused from `mission_lunisolar_closure`)
- **`run_perturbative_scaling.py`** — orchestrator for the (λ_J2, λ_3body) sweep
- **`run_reduced_models.py`** — idealized circular-Moon and synthetic-Moon propagations
- **`make_figures.py`** — publication-quality figures from results.json
- **`tests/`** — 12+ tests covering:
  - Force-level identity at machine precision
  - Synthetic oracle on cross-term estimator
  - Perturbative-scaling fit invariants
  - Frame-convention checks (no `_rot3` transpose)
  - Sign-convention tests
  - Decision rule post-conditions
  - Bug-class mutants from 015–020 history

---

## 9. Status at session open (2026-09-03)

- Mission card (this README): WRITTEN
- FET evaluation: PASSED (this document)
- Literature review: COMPLETE (logged in plan.md)
- Plan-of-record: in `plan.md`
- Reduced-model experiments: in progress
- Full 18.6-yr discriminating test: pending reduced-model completion

---

## 10. Recommended next action (after this mission)

If **H1-SUPPORTED**: graduate the perturbative-scaling estimator as a reusable methodology in `src/lab_utils/`; the J2 × Lunisolar coupling term becomes a new closed-form secular component to add to `src/lab_utils/orbits.py`.

If **H1-FALSIFIED**: the mission's contribution is to **bound the actual J2 × Lunisolar coupling** and identify the most likely alternative mechanism. The forced-secular lunar nodal mode becomes the next priority candidate.

If **H1-PARTIALLY-SUPPORTED** or **UNRESOLVED**: spawn a follow-on mission investigating the forced-secular lunar nodal mode explicitly, with the J2 × Lunisolar coupling bound as a known minor contribution.

---

## 11. Reference / supersession record

- **`mission_lunisolar_closure` (2026-09-03)**: established the 18.6-yr finding; PARTIALLY-VERIFIED-WITH-OPEN-QUESTION; recommended "spawn a follow-on investigation: derive the J2 × Lunisolar secular cross-coupling term and test whether it accounts for the sign disagreement at i_sso." This mission is that follow-on.
- **Exp 017–020 + audits**: the corrected formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i` (Convention B / Murray & Dermott) is the leading-order result; this mission tests whether a NEXT-ORDER term (J2 × Lunisolar) is needed.
- **Audit-020 Track 1 (2026-08-31)**: established the Convention-B sign convention (the lab's `+` sign is correct and matches data at all three inclinations at the 1-yr arc). This mission operates on the same Convention-B formula.

---
