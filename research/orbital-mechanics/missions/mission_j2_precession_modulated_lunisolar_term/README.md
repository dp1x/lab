# mission_j2_precession_modulated_lunisolar_term — Can J2-driven plane precession quantitatively explain the observed J2×lunisolar interaction without empirical fitting?

**Mission type:** Discrepancy mission (LAB_CONSTITUTION.md §2.5), follow-on to `mission_j2_lunisolar_coupling` and `mission_lunisolar_closure`.
**Status:** COMPLETE (2026-09-22) — **H-mod FALSIFIED** per the pre-registered decision rule (§6); verdict + evidence in §10, report in `localdocs/reports/mission-j2-precession-modulated-lunisolar-term-2026-09-22.md`, knowledge note in `localdocs/knowledge/j2-precession-modulated-lunisolar-term.md`.
**Constitutional authority:** LAB_CONSTITUTION.md §§2.3, 9, 12–13.
**Date established:** 2026-09-05. **Date completed:** 2026-09-22.
**Budget:** ≤10 hr single-core, ≤1 GB RAM, local commodity compute only. Reduced models seconds–minutes; short-arc RK4 sweeps ~20 min parallel; final 18.6-yr validation ~30–70 min parallel (6 propagations).

---

## 1. Question

Can the strong non-additive J2×lunisolar interaction reported by `mission_j2_lunisolar_coupling` — perturbative cross coefficient `a11 ≈ -7.85e-4 deg/day` (SNR ≈ 6.9) at 90 d, and 1-yr force-mode residual `R_J2x3b` representing 74–92% of the combined lunisolar contribution — be derived from a controlled perturbation/averaging treatment as an explicit J2-precession-modulated lunisolar term that quantitatively predicts the measured interaction **without fitting an empirical coefficient**, and what is its regime of validity?

`mission_lunisolar_closure` (18.6-yr DE441, h=600 km): numerical lunisolar (full−J2) −2.29e-2 deg/day at i_sso vs corrected doubly-averaged quadrupole +1.35e-4 (sign wrong, −170×); +4.70e-3 at i=90 vs +1.74e-4 (+27×); −3.47e-4 at i=30 vs +4.55e-5 (−7.6×). `mission_j2_lunisolar_coupling` (1-yr, corrected isolation): `R_J2x3b = (full−J2) − (Sun+Moon−2body)` = −1.16e-3 (i_sso), −4.21e-3 (i=90), −2.42e-4 (i=30) deg/day. Direct second-order Lie cross term scales ~J2·(n3/n)² ~1e-9, far too small. The phrase "J2-precession-modulated lunisolar coupling" is a **hypothesis**, not a result. This mission derives it falsifiably or refutes it.

---

## 2. FET verdict (passed before execution, all 7 gates binary)

| Gate | Verdict | Reason |
|---|---|---|
| 1. Reasoning/compute ratio | PASS | Constrained by Kaula/averaging derivation, resonance/estimator theory, adversarial review; compute is cheap reduced models + short-arc sweeps + one 18.6-yr validation (~4 hr single-core equivalent). |
| 2. Independent validation | PASS | Byte-pinned DE441 Sun+Moon (sha256 f2c4f048… sun, aee85099… moon), closed-form theory, synthetic oracles, force-level identity, independent derivation + independent code path. |
| 3. Durable knowledge | PASS | Bounds/augments or refutes the lab secular lunisolar canon at LEO; reusable prescribed-precession + scaling methodology; cited by station-keeping/estimation work. |
| 4. Hypothesis-distinguishing | PASS | Distinguishes direct cross-order vs precession-modulation of periodic term vs estimator aliasing via inclination/altitude/frequency/window signatures (see §4). |
| 5. Adversarial-survivable | PASS | Pre-registered quantitative decision rule (§6); hostile mutants; phase/window/cadence tests; audit-grade bug guards. |
| 6. Capability-advancing | PASS | Prescribed-precession control, perturbative-scaling hardening, synthetic-estimator harness as reusable methodology. |
| 7. Deterministic modest resources | PASS | Streaming RK4 (no full-trajectory storage), fixed dt=60 s, daily DE441, ≤10 hr single-core, ≤1 GB RAM, local only. |

**Selection rationale:** Highest-ROI FET-passing candidate addressing the open question left by two prior missions. Backlog alternatives (Estimation Doctrine Graduation, Repeat-Ground-Track, NRHO) have lower immediate information gain for this discrepancy. No numeric FET score is published (constitution §9.2). ROI tuple is a private ranking input, not a published score.

---

## 3. Hypotheses (pre-registered)

- **H-mod (J2-precession modulation):** The dominant interaction is the adiabatic modulation of the singly-averaged lunisolar node response by the J2-driven plane precession `Omega(t)=Omega0+Omega_dot_J2·t`. The instantaneous lunisolar node rate `Omega_dot_3b(Omega(t), Sun/Moon geometry)` contains long-period harmonics in `Omega−Omega3`, `2(Omega−Omega3)`, `Omega−lambda_sun`, etc. Over a finite window W the OLS secular estimator returns the true doubly-averaged secular **plus** the finite-window bias of those harmonics. For a harmonic `A·cos(omega·t+phi)` with `omega` a beat involving `Omega_dot_J2`, the OLS slope bias is `≈ −A·omega·sin(phi+omega·W/2)·sinc(omega·W/2)`-type (exact formula in derivation doc; scales as amplitude × frequency for `|omega|·W << 1`). Since `A ∝ 3body strength` and `omega ∝ J2 precession` (or beat with solar/lunar rates), the apparent extra slope scales as `lambda_J2 × lambda_3body` with no true second-order secular term. At sun-synchronous resonance (`Omega_dot_J2 ≈ n_sun`) the solar beat `omega → 0`, so the bias persists even at 18.6-yr and masquerades as secular. The `W→∞` asymptotic secular remains the doubly-averaged quadrupole; the modulation is long-period, not a new secular canon term, unless a zero-frequency beat exists.
- **H-x (direct cross-order):** A genuine second-order secular term ∝ J2·(mu3/mu)·(a/a3)³ exists with magnitude comparable to the observed R. Requires Lie-transform magnitude ~1e-9 relative, so predicts |R| ≪ observed. H-x is the null against which H-mod is tested.
- **H-alias (pure estimator/mean-osculating artifact without J2-frequency dependence):** The residual is OLS aliasing of evection/variation/annual terms whose frequencies do **not** involve `Omega_dot_J2`. Predicts R independent of prescribed precession rate and surviving frozen-plane control. Distinguished from H-mod by the prescribed-precession and frozen-plane tests.
- **H-other:** Octupole, lunar e/i variation, solar-lunar cross, or frame/time-scale error dominates. Each has distinct altitude/distance-power/inclination/frequency signatures (§5.4).

---

## 4. Candidate theory (falsifiable, no empirical coefficient)

Full derivation in `DERIVATION.md` (lead) + independent subagent derivation (returned text, reconciled in report). Summary:

- Start from `R = R_J2 + R_3b`, Lagrange `dOmega/dt = [n·a²·sqrt(1−e²)·sin i]⁻¹ · ∂R/∂i`.
- Singly-averaged (over satellite M) lunisolar rate: `Omega_dot_3b(t) = S0(i,i3,e3,a/a3) + Σ_k A_k(i,i3,…)·cos(omega_k·t + phi_k)`, where `S0` is the doubly-averaged secular `(3/8)n(mu3/mu)(a/a3)³·sin2(i−i3)/sin i` (Convention B, prograde at tested inclinations), and `A_k`, `omega_k`, `phi_k` are the long-period amplitudes/frequencies/phases. Leading long-period arguments involve `Omega(t)−Omega3(t)`, `2(Omega(t)−Omega3(t))`, `Omega(t)−lambda_sun(t)`, lunar nodal `Omega3(t)` (18.6-yr), annual solar, monthly lunar terms. With `Omega(t)=Omega0+Omega_dot_J2·t`, each `omega_k = m·Omega_dot_J2 + n·n_sun + p·n_moon + q·Omega_dot_3 + …` for small integers.
- Finite-window OLS slope over W: `s_hat = S0 + S_J2 + Σ_k B_k`, `B_k = A_k·omega_k·W·w(omega_k·W,phi_k)` with `|B_k| ≤ |A_k·omega_k|` for short windows and `B_k → O(A_k/(W²·omega_k))` or `O(A_k/W)` for fast harmonics at long W, and `B_k → quasi-constant` for resonant `omega_k→0`. Full `w()` derived in DERIVATION.md and tested against synthetic oracles. No fitted constant enters `S0`, `A_k`, or `omega_k`; all follow from `n, a/a3, mu3/mu, J2(R/p)², i/i3`, and ephemeris frequencies.
- Perturbation order: first-order in 3body for `A_k` and `S0`; J2 enters only through `Omega(t)` (adiabatic/resummation, **not** a conventional second-order Lie cross term). Taylor expansion for small `Omega_dot_J2·W` gives apparent mixed scaling `∝ lambda_J2·lambda_3body`, explaining the Phase-B `a11` without a true `J2·3b` secular term.
- Predicted scalings: `S0 ∝ n·(a/a3)³`; `A_k ∝ n·(a/a3)³·F_k(i,i3)`; `omega_k ∝ Omega_dot_J2 ∝ n·J2·(R/p)²·cos i` plus solar/lunar rates; `B_k ∝ A_k·omega_k ∝ n²·J2·(mu3/mu)·(R/p)²·(a/a3)³` for short-window non-resonant beats. Altitude scaling is therefore **steeper** than `S0` alone (extra `n·J2·(R/p)²` factor). Inclination dependence follows `F_k(i,i3)·cos i` (non-resonant) plus resonant solar term peaking where `Omega_dot_J2 ≈ n_sun` (near SSO). Sign depends on `phi_k` (phase) and `d(Omega_dot_3b)/dOmega` sign; the mission predicts the **1-yr R sign pattern** from the derived `F_k` + 2026 phase (lunar node near descending, solar longitude at epoch) and tests it — sign agreement at all three inclinations is required (see §6).
- Asymptotic claim: H-mod predicts the `W→∞` limit is `S0` (plus small true second-order ~1e-9 relative), **not** a new secular canon term. A persistent `W→∞` offset after full beat averaging would falsify H-mod in its pure form and require H-x or H-other. Sun-sync resonance is the exception where the solar beat period →∞ and the bias remains quasi-secular on decadal windows; this is still long-period, not asymptotic secular, and must be labeled as such.

---

## 5. Frozen protocol

| Item | Value | Justification |
|---|---|---|
| Frame / time | ECI mean-of-date; Sun/Moon ICRF/J2000→MOD via FIXED IAU-1976 precession; TDB daily DE441 | Continuity with closure/coupling missions; convention pinned by regression |
| Integrator | RK4 fixed-step dt=60 s (headline); dt∈{120,60,30} ladder | Verified design order p≈4.5; ladder tests integration artifact |
| Modes | kepler_only, j2_only, sun_only, moon_only, sun_moon, sun_moon_j2 (6 modes). Corrected isolation `use_j2 = mode in ("j2_only","sun_moon_j2")` guarded by regression | Remediation from coupling mission; bug-class regression required |
| Third-body accel | Geocentric direct+indirect `mu3[(r−r3)/\|r−r3\|³ − r3/\|r3\|³]`; force-level identity at 50 random states to machine precision | Prior identity standard |
| Ephemeris | Byte-pinned 19-yr DE441 Sun+Moon 2026-01-01→2045-01-01 daily (sha256 sun f2c4f048…, moon aee85099…); linear interpolation; synthetic circular-Moon control (fixed i3, e3=0) | Provenance verified at startup; synthetic isolates lunar e/i |
| Epoch / ICs | T0 = JD_J2000+820476800 s; circular LEO, J2-canonical mean elements at each (h,i); single anomalistic phase + 4-phase ensemble for headline | Continuity; phase test per estimator doctrine |
| Altitudes | h∈{500,600,800} km (LEO-like; 600 headline) | Altitude power-law discriminator |
| Inclinations | i∈{97.79 (SSO), 90, 30} headline + i∈{60,0} discriminators (prograde/equatorial structure) | Sign-reversal + inclination-function test |
| Arcs | Reduced/analytic seconds; 90-d perturbative + prescribed-precession sweeps; 365-d force-mode headline; 18.6-yr final validation (2 modes×3 i = 6 propagations) | Cheap-first, full-arc only for final discrimination |
| Scaling | lambda_J2∈{0,0.25,0.5,0.75,1.0,1.5,2.0}, lambda_3body∈{0,0.5,1.0,2.0} (90-d, i_sso) + small-scale symmetry set {±0.25,±0.5} for detection-limit test | Repeats Mission-2 methodology with stability + symmetry checks |
| Prescribed-precession control | (a) Frozen-plane: counter-rotate to hold Omega fixed (kills modulation, keeps forces); (b) Kinematic precession without J2 force: rigid z-rotation at prescribed `Omega_dot_p ∈ {−2,−1,−0.5,0,+0.5,+1,+2}×Omega_dot_J2(600,i_sso)` with 3body only | Direct frequency-dependence test; H-mod must vanish in (a) and track `Omega_dot_p` in (b) |
| Estimators | Direct OLS, secant, harmonic regression (f), node-vector OLS (n) at ascending-node crossings + phase-locked 2-window cross-check; output cadence every node (headline) vs every 4th node | Four-estimator agreement + cadence invariance required |
| Synthetic harness | `Omega(t)=s·t+Σ A_k sin(omega_k t+phi_k)` with known s; tests OLS/phase-locked bias vs analytic `B_k`, resonance, and no-secular masquerade cases | Estimator validation independent of propagator |
| Constants | MU_EARTH, J2, R_EARTH from src/lab_utils; MU_SUN=132712440018.0, MU_MOON=4902.8001 km³/s², AU=149597870.7 km; documented in results provenance | Reproducibility |

---

## 6. Pre-registered decision rule (quantitative, evaluated before any 18.6-yr reinterpretation)

Headline observable: 1-yr non-additive residual `R = (full−j2_only) − (sun_only+moon_only−2·kepler_only)`-equivalent (exact subtraction in code, documented) via direct-OLS at ascending nodes, dt=60 s, h=600 km. Prior (to be reproduced, not assumed): R≈−1.16e-3 (SSO), −4.21e-3 (90), −2.42e-4 (30) deg/day.

**H-mod is SUPPORTED iff all hold (no empirical coefficient in theory):**
- (a) Mathematical existence: closed derivation from `R_3b` + `Omega(t)` gives explicit `S0 + Σ B_k` with stated order, dimensions, limits; `B_k→0` for frozen `Omega` and for `J2→0` at fixed W (analytic + synthetic-oracle agreement to <1% on bias formula).
- (b) Scaling: 90-d fit cross term `a11` >3σ with H-mod sign; small-scale symmetry set reproduces `a11` within 2σ; prescribed-precession slope `dR/dOmega_dot_p` at fixed 3body is nonzero (>3σ) with H-mod sign near `Omega_dot_p=0`, and `R→0` within detection limit when `Omega_dot_p=0` AND `lambda_J2=0` (detection limit quantified from `lambda=0` scatter).
- (c) Inclination: predicted R sign matches numerical R sign at **all three** headline inclinations (SSO, 90, 30) at 1-yr; at least one extra discriminator (60 or 0 deg) also matches sign.
- (d) Magnitude: predicted |R| within factor 2 of measured |R| at ≥2 of 3 headline inclinations at 1-yr **with no fitted scale** (all constants from model/ephemeris); altitude ratio R(500)/R(800) predicted within factor 2 (tests the extra `n·J2·(R/p)²` steepening vs pure `(a/a3)³`).
- (e) Robustness: dt 60→30 s changes |R| <20%; 4 estimators agree within 30% at each headline i; 4-phase spread <50% of mean |R|; decimated cadence <10%; 90-d vs 365-d R follows predicted `w(omega·W)` window function within factor 2 (not an arbitrary drift).
- (f) Alternatives: octupole bound <10% of |R| at all headline i; frozen-plane R is <25% of modulated R (modulation killed); synthetic-Moon vs real-Moon difference accounted for by `A_k(i3(t),e3)` variation, not by a missing secular; frame/sign mutants (J2 sign, mu3 sign, distance power ±1, rot direction, inclination-factor swap, precession-sign flip) each caught by tests and each moves R as H-mod predicts (sign/scale), not arbitrarily.
- (g) Independence: independent Kaula-track derivation agrees on functional form + scaling + resonance condition (up to O(1) geometry derived, not fitted); independent code path (averaged-element propagator, materially different algorithm) reproduces R sign and order-of-magnitude at headline i.

**Partial / falsified:**
- PARTIAL if (a)+(b)+(g) hold but (c) or (d) holds at only 1 inclination or within factor 3–5: mechanism exists but explains only part; report fraction explained per inclination/altitude.
- FALSIFIED if wrong sign at ≥2 headline i, or magnitude off by >5× at ≥2 i, or (b) fails (no `Omega_dot_p` dependence / no vanish at J2=0), or frozen-plane does not kill R, or an alternative (e.g. mean/osculating + forced nodal mode without J2-frequency dependence) predicts R better (factor-2 vs >5×). Then report strongest surviving mechanism, not H-mod.
- Tolerances are set from actual 1-yr scatter (estimator ~4%, phase/detection ~10–20%, synthesis margin → factor-2 headline). 18.6-yr validation uses the same theory `w(omega·W)` prediction (resonant solar bias persists, non-resonant attenuates); it is a **consistency check**, not a refit. No "new secular canon" claim without `W→∞` evidence; H-mod as derived is long-period.

---

## 7. Implementation map

- `DERIVATION.md` — lead derivation (disturbing function → singly-averaged rate → finite-window bias → scalings → resonance → falsifiers). Symbolic record where practical + dimensional/limit checks.
- `theory_modulation.py` — closed-form `S0`, `A_k/omega_k/B_k`, `predicted_R(h,i,W,phase)` with **no fitted constants**; used by tests + campaign analysis.
- `mission_experiment.py` — streaming RK4 + 6-mode isolation + lambda multipliers + prescribed-precession/frozen-plane controls + 4 estimators + synthetic harness (reuses closure/coupling machinery; corrected `use_j2` guarded).
- `averaged_propagator.py` — independent averaged-element path (secular J2 + singly-averaged lunisolar with prescribed `Omega_dot`); materially different algorithm for (g).
- `run_reduced.py` / `run_scaling.py` / `run_headline.py` / `run_final_18yr.py` — staged campaigns (cheap-first).
- `tests/test_mission_j2_precession_modulated_lunisolar_term.py` + `tests/test_force_mode_isolation.py` — focused + regression tests (§6 discriminators, mutants, oracles, determinism).
- `results/` — `reduced.json`, `scaling.json`, `headline_365d.json`, `final_18yr.json` (+ analyses), `results.json` synthesis, figures.
- Report: `localdocs/reports/mission-j2-precession-modulated-lunisolar-term-2026-09-05.md`; knowledge: `localdocs/knowledge/j2-precession-modulated-lunisolar-term.md` (supersession explicit).

---

## 8. Limitations / non-claims

- Point-mass Sun/Moon, no SRP/drag, no J3+, linear daily-DE441 interpolation; circular LEO ICs; single epoch + phase ensemble (not full epoch sweep).
- Theory is singly-averaged + adiabatic `Omega(t)`; evection/variation short-period content enters via `A_k` bounds and synthetic tests, not full Brown-theory.
- 18.6-yr final is one window (2026–2045); multi-window `W→∞` extrapolation is diagnostic, not a new limit claim.
- No new dependencies; numpy/matplotlib/pytest only.

---

## 9. Reference / supersession

- Exp 018 corrected quadrupole (leading-order, Convention B) preserved; NOT retired. If H-mod holds, canon is **augmented** by an explicitly long-period (not secular) modulation correction with stated domain; if H-mod fails, canon stands with 18.6-yr discrepancy still open and strongest alternative named.
- `mission_lunisolar_closure` PARTIALLY-VERIFIED-WITH-OPEN-QUESTION and `mission_j2_lunisolar_coupling` H1-PARTIALLY-SUPPORTED are preserved; this mission tests the coupling-mechanism refinement those missions recommended.
- Prior numeric figures (a11, R, 18.6-yr rates) are **priors to be reproduced**, not premises.

---

## 10. VERDICT (2026-09-22): H-mod FALSIFIED

**Decision-rule score (pre-registered §6):**

| Gate | Result |
|---|---|
| (a) Mathematical existence, `B_k→0` for frozen Ω / `J2→0` | **PASS** — closed derivation + synthetic-oracle OLS-bias recovery <1% |
| (b) `a11` >3σ with predicted sign | **FAIL** — `a11 = −1.49e-4 ± 0.83e-4` deg/day, SNR **1.79** (prior mission reported −7.85e-4, SNR 6.9; not reproduced with corrected interp) |
| (b) `dR/dΩ̇_p` >3σ at 0 | **FAIL** — local slope **2.42σ** |
| (b) `R→0` at `λ_J2=0` and `Ω̇_p=0` | PASS — ~8e-5 deg/day detection scatter |
| (c) Predicted R sign at all 3 headline i | **FAIL** — SSO wrong, 30° wrong (predicted prograde; measured −2.96e-4 at 30°) |
| (d) Magnitude within 2× at ≥2 i, no fitted scale | **FAIL** — measured SSO R +1.296e-3 vs theory ceiling ≤2.4e-4 (≥5.4×); 30° sign-flip |
| (e) Robustness (dt ladder) | PASS — 60→30 s changes R by 0.9% (<20%) |
| (f) Octupole <10% of R | PASS |
| (f) Frozen-plane R <25% of modulated R | **FAIL at SSO** — frozen R = **52%** of modulated (survives control); 90° identical (trivially, Ω̇_J2=0); 30° 21% (<25%) |
| (g) Independent derivation agreement | PASS — functional form, scalings, resonance condition |
| (g) Independent code path reproduces R sign+order | **FAIL at SSO/30°** — averaged path agrees with closed-form THEORY (SSO +1.77e-4 ≈ S0 1.35e-4 + bias; 30° +4.52e-5 ≈ S0 4.55e-5) but does **not** reproduce the Cowell R (sign flip at 30°, ≥7× at SSO) |

**FALSIFICATION triggers met (any one suffices):** magnitude >5× at ≥2 i (d); (b) SNR failures; frozen-plane does not kill R at SSO.

**Final 18.6-yr validation (consistency, not refit):** `full − j2_only` at h=600 km reproduces `mission_lunisolar_closure`: SSO **−2.288e-2** (closure −2.29e-2), 90° **+5.027e-3** (+4.70e-3), 30° **−4.021e-4** (−3.47e-4) — same signs, magnitudes within 7%/16% (corrected interpolation + streaming implementation). The corrected doubly-averaged quadrupole remains NOT a valid predictor of the osculating OLS rate at LEO.

**Strongest surviving mechanism (recorded, not claimed as canon):** **J2-baseline mean-element leakage** — the paired subtraction `(full−j2) − (isolated)` compares a J2-precessing plane sampled at ascending nodes against stationary-plane runs; the residual is degenerate in λ-scaling and independent of Ω̇_J2 frequency, consistent with every failed discriminator ((b), (c)-30°, frozen 52%). Secondary contributor: quadrupole amplitudes `A_k` from mean lunar geometry (i3=28.58°, e3=0) underestimate the real 2026 lunar field (i3≈18.3°, e3=0.055; cf. audit-020 Track 4 and the coupling mission's synthetic-vs-real Moon ~2× finding).

**Canon impact:** Exp 018 corrected quadrupole stands as the leading-order secular term and is now also bounded as an *osculating-OLS predictor*: the gap between it and Cowell R is real and NOT explained by J2-precession modulation. `mission_j2_lunisolar_coupling`'s "J2-precession-modulated coupling" phrasing is **retracted as a mechanism explanation**; its measured a11/R values remain valid priors (a11 not reproduced at 90 d with corrected interp — flag for follow-up).

**Limitations (honest):** 4-estimator-spread and 4-phase-spread gates (§6e) not executed by the run scripts (estimator theory covered by synthetic OLS-bias oracle instead); adversarial mutant battery not executed in-session (2 hostile-review subagents failed on provider errors); single epoch (2026); circular ICs; quadrupole + point-mass physics.

**Recommended next:** (1) quantify the J2-baseline mean-element/geometry-sampling leakage with a mean-element (Brouwer-style) reference propagator — this is now the leading explanation and is directly testable; (2) re-derive `A_k` with real lunar e3/i3(t) from the byte-pinned DE441 Moon (not secular means); (3) reconcile the 90-d a11 (−7.85e-4 SNR 6.9 prior vs −1.49e-4 SNR 1.8 corrected-interp present mission).

