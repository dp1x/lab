# mission_mean_element_leakage — Is the "J2 × lunisolar mystery interaction" a measurement artifact of a breathing J2 baseline and node-sampling, or missing physics?

**Mission type:** Discrepancy mission (LAB_CONSTITUTION.md §2.5), follow-on to `mission_j2_precession_modulated_lunisolar_term` (H-mod FALSIFIED), `mission_j2_lunisolar_coupling`, and `mission_lunisolar_closure`.
**Status:** ACTIVE — card + pre-registered decision rule written BEFORE Phase 0 results.
**Constitutional authority:** LAB_CONSTITUTION.md §§2.3, 3, 9, 12–13.
**Date established:** 2026-09-24.
**Budget:** ≤4 hr wall on 4 workers, ≤1 GB RAM, local commodity compute only. Phase 0 ≈ 5 min (1-yr) / ≈ 40 min (18.6-yr); estimator ladder minutes; no GPU.

---

## 1. Question

Every mission in the lunisolar chain measures the orbit's node angle **only at equator crossings** and subtracts a **constant** J2 baseline taken from a separate J2-only run:

```
R = (sun_moon_j2 − j2_only) − [(sun_only − kepler) + (moon_only − kepler)]
```

That single number `R` is the entire "J2 × lunisolar interaction" claim (1-yr: `+1.296e-3` deg/day at SSO; 18.6-yr: SSO `−2.288e-2`). Two mechanisms can create `R` without any new physics:

1. **H-baseline (breathing baseline).** The J2 drift rate depends strongly on the orbit's mean tilt and size: `dΩ̇_J2/di = +0.1258 deg/day per degree`, `dΩ̇_J2/da = −4.94e-4 deg/day per km` (h=600 km, i=97.79°). Sun/Moon slowly change the mean `i` and `a`, so the real run's J2 baseline *differs from the constant one we subtract*. Required drift to explain the observed numbers: **0.010° (1-yr R)** and **0.18° (18.6-yr SSO)** — both plausible for lunisolar secular drift at LEO.
2. **H-sampling (node-crossing estimator bias).** The node is sampled only when the satellite crosses the equator. For a wobbling, precessing orbit this is a biased estimator of the mean drift (audit-019 Track F: bias `1–3e-4 deg/day`, "enough to explain a factor of 2–10").

If either holds, the chain's headline finding is a **measurement artifact**, not a J2×lunisolar force interaction.

## 2. Discriminating predictions (written before results)

| Observable | H-baseline predicts | H-sampling predicts | Falsified plane-swivel (H-mod) predicted |
|---|---|---|---|
| Inclination dependence | ∝ `sin i` — **max at i=90°**, zero at i=0° | ∝ estimator geometry, non-zero at 90° | ∝ `cos i` — max at SSO, **zero at 90°** |
| Altitude dependence | ∝ `a^-3.5` (baseline sensitivity) | weak | ∝ third-body scaling |
| Time-varying baseline removes R? | yes (≥70% at all three `i`) | only partly | no (52% survived the frozen-plane control) |
| Mean-element estimators (h-vector, mean integration) | agree with baseline-corrected value | disagree with node-crossing value by ≤25% | agree with nothing |

The `sin i` vs `cos i` split is decisive: we already own data at `i=97.79°` and `i=90°`.

## 3. Hypotheses (pre-registered)

- **H-baseline:** `R` is the constant-baseline subtraction error `⟨Ω̇_J2(a_f, i_f) − Ω̇_J2(a_j, i_j)⟩`.
- **H-sampling:** `R` is the osculating node-crossing estimator bias relative to the mean-element rate.
- **H-phys:** genuine unmodelled physics (real lunar e3/i3(t) geometry, forced nodal modes) — the arm that must be re-derived with time-varying lunar elements if both H-baseline and H-sampling fail.
- **H-null:** none of the above; report OPEN with the strongest surviving candidate.

## 4. Pre-registered decision rule (binding)

Measured against the 1-yr campaign values (`R` = +1.296e-3 / +6.703e-4 / −2.962e-4 deg/day at SSO / 90° / 30°) and the 18.6-yr SSO value (−2.288e-2):

- **H-baseline CONFIRMED** if the time-varying-baseline correction removes **≥70% of `R` at all three inclinations** (1-yr) **and** the residual 18.6-yr SSO magnitude falls below **25% of 2.288e-2** (`< 5.7e-3`).
- **H-sampling CONFIRMED** if mean-element estimators (A: angular-momentum-vector; C: mean-element integration) give **<25% of the node-crossing `R`** at all three inclinations **while the baseline correction alone does not**.
- **PARTIAL (both):** report the split, per inclination and altitude; the two contributions must be separately quantified and must not double-count.
- **H-phys required** if the joint correction leaves **>50% of `R`** at any inclination; then re-derive the lunisolar amplitudes with real lunar `e3/i3(t)` from the byte-pinned DE441 Moon before any canon claim.
- **FALSIFIED (all artifact arms)** if no combination of baseline + sampling correction removes **≥50% of `R`** at the two `i` where H-mod already failed (SSO and 30°).
- **Gate before any real-data claim:** the estimator ladder must recover a *known* injected drift from synthetic data to <1% (audit-020 Track 3 style oracle). If the gate fails, no real-data result is reported as evidence.

## 5. Phases

0. **Baseline-breathing pre-check (this stage).** 1-yr, 3 inclinations, `j2_only` + `sun_moon_j2`, recording osculating `a, e, i, Ω` at every ascending-node crossing; compute the baseline-difference rate and compare with `R`. Then the 18.6-yr version if the 1-yr arm is informative.
1. **Estimator ladder.** (A) angular-momentum-vector drift (no node detection); (C) mean-element integration (upgrade the predecessor's 1D averaged propagator); bridge estimators: FFT subtraction + multi-window joint fit (audit-020 Track 2 §8.4). Synthetic oracle gate.
2. **Mean-element re-test of the mystery.** Same 6 force modes × 3 inclinations, 1-yr and 18.6-yr, every estimator, with and without the breathing-baseline correction.
3. **Real-geometry arm (only if needed).** Re-derive lunisolar amplitudes with time-varying lunar `i3(t), e3(t)` from the pinned DE441 Moon.
4. **Docs + staged commits:** card → Phase 0 (+tests, results) → ladder → re-test → report/knowledge/roadmap/AGENTS.

## 6. Implementation map

- `phase0_baseline_breathing.py` — streaming RK4 sharing the predecessor's force model, DE441 snapshots and corrected interpolation (imported, so the interpolation fix is not duplicated); records elements at crossings; computes the baseline-difference rate and the measured `R`; `--years` selects the arc.
- `tests/test_mean_element_leakage.py` — propagator agreement vs predecessor `propagate()` (guards divergence), analytic-vs-numeric J2 sensitivity check, synthetic closure of the baseline formula.
- `results/` — `phase0_1yr.json`, `phase0_18yr.json` (provenance: snapshot sha256, code hashes; `wall_s` deliberately NOT recorded).
- Report `localdocs/reports/mission-mean-element-leakage-<date>.md`; knowledge note `localdocs/knowledge/mean-element-leakage.md`.

## 7. Limitations / non-claims

- Phase 0 uses osculating elements at crossings (not a full Brouwer mean-element transform) — a *lower bound* on the effect; the audit trail (019 Track F, 020 Track 2/5) says the full transform is the wrong tool for this problem.
- No new physics is added; a null result on the artifact arms does NOT imply a force-level J2×lunisolar coupling.
- Single epoch (2026), circular ICs, point-mass Sun/Moon, quadrupole-level lunisolar.
- H-mod (plane swivel) remains falsified; this mission does not reopen it.

## 8. References / supersession

- Predecessor: `mission_j2_precession_modulated_lunisolar_term` (H-mod FALSIFIED, 2026-09-22, commit `3389949`) — retracts the coupling mission's mechanism claim.
- `localdocs/reports/audit-019-track-F-mean-vs-osculating.md` (bias theory + verdict), `audit-020-track-2-periodic-terms-and-bias.md` (§7.3 Brouwer subtraction is the wrong tool; §8.4 FFT + multi-window bridge), `audit-020-track-5-independent-estimator.md` (§3 estimator ladder C + A; §4 code design).
- If H-baseline or H-sampling is confirmed, `localdocs/knowledge/lunisolar-closure-021.md`'s 18.6-yr "sign disagreement" is reclassified as a measurement artifact and the roadmap row's "Brouwer-style propagator" wording is corrected.

---

## 9. PHASE-0 VERDICT (appended 2026-10-03; §4 above is unchanged)

**Verdict: H-baseline NOT CONFIRMED. The mission's central question is not answerable
on an 18.6-year arc.** Full evidence: `localdocs/reports/mission-mean-element-leakage-2026-10-03.md`.

§4 was evaluated clause-by-clause without retuning any threshold:

| §4 clause | Verdict |
|---|---|
| 1a (≥70% of `R` at all three inclinations, 1 yr) | **FAIL** — 30° gives 9.8% (smoothed) / 27.5% (raw) |
| 1b (18.6-yr SSO residual < 5.72e-3) | **PASS**, but §4 never defines "residual" |
| H-baseline CONFIRMED (1a ∧ 1b) | **FAIL** |
| H-sampling CONFIRMED | **UNDEFINED** — estimator data did not exist at write time |
| PARTIAL / H-phys / FALSIFIED | **UNDEFINED** at 18.6 yr (2-mode arm; `R_measured` correctly omitted) |
| Oracle gate (<1%) | was **UNSATISFIED** when Phase 0 was written; satisfied 2026-10-03 for EST-1/2/1b |

**The blocking finding.** The residual is **not a window-independent quantity**. Measured
on the same propagation with the same code, only the analysis window varies:

- residual ranges **+1.69e-3 … −1.90e-3 deg/day** and **changes sign** across sliding
  windows; its peak-to-peak swing is 1.4–4.2× its own central value;
- the committed 18.6-yr value (−8.56e-4) is **smaller than its own window-to-window
  standard deviation**;
- the detrended node-difference series is **99.4%** explained by one ~35 deg sinusoid at
  a ~9000-day period; the arc spans only **0.894 lunar nodal cycles**, over which the
  secular term and the long-period term are **not separable** (the joint fit is degenerate).

A secular rate cannot flip sign with the window. The residual is a **single-phase sample
of long-period forced content**, not evidence for a dynamical interaction.

**Falsified this session** (all independently measured):
- the "165 km semi-major-axis decay" — **RK4 step-size artifact**, 5th order in dt
  (p = 5.001/5.000/5.000), −143.0 km predicted from `kepler_only` vs −164.8 km committed;
- `a`-drift contamination of the residual — **common-mode, cancels** (the pure-`a`
  channel is 0.034% of `D` at SSO and 0.000% at 90°);
- mean-element / Jensen-gap inadequacy — **bounded at 0.17% of the residual** at SSO;
  the short-period `Ω` term is **exactly zero at the ascending node** (`u = 0`) and
  aliases to zero in a per-orbit OLS slope;
- node-crossing sampling bias — **exonerated**: a node-detection-free estimator reproduces
  the committed `combined` to **5.8e-9 deg/day** (oracle-gated at 0.083–0.257%);
- a claimed OLS-conditioning artifact of +2.17e-3 deg/day — **does not reproduce**
  (4.7e-18 on the real grid; re-scoring moves the residual by 8.9e-16).

**Status: still ACTIVE**, for one reason only: the arc cannot contain the signal. The
evidence-indicated next step is a **multi-cycle, phase-decorrelated determination of the
long-period nodal term**; formal selection is a human decision (§4.3).

