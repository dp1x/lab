# mission_j2_precession_modulated_lunisolar_term — Final Scientific Report (2026-09-22)

**Verdict:** H-mod (J2-precession-modulated lunisolar coupling) **FALSIFIED** per the pre-registered decision rule. Mission card: `research/orbital-mechanics/missions/mission_j2_precession_modulated_lunisolar_term/README.md` (§10). Knowledge: `localdocs/knowledge/j2-precession-modulated-lunisolar-term.md`.

## 1. Question and prior state

`mission_lunisolar_closure` (18.6-yr DE441) found the corrected doubly-averaged quadrupole `(3/8) n (mu3/mu)(a/a3)^3 sin2(i-i3)/sin i` disagrees with the numerically integrated lunisolar node rate at LEO (SSO: +1.35e-4 predicted vs −2.29e-2 measured, sign wrong). `mission_j2_lunisolar_coupling` measured a non-additive J2×lunisolar residual (74–92% of the combined lunisolar contribution at 1-yr) and named "J2-precession-modulated coupling" as the candidate mechanism — explicitly a hypothesis, not a result. This mission derived the mechanism falsifiably: singly-averaged quadrupole harmonics `A_k` (quadrature, no fitted constants), beat frequencies `omega_k` containing `Omega_dot_J2`, exact finite-window OLS-bias kernel `g(omega_k W, phi)`, and pre-registered gates (README §6) written before full-scale validation.

## 2. Method

- **Theory:** DERIVATION.md (lead) + independent Kaula-track derivation (read-only subagent, returned text) + `theory_modulation.py` (S0, quadrature `A_k`, exact `g()`, octupole bound).
- **Numerics:** streaming RK4 Cowell (`mission_experiment.py`), 6-mode isolation (`use_j2 = mode in ("j2_only","sun_moon_j2")`), corrected DE441 linear interpolation (regression-pinned midpoint test; the closure/coupling missions' sign-reversed-weights bug quantified at ~15% on 90-d R), byte-pinned DE441 Sun+Moon (sha256-pinned provenance in every JSON), IAU-1976 precession, dt=60 s, 4-worker parallel (BLAS threads pinned to 1 after a host commit-charge failure with threaded OpenBLAS — WinError 1455).
- **Controls:** prescribed-precession (active rigid rotation of the full Cartesian state), frozen-plane open-loop counter-rotation at the analytic J2 rate, with a bias-invariance regression guard (`freeze+δ` shifts both raw J2-mode rates by `−δ`, paired R invariant — proving feed-forward subtraction, not closed-loop cancellation).
- **Independent path:** 1D averaged-element propagator (`averaged_propagator.py`) integrating `dOmega/dt = Omega_dot_J2 + S0 + Σ A_k cos(theta_k)` — materially different algorithm, used as a theory-consistency check.

**Remediations during execution (documented in code):** (1) the interpolation bug above; (2) a passive rotating-frame control formulation (Coriolis/centrifugal + rotating-frame observable) was attempted 2026-09-12, never went green (observable sign defect), and was reverted to active semantics with the frozen control re-interpreted as an open-loop common rotation — all campaign results reported here use the active implementation, with unfrozen rows bit-identical across both formulations (verified against the pre-change campaign).

## 3. Results

### 3.1 Headline 365-d force-mode campaign (h=600 km; correct interp)

| i (deg) | combined (full−j2) | isolated (sun+moon−kep) | R = combined−isolated |
|---|---|---|---|
| 97.7876 (SSO) | +1.165e-3 | −1.307e-4 | **+1.296e-3** |
| 90 | +5.007e-4 | −1.697e-4 | **+6.703e-4** |
| 30 | −3.457e-4 | −4.955e-5 | **−2.962e-4** |

Altitude ladder (SSO): R(500 km) = +7.752e-4, R(800 km) = +2.007e-3 deg/day (rises with altitude, qualitatively consistent with third-body scaling).

### 3.2 Scaling and prescribed-rate discriminators (90-d, SSO)

- `a11 = −1.489e-4 ± 0.831e-4` deg/day, **SNR 1.79** (prior mission: −7.85e-4, SNR 6.9 — NOT reproduced with corrected interpolation).
- Prescribed-precession sweep: response-vs-command removes the commanded rate; `R_vs_stationary` vs `Ω̇_p` is non-monotone; local slope at `Ω̇_p=0` = **2.42σ** (needs >3).
- dt ladder 120→60→30 s: R changes 0.9% (<20% gate). Detection scatter ~8e-5 deg/day.

### 3.3 Frozen-plane control (365-d)

| i | frozen R | modulated R | ratio |
|---|---|---|---|
| SSO | +6.71e-4 | +1.296e-3 | **52%** (kill criterion <25% — FAIL) |
| 90 | identical (Ω̇_J2=0, control inert) | +6.703e-4 | — |
| 30 | +6.10e-5 | −2.962e-4 | 21% (<25%) |

### 3.4 Theory vs measurement (the core falsification)

- Theory ceiling (S0 + constructive harmonic bias) ≈ **2.4e-4** deg/day at h=600 km. Measured R = 1.296e-3 (SSO) / 6.70e-4 (90°) / −2.96e-4 (30°): **≥5.4× too large at SSO, sign-flipped at 30°**.
- Independent averaged path reproduces the THEORY (SSO +1.770e-4 vs S0 1.348e-4; 30° +4.522e-5 vs S0 4.547e-5; 90° ≈+3e-6) — internally consistent — but NOT the Cowell R. The discrepancy lives in the osculating-sampling layer.

### 3.5 Final 18.6-yr validation (consistency check, no refit)

`full − j2_only` over one full lunar nodal cycle (~103k node samples per propagation, dt=60 s, byte-pinned DE441): SSO **−2.288e-2**, 90° **+5.027e-3**, 30° **−4.021e-4** deg/day — reproducing `mission_lunisolar_closure` (−2.29e-2 / +4.70e-3 / −3.47e-4) in sign and within 7%/16% in magnitude (corrected interpolation + streaming implementation). Quadrupole-vs-Cowell discrepancy at W=18.6 yr STANDS.

## 4. Decision-rule score

(a) PASS; (b) FAIL×2 (a11 SNR 1.79; prescribed slope 2.42σ) with R→0 PASS; (c) FAIL (SSO, 30°); (d) FAIL; (e) PASS (dt); (f) octupole PASS, frozen-plane FAIL at SSO; (g) derivation PASS, independent-path FAIL at SSO/30°. **FALSIFIED** — multiple independent triggers.

## 5. Interpretation and surviving mechanism

The modulation is real but explains neither the magnitude nor the sign of R, and fails every J2-frequency discriminator. The residual is dominated by **J2-baseline mean-element leakage / geometry-sampling on a precessing plane**: paired subtraction compares node-sampled osculating rates on a precessing plane against stationary-plane runs, producing a J2-frequency-independent offset degenerate in λ-scaling — consistent with all failed gates. Secondary contributor: quadrupole amplitudes from secular-mean lunar geometry underestimate the real lunar field (2026 i3 ≈ 18.3°, e3 = 0.055; cf. audit-020 Track 4; coupling mission's synthetic-vs-real Moon ~2×).

## 6. Canon impact and supersession

- `lunisolar-perturbation-018.md`: corrected quadrupole stands; additionally bounded as an osculating-OLS predictor (gap to Cowell R is real, NOT modulation).
- `j2-lunisolar-coupling.md`: mechanism claim "J2-precession-modulated coupling" RETRACTED; measured a11/R kept as priors, with the a11 value flagged stale (interp fix moved it ~5×).
- `lunisolar-closure-021.md`: open question (18.6-yr SSO discrepancy) remains OPEN.

## 7. Limitations

4-estimator-spread and 4-phase-spread gates not executed by run scripts (estimator-bias physics covered by synthetic OLS-bias oracle); adversarial mutant battery not executed (two hostile-review subagents failed on provider errors); single 2026 epoch; circular ICs; quadrupole point-mass physics; frozen-plane control is an open-loop common rotation (alters Sun/Moon geometry; documented in code).

## 8. Recommended next

1. **Quantify the mean-element/geometry-sampling leakage** with a Brouwer-style mean-element reference propagator (directly testable; now the leading explanation).
2. **Re-derive `A_k` with real lunar e3/i3(t)** from the byte-pinned DE441 Moon.
3. **Reconcile the 90-d a11** (−7.85e-4 SNR 6.9 prior vs −1.49e-4 SNR 1.8 here) — determines how much of the coupling mission's Phase-B sensitivity survives the interpolation fix.

