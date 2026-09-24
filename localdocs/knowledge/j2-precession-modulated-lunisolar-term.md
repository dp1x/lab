# J2-precession-modulated lunisolar term at LEO — FALSIFIED as the coupling mechanism

**Date:** 2026-09-22. **Status:** FINAL (mission `mission_j2_precession_modulated_lunisolar_term` COMPLETE).
**Supersession:** this note RETRACTS the mechanism claim in `j2-lunisolar-coupling.md` ("J2-precession-modulated Lunisolar coupling" as the explanation of the non-additive residual) and AUGMENTS `lunisolar-closure-021.md` (the 18.6-yr discrepancy it left open is NOT closed by this mechanism). The corrected doubly-averaged quadrupole formula (`lunisolar-perturbation-018.md`) remains the leading-order secular canon and is now additionally bounded as an osculating-OLS predictor.

## Question

Can J2-driven nodal precession `Omega_dot_J2 = -(3/2) n J2 (R/p)^2 cos i`, entering the singly-averaged lunisolar node response as beat frequencies `omega_k` (with finite-window OLS bias `B_k = A_k g(omega_k W, phi)`), quantitatively explain the J2×lunisolar residual measured by the coupling mission — without any fitted coefficient?

## Verdict: FALSIFIED (pre-registered rule, mission README §6)

- **Exists but too small:** the modulation is a real adiabatic long-period effect (lead derivation + independent Kaula track + independent averaged-element propagator all agree on functional form, scalings, sun-sync resonance). Magnitude ceiling ~2.4e-4 deg/day at h=600 km; measured residuals 2.9e-4…1.3e-3 deg/day — 5–17× too large, sign wrong at 30° (measured retrograde, predicted prograde).
- **Frequency-dependence discriminators fail:** perturbative cross `a11` SNR 1.79 (needs >3); prescribed-precession local slope 2.42σ (needs >3). The residual does NOT scale with the J2-precession frequency the way H-mod requires.
- **Frozen-plane control fails at SSO:** with J2 precession counter-rotated (open-loop, bias-invariance-guarded), 52% of the modulated R survives (pre-registered kill criterion <25%). Most of R is J2-frequency-INDEPENDENT.
- **Independent path confirms the theory, not the measurement:** the averaged-element propagator reproduces S0+bias (SSO +1.77e-4, 30° +4.52e-5) — i.e. the theory is internally consistent — while the Cowell osculating-node campaign measures R 5–17× larger with a 30° sign flip. The gap is in the mean-vs-osculating / sampling layer, not in the averaged theory.
- **18.6-yr consistency:** final validation reproduces `mission_lunisolar_closure` (SSO −2.288e-2, 90° +5.027e-3, 30° −4.021e-4 deg/day) with corrected interpolation. The quadrupole-vs-Cowell discrepancy at the full lunar nodal cycle STANDS, unexplained by modulation.

## Durable knowledge

1. **Osculating-OLS sampling on a J2-precessing plane is NOT a small correction.** Paired force-mode subtraction on an osculating, ascending-node-sampled observable carries a geometry-dependent offset (the leading surviving explanation for R) that is invisible to averaged theory and degenerate in lambda-scaling. Any future LEO node-rate campaign must include a mean-element (Brouwer-style) reference before attributing residuals to physics.
2. **Prior a11 = −7.85e-4 (SNR 6.9, 90 d) did not reproduce** with corrected DE441 interpolation (−1.49e-4, SNR 1.79). The shared-lineage interpolation bug fix (`r_lo + frac*(r_hi − r_lo)`) changed the 90-d estimate by ~5×; treat Phase-B coupling-sensitivity numbers from pre-fix missions as stale.
3. **Methodology (reusable):** open-loop counter-rotation controls with bias-invariance regression (`freeze + delta` shifts raw rates by `−delta`, paired R invariant); exact OLS-bias kernel `g(omega W, phi)` with synthetic oracle; quadrature-derived quadrupole amplitudes (no fitted constants); independent 1D averaged-element path as a theory-consistency check (NOT a measurement proxy).
4. **Open question carried forward:** what produces the 30° sign flip and the 18.6-yr SSO retrograde −2.29e-2? Leading candidate: mean-element leakage + real lunar e3/i3(t) geometry (2026 lunar i3 ≈ 18.3°, not the 28.58° secular mean). Quantifying this is the recommended next mission.

## What to cite

- `research/orbital-mechanics/missions/mission_j2_precession_modulated_lunisolar_term/` (README §10 verdict; DERIVATION.md; theory_modulation.py; mission_experiment.py; averaged_propagator.py; results/*.json incl. final_18yr.json + independent_check.json)
- `localdocs/reports/mission-j2-precession-modulated-lunisolar-term-2026-09-22.md`
- Prior chain: `lunisolar-closure-021.md`, `j2-lunisolar-coupling.md` (mechanism claim retracted; measurements kept as priors), `lunisolar-perturbation-018.md`, audits 018–020.

## Limitations

- 4-estimator-spread and 4-phase-spread gates not executed by run scripts (synthetic OLS-bias oracle covers estimator-bias physics); adversarial mutant battery not executed (subagent provider failures); single 2026 epoch; circular ICs; quadrupole point-mass physics; frozen-plane control is an open-loop common rotation (alters Sun/Moon geometry; documented in code).
