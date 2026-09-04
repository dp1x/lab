# Knowledge Note — J2 × Lunisolar Coupling at LEO

**Date**: 2026-09-03
**Status**: H1-PARTIALLY-SUPPORTED — J2 × Lunisolar coupling is real and dominant at LEO, but the mechanism is **J2-precession-modulated Lunisolar coupling**, not the direct Lie-transform cross term predicted by naive perturbation theory.

## Key findings

### 1. The direct Lie-transform J2 × 3b cross term is too small
Literature (Murray & Dermott §2.10, Brouwer & Clemence §11/17, Cook lunisolar papers) shows the direct cross term scales as `J2 × (n₃/n)²`. At LEO SSO this is ~10⁻⁹ — far too small to explain the 018 / 020 / 021 discrepancies.

### 2. The kinematic J2-precession modulation IS large
The 1-yr arc force-mode decomposition (Phase C, post-remediation) shows:
- Luni_combined (full − J2) at i_sso: -1.44e-3 deg/day
- Luni_isolated (Sun-only + Moon-only): -2.76e-4 deg/day
- **R_J2x3b (cross-coupling residual): -1.16e-3 deg/day = 80.8% of Luni_combined**

The majority of what is conventionally called "Lunisolar RAAN drift" at LEO is actually the **J2-precession-modulated Lunisolar coupling**: J2 secular Ω drift (~1 deg/day at h=600 km i_sso) modulates the orbit-plane orientation in the Sun/Moon field. When Sun and Moon are propagated WITHOUT J2, the orbit plane is stationary and this modulation does not occur — hence the isolated modes give much smaller Lunisolar rates.

### 3. The perturbative scaling confirms the cross-coupling
Phase B (90-d, i_sso) 2-D polynomial fit:
- Cross-term a11 = -7.85e-4 deg/day, SNR = 6.89 → **statistically significant**
- J2² coefficient a20 = +5.21e-3, SNR = 36.6 → J2² second-order is real
- 3b² coefficient a02 = -8.8e-7, SNR = 0.01 → 3b is linear (no self-coupling)

The cross term scaling as λ_J2 × λ_3body is unambiguous evidence of a genuine coupling, not an estimator artifact.

### 4. The corrected formula is missing this term
The doubly-averaged quadrupole formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i` averages over the satellite's mean anomaly AND assumes a fixed orbit plane. The J2-precession modulation is a **mean-osculating coupling** that the formula cannot capture.

## Why this matters

The lab's corrected Lunisolar secular canon is INCOMPLETE. The corrected formula is the leading-order doubly-averaged term; it must be augmented by the J2-precession-modulated term to be a valid asymptotic predictor of the osculating-element secular rate at LEO under real DE441 ephemerides.

The mission's recommended next action is to derive this J2-precession-modulated term explicitly, parameterized by the J2 precession rate and the orbital geometry. This is a NEW term to add to the lab's secular canon.

## Implementation artifacts

- `research/orbital-mechanics/missions/mission_j2_lunisolar_coupling/`: streaming RK4 + mode isolation + perturbative multipliers + 19 tests
- `localdocs/reports/mission-j2-lunisolar-coupling-2026-09-03.md`: full scientific report
- `localdocs/knowledge/lunisolar-closure-021.md`: prior mission's knowledge note (supersession record)

## Bug class lesson (audit-grade)

The mission started with a critical implementation bug (`use_j2 = mode != "kepler_only"` instead of `use_j2 = mode in ("j2_only", "sun_moon_j2")`). This caused every "non-Kepler" mode to silently include J2, contaminating the force-mode decomposition with catastrophic (1 deg/day) residuals. The analysis script's R_J2x3b sanity check caught it. The remediation commit documents the bug, the fix, and the corrected numerical results. This is consistent with the audit-018 / 019 / 020 doctrine that adversarial-style sanity checks catch implementation errors before they propagate.

## Connection to prior chain

- Exp 015–020: established the Lunisolar discrepancy chain
- mission_lunisolar_closure: 18.6-yr arc refutes the corrected formula's asymptotic validity
- **THIS MISSION**: identifies the J2-precession-modulated Lunisolar coupling as the dominant missing physics; recommends deriving it explicitly
- Next mission (recommended): derivation of the explicit J2-modulation term + 18.6-yr re-run with corrected mode isolation
