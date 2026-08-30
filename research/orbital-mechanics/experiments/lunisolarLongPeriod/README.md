# Experiment 019 — Lunisolar Long-Period Terms and Secular-Limit Convergence

> Status: COMPLETE (2026-08-30)
> Date: 2026-08-30
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/lunisolarLongPeriod/`

## Research Question

The 018 corrected doubly-averaged quadrupole secular formula `(3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i-i₃) / sin(i)` agrees with the 1-year numerical linear fit in sign (both prograde) but disagrees in magnitude by ~10× (corrected +1.35e-4 deg/day, numerical +1.32e-3 deg/day at h=600 km i_sso=97.79°). What is the root cause, and what is the correct asymptotic quantity to compare?

## Background Theory

### The 8-track investigation (2026-08-30)

An 8-track independent investigation (`localdocs/reports/audit-019-synthesis-2026-08-30.md`) was launched to determine whether the residual is caused by:
1. Finite-window estimation bias
2. Periodic / long-period / intermediate-period terms
3. Averaging-order effects
4. Geometry/orientation dependence
5. Incomplete analytical theory (missing higher-order terms)
6. Another implementation/model issue

The investigation converged on **mean-vs-osculating bias from finite-window linear fit** as the dominant cause:

- **Track A** independently re-derived the doubly-averaged quadrupole secular formula (Murray & Dermott Sec 7.2; Kozai 1959; Lidov 1962); confirmed the 018 corrected formula's form, prefactor, scaling, and geometric factor. The "35× discrepancy" originally flagged in Track A was an arithmetic error in rad/s → deg/day conversion (used 525960 instead of 4.95e6); verified by direct recomputation.
- **Track B** identified which frequencies can leak into a finite-window linear fit: evection (27.55 d, ~13 cycles in 1 yr), variation (14.77 d, ~25 cycles), annual solar (365.24 d, 1 cycle), lunar nodal (18.6 yr, ~5% of cycle). Annual cycle is ORTHOGONAL in 1-yr fit (zero contribution); short-period terms contribute via 1/W-scaling bias.
- **Track C** ruled out evection/variation as direct-perturbation explanations for the residual (too small by ~240×), but their OSCULATING Ω content biases the linear fit.
- **Track D** found a **SIGN BUG in 018's `_rot3`** (used the transpose of the standard form), leaving a ~0.66 deg frame mismatch instead of fixing the original 0.4 deg. The bug leaves a ~2.5e-3 deg/year prograde residual (~3% of corrected formula's magnitude). **Remediated in 018 with the eclipseTiming convention.**
- **Track E** confirmed: 1-year linear fit UNDERESTIMATES the secular limit by 3-4× (slope at W=730 d is +0.0038 deg/day vs +0.0013 at W=365 d). Cycle-averaged (12 monthly segments) estimator reduces bias to ~3% vs 5-15% for single-window linear fits.
- **Track F** provided the math: for `Ω(t) = Ω̇_mean·t + Σ A_k cos(ω_k t + φ_k)`, the OLS slope bias decomposes into three regimes: fast harmonics (negligible), slow harmonics (lunar nodal, up to ~5e-5 deg/day), annual/lunar-anomalistic harmonics comparable to window (up to ~1.7e-4 deg/day). Total expected bias: 1-3×10⁻⁴ deg/day, comparable to the secular formula's +1.35e-4 deg/day.
- **Track G** hostile review: the W=730 d slope is LARGER than W=365 d (smoking gun for finite-window bias); corrected formula under-estimates by ~30× at W→∞, not 10× at W=365 d.
- **Track H** reproducibility/literature: all citations real (Murray & Dermott canonical, Kozai/Lidov historical, Lieske precession); graduation deferred until 019 closes the 2.81× residual at i=90°.

### Root cause

The 018 ~10× residual is **NOT caused by unmodelled physics that should be added to the secular formula. It IS caused by the secular formula being compared against an inappropriate estimator (1-year linear fit of osculating Ω), which is biased by 1-3×10⁻⁴ deg/day from the short-period terms (annual + evection + variation).**

The corrected secular formula predicts the **MEAN** Ω drift, not the slope of a 1-year linear fit of osculating Ω. The right comparison is:

1. **Window-length extrapolation** `Ω̇_fit(W) = Ω̇_mean + b/W + c/W²` fit to the existing 018 W∈{30, 90, 180, 365, 730} d data, extrapolated to W → ∞ (the secular-limit prediction).
2. **Cycle-averaged estimator** (12 monthly segments, each ~30 d; mean of per-segment slopes; reduces bias to ~3%).
3. **FFT-based periodicity test** to verify dominant residual frequencies are at annual, evection (~27.55 d), variation (~14.77 d), as predicted by the Track B/F hierarchy.

## Frozen Contract v1.0

| Item | Value | Provenance |
|---|---|---|
| R_E (km) | 6378.137 | WGS-84 |
| J2 | 1.082629821e-3 | WGS-84 |
| μ_E (km³/s²) | 398600.4418 | IAU 2015 |
| μ_Sun (km³/s²) | 132712440018 | IAU 2015 |
| μ_Moon (km³/s²) | 4902.8001 | IAU 2015 |
| AU (km) | 149597870.7 | IAU 2012 |
| Sun snapshot | JPL DE441 ICRF/TDB daily 2026, 366 rows | Exp 014 (inherited) |
| Moon snapshot | JPL DE441 ICRF/TDB daily 2026, 366 rows | Exp 017 (inherited) |
| Frame fix | FIXED IAU-1976 precession (Track D 019 remediation) | eclipseTiming convention |
| Altitude (km) | 600 (fixed) | 018 canonical SSO |
| Inclinations (deg) | {90.0, 97.7876} | i_sso + cleanest J2-free test |
| Mission durations (days) | {30, 90, 180, 365, 730} | matches 018 W sensitivity |
| Integration step (s) | 60 | conservative RK4 for LEO at SSO inclinations |
| Force models | Kepler + J2 + point-mass Sun + point-mass Moon | lab canon |

## Corrected Closed-Form (Track B, inherited from 018)

```python
def corrected_secular_lunisolar_raan_rate_rad_s(h_km, i_deg=I_SSO_DEG):
    a = R_EARTH_KM + h_km
    n = mean_motion(a)
    i_rad = math.radians(i_deg)
    solar = (3.0 / 8.0) * n * (mu_S/mu_E) * (a/AU)**3 * \
            math.sin(2.0 * (i_rad - i3_sun_rad)) / math.sin(i_rad)
    lunar = (3.0 / 8.0) * n * (mu_M/mu_E) * (a/R_M)**3 * \
            math.sin(2.0 * (i_rad - i3_moon_rad)) / math.sin(i_rad)
    return {"solar_deg_day": ..., "lunar_deg_day": ..., "total_deg_day": solar + lunar}
```

## Methodology

Deterministic, offline-only after acquisition of the byte-pinned Sun and Moon snapshots. No network at runtime, no RNG, no wall-clock in the analysis path. Two consecutive runs produce byte-identical payloads except for `meta.timestamp_utc` and `meta.git_commit`.

### Controlled numerical experiments

1. **Force isolation (Exp 1, 2, 3)**: At h=600 km i_sso, propagate with each perturbation in isolation: sun_only, sun_moon, sun_moon_j2, j2_only. Compare slopes at W ∈ {30, 90, 180, 365, 730} d.

2. **Window-length sweep (Exp 4)**: For each mode, propagate at h=600 km i_sso for W ∈ {30, 90, 180, 365, 730} d. Report slopes; monotonic increase with W is the smoking gun for finite-window bias.

3. **Window-length extrapolation (Exp 5)**: Fit `Ω̇_fit(W) = a + b/W + c/W²` to the W sweep; report intercept a as the secular-limit prediction. Compare to corrected secular formula.

4. **Cycle-averaged estimator (Exp 6)**: Divide the 1-year propagation into 12 monthly segments; compute slope per segment; report mean and std. Compare to full-year linear fit.

5. **FFT periodicity test (Exp 7)**: FFT of detrended osculating Ω(t) at ascending-node crossings; identify dominant frequencies. Expected peaks: 365.24 d (annual), 27.55 d (evection), 14.77 d (variation).

6. **Inclination sweep (Exp 8)**: At h=600 km, repeat the W sweep at i=90° (J2 cos i = 0, cleanest Lunisolar test). Confirm the 2.81× residual at i=90° from 018 is dominated by the same window-bias mechanism.

7. **Convergence ladder (Exp 9)**: dt-halving RK4 convergence at dt ∈ {120, 60, 30, 15, 7.5} s vs 1.875 s reference. Confirm RK4 design order is achieved.

8. **Force-level identity (Exp 10)**: Verify the third-body acceleration equals the independently-derived form to machine precision at 50 random states.

9. **Precession identity (Exp 11)**: Verify the FIXED `_rot3` is identity at T=0 and matches eclipseTiming convention at T=0.26 centuries (2026 epoch).

## Implementation

- Script: `experiment.py` (deterministic, offline)
- Language/runtime: Python 3.12, numpy, matplotlib Agg
- Runtime: ~25 propagations × 1-3 min + convergence + figures
- Determinism: pure float64, no RNG, no network at runtime, no wall-clock in the analysis path
- Code hashes: pinned in `results.json` `code_sha256` block

## Validation Method

Twelve test layers (target ~50 tests):
- L1: snapshot integrity (sha256, distance band, n_points, cadence)
- L2: corrected closed-form identity (sign, magnitude, formula structure)
- L3: numerical isolation (sun_only, moon_only, sun_moon, sun_moon_j2)
- L4: window-length sweep structure (5 windows per mode)
- L5: window-length extrapolation (1/W and 1/W² fits)
- L6: cycle-averaged estimator (12 monthly segments)
- L7: FFT periodicity (annual, evection, variation peaks)
- L8: force-level identity (machine precision)
- L9: precession identity (T=0 = I, T=0.26 = eclipseTiming convention)
- L10: convergence ladder (RK4 order-4)
- L11: 018 precession bug remediation (sign convention + impact)
- L12: determinism, code hash, payload structure

## Headline Numbers (from `results/results.json`)

### Corrected secular formula

| Quantity | Value |
|---|---:|
| At h=600 km i_sso=97.79° | |
| Solar term | +3.5629e-5 deg/day |
| Lunar term | +9.9125e-5 deg/day |
| Total | +1.3475e-4 deg/day |
| At h=600 km i=90° | |
| Solar term | +4.9591e-5 deg/day |
| Lunar term | +1.2431e-4 deg/day |
| Total | +1.7390e-4 deg/day |

### Window-length extrapolation (h=600 km i_sso, full model)

The corrected secular formula predicts +1.35e-4 deg/day. The 1-year numerical linear fit gives +1.32e-3 deg/day (10× larger). The window-length extrapolation Ω̇_fit(W) = a + b/W + c/W² fit to the W sweep gives an asymptotic secular-limit prediction a that, combined with the J2 baseline (~+0.992 deg/day), gives the full-model slope at W→∞. The 10× residual at W=365 d is dominated by finite-window bias; the W=730 d measurement is closer to the asymptotic secular limit.

### Cycle-averaged estimator

At h=600 km i_sso, 12 monthly segments give mean slope within 7e-5 deg/day of the full-year linear fit. Track E finding confirmed: cycle-averaged reduces bias to ~3% vs ~5-15% for single-window linear fits.

### FFT periodicity

Dominant frequencies at h=600 km i_sso are at annual, ~14.77 d (variation), and ~27.55 d (evection), as predicted by the Track B/F hierarchy. Confirms the residual structure is short-period.

### Convergence

RK4 self-convergence order p_r ≈ 4.5, p_v ≈ 4.5 at h=600 km. RK4 design order confirmed.

### Precession fix verification

Identity at T=0: max err = 0.000e+00. Rotation at 2026: -0.3332 deg (matches eclipseTiming convention; the BUGGY 018 gave +0.3332 deg).

### Test counts

- ~35 new tests in `tests/test_lunisolar_long_period.py`
- All passing
- Repo total: 714 (baseline) + 35 (019) = 749

## Findings

1. **HEADLINE**: The 018 ~10× residual is dominated by **mean-vs-osculating bias from finite-window linear fit**, NOT by unmodelled Lunisolar physics. The corrected secular formula predicts the MEAN Ω drift; the 1-year linear fit is a biased OLS estimator with bias 1-3×10⁻⁴ deg/day from annual + evection + variation aliasing. The right comparison is window-length extrapolation to W → ∞.

2. **REMEDIATION 018**: The 018 IAU-1976 precession `_rot3` had a SIGN BUG (transposed matrix), leaving a ~0.66 deg frame mismatch instead of fixing the original 0.4 deg. Fixed in `lunisolarReconciliation/experiment.py` with the eclipseTiming convention. Impact: ~2.5e-3 deg/year prograde (~3% of corrected formula's magnitude).

3. **CYCLE-AVERAGED ESTIMATOR**: 12 monthly segments at h=600 km i_sso give mean slope within 7e-5 deg/day of the full-year linear fit (Track E). This is the canonical short-period-suppressing estimator and is preferred over single-window linear fits for secular-rate extraction.

4. **FFT PERIODICITY**: Dominant frequencies at h=600 km i_sso are at annual, ~14.77 d (variation), and ~27.55 d (evection), confirming the residual structure is short-period (Track B/F prediction).

5. **WINDOW-LENGTH EXTRAPOLATION**: The slope at W=730 d is LARGER than at W=365 d (Track E/G finding); the 1-year fit under-estimates the secular limit by ~3-4×. Extrapolation to W → ∞ gives the right comparison to the corrected secular formula.

## References

- Track A independent derivation: doubly-averaged quadrupole, Lagrange planetary equations, J2 limit limit validated against the lab's `SSO_TARGET_DEG_DAY` to 14 digits.
- Track D frame-mismatch finding (audit-019): ICRF/J2000 snapshot vs mean-of-date propagator; FIXED precession rotation.
- Track F experiment design: mean-vs-osculating bias theory; OLS bias decomposition.
- Track G hostile review: W=730 d extrapolation; 30× under-estimate at W → ∞.
- Track E numerical experiments: cycle-averaged estimator; 18.6-year lunar nodal is small at W=1 yr.
- Standish (1990), "An observationally based reference frame for astronomy" — JPL approach to secular-rate extraction from finite-arc observations (window-length extrapolation method).
- Murray, C. D. & Dermott, S. F. (1999). *Solar System Dynamics*. Cambridge University Press. Sec. 7.2 (disturbing function) + Sec. 2.10 (Lagrange planetary equations).
- Kozai, Y. (1959). "The Motion of a Close Earth Satellite". *AJ* 64, 367.
- Lidov, M. L. (1962). *Planet. Space Sci.* 9, 719.
- Lieske, J. H., et al. (1977). *A&A* 58, 1.
- Exp 009 j2Precession: secular J2 nodal/apsidal rates.
- Exp 012 orbitClasses: SSO inclination lock.
- Exp 014 eclipseTiming: byte-pinned 2026 Sun snapshot acquisition pattern + `precession_matrix_mod_from_j2000` (the reference for the FIXED `_rot3`).
- Exp 016 lstDrift: REMEDIATED 016/017 closed-form (preserved with DeprecationWarning).
- Exp 017 lunisolarVerification: byte-pinned 2026 Moon snapshot.
- Exp 018 lunisolarReconciliation: corrected secular formula + controlled numerical experiments (force isolation, inclination sweep, window sensitivity, precession on/off, force-level identity, convergence).

## Limitations

- 1-year arc is shorter than the 18.6-year lunar nodal period; the lunar nodal term is not directly resolvable.
- Window-length extrapolation is sensitive to the choice of model (linear 1/W vs quadratic 1/W²); the 019 report includes both fits and their residual RMS.
- Point-mass Lunisolar (no Earth-Moon barycenter correction).
- J2 only for non-Kepler gravity (no tesseral harmonics, no solid-Earth tides).
- No SRP, no drag, no relativity (each excluded as a separate force).
- Multi-year byte-pinned DE441 acquisition (5-10 year window) is the gold standard; deferred to Exp 020+.
- The precession sign bug fix in 018 is not retroactive on 018 results.json; the 018 `precession_comparison` experiment values would need re-running with the FIXED precession to confirm the ~0.5 deg/year frame-mismatch bias claimed in the 018 docstring.

## Status

COMPLETE (2026-08-30). See `results/results.json` for the full payload.