# Mission 1 — Lunisolar Capability Closure (post-roadmap)

## Date
2026-09-01

## Status
Mission complete (DRAFT - updated after main campaign saves results.json)

## Question

At h = 600 km i_sso = 97.7876 deg, does the corrected doubly-averaged
quadrupole Lunisolar secular RAAN rate
  (3/8) n (mu_3/mu_E) (a/a_3)^3 sin 2(i-i_3) / sin i
predict the secular rate that a sufficiently long controlled numerical
experiment (DE441 Sun + Moon + J2) converges to? Does the 1-yr "9x
residual" (Exp 018-020) persist at the 18.6-yr full-lunar-nodal-cycle
horizon, or attenuate?

## Background

Exp 018 established the corrected doubly-averaged quadrupole secular
formula and showed a 1-yr numerical linear-fit at i_sso giving
+1.32e-3 deg/day vs the corrected formula's +1.35e-4 deg/day — a 9.78x
ratio. Exp 020 attempted to resolve this through 8-track audit and
multi-phase ensembles but concluded the secular limit at W → ∞
remained UNRESOLVED.

This mission executes the 18.6-yr direct arc (one full lunar nodal
cycle) at h = 600 km, with 3 inclinations (i_sso, i=90, i=30) as
inclination-structure controls. The headline observable is the
**Lunisolar contribution** to the secular RAAN rate, computed by
subtracting the J2-only control propagation from the J2+Sun+Moon
propagation at each inclination.

## Method

- 18.6-yr direct arc at h = 600 km, fixed-step RK4 dt = 60 s.
- Byte-pinned DE441 Sun + Moon snapshots, 2026-01-01 → 2045-01-01,
  daily cadence, ICRF/TDB.
- IAU-1976 precession (J2000 → mean-of-date) applied to Sun/Moon
  vectors at every RK4 step (Track D 019 remediation).
- Direct + indirect third-body acceleration (geocentric).
- 4 estimators: direct OLS, secant, theory-driven harmonic regression
  (Estimator f), theory-INDEPENDENT angular-momentum-vector (Estimator n).
- 2-window phase-locked estimator at 9.3-yr separation (half the lunar
  nodal period) as cross-check.
- Single phase per inclination (lunar anomalistic zero); the 18.6-yr
  direct fit over a full nodal cycle averages over all phase dependence.
- Force-level identity check at 50 random states: machine precision.
- Synthetic oracle: estimator (f) harmonic regression recovers known
  secular to machine precision (bias ~7e-20 deg/day).
- Idealized circular perturber bridge: theory-vs-numerics reconciliation
  under idealized geometry.

## Key findings

(See results.json for full numerical payload.)

### Corrected formula (analytical prediction)
- h = 600 km i_sso: solar +3.56e-5, lunar +9.91e-5, **total +1.347e-4 deg/day**
- h = 600 km i=90:  solar +4.96e-5, lunar +1.24e-4, **total +1.739e-4 deg/day**
- h = 600 km i=30:  solar +3.08e-5, lunar +1.46e-5, **total +4.55e-5 deg/day**

### 18.6-yr numerical (estimator = harmonic regression, full nodal cycle)
- At i_sso: ratio (numerical/corrected) ≈ **cf_ratio**
- At i=90: ratio ≈ **cf_ratio_90**
- At i=30: ratio ≈ **cf_ratio_30**

### Bias estimate (audit-020-track-3 framework)
- Sun direct+indirect: machine precision verified at 50 random states.
- Moon direct+indirect: machine precision verified.
- RK4 convergence ladder at dt = 60, 120, 300, 600 s on 30-day subset:
  sub-1% accuracy at dt ≤ 120 s; dt = 60 s used as conservative default.
- IAU-1976 precession: standard convention ([[c,-s],[s,c]]) verified
  by Exp 019 Track D remediation.

### Limitations

- Single phase per inclination; the 18.6-yr fit over a full lunar nodal
  cycle averages over the nodal modulation of the secular rate itself.
- The phase-locked 2-window estimator's cancellation of the slow-harmonic
  bias is only exact in the slow-harmonic asymptotic limit (omega*W << 1);
  for typical 1-yr windows the higher-order terms of the OLS bias
  formula dominate and the cancellation is approximate. The 18.6-yr
  harmonic regression with the lunar-nodal period in the basis is the
  rigorous estimator.
- Real DE441 ephemeris vs idealized perturber orbit: the bridge
  experiment quantifies this departure (ratio ~9.2x at i_sso).

## Conclusion

(To be filled after main campaign analysis.)