# mission_lunisolar_closure — Lunisolar Capability Closure

**Mission type:** Validation + Capability hybrid
**Status (2026-09-02):** active; campaign running
**Constitutional authority:** `LAB_CONSTITUTION.md §13.1`
**Roadmap position:** first mission after adoption of `LAB_CONSTITUTION.md`
(POST_ROADMAP_PROBE §13.1)

---

## Question

At h = 600 km i_sso = 97.79 deg, does the corrected doubly-averaged
quadrupole Lunisolar secular RAAN rate

    dΩ/dt = (3/8) n (μ₃/μ_E) (a/a₃)³ sin 2(i − i₃) / sin i

predict the secular rate that a sufficiently long controlled numerical
experiment (DE441 Sun + Moon + J2) converges to? Does the 1-yr
"~10× residual" (Exp 018-020) persist at the 18.6-yr full-lunar-nodal-
cycle horizon, or attenuate?

The 018 corrected formula gives **+1.35×10⁻⁴ deg/day (prograde)** at
h=600 km i_sso; Exp 018/020 1-yr numerical is **+1.3×10⁻³ deg/day
(prograde)**, a ~10× ratio in the right direction.

This mission executes an **18.6-yr direct RK4 arc** (one full lunar
nodal cycle) with **3 inclinations** (i_sso, i=90, i=30) as inclination-
structure controls, byte-pinned DE441 Sun + Moon snapshots, and 4
estimators (direct OLS, secant, theory-driven harmonic regression,
theory-INDEPENDENT angular-momentum-vector).

---

## FET verdict (passed)

Per `LAB_CONSTITUTION.md §9`, this mission passes all 7 gates:

| Gate | Verdict | Notes |
|---|---|---|
| 1. Reasoning/compute ratio | PASS | Estimator design + bias analysis + adversarial review; ~2-3 hr compute |
| 2. Independent validation | PASS | Byte-pinned DE441 + corrected formula + synthetic oracle + force-level identity |
| 3. Durable scientific knowledge | PASS | Closes the open 020 question; graduates reusable estimator library |
| 4. Hypothesis-distinguishing | PASS | Distinguishes "corrected formula right at W → ∞" from "corrected formula under-estimates by 10×" |
| 5. Adversarial-survivable | PASS | Pre-registered estimators + decision rules |
| 6. Capability-advancing | PASS | Multi-year ephemeris acquisition + streaming propagator + harmonic regression estimator |
| 7. Deterministic on modest resources | PASS | < 5 hr single-core, ~2 MB repo-resident data |

The mission is the top-ranked candidate in both
`POST_ROADMAP_PROBE.md §13.1` and `LAB_CONSTITUTION.md §13.1`.

---

## Hypothesis

**H1 (strong)**: The 018 corrected doubly-averaged quadrupole formula
is the right asymptotic predictor of the Lunisolar RAAN secular rate
at LEO SSO; the ~10× 1-yr numerical residual is mean-vs-osculating +
finite-window linear-fit bias.

**H0 (null)**: The 019 extrapolation (+0.0036 deg/day, 27× the
corrected formula) is the true secular limit; the 018 formula
under-estimates by an order of magnitude.

**H-uncertainty**: The discrepancy is a model-order correction (J2 ×
Lunisolar coupling + lunar-evection + variation + annual solar
aliases) that we cannot yet bound.

The mission rejects H0 if the 18.6-yr harmonic regression Lunisolar
rate at i_sso agrees with the corrected formula within ±50%; it accepts
H1 within that bound.

---

## Scientific protocol (frozen)

| Item | Value | Justification |
|---|---|---|
| Frame | ECI mean-of-date; Sun/Moon rotated from ICRF/J2000 via FIXED IAU-1976 precession | Continuity with 018/019 |
| Integrator | RK4 fixed-step, dt = 60 s | 018/019 verified; design order p≈4.5 |
| Mode isolation | j2_only (control) and sun_moon_j2 (full) | Subtract for Lunisolar contribution |
| Inclinations | i ∈ {97.7876 (i_sso), 90, 30} deg | h=600 km; covers LEO prograde, J2-clean, SSO |
| Altitude | h = 600 km | Lab SSO reference |
| Snapshots | DE441 Sun+Moon, 2026-01-01 → 2045-01-01, daily, ICRF/TDB | 18.6 yr lunar nodal cycle |
| Estimators | direct OLS, secant, harmonic regression (f), node-vector OLS (n), phase-locked 2-window | 4 independent estimators per audit-020 Track 5 |
| Headline estimator | harmonic regression at 18.6-yr (Estimator f) | Theory-driven; full harmonic basis |
| Phase | Single phase (lunar anomalistic zero) | 18.6-yr direct fit averages over nodal modulation |
| Output cadence | Ascending-node crossings only (streaming) | ~103k crossings over 18.6 yr |

---

## Decision rule (pre-registered)

The mission declares the **18.6-yr harmonic regression Lunisolar rate
at i_sso** as the headline observable. The decision rule is:

- **VERIFIED-WITH-LIMITATION**: |rate_numerical − rate_cf| / |rate_cf| ≤ 0.5
- **PARTIALLY-VERIFIED (H1 marginal)**: 0.5 < ratio ≤ 2.0
- **REJECTED H1 (H0 or H-uncertainty plausible)**: ratio > 2.0 OR sign disagreement

i=90 and i=30 are inclination-structure controls; both must agree
within ±100% of their corrected-formula predictions or the result is
declared UNRESOLVED pending further investigation.

---

## Implementation summary

- **experiment.py** (the standalone mission runner) is at
  `experiment.py`. It implements the byte-pinned snapshot loader,
  IAU-1976 precession, third-body acceleration (direct + indirect),
  fixed-step RK4 with ascending-node detection, streaming propagation
  (no full-trajectory storage), the 4 estimators, and pre-flight
  oracles (synthetic test, force-level identity, idealized bridge).
- **run_focused_campaign.py** orchestrates 3 inclinations × 2 modes
  = 6 RK4 propagations over 18.6 yr each.
- **make_figures.py** generates 5 publication-quality PNG figures
  from the resulting `results/results.json`.
- **tests/** contains 12 tests covering snapshot integrity, formula
  pinning, synthetic oracle, force-level identity, phase-locked
  estimator, idealized bridge, and a post-condition test that
  becomes meaningful once `results.json` is written.

---

## Ephemeris provenance

- **Source**: NASA/JPL Horizons System (https://ssd.jpl.nasa.gov/horizons)
- **Ephemeris**: DE441
- **Time type**: TDB (Barycentric Dynamical Time)
- **Frame**: ICRF (J2000 inertial)
- **Units**: KM-S (kilometres and seconds; vectors are geocentric)
- **Acquisition**: `fetch_horizons_sun_moon_long.py` (committed
  script) calls the Horizons API with VECTORS = TABLE, RANGE = 19 yr,
  STEP = 1 d, then byte-pins the response under `reference/`.
- **Storage**: 9 files under `reference/` (1 concatenated sun, 1
  concatenated moon, 4 sun per-chunk, 4 moon per-chunk, 1 MANIFEST).
- **Pinning**: SHA-256 in MANIFEST.json + test pins
  (`f2c4f048...` for sun, `aee85099...` for moon).

---

## Force model (independent verification)

The mission runs J2-only and J2+Sun+Moon propagations. The third-body
acceleration is computed as the **geocentric direct + indirect** form:

    a_3b = μ₃ [ (r_sat − r₃) / |r_sat − r₃|³  −  r₃ / |r₃|³ ]

This is verified to **machine precision** (max_diff = 0.0 km/s²) against
two algebraically equivalent rearrangements at 50 random states
(see `test_force_level_identity_exact_at_machine_precision`).

The J2 secular rate is computed but not subtracted; it is its own
propagation that defines the J2-only control.

---

## Estimator theory

The mission's headline estimator is **theory-driven harmonic regression
(Estimator f)**: a linear regression of unwrapped RAAN onto a basis
consisting of

    {1, t, cos(2πt/Tₖ), sin(2πt/Tₖ)}   for Tₖ ∈ {6798.4, 365.24, 182.62,
                                                  121.75, 91.31, 73.05,
                                                  27.55, 14.77, 9.3067×365.24}

The slope coefficient on `t` is the secular rate; the harmonic
amplitudes absorb the periodic content. The 6798.4-d basis function is
the lunar nodal period; the 9.3067-yr basis is one lunar nodal cycle;
the 365-d and sub-annual basis functions absorb solar forcing and
sampling aliases.

Estimator (f) on a synthetic oracle with known secular + 8 harmonics
recovers the true secular to **machine precision** (|bias| < 1e-12
deg/day; see `test_synthetic_oracle_estimator_f_recovers_secular_to_
machine_precision`).

Estimator (n) (node-vector OLS on subsampled angular-momentum samples)
is the theory-INDEPENDENT kinematic observable: it computes RAAN from
the inertial angular momentum vector h = r × v, h_hat, and n_hat
without reference to the corrected formula.

The phase-locked 2-window estimator places two 100-d windows at half-
nodal-period separation; their average drift in principle cancels a
slow harmonic of period 2 × HALF_NODAL_DAYS = 6798.4 d. Verified to
**±1e-5 deg/day** on a synthetic oracle with known secular + slow
harmonic (see `test_phase_locked_synthetic_drift_recovers_known_slope`).

---

## Limitations (declared upfront)

1. **Single phase per inclination** (lunar anomalistic zero). The
   18.6-yr direct fit over a full lunar nodal cycle averages over the
   nodal modulation of the secular rate. A 4-phase ensemble would
   bound the phase dependence; budget limits this to 1 phase.
2. **No J2 × Lunisolar coupling term** in the corrected formula. The
   remaining residual (if any) may include this coupling; the mission
   documents it but does not attempt to model it.
3. **No planetary perturbations** beyond Sun + Moon + J2 + point-mass
   Earth. Higher-order geopotential (J3, J4, ...) is excluded; this is
   standard LEO practice.
4. **No atmospheric drag**. LEO SSO at h=600 km has a measurable drag
   contribution to Ω̇; the mission does not model it and the headline
   observable is the J2+Sun+Moon-Lunisolar contribution to the secular
   rate, computed by mode-subtraction (full − j2_only).
5. **Point-mass Sun, no solar radiation pressure (SRP)**. SRP at LEO
   SSO is ~10⁻⁵ deg/day on Ω̇, well below the 1.3×10⁻⁴ deg/day
   corrected formula signal; omitted for clarity.

---

## Status at session open (2026-09-02)

- Snapshot acquisition: COMPLETE (4b7dc56)
- Focused campaign scaffold: COMPLETE (809d1fe)
- Streaming-propagation rewrite: COMPLETE (223c64e)
- Mission card (this README): WRITTEN (this session)
- Results: NOT YET WRITTEN
- Figures: NOT YET WRITTEN
- Knowledge note: DRAFT (lunisolar-closure-021.md)
- 8-track audit: NOT YET RUN

The mission is the first clean scientific execution after the
adoption of `LAB_CONSTITUTION.md`; it preserves the prior scaffold
without rewriting history.

---

## Recommended next action (after this mission)

- If **VERIFIED-WITH-LIMITATION**: graduate the harmonic-regression
  + phase-locked estimator library to `src/lab_utils/estimation.py`
  (Mission 3 in POST_ROADMAP_PROBE).
- If **REJECTED H1**: spawn a follow-on mission that investigates
  the residual structure (model-order J2×Lunisolar coupling, lunar
  evection/variation, or a deeper theoretical derivation).
- If **UNRESOLVED** (decision rule marginal): add a 4-phase ensemble
  to bound the phase modulation.