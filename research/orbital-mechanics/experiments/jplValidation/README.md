# Experiment Card: JPL Ephemeris Validation (ISS / Horizons ICRF states)

> Status: complete
> Date: 2026-08-24
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/jplValidation/`

## Research Question

How well do progressively enriched deterministic laboratory force models —
M1 two-body, M2 two-body+J2, M3 two-body+J2+drag — reproduce an authoritative
external ephemeris (NASA JPL Horizons geometric ICRF vector states of the
International Space Station), and how can the resulting discrepancy be
decomposed into numerical integration error, reference sampling/interpolation
effects, time-system bookkeeping, initialization, constants/conventions, and a
jointly-attributed remainder?

**Framing commitment (pre-registered):** the reference is an SGP4/TLE-provenance
product whose own error grows ~1–3 km/day from its trajectory revision date
(disclosed in the Horizons object sheet and manual). It is an authoritative
external reference, **not** metaphysical ground truth; residuals therefore
measure model-vs-reference divergence *jointly* with reference uncertainty,
which are never separated. The experiment validates (a) internal consistency of
the model hierarchy and (b) upper bounds on combined model+reference
discrepancy — never absolute accuracy.

## Background Theory

- Two-body + $J_2$ Cowell dynamics in a geocentric pseudo-inertial frame
  (identical to Exps 006/009/012 canon); exponential-atmosphere drag with
  ballistic coefficient $\beta = m/(C_D A)$ and co-rotating atmosphere (Exp 010).
- Reference-built local orbital frame (RTN/RSW): $\hat r$, $\hat c = \widehat{r\times v}$,
  $\hat t = \hat c\times\hat r$; residuals projected at matched epochs. The transverse
  axis is perpendicular to $\hat r$, NOT velocity-aligned (differs up to 47.7° for
  e=0.74 orbits — irrelevant here but convention declared).
- Time: TDB as uniform dynamical argument on both sides (leap-free; TT−TDB ≤ 2 ms
  never mixed). The `JDTDB` column is the sole clock; calendar fields are validated
  against their row's Julian day via two independent Gregorian-ordinal pipelines
  and then discarded.
- Exact-grid alignment: our integrator output grid lands exactly on snapshot
  epochs (internal RK4 substeps between them), so time-reconstruction error of
  the reference is identically zero in all headline numbers.

## References

Real references only:

- NASA/JPL Horizons API v1.3 documentation (`ssd-api.jpl.nasa.gov/doc/horizons.html`)
  and Horizons system manual v4.98d (`ssd.jpl.nasa.gov/horizons/manual.html`) —
  query semantics, frame/time-system declarations, spacecraft-trajectory
  provenance and degradation guidance ("10s or 100's of km" beyond a few days).
- Horizons object sheet for ISS (-125544), retrieved 2026-08-23/24 UTC:
  "Trajectory is TLE-based… low accuracy more than a few days past [revision]."
  Revision at acquisition: Aug 23, 2026. Pinned verbatim under `reference/`.
- Vallado, *Fundamentals of Astrodynamics and Applications* — exponential
  atmosphere table (ATMOSEXP.DAT via CelesTrak, pinned by Exp 010), TEME/SGP4
  accuracy discussion, RTN decomposition conventions.
- Montenbruck & Gill, *Satellite Orbits* (2000) — residual decomposition practice.
- IERS leap-second table (Bulletin C-72, July 2026): no leap second through
  mid-2027; TT−UTC = 69.184 s throughout the window (context only; window is TDB).
- Laboratory canon: `src/lab_utils/orbits.py`, `src/lab_utils/integrators.py`
  (equivalence-pinned to Exps 002/006/008/009/012), Exp 009 first-order nodal
  rates, Exp 010 drag machinery and doctrine.

## Assumptions

- Verified: snapshot identity/frame/units/time-system from response header
  echoes (asserted at load); |r|,|v| plausibility; motion continuity; h-vector
  alignment; osculating inclination ≈ 51.63°; nodal regression ≈ −4.96 deg/day
  matching the Exp 009 first-order anchor within 5%.
- Plausible (declared, not tuned): ballistic coefficient band {50, 100, 200, 400}
  kg/m² from the Exp 010 sweep; primary β=100 implies $A_{eff}=m/(C_D\beta)\approx$
  1906 m² using public ISS mass ≈ 419,725 kg, C_D = 2.2 — plausible mid-range
  projected area for the ISS geometry.
- Idealization: constant-β drag with static exponential atmosphere ignores
  attitude-dependent area and thermospheric weather; J2-only gravity ignores
  third body, tesseral/zonal harmonics beyond J2, SRP, relativity (all bounded
  qualitatively in Limitations; none implemented — clean exclusion doctrine).

## Methodology

1. **One-time online acquisition** (`fetch_horizons_snapshot.py`): two serial
   GET requests ≥3 s apart (object-data sheet; then geometric ICRF/TDB state
   vectors, CENTER='500@399', REF_PLANE='FRAME', REF_SYSTEM='ICRF',
   VEC_TABLE='2', VEC_CORR='NONE', OUT_UNITS='KM-S', CSV_FORMAT='YES',
   CAL_FORMAT='BOTH', TIME_DIGITS='FRACSEC'), window 2026-Aug-24 00:00 → Aug 27
   00:00 TDB at 120 s cadence (2161 rows, 3 days starting one day after the
   trajectory revision date). Raw bytes + MANIFEST.json (SHA-256, exact query
   params, retrieval timestamps, header metadata) pinned under `reference/`;
   refuse-to-overwrite guard. All numerics afterwards run fully offline.
2. **Load-time verification** (every run): SHA-256 recomputation against the
   manifest before ANY parsing (hard fail); structural checks (SOE/EOE, exactly
   2161 rows, trailing-comma fingerprint); dual independent parse pipelines
   (strict split vs whole-row regex) compared cell-by-cell; per-column magnitude
   gates; chord/speed continuity; h-vector alignment; inclination gate;
   calendar↔JDN consistency via independent ordinal pipelines.
3. **Initialization (declared)**: x(t0) = reference state at the first snapshot
   epoch, taken verbatim for every model. No parameter is fitted against any
   t > t0 data. Row-0 residuals are identically zero by construction.
4. **Model hierarchy** on dense substep grids landing exactly on snapshot epochs:
   M1/M2 via canonical `j2_rhs` + `rk4_propagate`; M3 via the Exp 010 donor
   `propagate_3d_rk4_drag` (Vallado atmosphere verbatim, ω_atm = ω_E);
   headline dt = 120 s / 8 = 15 s; μ sensitivity variants (SGP4-heritage
   398600.8, DE440 398600.435507) run through identical code paths.
5. **Pre-registered analysis commitments** (fixed before any residual was
   computed): window/cadence as acquired; residuals ONLY at snapshot epochs;
   primary metric = RIC along-track RMS; dt-ladder BEFORE residual inspection;
   frozen constants; β band with primary member chosen by the stated area rule;
   contamination gate: second-difference jump > 100 m between consecutive
   along-track residuals rejects the window before any comparison; decision
   rule: an improvement claim requires exceeding BOTH the seeded block-bootstrap
   95% CI (1-day blocks, B=200, seed 137) AND the declared reference envelope
   (3 km/day × days-since-revision); otherwise the automatic label is
   "indistinguishable given reference uncertainty". Deviations may appear only
   as dated addenda — none occurred.

## Implementation

- Script: `experiment.py` (analysis, deterministic, offline) +
  `fetch_horizons_snapshot.py` (one-time acquisition, refuses overwrite)
- Language/runtime: Python 3.12, numpy 2.5.1, matplotlib 3.11.1 (Agg)
- Runtime: `$REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/jplValidation/experiment.py` (~25 s single core)
- Determinism: pure float64; no wall-clock input to any number; bootstrap RNG
  seeded (137); figures regenerated deterministically from recorded series
  (fixed layouts, dpi=150). `results.json` meta carries the designed-in
  `timestamp_utc` — the only field allowed to differ between runs.
- Dependencies: numpy, matplotlib only (already pinned in uv.lock); acquisition
  uses stdlib urllib/hashlib — zero dependency changes.
- Reuse: `lab_utils.orbits` / `lab_utils.integrators` canon (bit-exact-pinned);
  Exp 010 drag propagator via single-hop importlib; Exp 006 propagator loaded
  for bit-exact null testing. No scaffolding rebuilt.

## Validation Method

(`tests/test_jpl_validation.py`, 46 tests; layers L0–L14)

- L0 committed-artifact integrity: SHA-256 of both raw snapshots re-derived
  independently in pytest; manifest schema; `.gitattributes` `-text` protection
  (guards the byte-hash scheme against `core.autocrlf` corruption on fresh clones).
- L1 header identity assertions (target -125544, Earth(399) center, ICRF,
  KM-S, GEOMETRIC, TDB); wrong-object and missing-frame-declaration mutants rejected.
- L2 parser exactness incl. malformed-table mutants (extra column, truncation,
  bad floats); regex-vs-split pipeline equality over all 2161 rows.
- L3 epoch alignment: independent ordinal pipelines vs hand-derived anchors
  (1970/2000 epochs, century spans); calendar↔JDN day-level + time-of-day
  consistency; corrupted-calendar mutants rejected; monotone 120 s spacing.
- L4 plausibility gates: km→m, s→day, state-order swap, velocity sign flip,
  ecliptic-rotation (23.44° obliquity) mutants each caught by named gates.
- L5 physics anchors: osculating inclination 51.63°±0.75°; nodal regression
  −4.9613 deg/day vs first-order J2 analytic within 5% and inside the
  [−5.2, −4.7] ISS anchor band (Exp 009 lineage).
- L6 model nulls: J2-off path bit-exact vs Exp 006 donor propagator; drag β=0
  path bit-exact vs the J2 canon; β<0 rejected by the donor.
- L7 RIC construction: orthonormality, right-handedness, pure-along-track known
  answer, predicted-frame diagnostic < 2% of signal.
- L8 integration self-convergence order > 3 (same-physics differences only —
  total-vs-reference residuals plateau at the model floor and carry no order).
- L9 Hermite reconstruction within analytic bound A(ωh)⁴/384; linear-interpolant
  mutant flagged as bound violation.
- L10 predictor isolation: mutating any FUTURE reference row leaves predictions
  bit-identical (no interpolation peeking); truncation invariance (no hidden
  consumption / NaN fallback).
- L11 numeric determinism (bit-equal double runs, seeded bootstrap reproducible).
- L12 adversarial mutants: UTC-shift magnitude documented (~69.184 s × 7.66 km/s
  ≈ 530 km — labeling catastrophe quantified); tampered-snapshot hash failure;
  stale byte-count failure; pre-registration thresholds pinned.
- L13 stale-results guard: results.json carries SHA-256 of every source file it
  depends on; test asserts they match the working tree.
- L14 offline doctrine: no network imports anywhere in the analysis path.

## Results

Headline (window 2026-Aug-24 → -27 TDB, 3 days, dt=15 s, RIC from reference):

| Model | RMS 3-D [km] | Transverse RMS [km] | Trend c1 [km/day] | Verdict |
|---|---|---|---|---|
| M1 two-body | 1346.52 | 1210.31 | +690.02 | baseline |
| M2 +J2 | 8.216 | 8.068 | −7.08 | skill 0.9933 (CI 0.9933–0.9979) → improvement exceeds CI & envelope |
| M3 β=100 (primary) | 26.506 | 26.460 | +18.24 | skill −2.279 (CI −6.42…−2.32) → **worse than previous model** |
| M3 β=200 | 9.808 | 9.685 | +5.57 | worse than M2 |
| M3 β=400 | 3.133 | 2.723 | −0.76 | better than M2 (band-edge observation, see below) |
| M3 β=50 | 60.496 | 60.471 | +43.64 | much worse |
| M2 with SGP4-heritage GM | 6.612 | 6.427 | −5.89 | constants sensitivity, not adopted |

Findings (full labeled set in `results/results.json`):

1. **FINDING**: J2 removes 99.33% of the residual RMS. M1's +690 km/day
   in-track drift reflects the purely-Keplerian model's constant mean-motion
   mismatch — the mean-vs-osculating semimajor-axis offset at the
   initialization epoch (J2 short-period amplitude ~9 km class) plus the
   absent J2 secular rates — both removed by the M1→M2 enrichment.
2. **FINDING (negative, reported verbatim)**: with the declared primary β=100
   and the untuned Vallado atmosphere, drag WORSENS agreement (robust: CI
   excludes zero on the negative side). No parameter was retuned to hide this.
3. **INFERENCE**: the β-band responds monotonically (c1: +43.6 → +18.2 → +5.6 →
   −0.76 km/day across 50→400), crossing zero only at the band edge (β=400 →
   RMS 3.13 km, 2.6× better than M2). The drag signature IS detectable in the
   reference, but the compatible effective β lies at/beyond the pre-declared
   band edge. Extending the band or selecting the nicest member post hoc would
   be tuning and is deliberately not done; a refined-β study must be a separate
   declared experiment.
4. **INFERENCE**: M2's leftover −7.08 km/day drift reduces to −5.89 km/day
   under SGP4-heritage GM (constants absorption) and corresponds to the
   reference behaving as if its mean semilatus decays ~25–40 m/day relative to
   our non-decaying model — consistent with TLE-B*-absorbed drag plus
   density weather; not independently separable from reference uncertainty.
5. Error budget (all bounded BEFORE interpretation): integration ≤ 19.6 m
   (self-convergence order p≈4, Richardson-corrected bound smaller);
   reference interpolation: diagnostic only, 58.03 km measured vs 59.20 km
   analytic bound at worst subsampled cadence (20 min) — exact-grid alignment
   makes it identically zero in headlines; epoch-tag uncertainty ≤ 43 μs
   → ≤ 0.33 mm along-track charge; initialization round-trip 13 nm;
   frame-origin diagnostic 0.93 m; GM variants ±0.03 km RMS effect
   (DE440) / −1.6 km (WGS72). Remainder attributed jointly to external-reference
   uncertainty (TLE fit, SGP4 force truncation, maneuvers) + unmodelled physics;
   NOT separated — the reference's own error is not independently characterized.
6. Contamination gate cleared: max second-difference 41.7 m < 100 m threshold
   (no maneuver/TLE-handover flag inside the window).

Figures (`results/figures/`, one claim each):

- `f1_residual_hierarchy.png` — |Δr|(t) for M1/M2/M3(primary): hierarchy collapse.
- `f2_ric_components.png` — M2 residual structure in the reference-built RTN frame.
- `f3_dt_convergence.png` — self-convergence ladder with measured order.
- `f4_interpolation_bound.png` — Hermite reconstruction vs analytic bound.
- `f5_beta_gm_sensitivity.png` — β band + GM variant trends vs declared envelope.

## Limitations

- Single target (ISS), single 3-day window immediately after one revision date;
  conclusions are window-scoped. A Molniya-class secondary anchor was planned
  but Horizons holds NO stored Molniya trajectories (major-body catalog checked);
  the user-TLE ingest path would change the provenance chain to a differential-
  SGP4 experiment — deferred as future work rather than half-implemented.
- The reference is TLE/SGP4-provenance with kilometer-per-day-scale uncertainty;
  nothing here measures absolute propagator accuracy, only divergence from this
  specific external product.
- Constant-β drag + static exponential atmosphere idealize attitude-dependent
  area and thermospheric weather; the β=400 band-edge observation is NOT a
  calibrated ballistic coefficient.
- Unmodelled physics (third body, harmonics beyond J2, SRP, relativity) remains
  jointly attributed with reference uncertainty.
- TEME→ICRF conversion inside Horizons is undocumented; bounded only as part of
  the remainder.
- Bit-exactness claims are same-platform/same-version (numpy 2.5.1, Python
  3.12.13); cross-platform equivalence should use ~1e-12 relative tolerance.

## Future Improvements

- Declared follow-up: refined effective-drag study with published area models
  and space-weather covariates (separate experiment, disjoint window).
- Differential arm: pin CelesTrak TLEs, request Horizons COMMAND='TLE' vectors
  (their SGP4/SDP4) and compare against our own SGP4-free models on identical
  elements — isolates the reference generator from the reference data.
- Molniya-class anchor via the user-TLE path (third-body-dominated regime).
- Lunisolar third-body term as M4 once a second consumer justifies graduating it.

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command: `$REPO_ROOT\.venv\Scripts\python.exe -m pytest && $REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/jplValidation/experiment.py`
  (generic `uv sync && uv run pytest && uv run python …` when uv is on PATH).
- Acquisition is one-time and idempotence-guarded; rerunning the experiment
  never touches the network. Byte-reproduction is defined by the COMMITTED
  snapshot files (`-text` in `.gitattributes`; hashes pinned in MANIFEST.json),
  not by re-querying Horizons (whose TLE-based trajectories change continuously).
- Deterministic regeneration verified: two independent runs produce identical
  `results.results` payloads (only `meta.timestamp_utc` differs); figure MD5s
  identical across runs. results.json additionally pins SHA-256 of every source
  file the result depends on (stale-run guard, asserted by L13).
- Provenance: `results/results.json:provenance.reference_snapshot_files` echoes
  the pinned file hashes; `code_sha256` binds result to code.
