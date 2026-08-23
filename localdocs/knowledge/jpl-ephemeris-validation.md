---
tags: [orbital-mechanics, validation, ephemeris, JPL, horizons, reference-data]
date: 2026-08-24
aliases: [JPL Validation, Exp 013, Horizons Validation]
links: [[j2-precession]], [[orbit-decay]], [[orbit-classes]], [[ground-tracks]]
---

# JPL Ephemeris Validation (Exp 013)

## Summary

First external-answer-key experiment of the flagship: the deterministic
two-body / +J2 / +J2+drag hierarchy was compared against pinned NASA JPL
Horizons geometric ICRF/TDB vector states of the ISS (-125544) over a
3-day window starting one day after the trajectory revision date. J2 removes
99.33% of the residual RMS (skill CI excludes zero); the declared drag tier
with untuned parameters **worsens** agreement at primary β=100 kg/m² — reported
verbatim — while the pre-declared β-band responds monotonically and crosses
zero only at the band edge (β=400 → 3.13 km RMS vs M2's 8.22 km). The band-edge
observation is documented without tuning: refining β would be a separate
declared experiment.

## Content

- **Reference provenance discipline**: Horizons satellite trajectories are
  TLE/SGP4-based predicts whose error grows ~1–3 km/day from the trajectory
  revision date (disclosed in the object sheet itself). The snapshot was pinned
  byte-for-byte with SHA-256 inside the repository (`reference/` +
  `MANIFEST.json`), marked `-text` in `.gitattributes` because `core.autocrlf`
  would otherwise corrupt the hashes on fresh clones. All numerics run offline;
  re-querying Horizons does NOT reproduce the data (trajectories update
  continuously) — reproduction means the committed bytes.
- **Exact-grid alignment**: propagator substeps land exactly on snapshot
  epochs, making reference time-reconstruction error identically zero in all
  headline numbers. Hermite interpolation through published states+velocities
  stays within A(ωh)⁴/384 even at 20-min cadence (~59 km bound) — useful for
  figures, catastrophic if fed into residuals.
- **Time systems**: vector tables accept TDB (used here; leap-free uniform
  dynamical scale); `JDTDB` is the sole clock; calendar columns validated via
  independent Gregorian-ordinal pipelines then discarded. TT−UTC = 69.184 s
  would fake ~530 km along-track error if scales were mixed.
- **Error decomposition before interpretation**: integration ≤19.6 m
  (self-convergence order ≈4; total-vs-reference residuals plateau at the model
  floor and carry NO order information — Exp 009 doctrine rediscovered),
  initialization round-trip 13 nm, frame-origin diagnostic 0.93 m, epoch-tag
  charge ≤0.33 mm; remainder attributed jointly to reference uncertainty +
  unmodelled physics (never separated).
- **Negative result doctrine**: the drag tier's worsening (robust CI) is a
  finding, not an embarrassment: M2's leftover −7 km/day drift behaves like the
  reference decaying ~25–40 m/day relative to our model (TLE-B*-absorbed drag);
  the declared Vallado atmosphere at β=100 over-decays (~120 m/day class),
  overshooting forward. Monotone band response localizes compatibility near
  β≈300–500 kg/m² without any post-hoc selection.
- **Mean-vs-osculating trap**: seeding pure Kepler from the osculating state
  bakes in a constant mean-motion offset (short-period a amplitude ~9 km class)
  that shows up as a huge linear transverse trend (+690 km/day observed).
  Enriched models integrating the actual perturbation dynamics don't suffer
  this; element-seeding experiments must account for it.
- **Anti-overfitting mechanics that worked**: pre-registered window/cadence/
  metric/dt/constants/decision-rule; contamination gate (100 m second-difference
  jump — max observed 41.7 m, window clean); seeded block bootstrap for skill
  CIs; code-hash binding of results.json to source (stale-run guard caught a
  real mid-experiment edit during development).

## Source Experiments

* `research/orbital-mechanics/experiments/jplValidation/` — Horizons-pinned
  validation hierarchy with decomposition budget; runnable:
  `$REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/jplValidation/experiment.py`
* [[j2-precession]] (Exp 009) — first-order nodal-rate anchor reused as an
  independent physics check on the reference itself (measured −4.9613 deg/day).
* [[orbit-decay]] (Exp 010) — drag machinery, atmosphere table, β-band sweep.
* [[orbit-classes]] (Exp 012) — graduated canon consumed via `src/lab_utils`.
* [[ground-tracks]] (Exp 008) — donor lineage of the Kepler/element machinery.

## Key Takeaways

1. External answer keys must be machine-pinned bytes, never re-fetchable URLs.
2. An authoritative reference with known degradation bounds what can be claimed:
   build the envelope into the decision rule BEFORE looking at residuals.
3. Negative and band-edge results are results; pre-registration converts them
   from embarrassments into measurements.
4. Self-convergence measures integrator order; residual-vs-reference plateaus do not.

## See Also

- `localdocs/roadmap.md` Phase 2 row 013; next: 014 eclipse timing / launch windows.
- `research/orbital-mechanics/experiments/jplValidation/README.md` (full card,
  pre-registration statement, budget tables).
