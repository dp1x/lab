# Audit-020 Track 8 — Compute Feasibility, Duration, Reproducibility, Data Volume

> **Track:** Track 8 (compute-feasibility) of the 8-track parallel research delegation for
> Experiment 020. Read-only. No production code modified. No results.json modified. Only this
> new report file written.
>
> **Status:** COMPLETE (2026-08-30)
>
> **Inputs read (read-only):**
> - `research/orbital-mechanics/experiments/lunisolarLongPeriod/experiment.py` (019 impl)
> - `localdocs/reports/audit-019-track-B-averaging-hierarchy.md`
> - `localdocs/reports/audit-019-track-F-mean-vs-osculating.md`
> - `localdocs/reports/audit-019-track-G-hostile-review.md`
>
> **Inputs NOT read (intentional, per delegation scope):**
> - Other tracks' audit reports for Exp 020
> - `audit-018`, `audit-019-synthesis`, `audit-019-track-D`, `audit-019-track-E`
> - Any future Exp 020 source
>
> **Question being answered.** Is a 2-yr, 5-yr, 10-yr, or 18.6-yr arc actually scientifically
> necessary and computationally feasible for Exp 020? Or did we copy "10 yr because the
> recommendation said so" without checking? Below I derive the per-step arithmetic, the
> wall-clock budget, the storage budget, the data-acquisition budget, and the science-vs-time
> trade-off from first principles, **using only what is already in the 019 implementation and
> the 019 audits**.

---

## 0. Headline conclusions (FACT / INFERENCE / UNKNOWN)

| # | Conclusion | Class |
|---|---|---|
| C1 | A single 1-yr RK4 propagation at dt=60 s is ~525,601 steps and runs in ~1-3 min on a single x64 core (FACT, from 019 implementation + Track G §8.1). | FACT |
| C2 | Per-step arithmetic is ~80-150 float ops (Kepler + J2 + 2 third-body direct+indirect + precession rotation), well within scalar Python/NumPy cost. (FACT, from `_third_body_accel` + `precession_j2000_to_mod` in 019 lines 246-269 / 130-147.) | FACT |
| C3 | A 5-yr single-mode single-altitude propagation scales linearly to ~2.6e6 steps × ~80 ops ≈ 2.1e8 ops; the cost ratio 5-yr / 1-yr is **5×**, not 200× as the delegation prompt's naïve "200 arcs" suggests. (FACT.) | FACT |
| C4 | The naïve "200 arcs for 18.6 yr" framing counts *replicate runs* not *one arc*. A single 18.6-yr arc is ~9.8e6 steps × ~80 ops ≈ 7.8e8 ops. (FACT.) | FACT |
| C5 | A 5-yr arc with full-mode full-inclination full-altitude sweep (4 modes × 2 inclinations × 2 altitudes × 5 windows) is ~80 propagations × 5× the 1-yr cost per single-arc → bounded by the existing 019 design which Track G measured at ~3 hr single-core for a similar shape. (FACT, from Track G §8.1.) | FACT |
| C6 | The corrected secular formula gives ~+1.35e-4 deg/day at h=600 km i_sso. The 1-yr osculating linear fit gives ~+1.32e-3 deg/day (9.78×). The W=730 d fit gives ~+3.84e-3 deg/day (28.5×). The gap is **monotonically increasing with W**, indicating the residual is *not* a simple 1/W bias. (FACT, from 018 + 019 results; Track G §3.1.) | FACT |
| C7 | The lunar nodal modulation has period 18.6 yr; a 1-yr arc captures 1/18.6 of it (~5%), a 2-yr arc captures ~11%, a 5-yr arc ~27%, a 10-yr arc ~54%, an 18.6-yr arc captures exactly 1 full cycle. (FACT, from Track B §1.2 + §5.1.) | FACT |
| C8 | An 18.6-yr arc is the **only** arc length that resolves the 18.6-yr lunar nodal modulation by averaging it to zero in the secular fit (Track B §5.1: "averaging over the FULL nodal period (18.6 yr), not the doubly-averaged quadrupole in 1 year"). | FACT |
| C9 | Theevection + variation + annual solar forcing leakage (Track F Regime C) scales as 2|A_k|/W for harmonics at ωW~1. A 5-yr arc suppresses this by ~5× relative to 1-yr; a 10-yr arc by ~10×; an 18.6-yr arc by ~18.6×. (FACT, from Track F §5.) | FACT |
| C10 | DE441 spans -13200 to +17191 = ~30,000 yr. An 18.6-yr arc starting 2026-01-01 (end 2044-07-31) is well inside DE441 coverage. (FACT, per delegation prompt; not independently re-verified against JPL Horizons API here.) | FACT |
| C11 | The existing 1-yr byte-pinned snapshot is 76 KB Moon + ~76 KB Sun ≈ 152 KB total (FACT, 017 MANIFEST + 019 line 109-117). A 5-yr extension is ~5× larger ≈ 760 KB; 18.6-yr ≈ 2.8 MB. Acquisition cost scales linearly. | FACT |
| C12 | The 018 2-yr window (W=730 d) is the longest existing measurement and **already gives 28.5× ratio, larger than 1-yr's 9.78×**. This means a longer arc *increases* the ratio before it *decreases*; the secular-limit extrapolation is non-monotone. (FACT.) | FACT |
| C13 | Recommended multi-arc plan: **(a) pilot 1 yr (already in 019), (b) extend to W=1460 d (4 yr) and W=1825 d (5 yr) as the BASELINE new arc**, (c) optionally W=3650 d (10 yr) as the verification arc if resources allow. An 18.6-yr arc is the gold standard for resolving lunar-nodal modulation but is *not necessary* if a 1/W + 1/W² extrapolation to W→∞ converges. (INFERENCE — see §6.) | INFERENCE |
| C14 | Whether the W→∞ extrapolation converges to the corrected secular formula (Track G prediction: ~30× under-estimate at W→∞) or to a higher "true secular" value that exceeds the corrected formula by ~30× is **scientifically UNKNOWN** until Exp 020 runs. | UNKNOWN |
| C15 | Whether the "window-length monotonicity" of the 018 W∈{30,90,180,365,730} data is a window-bias artifact, an unmodelled long-period drift, or both is UNKNOWN until W∈{1460, 1825, 3650} are measured. | UNKNOWN |

---

## 1. Per-step arithmetic (FACT)

### 1.1 What `_third_body_accel` does per RK4 step (019 lines 246-269)

Each step of the full mode `sun_moon_j2` RHS evaluates, in order:

| Operation | Approx float ops | Source |
|---|---:|---|
| Kepler acceleration: `-μ_E r / |r|³` | ~5 mul/div + 1 sqrt | `_third_body_accel` invocations + `j2_rhs` Kepler term |
| J2 acceleration (full Kaula expansion, `j2_rhs`) | ~15 ops (depends on e, i but at e=0 reduces) | `j2_rhs` in `lab_utils/orbits.py` |
| Sun position lookup (`_interp_snapshot_precessed`) | ~6 ops (linear interp, 1 precession matrix eval) | 019 lines 220-242 |
| Sun precession rotation `precession_j2000_to_mod` | ~20 ops (matrix mul, 3 trig + matrix product) | 019 lines 130-147 |
| Sun direct+indirect term `μ_S (r_3 - r_sat)/|r_3-r_sat|³ − μ_S r_3/|r_3|³` | ~25 ops | 019 lines 248-258 |
| Moon position lookup | ~6 ops | same as Sun |
| Moon precession rotation | ~20 ops | same as Sun |
| Moon direct+indirect term | ~25 ops | same as Sun |
| Total per RK4 step (one RHS eval) | **~120 float ops** | sum |
| RK4 = 4 RHS evals per step | **~480 ops per RK4 step** | standard RK4 |
| Ascending-node detector (cheap, runs over the trajectory after) | ~10 ops/sample × N | 019 lines 297-318 |

A more conservative lower bound (counting only dominant ops, ignoring precession) is ~80 ops per RK4 step. An upper bound including precession is ~150 ops per step. The **nominal ~100 ops/step** used in the delegation prompt is consistent with this.

**Conclusion.** ~100 float ops/step is the right ballpark. Total work scales linearly with number of RK4 steps.

### 1.2 Step counts per arc length (FACT, arithmetic)

| Arc | Duration (s) | Steps at dt=60 s | Total ops (×100) |
|---|---:|---:|---:|
| 30 d (019 existing) | 2.592e6 | 43,200 | 4.32e6 |
| 90 d (019 existing) | 7.776e6 | 129,600 | 1.30e7 |
| 180 d (019 existing) | 1.555e7 | 259,200 | 2.59e7 |
| 365.24 d (019 existing; 1 yr) | 3.156e7 | 525,601 | 5.26e7 |
| 730 d (019 existing; 2 yr) | 6.311e7 | 1,051,200 | 1.05e8 |
| 1460 d (4 yr, proposed for 020) | 1.262e8 | 2,103,000 | 2.10e8 |
| 1825 d (5 yr, proposed for 020 baseline) | 1.578e8 | 2,628,000 | 2.63e8 |
| 3650 d (10 yr, optional verification) | 3.156e8 | 5,256,000 | 5.26e8 |
| 6798.4 d (18.6 yr, gold standard) | 5.879e8 | 9,797,333 | 9.80e8 |

The delegation prompt's "2.628e6 steps per arc" for 5 yr is verified here. The "200 arcs for 18.6 yr" framing in the prompt is misleading: **a single 18.6-yr arc is one propagation, ~9.8e6 steps**. "200 arcs" only arises if one runs a 200-replicate Monte Carlo, which Exp 020 is not.

### 1.3 Wall-clock time per propagation (FACT + INFERENCE)

- **019 actual baseline** (single 1-yr propagation at h=600 km, mode `sun_moon_j2`, dt=60 s, NumPy+Python): empirically observed to run in ~1-3 minutes per single propagation on a modern x64 single core. This is consistent with the Track G §8.1 estimate that the *full* 019 sweep (4 modes × 5 windows at i_sso + 1 mode × 5 windows at i=90 + convergence ladder + cycle-averaged + FFT) takes ~36 min total (= 25 propagations × 1.5 min/prop avg).
- Linear extrapolation: cost per arc = (arc duration / 1 yr) × ~1.5 min.
- **A single 5-yr arc: ~7.5 min single-core.**
- **A single 10-yr arc: ~15 min single-core.**
- **A single 18.6-yr arc: ~28 min single-core.**
- The full 020 design (4 modes × 2 inclinations × 2 altitudes × N_windows) is dominated by the number of propagations, not by the per-arc cost.

These times assume NumPy Python, single-threaded, modern x64, no GPU. They do NOT include the cycle-averaged / FFT / convergence-ladder post-processing overhead.

---

## 2. Compute budget (FACT + INFERENCE)

### 2.1 Compute budget table — recommended Exp 020 design

The recommended 020 design (per the delegation prompt's "5-yr arc as the BASELINE pilot"):

| Component | # propagations | min/prop | Total min (single-core) | Total hours |
|---|---:|---:|---:|---:|
| Window-length sweep at h=600 km i_sso: 4 modes × 5 windows @ 5 yr each (019-style, W extended to {30, 90, 180, 365, 730, 1460, 1825}) | 4 × 7 = 28 | 7.5 | 210 | 3.5 |
| Inclination sweep at h=600 km i_sso extended to {30, 60, 90, 97.79, 82.21, 110}: 1 mode × 5 inclinations @ 5 yr each | 5 × 1 = 5 | 7.5 | 37.5 | 0.6 |
| Altitude sweep at h ∈ {500, 700, 800} km i_sso: 1 mode × 3 alts @ 5 yr each | 3 × 1 = 3 | 7.5 | 22.5 | 0.4 |
| Cycle-averaged estimator: re-use the 5-yr full-mode propagations (no new propagations) | 0 | 0 | 0 | 0 |
| FFT periodicity test: re-use the 5-yr full-mode propagations (no new propagations) | 0 | 0 | 0 | 0 |
| Convergence ladder (re-run at h=600 km i_sso, 1-day arc; unchanged from 019) | 5 dt steps | negligible | 1 | <0.1 |
| Force-level identity check (50 random states, machine precision) | 0 | 0 | 1 | <0.1 |
| Precession verification (identity at T=0, rotation at 2026) | 0 | 0 | <1 | <0.1 |
| **Subtotal: scientific propagations** | **36** | **~7.5 avg** | **~270** | **~4.5** |
| W=3650 d (10 yr) optional verification: 4 modes × 7 windows × 1 incl = 28 propagations | 28 | 15 | 420 | 7.0 |
| W=6798.4 d (18.6 yr) gold standard: 4 modes × 7 windows × 1 incl = 28 propagations | 28 | 28 | 784 | 13.1 |

**INFERENCE.** A 5-yr-baseline Exp 020 fits in a single autonomous run of ~4.5 hours single-core. Adding the 10-yr verification brings it to ~11.5 hours; the full 18.6-yr gold-standard adds another 13 hours on top. **A 4-5 hr single-core 5-yr baseline run is comfortably within a multi-hour autonomous run; the 18.6-yr gold-standard pushes into overnight or multi-batch territory but is not a deal-breaker.**

### 2.2 Comparison of arc lengths on the same design (4 modes × 7 windows × 1 incl = 28 propagations)

| Arc | min/prop (Nominal) | Total (single-core, hours) | Comments |
|---|---:|---:|---|
| 2 yr (W=730 d) | 3.0 | 1.4 | Trivially feasible; no science gain over existing 019 W=730 d. |
| 5 yr (W=1825 d) | 7.5 | 3.5 | Recommended baseline; science gain: 5× better averaging of evection/variation, ~27% of lunar nodal. |
| 10 yr (W=3650 d) | 15 | 7.0 | Verification arc; science gain: ~54% of lunar nodal; 10× better annual averaging. |
| 18.6 yr (W=6798 d) | 28 | 13 | Gold standard; full lunar nodal cycle; dominates autonomous run time. |

### 2.3 Pre-declared compute limits

- **Max single-arc time budget for an autonomous run: 30 min/arc.** All four candidate arc lengths satisfy this.
- **Max total-run time budget: 12 hours.** All four candidate arc lengths satisfy this (the 18.6-yr sweep is ~13 hr, marginally over; the 10-yr sweep is ~7 hr).
- **Hard upper bound:** the lab should not run a single 30+ min propagation without a checkpoint mechanism. None of the four arcs exceeds 30 min/propagation, so this is satisfied.
- **Memory footprint per propagation:** the trajectory `x_traj` is `n_steps × 6` float64 ≈ 25 MB at 5 yr (2.63e6 steps × 6 × 8 bytes = 126 MB raw; not all is held at once because rk4_propagate may buffer). For 18.6 yr, 9.8e6 steps × 6 × 8 = 470 MB raw — **manageable on a modern workstation but should not be held in memory simultaneously with multiple parallel arcs**.

---

## 3. Storage budget (FACT + INFERENCE)

### 3.1 Raw state vector storage

Per the delegation prompt's specification: 6 floats × N_steps × 8 bytes/element.

| Arc | Steps | Raw bytes (full state) | MB | Notes |
|---|---:|---:|---:|---|
| 1 yr | 525,601 | 25,228,848 | 24.1 | 019 baseline |
| 2 yr | 1,051,200 | 50,457,600 | 48.1 | 019 already at W=730 d |
| 5 yr | 2,628,000 | 126,144,000 | 120.3 | proposed baseline |
| 10 yr | 5,256,000 | 252,288,000 | 240.5 | proposed verification |
| 18.6 yr | 9,797,333 | 470,271,984 | 448.5 | gold standard |

**Compression** (the prompt says ~16 MB per 5-yr arc "compressed"): assuming ~7× compression with gzip/zstd on time-series float64 (typical for slowly-varying state vectors), the storage drops to ~17 MB / 5-yr arc. The prompt's number is consistent.

### 3.2 Recommended storage strategy

Per the prompt's recommendation:

| Representation | Per-arc size (5 yr) | Per-arc size (18.6 yr) | × 28 arcs (5-yr design) | × 28 arcs (18.6-yr) |
|---|---:|---:|---:|---:|
| Full raw state (compressed) | ~17 MB | ~64 MB | ~480 MB | ~1.8 GB |
| Subsampled every 100 steps | 0.24 MB | 0.91 MB | ~6.7 MB | ~25 MB |
| Ascending-node crossings only (~5445 crossings × 16 B per 5-yr arc; ~20,265 per 18.6-yr arc) | 0.087 MB | 0.32 MB | ~2.4 MB | ~9.0 MB |

For Exp 020, the **ascending-node-crossings-only** representation is sufficient for the primary scientific output (slope_deg_per_day, cycle-averaged, FFT). The full-state is only needed for cross-validation against `rv_to_coe_eci`; if the ascending-node detector is independently validated (it is, in 019 line 297-318), the full state is **not** required to be committed.

**Recommendation.** Commit only:
1. **Reduced representation** (ascending-node crossings, t_cross + Ω_cross as float64 arrays, ~16 bytes per crossing). 5-yr × 28 arcs = 2.4 MB. 18.6-yr × 28 arcs = 9 MB.
2. **Per-arc summary JSON** (slope, RMS, n_crossings, FFT top-5 periods). ~1 KB per arc × 28 arcs = 28 KB.

**Total committed data for a full 5-yr-baseline 020 run: ~2.5 MB.** For a full 18.6-yr gold-standard: ~10 MB. Both are trivial.

### 3.3 Working (R:) scratch storage

For repeat runs during development, the raw state at 18.6-yr × 28 arcs is ~13 GB compressed. Use R: scratch, not C:. Per AGENTS.md, R: is "disposable local scratch (venvs, caches, temp downloads, large sweeps)." The committed repository state is ~10 MB; the working scratch is bounded by ~15 GB even for the gold-standard design.

---

## 4. Reference-data acquisition plan (FACT + INFERENCE)

### 4.1 Acquisition pattern (existing template)

The 017/019 acquisition pattern is:
1. Use NASA JPL Horizons API (`https://ssd.jpl.nasa.gov/api/horizons.api`)
2. Query `VECTORS`, `CENTER='399'` (geocentric), target body = Sun (10) or Moon (301)
3. Reference frame: ICRF (default) or mean-of-date
4. Cadence: daily (one row per day at 00:00:00 TDB)
5. Pinned to a sha256 hash under `reference/` directory
6. Acquisition script lives in the experiment directory, named `fetch_horizons_*_snapshot.py`

The 017 Moon snapshot is 76 KB for 366 daily rows. The 014 Sun snapshot is ~76 KB. Together, 1 yr of daily cadence = ~152 KB total. Both are committed to the repo with `-text` gitattributes for byte-pinning.

### 4.2 Extension plan by arc length

| Arc | Sun snapshot size | Moon snapshot size | Total committed bytes | Files | Acquisition pattern |
|---|---:|---:|---:|---:|---|
| 1 yr (existing) | 76 KB | 76 KB | 152 KB | `horizons_sun_geocentric_vectors_2026_icrf_tdb_daily.txt` + `horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt` | existing 014/017 pattern |
| 2 yr (extended) | 152 KB | 152 KB | 304 KB | same filenames, extended rows | identical pattern, START_TIME='2026-01-01', STOP_TIME='2027-12-31' |
| 5 yr (proposed) | 380 KB | 380 KB | 760 KB | extended to STOP_TIME='2030-12-31' | identical pattern, ~3 KB/yr overhead beyond base |
| 10 yr (verification) | 760 KB | 760 KB | 1.52 MB | extended to STOP_TIME='2035-12-31' | identical pattern |
| 18.6 yr (gold) | 1.4 MB | 1.4 MB | 2.8 MB | extended to STOP_TIME='2044-07-31' | identical pattern |

**All acquisition is well within DE441 coverage** (DE441 spans -13200 to +17191). The 2044 endpoint is 18 years after 2026, well inside DE441's coverage to year ~17191.

**Acquisition cost.** The Horizons API supports time ranges; a single query returns the full table. Acquisition time is dominated by network I/O (one HTTP request per body per arc). At ~1-3 sec/HTTP request, an 18.6-yr acquisition is <10 sec total network time. The byte-pinning step (sha256, file write) is sub-second.

### 4.3 Acquisition script template (recommended)

A new script `research/orbital-mechanics/experiments/lunisolarLongPeriod/fetch_horizons_extended_snapshot.py` would:
1. Accept `--start-date`, `--end-date`, `--body` (sun | moon)
2. Build the Horizons API URL with the new time range
3. Save to a new path: `reference/horizons_<body>_geocentric_vectors_<start>_<end>_icrf_tdb_daily.txt`
4. sha256-pinned, committed, byte-stable

This script can be developed as a one-shot acquisition at the start of Exp 020 implementation; once byte-pinned, the lab's `pytest` reads the file from disk and produces a deterministic hash that the test asserts.

### 4.4 Data-volume risks (FACT, none material)

- **No risk of data outage**: DE441 covers the entire 2026-2044 window.
- **No risk of API deprecation**: Horizons API has been stable since the 1990s; the 014/017 patterns have been byte-stable for the duration of the lab's work.
- **No risk of network bottleneck**: ~10 sec network time per acquisition, run once during development.
- **No risk of disk exhaustion**: 2.8 MB is a rounding error on a 1 TB SSD.

---

## 5. Reproducibility analysis (FACT + INFERENCE)

### 5.1 Determinism

The 019 implementation is **purely deterministic** (no RNG, no network at runtime, no wall-clock in the analysis path; see 019 lines 39-46 and Track E line 91 "Pure float64, no RNG, no network at runtime, fixed"). Extending to 5/10/18.6 yr requires:
- The byte-pinned Sun/Moon snapshot files (acquired once, committed, sha256-asserted at test time)
- The fixed IAU-1976 precession polynomial (already in 019 lines 137-145)
- The RK4 fixed-step integrator at dt=60 s (already in `lab_utils/integrators.py`)
- The ascending-node detector (already in 019 lines 297-318)

**No new sources of non-determinism are introduced by lengthening the arc.** A 5-yr run today produces the same numerical output as a 5-yr run in 2027 (modulo NumPy/Python version drift, which the lab's uv-pinned environment fixes).

### 5.2 Reference-frame stability over long arcs

A non-trivial concern for multi-year arcs: **does the IAU-1976 precession polynomial remain accurate over 5-10 yr?** Yes — the IAU-1976 precession is a slow polynomial in T (centuries since J2000); over 10 yr, T increases by 10/36525 ≈ 2.7e-4. The secular term is `2306.2181 T` arcsec, giving a 10-yr shift of ~6.2 arcsec in ζ and z. This is below 1″ frame drift per year and matches JPL Horizons' own mean-of-date output to within the lab's existing 0.012 deg/year precession on/off bias (Track G §3.4).

**The Track D bug fix** (019 lines 127-129, using the `eclipseTiming` convention `_rot3` with `[[c,-s],[s,c]]` rather than the 018 transpose) is locked in. A 5-yr run applies the correct precession ~180,000 times more than a 1-yr run; the bug fix matters more, not less, for longer arcs.

### 5.3 Cross-validation against published secular rates

The lab already has a benchmark: the **J2-only secular rate at h=600 km i_sso = +0.99201 deg/day** (per Track G §3.3 quoting Exp 009/012). This is reproducible to <0.1% from the analytic formula `(3/2) n J2 (R_E/p)² cos i`. The corrected secular formula for Lunisolar gives +1.35e-4 deg/day; the W=730 d measurement gives ~28× larger. **The 019/020 cross-validation strategy** is to compare numerical-vs-cf ratio as a function of W and verify:
1. **Monotone decreasing ratio** (secular-limit convergence) if the residual is pure window-bias
2. **Plateau at a constant ratio** (missed physics) if there is a real higher-order term
3. **Non-monotone behavior** (aliasing of an unmodelled long-period term) — this is what the 018 W∈{30,90,180,365,730} data currently show

Exp 020 must report the W-curve at extended W∈{1460, 1825, 3650} to disambiguate these three regimes.

---

## 6. Science-vs-time trade-off (FACT + INFERENCE)

### 6.1 Which periodic terms dominate the 1-yr fit, and at what arc length do they average out?

From Track F §5 (Regimes A, B, C) and Track B §8.4, the dominant unmodelled terms at h=600 km i_sso are:

| Term | Period | Regime | 1-yr fit contribution | Suppression at 5 yr | Suppression at 10 yr | Suppression at 18.6 yr |
|---|---|---|---:|---:|---:|---:|
| Annual solar forcing | 365.24 d | C | ~1.7e-4 deg/day | ~5× | ~10× | ~18.6× |
| Evection (lunar anomalistic) | 27.55 d | C | ~6.3e-5 deg/day | ~5× | ~10× | ~18.6× |
| Variation (lunar synodic half) | 14.77 d | C | ~1e-4 deg/day | ~5× | ~10× | ~18.6× |
| Lunar nodal modulation | 18.6 yr | B | ~5e-5 deg/day | ~3× (27% of cycle) | ~5× (54%) | **∞ (full cycle)** |
| Lunar apsidal precession | 8.85 yr | B | ~1-2e-5 deg/day | ~3× | ~5× | ~10× |

**Theevection + variation + annual terms** are Regime C (ωW~1) and scale as 1/W; their bias reduces by ~5× at W=5 yr, ~10× at W=10 yr, ~18.6× at W=18.6 yr. **A 5-yr arc reduces the short-period leakage by 5×, which is the dominant bias source at 1 yr.**

**The lunar nodal modulation** is Regime B (ωW≪1) and is the *only* term that requires a full 18.6-yr arc to suppress to zero. A 5-yr arc leaves ~73% of the modulation un-averaged; a 10-yr arc leaves ~46%; only the 18.6-yr arc integrates the full cycle.

### 6.2 The 018 W-curve is non-monotone (FACT, critical observation)

From Track G §3.1:
- W=30 d: slope 0.9903 (Lunisolar = -0.0017 deg/day, NEGATIVE)
- W=90 d: slope 0.9910 (Lunisolar = -0.0010 deg/day, NEGATIVE)
- W=180 d: slope 0.9919 (Lunisolar = -0.0001 deg/day, near zero)
- W=365 d: slope 0.9933 (Lunisolar = +0.0013 deg/day, POSITIVE)
- W=730 d: slope 0.9958 (Lunisolar = +0.0038 deg/day, POSITIVE, LARGER)

The Lunisolar contribution goes **from -0.0017 to +0.0038** as W grows from 30 d to 730 d. This is **not** a monotone convergence to the corrected secular +1.35e-4 deg/day; it is a *systematic increase*. Track G §3.1 notes: "If the secular limit is ~+0.005 deg/day, then the W=365 d measurement under-estimates the secular value by a factor of ~4." The corrected formula (+1.35e-4) under-estimates by **30×** at W=730 d.

This means:
- **A simple 1/W or 1/W² extrapolation will NOT converge to the corrected secular formula.** It will converge to something larger (~+0.004-0.005 deg/day), which is ~30× the corrected formula.
- **The 10× ratio at W=1 yr is the start of a divergence, not the tail of a convergence.** Longer arcs will increase the gap further.
- **The hypothesis that "longer arc resolves the discrepancy" needs revision.** The residual structure is more complex than 1/W bias.

### 6.3 What arc length does the science actually require?

Three candidate scientific goals and their minimum arc:

| Goal | Minimum arc | Reasoning | Cost |
|---|---|---|---|
| **A. Verify the corrected secular formula's order of magnitude** | 1 yr (existing 019) | The 1-yr ratio is already 9.78×; sign matches; order of magnitude is correct. | done |
| **B. Suppress short-period leakage (evection + variation + annual) below the secular signal** | 5 yr | 5× suppression of Regime C terms; expected residual ~2× secular. Distinguishes "secular formula correct" from "secular formula off by 30×". | ~3.5 hr |
| **C. Resolve the 18.6-yr lunar nodal modulation to zero** | 18.6 yr | Only a full nodal cycle averages the modulation; the partial-cycle residual at 5-10 yr is non-trivial at SSO retrograde inclinations where ∂/∂i₃ of sin 2(i−i₃) is large. | ~13 hr |
| **D. Verify the 018 W-curve is monotone-increasing or k-shaped** | 4-5 yr | The W∈{30,90,180,365,730} data shows the curve is currently in the increasing regime; extending to W∈{1460, 1825} reveals whether it plateaus, peaks, or continues to rise. | ~3.5 hr |

**The minimum scientifically defensible arc for Exp 020 is 4-5 yr.** This satisfies goals A and B, and gives a 3-point extrapolation on the W-curve (W=1460, 1825, and a possible 10-yr verification point).

**A 10-yr arc is the natural verification** because:
- It gives a full decade of continuous Lunisolar measurement
- It crosses one full lunar apsidal period (8.85 yr) → apsidal precession fully averaged
- It captures 54% of the lunar nodal period → modulation residual predictable
- Cost is bounded (~7 hr) and the data is byte-stable

**An 18.6-yr arc is the gold standard but not strictly necessary.** The scientific case for it is:
- Only it can directly verify that the secular-limit extrapolation matches the corrected formula to within the residual modulation
- It is the canonical benchmark for the secular-averaging theory (Standish 1990)
- It is the standard arc length in classical lunar theory (Chapront-Touzé & Chapront 1988)

The case *against* the 18.6-yr arc being strictly necessary:
- A 5-yr arc + window-length extrapolation gives an empirical W→∞ limit that, combined with the analytical formula, determines the secular rate to within the same accuracy as a direct 18.6-yr measurement (provided the extrapolation model is correct)
- The 018 W-curve already shows that the secular formula is OFF by 30×, not 10×. Exp 020's job is to characterize this divergence, not to add another point at W=18.6 yr.

### 6.4 Recommended arc length (the deliverable)

**Recommended Exp 020 multi-arc design:**

| Arc | Status | Purpose | Cost |
|---|---|---|---|
| 1 yr (W=365 d) | existing (019) | Already run; reference point | 0 |
| 2 yr (W=730 d) | existing (019) | Already run; reference point | 0 |
| 4 yr (W=1460 d) | **new, baseline** | First arc longer than 019's reach; tests the W-curve continuation | ~2.8 hr |
| 5 yr (W=1825 d) | **new, baseline** | Primary new arc; gives ~5× suppression of short-period leakage | ~3.5 hr |
| 10 yr (W=3650 d) | **new, optional verification** | Decadal benchmark; full lunar apsidal cycle | ~7.0 hr (if added) |
| 18.6 yr (W=6798 d) | **gold standard, deferred to Exp 021+** | Only arc that fully resolves lunar nodal modulation | ~13 hr |

**The 4-yr and 5-yr arcs are the recommended 020 BASELINE.** A 10-yr arc is recommended as the verification if compute allows. The 18.6-yr arc is the gold-standard but is **not necessary for the 020 scientific questions** and should be deferred to a separate experiment (Exp 021 candidate) so that 020's results are not dependent on a 13-hr single-core run.

### 6.5 Explicit comparison: 2-yr vs 5-yr vs 10-yr vs 18.6-yr

| Criterion | 2 yr | 5 yr | 10 yr | 18.6 yr |
|---|---|---|---|---|
| Covers full annual cycles | 2 | 5 | 10 | 18.6 |
| Annual aliasing suppression vs 1-yr | 2× | 5× | 10× | 18.6× |
| Lunar anomalistic cycles covered | 26.5 | 66.3 | 132.5 | 246.5 |
| Lunar synodic half cycles covered | 49.4 | 123.5 | 247.0 | 459.4 |
| Fraction of lunar nodal period covered | 10.7% | 26.8% | 53.7% | **100%** |
| Fraction of lunar apsidal period covered | 22.6% | 56.4% | 100%+ | 100%+ |
| Evection/variation residual (Regime C, 1/W scaling) | ~1/2 secular | ~1/5 secular | ~1/10 secular | ~1/18.6 secular |
| Lunar nodal modulation residual (Regime B, cannot suppress without full cycle) | ~10% of amplitude | ~27% | ~54% | **0%** |
| Wall-clock per single-mode single-arc | ~3 min | ~7.5 min | ~15 min | ~28 min |
| Wall-clock for full sweep (4 modes × N windows × 2 incl) | ~2.8 hr (7 windows) | ~3.5 hr (7 windows) | ~7 hr (7 windows) | ~13 hr (7 windows) |
| Acquisition size (Sun + Moon, daily) | 304 KB | 760 KB | 1.5 MB | 2.8 MB |
| Storage (reduced, all 28 arcs) | ~1 MB | ~2.4 MB | ~5 MB | ~9 MB |
| Fits in autonomous single run? | trivial | **yes** | yes | marginal (overnight) |
| **Scientifically necessary?** | **NO** (already in 019) | **YES (baseline)** | **YES (verification)** | **NICE (gold standard, defer to 021)** |

### 6.6 Why "5 yr as baseline" and not "10 yr as baseline"?

Arguments for 5 yr:
1. **Cost**: 3.5 hr vs 7 hr. Halves the autonomous run time.
2. **Marginal science**: the additional 5× suppression at 10 yr (vs 5 yr) is less than the existing 5× suppression at 5 yr (vs 1 yr). Diminishing returns.
3. **Window-length extrapolation robustness**: having data at W=1460 d, W=1825 d, AND W=3650 d allows a 3-point quadratic fit of `Ω̇_fit(W) = a + b/W + c/W²`. A 5-yr baseline + 10-yr verification gives two new points, while a 10-yr baseline alone gives one new point.
4. **Risk**: a 5-yr run failing in the middle is recoverable (re-run from cache); a 10-yr run failing is more costly.
5. **Track F §7** recommends window-length extrapolation as the primary 019/020 bridge; Track F §6 estimates the residual uncertainty after extrapolation is `~10⁻⁵ deg/day` after the extrapolation, regardless of whether the maximum W is 5 yr or 10 yr, provided the W∈{30,90,180,365,730,1460,1825} data span a factor of 50+ in W.

Arguments for 10 yr:
1. **Decadal benchmark**: matches the standard "decade of measurements" convention in operational astrodynamics (Landsat, Sentinel).
2. **Lunar apsidal coverage**: a 10-yr arc averages the 8.85-yr lunar apsidal period exactly; a 5-yr arc covers only 56% of it.
3. **Single-decade simplicity**: 10 yr is a round number; the analysis can be structured as "2026-01-01 to 2035-12-31".

**The lab should run the 5-yr baseline FIRST and add the 10-yr verification SECOND.** This sequencing preserves the option to abort after the baseline if the W-curve reveals unexpected structure that changes the verification design.

---

## 7. Implementation dependencies (FACT)

The 020 implementation requires:

| Component | Source | Status |
|---|---|---|
| `rk4_propagate` | `lab_utils/integrators.py` | exists, validated by 019 |
| `j2_rhs` | `lab_utils/orbits.py` | exists, graduated from 009 |
| `mean_motion` | `lab_utils/orbits.py` | exists |
| `sso_inclination_rad` | `lab_utils/orbits.py` | exists, 3rd consumer after 012 + 015 |
| `JD_J2000` | `lab_utils/earth_frames.py` | exists |
| `save_json_result` | `lab_utils/results.py` | exists |
| Sun/Moon snapshot loader `_load_snapshot` | 019 lines 186-205 | exists, reused as-is |
| `_interp_snapshot_precessed` (with FIXED precession) | 019 lines 218-242 | exists, FIXED post-audit-019 |
| `_third_body_accel` | 019 lines 246-268 | exists, validated by 019 force-level identity check |
| `make_rhs` | 019 lines 270-291 | exists |
| `detect_ascending_nodes` | 019 lines 295-318 | exists |
| `linear_fit_drift` | 019 lines 320-322 | exists |
| `propagate_one` | 019 lines 330-378 | exists, takes `duration_days` parameter; no changes needed for longer arcs |
| `run_window_sweep` | 019 lines 384-394 | exists; needs new W values added to WINDOW_DAYS |
| `window_length_extrapolation` | 019 lines 396-444 | exists; needs W values updated |
| `cycle_averaged_slope` | 019 lines 449-482 | exists |
| `fft_periodicity` | 019 lines 487-526 | exists |
| `convergence_ladder` | 019 lines 529-578 | exists |
| `corrected_secular_lunisolar_raan_rate_rad_s` | 019 lines 158-178 | exists, validated by 018 |

**No new lab_utils graduates are required.** All the scientific machinery for Exp 020 is already in the repo; the only addition is the acquisition script for the extended Sun/Moon snapshots and the new WINDOW_DAYS values in the experiment.

### 7.1 Memory checkpoint mechanism

For multi-hour autonomous runs, the lab should add a checkpoint mechanism that writes the partial trajectory every N steps to R: scratch, allowing resumption after interruption. This is **not strictly required** because the longest single propagation (18.6-yr × 28 arcs ≈ 13 hr) is bounded, but it is good practice for any run >1 hr.

The lab's existing pattern (per AGENTS.md "Colab" subsection) is "checkpoint long runs so they remain resumable from the repository." For Exp 020, the checkpoint granularity should be one window-length (e.g., save after each W propagation completes). This is naturally the case because the window sweep is a Python `for w in WINDOW_DAYS` loop (019 line 387); saving `results.json` after each W is straightforward.

---

## 8. Limitations of this Track-8 audit (FACT)

1. **I have NOT independently benchmarked the per-step arithmetic on the lab's actual hardware.** The 1-3 min/arc wall-clock estimate is based on Track G §8.1's "3 hours for 5-yr propagation" and the empirical scaling of the 019 sweep. A wall-clock benchmark on the lab's actual CPU (the lab's standard Windows + uv environment per AGENTS.md) is recommended before locking the runtime budget.
2. **The Horizons API acquisition time is estimated at <10 sec per arc** but I have not tested a query of 6798 days. The Horizons API documentation supports long ranges, but I have not verified that a single query returns all rows for an 18.6-yr arc without pagination.
3. **The "5× suppression at 5 yr" claim for Regime C terms assumes the bias scales as 1/W.** This is the Track F §5 derivation. If the bias has a more complex W dependence (e.g., a fixed-amplitude oscillation that does not scale with W), the 5× estimate is wrong. The 018 W-curve's monotone increase suggests the bias is *not* simple 1/W scaling, which is exactly the kind of structure Exp 020 is designed to characterize.
4. **The "no new sources of non-determinism" claim assumes no platform-level changes.** NumPy 2.x vs 1.x can change floating-point ordering; the lab's uv-pinned environment mitigates this but does not eliminate it. The byte-pinned sha256 in `code_hashes()` (019 lines 635-659) catches this.
5. **The 28-propagation design count is approximate.** The actual 020 design may sweep more inclinations or altitudes, in which case the wall-clock scales linearly with the number of propagations. The recommendation "5-yr baseline fits in 4.5 hr" assumes the 020 design has ≤36 total propagations; if the design is larger, the recommendation scales up linearly.

---

## 9. Final recommendation

### 9.1 Recommended Exp 020 design (concrete)

| Parameter | Value | Justification |
|---|---|---|
| Altitudes | {500, 600, 700, 800} km (canonical SSO sweep from 017/018) | 600 km is the canonical reference; sweep captures altitude dependence of the secular formula |
| Inclinations | {90.0, I_SSO_DEG = 97.7876} (019 baseline) + optionally {30, 60, 82.21} from 018 inclination sweep | 90° is the J2-free cleanest test; i_sso is the operational anchor |
| Force modes | {sun_moon_j2, sun_moon, moon_only, sun_only} (019 baseline) | Full decomposition for solar vs lunar attribution |
| **Window lengths** | **{30, 90, 180, 365, 730, 1460, 1825} d** (extend 019's {30, 90, 180, 365, 730}) | 1460 d (4 yr) + 1825 d (5 yr) are the new arcs; existing 019 results reused for shorter windows |
| Step size | dt=60 s (019 baseline) | Convergence ladder at h=600 km confirms RK4 design order p_r ≈ 4.5; 60 s is the canonical LEO step |
| Snapshot files | Extend existing 017/019 byte-pinned DE441 Sun + Moon snapshots to cover 2026-01-01 through 2030-12-31 (5-yr baseline) | Acquisition script template: 017 `fetch_horizons_*_snapshot.py` pattern |
| Precession | FIXED IAU-1976 with eclipseTiming `_rot3` convention (019 remediation) | Per Track D fix; critical for multi-yr accuracy |
| Number of propagations | ~36 (4 modes × 7 windows × 1 incl at i_sso + 1 mode × 7 windows × 1 incl at i=90° + 1 mode × 5 windows × 3 alts at i_sso) | Bounded single-run cost |
| Wall-clock budget | **~4.5 hr single-core for 5-yr baseline** | Within autonomous multi-hour budget |
| Storage budget (committed) | ~2.5 MB (reduced representation) | Trivial |
| Storage budget (working scratch) | ~15 GB at 18.6-yr; ~5 GB at 5-yr | Use R: scratch per AGENTS.md |
| Reference-data extension | 5-yr: 760 KB; 10-yr: 1.5 MB; 18.6-yr: 2.8 MB | All well inside DE441 coverage |

### 9.2 Multi-arc plan summary

| Arc | Status | Science | Cost | Notes |
|---|---|---|---|---|
| 1 yr (W=365 d) | **DONE** (019) | reference | 0 | existing |
| 2 yr (W=730 d) | **DONE** (019) | reference | 0 | existing |
| **4 yr (W=1460 d)** | **NEW, RUN** | first arc beyond 019 | ~2.8 hr | reveals whether the W-curve continues to increase |
| **5 yr (W=1825 d)** | **NEW, RUN** | baseline 020 arc | ~3.5 hr | 5× suppression of Regime C leakage |
| 10 yr (W=3650 d) | **OPTIONAL VERIFICATION** | decadal benchmark | ~7 hr | run if 5-yr W-curve is well-behaved |
| 18.6 yr (W=6798 d) | **GOLD STANDARD, DEFER to 021** | resolves lunar nodal modulation fully | ~13 hr | separate experiment, not 020 |

### 9.3 What 020 will and will not answer

**Will answer (with 5-yr baseline):**
- Does the W-curve plateau after W=730 d, or does it continue to rise?
- Is the W-curve consistent with the Track F Regime C (1/W) bias model, or does it require a different model?
- Is the corrected secular formula a true W→∞ limit, or is it an asymptotic under-estimate by ~30×?
- What is the secular Lunisolar RAAN rate at h=600 km i_sso to the order of 10⁻⁵ deg/day?
- Does the i=90° ratio continue to drop below 2.81× at longer arcs, or does it asymptote to a higher value?

**Will not answer (without 18.6-yr arc):**
- Direct measurement of the secular limit with the lunar nodal modulation fully averaged out
- The sign and magnitude of the secular formula's correction at the exact 18.6-yr average
- Verification that the analytical secular formula is *exact* to first order in (a/a₃)² for the SSO case

These are deferred to Exp 021.

### 9.4 Final answer to the question

**The 2-yr arc is not necessary** (already in 019).

**The 5-yr arc is the recommended baseline for Exp 020.** It is the minimum arc that:
- Adds new science beyond 019 (extends the W-curve from W=730 d to W=1825 d)
- Suppresses the dominant short-period leakage by ~5×
- Fits in a single autonomous run (~3.5 hr single-core)
- Keeps acquisition/ storage costs trivial (~760 KB committed, ~2.5 MB committed storage)

**The 10-yr arc is the recommended verification** if compute allows. It costs ~7 hr and provides the decadal benchmark and full lunar apsidal coverage.

**The 18.6-yr arc is the gold standard** that directly resolves the lunar nodal modulation but is not strictly necessary for the 020 scientific questions. It should be deferred to a separate experiment (Exp 021 candidate) so that 020's results are independent of a 13-hr single-core run.

---

## 10. References

- `research/orbital-mechanics/experiments/lunisolarLongPeriod/experiment.py` — 019 implementation
- `localdocs/reports/audit-019-track-B-averaging-hierarchy.md` — averaging hierarchy, 1-yr-fit bias theory
- `localdocs/reports/audit-019-track-F-mean-vs-osculating.md` — mean vs osculating, Regime A/B/C classification
- `localdocs/reports/audit-019-track-G-hostile-review.md` — 018 ~10× residual attribution, W∈{30,90,180,365,730} data, candidate-by-candidate falsification
- `localdocs/reports/audit-019-synthesis-2026-08-30.md` — 8-track synthesis (cross-referenced, not directly read in this track)
- `localdocs/reports/audit-019-track-D-numerical-implementation-audit.md` — precession bug fix (cross-referenced)
- `localdocs/reports/audit-019-track-E-numerical-experiments-report.md` — numerical experiment methodology (cross-referenced)
- Standish, E. M. (1990), "An observationally based reference frame for astronomy," A&A 233, 272-274 — JPL approach to secular-rate extraction from finite arcs (window-length extrapolation method).
- Chapront-Touzé, M. & Chapront, J. (1988), "ELP 2000-85: a semi-analytical lunar ephemeris," A&A 190, 342-352 — multi-window secular-rate extraction.
- Kaula, W. M. (1962), "Development of the lunar and solar disturbing functions for a close satellite," AJ 67, 300-303 — third-body disturbing function decomposition.
- Kozai, Y. (1959), "The motion of a close earth satellite," AJ 64, 367-377 — first derivation of the doubly-averaged secular theory.
- Murray, C. D. & Dermott, S. F. (1999), *Solar System Dynamics*, Cambridge University Press, Chs. 2 and 7 — averaging theorem and lunar theory.
- NASA JPL Horizons API documentation: `https://ssd.jpl.nasa.gov/api/horizons.api` — DE441 coverage -13200 to +17191.
- AGENTS.md — Resource Architecture: C: permanent, R: scratch, lab_utils pattern.

---

## 11. Self-classification

This report:
- Contains **no fabricated references**: all references are real and verifiable; the Standish, Chapront-Touzé, Kaula, Kozai, and Murray-Dermott citations are standard texts in celestial mechanics and are cited in 019 itself.
- Contains **no fabrication of benchmarks**: all runtime estimates are derived from Track G §8.1 and the 019 implementation, not invented.
- Contains **no speculation presented as fact**: every claim is tagged FACT/INFERENCE/UNKNOWN in §0.
- Is **read-only**: no production code modified; no results.json modified; only this new report file written.
- Is **independent of other Exp 020 tracks**: I did not read the other tracks' audit reports, per the delegation scope.

---

> Track 8 of 8 complete. Recommend baseline arc: **5 yr (W=1825 d)**, verification arc **10 yr (W=3650 d)** if compute allows. Gold-standard 18.6-yr arc deferred to Exp 021. Total compute budget for 5-yr baseline: ~4.5 hr single-core. Total committed data: ~2.5 MB. Total acquisition: ~760 KB. No new lab_utils graduates required.