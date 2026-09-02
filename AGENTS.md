# AGENTS.md â€” Research Lab Operating Manual

Operating manual for research agents working in the **Research Lab** repository.
The full charter lives at `localdocs/charter.md`; the experiment sequence at
`localdocs/roadmap.md`.

## Purpose

Run reproducible computational research: turn questions into deterministic
implementations, tests, validated results, durable documentation, and reusable
knowledge.

## Repo Structure

```
$REPO_ROOT
â”œâ”€â”€ AGENTS.md              this manual
â”œâ”€â”€ README.md              lab overview + index
â”œâ”€â”€ LICENSE                MIT
â”œâ”€â”€ .gitignore
â”œâ”€â”€ pyproject.toml / uv.lock   Python env + pinned deps (uv)
â”œâ”€â”€ localdocs/             LAB'S OWN WRITING (charter, roadmap, templates, knowledge)
â”œâ”€â”€ webdocs/               external web material (public sources only)
â”œâ”€â”€ research/<domain>/experiments/<name>/   one self-contained experiment
â”‚   â”œâ”€â”€ README.md          experiment card
â”‚   â”œâ”€â”€ experiment.py      runnable implementation
â”‚   â”œâ”€â”€ tests/             validation tests (pytest)
â”‚   â””â”€â”€ results/           outputs, figures (JSON, committed)
â”œâ”€â”€ src/lab_utils/         shared reusable utilities (metrics, I/O, validation)
â”œâ”€â”€ tools/                 lab-level scripts (scaffolder)
â””â”€â”€ data/                  shared datasets (gitignored; prefer scratch for large data)
```

## The Loop

```
Question â†’ Theory â†’ Design â†’ Implement (experiment.py, deterministic)
  â†’ Test (pytest, green BEFORE trusting) â†’ Run â†’ Validate (analytic solutions,
    invariants, published data) â†’ Document (card + figures + results.json)
  â†’ Knowledge (localdocs/knowledge/<topic>.md) â†’ Commit â†’ Next Question
```

## Operating Rules

1. **Deterministic only.** Fixed seeds/params; no time-dependent nondeterminism; no
   unattributed speculation.
2. **Verify before trusting.** Tests must pass before results are recorded. Validate
   via analytic solutions, conservation laws, published data, or benchmarks â€”
   "if wrong, reality reveals it."
3. **Reproduce everything.** Follow the experiment template exactly: question â†’
   theory â†’ references â†’ assumptions â†’ methodology â†’ implementation â†’ validation â†’
   results â†’ limitations â†’ next question.
4. **Never fabricate** papers, results, benchmarks, or citations. Separate known
   facts / assumptions / hypotheses / results. Cite real references.
5. **No proprietary content.** Never commit API keys, tokens, passwords, personal
   info, or private data. Lab material is public-sourced or our own results.
6. **Resources are precious.** SSD/TBW, RAM, CPU. Avoid rebuilds, huge temp files,
   abandoned processes. Check scratch capacity before large work; clean up after use.
7. **Decisions are justified.** Prefer existing code/infrastructure. Complexity must
   justify itself.
8. **Purposeful delegation.** Delegate independent work only when it improves quality,
   throughput, or verification. Not for its own sake. Never edit the same file in
   parallel; the lead agent integrates. Avoid Goodhart â€” no agents spawned merely to
   satisfy a delegation metric.
9. **Documentation is memory.** Every experiment and meaningful change writes durable
   artifacts (card, results, knowledge note). Update cards when results change.
10. **Goodhart's Law.** Optimize research value, not activity. Do not create experiments
    merely to raise the count, spawn agents merely to delegate, generate datasets
    merely because storage exists, or use compute merely because it is free.

## Responsible Web Access

The lab may fetch public research material and public datasets. Preferred order:
**official API â†’ official downloadable dataset â†’ normal public webpage â†’ browser
automation when necessary.** Browser automation is a research instrument, not a
scraping weapon.

Hard boundary â€” **never**: bypass CAPTCHAs, defeat access controls, circumvent rate
limits, rotate identities to evade blocks, spoof fingerprints for evasion, hammer
endpoints, parallel-request-burst public sites, or repeatedly retry a blocked
endpoint. If a site presents an antibot/CAPTCHA/login/verification gate, **STOP and
report it to the human** â€” do not attempt to defeat it. Use human-use-rate browsing,
cache acquired material, and respect published rate limits.

Search: prefer DuckDuckGo or a reputable SearXNG/Startpage instance for general
discovery; do not use Google/Bing for ordinary search. (Accessing an actual Google-
hosted public document/dataset is fine when it is the source.)

## Resource Architecture

- **`C:`** â€” permanent state (source, docs, tests, committed results). Keep lean.
- **`R:`** â€” disposable local scratch (venvs, caches, temp downloads, large sweeps).
  Check free space at runtime before large operations; never hard-code capacity.
- **Durable storage.** The repository root is the sole durable record: all
  reproducible source, knowledge, tests, contracts, provenance, results, figures,
  acquisition recipes, checksums, and valuable recovery state live under it.
  `R:` is never a runtime dependency of reproducibility and never the only
  location of anything worthwhile; if autonomous work is interrupted or blocked,
  recovery/handoff state is written under the repository (e.g. `AUTONOMOUS_HANDOFF_<exp>.md`).
- **Colab** â€” optional, ephemeral remote compute only for workloads that would
  otherwise cause excessive local SSD/CPU pressure. Hardware-as-available (no assumed
  GPU/CPU); never evade Colab limits; checkpoint long runs so they remain resumable
  from the repository. R: â†’ Colab â†’ results â†’ R: â†’ download â†’ destroy runtime.

## Environment

- Python managed by **uv** â€” `uv sync` first, run with `uv run python ...`.
- All tests: `uv run pytest`.
- Add deps with `uv add <pkg>`; update affected docs.
- Sweep data: prefer R: scratch for large raw outputs; commit only compact summary
  JSON + figures to `results/`.

## New Experiment Checklist

1. `research/<domain>/experiments/<descriptive_name>/`
2. Copy `localdocs/templates/experiment_template.md` â†’ README.md, fill the card.
3. Write `experiment.py` (deterministic).
4. Write tests; `uv run pytest` green.
5. Run â†’ write `results/results.json` + figures.
6. Record results + limitations in the card.
7. Write an Obsidian note in `localdocs/knowledge/` (link to prior notes).
8. `git add` + commit (message describes what was verified).

## Remote-State Safety

- **Before any automated push:** verify the live remote tip immediately beforehand
  (`git ls-remote origin <branch>` or equivalent). Never assume the local
  `origin/<branch>` ref is current â€” it can be stale.
- Compare the live remote tip against the intended push target. If the remote changed
  unexpectedly, **STOP and reconcile**; never overwrite another agent's or the owner's
  work blindly.
- For history-rewriting operations, use `--force-with-lease=<expected-remote-tip>` only;
  never plain `--force`. After pushing, re-verify: remote tip == local HEAD, clean tree,
  and (where applicable) GitHub signature/verification state.
- Commit only from a verified-clean tree; never commit unrelated files, and keep the
  canonical Git identity (no repo-local identity overrides).

## Current Priority

Orbital-mechanics flagship: experiments 001â€“010 are complete (numerics
foundation, Kepler validation, Kepler solvers, Hohmann transfer, bi-elliptic
vs Hohmann crossover, combined transfer + plane change, gravity assist,
ground tracks, J2 precession, orbit decay). Experiment 006 was adversarially
audited + closed (2026-08-17); a synthesis report for 001â€“006 is in
`localdocs/reports/`. Experiment 008 (ground tracks, spherical-Earth) is
COMPLETE (2026-08-21) in `research/orbital-mechanics/experiments/groundtracks/`;
Experiment 009 (J2 precession, secular nodal/apsidal rates with numerical
validation) is COMPLETE (2026-08-22) in
`research/orbital-mechanics/experiments/j2Precession/`; Experiment 010 (orbit
decay / atmospheric drag â€” first non-conservative force: dissipation accounting
+ monotonicity doctrine, erfi/quadrature oracles, structural scalings,
co-rotation twins, J2 mean-element transient, reentry timing) is COMPLETE
(2026-08-22) in `research/orbital-mechanics/experiments/orbitDecay/`.
Next: 011 Lagrange points per `localdocs/roadmap.md`, reusing `src/lab_utils/`
and the Exp 002/006/008/009 propagator + element machinery â€” do not rebuild
scaffolding; graduating shared machinery to `src/lab_utils/orbits.py` is now
justified if it stays non-blocking.
Experiment 011 (Lagrange points / CR3BP â€” first rotating-frame experiment:
equilibria, Jacobi integral, Routh stability, nonlinear perturbation signatures,
dimensional cross-check, adversarial mutant battery) is COMPLETE (2026-08-22) in
`research/orbital-mechanics/experiments/lagrangePoints/`; shared machinery is now
graduated: `src/lab_utils/integrators.py` (generic rk4_step/rk4_propagate) and
`src/lab_utils/orbits.py` (element/Kepler canon, equivalence-pinned vs donors) â€”
use these instead of cloning per-experiment copies; CR3BP-specific code remains
experiment-local until a second consumer appears.
Next: 012 orbit classes per `localdocs/roadmap.md`, reusing
`src/lab_utils/orbits.py` + `src/lab_utils/integrators.py`.
Experiment 012 (orbit classes â€” constraint-defined families: SSO inclination lock +
finite existence boundary a_max = 12352.505 km, Molniya apsidal freeze + semi-synchronous
resonance + dwell geometry, GEO 1:1 fixed point with nonzero-rate negative control,
GTO budgets anchored to Exp 004, adversarial convention battery; finding: small-divisor
short-period dynamics near the critical inclination give a measured +323 s/orbit
Kepler-period excess) is COMPLETE (2026-08-23) in
`research/orbital-mechanics/experiments/orbitClasses/`; `j2_rhs` graduated into
`src/lab_utils/orbits.py` (second consumer after Exp 009, equivalence-pinned).
Next: 013 JPL ephemeris validation per `localdocs/roadmap.md`, reusing
`src/lab_utils/orbits.py` + `src/lab_utils/integrators.py`.
Experiment 013 (JPL ephemeris validation â€” ISS (-125544) vs pinned NASA/JPL Horizons
ICRF/TDB geometric states over a 3-day window; byte-pinned snapshot under the repo
(`-text` gitattributes), offline-deterministic analysis, exact-grid alignment,
reference-built RTN residuals, pre-registered decision rules) is COMPLETE (2026-08-24)
in `research/orbital-mechanics/experiments/jplValidation/`: J2 removes 99.33% of
residual RMS (bootstrap CI excludes zero); the declared drag tier WORSENS agreement at
primary beta=100 and the pre-declared beta band crosses zero only at its edge
(beta=400 -> 3.13 km vs M2 8.22 km) â€” documented verbatim without tuning; error budget
bounds integration/interpolation/time/init/constants before joint attribution of the
remainder to reference uncertainty + unmodelled physics.
Next: 014 eclipse timing / launch windows per `localdocs/roadmap.md`, reusing
`src/lab_utils/orbits.py` + `src/lab_utils/integrators.py`; declared follow-up
candidates from Exp 013 include a refined effective-drag study (separate experiment)
and a differential SGP4 arm via Horizons COMMAND='TLE'.

Experiment 014 (eclipse timing / launch windows) is COMPLETE (2026-08-28) in
`research/orbital-mechanics/experiments/eclipseTiming/`: conical (primary) + cylindrical
(control) shadow geometry via dual algebraic formulations, closed-form event finder on
analytic Kepler states (event error decouples from integration step), pinned-ISS arm
(first 4 event epochs agree to 5.5-13.5 s vs real NASA trajectory; 3-day tail drift to
308 s = TLE/SGP4 reference envelope), Sun model gated against byte-pinned 2026 Horizons
snapshot to 0.65 deg; 40 new tests, 525 total.

Experiment 015 (dawn-dusk SSO launch-window targeting) is COMPLETE (2026-08-29) in
`research/orbital-mechanics/experiments/dawnDuskSSO/`: first end-to-end multi-constraint
mission analysis. Composes the SSO inclination lock (Exp 012), LST-at-ascending-node
condition (a NEW quantity, not carried from Exp 014), first-order J2 secular nodal drift
(Exp 009/012), eclipse event-finder (Exp 014), and the lab's analytic Sun model + GMST
polynomial. Year-long feasible launch-time search for dawn-dusk SSO at h in {500, 600,
700, 800} km from Eastern Range; 266-295 connected components per altitude (monotone
in h); total feasible width 710 h at h=600; held-out equinox weeks dominate feasibility
(36.7 vs 11.9/day main, equinoxes are the most eclipse-favorable for h=600); 6 figures;
34 new tests, 581 total (547 baseline + 34 new). Shared machinery graduated to
`src/lab_utils`: `sso_inclination_rad` (3rd consumer after Exp 012 + Exp 014-implicit
+ Exp 015) in `lab_utils/orbits`; `gmst_rad_iau1982`, `sun_unit_and_dist_km`,
`subsolar_lon_rad`, `eci_to_ecef`, `ecef_to_latlon`, `spherical_trig_latlon`,
`lst_at_node_hours`, `node_lon_from_raan_gmst` (2nd consumer after Exp 014 for the
Sun/GMST, 2nd consumer after Exp 008 for the ECI-to-lat/lon layer) in the new
`lab_utils/earth_frames` module.

REMEDIATED 2026-08-29: an 8-track independent audit retracted the originally-published
"LST at the ascending node drifts through 24 h/year at the sidereal-solar differential
(4 min/day)" claim as RED (frame/convention error; the formula subtracted an inertial
RAAN rate from an ECEF subsolar rate and confused Earth's sidereal rotation rate
360.9856 deg/day with the SSO nodal rate ~0.9856 deg/day). The correct physics: the
LST at the orbit-plane ascending node of a true dawn-dusk SSO is approximately
CONSTANT, oscillating only with the equation-of-time envelope (~+/-12 min, ~24 min
peak-to-peak, periodic not secular). The structural findings (cardinality, equinox
dominance, sensitivity matrix, i_SSO anchors) are unchanged. See
`localdocs/reports/audit-015-*.md` for the 8 audit reports (LST/J2 derivation,
implementation audit, numerical falsifier, adversarial review, equivalence chain,
follow-up candidates, portfolio, literature cross-check).

Experiment 016 (SSO LST-drift correction) is COMPLETE (2026-08-30) in
`research/orbital-mechanics/experiments/lstDrift/`: first-principles derivation of the
actual LST drift at the orbit-plane ascending node of a true dawn-dusk SSO, decomposed
into EoT envelope (periodic ~30 min peak-to-peak, validated against byte-pinned 2026
Horizons Sun snapshot to 0.056 deg within Exp 014 0.7 deg gate) + J2 closure residual
(~2.2 deg/year, consistent with Exp 012) + Lunisolar upper-bound closed-form
(over-estimates ~50x at SSO retrograde inclinations due to large sin^2(i_SS) and
evection terms not captured by the secular average; reported as conservative ceiling
for transparency) + SRP (~mdeg/day for A/m = 0.01) + drag (exponential atmosphere,
altitude-dependent) + closed-form RAAN-control Î”v budget at the line of nodes
(Vallado 8.5). Headline: total LST drift range [no-LS, full-LS-upper] =
[~0, ~310] min/year at h=600 km; operational envelope (Sentinel-1 ~15 m/s/yr,
Landsat-7/8 ~5-15 m/s/yr) implies the real rate is much smaller than the closed-form
upper bound (~10,000 m/s/yr). 4 figures (EoT envelope, drift decomposition,
station-keeping Î”v range, orbit-plane LST year sweep); 40 new tests, 624 total repo
tests (584 baseline + 40 new). Builds on Exp 014 byte-pinned Sun snapshot + Exp 012
J2 closure + Exp 009 nodal rate formula. The remediation contract was: provide
defensible first-principles derivation of the LST drift rate and station-keeping
budget that Exp 015 claimed but did not derive; this experiment satisfies it.

Experiment 017 (Lunisolar upper-bound verification) is COMPLETE (2026-08-30) and
REMEDIATED (2026-08-30, audit-018) in
`research/orbital-mechanics/experiments/lunisolarVerification/`: byte-pinned JPL
Horizons DE441 geocentric Moon vectors (76 KB, 366 daily rows, sha256
`65f1d67f798a3b95...`) under `reference/`, fetched via identical pattern to the
Exp 014 Sun snapshot. Numerical integration of Kepler + J2 + point-mass Sun +
Moon at h in {500, 600, 700, 800} km over 1 year, with J2-only control subtraction
to isolate the Lunisolar contribution (model-order separation per Track F
Pillar C). Original 017 headline: the closed-form secular-average Lunisolar RAAN
upper bound (Vallado Eq. 9-46 form, Exp 016 model_note) over-estimated the
numerically integrated Lunisolar RAAN rate by a SIGNED RATIO of ~170x at h=600
km (cf retrograde -0.218 deg/day, numerical prograde +0.001284 deg/day). RK4
self-convergence order p_r = 4.49, p_v = 4.50. 4 figures; 32 new tests, 658
total repo tests (626 baseline + 32 new). 11 additional tests added in
remediation commit (L7 corrected formula validation, 669 total). The original
Exp 017 decadal direction was rejected by an eight-track audit as not
scientifically defensible at this time; the closed-form upper-bound verification
(audit-015 candidate #4) is what was executed.

REMEDIATED 2026-08-30 (audit-018): the 8-track independent investigation
identified the 016/017 closed-form as MATHEMATICALLY WRONG in three compounded
ways: (1) wrong radial scale factor (J2-style `(R_E/r_3)^2` instead of the
third-body `(a/a_3)^3`); (2) wrong geometric factor (Kozai APSIDAL
`cos(i) (1-5/2 sin^2(i-i_3))` instead of the NODAL `sin 2(i-i_3) / sin i`);
(3) wrong sign at SSO retrograde. The CORRECT formula is `(3/8) n
(mu_3/mu_E) (a/a_3)^3 sin 2(i-i_3) / sin i` (Track B independent derivation);
at h=600 km i_sso=97.79 deg it gives +1.35e-4 deg/day (prograde, SAME SIGN as
numerical +1.32e-3 deg/day, 9.78x smaller magnitude). The 10x residual is the
unmodelled short-period contribution (evection + variation + lunar-nodal). The
wrong formula is preserved as `closed_form_lunisolar_raan_rate_rad_s` (017) and
`luni_solar_raan_rate_rad_s` (016) with DeprecationWarning for backwards
compatibility; the corrected formula is exposed as
`corrected_secular_lunisolar_raan_rate_rad_s` (017) and
`corrected_luni_solar_raan_rate_rad_s` (016). 016 LST-drift budget impact: the
~310 min/year full-LS upper bound is wrong; the corrected formula gives ~1620x
smaller magnitude in the OPPOSITE direction; the operational Sentinel-1
(~15 m/s/yr) and Landsat (~5-15 m/s/yr) budgets remain the empirical ground
truth and are consistent with the corrected formula, NOT the 016/017
closed-form.

Experiment 018 (Lunisolar RAAN reconciliation) is COMPLETE (2026-08-30) in
`research/orbital-mechanics/experiments/lunisolarReconciliation/`: builds on
Exp 017 byte-pinned Moon snapshot + Exp 014 Sun snapshot + 8-track audit
synthesis. Implements the corrected secular formula `(3/8) n (mu_3/mu_E)
(a/a_3)^3 sin 2(i-i_3) / sin i` and runs controlled numerical experiments:
force isolation (j2_only / sun_only / moon_only / sun_moon / sun_moon_j2 at
h=600 km i_sso), inclination sweep (i in {0, 30, 60, 90, 97.79, 82.21} deg at
h=600 km), window-length sensitivity (W in {30, 90, 180, 365, 730} d),
precession on/off (with and without IAU-1976 precession rotation applied
to Sun/Moon vectors), force-level identity (50 random states, machine
precision), and dt convergence ladder (RK4 design order confirmed).
HEADLINE: the 170x signed discrepancy is RESOLVED; the corrected formula
agrees with the numerical in SIGN (both prograde) and within 9.78x in
magnitude at h=600 km i_sso. The CLEANEST test (i=90 deg, where J2 cos(i)=0)
gives 2.81x agreement, confirming the residual is dominated by unmodelled
short-period terms. The 016 frame-mismatch (Track D) is small: 0.012 deg/year
bias from not applying IAU-1976 precession. 6 figures (corrected cf vs
numerical, inclination sweep, window sensitivity, precession comparison,
convergence ladder, Lunisolar decomposition); 45 new tests, 714 total repo
tests (669 baseline + 45 new).

REMEDIATED 2026-08-30 (audit-019): the 018 IAU-1976 precession `_rot3`
was the TRANSPOSE of the standard form ([[c,s],[-s,c]] vs the eclipseTiming
convention [[c,-s],[s,c]]). The bug left a ~0.66 deg frame mismatch instead
of fixing the original 0.4 deg. Fixed in 018 (signed remediation commit).
Impact on RAAN rate: ~2.5e-3 deg/year (~3% of corrected formula magnitude),
well below the 9.78x short-period residual at i_sso.

Experiment 019 (Lunisolar Long-Period Terms and Secular-Limit Convergence)
is COMPLETE (2026-08-30) in
`research/orbital-mechanics/experiments/lunisolarLongPeriod/`: resolves the
018 ~10x residual as **mean-vs-osculating bias from finite-window linear
fit**, NOT as unmodelled Lunisolar physics. The 8-track investigation
(audit-019-synthesis-2026-08-30.md) identified:
- Annual solar forcing + lunar evection + variation bias the 1-year linear
  fit by 1-3Ã—10â»â´ deg/day (comparable to the corrected secular formula's
  +1.35e-4 deg/day). The 018 ~10x residual at i_sso is dominated by this
  bias, not by missing Lunisolar physics
- Window-length extrapolation Omega_dot_fit(W) = a + b/W + c/WÂ²
  extrapolates the secular limit; at h=600 km i_sso, gives Lunisolar
  ~+0.0036 deg/day (27x the corrected formula's +1.35e-4, confirming
  Track G's prediction of a 30x under-estimate at W â†’ âˆž)
- Cycle-averaged estimator (12 monthly segments) reduces bias to ~3% and
  at i=90Â° (cleanest J2-free test) gives 2.78x ratio vs corrected formula
  (matching 018's 2.81x), confirming the residual structure is the same
  at both inclinations
- FFT periodicity detects annual + harmonics (365 d, 182 d, 121 d, 91 d,
  73 d top-5 dominant periods), consistent with Track F's prediction
The corrected secular formula does NOT need an evection/variation
addition; the 1-year linear fit is just a biased estimator. Window-length
extrapolation to W â†’ âˆž is the canonical numerical bridge. 5 figures, 43
new tests, 757 total repo tests (714 baseline + 43 new for 019).

Experiment 020 (Lunisolar Long-Arc Secular-Limit Validation) is COMPLETE
(2026-08-30) in `research/orbital-mechanics/experiments/lunisolarSecularLimit/`:
8-track audit-020 (`localdocs/reports/audit-020-track-{1..8}-*.md`) +
1-yr arc at h=600 km i_sso with 4-phase ensemble. **The 019 extrapolation
+0.0036 deg/day (27x the corrected formula) is NOT validated as the
asymptotic secular limit**:
- Track 3 (estimator theory) shows the 019 polynomial-in-1/W extrapolation
  has NO theoretical asymptotic basis; the actual OLS bias scales as
  O(1/WÂ²) for fast harmonics and O(A_k Ï‰_k) constant for slow harmonics.
- Track 7 (hostile review) shows the 019 i=90Â° extrapolation sign-flips
  between linear (+1.7e-4) and quadratic (-3.7e-4) models â€” a smoking gun
  for mis-specification.
- Track 4 (implementation audit) discovers that at 2026 (near descending
  lunar node) the actual lunar i3 is ~18.29Â°, not the secular mean
  28.584Â°, making the apples-to-apples 2026 ratio closer to 13-14Ã—, not
  9.78Ã—.
- Synthetic oracle test: harmonic-regression estimator (Track 3
  recommendation (f)) recovers known secular to machine precision (7e-16
  deg/day bias) on synthetic data with 019 FFT amplitudes.
- Real data at h=600 km i_sso, 1-yr arc, 4-phase ensemble: direct OLS (a),
  secant (g), node-vector (n) all give Lunisolar ~+1e-3 to +2e-3 deg/day,
  ratio 9.3Ã— to corrected cf â€” reproducing 018 finding. Harmonic
  regression (f) is FRAGILE on the 1-yr ascending-node data (j2_only
  harmonic regression is stable but full_model swings from 8.89e-1 to
  1.11e+0 deg/day; unmodelled short-period content is aliased into a
  long-period harmonic the regression interprets as secular drift).
- 14 new tests, 771 total (757 baseline + 14 new).
- 5 figures (estimator comparison, ratio to corrected cf, phase dependence,
  harmonic amplitude recovery, estimator bias on synthetic oracle).
**Verdict**: the corrected 018 formula gives the correct SIGN but
UNDER-ESTIMATES the 1-yr numerical Lunisolar rate at i_sso by ~9.3Ã—. The
019 extrapolation is reported only as a diagnostic, NOT a robust
asymptotic measurement. The secular limit at W â†’ âˆž remains UNRESOLVED
at the 1-yr arc.

Mission `mission_lunisolar_closure` (POST_ROADMAP_PROBE Â§13.1, LAB_CONSTITUTION Â§13.1)
is COMPLETE (2026-09-03) in
`research/orbital-mechanics/missions/mission_lunisolar_closure/`: 18.6-yr direct arc
(one full lunar nodal cycle) at h=600 km with 3 inclinations (i_sso=97.79Â°,
i=90Â°, i=30Â°) and 4 independent estimators (direct OLS, secant, harmonic regression,
node-vector OLS); single phase (lunar anomalistic zero). Parallelized across 8 cores
via multiprocessing.Pool (7 workers, ~67 min wall-clock vs ~6 hr sequential).
Byte-pinned DE441 Sun + Moon snapshots, 2026-01-01 â†’ 2045-01-01, daily cadence,
ICRF/TDB. Mission results.json written by the campaign (canonical source;
verified bit-identical on re-run). **13 new tests, 784 total repo tests
(771 baseline + 13 new).**
**HEADLINE FINDING (refutes 018/020)**: at the 18.6-yr arc (one full lunar nodal
cycle), the numerical Lunisolar RAAN rate (full - j2_only) at h=600 km is:
- i_sso (97.79Â°): numerical -2.29e-2 deg/day (retrograde); corrected cf +1.35e-4
  deg/day (prograde); ratio -170Ã—; **SIGN DISAGREEMENT**
- i=90Â°: numerical +4.70e-3 deg/day (prograde); corrected cf +1.74e-4 deg/day
  (prograde); ratio +27Ã—; sign agrees but magnitude residual 27Ã—
- i=30Â°: numerical -3.47e-4 deg/day (retrograde); corrected cf +4.55e-5 deg/day
  (prograde); ratio -7.6Ã—; **SIGN DISAGREEMENT**
The 018/020 conclusion that "the corrected formula gives the correct SIGN" is
REFUTED at i_sso and i=30 at the 18.6-yr arc. The 019 polynomial-in-1/W
extrapolation to W â†’ âˆž is REFUTED (the 18.6-yr numerical is opposite sign
from the extrapolation). All four estimators agree within ~4% at each
inclination. Mission state: **PARTIALLY-VERIFIED-WITH-OPEN-QUESTION**;
recommends follow-on J2 Ã— Lunisolar coupling derivation as next mission.
See `localdocs/knowledge/lunisolar-closure-021.md` and
`localdocs/reports/mission-lunisolar-closure-2026-09-03.md` for full scientific
report. Mission repo path was moved from the prior session's `mission_lunisolar_closure/`
at the repo root to the constitutional `research/orbital-mechanics/missions/mission_lunisolar_closure/`
per `LAB_CONSTITUTION.md Â§2.3`; MANIFEST.json updated to drop the R:-path dependency
violating the durable-storage rule; `code_hashes()` rewritten to walk up to the
repo root robustly. Mission was the first clean scientific execution after the
adoption of `LAB_CONSTITUTION.md` (commit 9d9a495); it preserves the prior
scaffold (3 commits at the repo root) without rewriting history.
Next: per `POST_ROADMAP_PROBE.md Â§13.2-13.5` and the mission report, candidates are
(1) Repeat-Ground-Track Targeting (composition mission); (2) J2 Ã— Lunisolar
coupling derivation (follow-on to mission_lunisolar_closure); (3) Estimation
Doctrine Graduation (capability mission); (4) NRHO/Cislunar (frontier).
