# AGENTS.md — Research Lab Operating Manual

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
├── AGENTS.md              this manual
├── README.md              lab overview + index
├── LICENSE                MIT
├── .gitignore
├── pyproject.toml / uv.lock   Python env + pinned deps (uv)
├── localdocs/             LAB'S OWN WRITING (charter, roadmap, templates, knowledge)
├── webdocs/               external web material (public sources only)
├── research/<domain>/experiments/<name>/   one self-contained experiment
│   ├── README.md          experiment card
│   ├── experiment.py      runnable implementation
│   ├── tests/             validation tests (pytest)
│   └── results/           outputs, figures (JSON, committed)
├── src/lab_utils/         shared reusable utilities (metrics, I/O, validation)
├── tools/                 lab-level scripts (scaffolder)
└── data/                  shared datasets (gitignored; prefer scratch for large data)
```

## The Loop

```
Question → Theory → Design → Implement (experiment.py, deterministic)
  → Test (pytest, green BEFORE trusting) → Run → Validate (analytic solutions,
    invariants, published data) → Document (card + figures + results.json)
  → Knowledge (localdocs/knowledge/<topic>.md) → Commit → Next Question
```

## Operating Rules

1. **Deterministic only.** Fixed seeds/params; no time-dependent nondeterminism; no
   unattributed speculation.
2. **Verify before trusting.** Tests must pass before results are recorded. Validate
   via analytic solutions, conservation laws, published data, or benchmarks —
   "if wrong, reality reveals it."
3. **Reproduce everything.** Follow the experiment template exactly: question →
   theory → references → assumptions → methodology → implementation → validation →
   results → limitations → next question.
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
   parallel; the lead agent integrates. Avoid Goodhart — no agents spawned merely to
   satisfy a delegation metric.
9. **Documentation is memory.** Every experiment and meaningful change writes durable
   artifacts (card, results, knowledge note). Update cards when results change.
10. **Goodhart's Law.** Optimize research value, not activity. Do not create experiments
    merely to raise the count, spawn agents merely to delegate, generate datasets
    merely because storage exists, or use compute merely because it is free.

## Responsible Web Access

The lab may fetch public research material and public datasets. Preferred order:
**official API → official downloadable dataset → normal public webpage → browser
automation when necessary.** Browser automation is a research instrument, not a
scraping weapon.

Hard boundary — **never**: bypass CAPTCHAs, defeat access controls, circumvent rate
limits, rotate identities to evade blocks, spoof fingerprints for evasion, hammer
endpoints, parallel-request-burst public sites, or repeatedly retry a blocked
endpoint. If a site presents an antibot/CAPTCHA/login/verification gate, **STOP and
report it to the human** — do not attempt to defeat it. Use human-use-rate browsing,
cache acquired material, and respect published rate limits.

Search: prefer DuckDuckGo or a reputable SearXNG/Startpage instance for general
discovery; do not use Google/Bing for ordinary search. (Accessing an actual Google-
hosted public document/dataset is fine when it is the source.)

## Resource Architecture

- **`C:`** — permanent state (source, docs, tests, committed results). Keep lean.
- **`R:`** — disposable local scratch (venvs, caches, temp downloads, large sweeps).
  Check free space at runtime before large operations; never hard-code capacity.
- **Durable storage.** The repository root is the sole durable record: all
  reproducible source, knowledge, tests, contracts, provenance, results, figures,
  acquisition recipes, checksums, and valuable recovery state live under it.
  `R:` is never a runtime dependency of reproducibility and never the only
  location of anything worthwhile; if autonomous work is interrupted or blocked,
  recovery/handoff state is written under the repository (e.g. `AUTONOMOUS_HANDOFF_<exp>.md`).
- **Colab** — optional, ephemeral remote compute only for workloads that would
  otherwise cause excessive local SSD/CPU pressure. Hardware-as-available (no assumed
  GPU/CPU); never evade Colab limits; checkpoint long runs so they remain resumable
  from the repository. R: → Colab → results → R: → download → destroy runtime.

## Environment

- Python managed by **uv** — `uv sync` first, run with `uv run python ...`.
- All tests: `uv run pytest`.
- Add deps with `uv add <pkg>`; update affected docs.
- Sweep data: prefer R: scratch for large raw outputs; commit only compact summary
  JSON + figures to `results/`.

## New Experiment Checklist

1. `research/<domain>/experiments/<descriptive_name>/`
2. Copy `localdocs/templates/experiment_template.md` → README.md, fill the card.
3. Write `experiment.py` (deterministic).
4. Write tests; `uv run pytest` green.
5. Run → write `results/results.json` + figures.
6. Record results + limitations in the card.
7. Write an Obsidian note in `localdocs/knowledge/` (link to prior notes).
8. `git add` + commit (message describes what was verified).

## Remote-State Safety

- **Before any automated push:** verify the live remote tip immediately beforehand
  (`git ls-remote origin <branch>` or equivalent). Never assume the local
  `origin/<branch>` ref is current — it can be stale.
- Compare the live remote tip against the intended push target. If the remote changed
  unexpectedly, **STOP and reconcile**; never overwrite another agent's or the owner's
  work blindly.
- For history-rewriting operations, use `--force-with-lease=<expected-remote-tip>` only;
  never plain `--force`. After pushing, re-verify: remote tip == local HEAD, clean tree,
  and (where applicable) GitHub signature/verification state.
- Commit only from a verified-clean tree; never commit unrelated files, and keep the
  canonical Git identity (no repo-local identity overrides).

## Current Priority

Orbital-mechanics flagship: experiments 001–010 are complete (numerics
foundation, Kepler validation, Kepler solvers, Hohmann transfer, bi-elliptic
vs Hohmann crossover, combined transfer + plane change, gravity assist,
ground tracks, J2 precession, orbit decay). Experiment 006 was adversarially
audited + closed (2026-08-17); a synthesis report for 001–006 is in
`localdocs/reports/`. Experiment 008 (ground tracks, spherical-Earth) is
COMPLETE (2026-08-21) in `research/orbital-mechanics/experiments/groundtracks/`;
Experiment 009 (J2 precession, secular nodal/apsidal rates with numerical
validation) is COMPLETE (2026-08-22) in
`research/orbital-mechanics/experiments/j2Precession/`; Experiment 010 (orbit
decay / atmospheric drag — first non-conservative force: dissipation accounting
+ monotonicity doctrine, erfi/quadrature oracles, structural scalings,
co-rotation twins, J2 mean-element transient, reentry timing) is COMPLETE
(2026-08-22) in `research/orbital-mechanics/experiments/orbitDecay/`.
Next: 011 Lagrange points per `localdocs/roadmap.md`, reusing `src/lab_utils/`
and the Exp 002/006/008/009 propagator + element machinery — do not rebuild
scaffolding; graduating shared machinery to `src/lab_utils/orbits.py` is now
justified if it stays non-blocking.
Experiment 011 (Lagrange points / CR3BP — first rotating-frame experiment:
equilibria, Jacobi integral, Routh stability, nonlinear perturbation signatures,
dimensional cross-check, adversarial mutant battery) is COMPLETE (2026-08-22) in
`research/orbital-mechanics/experiments/lagrangePoints/`; shared machinery is now
graduated: `src/lab_utils/integrators.py` (generic rk4_step/rk4_propagate) and
`src/lab_utils/orbits.py` (element/Kepler canon, equivalence-pinned vs donors) —
use these instead of cloning per-experiment copies; CR3BP-specific code remains
experiment-local until a second consumer appears.
Next: 012 orbit classes per `localdocs/roadmap.md`, reusing
`src/lab_utils/orbits.py` + `src/lab_utils/integrators.py`.
Experiment 012 (orbit classes — constraint-defined families: SSO inclination lock +
finite existence boundary a_max = 12352.505 km, Molniya apsidal freeze + semi-synchronous
resonance + dwell geometry, GEO 1:1 fixed point with nonzero-rate negative control,
GTO budgets anchored to Exp 004, adversarial convention battery; finding: small-divisor
short-period dynamics near the critical inclination give a measured +323 s/orbit
Kepler-period excess) is COMPLETE (2026-08-23) in
`research/orbital-mechanics/experiments/orbitClasses/`; `j2_rhs` graduated into
`src/lab_utils/orbits.py` (second consumer after Exp 009, equivalence-pinned).
Next: 013 JPL ephemeris validation per `localdocs/roadmap.md`, reusing
`src/lab_utils/orbits.py` + `src/lab_utils/integrators.py`.
Experiment 013 (JPL ephemeris validation — ISS (-125544) vs pinned NASA/JPL Horizons
ICRF/TDB geometric states over a 3-day window; byte-pinned snapshot under the repo
(`-text` gitattributes), offline-deterministic analysis, exact-grid alignment,
reference-built RTN residuals, pre-registered decision rules) is COMPLETE (2026-08-24)
in `research/orbital-mechanics/experiments/jplValidation/`: J2 removes 99.33% of
residual RMS (bootstrap CI excludes zero); the declared drag tier WORSENS agreement at
primary beta=100 and the pre-declared beta band crosses zero only at its edge
(beta=400 -> 3.13 km vs M2 8.22 km) — documented verbatim without tuning; error budget
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
(36.7 vs 11.9/day main, equinoxes are the most eclipse-favorable for h=600); LST at the
ascending node drifts through 24 h/year at the sidereal-solar differential (4 min/day) —
the "LST-constant" intuition from the host research track was wrong; 6 figures; 34 new
tests, 581 total (547 baseline + 34 new). Shared machinery graduated to `src/lab_utils`:
`sso_inclination_rad` (3rd consumer after Exp 012 + Exp 014-implicit + Exp 015) in
`lab_utils/orbits`; `gmst_rad_iau1982`, `sun_unit_and_dist_km`, `subsolar_lon_rad`,
`eci_to_ecef`, `ecef_to_latlon`, `spherical_trig_latlon`, `lst_at_node_hours`,
`node_lon_from_raan_gmst` (2nd consumer after Exp 014 for the Sun/GMST, 2nd consumer
after Exp 008 for the ECI-to-lat/lon layer) in the new `lab_utils/earth_frames` module.
Next: 016 eclipse-aware station-keeping for dawn-dusk SSOs per `localdocs/roadmap.md`,
reusing the Exp 015 feasible-set table; declared follow-up candidates include a refined
J2 mean-vs-osculating coupling study and a multi-year LST-drift compensation budget.
