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

## Current Priority

Orbital-mechanics flagship: experiments 001–009 are complete (numerics
foundation, Kepler validation, Kepler solvers, Hohmann transfer, bi-elliptic
vs Hohmann crossover, combined transfer + plane change, gravity assist,
ground tracks, J2 precession). Experiment 006 was adversarially audited +
closed (2026-08-17); a synthesis report for 001–006 is in `localdocs/reports/`.
Experiment 008 (ground tracks, spherical-Earth) is COMPLETE (2026-08-21) in
`research/orbital-mechanics/experiments/groundtracks/`; Experiment 009
(J2 precession, secular nodal/apsidal rates with numerical validation) is
COMPLETE (2026-08-22) in `research/orbital-mechanics/experiments/j2Precession/`.
Next: 010 orbit decay, reusing `src/lab_utils/`, the Exp 002 RK4 machinery,
and the Exp 004–009 closed-form machinery — do not rebuild scaffolding. See
`localdocs/roadmap.md`.
