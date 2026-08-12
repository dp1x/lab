# Research Lab

A small, reproducible computational research environment for structured investigations
across mathematically and computationally tractable domains.

## Mission

The laboratory conducts computational investigations, validates them against theory
or external reality, and accumulates the results, code, and documentation into a
coherent, reproducible research record.

The laboratory produces:

- reproducible computational experiments with traceable results
- validated computational models
- reusable software infrastructure
- technical documentation and an Obsidian-linked knowledge base
- engineering-quality, open-source-style artifacts

This is a monorepo laboratory, operated as a small independent research organization.

## Core Principles

1. **Deterministic Engineering.** Prioritize work where the computer can judge whether
   the result is right — simulations, mathematics, algorithms, optimization.
   *"If the result is wrong, reality should reveal it."*
2. **Validation Before Trust.** Tests must pass and results must agree with analytic
   solutions, conservation laws, published data, or benchmarks before they are recorded.
3. **Reproducibility.** Every experiment contains: question → theory → references →
   assumptions → methodology → implementation → validation → results → limitations →
   next question. The full run is `uv sync && uv run pytest && uv run python experiment.py`.
4. **Scientific Integrity.** Never fabricate papers, results, benchmarks, or citations.
   Separate known facts, assumptions, hypotheses, and results. Cite real references.
5. **Self-Improving Laboratory.** Improve templates, tooling, and workflow only with
   measurable benefit; research output comes before framework rewrites.
6. **Documentation Is Memory.** Every meaningful activity creates durable artifacts:
   experiment cards, figures, results, and Obsidian-compatible knowledge notes.
7. **Purposeful Delegation.** Subagents may be used when doing so materially improves
   research quality, throughput, or verification. Delegation is a means, not a quota.
8. **Goodhart's Law.** Do not optimize for activity metrics instead of research value —
   no experiments merely to increase count, no agents merely to fill quotas, no
   compute/datasets generated merely because resources are available.

## Repository Layout

```
├── AGENTS.md              operating manual for research agents
├── README.md              this file
├── LICENSE                MIT
├── pyproject.toml / uv.lock  Python env + pinned dependencies (uv)
├── localdocs/             THE LAB'S OWN WRITING
│   ├── charter.md         source of truth (mission, philosophy, standards)
│   ├── roadmap.md         experiment sequence (orbital-mechanics flagship)
│   ├── templates/         experiment card + knowledge note templates
│   ├── knowledge/         Obsidian-compatible knowledge base (compounding)
│   └── reports/           multi-experiment syntheses (every ~5 experiments)
├── webdocs/               EXTERNAL web material fetched for reference (public only)
├── research/              one dir per domain, experiments self-contained
│   ├── <domain>/experiments/<name>/
│   │   ├── README.md      experiment card
│   │   ├── experiment.py  runnable implementation
│   │   ├── tests/         validation tests (pytest)
│   │   └── results/       outputs, figures (JSON preferred, committed)
│   └── <domain>/README.md  domain index
├── src/lab_utils/         shared reusable utilities (metrics, I/O, validation)
├── tools/                 lab-level scripts (scaffolder, ...)
└── data/                  shared datasets (gitignored; large ephemeral data preferred on scratch)
```

## Quickstart

Requires [uv](https://docs.astral.sh/uv/) (Python 3.12):

```bash
uv sync          # create venv, install dependencies
uv run pytest    # run the full test suite
uv run python research/<domain>/experiments/<name>/experiment.py
```

## Research Workflow

```
Question → Theory → Design → Implement → Test (pytest) → Run →
Validate → Document (card + figures + results) → Knowledge note → Commit → Next
```

Every experiment: card → deterministic code → tests green → run → validate against
theory/real data → results → knowledge note → commit.

## Research Domains

1. **Orbital Mechanics — flagship (next)** — orbits, transfers, perturbations; validated
   against Kepler's laws and, where feasible, NASA/JPL Horizons ephemeris.
2. **Numerics — foundation (001 complete)** — verified numerical methods (integrators,
   error analysis) that applied domains build on.
3. Energy Systems (planned) — power flow, batteries, grid modelling.
4. Computer Architecture (planned) — pipelines, caches, scheduling.
5. Cybersecurity (planned) — cryptography, protocols, attack modelling.

## Resource Architecture

- **`$REPO_ROOT`** — permanent laboratory state: source, docs, tests, committed
  results. Kept lean; avoid unnecessary writes/rebuilds.
- **`R:`** — disposable local scratch (virtualenvs, caches, downloaded temp files, large
  sweeps). Check free space at runtime before large operations; do not hard-code capacity.
- **Google Colab** — optional, *ephemeral* remote compute used only when a workload would
  otherwise cause excessive local SSD/CPU pressure. Hardware is dynamically allocated (no
  assumed GPU/CPU); never attempt to circumvent Colab limits.

## Operating Reference

- [Laboratory charter (source of truth)](localdocs/charter.md)
- [Roadmap & experiment sequence](localdocs/roadmap.md)
- [Laboratory operating manual](AGENTS.md)
- [Experiment template](localdocs/templates/experiment_template.md)
- [Knowledge base](localdocs/knowledge/)
- [Research overview](research/)
