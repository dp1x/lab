# Computational Research Laboratory

A long-term computational research environment that converts temporary abundant AI
inference compute into permanent technical value.

## Mission

The laboratory produces:

- reproducible scientific experiments
- validated computational models
- engineering simulations
- technical documentation
- reusable software infrastructure
- a growing scientific knowledge base
- open-source-quality engineering artifacts

This is a monorepo laboratory, operated as a small independent research organization.

## Core Principles

1. **Deterministic Engineering** — prioritize work where reality provides feedback
   (simulations, mathematics, physics, algorithms, optimization). *"If the result is
   wrong, reality should reveal it."*
2. **Reproducibility Above Everything** — every experiment contains research question,
   theory, references, assumptions, methodology, implementation, validation, results,
   limitations, and future improvements. Another person can reproduce it.
3. **Scientific Integrity** — never fabricate papers, results, benchmarks, or outcomes.
   Separate known facts / assumptions / hypotheses / results.
4. **Self-Improving Laboratory** — improve templates, tooling, and workflow only with
   measurable benefit; research output comes before framework rewrites.
5. **Documentation Is Memory** — every meaningful activity creates durable artifacts:
   experiment cards, reports, notes, results, Obsidian-compatible knowledge entries.
6. **Aggressive Delegation and Parallelism** — heavy subagent use for independent
   work; never concurrent edits to the same file.

## Repository Layout

```
├── AGENTS.md                operating manual for AI agents
├── README.md                this file
├── pyproject.toml / uv.lock   Python environment + pinned dependencies (uv)
├── localdocs/               THE LAB'S OWN WRITING
│   ├── charter.md          source of truth (16-section master document)
│   ├── roadmap.md          experiment sequence (orbital-mechanics flagship)
│   ├── templates/          experiment card + note templates
│   ├── knowledge/          Obsidian-compatible knowledge base (compounding)
│   └── reports/            multi-experiment syntheses
├── webdocs/                 EXTERNAL web material (public sources only)
├── research/
│   ├── <domain>/experiments/<name>/   self-contained experiment
│   │   ├── README.md                 experiment card
│   │   ├── experiment.py             runnable implementation
│   │   ├── tests/                    validation tests (pytest)
│   │   └── results/                  outputs, figures (JSON, committed)
│   └── <domain>/README.md           domain index
├── src/lab_utils/          shared reusable utilities (metrics, I/O, validation)
├── tools/                  lab-level scripts (scaffolder, ...)
└── data/                   shared datasets (gitignored)
```

## Quickstart

Requires [uv](https://docs.astral.sh/uv/) (Python 3.12 managed by uv):

```bash
uv sync            # create venv, install dependencies
uv run pytest      # run all tests
uv run python research/<domain>/experiments/<name>/experiment.py
```

## Research Workflow

```
Research → Question → Theory → Design → Implementation → Testing → Simulation
    → Validation → Documentation → Knowledge Update → Improvement → Next
```

Every experiment: card → code → tests green → run → validate → results →
knowledge note (Obsidian) → commit.

## Operating Reference

- [Laboratory charter (source of truth)](localdocs/charter.md)
- [Roadmap & experiment sequence](localdocs/roadmap.md)
- [Laboratory operating manual](AGENTS.md)
- [Experiment template](localdocs/templates/experiment_template.md)
- [Knowledge base](localdocs/knowledge/)
- [Research overview](research/)

## Research Domains

1. **Orbital Mechanics — flagship (planned next)** (orbits, transfers,
   perturbations; validation against NASA/JPL Horizons is a target, not yet achieved)
2. Numerics (foundation — 001 complete)
3. Energy Systems (planned — power flow, batteries, grid)
4. Computer Architecture (planned — pipelines, caches, scheduling)
5. Cybersecurity (planned — cryptography, protocols)