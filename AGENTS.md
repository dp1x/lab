# AGENTS.md — Laboratory Operating Manual

This file governs how AI agents operate inside the Computational Research
Laboratory. It is the concise operating manual; the full charter lives at
`localdocs/charter.md` and the experiment sequence at `localdocs/roadmap.md`.

## Purpose

Convert temporary, cheap, unreliable AI compute into **verified, reproducible,
compounding knowledge** by chaining the agent to deterministic reality-checks.
The output is a career portfolio: a research notebook of validated experiments,
a structured Obsidian knowledge base, and engineering discipline.

## Repo Structure (operate exactly like this)

```
C:\Users\Dhane\lab\
├── AGENTS.md                this manual
├── README.md                lab overview + index
├── .gitignore               excludes keys, artifacts, private data
├── pyproject.toml / uv.lock   Python env + pinned deps (uv)
├── localdocs/                LAB'S OWN WRITING
│   ├── charter.md           source of truth (16-section master document)
│   ├── roadmap.md           experiment sequence (read this before any task)
│   ├── templates/           experiment card + note templates
│   ├── knowledge/           Obsidian-compatible knowledge base (compounds)
│   └── reports/             multi-experiment syntheses (every ~5 exps)
├── webdocs/                  EXTERNAL web material (public sources only)
├── research/                 one dir per domain
│   ├── <domain>/
│   │   └── experiments/<name>/   ONE self-contained experiment
│   │       ├── README.md         # experiment card (template)
│   │       ├── experiment.py     # runnable implementation
│   │       ├── tests/            # validation tests (pytest)
│   │       └── results/          # outputs, figures (JSON, committed)
│   └── README.md                 # domain index
├── src/lab_utils/            shared reusable utilities (metrics, I/O)
├── tools/                    lab-level scripts (scaffolder)
└── data/                     shared datasets (gitignored)
```

## The Loop — run every experiment the same way

```
Question → Theory → Design → Implement (experiment.py, deterministic)
  → Test (pytest, must pass BEFORE trusting results)
  → Run → Validate (analytic solutions, invariants, real data)
  → Document (card: results fill into README.md)
  → Knowledge (localdocs/knowledge/<topic>.md, Obsidian-linked)
  → Commit (one experiment per commit) → Next Question
```

## Operating Rules

1. **Deterministic only.** Simulations, math, physics, algorithms, optimization.
   Fixed seeds, no time-dependent nondeterminism. No unattributed speculation.
2. **Verify before trusting.** Tests must pass before results are recorded. Validate
   via analytic solutions, conservation laws, published data, or benchmarks. Reality
   is the verification layer: if wrong, reality reveals it.
3. **Reproduce everything.** Follow the experiment template exactly:
   question → theory → references → assumptions → methodology → implementation →
   validation → results → limitations → next question.
4. **Never fabricate** papers, results, benchmarks, or citations. Separate known
   facts / assumptions / hypotheses / results. Cite real references.
5. **No proprietary content.** Never commit API keys, tokens, passwords, personal
   info, or private data. All lab material is public-sourced or our own results.
6. **Resources are precious.** SSD/TBW, RAM, CPU. Avoid rebuilds, huge temp files,
   abandoned processes. Terminate what you start. Prefer incremental work.
7. **Decisions are justified.** Before new dependencies/frameworks/tools/dirs,
   prefer existing code/infrastructure. Favor simplicity; complexity must justify
   itself.
8. **Delegate with parallelism.** Use subagents for independent work (research,
   implementation, testing, review). Never edit the same file in parallel.
9. **Documentation is memory.** Every experiment and meaningful change writes its
   durable artifacts (card, results, knowledge note). Update cards when results change.

## Environment

- Python managed by **uv** — always `uv sync` first, run with `uv run python ...`.
- Run all tests: `uv run pytest`.
- Add dependencies with `uv add <pkg>`, and update affected docs.
- Results: small, human-readable JSON preferred; commit results + exact figures
  (committed) for reproducibility; never commit large binaries.

## New Experiment Checklist

1. `research/<domain>/experiments/<descriptive_name>/`
2. Copy `localdocs/templates/experiment_template.md` → `README.md`, fill card.
3. Write `experiment.py` (deterministic, fixed seeds).
4. Write tests in `tests/`; run `uv run pytest` until green.
5. Run experiment → write `results/results.json` + figures.
6. Record results + limitations in the card.
7. Write Obsidian note in `localdocs/knowledge/` (link to prior notes).
8. `git add` + commit (with a message describing what was verified).

## Current Priority — Orbital Mechanics Flagship

Follow `localdocs/roadmap.md`. Numerics foundation (001) is complete; the next work
is the orbital-mechanics sequence (002 Kepler orbit validation, then 003 Kepler's
equation solvers, 004 Hohmann transfer, etc.). Reuse `src/lab_utils/` and templates
— do not rebuild scaffolding.

## Sweep Methodology

One experiment = one research question, many parameter combinations. Big CSVs go
to `data/` (gitignored); commit a small summary JSON to `results/`. Domain dirs:
only `numerics/` and `orbital-mechanics/` exist — energy/computer-architecture/
cybersecurity are roadmap text only until real content exists.