# Computational Research Laboratory — Charter

Master reference document for the laboratory. This document defines the mission,
philosophy, standards, structure, and operating rules. AGENTS.md is the concise
operating manual derived from this document; `charter.md` is the source of truth.

---

## 1. Mission

Convert temporary abundant AI inference capability into permanent technical value.

The laboratory is a personal research environment producing:

- reproducible scientific experiments
- validated computational models
- engineering simulations
- technical documentation
- reusable software infrastructure
- a growing scientific knowledge base
- open-source-quality engineering artifacts

The goal is not code generation. The goal is a continuously improving research
system:

```
Research → Experiment → Validation → Knowledge → Improvement → New Research
```

Assets must remain valuable even if AI tools disappear tomorrow.

## 2. Core Philosophy

Reality is the verification layer. The laboratory focuses on deterministic domains
where results can be tested:

- mathematics
- computational physics
- engineering simulation
- optimization
- algorithms
- software systems
- hardware analysis
- cybersecurity research

Avoid unsupported speculation. The standard:

> "If the result is wrong, reality should reveal it."

## 3. Scientific Integrity

Never fabricate:

- papers
- citations
- benchmarks
- experimental results
- equations
- measurements

Separate all claims into:

- **Known** — established facts (textbooks, papers, official docs, verified sources).
- **Assumptions** — simplifications or modelling choices.
- **Hypotheses** — ideas being tested.
- **Results** — actual experiment outputs.

## 4. Reproducibility Standard

Every experiment must be reproducible and must include:

1. **Research Question** — what is being investigated.
2. **Background Theory** — equations, principles, assumptions, expected behavior.
3. **References** — real academic papers, textbooks, official documentation.
4. **Methodology** — algorithms, parameters, implementation choices.
5. **Implementation** — source code, dependencies, configuration, run instructions.
6. **Validation** — analytic solutions, known results, benchmarks, conservation laws.
7. **Results** — numerical outputs, plots, tables, logs.
8. **Limitations** — missing physics, assumptions, numerical caveats.
9. **Next Research Question** — every experiment seeds the next one.

See `localdocs/templates/experiment_template.md` for the canonical card format.

## 5. Laboratory Operating Loop

```
Research
    ↓
Question
    ↓
Theory
    ↓
Design
    ↓
Implementation
    ↓
Testing
    ↓
Simulation
    ↓
Validation
    ↓
Documentation
    ↓
Knowledge Base Update
    ↓
Framework Improvement
    ↓
Next Experiment
```

Knowledge compounds over time. Each completed experiment produces a knowledge
note that later experiments build on.

## 6. Repository Structure

Monorepo — all domains share scientific tools, documentation standards, templates,
the Python environment, validation methods, and workflow.

```
Computational-Research-Lab/
├── AGENTS.md                # operating manual for AI agents (concise)
├── README.md                # lab overview + index
├── .gitignore                # excludes artifacts, keys, private data
├── pyproject.toml / uv.lock   # Python env + pinned dependencies (uv)
├── localdocs/                # LAB'S OWN WRITING (charter, roadmap, templates, notes)
│   ├── charter.md           # this document (source of truth)
│   ├── roadmap.md           # domain roadmap and experiment sequence
│   ├── templates/           # experiment card + note templates
│   ├── knowledge/           # Obsidian-compatible knowledge base (compounding)
│   └── reports/             # syntheses spanning multiple experiments
├── webdocs/                 # EXTERNAL material from the web (public sources)
│   └── README.md           # policy + index
├── research/                 # one dir per domain, experiments self-contained
│   ├── README.md           # research overview + index
│   ├── numerics/experiments/         # foundation domain (verified numerical methods)
│   ├── orbital-mechanics/experiments/   # FLAGSHIP domain
│   └── <domain>/experiments/        # created only when it has real content
├── src/lab_utils/         # shared reusable utilities (metrics, I/O, validation)
├── tools/                 # lab-level scripts (scaffolder, etc.)
└── data/                  # shared datasets (gitignored)
```

Each experiment directory is self-contained:

```
research/<domain>/experiments/<name>/
├── README.md       # experiment card
├── experiment.py   # runnable implementation
├── tests/          # validation tests (pytest)
└── results/        # numeric outputs, figures (JSON preferred, committed)
```

## 7. Technology Stack

Primary language: Python. Scientific ecosystem: NumPy, SciPy, Matplotlib, pandas,
SymPy, Jupyter, pytest.

Later additions by phase:

- Aerospace engineering: Astropy, poliastro, SPICE tools.
- Energy systems: PyPSA, pandapower, oemof.
- Visualization: Plotly, PyVista.

Add dependencies with `uv add <pkg>` and update relevant docs. Avoid unnecessary
dependencies.

## 8. Shared Infrastructure

Maintain reusable utilities under `src/lab_utils/`:

- **Metrics** — error measurement, convergence analysis, comparison, statistics.
- **Results** — experiment metadata, parameters, outputs, version info (JSON).
- **Templates** — every experiment uses the standard card + note templates.

Favor existing infrastructure over new ad-hoc solutions. Improve utilities only
when there is measurable benefit.

## 9. Knowledge Base

Documentation is memory. Every experiment creates:

- Markdown report (experiment card)
- Obsidian-compatible note in `localdocs/knowledge/`
- results, logs, references

Knowledge compounds:

```
Experiment 001 → knowledge note → Experiment 002 builds on it → improved methodology
```

## 10. Research Roadmap

### Phase 1 — Numerics (foundation)

- **001 Numerical Integrator Study** — COMPLETE. Convergence order and energy
  preservation of integrators (Euler, RK2, RK4, symplectic Euler, velocity Verlet).
- Foundation name `physics/` → `numerics/`: the substrate is verified numerical
  methods, not "physics" broadly.
- Duffing oscillator (nonlinear dynamics, chaos, Poincaré sections) deferred —
  orbital mechanics launched first. Continue foundation work only where it
  directly serves the flagship.

### Phase 2 — Orbital Mechanics (flagship domain, renamed from aerospace)

How objects move under gravity in space. Compact physics, closed-form checks,
real-world data (NASA/JPL Horizons) as the answer key. Sequence:

- Kepler orbit validation (002)
- Kepler's equation solvers (003)
- Hohmann transfer (004)
- Bi-elliptic vs Hohmann (005)
- Plane-change maneuvers (006)
- Gravity assist / slingshot (007)
- Ground tracks (008)
- J2 precession (009)
- Orbit decay (010)
- Lagrange points (011)
- Orbit classes (012)
- JPL ephemeris validation (013)
- Eclipse timing, launch windows, trajectory optimization (014+)

Sweep methodology: one experiment, one question, many parameter combinations;
CSVs → `data/` (gitignored), summary JSON → `results/` (committed).

### Phase 3 — Energy Systems (second pillar)

Renewable generation, battery systems, grid optimization, energy modelling:

- solar forecasting, battery degradation modelling, grid simulations

### Phase 4 — Computer Architecture systems engineering simulation

- CPU pipeline simulator, cache simulator, scheduling algorithms.

### Phase 5 — Cybersecurity — protocol simulation & analysis

Cryptographic analysis, secure protocol modelling, vulnerability testing frameworks.

## 11. Delegation and Parallelism

Subagents are encouraged for independent work: literature research, implementation,
testing, verification, documentation, review. Run independent work in parallel.
Never let multiple agents edit the same file simultaneously. The lead agent integrates
outputs. Parallelism must increase quality, not create confusion.

## 12. Operating Rules for AI Agents

Agents behave as a research organization — determine the question, design the
experiment, implement, validate, document, improve. See AGENTS.md for the concise
operating manual.

## 13. Self-Improvement

Improve templates, utilities, testing, workflow, and documentation through measured
gains. Do not endlessly redesign infrastructure. Infrastructure exists to enable
research — research output comes first.

## 14. Resource Efficiency

Constraints: SSD lifetime (TBW), storage, RAM, CPU. Avoid: unnecessary
compilation, giant temporary files, repeated rebuilds, duplicated datasets,
abandoned processes. Clean temporary files; stop unnecessary processes; free
resources. Prefer incremental changes, lightweight tools, efficient workflows.

## 15. Decision Making

Before adding dependencies/frameworks/tools/architecture ask: Is it necessary?
Does an existing solution exist? Does it improve research output? Prefer
simplicity. Complexity most yet justify itself.

## 16. Definition of Success

The laboratory succeeds when it produces over time:

- **Technical assets** — validated simulations, software tools, datasets, docs.
- **Knowledge assets** — scientific notes, experiment history, lessons learned.
- **Portfolio assets** — reproducible engineering projects, open-source, proven skill.
- **Personal growth** — deeper understanding of mathematics, physics, engineering, CS.

---

## Final Principle

The laboratory is not an AI code generator. It is a computational research system.

AI provides acceleration. The permanent value comes from validated experiments,
documented knowledge, engineering discipline, accumulated understanding, and
reproducible work.

> Build continuously. Validate everything. Document everything. Improve the system
> over time. Prioritize science and engineering that reality can verify.