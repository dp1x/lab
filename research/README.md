# research/ — Experiments by Domain

Every domain holds self-contained experiments under `experiments/`. Each experiment
has an experiment card (README.md), code, tests, and results. See
`localdocs/templates/experiment_template.md` for the format.

## Domains

| Domain | Status | Focus |
|--------|--------|-------|
| [numerics](numerics/) | active (foundation) | numerical methods, error analysis — feeds all others |
| [orbital-mechanics](orbital-mechanics/) | **FLAGSHIP** | orbits, transfers, perturbations, satellite trajectory |
| energy | planned (roadmap text only) | power flow, batteries, grid modelling |
| computer-architecture | planned (roadmap text only) | pipelines, caches, scheduling |
| cybersecurity | planned (roadmap text only) | cryptography, protocols, attack modelling |

Only domains with real content get a directory. Planned domains stay text in
`localdocs/roadmap.md` until an experiment exists. See the roadmap for the full
sequence.