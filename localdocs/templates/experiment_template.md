# Experiment Card: <Title>

> Status: planning | running | complete
> Date: YYYY-MM-DD
> Domain: <domain>
> Experiment dir: `research/<domain>/<name>/`

## Research Question

_What question does this experiment answer? State it precisely._

## Background Theory

_Theory behind the problem. Equations, models, and the physics/math that ground the
experiment._

## References

- Author, Title, Year, source. (Real references only — no fabricated citations.)

## Assumptions

- _List all assumptions made. Mark each as verified / plausible / idealization._

## Methodology

_How the experiment is performed: setup, parameters, algorithms, sample sizes. Include
enough detail to reproduce._

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, scipy, matplotlib
- Runtime: `uv run python experiment.py`
- Determinism: fixed seeds / no RNG (state how)

## Validation Method

_How are results checked? Analytic solutions, conservation laws, cross-checks,
convergence tests, unit tests (list test file)._

## Results

_Summary tables, key numbers, figures (stored in `results/`)._

## Limitations

_What the results do NOT cover; numerical and modeling caveats._

## Future Improvements

_Next steps and open extensions._

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command to reproduce: `uv sync && uv run pytest && uv run python experiment.py`
