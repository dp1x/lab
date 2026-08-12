"""Scaffolder: create a new self-contained experiment from the templates.

Usage:
    uv run python tools/new_experiment.py <domain> <experiment-name>

Creates:
    research/<domain>/experiments/<name>/
    ├── README.md            # experiment card from template
    ├── experiment.py        # runnable skeleton
    ├── tests/test_<name>.py # validation test skeleton
    └── results/             # outputs, figures

Test module basenames must be unique per experiment, tests/ must NOT have an
__init__.py, and tests load experiment.py via importlib from its explicit
path: two experiments both named "experiment.py" (and two tests/test_*.py
files) otherwise collide in pytest's module registry and sys.modules.

Domain must be an existing dir under research/ or a new one (auto-created).
Update research/README.md when adding a domain (docs follow structure).
"""

from __future__ import annotations

import re
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research"
TEMPLATE_DIR = ROOT / "localdocs" / "templates"
CARD_TEMPLATE = TEMPLATE_DIR / "experiment_template.md"

SKELETON_PY = '''"""<Title> — experiment implementation.

Laboratory loop: implement → test (pytest green) → run → validate → document.
Determinism rule: fixed seeds / no RNG.
"""

from __future__ import annotations

import numpy as np


def answer_skeleton(a: float, b: float) -> float:
    """Placeholder: replace with the real computation."""
    return a + b


def main() -> None:
    print("Experiment scaffold ready. Implement the question, then fill the card.")


if __name__ == "__main__":
    main()
'''

TEST_SKELETON = '''"""Validation tests for {title}.

These must pass BEFORE any results are trusted (laboratory rule: verify before
trust). Replace the placeholder with the experiment's real invariants.

The experiment module is loaded via importlib from its explicit path (see
tools/new_experiment.py) so that multiple experiments with an "experiment.py"
module never collide in pytest/sys.modules.
"""

import importlib.util
from pathlib import Path

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "{module_name}", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)


def test_skeleton():
    assert experiment.answer_skeleton(2, 3) == 5
'''


def slugify(name: str) -> str:
    return re.sub(r"[^a-z0-9-]", "-", name.strip().lower())


def fill_card(template: str, domain: str, name: str) -> str:
    today = date.today().isoformat()
    return (
        template.replace("<Title>", name.replace("-", " ").title())
        .replace("<domain>", domain)
        .replace("YYYY-MM-DD", today)
    )


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    domain, name = argv
    domain = slugify(domain)
    name = slugify(name)
    if not name or not domain:
        print(f"error: invalid domain or name: {argv}", file=sys.stderr)
        return 2

    domain_dir = RESEARCH / domain
    domain_dir.mkdir(parents=True, exist_ok=True)
    if not (domain_dir / "README.md").exists() and not (domain_dir / "experiments").exists():
        (domain_dir / "README.md").write_text(
            f"# {domain.title()} — Experiments\n\nNew domain. See `localdocs/roadmap.md`.\n",
            encoding="utf-8",
        )
    exp_dir = domain_dir / "experiments" / name
    if exp_dir.exists():
        print(f"error: {exp_dir} already exists", file=sys.stderr)
        return 1

    exp_dir.mkdir(parents=True)
    tests_dir = exp_dir / "tests"
    tests_dir.mkdir()
    (exp_dir / "results").mkdir()

    card = CARD_TEMPLATE.read_text(encoding="utf-8")
    (exp_dir / "README.md").write_text(fill_card(card, domain, name), encoding="utf-8")
    (exp_dir / "experiment.py").write_text(SKELETON_PY, encoding="utf-8")
    module_name = f"{name.replace('-', '_')}_experiment"
    (tests_dir / f"test_{name.replace('-', '_')}.py").write_text(
        TEST_SKELETON.format(
            title=name.replace("-", " ").title(), module_name=module_name
        ),
        encoding="utf-8",
    )

    print(f"created {exp_dir}")
    print("next: fill README card, implement experiment.py, write real tests, run")
    print("  uv run pytest")
    print("  uv run python experiment.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))