"""Regression guard: every byte-pinned external snapshot must be `-text`-protected.

Why this test exists (2026-10-03 provenance audit)
---------------------------------------------------
The repository sets ``core.autocrlf=true`` at *system* level on this host. Git therefore
converts LF -> CRLF on checkout unless a path is marked ``-text`` (never convert) or given
``eol=lf``.

``mission_lunisolar_closure/.gitattributes`` attempted to protect the 2.7 MB of byte-pinned
DE441 Sun/Moon vectors, but **a nested ``.gitattributes`` resolves its patterns relative
to its own directory**. Its patterns were written repo-root-relative, so they resolved to
``mission_lunisolar_closure/research/orbital-mechanics/...`` -- a path that does not exist
-- and all four rules were silently inert. ``git check-attr text`` reported ``unspecified``
for every DE441 file. A contributor cloning on Linux would have received LF files whose
SHA-256 did not match the pinned ``f2c4f048...`` / ``aee85099...`` values, with no
explanation and no failing test.

The failure mode is invisible by construction: an attribute rule that matches nothing is
indistinguishable from an absent rule. This test makes it visible.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]

#: Byte-pinned external snapshots whose identity is a recorded SHA-256.
PINNED = [
    "research/orbital-mechanics/experiments/jplValidation/reference",
    "research/orbital-mechanics/experiments/eclipseTiming/reference",
    "research/orbital-mechanics/experiments/lunisolarVerification/reference",
    "research/orbital-mechanics/missions/mission_lunisolar_closure/reference",
]


def _check_attr(path: str) -> str:
    out = subprocess.run(
        ["git", "check-attr", "text", "--", path],
        cwd=ROOT, capture_output=True, text=True, timeout=60,
    )
    if out.returncode != 0:
        pytest.skip(f"git check-attr unavailable: {out.stderr.strip()}")
    line = out.stdout.strip()
    # Format: "<path>: text: <value>"
    return line.rsplit(":", 1)[-1].strip()


@pytest.mark.parametrize("rel_dir", PINNED)
def test_pinned_snapshot_directories_exist(rel_dir: str) -> None:
    assert (ROOT / rel_dir).is_dir(), f"missing pinned-reference directory: {rel_dir}"


def test_every_pinned_snapshot_is_marked_no_text() -> None:
    """Every ``.txt``/MANIFEST under a pinned dir must resolve to ``text: unset``."""
    offenders = []
    checked = 0
    for rel_dir in PINNED:
        d = ROOT / rel_dir
        for p in sorted(d.rglob("*")):
            if not p.is_file() or p.suffix not in (".txt", ".json"):
                continue
            checked += 1
            val = _check_attr(p.relative_to(ROOT).as_posix())
            if val != "unset":
                offenders.append((p.relative_to(ROOT).as_posix(), val))
    assert checked > 0, "no pinned snapshot files were found to check"
    assert not offenders, (
        "byte-pinned snapshots are NOT protected from line-ending conversion "
        f"(text: expected 'unset', got other). core.autocrlf=true would corrupt them "
        f"on clone. Offenders: {offenders}"
    )


def test_source_files_are_pinned_to_lf() -> None:
    """``*.py`` must be checked out LF so worktree bytes match git blobs.

    Result artifacts fingerprint source files with ``sha256(path.read_bytes())``. With
    ``core.autocrlf=true`` and no ``eol`` attribute those worktree bytes differ from the
    blob, so a recorded fingerprint cannot be verified from a fresh clone.

    ``check-attr text`` reports ``set`` for a file that carries an explicit ``text``
    attribute, which is exactly what ``*.py text eol=lf`` produces. The test asserts the
    attribute is PRESENT, not that it equals any particular value.
    """
    unset = []
    for rel in ("src/lab_utils/estimation.py", "src/lab_utils/ephemeris.py",
                "research/orbital-mechanics/missions/mission_mean_element_leakage/"
                "phase0_baseline_breathing.py"):
        p = ROOT / rel
        assert p.is_file(), rel
        if _check_attr(rel) == "unspecified":
            unset.append(rel)
    assert not unset, (
        f"source files carry no text/eol attribute; with core.autocrlf=true their "
        f"worktree bytes will not match the git blob: {unset}"
    )


def test_gitattributes_patterns_are_root_relative_or_dir_relative() -> None:
    """Nested .gitattributes must not carry repo-root-relative patterns.

    Such a pattern resolves inside the nested file's own directory and silently matches
    nothing. Any ``*.txt``-style pattern here is a bug; ``reference/...`` is correct.
    """
    nested = (ROOT / "research/orbital-mechanics/missions/mission_lunisolar_closure"
              / ".gitattributes")
    if not nested.is_file():
        pytest.skip("no nested .gitattributes present")
    bad = []
    for ln in nested.read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        pattern = s.split()[0]
        if pattern.startswith("research/") or pattern.startswith("src/"):
            bad.append(pattern)
    assert not bad, (
        f"nested .gitattributes patterns are resolved relative to its own directory; "
        f"repo-root-relative patterns match nothing: {bad}"
    )