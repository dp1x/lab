"""Results I/O: deterministic, human-readable JSON result storage.

Result files carry metadata (package versions, timestamp, git commit if available)
so every recorded result is traceable to the code that produced it.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from typing import Any

import numpy as np

__all__ = ["save_json_result"]


def _git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=os.getcwd(),
        )
        return out.stdout.strip() if out.returncode == 0 else "not-a-git-repo"
    except Exception:
        return "unknown"


def _numpy_to_python(obj: Any) -> Any:
    """Recursively convert numpy scalars/arrays to JSON-serializable Python types."""
    if isinstance(obj, dict):
        return {str(k): _numpy_to_python(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_numpy_to_python(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float):
        return round(obj, 12)  # avoid float repr noise in saved JSON
    return obj


def save_json_result(
    path: str, data: dict[str, Any], *, name: str, description: str = ""
) -> str:
    """Save experiment results to ``path`` as JSON with traceability metadata.

    Returns the absolute path written.
    """
    payload: dict[str, Any] = {
        "meta": {
            "name": name,
            "description": description,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "git_commit": _git_commit(),
            "python_platform": platform.platform(),
            "python_version": platform.python_version(),
        },
        "results": _numpy_to_python(data),
    }
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    return path
