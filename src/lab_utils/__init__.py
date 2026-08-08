"""lab_utils — shared utilities for the Computational Research Laboratory.

Contains only reusable code: metrics, numerical helpers, and result I/O.
Each module is self-contained and covered by tests in `src/lab_utils/tests/`.
"""

from lab_utils.metrics import (
    convergence_rate,
    l2_norm_error,
    max_abs_error,
    relative_l2_error,
)
from lab_utils.results import save_json_result

__all__ = [
    "l2_norm_error",
    "max_abs_error",
    "relative_l2_error",
    "convergence_rate",
    "save_json_result",
]
