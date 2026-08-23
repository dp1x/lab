"""lab_utils — shared utilities for the Research Lab.

Contains only reusable code: metrics, integrators, orbital canon, and result I/O.
Each module is self-contained and covered by tests in `src/lab_utils/tests/`.
"""

from lab_utils.integrators import rk4_propagate, rk4_step
from lab_utils.metrics import (
    convergence_rate,
    l2_norm_error,
    max_abs_error,
    relative_l2_error,
)
from lab_utils.orbits import (
    ECC_GUARD_ABS,
    J2_EARTH,
    MU_EARTH_KM3S2,
    NODE_GUARD_REL,
    OMEGA_EARTH_RAD_S,
    R_EARTH_KM,
    coe_to_rv_eci,
    mean_motion,
    orbital_period,
    rotation_matrix_313,
    rv_to_coe_eci,
    seed_state,
    solve_kepler,
    steps_per_orbit,
    true_anomaly_from_E,
    j2_rhs,
)
from lab_utils.results import save_json_result

__all__ = [
    "l2_norm_error",
    "max_abs_error",
    "relative_l2_error",
    "convergence_rate",
    "save_json_result",
    "rk4_step",
    "rk4_propagate",
    "MU_EARTH_KM3S2",
    "R_EARTH_KM",
    "OMEGA_EARTH_RAD_S",
    "J2_EARTH",
    "NODE_GUARD_REL",
    "ECC_GUARD_ABS",
    "solve_kepler",
    "true_anomaly_from_E",
    "orbital_period",
    "mean_motion",
    "rotation_matrix_313",
    "coe_to_rv_eci",
    "rv_to_coe_eci",
    "seed_state",
    "steps_per_orbit",
    "j2_rhs",
]
