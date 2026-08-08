"""Comparative numerical integration of ODEs.

Study on the simple harmonic oscillator (x'' + omega^2 x = 0), which has a
closed-form analytic solution. Measures:

  1. Global convergence order of five integrators (error vs stepsize).
  2. Long-term energy preservation (drift of total mechanical energy).

Reference: J. C. Butcher, "Numerical Methods for Ordinary Differential
Equations", 3rd ed., Wiley, 2016. Energy behaviour of symplectic integrators:
E. Hairer, C. Lubich, G. Wanner, "Geometric Numerical Integration", 2nd ed.,
Springer, 2006.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lab_utils.metrics import convergence_rate, max_abs_error
from lab_utils.results import save_json_result

# --- Model definition -----------------------------------------------------

OMEGA = 1.0  # angular frequency [rad/s]
X0, V0 = 1.0, 0.0  # initial conditions

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"


def acceleration(x: float) -> float:
    """Acceleration of the harmonic oscillator."""
    return -OMEGA**2 * x


def analytic_solution(t: np.ndarray, x0: float = X0, v0: float = V0) -> np.ndarray:
    """Exact position x(t) = x0 cos(wt) + (v0/w) sin(wt)."""
    return x0 * np.cos(OMEGA * t) + (v0 / OMEGA) * np.sin(OMEGA * t)


def total_energy(x: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Mechanical energy E = 1/2 v^2 + 1/2 omega^2 x^2 (unit mass)."""
    return 0.5 * v**2 + 0.5 * OMEGA**2 * x**2


# --- Integrators ----------------------------------------------------------
# All integrators advance the state y = [x, v] from t0 with step h over
# n_steps, returning an array of shape (n_steps + 1, 2).


def forward_euler(t0: float, y0: np.ndarray, h: float, n_steps: int) -> np.ndarray:
    y = y0.astype(float).copy()
    states = np.empty((n_steps + 1, 2))
    states[0] = y
    t = t0
    for i in range(n_steps):
        y = y + h * np.array([y[1], acceleration(y[0])])
        states[i + 1] = y
        t += h
    return states


def rk2_midpoint(t0: float, y0: np.ndarray, h: float, n_steps: int) -> np.ndarray:
    y = y0.astype(float).copy()
    states = np.empty((n_steps + 1, 2))
    states[0] = y
    t = t0
    for i in range(n_steps):
        k1 = np.array([y[1], acceleration(y[0])])
        k2 = np.array(
            [y[1] + h / 2 * k1[1], acceleration(y[0] + h / 2 * k1[0])]
        )
        y = y + h * k2
        states[i + 1] = y
        t += h
    return states


def rk4(t0: float, y0: np.ndarray, h: float, n_steps: int) -> np.ndarray:
    y = y0.astype(float).copy()
    states = np.empty((n_steps + 1, 2))
    states[0] = y
    t = t0
    for i in range(n_steps):
        k1 = np.array([y[1], acceleration(y[0])])
        k2 = np.array([y[1] + h / 2 * k1[1], acceleration(y[0] + h / 2 * k1[0])])
        k3 = np.array([y[1] + h / 2 * k2[1], acceleration(y[0] + h / 2 * k2[0])])
        k4 = np.array([y[1] + h * k3[1], acceleration(y[0] + h * k3[0])])
        y = y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        states[i + 1] = y
        t += h
    return states


def symplectic_euler(t0: float, y0: np.ndarray, h: float, n_steps: int) -> np.ndarray:
    """Semi-implicit Euler: v updated with new x — symplectic, 1st order."""
    x, v = float(y0[0]), float(y0[1])
    states = np.empty((n_steps + 1, 2))
    states[0] = (x, v)
    for _ in range(n_steps):
        v = v + h * acceleration(x)
        x = x + h * v
        states[_ + 1] = (x, v)
    return states


def velocity_verlet(t0: float, y0: np.ndarray, h: float, n_steps: int) -> np.ndarray:
    """Velocity Verlet (2nd order, symplectic for separable Hamiltonians)."""
    x, v = float(y0[0]), float(y0[1])
    states = np.empty((n_steps + 1, 2))
    states[0] = (x, v)
    for i in range(n_steps):
        a_old = acceleration(x)
        x_new = x + v * h + 0.5 * a_old * h**2
        a_new = acceleration(x_new)
        v_new = v + 0.5 * (a_old + a_new) * h
        x, v = x_new, v_new
        states[i + 1] = (x, v)
    return states


METHODS: dict[str, callable] = {
    "euler": forward_euler,
    "rk2_midpoint": rk2_midpoint,
    "rk4": rk4,
    "symplectic_euler": symplectic_euler,
    "velocity_verlet": velocity_verlet,
}

# Theoretical global error orders.
THEORETICAL_ORDER = {
    "euler": 1.0,
    "rk2_midpoint": 2.0,
    "rk4": 4.0,
    "symplectic_euler": 1.0,
    "velocity_verlet": 2.0,
}

STEPSIZES = [0.1, 0.05, 0.025, 0.0125]
T_PERIOD = 2 * np.pi  # one full period
T_LONG = 200 * np.pi  # energy-drift horizon
H_ENERGY = 0.05  # stepsize for the energy study


def integrate(method: str, h: float, t_end: float) -> np.ndarray:
    n_steps = int(round(t_end / h))
    return METHODS[method](0.0, np.array([X0, V0]), h, n_steps)


def convergence_table() -> dict:
    """Global error at T_PERIOD vs analytic solution.

    The analytic reference is evaluated exactly at the integrator's time grid
    ``t_i = i*h`` (``t_0 = 0``), so no time-argument mismatch is measured.
    """
    errors: dict[str, list[float]] = {m: [] for m in METHODS}
    for h in STEPSIZES:
        for m in METHODS:
            states = integrate(m, h, T_PERIOD)
            t_num = np.arange(len(states)) * h
            errors[m].append(max_abs_error(states[:, 0], analytic_solution(t_num)))
    orders = {
        m: convergence_rate(np.array(errors[m]), np.array(STEPSIZES)).tolist()
        for m in METHODS
    }
    return {"stepsizes": STEPSIZES, "errors": errors, "measured_orders": orders}


def energy_study() -> dict:
    """Energy deviation |E - E0| over a long horizon for each method."""
    e0 = total_energy(X0, V0)
    t = np.arange(0.0, T_LONG + H_ENERGY, H_ENERGY)
    result = {}
    for m in METHODS:
        states = integrate(m, H_ENERGY, T_LONG)
        e = total_energy(states[:, 0], states[:, 1])
        dev = np.abs(e - e0)
        result[m] = {
            "max_deviation": float(np.max(dev)),
            "final_deviation": float(np.abs(e[-1] - e0)),
        }
    return {
        "initial_energy": float(e0),
        "horizon": float(T_LONG),
        "stepsize": H_ENERGY,
        "energy_deviations": result,
    }


def make_figures() -> list[str]:
    """Write convergence and energy figures; return their paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    for m in METHODS:
        errs = []
        for h in STEPSIZES:
            states = integrate(m, h, T_PERIOD)
            t_num = np.arange(len(states)) * h
            errs.append(max_abs_error(states[:, 0], analytic_solution(t_num)))
        axes[0].loglog(STEPSIZES, errs, marker="o", label=m)
    axes[0].set_xlabel("stepsize h")
    axes[0].set_ylabel("global error at t = T")
    axes[0].set_title("Convergence")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, which="both", alpha=0.3)

    e0 = total_energy(X0, V0)
    for m in ["euler", "rk4", "symplectic_euler", "velocity_verlet"]:
        states = integrate(m, H_ENERGY, T_LONG)
        e = total_energy(states[:, 0], states[:, 1])
        t_num = np.arange(len(e)) * H_ENERGY
        axes[1].semilogy(t_num, np.abs(e - e0) + 1e-16, label=m)
    axes[1].set_xlabel("time t")
    axes[1].set_ylabel("|E(t) - E(0)|")
    axes[1].set_title("Energy deviation, h = 0.05")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    (RESULTS_DIR / "figures").mkdir(parents=True, exist_ok=True)
    p1 = RESULTS_DIR / "figures" / "convergence.png"
    p2 = RESULTS_DIR / "figures" / "energy_deviation.png"
    fig.savefig(p1, dpi=150)
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    return [str(p1), str(p2)]


def main() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    conv = convergence_table()
    energy = energy_study()
    figures = make_figures()

    print("=== Convergence: global error at t = 2pi ===")
    print(f"{'method':<16}", *[f"h={h:<8}" for h in STEPSIZES], "order")
    for m in METHODS:
        avg_order = np.mean(conv["measured_orders"][m])
        print(
            f"{m:<16}",
            *[f"{e:.2e}" for e in conv["errors"][m]],
            f"{avg_order:.2f} (theory {THEORETICAL_ORDER[m]:.0f})",
        )
    print("\n=== Energy deviation over t in [0, 200pi], h = 0.05 ===")
    for m, d in energy["energy_deviations"].items():
        print(f"{m:<16} max={d['max_deviation']:.3e} final={d['final_deviation']:.3e}")

    result = {
        "convergence": conv,
        "energy_study": energy,
        "theoretical_orders": THEORETICAL_ORDER,
        "figures": [Path(f).name for f in figures],
    }
    path = save_json_result(
        RESULTS_DIR / "results.json",
        result,
        name="ode_integrator_comparison",
        description=(
            "Convergence order and energy preservation of five ODE integrators "
            "on the harmonic oscillator."
        ),
    )
    print(f"\nSaved results -> {path}")
    return result


if __name__ == "__main__":
    main()
