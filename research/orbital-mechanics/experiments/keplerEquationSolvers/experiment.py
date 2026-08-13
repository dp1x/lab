"""Kepler's equation solvers: Newton, bisection, and the Fourier-Bessel series.

Solves M = E - e sin E for the eccentric anomaly E, M in [0, 2*pi), e in [0,1):

  1. Newton iteration (starters M and M + e sin M) -- quadratic convergence.
  2. Bisection on the bracket [M, pi] (M <= pi) -- linear, guaranteed.
  3. Fourier-Bessel (Lagrange/Bessel) series E - M = sum (2/n) J_n(n e) sin(nM),
     whose per-parameter decay ratio is asymptotically
     q(e) = e * exp(sqrt(1-e^2)) / (1 + sqrt(1-e^2))  (Watson, Treatise on the
     Theory of Bessel Functions, Sec. 8.4 asymptotics) -- geometric but
     catastrophically slow as e -> 1.

The study measures convergence orders and rates, iteration counts vs
eccentricity, the series decay ratio vs the theoretical q(e), and agreement
among the three solvers. Independent cross-check: fixed-point iteration
E <- M + e sin E (linear convergence with asymptotic rate e cos E*, worst
case e at M = 0) as a fourth, independent reference.

References: R. Borghi, "The Kepler equation...", Mathematics 12(1):154, 2024
(arXiv:2312.01437); Philcox, Goodman & Slepian, "Kepler's Goat Herd", 2021
(arXiv:2103.15829); P. Colwell, "Solving Kepler's Equation over Three
Centuries", Willmann-Bell, 1993; G. N. Watson, "A Treatise on the Theory of
Bessel Functions", 2nd ed., Cambridge UP, 1944, Sec. 8.4; H. D. Curtis,
"Orbital Mechanics for Engineering Students", 4th ed., Elsevier, 2021, Ch. 3.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

from lab_utils.results import save_json_result

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

TOL = 1e-14  # residual tolerance for the root solvers
MAX_ITER_NEWTON = 50
MAX_ITER_BISECT = 80
MAX_TERMS = 2048  # series truncation cap for the study

# Eccentricities for the solver studies.
E_SWEEP = [0.10, 0.30, 0.60, 0.85, 0.90, 0.95, 0.99]
E_SERIES = [0.50, 0.70, 0.85, 0.90, 0.95]
N_SERIES = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
# Test points for the M grid (covers both half-planes, incl. near pi).
M_GRID = np.array(
    [0.0, 0.1, 0.3, 1.0, np.pi / 2, np.pi - 0.1, np.pi, np.pi + 0.1,
     3.0, 5.0, 2 * np.pi - 0.1, 2 * np.pi - 1e-9, 2 * np.pi]
)


def kepler_residual(E: np.ndarray, M: np.ndarray, e: float) -> np.ndarray:
    """Residual f(E) = E - e sin E - M (asserts the equation)."""
    return E - e * np.sin(E) - M


# --- Root solvers ----------------------------------------------------------


def solve_newton(M: float, e: float, starter: str = "msin") -> dict:
    """Newton iteration for M = E - e sin E.

    starter: "m" uses E0 = M; "msin" uses E0 = M + e sin M. f' = 1 - e cos E
    is positive for e < 1, so f is strictly increasing and Newton converges
    from any starter on either side of the root. The msin starter starts
    above the root for a majority of M in [0, pi] (and the m starter for all
    M > pi); convergence is unaffected.
    Returns E, iterations, function evaluations, residual history.
    """
    E = M + (e * np.sin(M) if starter == "msin" else 0.0)
    history: list[float] = []
    evals = 0
    for it in range(1, MAX_ITER_NEWTON + 1):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        history.append(abs(f))
        evals += 2
        if abs(f) < TOL:
            return {"E": E, "iterations": it, "evals": evals,
                    "history": history, "converged": True}
        E = E - f / fp
    raise RuntimeError(f"Newton did not converge (e={e}, M={M})")


def solve_bisection(M: float, e: float) -> dict:
    """Bisection on the canonical bracket.

    For M in [0, pi]: f(M) = -e sin M <= 0 and f(pi) = pi - M >= 0, so the
    root lies in [M, pi]; for M in (pi, 2pi) the symmetry E(2pi - M) =
    2pi - E(M) reduces to the same case. M = 0 -> E = 0, M = pi -> E = pi
    exactly. Returns E, iterations, evals, bracket-width history.
    """
    if M == 0.0:
        return {"E": 0.0, "iterations": 0, "evals": 0,
                "history": [], "converged": True}
    if M == np.pi:
        return {"E": np.pi, "iterations": 0, "evals": 0,
                "history": [], "converged": True}
    flipped = M > np.pi
    m = 2.0 * np.pi - M if flipped else M
    lo, hi = m, np.pi
    history: list[float] = []
    evals = 0
    for it in range(1, MAX_ITER_BISECT + 1):
        mid = 0.5 * (lo + hi)
        history.append(hi - lo)
        evals += 1
        if hi - lo < TOL:
            E = 0.5 * (lo + hi)
            if flipped:
                E = 2.0 * np.pi - E
            return {"E": E, "iterations": it, "evals": evals,
                    "history": history, "converged": True}
        if mid - e * np.sin(mid) - m > 0.0:
            hi = mid
        else:
            lo = mid
    raise RuntimeError(f"Bisection did not converge (e={e}, M={M})")


def solve_fixed_point(M: float, e: float, max_iter: int = 100000) -> dict:
    """Fixed-point iteration E <- M + e sin E (linear, rate e cos E* <= e).

    Independent reference solver: cheap, different convergence mechanism.
    The derivative of the iteration map at the root is e cos E*; the worst
    case (M = 0) gives rate e. None of the accuracy results depend on it.
    """
    E = M
    history: list[float] = []
    for it in range(1, max_iter + 1):
        E_new = M + e * np.sin(E)
        history.append(abs(E_new - E))
        if abs(E_new - E) < TOL:
            return {"E": E_new, "iterations": it, "evals": it,
                    "history": history, "converged": True}
        E = E_new
    raise RuntimeError(f"Fixed-point did not converge (e={e}, M={M})")


# --- Fourier-Bessel series -------------------------------------------------


def jn_miller(n: int, z: float) -> float:
    """J_n(z) by Miller backward recurrence (stable for z < n).

    The alternating power series J_n(z) = sum_k (-1)^k (z/2)^(n+2k)/(k!(n+k)!)
    suffers catastrophic cancellation at z ~ n (largest term ~ e^(c n) while
    J_n(n e) ~ e^(-c' n)); for z < n the backward recurrence
    J_{k-1} = (2k/z) J_k - J_{k+1} evaluated from a large index downward is
    numerically stable (DLMF 10.74(iii); Gil, Segura & Temme, Numerical
    Methods for Special Functions, 2007). The seed sequence is normalized
    with the exact sum identity J_0(z) + 2 sum_{k>=1} J_{2k}(z) = 1 (only
    even orders contribute, since J_{-k} = (-1)^k J_k). Values are rescaled
    in-flight so the raw peak (huge for z >> 60) cannot overflow.
    """
    M = n + 60 + int(2.0 * np.sqrt(n))  # depth of the recurrence tail
    jp = 0.0  # r[k+1]
    j = 1e-300  # r[k]
    s = 0.0  # accumulated normalized sum (same scale as the window)
    rn = 0.0  # r[n], captured when k - 1 == n
    for k in range(M, 0, -1):
        jm = (2.0 * k / z) * j - jp  # r[k-1]
        if abs(jm) > 1e200 or abs(j) > 1e200:
            f = 1e-200
            jm *= f
            j *= f
            jp *= f
            s *= f
            rn *= f
        idx = k - 1
        if idx == 0:
            s += jm
        elif idx % 2 == 0:
            s += 2.0 * jm
        if idx == n:
            rn = jm
        jp, j = j, jm
    if rn == 0.0:
        # True value is below the float64 range (captured raw value was
        # rescaled to zero); exact arithmetic would give < 1e-308 for z = ne
        # when e is small and n large, so 0.0 is the correct rounding.
        return 0.0
    return rn / s


@lru_cache(maxsize=256)
def series_coefficients(e: float, n_max: int) -> np.ndarray:
    """Fourier-Bessel coefficients c_n = (2/n) J_n(n e) for n = 1..n_max."""
    if e == 0.0:
        return np.zeros(n_max)
    coeffs = np.empty(n_max)
    for n in range(1, n_max + 1):
        coeffs[n - 1] = 2.0 * jn_miller(n, n * e) / n
    return coeffs


def solve_series(M: float, e: float, n_terms: int) -> dict:
    """E = M + sum_{n=1..n_terms} c_n sin(n M). Returns E and residual."""
    c = series_coefficients(e, n_terms)
    n = np.arange(1, n_terms + 1)
    delta = float(np.sum(c * np.sin(n * M)))
    E = M + delta
    return {"E": E, "n_terms": n_terms,
            "residual": abs(kepler_residual(E, M, e))}


def series_residuals(M_grid: np.ndarray, e: float, n_terms: int) -> np.ndarray:
    """Max absolute residual over the M grid for a fixed truncation."""
    c = series_coefficients(e, n_terms)
    n = np.arange(1, n_terms + 1)
    worst = 0.0
    for m in M_grid:
        E = m + float(np.sum(c * np.sin(n * m)))
        worst = max(worst, abs(kepler_residual(E, m, e)))
    return worst


def q_theory(e: float) -> float:
    """Asymptotic per-term decay ratio of the Fourier-Bessel series.

    q(e) = e * exp(sqrt(1-e^2)) / (1 + sqrt(1-e^2)), from J_n(n e) ~
    (2 pi sqrt(1-e^2) n)^(-1/2) q^n (Watson 1944, Sec. 8.4).
    """
    chi = np.sqrt(1.0 - e**2)
    return e * np.exp(chi) / (1.0 + chi)


# --- Studies ---------------------------------------------------------------


def _local_order(hist: list[float]) -> float:
    """Convergence order from residual steps clear of the round-off plateau.

    For order-p convergence |f_{k+1}| ~ C |f_k|^p, so p = log(r_{k+1}/r_k) /
    log(r_k/r_{k-1}). The very last residual usually sits on the round-off
    plateau (~1e-16 noise), so the order is read from the preceding step.
    Residuals that hit exact 0.0 are excluded.
    """
    pos = np.asarray(hist, dtype=float)
    pos = pos[pos > 0.0]
    if len(pos) < 3:
        return float("nan")
    l = np.log(pos)
    if len(pos) >= 4:
        return float((l[-2] - l[-3]) / (l[-3] - l[-4]))
    return float((l[-1] - l[-2]) / (l[-2] - l[-3]))


def newton_convergence_study() -> dict:
    """Residual history and measured convergence order per (e, starter)."""
    out = {}
    for e in (0.3, 0.6, 0.9):
        for starter in ("m", "msin"):
            res = solve_newton(1.0, e, starter)
            out[f"e={e}/starter={starter}"] = {
                "iterations": res["iterations"],
                "evals": res["evals"],
                "measured_order": _local_order(res["history"]),
                "history": [float(x) for x in res["history"]],
            }
    return out


def bisection_convergence_study() -> dict:
    """Bracket-width history: measured halving factor per iteration."""
    out = {}
    for e in (0.3, 0.6, 0.9):
        res = solve_bisection(1.0, e)
        widths = np.asarray(res["history"], dtype=float)
        ratios = widths[1:] / widths[:-1]
        out[f"e={e}"] = {
            "iterations": res["iterations"],
            "evals": res["evals"],
            "mean_halving_factor": float(np.mean(ratios)),
            "min_halving_factor": float(np.min(ratios)),
            "max_halving_factor": float(np.max(ratios)),
            "history": [float(x) for x in res["history"]],
        }
    return out


def eccentricity_sweep() -> dict:
    """Iterations to TOL across e, over the full M grid, all three solvers."""
    cells = []
    for e in E_SWEEP:
        row = {"e": e}
        max_iters = {}
        for starter in ("m", "msin"):
            worst = 0
            for m in M_GRID:
                worst = max(worst, solve_newton(m, e, starter)["iterations"])
            max_iters[f"newton_{starter}"] = worst
        worst_bis = max(solve_bisection(m, e)["iterations"] for m in M_GRID)
        max_iters["bisection"] = worst_bis
        row["worst_case_iterations"] = max_iters
        # Newton evaluates f and f' per iteration; bisection evaluates f once.
        row["worst_case_evals"] = {
            "newton_m": 2 * max_iters["newton_m"],
            "newton_msin": 2 * max_iters["newton_msin"],
            "bisection": max_iters["bisection"],
        }
        cells.append(row)
    return {"e_sweep": E_SWEEP, "m_grid": M_GRID.tolist(), "tol": TOL,
            "cells": cells}


def series_study() -> dict:
    """Residual vs truncation; measured coefficient-decay ratio vs q_theory."""
    cells = []
    for e in E_SERIES:
        residuals = [float(series_residuals(M_GRID, e, n)) for n in N_SERIES]
        # Watson's q is a property of the coefficients: |c_n| ~ C q^n (2 pi chi n)^-1/2.
        # The sqrt(n) prefactor is gone by n ~ 512; low-e tails underflow to
        # exact 0.0 beyond some n, so only strictly positive coefficients are
        # used in the geometric mean of consecutive ratios.
        c = series_coefficients(e, MAX_TERMS)
        tail = np.abs(c[511:])
        tail = tail[tail > 0.0]
        if len(tail) > 2:
            q_measured = float(np.exp(np.mean(np.log(tail[1:] / tail[:-1]))))
        else:
            q_measured = float("nan")
        cells.append({
            "e": e,
            "q_theory": q_theory(e),
            "q_measured": q_measured,
            "residuals": residuals,
            "n_terms": N_SERIES,
        })
    return {"n_terms": N_SERIES, "cells": cells}


def solver_agreement() -> dict:
    """Max disagreement between solvers over the M grid (per e)."""
    out = {}
    for e in (0.1, 0.3, 0.6, 0.85, 0.95):
        d_newton_bisect, d_series_newton, d_fixed_newton = 0.0, 0.0, 0.0
        for m in M_GRID:
            e_newton = solve_newton(m, e, "msin")["E"]
            e_bisect = solve_bisection(m, e)["E"]
            d_newton_bisect = max(d_newton_bisect, abs(e_newton - e_bisect))
            e_series = solve_series(m, e, 2048)["E"]
            d_series_newton = max(d_series_newton, abs(e_series - e_newton))
            e_fixed = solve_fixed_point(m, e)["E"]
            d_fixed_newton = max(d_fixed_newton, abs(e_fixed - e_newton))
        out[f"e={e}"] = {
            "max_|newton-bisection|": d_newton_bisect,
            "max_|series2048-newton|": d_series_newton,
            "max_|fixedpoint-newton|": d_fixed_newton,
        }
    return out


def special_values() -> dict:
    """Closed-form anchors: E(0) = 0, E(pi) = pi, E(2pi) = 2pi for all e."""
    out = {}
    for e in E_SWEEP:
        out[f"e={e}"] = {
            "E(0.0)": solve_newton(0.0, e, "msin")["E"],
            "E(pi)": solve_newton(np.pi, e, "msin")["E"],
            "E(2pi)": solve_newton(2.0 * np.pi, e, "msin")["E"],
        }
    return out


# --- Figures ---------------------------------------------------------------


def make_figures(newton: dict, bisect: dict, series: dict) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (RESULTS_DIR / "figures").mkdir(parents=True, exist_ok=True)
    paths = []

    # 1. Newton residual histories (quadratic convergence).
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for key, data in newton.items():
        hist = np.clip(np.asarray(data["history"]), 1e-300, None)
        ax.semilogy(np.arange(1, len(hist) + 1), hist, "o-", ms=3, label=key)
    ax.set_xlabel("iteration k")
    ax.set_ylabel("|E_k - e sin E_k - M|")
    ax.set_title("Newton: quadratic convergence (M = 1.0 rad)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "newton_convergence.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 2. Bisection bracket-width histories (linear halving).
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for key, data in bisect.items():
        widths = np.asarray(data["history"])
        ax.semilogy(np.arange(1, len(widths) + 1), widths, "-", lw=1.2,
                    label=f"{key} (factor {data['mean_halving_factor']:.4f})")
    ax.axhline(1e-14, color="k", ls=":", lw=0.8)
    ax.set_xlabel("iteration k")
    ax.set_ylabel("bracket width")
    ax.set_title("Bisection: linear convergence (M = 1.0 rad)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "bisection_convergence.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    # 3. Series residual vs truncation and q(e) fit vs theory.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.5))
    for cell in series["cells"]:
        e = cell["e"]
        ax1.loglog(N_SERIES, cell["residuals"], "o-", ms=3, label=f"e = {e}")
    ax1.set_xlabel("number of series terms N")
    ax1.set_ylabel("max |residual| over M grid")
    ax1.set_title("Fourier-Bessel series: slow geometric decay")
    ax1.grid(True, which="both", alpha=0.3)
    es = [c["e"] for c in series["cells"]]
    qm = [c["q_measured"] for c in series["cells"]]
    qt = [c["q_theory"] for c in series["cells"]]
    ax2.plot(es, qt, "k--", label="q_theory(e) = e e^{chi}/(1+chi)")
    ax2.plot(es, qm, "ro-", ms=4, label="q_measured (coefficient tail mean)")
    ax2.set_xlabel("eccentricity e")
    ax2.set_ylabel("series decay ratio")
    ax2.set_title("Measured decay ratio vs Watson Sec. 8.4 asymptotics")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    p = RESULTS_DIR / "figures" / "series_convergence.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p))

    return paths


# --- Main ------------------------------------------------------------------


def main() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    newton = newton_convergence_study()
    bisect = bisection_convergence_study()
    sweep = eccentricity_sweep()
    series = series_study()
    agree = solver_agreement()
    special = special_values()
    figures = make_figures(newton, bisect, series)

    print("=== Newton convergence (M = 1.0 rad) ===")
    for key, data in newton.items():
        print(f"{key:<24} iters={data['iterations']:>2}  order={data['measured_order']:.3f}")
    print("=== Bisection convergence (M = 1.0 rad) ===")
    for key, data in bisect.items():
        print(f"{key:<12} iters={data['iterations']:>2}  halving={data['mean_halving_factor']:.5f}")
    print("=== Worst-case iterations over M grid ===")
    for row in sweep["cells"]:
        w = row["worst_case_iterations"]
        print(
            f"e={row['e']:<5} newton_m={w['newton_m']:>2} newton_msin={w['newton_msin']:>2}"
            f" bisection={w['bisection']:>2}"
        )
    print("=== Series decay ratio vs theory ===")
    for cell in series["cells"]:
        qm = cell["q_measured"]
        qs = f"{qm:.6f}" if qm is not None else "n/a"
        print(f"e={cell['e']:<5} q_meas={qs}  q_theory={cell['q_theory']:.6f}")
    print("=== Solver agreement over M grid ===")
    for key, d in agree.items():
        print(
            f"{key:<6} |newton-bisection|={d['max_|newton-bisection|']:.2e}"
            f" |series-newton|={d['max_|series2048-newton|']:.2e}"
            f" |fixed-newton|={d['max_|fixedpoint-newton|']:.2e}"
        )
    print("=== Special values (E(0), E(pi), E(2pi)) ===")
    for key, d in special.items():
        print(f"{key:<8} {d['E(0.0)']:.6f} {d['E(pi)']:.6f} {d['E(2pi)']:.6f}")

    result = {
        "tol": TOL,
        "newton_study": newton,
        "bisection_study": bisect,
        "eccentricity_sweep": sweep,
        "series_study": series,
        "solver_agreement": agree,
        "special_values": special,
        "figures": [Path(p).name for p in figures],
    }
    path = save_json_result(
        RESULTS_DIR / "results.json",
        result,
        name="kepler_equation_solvers",
        description=(
            "Newton, bisection and Fourier-Bessel series solvers for Kepler's "
            "equation: convergence orders, iteration counts, series decay "
            "ratio vs theory, cross-solver agreement."
        ),
    )
    print(f"\nSaved results -> {path}")
    return result


if __name__ == "__main__":
    main()