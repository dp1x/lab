"""Validation tests for the Kepler's equation solvers experiment.

These must pass BEFORE any results are trusted (laboratory rule: verify before
trust). The experiment module is loaded via importlib from its explicit path
so multiple experiments with an "experiment.py" module never collide in
pytest/sys.modules (see tools/new_experiment.py).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "kepler_equation_solvers_experiment", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)

solve_newton = experiment.solve_newton
solve_bisection = experiment.solve_bisection
solve_fixed_point = experiment.solve_fixed_point
solve_series = experiment.solve_series
series_coefficients = experiment.series_coefficients
jn_miller = experiment.jn_miller
q_theory = experiment.q_theory
kepler_residual = experiment.kepler_residual
M_GRID = experiment.M_GRID


# --- Closed-form anchors ---------------------------------------------------


def test_special_values():
    """E(0) = 0, E(pi) = pi, E(2pi) = 2pi exactly for all e."""
    for e in (0.0, 0.3, 0.6, 0.9, 0.99):
        for m, expected in ((0.0, 0.0), (np.pi, np.pi), (2 * np.pi, 2 * np.pi)):
            E = solve_newton(m, e, "msin")["E"]
            assert abs(E - expected) < 1e-14, f"e={e}, M={m}: E={E}"


@pytest.mark.parametrize("e", [0.1, 0.3, 0.6, 0.85, 0.95])
@pytest.mark.parametrize("m", [0.1, 1.0, np.pi / 2, np.pi - 0.1, 3.0, 5.5])
def test_all_solvers_reach_machine_precision(e, m):
    """Every solver leaves |E - e sin E - M| below 1e-13."""
    for solver in (lambda: solve_newton(m, e, "m"),
                   lambda: solve_newton(m, e, "msin"),
                   lambda: solve_bisection(m, e),
                   lambda: solve_fixed_point(m, e)):
        assert solver()["converged"]
        E = solver()["E"]
        assert abs(kepler_residual(E, m, e)) < 1e-13, f"e={e}, M={m}, E={E}"


def test_newton_quadratic_order():
    """Newton's residual steps clear of the plateau must show p ~ 2."""
    for e in (0.3, 0.6, 0.9):
        res = solve_newton(1.0, e, "msin")
        p = experiment._local_order(res["history"])
        assert np.isfinite(p), f"e={e}: unmeasurable"
        assert abs(p - 2.0) < 0.6, f"e={e}: measured order {p:.3f}"


def test_bisection_halving_factor():
    """Bisection bracket width must halve each iteration (factor 0.5)."""
    for e in (0.3, 0.6, 0.9):
        res = solve_bisection(1.0, e)
        widths = np.asarray(res["history"])
        ratios = widths[1:] / widths[:-1]
        assert np.all(np.abs(ratios - 0.5) < 0.01), f"e={e}: {ratios}"


def test_bisection_iterations_match_halving_theory():
    """Width pi must shrink below 1e-14 in ceil(log2(pi/1e-14)) ~ 49 steps."""
    for e in (0.3, 0.6, 0.9):
        it = solve_bisection(1.0, e)["iterations"]
        assert 44 <= it <= 60, f"e={e}: iterations={it}"


def test_fixed_point_rate_is_e_times_cos_E():
    """Fixed-point rate is e cos E* (linear); worst case is e at M = 0."""
    e, m = 0.5, 1.0
    res = solve_fixed_point(m, e, max_iter=20000)
    h = np.asarray(res["history"])
    ratios = h[10:] / h[9:-1]
    estar = solve_newton(m, e, "msin")["E"]
    expected = e * np.cos(estar)
    measured = float(np.mean(ratios))
    assert measured == pytest.approx(expected, rel=0.05)
    assert measured < e


# --- Fourier-Bessel series -------------------------------------------------


def test_series_first_coefficients():
    """c_1 = 2 J_1(e)/1, and for e -> 0 the series reduces to E ~ M + e sin M."""
    for e in (0.1, 0.5, 0.9):
        c = series_coefficients(e, 3)
        assert c[0] == pytest.approx(2.0 * _j1(e), rel=1e-10)


def _j1(z: float) -> float:
    """Independent J_1 via its power series (straightforward floats)."""
    s = 0.0
    for k in range(60):
        term = ((-1) ** k) * (z / 2) ** (1 + 2 * k) / (
            math.factorial(k) * math.factorial(k + 1)
        )
        s += term
        if abs(term) < 1e-18 * abs(s):
            break
    return s


def test_series_coefficients_match_published_values():
    """c_n vs published J_n constants (DLMF 10.2 / Abramowitz-Stegun)."""
    c = series_coefficients(0.5, 2)
    # c_1 = 2 J_1(0.5), c_2 = J_2(1.0); values from DLMF Table 10.2.
    assert c[0] == pytest.approx(2.0 * 0.24226845767487388, rel=1e-9)
    assert c[1] == pytest.approx(0.11490348493190048, rel=1e-9)
    j0 = jn_miller(0, 1.0)
    assert j0 == pytest.approx(0.76519768655796655, rel=1e-9)


@pytest.mark.parametrize("e", [0.3, 0.6, 0.85, 0.9])
def test_series_matches_newton_at_large_truncation(e):
    for m in (0.1, 1.0, np.pi - 0.1, 5.0):
        E_s = solve_series(m, e, 2048)["E"]
        E_n = solve_newton(m, e, "msin")["E"]
        assert abs(E_s - E_n) < 1e-9, f"e={e}, M={m}: {abs(E_s - E_n):.2e}"


def test_series_decay_matches_theory():
    """Measured per-term decay ratio equals q_theory within 15%."""
    for e in (0.5, 0.7, 0.85):
        cells = experiment.series_study()["cells"]
        cell = next(c for c in cells if c["e"] == e)
        assert cell["q_measured"] is not None
        assert abs(cell["q_measured"] - cell["q_theory"]) / cell["q_theory"] < 0.15, (
            f"e={e}: measured {cell['q_measured']:.5f}, theory {cell['q_theory']:.5f}"
        )


def test_q_theory_trend():
    """q(e) must increase toward 1 as e -> 1 (slow series near parabolic)."""
    qs = [q_theory(e) for e in (0.5, 0.7, 0.85, 0.9, 0.95, 0.99)]
    assert all(b > a for a, b in zip(qs, qs[1:]))
    assert qs[-1] > 0.998


def test_series_needs_many_terms_at_high_e():
    """At e = 0.9 the series still has ~1e-3 residual at N = 64 (slow)."""
    res_64 = experiment.series_residuals(M_GRID, 0.9, 64)
    res_2048 = experiment.series_residuals(M_GRID, 0.9, 2048)
    assert res_64 > 1e-5  # genuinely slow at high e
    assert res_2048 < 1e-9  # but eventually converges


# --- Agreement across solvers ---------------------------------------------


def test_cross_solver_agreement_full_grid():
    agree = experiment.solver_agreement()
    for key, d in agree.items():
        assert d["max_|newton-bisection|"] < 1e-12, key
        assert d["max_|fixedpoint-newton|"] < 1e-11, key
    # Series at e = 0.95 and N = 2048 is residual-limited (~1e-8); looser.
    assert agree["e=0.95"]["max_|series2048-newton|"] < 1e-7


# --- Determinism -----------------------------------------------------------


def test_determinism_across_processes():
    """A fresh interpreter must produce bit-identical numerical output."""
    import json
    import subprocess
    import sys

    exp_dir = Path(__file__).resolve().parents[1]
    script = (
        "import importlib.util, json\n"
        f"spec = importlib.util.spec_from_file_location('exp3', "
        f"{str(exp_dir / 'experiment.py')!r})\n"
        "m = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(m)\n"
        "payload = {\n"
        "    'n': m.solve_newton(1.0, 0.85, 'msin'),\n"
        "    's': m.series_study(),\n"
        "    'b': m.solve_bisection(2.9, 0.9),\n"
        "}\n"
        "print(json.dumps(payload, sort_keys=True, default=float))\n"
    )
    out = subprocess.check_output([sys.executable, "-c", script], text=True).strip()
    here = {
        "n": solve_newton(1.0, 0.85, "msin"),
        "s": experiment.series_study(),
        "b": solve_bisection(2.9, 0.9),
    }
    expected = json.dumps(here, sort_keys=True, default=float)
    assert out == expected, "results differ between interpreter processes"