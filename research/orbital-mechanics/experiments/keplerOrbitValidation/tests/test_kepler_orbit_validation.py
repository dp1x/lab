"""Validation tests for the Kepler orbit validation experiment.

These must pass BEFORE any results are trusted (laboratory rule: verify before
trust). The tests check the propagator and the closed-form machinery against
Kepler's laws, conservation invariants, and the IAU real-units anchor, using
tolerances wide enough to be robust but narrow enough to catch real bugs.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_EXP_DIR = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "kepler_orbit_validation_experiment", _EXP_DIR / "experiment.py"
)
experiment = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(experiment)

A_DEFAULT = experiment.A_DEFAULT
AU_KM = experiment.AU_KM
EARTH_E = experiment.EARTH_E
K2_E = experiment.K2_E
K3_A = experiment.K3_A
K3_E = experiment.K3_E
MU = experiment.MU
MU_SUN_KM3_S2 = experiment.MU_SUN_KM3_S2
ORBITS_K1 = experiment.ORBITS_K1
SIDEREAL_YEAR_DAYS = experiment.SIDEREAL_YEAR_DAYS
STEPS_PER_ORBIT = experiment.STEPS_PER_ORBIT
acceleration = experiment.acceleration
circular_orbit_check = experiment.circular_orbit_check
conservation_study = experiment.conservation_study
earth_case = experiment.earth_case
initial_state = experiment.initial_state
kepler1_shape_check = experiment.kepler1_shape_check
kepler2_equal_areas = experiment.kepler2_equal_areas
kepler3_harmonic_law = experiment.kepler3_harmonic_law
kepler_solution = experiment.kepler_solution
measure_period = experiment.measure_period
orbital_elements = experiment.orbital_elements
propagate_rk4 = experiment.propagate_rk4
run_case = experiment.run_case
solve_kepler = experiment.solve_kepler

# --- Closed-form machinery ------------------------------------------------


def test_kepler_solution_satisfies_the_ode():
    """d^2 r / dt^2 = -mu r / r^3 on the analytic solution (central force).

    Fine grid: central-difference truncation error O(dt^2 * |a''|) would
    otherwise swamp the check near periapsis, where acceleration varies fast.
    """
    a, e, mu = 1.0, 0.6, MU
    t = np.linspace(0.0, 2 * np.pi, 60000)
    states = kepler_solution(a, e, mu, t)
    dt = t[1] - t[0]
    for i in range(2, len(states) - 2, 8):
        r = states[i, :2]
        acc_num = (states[i + 1, 2:] - states[i - 1, 2:]) / (2 * dt)
        acc_true = acceleration(r, mu)
        assert np.allclose(acc_num, acc_true, rtol=1e-3, atol=1e-8), f"i={i}"


def test_kepler_solution_periapsis_ic():
    a, e, mu = 1.0, 0.6, MU
    state = kepler_solution(a, e, mu, np.array([0.0]))[0]
    r0, v0 = state[:2], state[2:]
    assert r0[0] == pytest.approx(a * (1 - e), rel=1e-12)
    assert r0[1] == pytest.approx(0.0, abs=1e-15)
    v_p = np.sqrt(mu * (1 + e) / (a * (1 - e)))
    assert v0[0] == pytest.approx(0.0, abs=1e-15)
    assert v0[1] == pytest.approx(v_p, rel=1e-12)


def test_kepler_solution_conserves_specific_energy():
    a, e, mu = 1.0, 0.6, MU
    t = np.linspace(0.0, 10 * np.pi, 2000)
    states = kepler_solution(a, e, mu, t)
    r = np.hypot(states[:, 0], states[:, 1])
    eps = 0.5 * (states[:, 2] ** 2 + states[:, 3] ** 2) - mu / r
    assert np.allclose(eps, -mu / (2 * a), atol=1e-12)


def test_solve_kepler_reproduces_mean_anomaly():
    """KEPLER-LAW CHECK: E -> M = E - e sin E round trip."""
    e = 0.85
    E = np.linspace(0.0, 2 * np.pi, 500)
    M = E - e * np.sin(E)
    E_back = solve_kepler(M, e)
    assert np.allclose(E_back, E, atol=1e-12)


def test_initial_state_recovers_elements():
    for e in (0.0, 0.3, 0.6, 0.85):
        r0, v0 = initial_state(1.0, e, MU)
        els = orbital_elements(r0, v0, MU)
        assert els["eccentricity"] == pytest.approx(e, abs=1e-12)
        assert els["semi_major_axis"] == pytest.approx(1.0, rel=1e-12)
        assert els["energy"] == pytest.approx(-MU / 2.0, rel=1e-12)
        assert els["angular_momentum"] == pytest.approx(
            np.sqrt(MU * 1.0 * (1.0 - e**2)), rel=1e-12
        )


def test_vis_viva_matches_energy_at_periapsis():
    # Specific energy from the orbit's elements must match vis-viva at r_p.
    a, e, mu = 2.0, 0.5, MU
    r0, v0 = initial_state(a, e, mu)
    r_p = a * (1 - e)
    v_sq = np.hypot(*v0) ** 2
    eps = 0.5 * v_sq - mu / r_p
    assert eps == pytest.approx(-mu / (2 * a), rel=1e-12)


# --- Kepler I --------------------------------------------------------------


def test_rk4_matches_kepler_solution():
    """The propagator reproduces the analytic conic pointwise over 5 orbits."""
    t, states = run_case(1.0, 0.6, MU, ORBITS_K1, STEPS_PER_ORBIT)
    ana = kepler_solution(1.0, 0.6, MU, t)
    err = np.hypot(states[:, 0] - ana[:, 0], states[:, 1] - ana[:, 1])
    ana_r = np.hypot(ana[:, 0], ana[:, 1])
    assert float(np.max(err / ana_r)) < 1e-3


def test_conic_equation_holds_along_trajectory():
    """r = p / (1 + e cos nu) with (p, e) from the t=0 state (all e)."""
    for e in (0.3, 0.6, 0.85):
        case = kepler1_shape_check(1.0, e, MU)
        assert case["max_rel_conic_error"] < 1e-5, f"e={e}"
        assert case["e_from_eccentricity_vector"] == pytest.approx(e, abs=1e-6)


def test_circular_orbit_is_circle():
    case = circular_orbit_check(1.0, MU)
    assert case["max_rel_radius_variation"] < 1e-6
    assert case["measured_eccentricity"] < 1e-8


# --- Kepler II -------------------------------------------------------------


def test_equal_areas_in_equal_times():
    k = kepler2_equal_areas(1.0, K2_E, MU)
    assert k["max_rel_interval_deviation"] < 1e-3


def test_full_orbit_sector_area_is_pi_ab():
    k = kepler2_equal_areas(1.0, K2_E, MU)
    assert k["full_orbit_area_rel_error"] < 1e-3


def test_areal_velocity_is_h_over_2():
    k = kepler2_equal_areas(1.0, K2_E, MU)
    assert k["areal_velocity_rel_error"] < 1e-3


# --- Kepler III ------------------------------------------------------------


def test_kepler3_loglog_slope_is_3_over_2():
    k = kepler3_harmonic_law(MU)
    assert abs(k["loglog_slope"] - 1.5) < 0.01
    assert k["T2_over_a3_max_rel_err"] < 1e-3


def test_kepler3_period_matches_theory_all_cells():
    k = kepler3_harmonic_law(MU)
    for cell in k["cells"]:
        assert cell["period_rel_error"] < 1e-3, (
            f"a={cell['a']} e={cell['e']}: measured {cell['period_measured']}, "
            f"theory {cell['period_theory']}"
        )


def test_kepler3_subset_slope():
    """Cheap slope check on a subset (kept as a fast independent probe)."""
    a_s = [0.5, 2.0, 8.0]
    e_s = [0.10, 0.85]
    log_a, log_T = [], []
    for a in a_s:
        for e in e_s:
            t, states = run_case(a, e, MU, 2.05)
            t_meas = measure_period(t, states, e)
            log_a.append(np.log10(a))
            log_T.append(np.log10(t_meas))
    slope = np.polyfit(log_a, log_T, 1)[0]
    assert abs(slope - 1.5) < 0.01


# --- Conservation and convergence -----------------------------------------


def test_invariants_conserved_over_10_orbits():
    cons = conservation_study(1.0, 0.6, MU)
    assert cons["max_rel_energy_drift"] < 1e-6
    assert cons["max_rel_angular_momentum_drift"] < 1e-8


def test_propagator_convergence_order_is_four():
    conv = experiment.propagator_convergence(1.0, 0.6, MU)
    slopes = conv["measured_order"]
    # Pre-asymptotic slopes must approach 4; coarsest cells still feel O(h^5)
    # terms (measured 4.39, 4.24, 4.13, 4.07).
    assert all(abs(s - 4.0) < 0.5 for s in slopes), f"slopes: {slopes}"
    assert abs(slopes[-1] - 4.0) < 0.2


# --- Real-units Earth anchor -----------------------------------------------


def test_earth_static_constants():
    """Sanity on the verified constants themselves."""
    assert abs(MU_SUN_KM3_S2 - 1.3271244e11) < 1e3
    assert abs(AU_KM - 1.495978707e8) < 1.0
    assert 0.016 < EARTH_E < 0.018


def test_earth_period_matches_ia_u_constants():
    est = earth_case()
    assert est["period_rel_error"] < 1e-4  # numeric vs closed form
    assert abs(est["e_measured"] - EARTH_E) < 1e-3
    # The two-body prediction sits ~1.5e-6 above the sidereal year (Earth-Moon
    # system mass + nominal-vs-actual a). Verify the anchor lands in that band.
    diff = est["rel_diff_vs_sidereal_year"]
    assert 1.0e-6 < diff < 2.0e-6, f"diff {diff:.3e}"
    # Independent check of the closed form: T(days) ~= 2 pi sqrt(a^3/mu).
    pred = 2.0 * np.pi * np.sqrt(AU_KM**3 / MU_SUN_KM3_S2) / 86400.0
    assert pred == pytest.approx(est["period_theory_days"], rel=1e-12)
    assert pred == pytest.approx(365.256898, rel=1e-5)
    assert SIDEREAL_YEAR_DAYS == pytest.approx(365.256363, rel=1e-6)


# --- Determinism -----------------------------------------------------------


def test_determinism_across_processes():
    """A fresh interpreter must produce bit-identical numerical output.

    The subprocess loads the experiment module from its explicit path, so the
    check is independent of sys.path and of other experiments in the repo.
    """
    import json
    import subprocess
    import sys

    exp_dir = Path(__file__).resolve().parents[1]
    script = (
        "import importlib.util, json, sys\n"
        f"spec = importlib.util.spec_from_file_location('exp2', "
        f"{str(exp_dir / 'experiment.py')!r})\n"
        "m = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(m)\n"
        "payload = {\n"
        "    'k1': m.kepler1_shape_check(1.0, 0.6, 1.0),\n"
        "    'k2': m.kepler2_equal_areas(1.0, 0.6, 1.0),\n"
        "    't': m.measure_period(*m.run_case(1.0, 0.3, 1.0, 2.05), 0.3),\n"
        "}\n"
        "print(json.dumps(payload, sort_keys=True, default=float))\n"
    )
    out = subprocess.check_output([sys.executable, "-c", script], text=True).strip()
    here = {
        "k1": kepler1_shape_check(1.0, 0.6, 1.0),
        "k2": kepler2_equal_areas(1.0, 0.6, 1.0),
        "t": measure_period(*run_case(1.0, 0.3, 1.0, 2.05), 0.3),
    }
    expected = json.dumps(here, sort_keys=True, default=float)
    assert out == expected, "results differ between interpreter processes"