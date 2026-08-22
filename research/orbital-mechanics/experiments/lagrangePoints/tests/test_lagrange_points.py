"""Tests for Experiment 011: Lagrange points / CR3BP (rotating-frame dynamics).

Doctrine follows Exp 009/010: theory constants and oracle formulas are duplicated
inline (NEVER imported from the experiment module), discriminators are reimplemented
from first principles where the anti-shared-algebra rule demands it, and heavy runs
are cached at module level. Banners L1..L7 mirror the pre-registered failure catalog
(Track E): L1 frames/identities, L2 equilibria/pinned values, L3 propagation-vs-
reference, L4 conservation/invariants, L5 stability + nonlinear perturbation,
L6 dimensional cross-check/determinism, L7 adversarial mutants.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Module loading (experiment under test)
# --------------------------------------------------------------------------- #
_EXP_PATH = Path(__file__).resolve().parents[1] / "experiment.py"
_spec = importlib.util.spec_from_file_location("lagrange_points_experiment", _EXP_PATH)
assert _spec is not None and _spec.loader is not None
experiment = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(experiment)

# --------------------------------------------------------------------------- #
# Inline independent theory (duplicated on purpose -- never imported)
# --------------------------------------------------------------------------- #
GM_E = 398600.4418
GM_M = 4902.800066
MU_EM = GM_M / (GM_E + GM_M)  # bit-exact MU-firewall reference
A_EM_KM = 384400.0
AU_KM = 149597870.7
MU_ROUTH = (9.0 - math.sqrt(69.0)) / 18.0


def grad_om_inline(x, y, z, mu):
    """Independent transcription of grad(omega_eff)."""
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    return (
        x - (1 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3,
        y - (1 - mu) * y / r1**3 - mu * y / r2**3,
        -(1 - mu) * z / r1**3 - mu * z / r2**3,
    )


def jacobi_inline(state, mu):
    x, y, z, vx, vy, vz = state
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    om = (1 - mu) / r1 + mu / r2 + 0.5 * (x * x + y * y)
    return 2.0 * om - (vx * vx + vy * vy + vz * vz)


def rot_z_inline(th):
    c, s = math.cos(th), math.sin(th)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rhs_rot_inline(mu):
    def f(t, s):
        gx, gy, gz = grad_om_inline(s[0], s[1], s[2], mu)
        return np.array([s[3], s[4], s[5], 2.0 * s[4] + gx, -2.0 * s[3] + gy, gz])

    return f


def rhs_inert_inline(mu):
    def f(t, s):
        R = rot_z_inline(t)
        p1 = R @ np.array([-mu, 0.0, 0.0])
        p2 = R @ np.array([1.0 - mu, 0.0, 0.0])
        d1 = s[:3] - p1
        d2 = s[:3] - p2
        a = -(1 - mu) * d1 / np.linalg.norm(d1) ** 3 - mu * d2 / np.linalg.norm(d2) ** 3
        return np.concatenate([s[3:], a])

    return f


def quintic_gamma(name, x, mu):
    if name == "L1":
        g = 1.0 - mu - x
        return g**5 - (3 - mu) * g**4 + (3 - 2 * mu) * g**3 - mu * g**2 + 2 * mu * g - mu
    if name == "L2":
        g = x - (1.0 - mu)
        return g**5 + (3 - mu) * g**4 + (3 - 2 * mu) * g**3 - mu * g**2 - 2 * mu * g - mu
    g = -mu - x
    return g**5 + (2 + mu) * g**4 + (1 + 2 * mu) * g**3 - (1 - mu) * g**2 - 2 * (1 - mu) * g - (1 - mu)


def bisect_inline(f, a, b, n=300):
    fa = f(a)
    for _ in range(n):
        m = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0 or (b - m) < 1e-18:
            return m
        if (fa < 0.0) != (fm < 0.0):
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


# --------------------------------------------------------------------------- #
# Module-level cached heavy runs (computed once per pytest session)
# --------------------------------------------------------------------------- #
_CACHE: dict = {}


def _cached(key, builder):
    if key not in _CACHE:
        _CACHE[key] = builder()
    return _CACHE[key]


def _mut_study():
    return _cached("mut", experiment.study_inertial_consistency_and_mutants)


def _jac_study():
    return _cached("jac", lambda: experiment.study_jacobi_drift()[0])


def _pert_study():
    return _cached("pert", experiment.study_perturbations)


def _dim_study():
    return _cached("dim", experiment.study_dimensional_crosscheck)


def _brid_study():
    return _cached("brid", experiment.study_bridge_and_drift_law)


def _sym_study():
    return _cached("sym", experiment.study_symmetries_and_limits)


def _eq_study():
    return _cached("eq", experiment.study_equilibria)


# =========================================================================== #
# L1 -- frame transformations & physical identities
# =========================================================================== #
class TestL1FramesAndIdentities:
    def test_round_trip_identity_several_angles(self):
        rng_state = [0.31, -0.52, 0.11, 0.07, 0.23, -0.05]
        worst = 0.0
        for th in (0.0, 0.3, 1.0, 2.5, math.pi, 4.7):
            s = np.array(rng_state)
            back = experiment.inert_to_rot(experiment.rot_to_inert(s, th), th)
            worst = max(worst, float(np.max(np.abs(back - s))))
        assert worst < 1e-14

    def test_rotation_matrix_handedness(self):
        R = rot_z_inline(0.7)
        assert abs(np.linalg.det(R) - 1.0) < 1e-15
        assert np.max(np.abs(R.T @ R - np.eye(3))) < 1e-15

    def test_angular_velocity_of_frame_consistent(self):
        """A point fixed in the rotating frame moves inertially as d/dt[R(th)r]
        = R(th) (omega x r): finite difference must match the rotated omega x r."""
        r_fixed = np.array([0.6, 0.25, 0.0])
        th, h = 0.9, 1e-6
        v_fd = (rot_z_inline(th + h) @ r_fixed - rot_z_inline(th - h) @ r_fixed) / (2 * h)
        v_exact = rot_z_inline(th) @ np.cross([0.0, 0.0, 1.0], r_fixed)
        assert np.max(np.abs(v_fd - v_exact)) < 1e-8

    def test_coriolis_sign_isolated(self):
        """Acceleration difference between moving and rest particle isolates
        -2 omega x v: moving +x must give a_y = -2u exactly."""
        mu = MU_EM
        p = [0.7, 0.3, 0.1]
        u = 0.137
        fr = experiment.rhs_rotating(mu)
        a_v = fr(0.0, np.array([*p, u, 0.0, 0.0]))[3:]
        a_0 = fr(0.0, np.array([*p, 0.0, 0.0, 0.0]))[3:]
        diff = a_v - a_0
        assert abs(diff[0]) < 1e-14 and abs(diff[2]) < 1e-14
        assert abs(diff[1] - (-2.0 * u)) < 1e-13

    def test_centrifugal_outward_at_rest(self):
        mu = 1e-6
        g = experiment.grad_omega(10.0, 0.0, 0.0, mu)
        assert g[0] > 9.9  # outward ~ +10

    def test_L4_L5_exact_coordinates_and_residual(self):
        pts = experiment.lagrange_points(MU_EM)
        assert abs(pts["L4"][0] - (0.5 - MU_EM)) < 1e-17
        assert abs(pts["L5"][1] + math.sqrt(3) / 2) < 1e-16
        for name in ("L4", "L5"):
            gx, gy, gz = grad_om_inline(pts[name][0], pts[name][1], pts[name][2], MU_EM)
            assert max(abs(gx), abs(gy), abs(gz)) < 1e-14


# =========================================================================== #
# L2 -- equilibria: brackets, ordering, pinned values, anchors
# =========================================================================== #
class TestL2Equilibria:
    def test_proven_brackets_and_ordering_multiple_mu(self):
        for mu in (MU_EM, 3.040423452319562e-06, 0.03852, 0.05, 1e-3, 1e-6):
            r = experiment.solve_collinear(mu)
            assert 0.5 - mu <= r["L1"] <= 1.0 - mu
            assert 1.0 - mu <= r["L2"] <= 2.0 - mu
            assert -1.0 - mu <= r["L3"] <= -1.0
            assert r["L3"] < -mu < r["L1"] < 1.0 - mu < r["L2"]

    def test_quintic_residuals_independent_family(self):
        for mu in (MU_EM, 0.05, 0.3):
            r = experiment.solve_collinear(mu)
            scale = {"L1": 0.03, "L2": 0.04, "L3": 8.0}
            for name in ("L1", "L2", "L3"):
                q = quintic_gamma(name, r[name], mu)
                assert abs(q) < 1e-13 * scale[name], f"{name} mu={mu}: {q}"

    def test_vector_residual_all_points_all_cases(self):
        for case, mu in experiment.CASES.items():
            pts = experiment.lagrange_points(mu)
            for name, pos in pts.items():
                gx, gy, gz = grad_om_inline(pos[0], pos[1], pos[2], mu)
                assert max(abs(gx), abs(gy), abs(gz)) < 2e-14, f"{case}/{name}"

    def test_mpmath_anchor_agreement_em(self):
        eq = _eq_study()
        dev = eq["mpmath_anchor_dev"]
        for name, d in dev.items():
            assert d < 5e-16, f"{name}: {d:.3e}"

    def test_jacobi_critical_value_ordering(self):
        assert _eq_study()["c_ordering_ok"]

    def test_mu_firewall_bit_exact(self):
        """mu must equal GM_M/(GM_E+GM_M) bit-exactly; L4 x must be 0.5-mu exactly."""
        assert experiment.MU_EM == MU_EM
        pts = experiment.lagrange_points(MU_EM)
        assert pts["L4"][0] == 0.5 - MU_EM
        assert pts["L5"][0] == 0.5 - MU_EM

    def test_mission_anchors_percent_level(self):
        an = _eq_study()["anchors"]
        assert abs(an["em_L1_km_from_moon"] - 58000.0) / 58000.0 < 0.02
        assert abs(an["em_L2_km_from_moon"] - 64500.0) / 64500.0 < 0.02
        assert abs(an["em_L3_km_from_earth"] - 381700.0) / 381700.0 < 0.001
        assert abs(an["sem_L1_km_from_earth"] - 1.5e6) / 1.5e6 < 0.03
        assert abs(an["sem_L2_km_from_earth"] - 1.5e6) / 1.5e6 < 0.03

    def test_hill_asymptotics_and_L3_offset_law(self):
        sym = _sym_study()
        for key in ("mu_0.001", "mu_1e-06"):
            lim = sym["singular_limits"][key]
            assert abs(lim["g1_over_series"] - 1.0) < 1e-3
            assert abs(lim["g2_over_series"] - 1.0) < 1e-3
            assert abs(lim["xL3_plus_1_over_mu"] - (-5.0 / 12.0)) < 1e-6
        slopes = sym["hill_slope_convergence"]["pairwise_orders"]
        # local exponent approaches 1/3 from below; along ASCENDING mu the series
        # correction gamma1 = alpha(1 - alpha/3 - ...) pushes the exponent down
        assert abs(slopes[0] - 1.0 / 3.0) < 0.005
        assert all(slopes[i + 1] < slopes[i] for i in range(len(slopes) - 1))

    def test_equal_mass_symmetry(self):
        sym = _sym_study()["equal_mass"]
        assert sym["x_L1_exact_zero"] == 0.0
        assert abs(sym["sym_x2_plus_x3"]) < 1e-12
        assert abs(sym["L4_x"]) < 1e-16


# =========================================================================== #
# L3 -- propagation vs independent references
# =========================================================================== #
class TestL3PropagationReferences:
    def test_inertial_consistency_round_trip_clean(self):
        """THE centerpiece discriminator: rotating-frame propagation mapped to the
        inertial frame must satisfy Newton's law with moving primaries."""
        res = _mut_study()
        assert res["clean_relative_residual"] < 1e-10

    def test_kinematic_closed_form_order(self):
        order = experiment.study_rk4_order_verification()
        assert all(o > 3.7 for o in order["pairwise_orders"])

    def test_mirror_field_law(self):
        """f(Mw) = -M f(w) with M = flip(y, vx, vz) must hold to FP noise."""
        assert _sym_study()["mirror_field_law_max_violation"] < 1e-13

    def test_half_mass_pi_rotation_equivariance(self):
        assert _sym_study()["half_mass_pi_rotation_max_violation"] < 1e-13


# =========================================================================== #
# L4 -- conservation laws & invariants
# =========================================================================== #
class TestL4Conservation:
    def test_capture_class_drift_orders_above_floor(self):
        jac = _jac_study()
        cap = jac["class_capture_C2minus1e-3"]
        drifts = list(cap["drifts"].values())
        floor_ref = 20 * np.finfo(float).eps * abs(3.19)
        pairs = [
            o
            for o, d_next in zip(cap["orders"], drifts[1:])
            if d_next > floor_ref
        ]
        assert len(pairs) >= 3
        assert all(3.5 < o < 5.0 for o in pairs)

    def test_retrograde_floor_approach(self):
        fp = _jac_study()["class_retrograde_rho0p30"]["floor_probe"]
        vals = list(fp.values())
        assert max(vals) < 1e-11 and min(vals) > 0.0

    def test_bounded_orbit_reaches_quantization_floor_plateau(self):
        plat = _jac_study()["class_L4_LP_floor_plateau"]
        drifts = list(plat["drifts"].values())
        # quantization floor lives on the finest rungs: equal to within 1 ulp of C
        assert max(drifts[-2:]) < 5e-15
        assert abs(drifts[-1] - drifts[-2]) < 1e-15
        assert plat["plateau_max_over_min"] > 1.0

    def test_bridge_identity_and_energy_frame_law(self):
        brid = _brid_study()["n_16384"]
        assert brid["bridge_identity_max"] < 1e-13
        assert brid["identity_E_plus_C2_minus_Lz_max"] < 1e-14
        # E_I itself is NOT conserved (moving primaries do net work)
        assert brid["EI_drift_measured"] != 0.0

    def test_ei_drift_law_matches_integral(self):
        brid = _brid_study()
        assert brid["n_4096"]["drift_law_rel_mismatch"] < 1e-4
        ratio = brid["drift_law_mismatch_reduction_ratio"]
        assert ratio > 8.0  # ~h^4 defect under dt/4 refinement

    def test_spatial_jacobi_true_vs_frozen_evaluator(self):
        """Kills the planar-frozen-Jacobi class: true spatial conservation holds
        while a z-frozen evaluator leaks orders of magnitude more."""
        mut = _mut_study()
        leak = mut["mutant_signals"]["planar_frozen_jacobi"]
        true_drift = mut["spatial_true_drift"]
        assert true_drift < 1e-12
        assert leak > 1e-6


# =========================================================================== #
# L5 -- linear stability + nonlinear perturbation signatures
# =========================================================================== #
class TestL5Stability:
    def test_collinear_closed_form_vs_numeric_em(self):
        stab = experiment.study_stability()["cases"]["earth_moon"]
        for name in ("L1", "L2", "L3"):
            entry = stab[name]
            assert entry["closed_vs_numeric_rel"] < 1e-12
            assert entry["max_real_eigenvalue"] > 1e-6  # unstable for every mu

    def test_triangular_rates_vs_numeric_em(self):
        stab = experiment.study_stability()["cases"]["earth_moon"]["L4"]
        imag = sorted(abs(x) for x in stab["nu_numeric_imag"])
        # spectrum = {+-i*nu_long (x2), +-i*nu_short (x2), +-i (vertical, x2)}
        assert abs(imag[0] - stab["nu_long_closed"]) / stab["nu_long_closed"] < 1e-8
        assert abs(imag[2] - stab["nu_short_closed"]) / stab["nu_short_closed"] < 1e-10
        assert any(abs(i - 1.0) < 1e-9 for i in imag)  # vertical pair exactly omega_z=1
        assert stab["routh_stable"]

    def test_vertical_frequency_exactly_one_at_triangular(self):
        assert experiment.vertical_frequency(0.5 - MU_EM, math.sqrt(3) / 2, MU_EM) == 1.0

    def test_routh_criterion_across_threshold_grid(self):
        for mu, expect_stable in ((0.01, True), (0.03, True), (0.05, False), (0.1, False)):
            _, _, _, stable, qa, qb = experiment.triangular_closed_form_rates(mu)
            assert bool(stable) == expect_stable
            if not expect_stable:
                assert qa > 0.0 and qb > 0.0  # complex quartet rates exist

    def test_boundary_degeneracy_at_mu_routh(self):
        stab = experiment.study_stability()
        g = stab["boundary_gamma_minus_quarter_50dps"]
        assert abs(g) < 1e-30  # exact degeneracy in high precision
        case = stab["cases"]["routh_boundary"]
        assert abs(case["routh_product"] - 1.0) < 1e-15

    def test_unstable_growth_rate_matches_linear_sigma(self):
        pert = _pert_study()["L1_unstable_growth"]
        e4 = pert["eps_0.0001"]
        assert e4["rel_error"] < 5e-3
        ratios = pert["bias_scaling_ratios"]
        for rr in ratios:
            assert 50.0 < rr < 200.0  # bias scales ~ proportionally to eps

    def test_LP_frequency_matches_linear_theory(self):
        pert = _pert_study()["L4_longperiod_mode"]
        e4 = pert["amp_0.0001"]
        assert e4["rel_error"] < 1e-6
        assert e4["cj_drift"] < 1e-14


# =========================================================================== #
# L6 -- dimensional/nondimensional cross-check + determinism
# =========================================================================== #
class TestL6UnitsAndDeterminism:
    def test_dimensional_equilibria_match_mapping(self):
        dim = _dim_study()
        for name in ("L1", "L2", "L3"):
            assert dim["equilibria"][name + "_rel_to_L"] < 1e-12

    def test_jacobi_scaling_exact(self):
        scal = _dim_study()["jacobi_scaling_rel"]
        for name, v in scal.items():
            assert v < 1e-15, f"{name}: {v:.3e}"

    def test_trajectory_correspondence_90day(self):
        tr = _dim_study()["trajectory_90day"]
        assert tr["final_pos_rel"] < 1e-12
        assert tr["final_vel_rel"] < 1e-12

    def test_jacobi_along_mapped_path(self):
        assert _dim_study()["jacobi_along_path_worst_rel"] < 1e-12

    def test_determinism_equilibria_rerun_bitwise(self):
        first = experiment.lagrange_points(MU_EM)
        second = experiment.lagrange_points(MU_EM)
        for name in first:
            assert np.array_equal(first[name], second[name])


# =========================================================================== #
# L7 -- adversarial mutants (inline reimplementations; each MUST be caught)
# =========================================================================== #
class TestL7Mutants:
    def test_coriolis_flip_caught_by_trajectory_discriminator(self):
        """Coriolis flip is invisible to spectra AND Jacobi; only this catches it."""
        mut = _mut_study()
        assert mut["clean_relative_residual"] < 1e-10
        assert mut["mutant_signals"]["coriolis_flip"] > 1e-2

    def test_spectra_blind_to_coriolis_flip_documented(self):
        shift = _mut_study()["coriolis_flip_max_eigen_shift"]
        assert shift < 1e-13  # why trajectory-level tests are mandatory

    def test_centrifugal_drop_caught_by_equilibrium_residual(self):
        """Gravity-only gradient cannot balance: residual at true L1 is O(1)."""
        mu = MU_EM
        pts = experiment.lagrange_points(mu)

        def g_grav_only(x, y, z):
            r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
            r2 = math.sqrt((x - 1 + mu) ** 2 + y * y + z * z)
            return (
                -(1 - mu) * (x + mu) / r1**3 - mu * (x - 1 + mu) / r2**3,
                -(1 - mu) * y / r1**3 - mu * y / r2**3,
                -(1 - mu) * z / r1**3 - mu * z / r2**3,
            )

        gx, gy, gz = g_grav_only(pts["L1"][0], 0.0, 0.0)
        assert abs(gx) > 0.1  # vs production gate 2e-14

    def test_mapping_sign_flip_caught(self):
        assert _mut_study()["mutant_signals"]["mapping_sign_flip"] > 1e-2

    def test_mu_convention_caught_by_firewall_grade_signal(self):
        sig = _mut_study()["mutant_signals"]["mu_convention_gamma_rel_shift"]
        assert sig > 1e-12  # firewall bound (percent-level anchors CANNOT catch it)

    def test_body_swap_domain_rejected(self):
        assert _mut_study()["mutant_signals"]["body_swap_domain_rejected"] == 1.0

    def test_all_registered_mutants_caught(self):
        assert _mut_study()["all_mutants_caught"]

    def test_jacobi_sign_flip_detected_by_inline_bridge(self):
        """Flipping the C convention breaks the bridge identity by O(1)."""
        mu = MU_EM
        pts = experiment.lagrange_points(mu)
        s = np.concatenate([pts["L4"], [0.05, -0.02, 0.01]])
        c_flipped = -jacobi_inline(s, mu)  # mutant convention v^2 - 2*Omega
        si = experiment.rot_to_inert(s, 0.7)
        energy, hz = experiment.inertial_state_quantities(si, mu, 0.7)
        assert abs(c_flipped - 2.0 * (hz - energy)) > 1.0

    def test_root_branch_swap_detected_by_containment(self):
        """Swapping L1/L2 labels violates the proven bracket containment."""
        r = experiment.solve_collinear(MU_EM)
        assert not (0.5 - MU_EM <= r["L2"] <= 1.0 - MU_EM)  # true L2 is not in L1's bracket
        assert not (1.0 - MU_EM <= r["L1"] <= 2.0 - MU_EM)
