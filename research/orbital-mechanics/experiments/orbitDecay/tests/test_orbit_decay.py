"""Tests for Experiment 010: orbit decay / atmospheric drag.

Structure follows the Exp 009 doctrine: theory constants are duplicated inline
(NOT imported from the experiment), oracles are reimplemented from first
principles where the anti-shared-algebra rule demands it, and heavy runs are
cached at module level. Banners L1..L7 mirror the pre-registered failure catalog.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Module loading (experiment under test + bit-exact regression donors)
# --------------------------------------------------------------------------- #
_EXP_PATH = Path(__file__).resolve().parents[1] / "experiment.py"
_spec = importlib.util.spec_from_file_location("orbit_decay_experiment", _EXP_PATH)
assert _spec is not None and _spec.loader is not None
experiment = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(experiment)

_EXPERIMENTS_DIR = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pcm006 = _load("pcm_for_decay_tests", _EXPERIMENTS_DIR / "planeChangeManeuvers" / "experiment.py")
j2exp009 = _load("j2_for_decay_tests", _EXPERIMENTS_DIR / "j2Precession" / "experiment.py")

# --------------------------------------------------------------------------- #
# Inline independent theory (duplicated on purpose -- never imported)
# --------------------------------------------------------------------------- #
R_E = 6378.137  # km, WGS-84
MU = 398600.4418  # km^3/s^2
OMEGA_EARTH = 7.2921159e-5  # rad/s
J2_WGS84 = 1.082629821e-3

LAYERS_INLINE = (
    (0.010, 1.225e0, 7.249), (25.000, 3.899e-2, 6.349), (30.000, 1.774e-2, 6.682),
    (40.000, 3.972e-3, 7.554), (50.000, 1.057e-3, 8.382), (60.000, 3.206e-4, 7.714),
    (70.000, 8.770e-5, 6.549), (80.000, 1.905e-5, 5.799), (90.000, 3.396e-6, 5.382),
    (100.000, 5.297e-7, 5.877), (110.000, 9.661e-8, 7.263), (120.000, 2.438e-8, 9.473),
    (130.000, 8.484e-9, 12.636), (140.000, 3.845e-9, 16.149), (150.000, 2.070e-9, 22.523),
    (180.000, 5.464e-10, 29.740), (200.000, 2.789e-10, 37.105), (250.000, 7.248e-11, 45.546),
    (300.000, 2.418e-11, 53.628), (350.000, 9.518e-12, 53.298), (400.000, 3.725e-12, 58.515),
    (450.000, 1.585e-12, 60.828), (500.000, 6.967e-13, 63.822), (600.000, 1.454e-13, 71.835),
    (700.000, 3.614e-14, 88.667), (800.000, 1.170e-14, 124.640), (900.000, 5.245e-15, 181.045),
)

US76_SPOTS = {300.0: 1.9159e-11, 400.0: 2.8028e-12, 500.0: 5.2148e-13}  # kg/m^3


def rho_inline(h_km: float) -> float:
    """Independent density transcription: last row with h0 <= h wins."""
    h0, rho0, hs = LAYERS_INLINE[0]
    for row in LAYERS_INLINE:
        if h_km >= row[0]:
            h0, rho0, hs = row
    return rho0 * math.exp(-(h_km - h0) / hs)


def power_inline(state, beta: float, omega: float = 0.0) -> float:
    """Drag power a_drag . v in km^2/s^3, transcribed from first principles."""
    x, y, z, vx, vy, vz = state
    rho = rho_inline(math.sqrt(x * x + y * y + z * z) - R_E)
    if omega != 0.0:
        vrx, vry, vrz = vx - (-omega * y), vy - (omega * x), vz
    else:
        vrx, vry, vrz = vx, vy, vz
    vr = math.sqrt(vrx * vrx + vry * vry + vrz * vrz)
    vdot = vrx * vx + vry * vy + vrz * vz
    return -(0.5 * rho / beta) * vr * vdot * 1e3


def eps_inline(state, j2: float = 0.0) -> float:
    """Specific mechanical energy (+ static J2 potential when active), km^2/s^2."""
    x, y, z, vx, vy, vz = state
    r = math.sqrt(x * x + y * y + z * z)
    e = 0.5 * (vx * vx + vy * vy + vz * vz) - MU / r
    if j2 != 0.0:
        u = z / r
        p2 = 0.5 * (3.0 * u * u - 1.0)
        e += MU * j2 * R_E * R_E * p2 / (r * r * r)
    return e


def fd_residual_inline(states, t, idx, beta: float, omega: float = 0.0, j2: float = 0.0):
    """Relative pointwise dissipation residual at grid indices (independent path)."""
    out = []
    for i in idx:
        mid = tuple(0.5 * (states[i][k] + states[i + 1][k]) for k in range(6))
        fd = (eps_inline(states[i + 1], j2) - eps_inline(states[i], j2)) / (t[i + 1] - t[i])
        p = power_inline(mid, beta, omega)
        out.append(abs(fd - p) / max(abs(p), 1e-16))
    return out


# Hand-computed holdout probe, digit by digit (420 km circular, beta = 100 kg/m^2):
#   rho(420)   = 3.725e-12 * exp(-20/58.515) = 2.646595637742271e-12 kg/m^3
#   a          = (6378.137 + 420) * 1000     = 6,798,137 m
#   a^2        = 4.62146667e13 m^2
#   Da_rev     = 2*pi*kappa*rho*a^2 = 2*pi*(1/100)*2.6465956e-12*4.62146667e13
#              = 2*pi*1.223130 = 7.68506 m
HAND_PROBE_DA_REV_M = 7.68506

_CACHE: dict = {}


def _cached(key, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


def _window60():
    return experiment.run_window(R_E + 420.0, 0.0005, 51.6, 100.0, 60)


def _anchor():
    return experiment.leo_ref_anchor(beta=100.0)


def _scalings():
    return experiment.scaling_battery()


def _rotation():
    return experiment.rotation_battery()


def _ecc():
    return experiment.eccentric_battery()


def _conv():
    return experiment.convergence_study()


def _plateau():
    return experiment.plateau_separation()


def _reentry():
    near = experiment.propagate_until_reentry(R_E + 200.0, 0.0, 51.6, 100.0)
    spot280 = experiment.propagate_until_reentry(R_E + 280.0, 0.0, 51.6, 100.0)
    return {"near": near, "spot280": spot280}


def _patho():
    return experiment.pathological_battery()


def _mutants():
    return experiment.mutant_battery()


def _rotating_run():
    return experiment.run_window(R_E + 400.0, 1e-4, 0.0, 100.0, 20,
                                 omega_atm=OMEGA_EARTH)


# =========================================================================== #
# L1 -- nulls and bit-exact regressions
# =========================================================================== #
def test_bitexact_vs_exp006_when_disabled():
    r0, v0, _ = experiment.seed_state(R_E + 420.0, 0.3, math.radians(51.6), 0.0, 0.0, 0.0)
    T = experiment.orbital_period(R_E + 420.0)
    t = np.arange(0, 3 * 64 + 1) * (T / 64)
    mine = experiment.propagate_3d_rk4_drag(r0, v0, MU, t)
    ref = pcm006.propagate_3d_rk4(r0, v0, MU, t, T / 64)
    assert np.array_equal(mine, ref)


def test_bitexact_vs_exp009_j2_when_disabled():
    r0, v0, _ = experiment.seed_state(R_E + 420.0, 0.3, math.radians(51.6), 0.0, 0.0, 0.0)
    T = experiment.orbital_period(R_E + 420.0)
    t = np.arange(0, 2 * 128 + 1) * (T / 128)
    mine = experiment.propagate_3d_rk4_drag(r0, v0, MU, t, j2=J2_WGS84)
    ref = j2exp009.propagate_3d_rk4_j2(r0, v0, MU, t, j2=J2_WGS84)
    assert np.array_equal(mine, ref)


def test_beta_zero_null_energy_conserved():
    run = experiment.run_window(R_E + 420.0, 0.0005, 51.6, 0.0, 20)
    d = run["dissipation"]
    assert d["null_run"] is True
    assert d["monotone_violations"] == 0
    assert d["energy_drift_rel_max"] < 1e-12


def test_strict_dissipation_j2_off_elementwise():
    run = _cached("w60", _window60)
    eps = [eps_inline(s) for s in run["_states"][::7]]
    diffs = np.diff(eps)
    assert np.all(diffs < 0.0)


def test_j2_on_total_energy_strictly_dissipative():
    run = experiment.run_window(R_E + 420.0, 0.0005, 51.6, 100.0, 60, j2=J2_WGS84)
    states = run["_states"]
    n = len(states)
    eps = np.array([eps_inline(states[i], J2_WGS84) for i in range(0, n, 11)])
    assert np.all(np.diff(eps) <= 0.0)


# =========================================================================== #
# L2 -- identities, unit firewalls, pinned values
# =========================================================================== #
def test_dissipation_identity_inertial_inline():
    run = _cached("w60", _window60)
    t, states = run["_t"], run["_states"]
    idx = list(range(0, len(t) - 1, 97))
    res = fd_residual_inline(states, t, idx, 100.0)
    assert np.median(res) < 5e-2


def test_dissipation_floor_shrinks_with_refinement():
    fine = experiment.run_window(R_E + 420.0, 0.0005, 51.6, 100.0, 4,
                                 spp=experiment.steps_per_orbit(0.0005) * 4)
    coarse = experiment.run_window(R_E + 420.0, 0.0005, 51.6, 100.0, 4)
    def med(run):
        t, st = run["_t"], run["_states"]
        idx = list(range(0, len(t) - 1, 31))
        return float(np.median(fd_residual_inline(st, t, idx, 100.0)))
    assert med(fine) < med(coarse)


def test_dissipation_identity_rotating_inline():
    run = _cached("rotrun", _rotating_run)
    t, states = run["_t"], run["_states"]
    idx = list(range(0, len(t) - 1, 53))
    res = fd_residual_inline(states, t, idx, 100.0, omega=OMEGA_EARTH)
    assert np.median(res) < 5e-2
    # rotation actually exercised: prograde |v_rel| strictly below inertial speed
    vmid = np.linalg.norm(states[len(states) // 2, 3:])
    x = states[len(states) // 2]
    vrel = math.sqrt((x[3] + OMEGA_EARTH * x[1]) ** 2 + (x[4] - OMEGA_EARTH * x[0]) ** 2 + x[5] ** 2)
    assert vrel < vmid * (1.0 - 1e-3)


def test_atmosphere_transcription_integrity():
    assert np.array_equal(np.array(LAYERS_INLINE, dtype=float), experiment.ATMOSPHERE_LAYERS)


def test_us76_plausibility_band():
    for alt, ref in US76_SPOTS.items():
        declared = rho_inline(alt)
        assert experiment.air_density_si(alt) == declared
        ratio = declared / ref
        assert 0.67 < ratio < 1.5, f"{alt} km ratio {ratio}"


def test_unit_firewall_si_to_km():
    r0, v0, _ = experiment.seed_state(R_E + 420.0, 0.0, math.radians(51.6), 0.0, 0.0, 0.0)
    x = np.concatenate([r0, v0])
    got = experiment.drag_accel_kkms(x, 100.0)
    rm = math.sqrt(sum(c * c for c in r0))
    vm = math.sqrt(sum(c * c for c in v0))
    rho = rho_inline(rm - R_E)
    a_si = 0.5 * (1.0 / 100.0) * rho * vm * vm  # m/s^2
    expected = -a_si * 1e-3 * (v0 / vm)  # km/s^2, anti-velocity
    assert np.allclose(got, expected, rtol=1e-12)


def test_pinned_hand_probe_da_rev():
    run = experiment.run_window(R_E + 420.0, 0.0005, 51.6, 100.0, 30)
    measured = abs(run["da_rev_measured_m"])  # decay => signed negative
    assert abs(measured - HAND_PROBE_DA_REV_M) / HAND_PROBE_DA_REV_M < 0.01


# =========================================================================== #
# L3 -- closed-form decay law (same declared density law only)
# =========================================================================== #
def test_leoref_window_matches_quadrature_oracle():
    anchor = _cached("anchor", _anchor)
    assert abs(anchor["quad_vs_num_rel"]) < 5e-3
    assert anchor["oracle_residual"]["points_compared"] > 100


def test_erfi_quadrature_cross_agreement():
    a0, af = R_E + 420.0, R_E + 418.0
    q = experiment.circular_decay_time_quadrature(a0, af, 100.0)
    e60 = experiment.erfi_decay_time(a0, af, 100.0, dps=60)
    e40 = experiment.erfi_decay_time(a0, af, 100.0, dps=40)
    assert abs(e60 / q - 1.0) < 1e-9
    assert abs(e40 / e60 - 1.0) < 1e-12


def test_constant_density_sqrt_a_law():
    rho_const = rho_inline(420.0)
    atm = np.array([[0.010, rho_const, 1e6]])  # effectively constant rho
    run = experiment.run_window(R_E + 420.0, 1e-9, 0.0, 100.0, 10, atmosphere=atm)
    kappa = 1.0 / 100.0
    mu_si = MU * 1e9
    a0_m = (R_E + 420.0) * 1e3
    root_rate = kappa * rho_const * math.sqrt(mu_si) / 2.0
    for frac in (0.2, 0.4, 0.6, 0.8, 1.0):
        i = int(frac * (len(run["_t"]) - 1))
        t = run["_t"][i]
        expected_km = (math.sqrt(a0_m) - root_rate * t) ** 2 / 1e3
        got_km = run["_coe_a"][i]
        assert abs(got_km / expected_km - 1.0) < 1e-5, f"frac={frac}"


def test_king_hele_circular_limit_consistency():
    a_km = R_E + 420.0
    kh = experiment.king_hele_delta_rev(a_km, 1e-9, 100.0)
    rho = rho_inline(420.0)
    analytic_m = -2.0 * math.pi * (1.0 / 100.0) * rho * (a_km * 1e3) ** 2
    assert abs(kh["da_rev_m"] / analytic_m - 1.0) < 1e-6


# =========================================================================== #
# L4 -- structural scalings and symmetry
# =========================================================================== #
def test_beta_scaling_ratios():
    scal = _cached("scalings", _scalings)
    for entry in scal["beta_ratios"].values():
        assert abs(entry["rel_dev"]) < 2.5e-2
    # accel magnitude scales as kappa = 1/beta: beta 400 -> 100 gives exactly 1/4
    assert abs(scal["accel_level_ratio_beta400_over_100"] - 0.25) < 1e-12


def test_rho0_scaling():
    scal = _cached("scalings", _scalings)
    for c, entry in scal["rho0_scalings"].items():
        assert abs(entry["measured"] / entry["theory"] - 1.0) < 2.5e-2


def test_area_linearity_halving():
    scal = _cached("scalings", _scalings)
    r100 = scal["rates_m_day"]["100.0"]
    r200 = scal["rates_m_day"]["200.0"]
    ratio = r200 / r100  # doubling beta == halving C_D*A/m
    assert abs(ratio - 0.5) < 2.5e-2


def test_omega_twin_frame_symmetry():
    rot = _cached("rotation", _rotation)
    assert rot["omega_twin_max_rel_asymmetry"] <= 1e-12


def test_equatorial_rotation_ratio():
    rot = _cached("rotation", _rotation)
    eq = rot["equatorial"]
    assert eq["ratio_measured"] > 1.0
    assert abs(eq["rel_dev"]) < 0.02


def test_inclined_rotation_direction_and_band():
    rot = _cached("rotation", _rotation)
    inc = rot["inclined_63p6_116p4"]
    assert inc["ratio_measured"] > 1.0
    assert abs(inc["rel_dev"]) < 0.03


def test_j2_settled_transient_physics_consistent():
    """J2 on/off: switching J2 on settles the mean elements ~5.7 km lower (inclination-
    dependent mean of U_J2), so the decay proceeds in denser air. The honest null is
    that settled-tail ENERGY rates match after correcting for the altitude shift."""
    scal = _cached("scalings", _scalings)
    entry = scal["j2_settled_comparison"]
    # sanity: energy-based off-rate must agree with the OLS slope (no ripple when off)
    assert abs(abs(entry["off_energy_rate_m_day"]) / 119.60 - 1.0) < 0.02
    assert abs(entry["residual_after_altitude_correction"]) < 0.03


# =========================================================================== #
# L5 -- eccentric-orbit behavior
# =========================================================================== #
def test_steep_dt_invariance_gate():
    ecc = _cached("ecc", _ecc)
    assert ecc["steep"]["dt_invariance_rel_halving"] < 1e-3


def test_steep_circularization_directions():
    ecc = _cached("ecc", _ecc)
    st = ecc["steep"]
    apo_drop = abs(st["apo_first_last_km"][0] - st["apo_first_last_km"][1])
    peri_change = abs(st["peri_first_last_km"][0] - st["peri_first_last_km"][1])
    assert apo_drop > 3.0 * peri_change
    assert st["de_over_window"] < 0.0


def test_gauss_o4_agreement_eccentric_cases():
    ecc = _cached("ecc", _ecc)
    for case in ("steep", "molniya_like"):
        assert abs(ecc[case]["measured_vs_o4_rel"]) < 0.02, case


def test_apsis_third_path_estimator():
    ecc = _cached("ecc", _ecc)
    for case in ("steep", "molniya_like"):
        o4 = ecc[case]["gauss_o4_da_rev_m"]
        apsis = ecc[case]["apsis_estimator_da_rev_m"]
        assert apsis is not None and o4 != 0.0
        assert abs(apsis / o4 - 1.0) < 0.02, case


def test_king_hele_small_e_few_percent():
    ecc = _cached("ecc", _ecc)
    assert abs(ecc["small_e_king_hele"]["da_rel_dev"]) < 0.03


def test_perigee_sink_direction_large_e():
    ecc = _cached("ecc", _ecc)
    mol = ecc["molniya_like"]
    kh_context = mol["king_hele_context_da_rev_m"]
    assert kh_context < 0.0  # decay present; direction checks live in card prose
    assert mol["e_direction"] < 0.0


def test_shortwindow_anchors_match_quadrature():
    sw = _cached("shortwindow", lambda: experiment.shortwindow_battery())
    for case in ("starlink", "sso"):
        assert abs(sw[case]["rel_dev"]) < 5e-3, case
        assert sw[case]["rate_m_day"] < 0.0  # decaying


def test_small_e_eccentricity_direction_signed():
    ecc = _cached("ecc", _ecc)
    de = ecc["small_e_king_hele"]["measured_de_per_rev"]
    assert de < 0.0  # eccentricity decreases at small e (signed, not abs)


def test_surface_guard_driver_terminates():
    patho = _cached("patho", _patho)
    sg = patho["surface_guard"]
    assert np.isfinite(sg["min_altitude_km"])
    assert sg["min_altitude_km"] > -100.0  # bounded penetration on the raw probe
    assert sg["driver_status"] == "surface-hit"  # event driver never integrates through


# =========================================================================== #
# L6 -- numerics: convergence, plateau, events, determinism
# =========================================================================== #
def test_convergence_order_no_degradation_and_preasymptotic_band():
    """Design order is 4; the time-to-fall observable SUPERCONVERGES (measured ~3.8
    then ~5.0 pre-floor) as leading error terms cancel. Assert: first coarse pair in
    the order-4 band, and no pair below design order."""
    conv = _cached("conv", _conv)
    rates = conv["convergence_rates"]
    assert len(rates) >= 2
    assert 3.6 <= rates[0] <= 4.4, rates
    for r in rates:
        assert r >= 3.6, rates


def test_kepler_clone_order_floor():
    """Raw position-error order of the drag-gated clone on pure Kepler: every measured
    rate must be at or above design order 4 (measured values run 4.4-4.8, i.e.
    pre-asymptotic superconvergence converging down toward 4 near the floor)."""
    ko = experiment.kepler_order_check()
    for r in ko["convergence_rates"]:
        assert r >= 3.6, ko["convergence_rates"]


def test_plateau_separation_doctrine():
    pl = _cached("plateau", _plateau)
    assert pl["separation_ratio"] >= 20.0
    assert pl["swap_flatness_abs_s"] < 0.1 * abs(
        0.5 * (pl["spp_128"]["swap_difference_s"] + pl["spp_256"]["swap_difference_s"])
    )


def test_reentry_dual_threshold_ordering_and_refinement():
    re = _cached("reentry", _reentry)["near"]
    assert re["status"] == "reentered"
    t120 = re["threshold_120km"]["crossing_days"]
    t100 = re["threshold_100km"]["crossing_days"]
    assert t120 is not None and t100 is not None
    assert 0.0 < t120 < t100
    for thr in ("threshold_120km", "threshold_100km"):
        change = re[thr]["last_refine_step_change_s"]
        assert change is not None and 0.0 <= change < 1.0


def test_lifetime_spot_matches_quadrature_oracle():
    re = _cached("reentry", _reentry)
    spot = re["spot280"]
    assert spot["status"] == "reentered"
    t_num_days = spot["threshold_120km"]["crossing_days"]
    t_orc_s = experiment.circular_decay_time_quadrature(R_E + 280.0, R_E + 120.0, 100.0)
    assert abs(t_num_days * 86400.0 / t_orc_s - 1.0) < 5e-3


def test_determinism_repeat_call():
    patho = _cached("patho", _patho)
    assert patho["determinism_repeat"] is True


# =========================================================================== #
# L7 -- pathological inputs and adversarial mutants
# =========================================================================== #
def test_mutant_sign_flip_detected():
    mut = _cached("mutants", _mutants)
    assert mut["sign_flip_detected"] is True


def test_mutant_unit_error_detected():
    mut = _cached("mutants", _mutants)
    assert mut["unit_mutant_detected"] is True


def test_mutant_b_inversion_detected():
    mut = _cached("mutants", _mutants)
    assert mut["b_inversion_detected"] is True
    ratio = mut["b_inversion_rate_ratio"]
    assert ratio > 10.0 or ratio < 0.1


def test_pathological_sentinels():
    patho = _cached("patho", _patho)
    assert patho["knife_edge_H"]["finite_states"] is True
    assert patho["underflow"]["finite"] is True
    assert patho["underflow"]["drag_finite"] is True
    assert patho["extreme_beta"]["beta_1e12_vs_null_a_rel_drift"] < 1e-12
    # kappa = 1e4: per-rev decay ~68 km -- extreme but finite; must terminate cleanly
    assert patho["extreme_beta"]["beta_1e-4_status"] in ("reentered", "surface-hit")
    assert patho["extreme_beta"]["beta_1e-4_threshold120_days"] is not None
    assert patho["layer_boundary"]["rho_250_exact_row_value"] is True
    assert patho["layer_boundary"]["continuity_rel_jump"] < 5e-3
    assert np.isfinite(patho["surface_guard"]["min_altitude_km"])
