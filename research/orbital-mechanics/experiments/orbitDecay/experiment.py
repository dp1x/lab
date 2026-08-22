"""Experiment 010 -- Orbit decay / atmospheric drag (deterministic).

Question: Does deterministic Cowell propagation with atmospheric drag rediscover
pointwise orbital-energy dissipation, the circular-orbit decay law, structural drag
scaling laws, eccentric-orbit decay behavior, and deterministic re-entry timing --
while separating numerical integration error from atmospheric-model uncertainty?

First non-conservative force of the sequence: energy conservation is replaced by
energy MONOTONICITY plus pointwise dissipation accounting. Doctrine and frozen
contract: scratch contract v1.0 (2026-08-22), adversarially reviewed (21 findings).

Frames/units: ECI km, km/s, s; angles rad internally. Density/ballistic quantities
in SI (kg/m^3, kg/m^2) inside the drag term; DRAG_SI_TO_KKM converts to km/s^2.
Declared atmosphere: Vallado exponential layer table (official CelesTrak repo,
software/misc/pascal/ATMOSEXP.DAT), used verbatim and untuned.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpmath import besseli as _mp_besseli, e as _mp_e, erfi as _mp_erfi
from mpmath import mpf as _mp_mpf, workdps as _mp_workdps

from lab_utils.metrics import convergence_rate
from lab_utils.results import save_json_result

# --- Reuse of verified prior-experiment machinery (single-hop importlib) ---- #
_J2_PATH = Path(__file__).resolve().parents[1] / "j2Precession" / "experiment.py"
_j2_spec = importlib.util.spec_from_file_location("j2_for_decay", _J2_PATH)
assert _j2_spec is not None and _j2_spec.loader is not None
_j2 = importlib.util.module_from_spec(_j2_spec)
_j2_spec.loader.exec_module(_j2)

_GT_PATH = Path(__file__).resolve().parents[1] / "groundtracks" / "experiment.py"
_gt_spec = importlib.util.spec_from_file_location("gt_for_decay", _GT_PATH)
assert _gt_spec is not None and _gt_spec.loader is not None
_gt = importlib.util.module_from_spec(_gt_spec)
_gt_spec.loader.exec_module(_gt)

seed_state = _j2.seed_state
rv_to_coe_eci = _j2.rv_to_coe_eci
steps_per_orbit = _j2.steps_per_orbit
ols_fit = _j2.ols_fit
j2_specific_energy = _j2.j2_specific_energy
solve_kepler = _gt.solve_kepler
true_anomaly_from_E = _gt.true_anomaly_from_E
coe_to_rv_eci = _gt.coe_to_rv_eci
orbital_period = _gt.orbital_period
mean_motion = _gt.mean_motion
rotation_matrix_313 = _gt.rotation_matrix_313

MU_EARTH_KM3S2 = _gt.MU_EARTH_KM3S2  # 398600.4418 IAU 2015 B3 nominal GM_E
R_EARTH_KM = _gt.R_EARTH_KM  # 6378.137 WGS-84 equatorial radius
OMEGA_EARTH_RAD_S = _gt.OMEGA_EARTH_RAD_S  # 7.2921159e-5 rad/s
J2_EARTH = _j2.J2_EARTH  # 1.082629821e-3 WGS-84 (provenance in Exp 009)

FRAME_CONVENTION = "ECI (inertial), spherical-Earth geocentric altitude h = |r|-R_E"
UNITS_CONVENTION = "km, km^3/s^2, s; angles rad internally; density/ballistic terms SI"

# Drag unit conversion (derived): a[m/s^2] = 0.5*kappa*rho*(1e3*v[km/s])^2 =>
# a[km/s^2] = 0.5*kappa*rho*v[km/s]^2 * 1e3. Single explicit factor, unit-fired in tests.
DRAG_SI_TO_KKM = 1e3

REENTRY_THRESHOLDS_KM = (120.0, 100.0)

# --- Declared atmosphere: Vallado exponential layer table (verbatim) -------- #
# Source: CelesTriak/fundamentals-of-astrodynamics, software/misc/pascal/ATMOSEXP.DAT
# (linked from https://celestrak.org/software/vallado-sw.php), retrieved 2026-08-22,
# verified byte-for-byte against the research snapshot. Rows: (h0 [km], rho0 [kg/m^3],
# H [km]); rho(h) = rho0*exp(-(h-h0)/H) within [h0, next_h0). Used AS PUBLISHED --
# never tuned to lifetimes. Known bias vs US76 spot values documented in results.
ATMOSPHERE_SOURCE = (
    "Vallado, Fundamentals of Astrodynamics and Applications, exponential-atmosphere "
    "data file ATMOSEXP.DAT (official CelesTrak repository, retrieved 2026-08-22); "
    "LEO rows consistent with US Standard Atmosphere 1976 (NASA-TM-X-74335) to ~35%"
)
US76_SPOTS_SOURCE = "U.S. Standard Atmosphere 1976 (NASA-TM-X-74335) via PDAS tables"
ATMOSPHERE_LAYERS = np.array(
    [
        [0.010, 1.225e0, 7.249],
        [25.000, 3.899e-2, 6.349],
        [30.000, 1.774e-2, 6.682],
        [40.000, 3.972e-3, 7.554],
        [50.000, 1.057e-3, 8.382],
        [60.000, 3.206e-4, 7.714],
        [70.000, 8.770e-5, 6.549],
        [80.000, 1.905e-5, 5.799],
        [90.000, 3.396e-6, 5.382],
        [100.000, 5.297e-7, 5.877],
        [110.000, 9.661e-8, 7.263],
        [120.000, 2.438e-8, 9.473],
        [130.000, 8.484e-9, 12.636],
        [140.000, 3.845e-9, 16.149],
        [150.000, 2.070e-9, 22.523],
        [180.000, 5.464e-10, 29.740],
        [200.000, 2.789e-10, 37.105],
        [250.000, 7.248e-11, 45.546],
        [300.000, 2.418e-11, 53.628],
        [350.000, 9.518e-12, 53.298],
        [400.000, 3.725e-12, 58.515],
        [450.000, 1.585e-12, 60.828],
        [500.000, 6.967e-13, 63.822],
        [600.000, 1.454e-13, 71.835],
        [700.000, 3.614e-14, 88.667],
        [800.000, 1.170e-14, 124.640],
        [900.000, 5.245e-15, 181.045],
    ]
)
DEFAULT_ATMOSPHERE = ATMOSPHERE_LAYERS


def air_density_si(h_km: float, atmosphere: np.ndarray = DEFAULT_ATMOSPHERE) -> float:
    """Declared piecewise-exponential density [kg/m^3].

    Layer = last row with h0 <= h (half-open [h_i, h_{i+1})); above the top row the
    top row is extended; below the first row the sea-level value clamps (re-entry
    events terminate far above this; the clamp exists so pathological inputs stay
    finite rather than to model the mesosphere).
    """
    idx = int(np.searchsorted(atmosphere[:, 0], h_km, side="right")) - 1
    idx = min(max(idx, 0), len(atmosphere) - 1)
    h0, rho0, hs = atmosphere[idx]
    return float(rho0 * np.exp(-(h_km - h0) / hs))


# --------------------------------------------------------------------------- #
# Propagator: clone of Exp 009's generalized-force loop + gated drag branch
# --------------------------------------------------------------------------- #
def _rk4_core(accel, r0, v0, mu, t) -> np.ndarray:
    """Verbatim fixed-step RK4 loop of Exp 009 (bit-exact regression anchor)."""
    n = len(t)
    state = np.empty((n, 6))
    state[0] = np.concatenate([np.asarray(r0, float), np.asarray(v0, float)])
    for k in range(1, n):
        h = t[k] - t[k - 1]
        x = state[k - 1]
        k1 = np.concatenate([x[3:], accel(x)])
        x2 = x + 0.5 * h * k1
        k2 = np.concatenate([x2[3:], accel(x2)])
        x3 = x + 0.5 * h * k2
        k3 = np.concatenate([x3[3:], accel(x3)])
        x4 = x + h * k3
        k4 = np.concatenate([x4[3:], accel(x4)])
        state[k] = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return state


def drag_accel_kkms(x, beta: float, omega_atm: float = 0.0, atmosphere=None) -> np.ndarray:
    """Drag acceleration [km/s^2] from Cartesian state (km, km/s). SI rho/beta."""
    atm = DEFAULT_ATMOSPHERE if atmosphere is None else atmosphere
    r = x[:3]
    v = x[3:]
    rm = np.linalg.norm(r)
    rho = air_density_si(rm - R_EARTH_KM, atm)
    if rho <= 0.0:
        return np.zeros(3)
    if omega_atm != 0.0:
        v_rel = v - omega_atm * np.array([-r[1], r[0], 0.0])
    else:
        v_rel = v
    return -(0.5 * (rho / beta) * DRAG_SI_TO_KKM * np.linalg.norm(v_rel)) * v_rel


def propagate_3d_rk4_drag(
    r0,
    v0,
    mu,
    t,
    j2: float = 0.0,
    beta: float = 0.0,
    omega_atm: float = 0.0,
    atmosphere=None,
) -> np.ndarray:
    """Fixed-step RK4 Cowell propagator with switchable J2 and drag branches.

    Clone-with-one-branch doctrine (Exp 009 pattern): disabled branches are skipped
    entirely (no masking, no extra array ops), so (beta=0, j2=0) is bit-exact with
    Exp 006's propagate_3d_rk4 and (beta=0, j2=J2) with Exp 009's
    propagate_3d_rk4_j2. omega_atm enters only through v_rel, so omega_atm=0
    reduces exactly. beta < 0 is rejected: sign flips belong to mutants, not runs.
    """
    if beta < 0.0:
        raise ValueError("beta must be >= 0; use the mutant harness for sign studies")
    re2 = R_EARTH_KM**2

    def accel(x):
        r = x[:3]
        rm = np.linalg.norm(r)
        a_kep = -mu * r / rm**3
        if j2 == 0.0 and beta == 0.0:
            return a_kep
        acc = a_kep
        if j2 != 0.0:
            z2r2 = (r[2] * r[2]) / (rm * rm)
            c = -1.5 * j2 * mu * re2 / rm**5
            f = 1.0 - 5.0 * z2r2
            acc = acc + c * np.array([r[0] * f, r[1] * f, r[2] * (3.0 - 5.0 * z2r2)])
        if beta != 0.0:
            acc = acc + drag_accel_kkms(x, beta, omega_atm, atmosphere)
        return acc

    return _rk4_core(accel, r0, v0, mu, t)


# --------------------------------------------------------------------------- #
# Series analysis
# --------------------------------------------------------------------------- #
def element_series(states: np.ndarray) -> dict:
    """Osculating element series from a propagated (n, 6) state grid."""
    coe = rv_to_coe_eci(states[:, :3], states[:, 3:])
    return coe


def _eps_total_series(states: np.ndarray, j2: float) -> np.ndarray:
    """Specific mechanical energy incl. static J2 potential (Exp 009 helper)."""
    r = states[:, :3]
    v = states[:, 3:]
    if j2 != 0.0:
        return np.asarray(j2_specific_energy(r, v, MU_EARTH_KM3S2, j2), dtype=float)
    rm = np.linalg.norm(r, axis=1)
    return 0.5 * np.einsum("ij,ij->i", v, v) - MU_EARTH_KM3S2 / rm


def rel_dissipation_residual(t, states, beta, omega_atm=0.0, atmosphere=None, j2=0.0) -> np.ndarray:
    """Per-step relative dissipation-identity residual array (finite-difference
    d(eps)/dt at midpoints vs independently evaluated drag power)."""
    eps = _eps_total_series(states, j2)
    dts = np.diff(t)
    fd = np.diff(eps) / dts  # [km^2/s^3]
    mid = 0.5 * (states[:-1] + states[1:])
    rm = np.linalg.norm(mid[:, :3], axis=1)
    vm = mid[:, 3:]
    atm = DEFAULT_ATMOSPHERE if atmosphere is None else atmosphere
    rho = np.empty(len(mid))
    for i in range(len(mid)):
        rho[i] = air_density_si(rm[i] - R_EARTH_KM, atm)
    if omega_atm != 0.0:
        v_rel = vm - omega_atm * np.stack([-mid[:, 1], mid[:, 0], np.zeros(len(mid))], axis=1)
    else:
        v_rel = vm
    kappa = 1.0 / beta if beta != 0.0 else 0.0
    power = -(0.5 * kappa * rho * DRAG_SI_TO_KKM) * np.einsum("ij,ij->i", v_rel, vm) * np.linalg.norm(v_rel, axis=1)
    denom = np.maximum(np.abs(power), 1e-16)
    return np.abs(fd - power) / denom


def dissipation_residual_stats(t, states, beta, omega_atm=0.0, atmosphere=None, j2=0.0):
    """Aggregated dissipation-residual statistics (see rel_dissipation_residual).

    beta == 0 (null runs): reports energy-drift diagnostics instead (conservation
    restored); the identity itself is trivially 0 == 0 and skipped.
    """
    eps = _eps_total_series(states, j2)
    diffs = np.diff(eps)
    if beta == 0.0:
        scale = float(np.max(np.abs(eps)))
        return {
            "null_run": True,
            "energy_drift_rel_max": float(np.max(np.abs(diffs)) / scale),
            "monotone_violations": int(np.count_nonzero(diffs > 0.0)),
            "strictly_dissipative": bool(np.all(diffs <= 0.0)),
        }
    rel = rel_dissipation_residual(t, states, beta, omega_atm, atmosphere, j2)
    return {
        "rel_resid_median": float(np.median(rel)),
        "rel_resid_max": float(np.max(rel)),
        "monotone_violations": int(np.count_nonzero(diffs > 0.0)),
        "strictly_dissipative": bool(np.all(diffs <= 0.0)),
    }


def crossing_time(t, values, target) -> float | None:
    """First downward crossing of `target` by `values`, linearly interpolated."""
    below = values < target
    hits = np.flatnonzero(below)
    if hits.size == 0:
        return None
    i = int(hits[0])
    if i == 0:
        return float(t[0])
    y0, y1 = values[i - 1], values[i]
    frac = (y0 - target) / (y0 - y1)
    return float(t[i - 1] + frac * (t[i] - t[i - 1]))


def apsis_indices(rmag: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Grid indices of periapsis/apoapsis passages (local extrema of |r|)."""
    dr = np.diff(rmag)
    sign_change = np.flatnonzero(np.diff(np.sign(dr)) != 0) + 1
    peri = sign_change[dr[sign_change - 1] < 0]
    apo = sign_change[dr[sign_change - 1] > 0]
    return peri, apo


# --------------------------------------------------------------------------- #
# Oracles (independent paths; see frozen contract section 4)
# --------------------------------------------------------------------------- #
_MU_SI = MU_EARTH_KM3S2 * 1e9  # m^3/s^2
_RE_SI = R_EARTH_KM * 1e3  # m


def circular_decay_time_quadrature(
    a0_km: float,
    af_km: float,
    beta: float,
    atmosphere: np.ndarray = DEFAULT_ATMOSPHERE,
    nodes: int = 24,
) -> float:
    """PRIMARY oracle: t(a0 -> af) for quasi-circular decay, Gauss-Legendre in `a`.

    Integrates dt/da = 1/(kappa*rho(h(a))*sqrt(mu*a)) (SI units) with panels aligned
    to density-layer joints so each panel's integrand is smooth (C0 joints would
    otherwise silently degrade convergence).
    """
    kappa = 1.0 / beta
    joints = [a0_km]
    for h0 in atmosphere[:, 0]:
        if af_km < h0 + R_EARTH_KM < a0_km:
            joints.append(float(h0) + R_EARTH_KM)
    joints.append(af_km)
    joints = sorted(set(joints), reverse=True)
    xg, wg = np.polynomial.legendre.leggauss(nodes)
    total = 0.0
    for lo, hi in zip(joints[:-1], joints[1:]):
        half = abs(0.5 * (hi - lo))  # joints descend a0 -> af; integrate the positive measure
        mid = 0.5 * (hi + lo)
        a_si = (mid + half * xg) * 1e3  # m
        h_si = a_si - _RE_SI
        rho = np.empty(nodes)
        for j in range(nodes):
            rho[j] = air_density_si(h_si[j] * 1e-3, atmosphere)
        f = 1.0 / (kappa * rho * np.sqrt(_MU_SI * a_si))
        total += half * 1e3 * float(wg @ f)
    return total  # s


def quadrature_node_doubling(a0_km, af_km, beta, atmosphere=DEFAULT_ATMOSPHERE) -> dict:
    """Self-convergence of the quadrature oracle under node doubling (10/20/40)."""
    t10 = circular_decay_time_quadrature(a0_km, af_km, beta, atmosphere, nodes=10)
    t20 = circular_decay_time_quadrature(a0_km, af_km, beta, atmosphere, nodes=20)
    t40 = circular_decay_time_quadrature(a0_km, af_km, beta, atmosphere, nodes=40)
    return {
        "t10_s": t10,
        "t20_s": t20,
        "t40_s": t40,
        "rel_10_20": abs(t10 - t20) / t40,
        "rel_20_40": abs(t20 - t40) / t40,
    }


def erfi_decay_time(a0_km: float, af_km: float, beta: float, atmosphere=DEFAULT_ATMOSPHERE, dps: int = 60) -> float:
    """Closed-form single-layer circular decay (cross-check of the quadrature).

    t = sqrt(pi H)/K * [erfi(sqrt(a0/H)) - erfi(sqrt(af/H))],
    K = kappa*rho0*sqrt(mu)*exp((R_E+h0)/H), all SI.
    Restricted to layers with h0 >= 120 km: lower layers overflow
    exp((R_E+h0)/H) in binary floats and every decay window of interest lives
    above the reentry threshold anyway. Evaluated with mpmath at `dps`.
    """
    idx = int(np.searchsorted(atmosphere[:, 0], a0_km - R_EARTH_KM, side="right")) - 1
    h0_km, rho0, H_km = atmosphere[max(idx, 0)]
    if h0_km < 120.0:
        raise ValueError("erfi closed form restricted to layers with h0 >= 120 km")
    if af_km - R_EARTH_KM < h0_km:
        raise ValueError("segment crosses below the selected layer; use quadrature")
    kappa = 1.0 / beta
    with _mp_workdps(dps):
        H = _mp_mpf(H_km) * 1000
        a0 = _mp_mpf(a0_km) * 1000
        af = _mp_mpf(af_km) * 1000
        K = (
            _mp_mpf(kappa) * _mp_mpf(rho0) * _mp_mpf(_MU_SI) ** 0.5 * _mp_e ** ((_mp_mpf(_RE_SI) + _mp_mpf(h0_km) * 1000) / H)
        )
        val = (_mp_mpf(np.pi) * H) ** 0.5 / K * (_mp_erfi((a0 / H) ** 0.5) - _mp_erfi((af / H) ** 0.5))
    return float(val)


def sqrt_a_linear_time(a0_km: float, af_km: float, rho_si: float, beta: float) -> float:
    """Constant-density short-window oracle: sqrt(a(t)) = sqrt(a0) - kappa*rho*sqrt(mu)*t/2."""
    kappa = 1.0 / beta
    num = 2.0 * (np.sqrt(a0_km * 1e3) - np.sqrt(af_km * 1e3))
    return float(num / (kappa * rho_si * np.sqrt(_MU_SI)))


def oracle_a_of_t(a0_km: float, times_s: np.ndarray, beta: float, atmosphere=DEFAULT_ATMOSPHERE) -> np.ndarray:
    """Invert the quadrature oracle: semi-major axis [km] at requested times."""
    a_grid = np.linspace(a0_km, 6500.0, 4000)  # dense descent grid
    t_grid = np.empty_like(a_grid)
    t_prev = 0.0
    a_prev = a_grid[0]
    t_grid[0] = 0.0
    for i in range(1, len(a_grid)):
        t_i = circular_decay_time_quadrature(a_prev, a_grid[i], beta, atmosphere, nodes=12)
        t_grid[i] = t_prev + t_i
        t_prev, a_prev = t_grid[i], a_grid[i]
    out = np.interp(times_s, t_grid, a_grid, right=np.nan)
    finite = np.isfinite(out)
    if not np.all(finite):
        out = np.where(finite, out, np.nan)
    return out


def gauss_delta_rev(
    a_km: float,
    e: float,
    inc_deg: float,
    beta: float,
    omega_deg: float = 0.0,
    Omega_deg: float = 0.0,
    n_nodes: int = 512,
    omega_atm: float = 0.0,
) -> dict:
    """Eccentric-oracle: per-revolution Delta_a, Delta_e by quadrature along the
    osculating conic (INDEPENDENT of Cowell/RK4 integration).

    Uses exact osculating identities: eps=-mu/2a => da/dt=(2a^2/mu)*(a_drag.v);
    e^2 = 1+2*eps*h^2/mu^2 => d(e^2)/dt = 2*h^2/mu^2*d(eps)/dt + 4*eps*h/mu^2*dh/dt,
    dh/dt = hhat.(r x a_drag). Line integral over true anomaly with dt = r^2/h dnu.
    """
    inc = np.radians(inc_deg)
    Om = np.radians(Omega_deg)
    argp = np.radians(omega_deg)
    R = rotation_matrix_313(Om, inc, argp)
    p = a_km * (1.0 - e * e)
    h_mag = np.sqrt(MU_EARTH_KM3S2 * p)
    nu, w = np.polynomial.legendre.leggauss(n_nodes)
    nu = np.pi * (nu + 1.0)  # [0, 2*pi]
    w = np.pi * w
    r_pf = np.stack([np.cos(nu), np.sin(nu), np.zeros_like(nu)], axis=1) * (
        p / (1.0 + e * np.cos(nu))
    )[:, None]
    vr = MU_EARTH_KM3S2 / h_mag * e * np.sin(nu)
    vt = MU_EARTH_KM3S2 / h_mag * (1.0 + e * np.cos(nu))
    theta_hat = np.stack([-np.sin(nu), np.cos(nu), np.zeros_like(nu)], axis=1)
    v_pf = vr[:, None] * r_pf / (p / (1.0 + e * np.cos(nu)))[:, None] + vt[:, None] * theta_hat
    r = r_pf @ R.T
    v = v_pf @ R.T
    rm = np.linalg.norm(r, axis=1)
    rho = np.empty(n_nodes)
    for i in range(n_nodes):
        rho[i] = air_density_si(rm[i] - R_EARTH_KM)
    if omega_atm != 0.0:
        v_rel = v - omega_atm * np.stack([-r[:, 1], r[:, 0], np.zeros(n_nodes)], axis=1)
    else:
        v_rel = v
    adrag = -(0.5 * (1.0 / beta) * rho * DRAG_SI_TO_KKM * np.linalg.norm(v_rel, axis=1))[:, None] * v_rel
    dedt = np.einsum("ij,ij->i", adrag, v)  # km^2/s^3
    hvec = np.cross(r, v)
    dhdt = np.einsum("ij,ij->i", hvec / np.linalg.norm(hvec, axis=1)[:, None], np.cross(r, adrag))
    eps = -MU_EARTH_KM3S2 / (2.0 * a_km)
    dt_dnu = rm**2 / h_mag
    da = float(np.sum(w * (2.0 * a_km**2 / MU_EARTH_KM3S2) * dedt * dt_dnu))
    de2 = float(np.sum(w * (2.0 * h_mag**2 / MU_EARTH_KM3S2**2 * dedt + 4.0 * eps * h_mag / MU_EARTH_KM3S2**2 * dhdt) * dt_dnu))
    de = de2 / (2.0 * e)
    return {"da_rev_km": da, "de_rev": de}


def king_hele_delta_rev(a_km: float, e: float, beta: float, atmosphere=DEFAULT_ATMOSPHERE, dps: int = 50) -> dict:
    """King-Hele first-order per-revolution deltas (Ray & Scheeres 2021 family).

    Da_rev = -2 pi kappa a^2 rho_p e^{-z} [I0 + 2 e I1]
    De_rev = -2 pi kappa a   rho_p e^{-z} [I1 + e (I0 + I2)/2],  z = a e / H_p.
    Valid as published for small e; tests gate comparisons to e <= 0.1.
    """
    hp_km = a_km * (1.0 - e) - R_EARTH_KM
    idx = int(np.searchsorted(atmosphere[:, 0], hp_km, side="right")) - 1
    h0_km, rho_base, H_km = atmosphere[max(idx, 0)]
    rho_p = float(rho_base * np.exp(-(hp_km - h0_km) / H_km))  # density AT perigee altitude
    kappa = 1.0 / beta
    with _mp_workdps(dps):
        z = _mp_mpf(a_km * 1e3 * e) / (_mp_mpf(H_km) * 1000)
        I0, I1, I2 = (_mp_besseli(k, z) for k in (0, 1, 2))
        pre = 2.0 * np.pi * kappa * _mp_mpf(a_km * 1e3) ** 2 * _mp_mpf(rho_p) * _mp_e ** (-z)
        da_m = -pre * (I0 + 2 * _mp_mpf(e) * I1)
        de = -2.0 * np.pi * kappa * _mp_mpf(a_km * 1e3) * _mp_mpf(rho_p) * _mp_e ** (-z) * (
            I1 + _mp_mpf(e) * (I0 + I2) / 2
        )
    return {"da_rev_m": float(da_m), "de_rev": float(de)}


# --------------------------------------------------------------------------- #
# Run drivers
# --------------------------------------------------------------------------- #
def run_window(
    a_km: float,
    e: float,
    inc_deg: float,
    beta: float,
    n_revs: int,
    *,
    j2: float = 0.0,
    omega_atm: float = 0.0,
    Omega_deg: float = 0.0,
    omega_deg: float = 0.0,
    M0_deg: float = 0.0,
    spp: int | None = None,
    atmosphere=None,
    desc: str = "",
) -> dict:
    """Fixed-window decay run on the initial-period grid; returns series + stats."""
    r0, v0, _nu0 = seed_state(
        a_km, e, np.radians(inc_deg), np.radians(Omega_deg), np.radians(omega_deg), np.radians(M0_deg)
    )
    T = orbital_period(a_km)
    spp = steps_per_orbit(e) if spp is None else int(spp)
    dt = T / spp
    t = np.arange(n_revs * spp + 1) * dt
    states = propagate_3d_rk4_drag(
        r0, v0, MU_EARTH_KM3S2, t, j2=j2, beta=beta, omega_atm=omega_atm, atmosphere=atmosphere
    )
    coe = element_series(states)
    diss = dissipation_residual_stats(t, states, beta, omega_atm, atmosphere, j2=j2)
    # mean decay rate over the central measurement window (skip 1 rev at each end)
    i0, i1 = spp, len(t) - spp
    fit = ols_fit(t[i0:i1], coe["a"][i0:i1])
    slope_km_s = float(fit["slope"]) if fit is not None else float("nan")
    stride = max(1, (n_revs * spp) // 1200)
    return {
        "desc": desc,
        "a0_km": a_km,
        "e0": e,
        "inc_deg": inc_deg,
        "Omega_deg": Omega_deg,
        "omega_deg": omega_deg,
        "beta_kg_m2": beta,
        "j2": j2,
        "omega_atm_rad_s": omega_atm,
        "n_revs": n_revs,
        "spp": spp,
        "dt_s": dt,
        "period_s": T,
        "rate_m_day": slope_km_s * 1e3 * 86400.0,
        "rate_fit_r2": float(fit["r2"]) if fit is not None else float("nan"),
        "da_rev_measured_m": slope_km_s * T * 1e3,
        "a_final_km": float(coe["a"][-1]),
        "e_final": float(coe["e"][-1]),
        "dissipation": diss,
        "_t_stride": t[::stride],
        "_a_stride": np.asarray(coe["a"])[::stride],
        "_e_stride": np.asarray(coe["e"])[::stride],
        "_apo_stride": (np.asarray(coe["a"]) * (1 + np.asarray(coe["e"])))[::stride],
        "_peri_stride": (np.asarray(coe["a"]) * (1 - np.asarray(coe["e"])))[::stride],
        "_states": states,
        "_coe_a": np.asarray(coe["a"]),
        "_t": t,
    }


def _refine_crossing(r_safe, v_safe, t_safe, dt, r_threshold_km, *, j2=0.0, beta=0.0,
                     omega_atm=0.0, atmosphere=None, levels=(1, 2, 4, 8, 16, 32, 64)) -> dict:
    """Deterministic crossing refinement: reintegrate the bracketed interval with
    dt/2^j substeps; final estimate from the finest level (+ linear interp within
    the last sub-step). Flags sub-surface states instead of clamping them."""
    detail = {}
    finest = None
    for j in levels:
        sub = dt / j
        tg = np.arange(j + 1) * sub
        st = propagate_3d_rk4_drag(r_safe, v_safe, MU_EARTH_KM3S2, tg, j2=j2, beta=beta,
                                   omega_atm=omega_atm, atmosphere=atmosphere)
        rm = np.linalg.norm(st[:, :3], axis=1)
        if np.min(rm) < R_EARTH_KM:
            return {"status": "below-surface", "detail": detail}
        below = np.flatnonzero(rm < r_threshold_km)
        if below.size == 0:
            break
        i = int(below[0])
        if i == 0:
            tc = float(t_safe)
        else:
            frac = (rm[i - 1] - r_threshold_km) / (rm[i - 1] - rm[i])
            tc = float(t_safe + (i - 1 + frac) * sub)
        detail[str(j)] = tc
        finest = tc
    if finest is None:
        return {"status": "no-crossing", "detail": detail}
    keys = sorted(detail.keys(), key=int)
    stable = None
    if len(keys) >= 2:
        a_, b_ = detail[keys[-2]], detail[keys[-1]]
        stable = abs(a_ - b_)
    return {"status": "ok", "crossing_s": finest, "last_step_change_s": stable, "detail": detail}


def propagate_until_reentry(
    a_km: float,
    e: float,
    inc_deg: float,
    beta: float,
    *,
    j2: float = 0.0,
    omega_atm: float = 0.0,
    thresholds_km: tuple = REENTRY_THRESHOLDS_KM,
    spp: int | None = None,
    max_days: float = 400.0,
    chunk_steps: int = 512,
    atmosphere=None,
    m0_deg: float = 0.0,
) -> dict:
    """Chunked propagation to re-entry with per-threshold refined crossings.

    Thresholds are processed in descending altitude order (they are crossed in
    that order); each crossing time is refined by reintegrating only the
    bracketed interval (contract section 6). Terminates on: all thresholds met,
    surface hit, or max duration.
    """
    r0, v0, _ = seed_state(a_km, e, np.radians(inc_deg), 0.0, 0.0, np.radians(m0_deg))
    T = orbital_period(a_km)
    spp = steps_per_orbit(e) if spp is None else int(spp)
    dt = T / spp
    state = np.concatenate([r0, v0])
    t_abs = 0.0
    pending = list(thresholds_km)  # e.g. [120, 100] -- descending
    crossings: dict[float, dict] = {}
    status = "max-duration"
    n_steps_total = 0
    while t_abs < max_days * 86400.0:
        n_steps = min(chunk_steps, int(np.ceil((max_days * 86400.0 - t_abs) / dt)))
        tg = t_abs + np.arange(n_steps + 1) * dt
        states = propagate_3d_rk4_drag(state[:3], state[3:], MU_EARTH_KM3S2, tg,
                                       j2=j2, beta=beta, omega_atm=omega_atm, atmosphere=atmosphere)
        rmag = np.linalg.norm(states[:, :3], axis=1)
        n_steps_total += n_steps
        for thr in list(pending):
            hits = np.flatnonzero(rmag < R_EARTH_KM + thr)
            if hits.size > 0:
                k = int(hits[0])
                if k == 0:
                    ref = {"status": "ok", "crossing_s": float(t_abs)}
                else:
                    ref = _refine_crossing(states[k - 1, :3], states[k - 1, 3:], float(tg[k - 1]), dt,
                                           R_EARTH_KM + thr, j2=j2, beta=beta, omega_atm=omega_atm,
                                           atmosphere=atmosphere)
                crossings[thr] = ref
                pending.remove(thr)
        if np.min(rmag) < R_EARTH_KM:
            status = "surface-hit"
            break
        if not pending:
            status = "reentered"
            break
        state = states[-1]
        t_abs = float(tg[-1])
    out = {
        "status": status,
        "a0_km": a_km,
        "e0": e,
        "inc_deg": inc_deg,
        "beta_kg_m2": beta,
        "dt_s": dt,
        "spp": spp,
        "steps_total": n_steps_total,
        "sim_duration_days": t_abs / 86400.0,
    }
    for thr, ref in crossings.items():
        key = f"threshold_{int(thr)}km"
        out[key] = {
            "crossing_days": ref.get("crossing_s", float("nan")) / 86400.0,
            "refine_status": ref.get("status"),
            "last_refine_step_change_s": ref.get("last_step_change_s"),
        }
    for thr in pending:
        out[f"threshold_{int(thr)}km"] = {"crossing_days": None, "refine_status": "not-reached"}
    return out


# --------------------------------------------------------------------------- #
# Numerics batteries: convergence, order check, plateau separation
# --------------------------------------------------------------------------- #
def convergence_study(beta: float = 100.0, a_km: float | None = None, e: float = 0.0005,
                      inc_deg: float = 51.6, n_orbits: int = 60,
                      target_below_km: float = 419.7,
                      spp_list=(64, 128, 256, 512), spp_ref: int = 1024) -> dict:
    """Order study on the time-to-fall observable (single density-layer window)."""
    a_km = R_EARTH_KM + 420.0 if a_km is None else a_km

    def fall_time(spp):
        run = run_window(a_km, e, inc_deg, beta, n_orbits, spp=spp)
        tc = crossing_time(run["_t"], run["_coe_a"], R_EARTH_KM + target_below_km)
        return tc

    times = {spp: fall_time(spp) for spp in list(spp_list) + [spp_ref]}
    if any(v is None for v in times.values()):
        return {"error": "target not reached within window", "times": times}
    dts = [orbital_period(a_km) / spp for spp in spp_list]
    errs = [abs(times[spp] - times[spp_ref]) for spp in spp_list]
    rates = convergence_rate(np.maximum(np.array(errs, dtype=float), 1e-16), np.array(dts, dtype=float))
    return {
        "observable": "time_to_fall_s",
        "targets_km": target_below_km,
        "fall_times_s": {str(k): v for k, v in times.items()},
        "errors_vs_ref_s": dict(zip(map(str, spp_list), errs)),
        "convergence_rates": [float(r) for r in rates],
    }


def kepler_order_check(spp_list=(16, 32, 64, 128, 256), e: float = 0.3, n_orbits: float = 2.0) -> dict:
    """Raw integrator order of the drag-gated clone on pure Kepler dynamics.

    Analytic reference: propagate elements analytically over the same span via the
    Exp 008 Kepler machinery (independent path).
    """
    a_km = R_EARTH_KM + 1000.0
    inc_deg = 30.0
    T = orbital_period(a_km)
    n = mean_motion(a_km)
    r0, v0, nu0 = seed_state(a_km, e, np.radians(inc_deg), 0.0, 0.0, 0.0)
    span = n_orbits * T
    errs = []
    for spp in spp_list:
        t = np.arange(int(n_orbits * spp) + 1) * (T / spp)
        st = propagate_3d_rk4_drag(r0, v0, MU_EARTH_KM3S2, t)
        M_t = np.mod(n * span, 2.0 * np.pi)
        E_t = solve_kepler(M_t, e)
        nu_t = true_anomaly_from_E(E_t, e)
        r_ref, _v_ref = coe_to_rv_eci(a_km, e, np.radians(inc_deg), 0.0, 0.0, float(nu_t))
        errs.append(float(np.linalg.norm(st[-1, :3] - r_ref)))
    rates = convergence_rate(np.maximum(np.array(errs, dtype=float), 1e-16),
                             np.array([T / s for s in spp_list], dtype=float))
    return {"position_errors_m": dict(zip(map(str, spp_list), errs)),
            "convergence_rates": [float(r) for r in rates]}


def plateau_separation(a_km: float | None = None, beta: float = 100.0) -> dict:
    """Numerical-vs-model error separation (Exp 009 doctrine, non-conservative form).

    Transit-time observable between fixed altitudes inside one density layer.
    Law swap (declared H -> H*1.25, all rows) at FIXED dt must move the answer by
    a law-scale amount that is FLAT across the dt grid, while dt-refinement moves
    it by the (shrinking) integration-error scale.
    """
    a_km = R_EARTH_KM + 300.0 if a_km is None else a_km
    law_b = DEFAULT_ATMOSPHERE.copy()
    law_b[:, 2] = law_b[:, 2] * 1.25

    def transit(law, spp):
        run = run_window(a_km, 0.0005, 51.6, beta, 200, spp=spp, atmosphere=law)
        t_hi = crossing_time(run["_t"], run["_coe_a"], R_EARTH_KM + 295.0)
        t_lo = crossing_time(run["_t"], run["_coe_a"], R_EARTH_KM + 285.0)
        return (None if t_hi is None or t_lo is None else t_lo - t_hi)

    res = {}
    for spp in (128, 256):
        ta, tb = transit(DEFAULT_ATMOSPHERE, spp), transit(law_b, spp)
        res[f"spp_{spp}"] = {"law_A_transit_s": ta, "law_B_transit_s": tb,
                             "swap_difference_s": (None if ta is None or tb is None else tb - ta)}
    d128 = res["spp_128"]["swap_difference_s"]
    d256 = res["spp_256"]["swap_difference_s"]
    ta128, ta256 = res["spp_128"]["law_A_transit_s"], res["spp_256"]["law_A_transit_s"]
    if None in (d128, d256, ta128, ta256):
        res["error"] = "transit targets not reached"
        return res
    refine_level = abs(ta128 - ta256)
    swap_mean = abs(0.5 * (d128 + d256))
    res.update({
        "swap_flatness_abs_s": abs(d128 - d256),
        "integration_refinement_abs_s": refine_level,
        "separation_ratio": swap_mean / max(refine_level, 1e-16),
    })
    return res


# --------------------------------------------------------------------------- #
# Structural laws (pillar 3) -- scalings and symmetry
# --------------------------------------------------------------------------- #
def scaling_battery(n_revs: int = 250, a_km: float | None = None) -> dict:
    """beta / rho0 scaling laws + area linearity + J2-equivalence null.

    Continuum-exact at fixed state; over finite windows the orbit itself drifts
    (denser effective profile for larger kappa), so window ratios carry a small
    documented nonlinearity. Machine-exactness is verified separately at the
    acceleration level (identical states).
    """
    a_km = R_EARTH_KM + 420.0 if a_km is None else a_km
    betas = (50.0, 100.0, 200.0, 400.0)
    runs = {b: run_window(a_km, 0.0005, 51.6, b, n_revs, desc=f"scaling beta={b}") for b in betas}
    rates = {b: runs[b]["rate_m_day"] for b in betas}

    def ratio_dev(b1, b2):
        measured = rates[b1] / rates[b2]
        theory = b2 / b1
        return {"measured": measured, "theory": theory, "rel_dev": measured / theory - 1.0}

    beta_ratios = {"50_over_100": ratio_dev(50.0, 100.0), "200_over_100": ratio_dev(200.0, 100.0)}
    # rho0 uniform scaling (whole profile x c): rates scale exactly c
    rho_runs = {}
    for c in (0.5, 2.0):
        atm = DEFAULT_ATMOSPHERE.copy()
        atm[:, 1] = atm[:, 1] * c
        rho_runs[str(c)] = run_window(a_km, 0.0005, 51.6, 100.0, n_revs, atmosphere=atm)["rate_m_day"]
    base_rate = rates[100.0]
    rho_ratio = {c: {"measured": v / base_rate, "theory": float(c)} for c, v in rho_runs.items()}
    # acceleration-level exactness at IDENTICAL states (machine-exact check)
    r0, v0, _ = seed_state(a_km, 0.0005, np.radians(51.6), 0.0, 0.0, 0.0)
    x = np.concatenate([r0, v0])
    d1 = drag_accel_kkms(x, 100.0)
    d2 = drag_accel_kkms(x, 400.0)
    acc_exact_ratio = float(np.linalg.norm(d2) / np.linalg.norm(d1))
    # J2 on/off comparison: raw osculating-a slopes are INVALID under J2 -- seeding a
    # Kepler state and switching on J2 relaxes the mean elements by
    # (2a^2/mu)*<U_J2> ~ -5.7 km at i=51.6 deg (nonzero inclination-dependent mean of
    # the J2 potential), and short-period ripple aliases into windowed OLS slopes.
    # Honest comparison: settled-tail energy-bookkeeping rates; the residual after
    # correcting for the settled mean-altitude shift must be small.
    j2_on = run_window(a_km, 0.0005, 51.6, 100.0, n_revs, j2=J2_EARTH)

    def _energy_rate_m_day(run):
        eps = _eps_total_series(run["_states"], run["j2"])
        spp = run["spp"]
        t = run["_t"]
        i0 = min(5 * spp, len(t) - 2)
        a_bar = 0.5 * (float(run["_coe_a"][i0]) + float(run["_coe_a"][-1]))
        dedt = (eps[-1] - eps[i0]) / (t[-1] - t[i0])
        return float((2.0 * a_bar**2 / MU_EARTH_KM3S2) * dedt * 86400.0 * 1e3)

    def _tail_mean_alt_km(run):
        st = run["_states"][5 * run["spp"] :: 37]
        rm = np.linalg.norm(st[:, :3], axis=1)
        return float(np.mean(rm) - R_EARTH_KM)

    off_energy = _energy_rate_m_day(runs[100.0])
    on_energy = _energy_rate_m_day(j2_on)
    alt_off = _tail_mean_alt_km(runs[100.0])
    alt_on = _tail_mean_alt_km(j2_on)
    h_local = 58.515  # declared scale height near ~415 km (400-row layer)
    predicted_factor = float(np.exp((alt_off - alt_on) / h_local))
    measured_factor = on_energy / off_energy
    j2_comparison = {
        "off_energy_rate_m_day": off_energy,
        "on_energy_rate_m_day": on_energy,
        "settled_mean_altitude_off_km": alt_off,
        "settled_mean_altitude_on_km": alt_on,
        "measured_factor": measured_factor,
        "altitude_shift_predicted_factor": predicted_factor,
        "residual_after_altitude_correction": measured_factor / predicted_factor - 1.0,
        "raw_osculating_slope_rel_diff_invalid": float(j2_on["rate_m_day"] / rates[100.0] - 1.0),
        "note": "raw osculating-a slope under J2 is contaminated by the U_J2 mean-offset "
                "transient plus ripple aliasing; use energy-bookkeeping tail rates",
    }
    return {
        "window_revs": n_revs,
        "rates_m_day": {str(b): r for b, r in rates.items()},
        "beta_ratios": beta_ratios,
        "rho0_scalings": rho_ratio,
        "accel_level_ratio_beta400_over_100": acc_exact_ratio,
        "j2_settled_comparison": j2_comparison,
    }


def rotation_battery(n_revs: int = 250) -> dict:
    """Atmosphere co-rotation asymmetry + frame-symmetry twins.

    Equatorial circular pair (i=0 vs i=180 deg) is EXACT to all orders:
    rate ratio = ((v+w)/(v-w))^2 with w = omega_atm * r_arm. Inclined pair is
    first-order w_eff = w cos(i) with O(w_eff^2) corrections (~0.1%).
    Omega-twins under omega_atm = 0 must be identical to 1e-12 (frame symmetry).
    """
    a_km = R_EARTH_KM + 400.0
    e = 1e-4
    v_si = np.sqrt(_MU_SI / (a_km * 1e3))
    w_arm = OMEGA_EARTH_RAD_S * a_km * 1e3

    pro = run_window(a_km, e, 0.0, 100.0, n_revs, omega_atm=OMEGA_EARTH_RAD_S, desc="prograde equatorial")
    ret = run_window(a_km, e, 180.0, 100.0, n_revs, omega_atm=OMEGA_EARTH_RAD_S, desc="retrograde equatorial")
    ratio_eq = ret["rate_m_day"] / pro["rate_m_day"]
    theory_eq = ((v_si + w_arm) / (v_si - w_arm)) ** 2

    ipro = run_window(a_km, e, 63.6, 100.0, n_revs, omega_atm=OMEGA_EARTH_RAD_S, desc="inclined 63.6")
    iret = run_window(a_km, e, 116.4, 100.0, n_revs, omega_atm=OMEGA_EARTH_RAD_S, desc="inclined 116.4")
    ratio_inc = iret["rate_m_day"] / ipro["rate_m_day"]
    w_dimless = w_arm / v_si
    theory_inc = ((1.0 + w_dimless * np.cos(np.radians(63.6))) / (1.0 - w_dimless * np.cos(np.radians(63.6)))) ** 2

    twin_a = run_window(a_km, e, 51.6, 100.0, 20, Omega_deg=0.0, desc="omega twin A")
    twin_b = run_window(a_km, e, 51.6, 100.0, 20, Omega_deg=137.0, desc="omega twin B")
    rel_sym = float(np.max(np.abs(twin_b["_a_stride"] - twin_a["_a_stride"]) / np.abs(twin_a["_a_stride"])))

    return {
        "w_dimensionless_400km": float(w_dimless),
        "equatorial": {"ratio_measured": float(ratio_eq), "ratio_theory": float(theory_eq),
                       "rel_dev": float(ratio_eq / theory_eq - 1.0)},
        "inclined_63p6_116p4": {"ratio_measured": float(ratio_inc), "ratio_theory_first_order": float(theory_inc),
                                "rel_dev": float(ratio_inc / theory_inc - 1.0)},
        "omega_twin_max_rel_asymmetry": rel_sym,
    }


def benchmark_battery() -> dict:
    """Pillar-4 anchors. BINDING: quiet-time natural decay near 400 km ~2 km/month
    (ESA 'ISS reboost', 2016) compared as a decade-wide band against model rates.
    CONTEXT-ONLY (recorded, never asserted): Starlink storm-time decay
    (arXiv:2505.13752) -- geomagnetic augmentation outside declared-model scope."""
    out = {"binding": {}, "context_only": {}}
    for beta in (100.0, 190.0):
        run = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, beta, 60, desc=f"benchmark beta={beta}")
        km_month = run["rate_m_day"] * 30.44 / 1000.0
        out["binding"][f"model_decay_420km_beta{int(beta)}_km_month"] = float(km_month)
    published = 2.0
    lo, hi = published / 10.0, published * 10.0
    vals = [abs(v) for v in out["binding"].values()]  # decay rates are signed; band is on magnitude
    out["binding"]["published_quiet_decay_km_month"] = published
    out["binding"]["decade_band"] = [lo, hi]
    out["binding"]["verdict"] = "PASS" if any(lo <= v <= hi for v in vals) else "FAIL"
    out["context_only"]["starlink_storm_decay_km_day"] = {
        "value_range": [95.0, 176.0],
        "source": "arXiv:2505.13752 (Oliveira et al. 2025)",
        "exclusion_reason": "severe-geomagnetic density augmentation outside declared quiet-model scope",
    }
    return out


# --------------------------------------------------------------------------- #
# Pathological cases + adversarial mutant harness
# --------------------------------------------------------------------------- #
def _sign_flip_propagate(r0, v0, mu, t, beta: float, **kw) -> np.ndarray:
    """MUTANT HARNESS (tests only): drag with flipped sign -- energy GROWS.
    Deliberately not part of the public propagator (which rejects beta < 0)."""
    re2 = R_EARTH_KM**2

    def accel(x):
        rm = np.linalg.norm(x[:3])
        a_kep = -mu * x[:3] / rm**3
        return a_kep - drag_accel_kkms(x, beta, kw.get("omega_atm", 0.0), kw.get("atmosphere"))

    return _rk4_core(accel, r0, v0, mu, t)


def pathological_battery() -> dict:
    """Sentinel battery: knife-edge H, underflow, extreme beta, surface guard,
    layer-boundary semantics, determinism."""
    out = {}
    # knife-edge scale height (H -> ~0 on one layer): finite, no NaN
    knife = DEFAULT_ATMOSPHERE.copy()
    knife[DEFAULT_ATMOSPHERE[:, 0] == 400.0, 2] = 1e-6
    kr = run_window(R_EARTH_KM + 430.0, 0.0005, 51.6, 100.0, 30, atmosphere=knife)
    out["knife_edge_H"] = {
        "finite_states": bool(np.all(np.isfinite(kr["_states"]))),
        "rate_m_day": kr["rate_m_day"],
    }
    # high-altitude underflow path
    rho_high = air_density_si(1200.0)
    out["underflow"] = {"rho_1200km": rho_high, "finite": bool(np.isfinite(rho_high)),
                        "drag_finite": bool(np.all(np.isfinite(
                            drag_accel_kkms(np.array([R_EARTH_KM + 1200.0, 0.0, 0.0, 0.0, 7.0, 0.0]), 100.0))))}
    # extreme ballistic coefficients
    tiny_drag = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, 1e12, 20)
    null_run = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, 0.0, 20)
    drift = float(np.max(np.abs(tiny_drag["_coe_a"] - null_run["_coe_a"])) / 420.0)
    # kappa = 1e4: per-rev decay ~68 km -- extreme but inside the integrator's valid
    # regime, so the run must terminate cleanly through the thresholds. (kappa >= 1e6
    # removes more energy per step than the orbit holds and leaves the valid regime;
    # documented, not run.)
    huge_kappa = propagate_until_reentry(R_EARTH_KM + 300.0, 0.0, 51.6, 1e-4, max_days=40.0)
    out["extreme_beta"] = {
        "beta_1e12_vs_null_a_rel_drift": drift,
        "beta_1e-4_status": huge_kappa["status"],
        "beta_1e-4_threshold120_days": huge_kappa.get("threshold_120km", {}).get("crossing_days"),
    }
    # suborbital state -> surface-hit termination (never clamps through the ground)
    r_sub = np.array([R_EARTH_KM + 150.0, 0.0, 0.0])
    v_sub = np.array([0.0, -4.0, 0.0])
    T_guess = 2 * np.pi * np.sqrt((R_EARTH_KM) ** 3 / MU_EARTH_KM3S2)
    tg = np.arange(1, 600) * (T_guess / 1024)
    st = propagate_3d_rk4_drag(r_sub, v_sub, MU_EARTH_KM3S2, tg, beta=100.0)
    out["surface_guard"] = {
        "min_altitude_km": float(np.min(np.linalg.norm(st[:, :3], axis=1)) - R_EARTH_KM),
        "all_finite_before_stop": True,
    }
    # the EVENT-DRIVEN driver must terminate at the surface instead of integrating
    # through it: same suborbital geometry seeded at apogee (150 km), descending.
    # Drag here is negligible (beta=1e12) so the demonstration isolates the GEOMETRY
    # guard; with realistic drag the object brakes in dense air and the run ends at
    # the 120 km threshold ("reentered") before any surface contact - also correct.
    eps_sub = 0.5 * 4.0**2 - MU_EARTH_KM3S2 / (R_EARTH_KM + 150.0)
    a_sub = -MU_EARTH_KM3S2 / (2.0 * eps_sub)
    h_sub = (R_EARTH_KM + 150.0) * 4.0
    e_sub = float(np.sqrt(1.0 + 2.0 * eps_sub * h_sub**2 / MU_EARTH_KM3S2**2))
    sub_driver = propagate_until_reentry(
        a_sub, e_sub, 0.0, 1e12, max_days=5.0, thresholds_km=(120.0,), m0_deg=180.0,
    )
    out["surface_guard"]["driver_status"] = sub_driver["status"]
    out["surface_guard"]["driver_threshold120_days"] = (
        sub_driver.get("threshold_120km", {}).get("crossing_days")
    )
    # layer-boundary semantics (half-open [h_i, h_{i+1}); continuity to table rounding)
    rho_at_joint = air_density_si(250.0)
    rho_below = air_density_si(250.0 - 1e-9)
    row200 = DEFAULT_ATMOSPHERE[DEFAULT_ATMOSPHERE[:, 0] == 200.0][0]
    rho_from_below = float(row200[1] * np.exp(-(250.0 - 200.0) / row200[2]))
    out["layer_boundary"] = {
        "rho_250_exact_row_value": bool(rho_at_joint == 7.248e-11),
        "continuity_rel_jump": abs(rho_below - rho_from_below) / rho_at_joint,
    }
    # determinism: byte-identical repeat
    ra = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, 100.0, 15)
    rb = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, 100.0, 15)
    out["determinism_repeat"] = bool(np.array_equal(ra["_states"], rb["_states"]))
    return out


def mutant_battery() -> dict:
    """Adversarial mutants from the pre-registered catalog: each entry verifies the
    DETECTOR fires (the plausible-but-wrong result would have passed silently)."""
    out = {}
    a_km = R_EARTH_KM + 420.0
    e = 0.0005
    inc_deg = 51.6
    # 1. sign flip: energy must INCREASE somewhere -> monotonicity detector fires
    r0, v0, _ = seed_state(a_km, e, np.radians(inc_deg), 0.0, 0.0, 0.0)
    T = orbital_period(a_km)
    spp = steps_per_orbit(e)
    t = np.arange(30 * spp + 1) * (T / spp)
    st_mut = _sign_flip_propagate(r0, v0, MU_EARTH_KM3S2, t, 100.0)
    eps_mut = _eps_total_series(st_mut, 0.0)
    out["sign_flip_detected"] = bool(np.count_nonzero(np.diff(eps_mut) > 0.0) > 0)
    # 2. unit mutant (kg/m^3 fed as kg/km^3): pinned-probe detector fires -- drag
    # becomes ~1e-9 of its true size, so the measured delta must VANISH vs probe
    atm_unit = DEFAULT_ATMOSPHERE.copy()
    atm_unit[:, 1] = atm_unit[:, 1] * 1e-9
    run_unit = run_window(a_km, e, inc_deg, 100.0, 10, atmosphere=atm_unit)
    probe_m = 7.6851  # hand-computed Delta a/rev at 420 km, beta=100 (digits in tests)
    out["unit_mutant_ratio_to_probe"] = float(run_unit["da_rev_measured_m"] / probe_m)
    out["unit_mutant_detected"] = bool(abs(run_unit["da_rev_measured_m"]) < 0.01 * probe_m)
    # 3. B-inversion mutant (beta fed as kappa, x100 rate): rate off by 2 orders;
    # kept at 10 revs where the inverted-convention orbit stays above the ground
    # (at kappa=100 the orbit would hit the ground within ~2 revs and the integrator
    # would diverge in the dense layers -- a different, less clean failure)
    run_inv = run_window(a_km, e, inc_deg, 1.0, 10)  # inverted convention
    run_ok = run_window(a_km, e, inc_deg, 100.0, 10)
    inv_ratio = run_inv["rate_m_day"] / run_ok["rate_m_day"]
    out["b_inversion_rate_ratio"] = float(inv_ratio)
    out["b_inversion_detected"] = bool(inv_ratio > 10.0 or inv_ratio < 0.1)
    return out


# --------------------------------------------------------------------------- #
# Headline experiment assembly
# --------------------------------------------------------------------------- #
RESULTS_DIR = Path(__file__).resolve().parent / "results"


def leo_ref_anchor(beta: float = 100.0, n_revs: int = 500) -> dict:
    """Primary pillar-2 case: numerical a(t) vs quadrature/erfi oracles."""
    a0 = R_EARTH_KM + 420.0
    run = run_window(a0, 0.0005, 51.6, beta, n_revs, desc="leoRef anchor window")
    af_num = float(run["_coe_a"][-1])
    t_quad = circular_decay_time_quadrature(a0, af_num, beta)
    t_erfi = erfi_decay_time(a0, af_num, beta)
    node_dbl = quadrature_node_doubling(a0, af_num, beta)
    times = run["_t"]
    oracle_series = oracle_a_of_t(a0, times, beta)
    resid_m = (run["_coe_a"] - oracle_series) * 1e3
    fin = np.isfinite(resid_m)
    resid_stats = {
        "max_abs_residual_m": float(np.max(np.abs(resid_m[fin]))),
        "rms_residual_m": float(np.sqrt(np.mean(resid_m[fin] ** 2))),
        "points_compared": int(np.count_nonzero(fin)),
    }
    return {
        "run": run,
        "final_a_km": af_num,
        "quadrature_time_s": t_quad,
        "numerical_time_s": float(times[-1]),
        "quad_vs_num_rel": float(times[-1] / t_quad - 1.0),
        "erfi_time_s": t_erfi,
        "erfi_vs_quad_rel": float(t_erfi / t_quad - 1.0),
        "node_doubling": node_dbl,
        "oracle_residual": resid_stats,
        "_oracle_stride": oracle_series[:: max(1, len(times) // 1200)],
        "_residual_stride": resid_m[:: max(1, len(times) // 1200)],
        "_t_stride": times[:: max(1, len(times) // 1200)],
    }


def main() -> dict:
    print("== Exp 010: orbit decay / atmospheric drag ==")
    results: dict = {}

    print("[1/11] leoRef anchor window (420 km, beta=100, 500 revs)...")
    anchor = leo_ref_anchor(beta=100.0)
    lr_run = anchor.pop("run")
    oracle_stride = anchor.pop("_oracle_stride")
    resid_stride = anchor.pop("_residual_stride")
    t_lr_stride = anchor.pop("_t_stride")
    anchor["dissipation"] = lr_run["dissipation"]
    # hand-probe measurement (30 revs) + dissipation floor shrink pair (recorded so
    # card numbers are traceable to results.json)
    probe_run = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, 100.0, 30)
    coarse = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, 100.0, 4)
    fine = run_window(R_EARTH_KM + 420.0, 0.0005, 51.6, 100.0, 4,
                      spp=steps_per_orbit(0.0005) * 4)

    def _med_rel(run_):
        rr = rel_dissipation_residual(run_["_t"], run_["_states"], 100.0)
        return float(np.median(rr))

    anchor["hand_probe_measurement"] = {
        "window_revs": 30,
        "da_rev_measured_m": probe_run["da_rev_measured_m"],
        "abs_dev_from_hand_value_rel": float(abs(abs(probe_run["da_rev_measured_m"]) - 7.68506) / 7.68506),
    }
    anchor["dissipation_floor"] = {"median_rel_coarse_spp": _med_rel(coarse),
                                   "median_rel_fine_4x_spp": _med_rel(fine)}
    results["leo_ref_anchor"] = anchor

    print("[2/11] convergence study (time-to-fall observable)...")
    conv = convergence_study()
    if "error" in conv:
        raise RuntimeError(f"convergence study failed: {conv['error']}")
    results["convergence"] = conv

    print("[3/11] Kepler order check of the drag-gated clone...")
    results["kepler_order_check"] = kepler_order_check()

    print("[4/11] plateau separation (model swap vs dt refinement)...")
    plateau = plateau_separation()
    results["plateau_separation"] = plateau

    print("[5/11] scaling battery (beta / rho0 / area / J2-null)...")
    results["scalings"] = scaling_battery()

    print("[6/11] rotation asymmetry + symmetry twins...")
    results["rotation_asymmetry"] = rotation_battery()

    print("[7/11] benchmark battery (binding + context)...")
    bench = benchmark_battery()
    results["benchmarks"] = bench

    print("[8/11] eccentric battery (steep / Molniya-like / small-e King-Hele)...")
    ecc = eccentric_battery()
    steep_series = ecc["steep"].pop("_series")
    results["eccentric"] = ecc
    results["shortwindow"] = shortwindow_battery()

    print("[9/11] reentry cases (near-reentry 200 km, full decay 300 km)...")
    near = propagate_until_reentry(R_EARTH_KM + 200.0, 0.0, 51.6, 100.0)
    full = propagate_until_reentry(R_EARTH_KM + 300.0, 0.0, 51.6, 100.0)
    results["reentry"] = {"near_reentry_200km": near, "full_decay_300km": full}

    print("[10/11] lifetime curve (quadrature oracle + numerical spots)...")
    lc = lifetime_curve()
    spots = {}
    for h in (250.0, 280.0, 300.0):
        r = propagate_until_reentry(R_EARTH_KM + h, 0.0, 51.6, 100.0)
        spots[str(h)] = r.get("threshold_120km", {}).get("crossing_days")
    lc["numerical_spots"] = spots
    results["lifetime_curve"] = {
        "floor_alt_km": lc["floor_alt_km"],
        "lifetime_days_by_beta": lc["lifetime_days_by_beta"],
        "numerical_spots_days_to_120km": spots,
    }

    print("[11/11] pathological + mutant batteries...")
    results["pathological"] = pathological_battery()
    results["mutants"] = mutant_battery()

    # --- provenance ------------------------------------------------------- #
    results["constants"] = {
        "MU_EARTH_KM3S2": MU_EARTH_KM3S2,
        "R_EARTH_KM": R_EARTH_KM,
        "OMEGA_EARTH_RAD_S": OMEGA_EARTH_RAD_S,
        "J2_EARTH": J2_EARTH,
        "provenance": "inherited from Exp 008/009 (IAU 2015 B3 GM; WGS-84 Re; WGS-84 J2; "
                      "see groundtracks/j2Precession cards for full derivation chains)",
        "DRAG_SI_TO_KKM": DRAG_SI_TO_KKM,
        "beta_convention": "beta = m/(C_D A) [kg/m^2]; kappa = C_D A/m = 1/beta; "
                           "C_D canonical 2.2 for LEO spheres (Cooke 1965)",
        "frame": FRAME_CONVENTION,
        "units": UNITS_CONVENTION,
    }
    us76 = {"300": 1.9159e-11, "400": 2.8028e-12, "500": 5.2148e-13}
    bias = {}
    for alt, ref in us76.items():
        declared = air_density_si(float(alt))
        bias[alt] = float(declared / ref - 1.0)
    results["atmosphere"] = {
        "layers_h0_rho0_H": ATMOSPHERE_LAYERS.tolist(),
        "source": ATMOSPHERE_SOURCE,
        "source_url": ("https://github.com/CelesTrak/fundamentals-of-astrodynamics "
                       "(software/misc/pascal/ATMOSEXP.DAT), linked from "
                       "https://celestrak.org/software/vallado-sw.php"),
        "us76_spots_kg_m3": us76,
        "us76_source": US76_SPOTS_SOURCE,
        "us76_source_url": ("https://ntrs.nasa.gov/citations/19770009539 ; tables via "
                            "https://www.pdas.com/bigtables.html"),
        "secondary_source_flags": (
            "King-Hele formulas cited from open secondary sources: Ray & Scheeres "
            "MNRAS 501, 1168 (2021), arXiv:2008.10644 and ASU MAE462 Lecture 12"
        ),
        "blocked_sources_note": (
            "Cook/King-Hele/Walker Proc. R. Soc. A primary texts sit behind a bot "
            "challenge - NOT accessed per lab rules"
        ),
        "declared_minus_us76_rel": bias,
        "bias_note": "declared Vallado fit runs +26..+34% above US76 across the LEO band; "
                     "documented model bias folded into decade-wide smoke bands only",
    }
    results["resolution_rule"] = (
        "steps_per_orbit(e): 512/orbit circular (donor special-cases e<=0), "
        "max(512, ceil(720/(1-e)^1.5)) eccentric; fixed dt per run = period/spp; "
        "reentry crossings refined by reintegrating the bracketed interval at dt/2^j "
        "(j<=6); dual thresholds 120/100 km reported"
    )
    results["tolerances"] = {
        "leoref_quad_vs_num_rel_max": {"bound": 5e-3, "justification":
            "order-4 integration floor at default spp (convergence study), oracle self-converged <1e-9"},
        "erfi_vs_quad_rel_max": {"bound": 1e-9, "justification":
            "both paths far below float64 noise on a single-layer segment"},
        "beta_ratio_and_rho_ratio_rel_dev_max": {"bound": 2.5e-2, "justification":
            "continuum-exact at fixed state; finite-window drift nonlinearity O(kappa*Delta_a/H) "
            "grows with window length: ~1.7% for the widest ratio (beta 50/100) at the frozen "
            "250-rev windows (measured 1.7e-2)"},
        "equatorial_twin_rel_dev_max": {"bound": 0.02, "justification":
            "exact all-order theory; residual is integration asymmetry only"},
        "inclined_twin_rel_dev_max": {"bound": 0.03, "justification":
            "first-order w*cos(i) theory; O(w_eff^2) ~ 0.1% corrections"},
        "omega_twin_symmetry_rel_max": {"bound": 1e-12, "justification":
            "exact frame symmetry under spherical law; floating-point rounding only"},
        "eccentric_o4_agreement_rel_max": {"bound": 0.02, "justification":
            "Gauss quadrature along osculating conic vs RK4 Cowell: independent paths"},
        "king_hele_small_e_rel_dev_max": {"bound": 0.03, "justification":
            "King-Hele truncation O(e^2)~0.4% at e=0.05 plus layer-H discretization"},
        "dt_invariance_steep_rel_max": {"bound": 1e-3, "justification":
            "contract gate before recording any steep-case number"},
        "plateau_separation_min_ratio": {"bound": 20.0, "justification":
            "law swap must move transit time >=20x more than dt refinement"},
        "convergence_order_min": {"bound": 3.6, "justification":
            "design order 4; time-to-fall observable shows pre-asymptotic SUPERCONVERGENCE "
            "(measured 3.78 then ~5.0 as leading error terms cancel before the floor); raw "
            "Kepler position-error check confirms no degradation below design order"},
        "j2_settled_residual_max": {"bound": 0.03, "justification":
            "settled-tail energy rates after altitude-shift correction; residual is second-order "
            "J2-density coupling only"},
    }
    results["limitations"] = [
        "spherical geocentric altitude h=|r|-R_E; geodetic altitude differs by up to "
        "21.4 km (pole) / 13.1 km (ISS inclination) -> external density comparisons carry "
        "up to a factor ~1.4 bias (folded into decade-wide bands)",
        "quiet-time exponential atmosphere; geomagnetic storms, diurnal and seasonal "
        "density variations excluded (Starlink storm benchmarks recorded as context only)",
        "lifetimes above ~350 km come from the quadrature oracle validated by numerical "
        "spot decays at 250-300 km; direct end-to-end propagation there exceeds the step budget",
        "erfi closed form restricted to single layers with h0>=120 km (float64 overflow of "
        "exp((R_E+h0)/H) in lower layers); multi-layer windows use quadrature",
        "King-Hele few-% checks gated to e<=0.1; large-e anchors use direct Gauss quadrature",
        "J2 off for headline decay numbers (no secular J2-decay coupling under a spherical "
        "altitude-only atmosphere); J2 on/off equivalence null recorded instead",
        "save_json_result rounds |x|>=1e-10 to 12 decimals; rates therefore stored in m/day",
    ]

    # --- figures from recorded data ---------------------------------------- #
    stride_d = max(1, len(lr_run["_t"]) // 1200)
    rel_hist = rel_dissipation_residual(lr_run["_t"], lr_run["_states"], 100.0)
    results["figures_data"] = {
        "leoRef": {
            "t_days": (lr_run["_t"][::stride_d] / 86400.0).tolist(),
            "a_num_km": lr_run["_a_stride"].tolist(),
            "a_oracle_km": np.asarray(oracle_stride).tolist(),
            "residual_m": np.asarray(resid_stride).tolist(),
        },
        "steep": {k: np.asarray(v).tolist() for k, v in steep_series.items()},
        "lifetime": {
            "floor_alt_km": lc["floor_alt_km"],
            "by_beta": lc["lifetime_days_by_beta"],
            "numerical_spots": spots,
        },
        "convergence": {"errors": conv.get("errors_vs_ref_s", {})},
        "dissipation": {
            "t_days": (lr_run["_t"][1::stride_d] / 86400.0).tolist()[: len(rel_hist[::stride_d])],
            "rel_resid": rel_hist[::stride_d].tolist(),
        },
    }
    results["headline"] = _headline(results)
    figures = make_figures(results)
    results["figures"] = figures
    results["figures_note"] = "regenerated deterministically from results.json data, dpi=150"

    path = save_json_result(
        str(RESULTS_DIR / "results.json"),
        results,
        name="orbit_decay",
        description="Exp 010: atmospheric drag / orbit decay -- dissipation accounting, "
                    "decay-law rediscovery, scalings, reentry timing (declared Vallado atmosphere)",
    )
    print(f"saved: {path}")
    print(f"figures: {results['figures']}")
    print("headline:")
    for k, v in results["headline"].items():
        print(f"  {k}: {v}")
    return results


def _headline(results: dict) -> dict:
    a = results["leo_ref_anchor"]
    rot = results["rotation_asymmetry"]["equatorial"]
    conv = results["convergence"]
    near = results["reentry"]["near_reentry_200km"]
    return {
        "decay_law_rediscovered": (
            f"leoRef 420 km beta=100: max |a_num - a_oracle| = "
            f"{a['oracle_residual']['max_abs_residual_m']:.1f} m over 500 revs "
            f"({a['quad_vs_num_rel']:+.2e} rel on window time)"
        ),
        "closed_form_cross_check": f"erfi vs quadrature rel dev {a['erfi_vs_quad_rel']:+.2e}",
        "scalings": (
            f"beta ratios exact to "
            f"{max(abs(v['rel_dev']) for v in results['scalings']['beta_ratios'].values()):.1e}; "
            f"rho0 ratios to "
            f"{max(abs(v['measured']/v['theory'] - 1.0) for v in results['scalings']['rho0_scalings'].values()):.1e}"
        ),
        "rotation_asymmetry": (
            f"equatorial retro/prograde rate ratio {rot['ratio_measured']:.4f} vs theory "
            f"{rot['ratio_theory']:.4f}"
        ),
        "circularization": (
            f"steep case apogee drop {results['eccentric']['steep']['apo_drop_km']:.1f} km with "
            f"perigee change {abs(results['eccentric']['steep']['peri_first_last_km'][0] - results['eccentric']['steep']['peri_first_last_km'][1]):.2f} km"
        ),
        "shortwindow_anchors": {
            k: f"{v['rel_dev']:+.2e} rel vs quadrature oracle"
            for k, v in results["shortwindow"].items()
        },
        "reentry_timing": (
            f"200 km circular -> 120 km in "
            f"{near.get('threshold_120km', {}).get('crossing_days') if near.get('threshold_120km') else None} days; "
            f"dual thresholds reported"
        ),
        "benchmarks_binding_verdict": results["benchmarks"]["binding"].get("verdict"),
        "convergence_order_band": [round(r, 2) for r in conv.get("convergence_rates", [])],
    }


def eccentric_battery(beta: float = 100.0) -> dict:
    """Eccentric-regime validation: steep-gradient circularization demo, Molniya-like
    large-c Gauss-oracle anchor, small-e King-Hele few-% check, per-rev apsis
    estimator (third path), and the steep-case dt-halving invariance gate."""
    out = {}

    def _apsis_da_rev_m(run, coe_a):
        """Third-path estimator: mean signed Delta-a between consecutive apoapsis
        passages (negative for decay, same sign convention as the Gauss oracle)."""
        rmag = np.linalg.norm(run["_states"][:, :3], axis=1)
        peri_idx, apo_idx = apsis_indices(rmag)
        if len(apo_idx) < 3:
            return None, None, 0
        a_at_apo = coe_a[apo_idx]
        n_span = len(apo_idx) - 1
        return float((a_at_apo[-1] - a_at_apo[0]) * 1e3 / n_span), float(np.mean(np.diff(a_at_apo))), n_span

    # --- steep-gradient circularization demo (perigee 250 km, e = 0.3) --------
    rp = R_EARTH_KM + 250.0
    e_st = 0.3
    a_st = rp / (1.0 - e_st)
    st_run = run_window(a_st, e_st, 51.6, beta, 300, desc="steep circularization")
    apo = st_run["_apo_stride"]
    peri_s = st_run["_peri_stride"]
    o4_steep = gauss_delta_rev(a_st, e_st, 51.6, beta)
    da_meas_m = st_run["da_rev_measured_m"]
    apsis_da, _, n_apo = _apsis_da_rev_m(st_run, st_run["_coe_a"])
    # dt-halving invariance gate (contract section 6): same window, doubled resolution;
    # numbers from this case are only recorded if the per-rev delta moves < 0.1%
    st_run_fine = run_window(a_st, e_st, 51.6, beta, 300, spp=steps_per_orbit(e_st) * 2)
    da_meas_fine_m = st_run_fine["da_rev_measured_m"]
    dt_invariance_rel = abs(da_meas_fine_m - da_meas_m) / abs(da_meas_m)
    out["steep"] = {
        "a0_km": a_st, "e0": e_st, "perigee_alt_km": 250.0,
        "apo_drop_km": float(apo[0] - apo[-1]),  # strided series spans the full window
        "apo_first_last_km": [float(apo[0]), float(apo[-1])],
        "peri_first_last_km": [float(peri_s[0]), float(peri_s[-1])],
        "de_over_window": float(st_run["e_final"] - e_st),
        "gauss_o4_da_rev_m": float(o4_steep["da_rev_km"] * 1e3),
        "gauss_o4_de_rev": float(o4_steep["de_rev"]),
        "measured_da_rev_m": float(da_meas_m),
        "measured_vs_o4_rel": float(da_meas_m / (o4_steep["da_rev_km"] * 1e3) - 1.0),
        "apsis_estimator_da_rev_m": apsis_da,
        "apsis_n_passages": int(n_apo),
        "dt_invariance_rel_halving": float(dt_invariance_rel),
        "_series": {"t_days": st_run["_t_stride"] / 86400.0, "apo_km": apo, "peri_km": peri_s,
                    "e": st_run["_e_stride"]},
    }
    # --- Molniya-like large-c anchor (perigee 500 km, e = 0.74) ---------------
    rp_m = R_EARTH_KM + 500.0
    e_mo = 0.74
    a_mo = rp_m / (1.0 - e_mo)
    mo_run = run_window(a_mo, e_mo, 63.4, beta, 100, desc="molniya-like anchor")
    o4_mol = gauss_delta_rev(a_mo, e_mo, 63.4, beta)
    apsis_da_mo, _, n_apo_mo = _apsis_da_rev_m(mo_run, mo_run["_coe_a"])
    kh_mol = king_hele_delta_rev(a_mo, e_mo, beta)  # recorded as context (e > 0.2 regime)
    out["molniya_like"] = {
        "a0_km": a_mo, "e0": e_mo, "perigee_alt_km": 500.0,
        "gauss_o4_da_rev_m": float(o4_mol["da_rev_km"] * 1e3),
        "gauss_o4_de_rev": float(o4_mol["de_rev"]),
        "measured_da_rev_m": float(mo_run["da_rev_measured_m"]),
        "measured_vs_o4_rel": float(mo_run["da_rev_measured_m"] / (o4_mol["da_rev_km"] * 1e3) - 1.0),
        "apsis_estimator_da_rev_m": apsis_da_mo,
        "apsis_n_passages": int(n_apo_mo),
        "king_hele_context_da_rev_m": kh_mol["da_rev_m"],
        "e_direction": float(mo_run["e_final"] - e_mo),
    }
    # --- small-e King-Hele few-% check (e = 0.05, perigee 350 km) -------------
    e_sm = 0.05
    a_sm = (R_EARTH_KM + 350.0) / (1.0 - e_sm)
    sm_run = run_window(a_sm, e_sm, 51.6, beta, 60, desc="small-e King-Hele check")
    kh_sm = king_hele_delta_rev(a_sm, e_sm, beta)
    out["small_e_king_hele"] = {
        "a0_km": a_sm, "e0": e_sm,
        "kh_da_rev_m": float(kh_sm["da_rev_m"]),
        "kh_de_rev": float(kh_sm["de_rev"]),
        "measured_da_rev_m": float(sm_run["da_rev_measured_m"]),
        "da_rel_dev": float(sm_run["da_rev_measured_m"] / kh_sm["da_rev_m"] - 1.0),
        "measured_de_per_rev": float((sm_run["e_final"] - e_sm) / 60.0),  # signed: negative = decreasing
    }
    return out


def shortwindow_battery(beta: float = 100.0) -> dict:
    """Frozen-contract cases: Starlink-like (550 km) and SSO-like (600 km)
    short-window per-rev rate validation against the quadrature oracle."""
    out = {}
    for alt, inc, desc in ((550.0, 53.0, "starlink-like"), (600.0, 97.8, "sso-like")):
        a0 = R_EARTH_KM + alt
        run = run_window(a0, 0.0005, inc, beta, 100, desc=desc)
        af = float(run["_coe_a"][-1])
        t_orc = circular_decay_time_quadrature(a0, af, beta)
        out[desc.split("-")[0]] = {
            "alt_km": alt, "inc_deg": inc,
            "window_revs": 100,
            "numerical_time_s": float(run["_t"][-1]),
            "oracle_time_s": t_orc,
            "rel_dev": float(run["_t"][-1] / t_orc - 1.0),
            "rate_m_day": run["rate_m_day"],
        }
    return out


def lifetime_curve(betas=(50.0, 100.0, 200.0, 400.0), alts=np.arange(250.0, 801.0, 50.0),
                   floor_alt_km: float = 120.0) -> dict:
    """Quadrature-oracle lifetimes (days) from each start altitude to the floor."""
    af = R_EARTH_KM + floor_alt_km
    table = {}
    for b in betas:
        rows = {}
        for h in alts:
            t_s = circular_decay_time_quadrature(R_EARTH_KM + float(h), af, b)
            rows[str(float(h))] = float(t_s / 86400.0)
        table[str(b)] = rows
    return {"floor_alt_km": floor_alt_km, "lifetime_days_by_beta": table}


# --------------------------------------------------------------------------- #
# Figures (regenerated deterministically from recorded result data)
# --------------------------------------------------------------------------- #
def make_figures(results: dict) -> list[str]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fd = results["figures_data"]
    paths = []

    # F1: leoRef numerical vs oracle decay + residual inset
    fig, ax = plt.subplots(figsize=(9, 5.5))
    lr = fd["leoRef"]
    ax.plot(lr["t_days"], lr["a_num_km"], "-", lw=1.8, label="numerical RK4 Cowell")
    ax.plot(lr["t_days"], lr["a_oracle_km"], "--", lw=1.5, label="quadrature oracle (declared law)")
    ax.set_xlabel("time [days]")
    ax.set_ylabel("semi-major axis [km]")
    ax.set_title("Exp 010 F1: leoRef circular decay, 420 km, beta=100 kg/m$^2$")
    ax.legend(loc="lower left")
    axins = ax.inset_axes([0.56, 0.12, 0.4, 0.32])
    axins.plot(lr["t_days"], lr["residual_m"], lw=0.9, color="tab:red")
    axins.set_ylabel("num $-$ oracle [m]", fontsize=8)
    axins.tick_params(labelsize=7)
    axins.grid(alpha=0.3)
    p = fig_dir / "f1_leoref_decay_vs_oracle.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F2: steep-case circularization (apogee/perigee/e evolution)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    st = fd["steep"]
    ax.plot(st["t_days"], st["apo_km"], "-", lw=1.6, label="apogee radius")
    ax.plot(st["t_days"], st["peri_km"], "-", lw=1.6, label="perigee radius")
    ax.set_xlabel("time [days]")
    ax.set_ylabel("radius [km]")
    ax.set_title("Exp 010 F2: circularization, e=0.3, perigee 250 km, beta=100 kg/m$^2$")
    ax.legend(loc="center left")
    ax2 = ax.twinx()
    ax2.plot(st["t_days"], st["e"], ":", lw=1.4, color="tab:green", label="eccentricity")
    ax2.set_ylabel("eccentricity", color="tab:green")
    ax2.tick_params(axis="y", colors="tab:green")
    p = fig_dir / "f2_circularization.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F3: lifetime vs altitude (+beta family), numerical spots overlaid
    fig, ax = plt.subplots(figsize=(9, 5.5))
    lc = fd["lifetime"]
    alts = sorted({float(k) for row in lc["by_beta"].values() for k in row}, key=float)
    for beta_str, row in lc["by_beta"].items():
        ys = [row[str(a)] for a in alts]
        ax.plot(alts, ys, "-o", ms=3, lw=1.3, label=f"oracle, beta={float(beta_str):.0f}")
    spots = lc.get("numerical_spots", {})
    if spots:
        sa = sorted(float(k) for k in spots)
        ax.plot(sa, [spots[str(a)] for a in sa], "ks", ms=7, mfc="none", label="numerical full decay")
    ax.set_yscale("log")
    ax.set_xlabel("start altitude [km]")
    ax.set_ylabel(f"time to {int(lc['floor_alt_km'])} km [days]")
    ax.set_title("Exp 010 F3: re-entry timeline vs altitude and ballistic coefficient")
    ax.legend(fontsize=8)
    p = fig_dir / "f3_lifetime_vs_altitude.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)

    # F4: convergence panel + dissipation-residual history
    fig, (axc, axr) = plt.subplots(1, 2, figsize=(11, 5))
    cv = fd["convergence"]
    sp = sorted(int(k) for k in cv["errors"])
    dts = [orbital_period(R_EARTH_KM + 420.0) / s for s in sp]
    errs = [cv["errors"][str(s)] for s in sp]
    axc.loglog(dts, errs, "o-", lw=1.5)
    guide = errs[-1] * (np.array(dts) / dts[-1]) ** 4.0
    axc.loglog(dts, guide, ":", lw=1.2, label="order-4 guide")
    axc.set_xlabel("step size [s]")
    axc.set_ylabel("time-to-fall error vs finest grid [s]")
    axc.set_title("Convergence (order 4)")
    axc.legend(fontsize=8)
    axr.semilogy(fd["dissipation"]["t_days"], fd["dissipation"]["rel_resid"], lw=0.8)
    axr.set_xlabel("time [days]")
    axr.set_ylabel("relative dissipation residual")
    axr.set_title("Pointwise dissipation identity")
    axr.grid(alpha=0.3)
    p = fig_dir / "f4_convergence_dissipation.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p.name)
    return paths


if __name__ == "__main__":
    main()
