"""Experiment 011 -- Lagrange points / circular restricted three-body problem.

Question: Does an independently derived, rotating-frame CR3BP model rediscover the
five classical Lagrange equilibrium positions, conserve the Jacobi integral under
RK4 propagation at the documented convergence order, reproduce the linear stability
classification (collinear unstable for every mass ratio; triangular stable iff
27*mu*(1-mu) < 1) together with its nonlinear perturbation signatures, and survive
a dimensional/nondimensional cross-check plus a pre-registered adversarial mutant
battery -- while transitioning the laboratory from single-primary dynamics into
rotating-frame multi-body machinery?

Frozen contract v1.0 (2026-08-22), six-track research panel + adversarial review.
Frames/units: nondimensional barycentric rotating frame; mu = m2/(m1+m2) in (0, 1/2];
+x from barycenter toward m2; omega = +z_hat (prograde/CCW from +z); primaries at
(-mu, 0, 0) and (1-mu, 0, 0); length unit = primary separation a, time unit = 1/n,
mass unit = m1+m2 (so n = 1, G(m1+m2) = 1). Effective potential (positive-style)
omega_eff = (1-mu)/r1 + mu/r2 + (x^2+y^2)/2 with r1 to the LARGER primary m1.
Jacobi constant C = 2*omega_eff - |v_rot|^2; exact scalings C_dim = n^2 a^2 C_nondim
and inertial bridge C = 2(n h_z - E_I). L4 = LEADING point = (1/2 - mu, +sqrt(3)/2).

Validation is layered per the anti-shared-algebra doctrine: scalar root-finding is
cross-checked against vector residuals AND two algebraically independent quintic
families AND mpmath 40-dps anchors; stability eigenvalues against closed forms;
conservation against the inertial bridge identity E_I + C/2 = L_z (E_I itself is
NOT conserved -- moving primaries do net work; drift law dE_I/dt =
mu(1-mu) y (r1^-3 - r2^-3) is verified instead); frame handling against an
inertial-consistency round trip that is the ONLY reliable killer of Coriolis-sign
mutants (spectra and Jacobi conservation are provably blind to it).
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lab_utils.integrators import rk4_propagate
from lab_utils.metrics import convergence_rate
from lab_utils.results import save_json_result

# --- Reuse of verified prior-experiment machinery (single-hop importlib) ---- #
_J2_PATH = Path(__file__).resolve().parents[1] / "j2Precession" / "experiment.py"
_j2_spec = importlib.util.spec_from_file_location("j2_for_lagrange", _J2_PATH)
assert _j2_spec is not None and _j2_spec.loader is not None
_j2 = importlib.util.module_from_spec(_j2_spec)
_j2_spec.loader.exec_module(_j2)

ols_fit = _j2.ols_fit

# --- Declared physical constants (portfolio canon, provenance Exp 008/009) --- #
GM_SUN_KM3S2 = 1.32712440018e11  # IAU 2015 nominal GM_Sun (km^3/s^2)
GM_EARTH_KM3S2 = 398600.4418  # IAU 2015 nominal GM_E (Exp 008 canon)
GM_MOON_KM3S2 = 4902.800066  # DE-series nominal GM_Moon
A_EM_KM = 384400.0  # classical mean Earth-Moon separation (idealization: circular)
AU_KM = 149597870.7  # IAU 2012 exact astronomical unit (km)

MU_EM = GM_MOON_KM3S2 / (GM_EARTH_KM3S2 + GM_MOON_KM3S2)  # 0.012150584077904827
MU_SEM = (GM_EARTH_KM3S2 + GM_MOON_KM3S2) / (
    GM_SUN_KM3S2 + GM_EARTH_KM3S2 + GM_MOON_KM3S2
)  # 3.040423452319562e-06 (Sun-(Earth+Moon); Sun-Earth-only ratio would be 3.0035e-6)
MU_ROUTH = (9.0 - math.sqrt(69.0)) / 18.0  # 0.0385208965045514 triangular threshold
N_EM_RAD_S = math.sqrt((GM_EARTH_KM3S2 + GM_MOON_KM3S2) / A_EM_KM**3)
T_EM_DAYS = 2.0 * math.pi / N_EM_RAD_S / 86400.0

MU_LITERATURE_NOTE = (
    "CR3BP literature often quotes mu_EM = 0.012150585609624 (JPL SOLAR SYSTEM table "
    "vintage); the GM-derived value used here differs by 1.53e-9 absolute "
    "(1.26e-7 relative). One convention pinned throughout so the dimensional "
    "cross-check is internally exact."
)

CASES = {
    "earth_moon": MU_EM,
    "sun_earth_moon": MU_SEM,
    "routh_boundary": MU_ROUTH,
    "above_threshold": 0.05,
    "equal_mass": 0.5,
    "mu_1e-3": 1e-3,
    "mu_1e-6": 1e-6,
}

FRAME_CONVENTION = (
    "nondim CR3BP: barycentric rotating frame, +x barycenter->m2, omega=+z_hat, "
    "m1>=m2 at (-mu,0,0)/(1-mu,0,0), mu=m2/(m1+m2), L=a, T=1/n, n=1"
)

TOLERANCES = {
    "equilibrium_residual_inf": 2e-14,
    "quintic_residual_rel": 1e-13,
    "mpmath_anchor_abs": 5e-16,
    "roundtrip_abs": 1e-14,
    "inertial_consistency_rel": 1e-10,
    "mutant_signal_min": 1e-2,
    "bridge_identity_abs": 1e-13,
    "drift_law_rel": 1e-4,
    "dim_crosscheck_rel": 1e-12,
    "jacobi_order_lo": 3.5,
    "jacobi_order_hi": 5.0,
    "growth_rate_rel_eps1e-4": 5e-3,
    "lp_frequency_rel": 1e-6,
    "zvc_bias_rel": 1e-2,
    "zvc_drift_guard": 1e-6,
}

SQRT3_2 = math.sqrt(3.0) / 2.0


# --------------------------------------------------------------------------- #
# Rotating-frame core physics
# --------------------------------------------------------------------------- #
def omega_eff(x, y, z, mu):
    """Effective potential (positive-style): gravity + centrifugal barrier."""
    r1 = np.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = np.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    return (1.0 - mu) / r1 + mu / r2 + 0.5 * (x * x + y * y)


def grad_omega(x, y, z, mu):
    """grad(omega_eff): equals the force side (gravity + outward centrifugal)."""
    r1 = np.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = np.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    return np.array(
        [
            x - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3,
            y - (1.0 - mu) * y / r1**3 - mu * y / r2**3,
            -(1.0 - mu) * z / r1**3 - mu * z / r2**3,
        ]
    )


def jacobi_constant(state, mu):
    """C = 2*omega_eff - |v_rot|^2 (state = [x,y,z,vx,vy,vz], nondimensional)."""
    x, y, z, vx, vy, vz = state
    return float(2.0 * omega_eff(x, y, z, mu) - (vx * vx + vy * vy + vz * vz))


def rhs_rotating(mu):
    """Autonomous rotating-frame RHS: r'' + 2 zhat x r' = grad(omega_eff)."""

    def f(t, s):
        g = grad_omega(s[0], s[1], s[2], mu)
        return np.array([s[3], s[4], s[5], 2.0 * s[4] + g[0], -2.0 * s[3] + g[1], g[2]])

    return f


def rot_z(th):
    c, s = math.cos(th), math.sin(th)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_to_inert(state, th):
    """Rotating -> inertial: r_I = R(th) r_R; v_I = R(th)(v_R + zhat x r_R)."""
    R = rot_z(th)
    r = R @ state[:3]
    v = R @ (state[3:] + np.array([-state[1], state[0], 0.0]))
    return np.concatenate([r, v])


def inert_to_rot(state, th):
    """Inertial -> rotating: r_R = R(-th) r_I; v_R = R(-th) v_I - zhat x r_R."""
    R = rot_z(-th)
    r = R @ state[:3]
    v = R @ state[3:] - np.array([-r[1], r[0], 0.0])
    return np.concatenate([r, v])


def rhs_inertial(mu):
    """Newton RHS in the INERTIAL barycentric frame with primaries on circular
    rails p_i(t) = R(t) p_i0. Time-dependent THROUGH the RK4 stages (freezing
    stage-time-dependent primaries silently reduces the scheme to order 1)."""

    def f(t, s):
        R = rot_z(t)
        p1 = R @ np.array([-mu, 0.0, 0.0])
        p2 = R @ np.array([1.0 - mu, 0.0, 0.0])
        d1 = s[:3] - p1
        d2 = s[:3] - p2
        a = -(1.0 - mu) * d1 / np.linalg.norm(d1) ** 3 - mu * d2 / np.linalg.norm(d2) ** 3
        return np.concatenate([s[3:], a])

    return f


def inertial_state_quantities(state_i, mu, th=0.0):
    """Specific energy E_I and angular momentum h_z in the inertial frame.
    Primary distances use the MOVING primaries p_i(th) = R(th) p_i0 -- the
    static rotating-axis formula is only valid at th = 0."""
    x, y, z, vx, vy, vz = state_i
    R = rot_z(th)
    p1 = R @ np.array([-mu, 0.0, 0.0])
    p2 = R @ np.array([1.0 - mu, 0.0, 0.0])
    r1 = float(np.linalg.norm(state_i[:3] - p1))
    r2 = float(np.linalg.norm(state_i[:3] - p2))
    energy = 0.5 * (vx * vx + vy * vy + vz * vz) - (1.0 - mu) / r1 - mu / r2
    h_z = x * vy - y * vx
    return energy, h_z


# --------------------------------------------------------------------------- #
# Equilibrium solvers
# --------------------------------------------------------------------------- #
def collinear_scalar(x, mu):
    """Scalar equilibrium equation on the axis: f(x) = 0 (piecewise-defined)."""
    d1 = abs(x + mu)
    d2 = abs(x - 1.0 + mu)
    return x - (1.0 - mu) * (x + mu) / d1**3 - mu * (x - 1.0 + mu) / d2**3


def collinear_scalar_prime(x, mu):
    """f'(x) = 1 + 2(1-mu)/|x+mu|^3 + 2 mu/|x-1+mu|^3 >= 1 (strict monotonicity)."""
    d1 = abs(x + mu)
    d2 = abs(x - 1.0 + mu)
    return 1.0 + 2.0 * (1.0 - mu) / d1**3 + 2.0 * mu / d2**3


def _bisect(f, a, b, iters=200):
    fa = f(a)
    if fa == 0.0:
        return a
    for _ in range(iters):
        m = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0 or (b - m) < 1e-17 * max(1.0, abs(m)):
            return m
        if (fa < 0.0) != (fm < 0.0):
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def _solve_bracketed(name, mu, lo_fn, hi_fn):
    """Bracket-secure collinear root: honors exact endpoint roots (mu=1/2 puts
    L1 exactly at 1/2-mu), shrinks the inner offset only while needed."""
    flo = collinear_scalar(lo_fn(), mu)
    if flo == 0.0:
        return lo_fn()
    d = (mu / 3.0) ** (1.0 / 3.0) / 1024.0
    for _ in range(64):
        hi = hi_fn(d)
        fhi = collinear_scalar(hi, mu)
        if fhi == 0.0:
            return hi
        if flo * fhi < 0.0:
            x = _bisect(lambda xx: collinear_scalar(xx, mu), lo_fn(), hi)
            for _ in range(50):
                fx = collinear_scalar(x, mu)
                step = fx / collinear_scalar_prime(x, mu)
                xn = x - step
                if not (lo_fn() <= xn <= hi) or xn == x:
                    break
                x = xn
            return min(max(x, lo_fn()), hi)
        if d < 1e-300:
            break
        d /= 16.0
    raise RuntimeError(f"{name} bracket failed for mu={mu}")


def solve_collinear(mu):
    """L1/L2/L3 via proven brackets [Track B theorem] + safeguarded Newton.

    Brackets: L1 in [1/2-mu, 1-mu], L2 in [1-mu, 2-mu], L3 in [-1-mu, -1];
    endpoint signs are exact for all mu in (0, 1/2]. Newton polish keeps iterates
    bracketed; monotonicity (f' >= 1) guarantees uniqueness per interval.
    At mu = 1/2 the L1 root sits exactly ON the lower endpoint (f(1/2-mu) = 0).
    """
    roots = {
        "L1": _solve_bracketed(
            "L1", mu, lambda: 0.5 - mu, lambda d: 1.0 - mu - d
        ),
        "L2": _solve_bracketed(
            "L2", mu, lambda: 1.0 - mu + (mu / 3.0) ** (1.0 / 3.0) / 1024.0, lambda d: 2.0 - mu
        ),
        "L3": _solve_bracketed(
            "L3", mu, lambda: -1.0 - mu, lambda d: -1.0
        ),
    }
    return roots


def lagrange_points(mu):
    """All five equilibria for mass ratio mu (frozen convention: L4 leading +y)."""
    col = solve_collinear(mu)
    return {
        "L1": np.array([col["L1"], 0.0, 0.0]),
        "L2": np.array([col["L2"], 0.0, 0.0]),
        "L3": np.array([col["L3"], 0.0, 0.0]),
        "L4": np.array([0.5 - mu, SQRT3_2, 0.0]),
        "L5": np.array([0.5 - mu, -SQRT3_2, 0.0]),
    }


def gamma_quintic_residual(name, x, mu):
    """Cleared-denominator quintic in the branch distance variable (independent
    algebra family #1; derived by substitution + multiplication through)."""
    if name == "L1":
        g = 1.0 - mu - x
        return g**5 - (3 - mu) * g**4 + (3 - 2 * mu) * g**3 - mu * g**2 + 2 * mu * g - mu
    if name == "L2":
        g = x - (1.0 - mu)
        return g**5 + (3 - mu) * g**4 + (3 - 2 * mu) * g**3 - mu * g**2 - 2 * mu * g - mu
    if name == "L3":
        g = -mu - x
        return g**5 + (2 + mu) * g**4 + (1 + 2 * mu) * g**3 - (1 - mu) * g**2 - 2 * (1 - mu) * g - (1 - mu)
    raise ValueError(name)


def direct_x_quintic_residual(name, x, mu):
    """Direct-x quintic P(x) = x(x+mu)^2(x-1+mu)^2 + s1(1-mu)(x-1+mu)^2 +
    s2 mu (x+mu)^2 (independent algebra family #2; sign pattern per branch)."""
    signs = {"L1": (-1.0, 1.0), "L2": (-1.0, -1.0), "L3": (1.0, 1.0)}
    s1, s2 = signs[name]
    return (
        x * (x + mu) ** 2 * (x - 1.0 + mu) ** 2
        + s1 * (1.0 - mu) * (x - 1.0 + mu) ** 2
        + s2 * mu * (x + mu) ** 2
    )


# --------------------------------------------------------------------------- #
# Linear stability machinery
# --------------------------------------------------------------------------- #
def hessian_entries(x0, y0, mu):
    """Second derivatives of omega_eff at (x0, y0, z=0) [closed forms]."""
    dx1, dx2 = x0 + mu, x0 - 1.0 + mu
    r1 = math.hypot(dx1, y0)
    r2 = math.hypot(dx2, y0)
    r13, r15, r23, r25 = r1**3, r1**5, r2**3, r2**5
    wxx = 1.0 - (1 - mu) * (1.0 / r13 - 3 * dx1 * dx1 / r15) - mu * (1.0 / r23 - 3 * dx2 * dx2 / r25)
    wyy = 1.0 - (1 - mu) * (1.0 / r13 - 3 * y0 * y0 / r15) - mu * (1.0 / r23 - 3 * y0 * y0 / r25)
    wxy = 3.0 * ((1 - mu) * dx1 * y0 / r15 + mu * dx2 * y0 / r25)
    wzz = -(1.0 - mu) / r13 - mu / r23
    return wxx, wyy, wxy, wzz


def stability_matrix(x0, y0, mu):
    """Planar 4x4 Jacobian at an equilibrium (accel block enters POSITIVE;
    Coriolis coupling antisymmetric: +2 at (vx-row, vy-col), -2 mirrored)."""
    wxx, wyy, wxy, _ = hessian_entries(x0, y0, mu)
    return np.array(
        [
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [wxx, wxy, 0.0, 2.0],
            [wxy, wyy, -2.0, 0.0],
        ]
    )


def vertical_frequency(x0, y0, mu):
    """omega_z = sqrt((1-mu)/r1^3 + mu/r2^3): vertical pair lambda = +/- i omega_z."""
    dx1, dx2 = x0 + mu, x0 - 1.0 + mu
    r1 = math.hypot(dx1, y0)
    r2 = math.hypot(dx2, y0)
    return math.sqrt((1.0 - mu) / r1**3 + mu / r2**3)


def collinear_closed_form_rates(x0, mu):
    """Closed-form sigma (real pair), nu (center pair), from A := (1-mu)/r1^3 +
    mu/r2^3 > 1: discriminant of the biquadratic is A(9A-8)."""
    dx1, dx2 = x0 + mu, x0 - 1.0 + mu
    big_a = (1.0 - mu) / abs(dx1) ** 3 + mu / abs(dx2) ** 3
    disc = big_a * (9.0 * big_a - 8.0)
    sigma = math.sqrt(((big_a - 2.0) + math.sqrt(disc)) / 2.0)
    nu = math.sqrt(((2.0 - big_a) + math.sqrt(disc)) / 2.0)
    omega_z = math.sqrt(big_a)
    return sigma, nu, omega_z, big_a


def triangular_closed_form_rates(mu):
    """lambda^4 + lambda^2 + gamma = 0 with gamma = (27/4) mu (1-mu).
    gamma <= 1/4 (Routh): purely imaginary spectrum, rates nu_short/nu_long.
    gamma >  1/4: complex quartet lambda = +alpha +/- i*beta, -alpha +/- i*beta
    (unstable); nu_* are NaN and the quartet rates are returned instead."""
    gamma = 0.25 * 27.0 * mu * (1.0 - mu)
    disc = 1.0 - 4.0 * gamma
    stable = 27.0 * mu * (1.0 - mu) < 1.0
    if disc >= 0.0:
        nu_short = math.sqrt((1.0 + math.sqrt(disc)) / 2.0)
        nu_long = math.sqrt((1.0 - math.sqrt(disc)) / 2.0)
        quartet_alpha = quartet_beta = float("nan")
    else:
        nu_short = nu_long = float("nan")
        sq_g = math.sqrt(gamma)
        quartet_alpha = math.sqrt((sq_g - 0.5) / 2.0)
        quartet_beta = math.sqrt((sq_g + 0.5) / 2.0)
    return nu_short, nu_long, gamma, stable, quartet_alpha, quartet_beta


def eigen_analysis(x0, y0, mu):
    """Numeric eigenvalues (planar 4x4 + decoupled vertical pair), classified by
    clustering real parts (ordering of eigvals output is NOT assumed)."""
    A = stability_matrix(x0, y0, mu)
    lam4 = np.linalg.eigvals(A)
    wz = vertical_frequency(x0, y0, mu)
    lam = np.concatenate([lam4, np.array([1j * wz, -1j * wz])])
    scale = max(1.0, float(np.max(np.abs(lam))))
    real_type = [l for l in lam if abs(l.imag) < 1e-8 * scale]
    imag_type = sorted(abs(float(l.imag)) for l in lam if l.imag != 0.0 and abs(l.real) < 1e-8 * scale)
    complex_type = [l for l in lam if abs(l.imag) >= 1e-8 * scale and abs(l.real) >= 1e-8 * scale]
    max_real = float(np.max(lam.real))
    return lam, max_real, sorted(float(l.real) for l in real_type), imag_type, complex_type


# --------------------------------------------------------------------------- #
# Independent validation oracles
# --------------------------------------------------------------------------- #
def inertial_bridge_drift(traj_rot, times, mu):
    """Vectorized bridge/oracle evaluation along a history.

    Returns (max|C - 2(hz - E_I)|, max|E_I + C/2 - L_z|) where both sides use
    states mapped to the inertial frame at each sample's own epoch.
    """
    worst_bridge = 0.0
    worst_identity = 0.0
    for s, th in zip(traj_rot, times):
        si = rot_to_inert(s, float(th))
        c = jacobi_constant(s, mu)
        energy, hz = inertial_state_quantities(si, mu, float(th))
        worst_bridge = max(worst_bridge, abs(c - 2.0 * (hz - energy)))
        worst_identity = max(worst_identity, abs(energy + 0.5 * c - hz))
    return worst_bridge, worst_identity


def inertial_energy_history(traj_rot, times, mu):
    """E_I along a mapped history + cumulative drift-law integral.
    dE_I/dt = mu(1-mu) y_R (r1^-3 - r2^-3) is evaluated in ROTATING components
    (derived: dE/dt = sum_i GM_i (r-p_i).pdot_i / r^3 with pdot_i = w x p_i,
    whose triple products reduce exactly to the rotating-frame form)."""
    e_hist = np.empty(len(times))
    power = np.empty(len(times))
    for k, (s, th) in enumerate(zip(traj_rot, times)):
        si = rot_to_inert(s, float(th))
        e_hist[k], _ = inertial_state_quantities(si, mu, float(th))
        x, y = s[0], s[1]
        r1 = math.hypot(x + mu, y)
        r2 = math.hypot(x - 1.0 + mu, y)
        power[k] = mu * (1.0 - mu) * y * (1.0 / r1**3 - 1.0 / r2**3)
    drift_integral = float(np.trapezoid(power, times)) if hasattr(np, "trapezoid") else float(np.trapz(power, times))
    return e_hist, drift_integral


def zvc_wall_position(c_target, mu, x_lo=None, x_hi=None):
    """Axis wall of the zero-velocity structure: unique x in (1-mu, x_L2) with
    2*Omega(x,0) = c_target (2*Omega strictly decreasing there). Bisection."""
    pts = lagrange_points(mu)
    if x_lo is None:
        x_lo = 1.0 - mu + 1e-9
    if x_hi is None:
        x_hi = float(pts["L2"][0]) - 1e-12

    def g(x):
        return 2.0 * omega_eff(x, 0.0, 0.0, mu) - c_target

    assert g(x_lo) > 0.0 > g(x_hi), "wall not bracketed"
    return _bisect(g, x_lo, x_hi)


def critical_escape_speed(x_s, c_target, mu):
    """Speed along +x at (x_s, 0) whose Jacobi constant equals c_target."""
    return math.sqrt(max(0.0, 2.0 * omega_eff(x_s, 0.0, 0.0, mu) - c_target))


def propagate_with_guard(f, t_grid, s0, mu, drift_guard):
    """Propagate; if Jacobi drift exceeds the guard, retry once at dt/4.
    Returns (traj, times_used, max_drift)."""
    traj = rk4_propagate(f, t_grid, s0)
    drift = max(abs(jacobi_constant(s, mu) - jacobi_constant(s0, mu)) for s in traj)
    if drift <= drift_guard or t_grid[1] - t_grid[0] < 1e-9:
        return traj, t_grid, drift
    fine = np.linspace(t_grid[0], t_grid[-1], 2 * (len(t_grid) - 1) + 1)
    traj_f = rk4_propagate(f, fine, s0)
    drift_f = max(abs(jacobi_constant(s, mu) - jacobi_constant(s0, mu)) for s in traj_f)
    return traj_f, fine, drift_f


# --------------------------------------------------------------------------- #
# Studies
# --------------------------------------------------------------------------- #
def long_period_ic(amp: float, mu: float) -> tuple[np.ndarray, float]:
    """Initial condition on the L4 long-period normal mode [Track C recipe]:
    u0 = X_L4 + amp * [Re w0, Re w1, Re w2, Re w3] / scale where w is the
    first-order eigenvector closest to lambda = -i*nu_long. Re(w) automatically
    carries the modal velocity layout (w[2] = -i*nu*w[0] => Re w2 = nu Im w0)."""
    pts = lagrange_points(mu)
    _, nu_l, _, _, _, _ = triangular_closed_form_rates(mu)
    A4 = stability_matrix(pts["L4"][0], pts["L4"][1], mu)
    w4, V4 = np.linalg.eig(A4)
    idx = int(np.argmin(np.abs(w4.imag - (-nu_l))))
    wl = w4[idx]
    w = V4[:, idx]
    scale = max(abs(w[0]), abs(w[1]))
    pos_off = np.array([w[0].real, w[1].real]) / scale
    vel_off = (-wl.imag) * np.array([w[0].imag, w[1].imag]) / scale
    s0 = np.concatenate([pts["L4"][:2] + amp * pos_off, [0.0], amp * vel_off, [0.0]])
    return s0, nu_l


def study_equilibria():
    out = {"cases": {}, "anchors": {}, "tolerances": TOLERANCES}
    for case, mu in CASES.items():
        pts = lagrange_points(mu)
        case_out = {}
        for name, pos in pts.items():
            g = grad_omega(pos[0], pos[1], pos[2], mu)
            entry = {
                "x": float(pos[0]),
                "y": float(pos[1]),
                "residual_inf": float(np.max(np.abs(g))),
            }
            if name in ("L1", "L2", "L3"):
                entry["gamma_quintic"] = float(gamma_quintic_residual(name, pos[0], mu))
                entry["direct_x_quintic"] = float(direct_x_quintic_residual(name, pos[0], mu))
                entry["bracket_ok"] = bool(
                    (name == "L1" and 0.5 - mu <= pos[0] <= 1.0 - mu)
                    or (name == "L2" and 1.0 - mu <= pos[0] <= 2.0 - mu)
                    or (name == "L3" and -1.0 - mu <= pos[0] <= -1.0)
                )
            case_out[name] = entry
        xs = [case_out[n]["x"] for n in ("L1", "L2", "L3")]
        case_out["ordering_ok"] = bool(-1.0 - mu < xs[2] < -mu < xs[0] < 1.0 - mu < xs[1] < 2.0 - mu)
        case_out["gamma2_gt_gamma1"] = bool((xs[1] - (1 - mu)) > ((1 - mu) - xs[0]))
        case_out["c_values"] = {
            name: jacobi_constant(np.concatenate([pos, np.zeros(3)]), mu) for name, pos in pts.items()
        }
        out["cases"][case] = case_out

    # ordering of critical Jacobi values C1 > C2 > C3 > C4 = C5.
    # At the equal-mass ENDPOINT the configuration is mirror-degenerate:
    # C2 = C3 exactly (L2/L3 mirror pair), so only that middle link relaxes.
    c_ord_ok = True
    for case in CASES:
        cv = out["cases"][case]["c_values"]
        if case == "equal_mass":
            c_ord_ok &= cv["L1"] > cv["L2"] and abs(cv["L2"] - cv["L3"]) < 1e-12 and cv["L3"] > cv["L4"]
        else:
            c_ord_ok &= cv["L1"] > cv["L2"] > cv["L3"] > cv["L4"]
        c_ord_ok &= abs(cv["L4"] - cv["L5"]) < 1e-12
    out["c_ordering_ok"] = bool(c_ord_ok)

    # mpmath 40-dps anchors for Earth-Moon
    from mpmath import mp, mpf, findroot

    mp.dps = 40
    mu_m = mpf(MU_EM)

    def f_ax(x):
        return x - (1 - mu_m) * (x + mu_m) / abs(x + mu_m) ** 3 - mu_m * (x - 1 + mu_m) / abs(x - 1 + mu_m) ** 3

    anchors = {}
    seeds = {"L1": mpf("0.84"), "L2": mpf("1.16"), "L3": mpf("-1.005")}
    for name, seed in seeds.items():
        xa = findroot(f_ax, seed)
        anchors[name] = float(xa)
    em = out["cases"]["earth_moon"]
    out["mpmath_anchor_dev"] = {
        name: abs(em[name]["x"] - anchors[name]) for name in ("L1", "L2", "L3")
    }

    # mission/solar-system anchors (scale+geometry, percent-level honesty)
    gamma1 = 1.0 - MU_EM - em["L1"]["x"]
    gamma2 = em["L2"]["x"] - (1.0 - MU_EM)
    gamma3 = -MU_EM - em["L3"]["x"]
    out["anchors"] = {
        "em_L1_km_from_moon": gamma1 * A_EM_KM,
        "em_L1_published_km": 58000.0,
        "em_L2_km_from_moon": gamma2 * A_EM_KM,
        "em_L2_published_km": 64500.0,
        "em_L3_km_from_earth": abs(em["L3"]["x"] + MU_EM) * A_EM_KM,
        "em_L3_published_km": 381700.0,
        "sem_L1_km_from_earth": (1.0 - MU_SEM - out["cases"]["sun_earth_moon"]["L1"]["x"]) * AU_KM,
        "sem_L2_km_from_earth": (out["cases"]["sun_earth_moon"]["L2"]["x"] - (1.0 - MU_SEM)) * AU_KM,
        "sel_published_km": 1.5e6,
        "trojan_leading_angle_deg_em": math.degrees(math.atan2(SQRT3_2, 0.5 - MU_EM)),
        "trojan_note": "L4 leads m2 by atan2(sqrt3/2, 1/2-mu): 60.607 deg at EM mu, -> 60 deg as mu->0 (Jupiter Trojans swarm +-60 deg)",
    }
    return out


def study_stability():
    out = {"cases": {}}
    for case, mu in CASES.items():
        pts = lagrange_points(mu)
        case_out = {}
        for name in ("L1", "L2", "L3"):
            x0 = float(pts[name][0])
            lam, max_real, real_pair, imag_pairs, _ = eigen_analysis(x0, 0.0, mu)
            sigma, nu, wz, big_a = collinear_closed_form_rates(x0, mu)
            num_sigma = max(real_pair) if real_pair else float("nan")
            case_out[name] = {
                "sigma_numeric": num_sigma,
                "sigma_closed": sigma,
                "nu_closed": nu,
                "nu_numeric_imag": imag_pairs[:2],
                "omega_z_closed": wz,
                "A_curvature": big_a,
                "max_real_eigenvalue": max_real,
                "classification": "unstable (saddle x center)" if max_real > 1e-12 else "stable",
                "closed_vs_numeric_rel": abs(num_sigma - sigma) / sigma,
            }
        ns, nl, gamma, stable, qa, qb = triangular_closed_form_rates(mu)
        x4, y4 = 0.5 - mu, SQRT3_2
        lam4_, max_real4, _, imag4, cx4 = eigen_analysis(x4, y4, mu)
        case_out["L4"] = {
            "nu_short_closed": ns,
            "nu_long_closed": nl,
            "quartet_alpha_closed": qa,
            "quartet_beta_closed": qb,
            "nu_numeric_imag": imag4,
            "numeric_quartet": sorted(
                ([float(l.real), float(l.imag)] for l in cx4), key=lambda p: (p[0], p[1])
            ),
            "gamma": gamma,
            "routh_stable": bool(stable),
            "max_real_eigenvalue": max_real4,
            "vertical_freq_numeric": 1.0 if any(abs(abs(float(l.imag)) - 1.0) < 1e-9 and l.real != 0.0 for l in lam4_) else None,
            "classification": (
                "stable (center x center)" if stable else "unstable (complex quartet)"
            ),
        }
        case_out["routh_threshold"] = MU_ROUTH
        case_out["routh_product"] = 27.0 * mu * (1.0 - mu)
        out["cases"][case] = case_out

    # boundary degeneracy demonstration (mpmath): gamma(mu_R) - 1/4
    from mpmath import mp, mpf

    mp.dps = 50
    mu_r = mpf(9) / 18 - mp.sqrt(69) / 18
    gam = mpf(27) / 4 * mu_r * (1 - mu_r)
    out["boundary_gamma_minus_quarter_50dps"] = float(gam - mpf(1) / 4)
    return out


def _ladder_orders(drifts, stepsizes):
    mask = np.asarray(drifts) > 0
    return convergence_rate(np.maximum(np.asarray(drifts)[mask], 1e-20), np.asarray(stepsizes)[mask])


def study_jacobi_drift():
    """Three trajectory classes, dt ladders into the roundoff floor [Track D]."""
    mu = MU_EM
    f = rhs_rotating(mu)
    pts = lagrange_points(mu)
    out = {}

    # class A: L4 long-period mode, amp=3e-3 (bounded quasi-periodic)
    s0_A, nu_lp = long_period_ic(3e-3, mu)
    span_A = 10.0 * 2 * math.pi / nu_lp
    ladder_A = {}
    for k in (256, 1024, 2048, 4096):
        tg = np.linspace(0.0, span_A, k + 1)
        tr = rk4_propagate(f, tg, s0_A)
        c0 = jacobi_constant(s0_A, mu)
        ladder_A[k] = max(abs(jacobi_constant(s, mu) - c0) for s in tr)
    ks_A = sorted(ladder_A)
    out["class_L4_LP_amp3e-3"] = {
        "span": span_A,
        "drifts": ladder_A,
        "orders": _ladder_orders(
            [ladder_A[k] for k in ks_A], [span_A / k for k in ks_A]
        ).tolist(),
    }

    # class A2: floor plateau at amp=1e-3 (drift quantized at ~2 ulp of C0,
    # constant across the ladder -- Track D class-A floor evidence)
    s0_A2, _ = long_period_ic(1e-3, mu)
    ladder_A2 = {}
    for k in (256, 1024, 4096, 16384):
        tg = np.linspace(0.0, span_A, k + 1)
        tr = rk4_propagate(f, tg, s0_A2)
        c0 = jacobi_constant(s0_A2, mu)
        ladder_A2[k] = max(abs(jacobi_constant(s, mu) - c0) for s in tr)
    out["class_L4_LP_floor_plateau"] = {
        "span": span_A,
        "drifts": ladder_A2,
        "plateau_max_over_min": max(ladder_A2.values()) / min(ladder_A2.values()),
    }

    # class B: temporary capture near L1, C just below C2
    c2 = jacobi_constant(np.concatenate([pts["L2"], np.zeros(3)]), mu)
    c_target = c2 - 1e-3
    x_s = 0.80
    vy = math.sqrt(2.0 * omega_eff(x_s, 0.0, 0.0, mu) - c_target)
    s0_B = np.array([x_s, 0.0, 0.0, 0.0, vy, 0.0])
    span_B = 30.0
    ladder_B = {}
    for n in (1024, 2048, 4096, 8192, 16384):
        tg = np.linspace(0.0, span_B, n + 1)
        tr = rk4_propagate(f, tg, s0_B)
        c0 = jacobi_constant(s0_B, mu)
        ladder_B[n] = max(abs(jacobi_constant(s, mu) - c0) for s in tr)
    ns_B = sorted(ladder_B)
    out["class_capture_C2minus1e-3"] = {
        "span": span_B,
        "c_target_minus_c2": -1e-3,
        "release_x": x_s,
        "vy": vy,
        "drifts": ladder_B,
        "orders": _ladder_orders(
            [ladder_B[k] for k in ns_B], [span_B / k for k in ns_B]
        ).tolist(),
    }

    # class C: retrograde circlet around m2, rho=0.30 (inside Hill sphere)
    rho = 0.30
    s0_C = np.array([1.0 - mu, rho, 0.0, math.sqrt(mu / rho) + rho, 0.0, 0.0])
    t_kep = 2.0 * math.pi * math.sqrt(rho**3 / mu)
    span_C = 5.0 * t_kep
    ladder_C = {}
    floor_run = {}
    for k in (512, 1024, 2048, 4096, 8192, 16384, 32768):
        tg = np.linspace(0.0, span_C, k + 1)
        tr = rk4_propagate(f, tg, s0_C)
        c0 = jacobi_constant(s0_C, mu)
        drift = max(abs(jacobi_constant(s, mu) - c0) for s in tr)
        if k <= 8192:
            ladder_C[k] = drift
        else:
            floor_run[k] = drift
    ks_C = sorted(ladder_C)
    out["class_retrograde_rho0p30"] = {
        "span": span_C,
        "t_keplerian": t_kep,
        "drifts": ladder_C,
        "orders": _ladder_orders(
            [ladder_C[k] for k in ks_C], [span_C / k for k in ks_C]
        ).tolist(),
        "floor_probe": floor_run,
    }
    return out, s0_A, s0_B, s0_C


def study_perturbations():
    """Analytic linear theory vs nonlinear RK4 perturbation experiments [Track C]."""
    mu = MU_EM
    f = rhs_rotating(mu)
    out = {}

    # --- unstable growth at L1 ---
    pts = lagrange_points(mu)
    X1 = pts["L1"]
    A_mat = stability_matrix(X1[0], 0.0, mu)
    w, V = np.linalg.eig(A_mat)
    idx = int(np.argmax(w.real))
    sigma = float(w[idx].real)
    v_hat = V[:, idx].real  # real eigenvector of the real pair (imag ~ 0)
    assert np.max(np.abs(V[:, idx].imag)) < 1e-10
    v_hat = v_hat / np.max(np.abs(v_hat))
    growth = {}
    for eps in (1e-4, 1e-6, 1e-8):
        # planar eigenvector (dx, dy, dvx, dvy) applied around the spatial point
        s0 = np.concatenate(
            [
                X1[:2] + eps * v_hat[:2],
                [0.0],
                eps * v_hat[2:],
                [0.0],
            ]
        )
        c_top = 100.0 * eps
        t_end_guess = math.log(100.0) / sigma
        tg = np.linspace(0.0, t_end_guess, 20001)
        tr = rk4_propagate(f, tg, s0)
        proj = (tr[:, :2] - X1[:2]) @ np.array([v_hat[0], v_hat[1]])
        inside = (proj >= 2.0 * eps) & (proj <= c_top)
        tt = tg[inside]
        cc = proj[inside]
        fit = ols_fit(tt, np.log(cc))
        slope = fit["slope"]
        growth[f"eps_{eps:g}"] = {
            "sigma_linear": sigma,
            "slope_measured": slope,
            "rel_error": abs(slope - sigma) / sigma,
            "bias_ratio_consistency": None,
            "r2": fit["r2"],
        }
    biases = [growth[k]["rel_error"] for k in growth]
    growth["bias_scaling_ratios"] = [
        biases[i] / biases[i + 1] for i in range(len(biases) - 1)
    ]
    out["L1_unstable_growth"] = growth

    # --- stable long-period oscillation at L4 ---
    nu_s, nu_l, _, _, _, _ = triangular_closed_form_rates(mu)
    freq_fit = {}
    span = 6.0 * 2.0 * math.pi / nu_l
    for amp in (1e-3, 1e-4, 1e-5):
        s0, _ = long_period_ic(amp, mu)
        n_steps = 12000
        tg = np.linspace(0.0, span, n_steps * 6 // 5 + 1)  # dt ~ T_long/12000
        tr = rk4_propagate(f, tg, s0)
        c0 = jacobi_constant(s0, mu)
        cj_drift = max(abs(jacobi_constant(s, mu) - c0) for s in tr)
        sig = tr[:, 0] - np.mean(tr[:, 0])

        def resid(fr):
            cos_t = np.cos(2 * math.pi * fr * tg)
            sin_t = np.sin(2 * math.pi * fr * tg)
            denom = np.dot(cos_t, cos_t) * np.dot(sin_t, sin_t) - np.dot(cos_t, sin_t) ** 2
            a_c = np.dot(sig, cos_t) / np.dot(cos_t, cos_t)
            b_s = np.dot(sig, sin_t) / np.dot(sin_t, sin_t)
            return float(np.sum((sig - a_c * cos_t - b_s * sin_t) ** 2))

        coarse = np.linspace(nu_l / (2 * math.pi) * 0.9, nu_l / (2 * math.pi) * 1.1, 201)
        best = min(coarse, key=resid)
        for _ in range(6):
            fine = np.linspace(best - (coarse[1] - coarse[0]), best + (coarse[1] - coarse[0]), 41)
            best = min(fine, key=resid)
            coarse = fine
        nu_fit = 2 * math.pi * best
        freq_fit[f"amp_{amp:g}"] = {
            "nu_long_linear": nu_l,
            "nu_measured": nu_fit,
            "rel_error": abs(nu_fit - nu_l) / nu_l,
            "cj_drift": cj_drift,
        }
    amps = [1e-3, 1e-4, 1e-5]
    errs = [freq_fit[f"amp_{a:g}"]["rel_error"] for a in amps]
    # quadratic-in-amplitude extrapolation to zero amplitude
    c_quad = errs[1] / amps[1] ** 2 if errs[1] > 0 else 0.0
    freq_fit["quadratic_coefficient_check"] = errs[0] / errs[1]
    freq_fit["extrapolated_rel_error"] = abs(errs[-1] - c_quad * amps[-1] ** 2) / nu_l
    out["L4_longperiod_mode"] = freq_fit
    return out


def study_dimensional_crosscheck():
    """Same physical Earth-Moon case represented dimensionally vs nondimensionally."""
    mu = MU_EM
    gm_tot = GM_EARTH_KM3S2 + GM_MOON_KM3S2
    n = N_EM_RAD_S
    L = A_EM_KM
    out = {"gm_total": gm_tot, "n_rad_s": n, "L_km": L, "T_days": T_EM_DAYS}

    # dimensional scalar equation (km): equilibrium directly in km units
    def f_dim(x_km):
        xm1 = -mu * L
        xm2 = (1.0 - mu) * L
        d1 = abs(x_km - xm1)
        d2 = abs(x_km - xm2)
        return (
            n**2 * x_km
            - GM_EARTH_KM3S2 * (x_km - xm1) / d1**3
            - GM_MOON_KM3S2 * (x_km - xm2) / d2**3
        ) / n**2

    pts = lagrange_points(mu)
    dim_eq = {}
    for name, br in (("L1", (0.5 * L, (1 - mu) * L)), ("L2", ((1 - mu) * L, 2.0 * L)), ("L3", (-2.0 * L, -mu * L))):
        dim_eq[name] = _bisect(f_dim, br[0] + 1.0, br[1] - 1.0)
        nondim_mapped = float(pts[name][0]) * L
        dim_eq[name + "_mapped"] = nondim_mapped
        dim_eq[name + "_abs_diff_km"] = abs(dim_eq[name] - nondim_mapped)
        dim_eq[name + "_rel_to_L"] = dim_eq[name + "_abs_diff_km"] / L

    # C scaling at equilibria
    def c_dim_axis(x_km, y_km=0.0, vx=0.0, vy=0.0):
        xm1, xm2 = -mu * L, (1.0 - mu) * L
        r1 = math.hypot(x_km - xm1, y_km)
        r2 = math.hypot(x_km - xm2, y_km)
        U_dim = GM_EARTH_KM3S2 / r1 + GM_MOON_KM3S2 / r2 + 0.5 * n**2 * (x_km**2 + y_km**2)
        return 2.0 * U_dim - ((vx) ** 2 + (vy) ** 2)

    scale = n**2 * L**2
    c_scal = {}
    for name in ("L1", "L2", "L3"):
        c_dim = c_dim_axis(dim_eq[name])
        c_nond = jacobi_constant(np.concatenate([pts[name], np.zeros(3)]), mu)
        c_scal[name] = abs(c_dim / scale - c_nond) / abs(c_nond)
    x4 = float(pts["L4"][0]) * L
    y4 = SQRT3_2 * L
    c_dim4 = c_dim_axis(x4, y4)
    c_nond4 = jacobi_constant(np.concatenate([pts["L4"], np.zeros(3)]), mu)
    c_scal["L4"] = abs(c_dim4 / scale - c_nond4) / abs(c_nond4)
    out["jacobi_scaling_rel"] = c_scal

    # 90-day trajectory correspondence around the L4 long-period mode
    s0_nd, _ = long_period_ic(1e-3, mu)
    span_nd = 90.0 * 86400.0 * n
    tg_nd = np.linspace(0.0, span_nd, 4001)
    tr_nd = rk4_propagate(rhs_rotating(mu), tg_nd, s0_nd)
    s0_dim = np.concatenate([s0_nd[:3] * L, s0_nd[3:] * n * L])

    # Dimensional ROTATING-frame EOM (same frame, SI-scaled): primaries FIXED at
    # (-mu*L, 0), ((1-mu)*L, 0); Coriolis 2n zhat x v; centrifugal +n^2(x,y,0).
    def f_dim_full(t_s, s):
        x, y, z, vx, vy, vz = s
        dx1, dy1 = x + mu * L, y
        dx2, dy2 = x - (1.0 - mu) * L, y
        r1 = math.sqrt(dx1 * dx1 + dy1 * dy1 + z * z)
        r2 = math.sqrt(dx2 * dx2 + dy2 * dy2 + z * z)
        gx = -GM_EARTH_KM3S2 * dx1 / r1**3 - GM_MOON_KM3S2 * dx2 / r2**3
        gy = -GM_EARTH_KM3S2 * dy1 / r1**3 - GM_MOON_KM3S2 * dy2 / r2**3
        gz = -GM_EARTH_KM3S2 * z / r1**3 - GM_MOON_KM3S2 * z / r2**3
        return np.array(
            [vx, vy, vz, 2.0 * n * vy + gx + n * n * x, -2.0 * n * vx + gy + n * n * y, gz]
        )

    tg_dim = tg_nd / n
    tr_dim = rk4_propagate(f_dim_full, tg_dim, s0_dim)
    final_err_pos = np.linalg.norm(tr_nd[-1, :3] * L - tr_dim[-1, :3])
    final_err_vel = np.linalg.norm(tr_nd[-1, 3:] * (n * L) - tr_dim[-1, 3:])
    v_scale = n * L
    out["trajectory_90day"] = {
        "final_pos_err_km": float(final_err_pos),
        "final_pos_rel": float(final_err_pos / L),
        "final_vel_rel": float(final_err_vel / v_scale),
    }

    # Jacobi along the whole mapped history
    worst = 0.0
    for s_nd, s_dim in zip(tr_nd[::40], tr_dim[::40]):
        c_nd = jacobi_constant(s_nd, mu)
        c_dm = c_dim_axis(s_dim[0], s_dim[1], s_dim[3], s_dim[4])
        worst = max(worst, abs(c_dm / scale - c_nd) / abs(c_nd))
    out["jacobi_along_path_worst_rel"] = worst
    out["equilibria"] = dim_eq
    return out


def study_zvc():
    """Zero-velocity-curve neck test at L2 with quantitative gates [Track D]."""
    mu = MU_EM
    f = rhs_rotating(mu)
    pts = lagrange_points(mu)
    x_L2 = float(pts["L2"][0])
    c2 = jacobi_constant(np.concatenate([pts["L2"], np.zeros(3)]), mu)
    delta = 1e-3
    out = {"c_L2": c2, "delta": delta, "x_L2": x_L2}

    # closed case: release at rest just inside the wall
    wall = zvc_wall_position(c2 + delta, mu)
    out["wall_position"] = wall
    s0 = np.array([wall - 1e-6, 0.0, 0.0, 0.0, 0.0, 0.0])
    tg = np.linspace(0.0, 150.0, 150001)
    traj, t_used, drift = propagate_with_guard(f, tg, s0, mu, TOLERANCES["zvc_drift_guard"])
    out["closed_case"] = {
        "max_x": float(np.max(traj[:, 0])),
        "never_past_wall": bool(np.max(traj[:, 0]) < wall),
        "never_past_L2": bool(np.max(traj[:, 0]) < x_L2),
        "cj_drift": drift,
    }

    # open case: launch at x=1.10 with speed from C = C2 - delta
    x_start = 1.10
    v_open = critical_escape_speed(x_start, c2 - delta, mu)
    s0 = np.array([x_start, 0.0, 0.0, v_open, 0.0, 0.0])
    tg = np.linspace(0.0, 60.0, 60001)
    traj, t_used, drift = propagate_with_guard(f, tg, s0, mu, TOLERANCES["zvc_drift_guard"])
    crossed_idx = np.argmax(traj[:, 0] > x_L2 + 1e-2) if np.any(traj[:, 0] > x_L2 + 1e-2) else -1
    out["open_case"] = {
        "launch_speed": v_open,
        "crossed_decisively": bool(crossed_idx >= 0),
        "crossing_time": float(t_used[crossed_idx]) if crossed_idx >= 0 else None,
        "cj_drift": drift,
    }

    # Escape sweep vs launch speed [redesigned per Track D fidelity findings]:
    # escape through the open neck is NECESSARILY confined to v >= v_crit
    # (energy side), but it is not sufficient on a fixed horizon -- temporary
    # capture delays crossings for a band of speeds, and near-threshold grazing
    # runs corrupt C at coarse dt. We therefore measure: (i) every crossing run
    # must satisfy v >= v_crit*(1 - 1e-2) after drift guarding (necessary-
    # condition direction); (ii) crossing time decreases as v grows beyond the
    # capture band. The old "critical speed bisection" premise (monotone
    # escape-in-horizon) was falsified by direct probing and removed.
    v_crit_analytic = critical_escape_speed(x_start, c2, mu)
    sweep = {}
    for k_rel in (1.02, 1.05, 1.10, 1.20, 1.35):
        v = v_crit_analytic * k_rel
        s0v = np.array([x_start, 0.0, 0.0, v, 0.0, 0.0])
        tg_v = np.linspace(0.0, 60.0, 120001)
        tr, _, dr = propagate_with_guard(f, tg_v, s0v, mu, TOLERANCES["zvc_drift_guard"])
        crossed = np.any(tr[:, 0] > x_L2 + 1e-2)
        t_cross = float(tg_v[int(np.argmax(tr[:, 0] > x_L2 + 1e-2))]) if crossed else None
        sweep[f"v_{k_rel:g}x_vcrit"] = {
            "v_over_vcrit": k_rel,
            "escaped_within_60tu": bool(crossed),
            "t_cross": t_cross,
            "cj_drift": dr,
            "trusted": bool(dr <= TOLERANCES["zvc_drift_guard"]),
        }
    escaping_trusted = [
        s["v_over_vcrit"] for s in sweep.values() if s["escaped_within_60tu"] and s["trusted"]
    ]
    sweep["analytic_v_crit"] = v_crit_analytic
    sweep["min_trusted_escape_v_over_vcrit"] = min(escaping_trusted) if escaping_trusted else None
    sweep["necessary_condition_honored"] = bool(
        escaping_trusted and min(escaping_trusted) >= 1.0 - 1e-2
    )
    out["escape_sweep"] = sweep
    return out


def study_inertial_consistency_and_mutants():
    """The centerpiece discriminator: rotating-vs-inertial propagation agreement,
    plus the pre-registered mutant battery (each mutant must be caught)."""
    mu = MU_EM
    span = 1.0
    n_steps = 1000
    t_grid = np.linspace(0.0, span, n_steps + 1)
    pts = lagrange_points(mu)

    # start near L1 displaced along +x with +y velocity (encounter-free)
    s0_rot = np.array([float(pts["L1"][0]) - 0.02, 0.05, 0.01, 0.03, 0.11, -0.02])
    tr_rot = rk4_propagate(rhs_rotating(mu), t_grid, s0_rot)
    s_end_via_rot = rot_to_inert(tr_rot[-1], span)

    tr_in = rk4_propagate(rhs_inertial(mu), t_grid, rot_to_inert(s0_rot, 0.0))
    scale = max(np.max(np.abs(s_end_via_rot[:3])), 1.0)
    clean_resid = float(np.linalg.norm(s_end_via_rot[:3] - tr_in[-1, :3]) / scale)

    results = {"clean_relative_residual": clean_resid, "horizon": span, "n_steps": n_steps}
    mutants = {}

    # coriolis sign flip
    def f_cor(t, s):
        g = grad_omega(s[0], s[1], s[2], mu)
        return np.array([s[3], s[4], s[5], -2.0 * s[4] + g[0], 2.0 * s[3] + g[1], g[2]])

    tr_mut = rk4_propagate(f_cor, t_grid, s0_rot)
    mutants["coriolis_flip"] = float(
        np.linalg.norm(rot_to_inert(tr_mut[-1], span)[:3] - tr_in[-1, :3]) / scale
    )

    # centrifugal drop (gravity-only gradient)
    def g_grav(x, y, z, mu):
        r1 = np.sqrt((x + mu) ** 2 + y * y + z * z)
        r2 = np.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
        return np.array(
            [
                -(1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3,
                -(1.0 - mu) * y / r1**3 - mu * y / r2**3,
                -(1.0 - mu) * z / r1**3 - mu * z / r2**3,
            ]
        )

    def f_cent(t, s):
        g = g_grav(s[0], s[1], s[2], mu)
        return np.array([s[3], s[4], s[5], 2.0 * s[4] + g[0], -2.0 * s[3] + g[1], g[2]])

    tr_mut = rk4_propagate(f_cent, t_grid, s0_rot)
    mutants["centrifugal_drop"] = float(
        np.linalg.norm(rot_to_inert(tr_mut[-1], span)[:3] - tr_in[-1, :3]) / scale
    )

    # rotation sense flipped in the MAPPING only (R(-t) applied to history)
    s_end_badmap = rot_to_inert(tr_rot[-1], -span)
    mutants["mapping_sign_flip"] = float(np.linalg.norm(s_end_badmap[:3] - tr_in[-1, :3]) / scale)

    # mu convention error (m2/m1 instead of m2/(m1+m2)): caught by the MU-firewall
    # (bit-exact GM re-derivation in tests) and the exact L4-x identity, NOT by
    # percent-level anchors [Track E finding]. Signal = relative gamma_L1 shift.
    mu_wrong = MU_EM / (1.0 - MU_EM)
    pts_wrong = lagrange_points(mu_wrong)
    gamma_true = 1.0 - mu - float(pts["L1"][0])
    gamma_wrong = 1.0 - mu_wrong - float(pts_wrong["L1"][0])
    mutants["mu_convention_gamma_rel_shift"] = abs(gamma_wrong - gamma_true) / gamma_true

    # body swap (treat m2 as m1 => mu' = 1-mu outside the frozen domain (0, 1/2]):
    # the proven-bracket construction must REJECT it (Track B midpoint sign flip).
    try:
        lagrange_points(1.0 - mu)
        rejected = False
    except Exception:
        rejected = True
    mutants["body_swap_domain_rejected"] = 1.0 if rejected else 0.0

    # planar-frozen Jacobi evaluator leak (spatial orbit, evaluator ignoring z)
    s0_sp = np.array([float(pts["L4"][0]), float(pts["L4"][1]), 0.05, 0.0, 0.0, 0.02])
    tg = np.linspace(0.0, 20.0, 20001)
    tr_sp = rk4_propagate(rhs_rotating(mu), tg, s0_sp)
    c0_true = jacobi_constant(s0_sp, mu)

    def jacobi_frozen(state, mu_):
        xy_pot = 2.0 * (
            (1.0 - mu_) / np.linalg.norm(state[:2] - np.array([-mu_, 0.0]))
            + mu_ / np.linalg.norm(state[:2] - np.array([1.0 - mu_, 0.0]))
            + 0.5 * 0.0  # frozen-z potential omits the true 3-D distances
        )
        return xy_pot - float(np.dot(state[3:], state[3:]))

    leak = max(abs(jacobi_frozen(s, mu) - jacobi_frozen(s0_sp, mu)) for s in tr_sp)
    true_drift = max(abs(jacobi_constant(s, mu) - c0_true) for s in tr_sp)
    mutants["planar_frozen_jacobi"] = float(leak)
    results["spatial_true_drift"] = float(true_drift)

    # spectra blindness to coriolis flip (documents why T-IC is required)
    A_clean = stability_matrix(float(pts["L1"][0]), 0.0, mu)
    Wxx, Wyy, Wxy, _ = hessian_entries(float(pts["L1"][0]), 0.0, mu)
    A_flip = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [Wxx, Wxy, 0, -2.0], [Wxy, Wyy, 2.0, 0.0]])
    ev_clean = np.sort_complex(np.linalg.eigvals(A_clean))
    ev_flip = np.sort_complex(np.linalg.eigvals(A_flip))
    results["coriolis_flip_max_eigen_shift"] = float(np.max(np.abs(ev_clean - ev_flip)))

    results["mutant_signals"] = mutants
    # per-mutant discriminator thresholds (each mutant must exceed its bound)
    thresholds = {
        "coriolis_flip": TOLERANCES["mutant_signal_min"],
        "centrifugal_drop": 1e-6,  # equilibrium residual blows past the 2e-14 gate
        "mapping_sign_flip": TOLERANCES["mutant_signal_min"],
        "mu_convention_gamma_rel_shift": 1e-12,  # firewall-grade exactness bound
        "body_swap_domain_rejected": 0.5,
        "planar_frozen_jacobi": TOLERANCES["zvc_drift_guard"],
    }
    caught = {name: bool(mutants[name] > thr) for name, thr in thresholds.items()}
    results["mutant_thresholds"] = thresholds
    results["mutant_caught"] = caught
    results["all_mutants_caught"] = bool(all(caught.values()))
    return results


def study_symmetries_and_limits():
    """Mirror/rotation symmetries + singular-limit asymptotics [Tracks D/B/E]."""
    out = {}
    mu = MU_EM
    rng = 0.0
    # mirror field law f(Mw) = -M f(w), M = diag(1,-1,1,-1,1,-1)
    worst = 0.0
    fr = rhs_rotating(mu)
    for x in np.linspace(-1.5, 1.8, 21):
        for y in np.linspace(-1.2, 1.2, 17):
            for vx in (-0.3, 0.0, 0.4):
                s = np.array([x, y, 0.1, vx, 0.2, 0.05])
                Ms = np.array([x, -y, 0.1, -vx, 0.2, -0.05])
                f_s = fr(rng, s)
                m_f_s = np.array([f_s[0], -f_s[1], f_s[2], -f_s[3], f_s[4], -f_s[5]])
                viol = np.max(np.abs(fr(rng, Ms) + m_f_s))
                worst = max(worst, viol)
    out["mirror_field_law_max_violation"] = float(worst)

    # equal-mass pi-rotation equivariance f(Pw) = P f(w)
    fr_half = rhs_rotating(0.5)
    worst = 0.0
    for x in np.linspace(-1.5, 1.5, 19):
        for y in np.linspace(-1.0, 1.0, 15):
            s = np.array([x, y, 0.05, 0.2, -0.1, 0.03])
            Ps = np.array([-x, -y, 0.05, -0.2, 0.1, 0.03])
            fP = fr_half(rng, Ps)
            Pf = np.array([-fr_half(rng, s)[0], -fr_half(rng, s)[1], fr_half(rng, s)[2], -fr_half(rng, s)[3], -fr_half(rng, s)[4], fr_half(rng, s)[5]])
            worst = max(worst, np.max(np.abs(fP - Pf)))
    out["half_mass_pi_rotation_max_violation"] = float(worst)

    # singular limits: Hill scaling ratios and L3 offset law
    lim = {}
    for mu_small in (1e-3, 1e-6):
        pts = lagrange_points(mu_small)
        alpha = (mu_small / 3.0) ** (1.0 / 3.0)
        g1 = (1.0 - mu_small) - float(pts["L1"][0])
        g2 = float(pts["L2"][0]) - (1.0 - mu_small)
        series1 = alpha * (1.0 - alpha / 3.0 - alpha**2 / 9.0)
        series2 = alpha * (1.0 + alpha / 3.0 - alpha**2 / 9.0)
        lim[f"mu_{mu_small:g}"] = {
            "g1_over_alpha": g1 / alpha,
            "g1_over_series": g1 / series1,
            "g2_over_series": g2 / series2,
            "xL3_plus_1_over_mu": (float(pts["L3"][0]) + 1.0) / mu_small,
            "xL3_offset_expected": -5.0 / 12.0,
        }
    out["singular_limits"] = lim

    # log-log slope dlog(gamma1)/dlog(mu) ~ 1/3 across decades
    mus = np.logspace(-6, -2, 9)
    g1s = []
    for m_ in mus:
        p = solve_collinear(float(m_))
        g1s.append((1.0 - float(m_)) - p["L1"])
    slopes = convergence_rate(np.array(g1s), mus)
    out["hill_slope_convergence"] = {
        "mus": mus.tolist(),
        "gamma1": g1s,
        "pairwise_orders": slopes.tolist(),
    }

    # equal-mass symmetry: x_L1 = 0 exactly, x_L3 = -x_L2, L4 = (0, sqrt3/2)
    half = CASES["equal_mass"]
    pts_h = lagrange_points(half)
    out["equal_mass"] = {
        "x_L1_exact_zero": float(pts_h["L1"][0]),
        "sym_x2_plus_x3": float(pts_h["L2"][0] + pts_h["L3"][0]),
        "L4_x": float(pts_h["L4"][0]),
    }
    return out


def study_rk4_order_verification():
    """Closed-form kinematic round trip: free particle seen from the rotating
    frame has exact solution r_R(t) = (t cos t, -t sin t, 0) [order-4 probe]."""
    t_span = 1.5
    s0 = inert_to_rot(np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]), 0.0)

    def f_kin(t, s):
        return np.array([s[3], s[4], s[5], 2.0 * s[4] + s[0], -2.0 * s[3] + s[1], 0.0])

    errors = {}
    for n in (125, 250, 500, 1000):
        tg = np.linspace(0.0, t_span, n + 1)
        tr = rk4_propagate(f_kin, tg, s0)
        exact = np.array(
            [
                t_span * math.cos(t_span),
                -t_span * math.sin(t_span),
                0.0,
                math.cos(t_span) - t_span * math.sin(t_span),
                -math.sin(t_span) - t_span * math.cos(t_span),
                0.0,
            ]
        )
        errors[n] = float(np.linalg.norm(rot_to_inert(tr[-1], t_span) - np.array([t_span, 0.0, 0.0, 1.0, 0.0, 0.0])))
    ns_arr = np.array(sorted(errors.keys()), dtype=float)
    errs_arr = np.array([errors[int(k)] for k in sorted(errors.keys())])
    orders = convergence_rate(errs_arr, t_span / ns_arr)
    return {
        "errors_by_n": {str(int(k)): v for k, v in errors.items()},
        "pairwise_orders": orders.tolist(),
        "expected_order": 4,
    }


def study_bridge_and_drift_law():
    """Bridge identity C = 2(hz - E_I), E_I+C/2=L_z, and the E_I drift law."""
    mu = MU_EM
    s0, _ = long_period_ic(1e-2, mu)
    span = 63.0  # ~3 long periods
    out = {}
    for n in (4096, 16384):
        tg = np.linspace(0.0, span, n + 1)
        tr = rk4_propagate(rhs_rotating(mu), tg, s0)
        bridge, ident = inertial_bridge_drift(tr, tg, mu)
        e_hist, drift_integral = inertial_energy_history(tr, tg, mu)
        e_drift_measured = float(e_hist[-1] - e_hist[0])
        rel_mismatch = abs(e_drift_measured - drift_integral) / max(abs(drift_integral), 1e-300)
        c_drift = max(abs(jacobi_constant(s, mu) - jacobi_constant(s0, mu)) for s in tr)
        out[f"n_{n}"] = {
            "bridge_identity_max": bridge,
            "identity_E_plus_C2_minus_Lz_max": ident,
            "c_drift": c_drift,
            "EI_drift_measured": e_drift_measured,
            "EI_drift_integral_predicted": drift_integral,
            "drift_law_rel_mismatch": rel_mismatch,
        }
    ratio = out["n_4096"]["drift_law_rel_mismatch"] / max(out["n_16384"]["drift_law_rel_mismatch"], 1e-300)
    out["drift_law_mismatch_reduction_ratio"] = ratio
    return out


# --------------------------------------------------------------------------- #
# Figures (deterministic, generated from recorded data)
# --------------------------------------------------------------------------- #
FIG_DIR = Path(__file__).resolve().parent / "results" / "figures"


def make_figures(eq_study, stab_study, jac_study, zvc_study):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures = []

    # F1: rotating-frame geometry with all five Lagrange points (Earth-Moon)
    fig, ax = plt.subplots(figsize=(9, 6.5))
    mu = MU_EM
    pts = lagrange_points(mu)
    ax.plot([-mu], [0], "o", color="tab:blue", markersize=14, label="Earth ($m_1$)")
    ax.plot([1 - mu], [0], "o", color="gray", markersize=5, label="Moon ($m_2$)")
    labels = {"L1": "L1", "L2": "L2", "L3": "L3", "L4": "L4 (leading)", "L5": "L5 (trailing)"}
    for name, pos in pts.items():
        ax.plot(pos[0], pos[1], marker="*", markersize=13, color="crimson")
        ax.annotate(labels[name], (pos[0], pos[1]), textcoords="offset points", xytext=(6, 4), fontsize=9)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th) - mu, np.sin(th), "--", lw=0.7, color="lightgray")
    ax.plot(np.cos(th) + 1 - mu, np.sin(th), "--", lw=0.7, color="lightgray")
    ax.set_aspect("equal")
    ax.set_xlabel("$x$ [-]")
    ax.set_ylabel("$y$ [-]")
    ax.set_title(f"Exp 011 F1: CR3BP geometry, Earth-Moon $\\mu$={mu:.7f}, barycentric rotating frame")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    p = FIG_DIR / "f1_geometry_em.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    figures.append(p.name)

    # F2: zero-velocity structure at the critical levels
    fig, ax = plt.subplots(figsize=(10, 7))
    xs = np.linspace(-1.6, 1.9, 700)
    ys = np.linspace(-1.35, 1.35, 520)
    XX, YY = np.meshgrid(xs, ys)
    with np.errstate(all="ignore"):
        ZZ = 2.0 * omega_eff(XX, YY, 0.0, mu)
    cv = eq_study["cases"]["earth_moon"]["c_values"]
    # increasing order required by contour; physics ordering is C1>C2>C3>C4=C5
    levels = [cv["L4"], cv["L3"], cv["L2"], cv["L1"]]
    cs = ax.contourf(XX, YY, ZZ, levels=[levels[0], 30], colors=["white", "#ffd9d9"], hatches=["", "//"])
    ax.contour(XX, YY, ZZ, levels=levels, colors=["C0", "C1", "C2", "C3"], linewidths=1.0)
    for name, pos in pts.items():
        ax.plot(pos[0], pos[1], marker="*", markersize=11, color="black")
    ax.plot([-mu], [0], "o", color="tab:blue", markersize=10)
    ax.plot([1 - mu], [0], "o", color="gray", markersize=4)
    ax.set_aspect("equal")
    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(ys[0], ys[-1])
    ax.set_title("Exp 011 F2: zero-velocity curves at critical Jacobi levels C3<C4=C5<C2<C1 (EM), forbidden region hatched")
    ax.set_xlabel("$x$ [-]")
    ax.set_ylabel("$y$ [-]")
    fig.tight_layout()
    p = FIG_DIR / "f2_zvc_levels_em.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    figures.append(p.name)

    # F3: eigenvalue spectrum in the complex plane (selected cases)
    fig, ax = plt.subplots(figsize=(9, 6))
    markers = {"L1": "o", "L2": "s", "L3": "^", "L4": "*"}
    for name in ("L1", "L2", "L3", "L4"):
        pos = lagrange_points(mu)[name]
        lam, _, _, _, _ = eigen_analysis(float(pos[0]), float(pos[1]), mu)
        ax.scatter(lam.real, lam.imag, marker=markers[name], s=42, facecolors="none", edgecolors=f"C{ord(name[1])-49}", label=f"{name} ($\\mu$={MU_EM:.5f})")
    lam05, _, _, _, _ = eigen_analysis(0.45, SQRT3_2, 0.05)
    ax.scatter(lam05.real, lam05.imag, marker="D", s=38, facecolors="none", edgecolors="black", label="L4 ($\\mu$=0.05, unstable quartet)")
    ax.axhline(0, color="gray", lw=0.6)
    ax.axvline(0, color="gray", lw=0.6)
    ax.set_xlabel("Re $\\lambda$")
    ax.set_ylabel("Im $\\lambda$")
    ax.set_title("Exp 011 F3: linear-stability spectra at the equilibria (EM system + above-threshold L4)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    p = FIG_DIR / "f3_eigen_spectrum.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    figures.append(p.name)

    # F4: Jacobi drift ladders + order guide + floor line
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))
    cap = jac_study["class_capture_C2minus1e-3"]
    nsB = sorted(int(k) for k in cap["drifts"])
    dB = np.array([cap["drifts"][n] for n in nsB])
    dtB = cap["span"] / np.array(nsB, dtype=float)
    ax1.loglog(dtB, dB, "o-", label="capture class B")
    guide = dB[0] * (dtB / dtB[0]) ** 4
    ax1.loglog(dtB, guide, "k:", label="$\\propto h^4$ guide")
    ret = jac_study["class_retrograde_rho0p30"]
    nsC = sorted(int(k) for k in ret["drifts"])
    dC = np.array([ret["drifts"][n] for n in nsC])
    dtC = ret["span"] / np.array(nsC, dtype=float)
    ax1.loglog(dtC, dC, "s-", label="retrograde class C")
    fl = ret["floor_probe"]
    nsF = sorted(int(k) for k in fl)
    dF = np.array([fl[n] for n in nsF])
    dtF = ret["span"] / np.array(nsF, dtype=float)
    ax1.loglog(dtF, dF, "^-", label="floor probe (turn-up)")
    ax1.set_xlabel("$h$ [-]")
    ax1.set_ylabel("max $|\\Delta C|$ [-]")
    ax1.set_title("Exp 011 F4a: Jacobi drift convergence (RK4, EM $\\mu$)")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.25, which="both")

    lp = jac_study["class_L4_LP_amp3e-3"]
    nsA = sorted(int(k) for k in lp["drifts"])
    dA = np.array([lp["drifts"][n] for n in nsA])
    dtA = lp["span"] / np.array(nsA, dtype=float)
    ax2.loglog(dtA, dA, "o-", label="L4 LP mode, $\\varepsilon$=1e-2")
    ax2.axhline(8.88e-16, color="red", ls="--", lw=0.8, label="2-ulp quantization floor")
    ax2.set_xlabel("$h$ [-]")
    ax2.set_ylabel("max $|\\Delta C|$ [-]")
    ax2.set_title("Exp 011 F4b: bounded-orbit drift reaches roundoff floor")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.25, which="both")
    fig.tight_layout()
    p = FIG_DIR / "f4_jacobi_convergence.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    figures.append(p.name)
    return figures


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> dict:
    print("[011] equilibria ...", flush=True)
    eq = study_equilibria()
    print("[011] stability ...", flush=True)
    stab = study_stability()
    print("[011] jacobi drift ...", flush=True)
    jac, *_ = study_jacobi_drift()
    print("[011] perturbations ...", flush=True)
    pert = study_perturbations()
    print("[011] dimensional cross-check ...", flush=True)
    dim = study_dimensional_crosscheck()
    print("[011] ZVC ...", flush=True)
    zvc = study_zvc()
    print("[011] inertial consistency + mutants ...", flush=True)
    mut = study_inertial_consistency_and_mutants()
    print("[011] symmetries + limits ...", flush=True)
    sym = study_symmetries_and_limits()
    print("[011] RK4 order verification ...", flush=True)
    order = study_rk4_order_verification()
    print("[011] bridge + drift law ...", flush=True)
    brid = study_bridge_and_drift_law()

    figures = make_figures(eq, stab, jac, zvc)

    headline = {
        "mu_EM": MU_EM,
        "x_L1": eq["cases"]["earth_moon"]["L1"]["x"],
        "x_L2": eq["cases"]["earth_moon"]["L2"]["x"],
        "x_L3": eq["cases"]["earth_moon"]["L3"]["x"],
        "em_L1_km_from_moon": eq["anchors"]["em_L1_km_from_moon"],
        "em_L2_km_from_moon": eq["anchors"]["em_L2_km_from_moon"],
        "sem_L1_km_from_earth": eq["anchors"]["sem_L1_km_from_earth"],
        "sem_L2_km_from_earth": eq["anchors"]["sem_L2_km_from_earth"],
        "routh_threshold": MU_ROUTH,
        "routh_boundary_gamma_minus_quarter": stab["boundary_gamma_minus_quarter_50dps"],
        "max_equilibrium_residual": max(
            eq["cases"][c][n]["residual_inf"] for c in CASES for n in ("L1", "L2", "L3", "L4", "L5")
        ),
        "inertial_consistency_clean_residual": mut["clean_relative_residual"],
        "mutant_catch_matrix": mut["mutant_caught"],
        "all_mutants_caught": mut["all_mutants_caught"],
        "capture_class_orders": jac["class_capture_C2minus1e-3"]["orders"],
        "retrograde_floor_turnup": jac["class_retrograde_rho0p30"]["floor_probe"],
        "lp_floor_plateau": jac["class_L4_LP_floor_plateau"],
        "L1_growth_rate_rel_err_at_1e-4": pert["L1_unstable_growth"]["eps_0.0001"]["rel_error"],
        "L4_LP_freq_rel_err_at_1e-4": pert["L4_longperiod_mode"]["amp_0.0001"]["rel_error"],
        "dim_traj_final_pos_rel": dim["trajectory_90day"]["final_pos_rel"],
        "dim_jacobi_along_path_worst_rel": dim["jacobi_along_path_worst_rel"],
        "zvc_closed_never_past_wall": zvc["closed_case"]["never_past_wall"],
        "zvc_open_crossed": zvc["open_case"]["crossed_decisively"],
        "zvc_necessary_condition_honored": zvc["escape_sweep"]["necessary_condition_honored"],
        "kinematic_roundtrip_orders": order["pairwise_orders"],
        "mirror_field_law_max_violation": sym["mirror_field_law_max_violation"],
    }

    results = {
        "constants": {
            "GM_SUN_KM3S2": GM_SUN_KM3S2,
            "GM_EARTH_KM3S2": GM_EARTH_KM3S2,
            "GM_MOON_KM3S2": GM_MOON_KM3S2,
            "GM_provenance": "IAU 2015 nominal GM_E/GM_Sun (Exp 008 canon); DE-series GM_Moon",
            "MU_EM": MU_EM,
            "MU_SEM": MU_SEM,
            "MU_ROUTH": MU_ROUTH,
            "MU_derivation": "mu = GM2/(GM1+GM2) from declared GMs (portfolio-consistent)",
            "MU_LITERATURE_NOTE": MU_LITERATURE_NOTE,
            "A_EM_KM": A_EM_KM,
            "AU_KM": AU_KM,
            "N_EM_RAD_S": N_EM_RAD_S,
            "T_EM_DAYS": T_EM_DAYS,
            "frame_convention": FRAME_CONVENTION,
            "jacobi_sign_convention": "C = 2*omega_eff - v_rot^2 ; C_dim = n^2 a^2 C",
        },
        "headline": headline,
        "equilibria": eq,
        "stability": stab,
        "jacobi_conservation": jac,
        "perturbation_validation": pert,
        "dimensional_crosscheck": dim,
        "zero_velocity_structure": zvc,
        "adversarial": mut,
        "symmetries_limits": sym,
        "rk4_order_verification": order,
        "bridge_drift_law": brid,
        "tolerances": TOLERANCES,
        "limitations": [
            "Classical circular CR3BP: primaries on exact circular rails; real Moon orbit is eccentric (e=0.055) and solar-perturbed.",
            "No ephemeris fidelity: anchors validate scale+geometry at percent level only.",
            "Collinear equilibria are linearly unstable: station-keeping costs and halo/Lyapunov orbit families are OUT of scope (future experiment).",
            "Fixed-step RK4: grazing trajectories (pericenter << 1e-2) need adaptive stepping or guards; ZVC runs carry a Jacobi-drift guard.",
            "Linear stability does not by itself prove nonlinear stability; KAM refinements (Arnold 1963; Deprit & Deprit-Bartholome 1967) cited, not re-derived.",
            "Triangular-point stability claim restricted to mu < mu_Routh excluding the degenerate boundary (1:1 resonance).",
        ],
        "figures": figures,
        "figures_note": "Deterministic matplotlib/Agg, dpi=150, generated from recorded result data.",
    }

    out_path = save_json_result(
        Path(__file__).resolve().parent / "results" / "results.json",
        results,
        name="lagrange_points_cr3bp",
        description="Experiment 011: CR3BP Lagrange points - equilibria, Jacobi integral, stability, adversarial validation",
    )
    print(f"[011] results -> {out_path}")
    return results


if __name__ == "__main__":
    main()
