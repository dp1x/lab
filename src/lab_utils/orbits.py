"""Shared orbital-mechanics canon: constants, Kepler solver, element machinery.

Canonical copies of machinery that already has multiple proven consumers chained
through per-experiment importlib hops (Exp 008/009/010 -> {008, 009, 006, 002}).
Bodies are transcribed verbatim from the donor experiments; experiments 001-010
retain frozen local versions pending opportunistic migration, and equivalence is
pinned by tests in ``src/lab_utils/tests/test_orbits_canon.py``.
``j2_rhs`` (J2 Cowell right-hand side) graduated at Exp 012 -- second consumer
after donor Exp 009; pinned against the donor propagator on J2-on and J2-off paths.

Conventions: ECI frame, km / km^3-s^2 / s units, angles in radians.
"""

from __future__ import annotations

import numpy as np

__all__ = [
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

# --------------------------------------------------------------------------- #
# Physical constants (provenance as pinned by Exp 008/009)
# --------------------------------------------------------------------------- #
MU_EARTH_KM3S2 = 398600.4418  # IAU 2015 nominal GM_E (km^3/s^2); JPL DE440 planet-only 398600.435507 differs 1.5e-8 rel
R_EARTH_KM = 6378.137  # WGS-84 equatorial radius (km)
OMEGA_EARTH_RAD_S = 7.2921159e-5  # WGS-84 / Vallado Table 3-1 (rad/s)
J2_EARTH = 1.082629821e-3  # WGS-84 J2 = sqrt(5)|C20_bar| (provenance Exp 009)

NODE_GUARD_REL = 1e-6  # |zhat x h| / |h| = sin(i) below this -> RAAN undefined
ECC_GUARD_ABS = 1e-8  # |e_vec| below this -> argument of periapsis undefined


# --------------------------------------------------------------------------- #
# Kepler solver and element helpers (donor: Exp 008 groundtracks)
# --------------------------------------------------------------------------- #
def solve_kepler(M: np.ndarray | float, e: float, tol: float = 1e-14, max_iter: int = 100) -> np.ndarray | float:
    """Newton solve M = E - e sin E for E. Vectorized, seed E0 = M + e sin M."""
    M_arr = np.asarray(M, dtype=float)
    scalar = M_arr.ndim == 0
    M_arr = np.atleast_1d(M_arr)
    E = M_arr + e * np.sin(M_arr)
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M_arr
        if np.max(np.abs(f)) < tol:
            break
        denom = 1.0 - e * np.cos(E)
        E = E - f / denom
    else:
        f = E - e * np.sin(E) - M_arr
        if np.max(np.abs(f)) >= tol:
            raise ValueError(f"Kepler did not converge e={e} max|f|={np.max(np.abs(f)):.2e}")
    if scalar:
        return float(E[0])
    return E


def true_anomaly_from_E(E: np.ndarray | float, e: float) -> np.ndarray | float:
    """Convert eccentric anomaly E to true anomaly nu (atan2, preserves quadrant)."""
    E_arr = np.asarray(E, dtype=float)
    scalar = E_arr.ndim == 0
    E_arr = np.atleast_1d(E_arr)
    cos_nu = (np.cos(E_arr) - e) / (1.0 - e * np.cos(E_arr))
    sin_nu = np.sqrt(1.0 - e * e) * np.sin(E_arr) / (1.0 - e * np.cos(E_arr))
    cos_nu = np.clip(cos_nu, -1.0, 1.0)
    nu = np.arctan2(sin_nu, cos_nu)
    if scalar:
        return float(nu[0])
    return nu


def orbital_period(a: float, mu: float = MU_EARTH_KM3S2) -> float:
    """Kepler period T = 2*pi*sqrt(a^3/mu)."""
    return 2.0 * np.pi * np.sqrt(a ** 3 / mu)


def mean_motion(a: float, mu: float = MU_EARTH_KM3S2) -> float:
    """Mean motion n = sqrt(mu/a^3) [rad/s]."""
    return np.sqrt(mu / a ** 3)


def rotation_matrix_313(Omega: float, inc: float, omega: float) -> np.ndarray:
    """Q = R_z(Omega)*R_x(i)*R_z(omega) (3-1-3). All rad. Returns 3x3."""
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(omega), np.sin(omega)
    Q = np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si, cw * si, ci],
    ], dtype=float)
    return Q


def coe_to_rv_eci(a: float, e: float, inc: float, Omega: float, omega: float, nu: float,
                  mu: float = MU_EARTH_KM3S2) -> tuple[np.ndarray, np.ndarray]:
    """Single-epoch perifocal -> ECI for state (r_eci, v_eci) at true anomaly nu (rad)."""
    p = a * (1.0 - e * e)
    r = p / (1.0 + e * np.cos(nu)) if abs(e) < 1.0 else a  # e==1 not used
    r_pf = np.array([r * np.cos(nu), r * np.sin(nu), 0.0], dtype=float)
    h = np.sqrt(mu * p)
    v_pf = np.array([-(mu / h) * np.sin(nu), (mu / h) * (e + np.cos(nu)), 0.0], dtype=float)
    Q = rotation_matrix_313(Omega, inc, omega)
    return Q @ r_pf, Q @ v_pf


# --------------------------------------------------------------------------- #
# State seeding + osculating element recovery (donor: Exp 009 j2Precession)
# --------------------------------------------------------------------------- #
def seed_state(a, e, inc, Om, om, M0, mu=MU_EARTH_KM3S2):
    """Classical elements + M0 -> epoch state (r0, v0, nu0). M0 -> E0 -> nu0 chain."""
    E0 = solve_kepler(np.mod(M0, 2.0 * np.pi), e)
    nu0 = true_anomaly_from_E(E0, e)
    r0, v0 = coe_to_rv_eci(a, e, inc, Om, om, float(nu0), mu)
    return r0, v0, float(nu0)


def rv_to_coe_eci(r, v, mu=MU_EARTH_KM3S2) -> dict:
    """Cartesian state(s) -> classical osculating elements with singular guards.

    h = r x v; node vector n = zhat x h; e_vec = (v x h)/mu - r/r.
    Guards: sin(i) < NODE_GUARD_REL -> Omega NaN (RAAN undefined);
    |e_vec| < ECC_GUARD_ABS -> omega NaN (periapsis direction undefined).
    Accepts (N,3)/(N,3) or single vectors; returns arrays of a,e,i,Omega,omega,nu.
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    single = r.ndim == 1
    if single:
        r = r[None, :]
        v = v[None, :]
    rn = np.linalg.norm(r, axis=1)
    hvec = np.cross(r, v)
    hn = np.linalg.norm(hvec, axis=1)
    nvec = np.column_stack([-hvec[:, 1], hvec[:, 0], np.zeros(hn.shape[0])])
    nn = np.linalg.norm(nvec, axis=1)
    evec = np.cross(v, hvec) / mu - r / rn[:, None]
    em = np.linalg.norm(evec, axis=1)
    energy = 0.5 * np.einsum("ij,ij->i", v, v) - mu / rn
    sma = -mu / (2.0 * energy)
    inc = np.arctan2(nn, hvec[:, 2])

    good_node = nn > NODE_GUARD_REL * hn
    good_ecc = em >= ECC_GUARD_ABS
    Om = np.full(hn.shape, np.nan)
    Om[good_node] = np.arctan2(nvec[good_node, 1], nvec[good_node, 0])
    om = np.full(hn.shape, np.nan)
    both = good_node & good_ecc
    if np.any(both):
        nhat = nvec[both] / nn[both, None]
        ehat = evec[both] / em[both, None]
        hhat = hvec[both] / hn[both, None]
        cosw = np.einsum("ij,ij->i", nhat, ehat)
        sinw = np.einsum("ij,ij->i", np.cross(nhat, ehat), hhat)
        om[both] = np.arctan2(sinw, cosw)
    nu = np.full(hn.shape, np.nan)
    gnu = em > ECC_GUARD_ABS
    if np.any(gnu):
        ehat_n = evec[gnu] / em[gnu, None]
        rhat = r[gnu] / rn[gnu, None]
        hhat_n = hvec[gnu] / hn[gnu, None]
        cosnu = np.einsum("ij,ij->i", ehat_n, rhat)
        sinnu = np.einsum("ij,ij->i", np.cross(ehat_n, rhat), hhat_n)
        nu[gnu] = np.arctan2(sinnu, cosnu)
    out = {
        "a": sma,
        "e": em,
        "inc": inc,
        "Omega": Om,
        "omega": om,
        "nu": nu,
        "energy": energy,
        "h_z": hvec[:, 2],
        "h_mag": hn,
    }
    if single:
        return {k: (float(val) if np.ndim(val) == 0 else val[0]) for k, val in out.items()}
    return out


def steps_per_orbit(e: float, base: int = 512, ecc_base: int = 720) -> int:
    """Documented resolution rule: >=512/orbit circular; eccentric scaled as
    ceil(ecc_base/(1-e)^{3/2}) (periapsis-resolution law, Exp 002/008 precedent)."""
    if e <= 0.0:
        return base
    return max(base, int(np.ceil(ecc_base / (1.0 - e) ** 1.5)))


def j2_rhs(mu: float, j2: float, R_eq_km: float = R_EARTH_KM):
    """Build f(t, x6) for r'' = -mu*r/|r|^3 + a_J2(r) in ECI (Z = spin axis).

    a_J2 = -(3 mu J2 R_eq^2)/(2 r^5) * [x(1-5z^2/r^2), y(1-5z^2/r^2), z(3-5z^2/r^2)],
    the gradient of U_J2 = mu J2 R_eq^2 P2(z/r)/r^3 (sign check pinned by tests).

    Graduated at Exp 012 (second consumer after donor Exp 009). Op order is
    transcribed verbatim from the donor's inline accel closure -- Kepler term
    first, j2 == 0 short-circuit preserved, c*(r_i*f) association preserved --
    so ``rk4_propagate(j2_rhs(...), t, x0)`` is arithmetic-compatible with the
    verified Exp 006/009 loop structure (equivalence pinned vs donors in
    ``src/lab_utils/tests/test_orbits_canon.py``). The RHS takes t explicitly
    (non-autonomous-safe signature required by rk4_propagate).
    """
    re2 = R_eq_km ** 2

    def f(_t: float, x: np.ndarray) -> np.ndarray:
        r = x[:3]
        rm = np.linalg.norm(r)
        a_kep = -mu * r / rm**3
        if j2 == 0.0:
            return np.concatenate([x[3:], a_kep])
        z2r2 = (r[2] * r[2]) / (rm * rm)
        c = -1.5 * j2 * mu * re2 / rm**5
        g = 1.0 - 5.0 * z2r2
        a_j2 = c * np.array([r[0] * g, r[1] * g, r[2] * (3.0 - 5.0 * z2r2)])
        return np.concatenate([x[3:], a_kep + a_j2])

    return f
