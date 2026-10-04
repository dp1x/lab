"""Instantaneous nodal-rate instrument and the fixed-frequency secular estimator.

This is the measurement apparatus for `mission_secular_identification`.

Why an instantaneous rate rather than a slope of the node angle
----------------------------------------------------------------
Every previously committed lunisolar number in this Lab is an ordinary
least-squares slope of an unwrapped node angle sampled at equator crossings.
That estimator has two structural defects, both quantified by the 2026-10-04
session (`localdocs/reports/mission-mean-element-leakage-2026-10-04.md`):

* it requires node detection and unwrapping, both of which alias near i = 90
  deg where the node line is ill-conditioned (1/cos^2 i -> infinity), and
* its window is the whole arc, so any forced term at a period comparable to
  the arc length leaks coherently into the slope.

The instantaneous rate

    Omega_dot = ( h_y * (-hdot_x) + h_x * hdot_y ) / (h_x^2 + h_y^2),
    h = r x v,  hdot = r x a

needs no angle, no unwrapping and no window. It is exact for a Newtonian
field because d/dt(r x v) = r x a. The remaining question is how to average it
into a *secular* rate without reintroducing a window bias; that is what
:func:`fit_secular_fixed_frequency` does.

The estimator here is deliberately the *continuous/state-based* instrument the
2026-10-04 handoff asked for. The discrete-angle OLS estimator is retained
only as a contrast in the results, never as a measurement.
"""
from __future__ import annotations

import math

import numpy as np

SEC_PER_DAY = 86400.0
DEG_PER_RAD = 180.0 / math.pi

#: Physical forcing periods relevant to the RAAN of a LEO satellite (days).
#: These are FIXED by nature, not fitted. LAB_CONSTITUTION.md section 10.1
#: forbids empirical frequency selection, so the estimator may only ever use
#: this list (or a subset of it).
PHYSICAL_PERIODS_D = {
    "lunar_node": 6798.383,      # lunar nodal regression period
    "annual": 365.256363,        # tropical year: Earth's heliocentric motion
    "sidereal_day_nodes": 1.0,   # placeholder slot, never used (see below)
    "anomalistic_month": 29.530589,
    "synodic_month": 29.530589,
    "evection_period": 14.765,   # half synodic: the evection term
    "variation_period": 14.730,
}
#: Periods the estimator must NOT carry: any period equal to or shorter than the
#: orbit period is either invisible in the nodal rate (short-period terms vanish
#: at the ascending node) or already removed by the instrument's own averaging.
DEFAULT_PERIODS = ("lunar_node", "annual", "anomalistic_month", "evection_period")


def instantaneous_node_rate(x: np.ndarray, accel: np.ndarray) -> float:
    """Exact instantaneous dOmega/dt (rad/s) from a state and its acceleration.

    ``x`` is the 6-vector [r, v]; ``accel`` is a(x, t). Returns rad/s.

    The node longitude is ``atan2(-h_x, h_y)``; differentiating the atan2
    identity gives the expression in the module docstring. The denominator
    ``h_x^2 + h_y^2 = |h|^2 cos^2 i`` is exactly why i = 90 deg is
    unmeasurable: at that inclination it is zero. :func:`node_conditioning`
    exposes that number so callers can check it before trusting a result.
    """
    r, v = np.asarray(x[:3], float), np.asarray(x[3:], float)
    h = np.cross(r, v)
    hd = np.cross(r, np.asarray(accel, float))
    denom = h[0] ** 2 + h[1] ** 2
    return (h[1] * (-hd[0]) + h[0] * hd[1]) / denom


def node_conditioning(inc_rad: float) -> float:
    """1/cos^2 i - the conditioning factor of the node-longitude measurement.

    ~54 at SSO, ~1.3e4 at 89.5 deg, undefined at 90 deg. Report it alongside
    every measured rate: a rate obtained at high conditioning is not comparable
    to one obtained at low conditioning without an uncertainty inflation.
    """
    c = math.cos(inc_rad)
    if abs(c) < 1e-12:
        return math.inf
    return 1.0 / c ** 2


def node_rate_to_deg_day(rate_rad_s: np.ndarray | float) -> np.ndarray | float:
    """rad/s -> deg/day."""
    return np.asarray(rate_rad_s) * SEC_PER_DAY * DEG_PER_RAD if np.ndim(rate_rad_s) else \
        float(rate_rad_s) * SEC_PER_DAY * DEG_PER_RAD


def resolve_periods(periods: tuple) -> tuple:
    """Accept either NAMES from :data:`PHYSICAL_PERIODS_D` or numeric days.

    Names are resolved through the fixed physical table, which is what keeps
    the estimator from ever accepting a fitted frequency.
    """
    out = []
    for p in periods:
        if isinstance(p, str):
            if p not in PHYSICAL_PERIODS_D:
                raise KeyError(f"unknown physical period {p!r}; "
                               f"known: {sorted(PHYSICAL_PERIODS_D)}")
            out.append(PHYSICAL_PERIODS_D[p])
        else:
            out.append(float(p))
    return tuple(out)


def design_matrix(t_days: np.ndarray, periods: tuple) -> np.ndarray:
    """[1, t, cos(w_p t), sin(w_p t)] for each fixed physical period."""
    t = np.asarray(t_days, float)
    cols = [np.ones_like(t), t]
    for p in resolve_periods(periods):
        w = 2.0 * math.pi / p
        cols.append(np.cos(w * t))
        cols.append(np.sin(w * t))
    return np.column_stack(cols)


def _vif_linear(t_days: np.ndarray, period_d: float) -> float:
    """Variance inflation of the secular column against one harmonic pair.

    1.0 means the harmonic does not bias the secular slope; large values mean
    the two are not separable on this arc. Pure linear algebra - no physics -
    so it is a property of the arc, not of the force model.
    """
    t = np.asarray(t_days, float)
    w = 2.0 * math.pi / period_d
    A = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    coef, *_ = np.linalg.lstsq(A, t, rcond=None)
    r2 = 1.0 - (t - A @ coef).var() / t.var()
    return 1.0 / (1.0 - r2)


def vif_table(t_days: np.ndarray, periods: tuple) -> dict:
    return {f"{p:g}": _vif_linear(t_days, p) for p in resolve_periods(periods)}


def fit_secular_fixed_frequency(t_days: np.ndarray, rate_deg_day: np.ndarray,
                                periods: tuple = DEFAULT_PERIODS,
                                sample_weights: np.ndarray | None = None) -> dict:
    """Joint fixed-frequency fit of a secular rate plus forced terms.

    Fits ``rate = a0 + a1*t + sum_p [b_p cos(w_p t) + c_p sin(w_p t)]`` with every
    ``w_p`` FIXED at a physical value, and returns the secular slope ``a1`` in
    deg/day.

    Why this is not HARKing: the frequencies are supplied by the caller from
    :data:`PHYSICAL_PERIODS_D`, which is derived from known celestial periods.
    Nothing here searches for a frequency that improves the fit.

    The returned ``slope_uncertainty`` is the classical regression standard
    error. For a deterministic propagator the residuals are model misfit, not
    measurement noise, so it quantifies "how well does this model describe the
    series", NOT "how well is the slope determined". Read ``vif_max``: that is
    the number that says whether the secular column is separable at all.
    """
    t = np.asarray(t_days, float)
    y = np.asarray(rate_deg_day, float)
    A = design_matrix(t, periods)
    if sample_weights is not None:
        w = np.sqrt(np.asarray(sample_weights, float))
        A, y = A * w[:, None], y * w
    # COLUMN EQUILIBRATION -- mandatory, not an optimisation.
    #
    # The raw design matrix is badly scaled: over a multi-year arc the `t`
    # column has norm ~2e4 while every harmonic column has norm ~140, and the
    # secular direction is also nearly parallel to the constant column. Measured
    # cond(A) = 1.76e4 on a 2-cycle arc, and plain `lstsq` applied to a PURELY
    # SECULAR signal returned 1e-24 instead of 1.35e-4 deg/day: the secular
    # coefficient was annihilated outright, not merely imprecise. Solving the
    # equilibrated system and rescaling the secular coefficient back is
    # algebraically identical and is what makes the estimator able to recover a
    # known secular rate at all.
    scale = np.linalg.norm(A, axis=0)
    scale[scale == 0.0] = 1.0
    An = A / scale
    coef_n, *_ = np.linalg.lstsq(An, y, rcond=None)
    coef = coef_n / scale
    resid = y - A @ coef
    n, k = A.shape
    dof = max(n - k, 1)
    sigma2 = float(np.sum(resid ** 2) / dof)
    cov = sigma2 * np.linalg.pinv(An.T @ An) / np.outer(scale, scale)
    slope = float(coef[1])
    vifs = vif_table(t, periods)
    return {
        "slope_deg_day": slope,
        "slope_stderr_deg_day": float(math.sqrt(cov[1, 1])),
        "intercept": float(coef[0]),
        "periods_d": tuple(resolve_periods(periods)),
        "n_samples": int(n),
        "n_params": int(k),
        "rms_resid": float(math.sqrt(np.mean(resid ** 2))),
        "max_abs_harmonic_deg_day": float(np.max(np.abs(coef[2:]))) if k > 2 else 0.0,
        "vif_max": max(vifs.values()),
        "vif_by_period": vifs,
        "cond_raw": float(np.linalg.cond(A)),
        "cond_equilibrated": float(np.linalg.cond(An)),
    }


def trapezoid_cycle_mean(t_days: np.ndarray, rate: np.ndarray) -> float:
    """Time-average of a rate over the sampled window, trapezoid rule.

    Exact for a sum of sinusoids whose period divides the window, which is why
    it is used as an independent cross-check at INTEGER nodal cycles. Note
    ``np.mean`` must not be used here: on ``linspace(0, W, N)`` it divides by N
    while the span is ``W/(N-1)``, an O(1/N) endpoint bias (measured 6.4e-4 ->
    1.6e-5 -> 1.6e-6 for N = 1e3, 4e4, 4e5) that can be mistaken for physics.
    """
    r = np.asarray(rate, float)
    t = np.asarray(t_days, float)
    return float(np.trapezoid(r, t) / (t[-1] - t[0]))


def ols_slope_of_rate(t_days: np.ndarray, rate: np.ndarray) -> float:
    """The estimator the lab used until 2026-10-04, kept as a contrast.

    Included so the results can state quantitatively how wrong it is on a
    near-one-cycle arc, rather than merely asserting that it is.
    """
    t = np.asarray(t_days, float)
    y = np.asarray(rate, float)
    return float(np.polyfit(t, y, 1)[0])