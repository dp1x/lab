"""Estimation doctrine - canonical Lab source for rate estimation and uncertainty.

Graduated 2026-10-02. This module exists because a 2026-10-02 repository audit
found **six independent ``ols_slope`` implementations across four files**, with
**inconsistent return types** (three returning ``(slope, sigma)`` tuples, three
returning bare floats), and **thirteen** separate campaign runners across four
missions. Every mission had re-derived its own covariance, which is why the
Mission-2 cross coefficient could be reported as "SNR 6.89" without any of the
significance machinery being auditable.

It is also constitution ``LAB_CONSTITUTION.md`` section 12.3 backlog candidate #1,
"Estimation Doctrine Graduation".

The caveat that matters
----------------------
``ols_slope_ci`` uses the classical i.i.d.-residual covariance. **For a
deterministic propagator there is no noise**: the residuals are *model misfit*,
not measurement error. The returned standard error therefore quantifies "how well
does a straight line fit this series", NOT "how well is the physical slope
determined". Mission 2 reported ``|a11| / sigma`` as an SNR and it was read as a
statistical detection; it was not one. Adopting scipy does **not** fix this -
``scipy.stats.linregress`` makes the identical assumption - and the
``localdocs/reports/scientific-stack-audit-2026-10-02.md`` audit records that no
library upgrade would have prevented any error the Lab has actually made.

What this module does fix is that the assumption is now stated in **one** place,
in **one** function, instead of being re-derived per mission.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import signal, stats

__all__ = [
    "bootstrap_ci",
    "cadence_uniformity",
    "dominant_period_days",
    "ols_slope",
    "ols_slope_ci",
    "rate_deg_day",
]

SEC_PER_DAY = 86400.0
DEG_PER_RAD = 180.0 / math.pi

#: Peak element count for one lombscargle block (~32 MB at float64). The host has
#: no pagefile, so a full (n_samples x n_freqs) broadcast would raise MemoryError
#: on any real node series.
FREQ_CHUNK_ELEMENTS = 4_000_000


def ols_slope(t_s: np.ndarray, y_rad: np.ndarray) -> float:
    """Ordinary-least-squares slope of ``y`` against ``t``.

    Faithful re-implementation of the six historical copies: the design matrix and
    ``lstsq`` call are identical, so results are bit-identical to the donors that
    Mission 2, Mission 3 and ``mission_mean_element_leakage`` currently import.
    ``tests/test_estimation.py`` pins that equivalence against the real 1-yr
    ascending-node series.

    Returns a **point estimate only**. Use :func:`ols_slope_ci` if an uncertainty
    is needed, and read that function's caveat first.
    """
    A = np.column_stack([np.ones_like(t_s), t_s])
    coef, *_ = np.linalg.lstsq(A, y_rad, rcond=None)
    return float(coef[1])


def ols_slope_ci(t_s: np.ndarray, y_rad: np.ndarray) -> tuple[float, float]:
    """Return ``(slope, standard_error)`` under the classical i.i.d. assumption.

    WARNING
    -------
    The standard error is computed from fit residuals, which is only a valid
    uncertainty when the residual scatter reflects random measurement noise.
    For deterministic RK4 output the residuals are *model misfit*. Report such a
    number as a goodness-of-fit diagnostic, never as a physical detection.

    Cross-checked against ``scipy.stats.linregress.stderr`` in the tests.
    """
    n = len(t_s)
    if n < 3:
        raise ValueError("ols_slope_ci needs at least 3 samples")
    slope = ols_slope(t_s, y_rad)
    A = np.column_stack([np.ones_like(t_s), t_s])
    resid = y_rad - A @ np.array([np.mean(y_rad) - slope * np.mean(t_s), slope])
    sigma2 = float(np.sum(resid ** 2) / (n - 2))
    cov = sigma2 * np.linalg.inv(A.T @ A)
    return slope, float(math.sqrt(cov[1, 1]))


def rate_deg_day(t_s: np.ndarray, y_rad: np.ndarray) -> float:
    """Angular rate in deg/day from epoch-seconds and radians.

    The canonical converter used by every lunisolar campaign. ``t_s`` is seconds
    since J2000 TDB; ``y_rad`` is an unwrapped angle in radians.
    """
    return ols_slope(t_s, y_rad) * SEC_PER_DAY * DEG_PER_RAD


def bootstrap_ci(data: np.ndarray, *, statistic=np.mean, n_resamples: int = 10000,
                 alpha: float = 0.05, seed: int = 0) -> tuple[float, float]:
    """Percentile bootstrap CI for ``statistic`` over ``data``.

    Wraps :func:`scipy.stats.bootstrap` with a **fixed seed** so results are
    reproducible, per AGENTS.md rule 1 ("deterministic only"). Pass ``statistic``
    as a 1-D-argument callable.
    """
    res = stats.bootstrap(
        (np.asarray(data, dtype=float),), statistic,
        n_resamples=n_resamples, confidence_level=1.0 - alpha,
        method="percentile", random_state=np.random.default_rng(seed),
    )
    lo, hi = res.confidence_interval
    return float(lo), float(hi)


def cadence_uniformity(t_s: np.ndarray) -> dict:
    """Describe the spacing irregularity of a sample series.

    Added because the Lab applies uniform-grid FFT machinery to ascending-node
    crossings, which are **not** uniformly spaced. Measured across the 1-yr Phase-0
    run the node count varied 5445 / 5449 / 5462 between force modes, i.e. ~0.3%
    spacing variation - small enough that the uniform-grid assumption in Exp 019's
    FFT arm is approximately fine, and *not* a bug. This function makes the
    assumption measurable rather than implicit.

    Returns relative standard deviation, coefficient of variation, and the peak
    fraction of samples whose spacing deviates from the mean by more than 1%.
    """
    t = np.asarray(t_s, dtype=float)
    if t.size < 3:
        raise ValueError("cadence_uniformity needs at least 3 samples")
    d = np.diff(t)
    mean_d = float(np.mean(d))
    return {
        "n": int(t.size),
        "mean_dt_s": mean_d,
        "std_dt_s": float(np.std(d)),
        "rel_std": float(np.std(d) / mean_d),
        "cv": float(np.std(d) / mean_d),
        "max_abs_frac_dev": float(np.max(np.abs(d - mean_d)) / mean_d),
        "frac_beyond_1pct": float(np.mean(np.abs(d - mean_d) / mean_d > 0.01)),
    }


def dominant_period_days(t_s: np.ndarray, y_rad: np.ndarray,
                         n_periods: int = 5,
                         n_grid: int = 20000,
                         min_separation: float = 1.10) -> list[tuple[float, float]]:
    """Lomb-Scargle periodogram of an unevenly sampled series.

    Returns up to ``n_periods`` ``(period_days, power)`` pairs, strongest first,
    de-duplicated so that sidelobes of one tone are not reported as separate
    periods. ``min_separation`` is the minimum ratio between two reported
    periods (1.10 = 10% apart), which suppresses sidelobes while still resolving
    physically distinct harmonics.

    ADDED, NOT YET APPLIED. Exp 019's periodicity arm used ``np.fft.rfft`` with a
    mean-spacing grid assumption. For unevenly sampled data the correct tool is
    Lomb-Scargle, which is what this function wraps. Migrating Exp 019 is a
    *scientific* change that would alter a committed result, so it is deliberately
    left for a mission with a pre-registered rule covering it.
    """
    t = np.asarray(t_s, dtype=float)
    y = np.asarray(y_rad, dtype=float)
    y = y - np.mean(y)          # remove DC (scipy >=1.17 deprecates precenter=)
    t_days = (t - t[0]) / SEC_PER_DAY
    span = max(float(t_days[-1] - t_days[0]), 1.0)
    # Log-spaced grid: a linear frequency grid cannot resolve long periods
    # (a uniform grid spanning 1/span..5 cyc/day resolves a 180 d period only to
    # ~4.5%). Peaks are then refined by a parabolic fit in log-period space.
    periods = np.geomspace(0.2, span, n_grid)
    freqs = 1.0 / periods
    # Chunk the frequency sweep. scipy's lombscargle broadcasts to (n_samples,
    # n_freqs), so a full 100k-node x 20k-grid sweep would want ~16 GB. The Lab
    # host has NO pagefile (see plan.md), so peak allocation is capped explicitly
    # rather than left to the OS.
    power = np.empty_like(freqs)
    chunk = max(1, FREQ_CHUNK_ELEMENTS // max(t_days.size, 1))
    for start in range(0, freqs.size, chunk):
        sl = slice(start, start + chunk)
        power[sl] = signal.lombscargle(t_days, y, 2.0 * np.pi * freqs[sl],
                                       normalize=True)

    def _refine(i: int) -> float:
        """Parabolic interpolation in log-power, returned in period space."""
        p = float(periods[i])
        if 0 < i < len(periods) - 1:
            y0, y1, y2 = np.log(power[i - 1]), np.log(power[i]), np.log(power[i + 1])
            denom = y0 - 2.0 * y1 + y2
            if denom > 0:                       # concave in log-power -> real peak
                delta = 0.5 * (y0 - y2) / denom
                if abs(delta) <= 1.0:
                    step = np.log(periods[i + 1]) - np.log(periods[i])
                    p = float(np.exp(np.log(p) + delta * step))
        return p

    # Select DISTINCT local maxima. Taking the raw top-N by power returns
    # sidelobes of a single strong tone (e.g. 40.00 / 39.98 / 40.01 d) and is
    # useless for the job this exists for: identifying separate harmonics.
    interior = np.arange(1, len(power) - 1)
    is_peak = (power[interior] > power[interior - 1]) & (power[interior] > power[interior + 1])
    peaks = interior[is_peak]
    if peaks.size == 0:
        return [(_refine(int(np.argmax(power))), float(np.max(power)))]
    peaks = peaks[np.argsort(power[peaks])[::-1]]

    chosen: list[tuple[float, float]] = []
    for i in peaks:
        p = _refine(int(i))
        if all(abs(np.log(p) - np.log(q)) > np.log(min_separation)
               for q, _ in chosen):
            chosen.append((p, float(power[i])))
            if len(chosen) >= n_periods:
                break
    return chosen