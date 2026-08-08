"""Numerical metrics: error measures and convergence-rate estimation."""

from __future__ import annotations

import numpy as np

__all__ = ["max_abs_error", "l2_norm_error", "relative_l2_error", "convergence_rate"]


def max_abs_error(approx: np.ndarray, exact: np.ndarray) -> float:
    """Maximum absolute error between two arrays of the same shape."""
    return float(np.max(np.abs(np.asarray(approx) - np.asarray(exact))))


def l2_norm_error(approx: np.ndarray, exact: np.ndarray) -> float:
    """L2 norm of the difference between two arrays of the same shape."""
    diff = np.asarray(approx) - np.asarray(exact)
    return float(np.sqrt(np.sum(diff**2)))


def relative_l2_error(approx: np.ndarray, exact: np.ndarray) -> float:
    """Relative L2 error; returns NaN if the reference is exactly zero."""
    exact = np.asarray(exact)
    denom = np.sqrt(np.sum(exact**2))
    if denom == 0:
        return float("nan")
    return float(l2_norm_error(approx, exact) / denom)


def convergence_rate(errors: np.ndarray, stepsizes: np.ndarray) -> np.ndarray:
    """Measured convergence order from error-vs-stepsize data.

    Estimates p in ``error ~ C * h**p`` by log-log slopes between consecutive
    (stepsizes, errors) pairs. Returns an array of length ``len(errors) - 1``.
    """
    errors = np.asarray(errors, dtype=float)
    stepsizes = np.asarray(stepsizes, dtype=float)
    if errors.shape != stepsizes.shape:
        raise ValueError("errors and stepsizes must have the same shape")
    if len(errors) < 2:
        raise ValueError("need at least two data points")
    if np.any(errors <= 0) or np.any(stepsizes <= 0):
        raise ValueError("errors and stepsizes must be strictly positive")
    return np.diff(np.log(errors)) / np.diff(np.log(stepsizes))
