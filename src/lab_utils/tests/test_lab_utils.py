"""Unit tests for lab_utils shared utilities."""

import os
import tempfile

import numpy as np
import pytest

from lab_utils.metrics import (
    convergence_rate,
    l2_norm_error,
    max_abs_error,
    relative_l2_error,
)
from lab_utils.results import save_json_result


def test_max_abs_error():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([1.1, 2.0, 2.7])
    assert max_abs_error(a, b) == pytest.approx(0.3)


def test_l2_norm_error():
    a = np.zeros(3)
    b = np.array([1.0, 2.0, 2.0])
    assert l2_norm_error(a, b) == pytest.approx(3.0)


def test_relative_l2_error():
    a = np.array([1.0, 2.0])
    b = np.array([2.0, 4.0])
    assert relative_l2_error(a, b) == pytest.approx(0.5)


def test_relative_l2_error_zero_reference():
    assert np.isnan(relative_l2_error(np.zeros(2), np.zeros(2)))


def test_convergence_rate_second_order():
    h = np.array([0.1, 0.05, 0.025])
    e = np.array([0.01, 0.0025, 0.000625])
    rates = convergence_rate(e, h)
    assert np.allclose(rates, 2.0, atol=1e-6)


def test_convergence_rate_validation():
    with pytest.raises(ValueError):
        convergence_rate(np.array([1.0]), np.array([0.1]))
    with pytest.raises(ValueError):
        convergence_rate(np.array([0.0, 1.0]), np.array([0.1, 0.05]))


def test_save_json_result_roundtrip():
    with tempfile.TemporaryDirectory() as tmp:
        path = save_json_result(
            os.path.join(tmp, "out.json"),
            {"metrics": {"error": 0.123456789012345, "arr": np.array([1, 2])}},
            name="test",
        )
        assert os.path.exists(path)
        import json

        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        assert data["meta"]["name"] == "test"
        assert data["results"]["metrics"]["arr"] == [1, 2]
        # Float rounded to 12 significant decimals to avoid repr noise.
        assert abs(data["results"]["metrics"]["error"] - 0.123456789012) < 1e-9
        assert "git_commit" in data["meta"]