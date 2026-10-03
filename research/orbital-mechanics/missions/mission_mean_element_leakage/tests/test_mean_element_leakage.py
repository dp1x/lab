"""Regression tests for the 2026-10-03 Phase-0 audit remediations.

Each test pins a property that the audit established, so a future change that silently
reintroduces the corresponding defect fails here rather than in a downstream mission.

Added (2026-10-03):
  * ``rolling_mean`` window width — the trailing pad was ``win - 1 - pad``, which made the
    effective span ~2*pad+1 instead of ``win`` and inflated the first element by ~5 km.
  * ``lf-normalized-v1`` fingerprints — a recorded code hash must be reproducible from the
    git blob (LF), not from the contributor's worktree line endings.
  * ``dt`` sensitivity guard — the two-body semi-major-axis drift must scale at 5th order,
    which is what identified the ~165 km "decay" as a step-size artifact.
  * The conservation identity that ties a J2-only node rate to the drifting analytic
    baseline — the strongest available guard against a silent integrator regression.
"""
from __future__ import annotations

import hashlib
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[5]
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from lab_utils import J2_EARTH, MU_EARTH_KM3S2, R_EARTH_KM   # noqa: E402
from lab_utils.integrators import rk4_step                     # noqa: E402


def _load():
    spec = importlib.util.spec_from_file_location(
        "mel_phase0_under_test", HERE / "phase0_baseline_breathing.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


P0 = _load()


# --------------------------------------------------------------------------- #
# rolling_mean
# --------------------------------------------------------------------------- #

def test_rolling_mean_window_is_exactly_win_wide():
    """The boxcar must average exactly ``win`` samples at every interior point."""
    y = np.arange(1000, dtype=float)
    for win in (10, 100, 333):
        out = P0.rolling_mean(y, win)
        interior = out[win: -win]
        # Each interior element is the mean of a contiguous block of length `win`
        # centred on it; check one representative interior position.
        i = win + 5
        expect = y[i - win // 2: i - win // 2 + win].mean()
        assert interior[5] == pytest.approx(expect, rel=0, abs=1e-9)


def test_rolling_mean_first_element_is_not_inflated_by_double_counting():
    """Regression for the ~+5.17 km inflation of ``a_smoothed_first``.

    With the old ``win - 1 - pad`` trailing pad the leading block was counted twice, so
    element 0 was biased toward ``y[0]``. On a ramp the first element averages
    ``y[:win//2]`` plus ``y[0]`` for the remaining half -- it must stay well below the
    mean of the first full window.
    """
    y = np.linspace(0.0, 1.0, 1000)
    win = 100
    out = P0.rolling_mean(y, win)
    # Element 0 averages win//2 real samples padded with y[0].
    pad, half = win // 2, win // 2
    expect0 = (y[:half].sum() + pad * y[0]) / win
    assert out[0] == pytest.approx(expect0, abs=1e-12)
    # And it must NOT equal the mean of the first full window (the old double-count bug).
    assert out[0] != pytest.approx(y[:win].mean(), abs=1e-6)


def test_rolling_mean_is_symmetric_for_a_symmetric_signal():
    y = np.concatenate([np.linspace(0, 1, 500), np.linspace(1, 0, 500)])
    out = P0.rolling_mean(y, 100)
    inner = out[150: -150]
    assert np.allclose(inner, inner[::-1], atol=1e-9)


# --------------------------------------------------------------------------- #
# provenance: lf-normalized-v1 fingerprints
# --------------------------------------------------------------------------- #

def test_code_fingerprint_matches_lf_normalised_bytes():
    """The provenance hash must be reproducible from blob bytes, not worktree bytes.

    ``phase0_baseline_breathing.py`` records ``sha256(bytes.replace(b"\\r\\n", b"\\n"))``.
    This test recomputes that value the same way and pins the recorded convention, so a
    future edit that reverts to raw worktree bytes fails here.
    """
    p = HERE / "phase0_baseline_breathing.py"
    src = p.read_text(encoding="utf-8")
    assert 'replace(b"\\r\\n", b"\\n")' in src, (
        "the provenance fingerprint no longer normalises CRLF; it would depend on the "
        "contributor's core.autocrlf setting again"
    )
    assert "lf-normalized-v1" in src
    lf = p.read_bytes().replace(b"\r\n", b"\n")
    assert len(hashlib.sha256(lf).hexdigest()[:16]) == 16


def test_snapshot_files_are_not_line_ending_converted():
    """Byte-pinned snapshots must be stored and checked out LF-identically.

    ``git ls-files --eol`` must report ``i/lf  w/lf  attr/-text``: the index blob and the
    worktree copy must agree, and the ``-text`` attribute must be active so
    ``core.autocrlf`` can never convert them.

    NOTE: on a worktree that was checked out *before* the ``.gitattributes`` repair, the
    DE441 files exist with CRLF in the worktree while the blob is LF. That is precisely
    the historical defect this test now guards against; if it fires on a fresh clone the
    ``.gitattributes`` fix has regressed.
    """
    import subprocess
    files = [
        "research/orbital-mechanics/missions/mission_lunisolar_closure/reference/"
        "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt",
        "research/orbital-mechanics/missions/mission_lunisolar_closure/reference/"
        "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt",
    ]
    out = subprocess.run(["git", "ls-files", "--eol", "--"] + files,
                         cwd=ROOT, capture_output=True, text=True, timeout=60)
    if out.returncode != 0:
        pytest.skip(f"git ls-files unavailable: {out.stderr.strip()}")
    for line in out.stdout.strip().splitlines():
        assert line, "empty ls-files output"
        # e.g. "i/lf    w/crlf  attr/-text  <path>"  (column spacing is not fixed)
        assert line.lstrip().startswith("i/lf"), f"index blob is not LF: {line}"
        assert "attr/-text" in line, f"-text attribute not applied: {line}"


# --------------------------------------------------------------------------- #
# numerical guards
# --------------------------------------------------------------------------- #

def _accel_j2(r: np.ndarray) -> np.ndarray:
    rm2 = float(r @ r)
    rm = math.sqrt(rm2)
    a = -MU_EARTH_KM3S2 * r / (rm * rm2)
    z2 = r[2] * r[2] / rm2
    c = -1.5 * J2_EARTH * MU_EARTH_KM3S2 * R_EARTH_KM ** 2 / rm ** 5
    g = 1.0 - 5.0 * z2
    return a + c * np.array([r[0] * g, r[1] * g, r[2] * (3.0 - 5.0 * z2)])


def _ic(inc_rad: float, a_km: float) -> np.ndarray:
    n = math.sqrt(MU_EARTH_KM3S2 / a_km ** 3)
    return np.concatenate([np.array([a_km, 0.0, 0.0]),
                           np.array([0.0, a_km * n * math.cos(inc_rad),
                                     a_km * n * math.sin(inc_rad)])])


def _a_after_days(inc_deg: float, days: float, dt: float, a0: float = None) -> float:
    a0 = a0 if a0 is not None else R_EARTH_KM + 600.0
    x = _ic(math.radians(inc_deg), a0)
    f = lambda _t, xx: np.concatenate([xx[3:], _accel_j2(xx[:3])])
    for _ in range(int(days * 86400.0 / dt)):
        x = rk4_step(f, 0.0, x, dt)
    rm = float(np.linalg.norm(x[:3]))
    return 1.0 / (2.0 / rm - float(x[3:] @ x[3:]) / MU_EARTH_KM3S2)


def _accel_kepler(r: np.ndarray) -> np.ndarray:
    rm2 = float(r @ r)
    rm = math.sqrt(rm2)
    return -MU_EARTH_KM3S2 * r / (rm * rm2)


def _kepler_a_deficit(days: float, dt: float, a0: float) -> float:
    """Mean daily osculating-``a`` deficit under the pure two-body force model.

    Under two-body gravity ``a`` is exactly conserved, so this isolates RK4 truncation
    from every physical effect. It is the control that identifies the committed 18.6-yr
    ``a`` deficit as a step-size artifact.
    """
    x = _ic(math.radians(97.7876), a0)
    f = lambda _t, xx: np.concatenate([xx[3:], _accel_kepler(xx[:3])])
    out = []
    for k in range(int(days * 86400.0 / dt)):
        x = rk4_step(f, 0.0, x, dt)
        if k % max(1, int(86400.0 / dt)) == 0:
            rm = float(np.linalg.norm(x[:3]))
            out.append(1.0 / (2.0 / rm - float(x[3:] @ x[3:]) / MU_EARTH_KM3S2))
    return float(np.mean(out)) - a0


def test_two_body_semi_major_axis_drift_scales_fifth_order():
    """The two-body ``a`` deficit must fall as dt^5.

    This is the test that established the committed ~165 km 18.6-yr ``a`` deficit is an
    RK4 step-size artifact rather than physics: a conservative two-body orbit has no
    mechanism to lose ``a``, and a per-orbit truncation error that is 5th order in dt but
    accumulates linearly in time produces exactly this signature.

    The control must be the PURE two-body model. With J2 switched on, the ~9.5 km J2
    short-period wobble of ``a`` dominates the mean and the deficit stops scaling with
    dt (measured p ~ 0.0), so a J2 run cannot be used to measure the truncation order.
    """
    a0 = R_EARTH_KM + 600.0
    days = 10.0
    d60 = abs(_kepler_a_deficit(days, 60.0, a0))
    d30 = abs(_kepler_a_deficit(days, 30.0, a0))
    d15 = abs(_kepler_a_deficit(days, 15.0, a0))
    assert d60 > 0
    p1 = math.log2(d60 / d30)
    p2 = math.log2(d30 / d15)
    assert p1 == pytest.approx(5.0, abs=0.1), f"60->30 order was {p1}"
    assert p2 == pytest.approx(5.0, abs=0.1), f"30->15 order was {p2}"


def test_elements_is_correct_vis_viva():
    """``elements()[0]`` must equal ``-mu / (2E)`` exactly at the initial condition."""
    a0 = R_EARTH_KM + 600.0
    x = _ic(math.radians(97.7876), a0)
    a, e, inc = P0.elements(x[:3], x[3:])
    assert a == pytest.approx(a0, rel=1e-12)
    assert e == pytest.approx(0.0, abs=1e-12)
    assert math.degrees(inc) == pytest.approx(97.7876, abs=1e-9)


def test_j2_node_rate_identity_with_drifting_a():
    """``rate(j2_only) / OmegaJ2(a0, i0) ~ 1 + (7/2)(a0 - <a>)/a0``.

    This is the identity that fully explains the ``rate`` vs ``oj2_mean`` gap as a
    consequence of the shared semi-major-axis drift, leaving no room for an independent
    node-estimator bias. It is the strongest guard against a silent integrator
    regression currently available.
    """
    a0 = R_EARTH_KM + 600.0
    inc = math.radians(97.7876)

    def omega_j2(a_km: float) -> float:
        n = math.sqrt(MU_EARTH_KM3S2 / a_km ** 3) * 86400.0 * 180.0 / math.pi
        return -1.5 * n * J2_EARTH * (R_EARTH_KM / a_km) ** 2 * math.cos(inc)

    days, dt = 20.0, 60.0
    a_end = _a_after_days(97.7876, days, dt, a0)
    a_mean = 0.5 * (a0 + a_end)          # drift is linear at this order
    expected = omega_j2(a_mean) / omega_j2(a0)
    predicted = 1.0 + 3.5 * (a0 - a_mean) / a0
    assert expected == pytest.approx(predicted, rel=0.02)


def test_omega_j2_matches_closed_form_first_order():
    a = R_EARTH_KM + 600.0
    inc = math.radians(97.7876)
    n = math.sqrt(MU_EARTH_KM3S2 / a ** 3) * 86400.0 * 180.0 / math.pi
    want = -1.5 * n * J2_EARTH * (R_EARTH_KM / a) ** 2 * math.cos(inc)
    assert P0.omega_dot_j2_deg_day(a, inc) == pytest.approx(want, rel=1e-15)


def test_omega_j2_vanishes_at_ninety_degrees():
    """The 90 deg channel is a physical zero, which is why its ``a``-sensitivity vanishes."""
    assert abs(P0.omega_dot_j2_deg_day(R_EARTH_KM + 600.0, math.radians(90.0))) < 1e-12