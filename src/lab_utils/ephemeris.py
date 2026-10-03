"""Ephemeris pipeline and third-body force model - canonical Lab source.

Graduated 2026-10-02 from
``research/orbital-mechanics/missions/mission_j2_precession_modulated_lunisolar_term/
mission_experiment.py`` (commit ``3389949``), which held the only copy after the
mission that introduced it completed.

Why this was graduated: the byte-pinned DE441 Sun/Moon loader, its interpolation,
the IAU-1976 precession rotation, and the point-mass third-body acceleration are
the single most load-bearing pieces of physics in the lunisolar chain - and every
mission was reaching into a *completed mission's folder* to import them. That is
exactly the shared-ancestry arrangement that let one DE441 interpolation defect
corrupt four missions at once. The fix for that class of failure is one canonical
home plus equivalence pinning, not four copies.

Equivalence contract: every function here is a **faithful move**, not a rewrite.
``tests/test_ephemeris.py`` pins each one bit-for-bit against the original donor
module. Do not "optimise", re-derive, or tidy these without re-pinning.

Known defect classes preserved here deliberately:

* ``interp_snapshot_buggy`` retains the sign-reversed interpolation weights that
  moved the Mission-2 cross coefficient ``a11`` by ~5x. It is kept so the bug's
  magnitude stays quantifiable and so the regression test can prove the fix holds.
* ``precession_j2000_to_mod`` is the corrected (non-transposed) form. The
  audit-019 transpose bug lived here; ``tests`` assert orthonormality.

Performance note (2026-10-02 profiling, NOT acted on): ``interp_snapshot`` calls
``precession_j2000_to_mod`` once per query, i.e. 8 times per RK4 step in the
``sun_moon_j2`` mode, and ``third_body_accel`` performs several ``np.linalg.norm``
calls on 3-element arrays. Together these are ~80% of the hot-path cost. A 4-6x
speedup is available by precessing the snapshot arrays once at load time and using
scalar math for 3-vectors - both algebraically equivalent. It was deliberately NOT
applied: it would change every committed result number in the repository and the
Lab's compute share of wall-clock is well under 1%. Revisit only when arc lengths
grow (cislunar/NRHO), where fixed-dt becomes the wall.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

from .earth_frames import JD_J2000
from .orbits import MU_EARTH_KM3S2

__all__ = [
    "LUNAR_GM",
    "SOLAR_GM",
    "circular_ic",
    "de441_snapshot_paths",
    "interp_snapshot",
    "interp_snapshot_buggy",
    "load_snapshot",
    "precession_j2000_to_mod",
    "third_body_accel",
]

SOLAR_GM = 132712440018.0
LUNAR_GM = 4902.8001

#: Sub-directory holding the byte-pinned DE441 geocentric Sun/Moon vectors
#: (2026-01-01 -> 2045-01-01, daily, ICRF/TDB), acquired by
#: ``mission_lunisolar_closure``.
DE441_REFERENCE_DIRNAME = Path(
    "research/orbital-mechanics/missions/mission_lunisolar_closure/reference"
)


def de441_snapshot_paths(repo_root: Path) -> tuple[Path, Path]:
    """Return ``(sun_path, moon_path)`` for the byte-pinned DE441 snapshots.

    The authoritative data live *inside the repository* (never on ``R:``), with
    checksums recorded in the sibling ``MANIFEST.json``.
    """
    ref = Path(repo_root) / DE441_REFERENCE_DIRNAME
    return (
        ref / "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt",
        ref / "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt",
    )


def _sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _rot3(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot2(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])


def precession_j2000_to_mod(t_s: float) -> np.ndarray:
    """IAU-1976 precession rotation J2000 -> mean equinox of date.

    Corrected form: ``_rot3`` uses ``[[c, -s], [s, c]]``. The audit-019 defect was a
    TRANSPOSE of this matrix, which left a ~0.66 deg frame mismatch; the fix is
    pinned by the orthonormality tests in ``tests/test_ephemeris.py``.
    """
    T = t_s / (86400.0 * 36525.0)
    sec = math.radians(1.0 / 3600.0)
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * sec
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * sec
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * sec
    return _rot3(-z) @ _rot2(theta) @ _rot3(-zeta)


def load_snapshot(path: Path) -> dict:
    """Parse a Horizons ``$$SOE``/``$$EOE`` vector table into a snapshot dict.

    Returns ``{"t_s", "r_eci_km", "sha256", "n_points"}``. The sha256 travels with
    the data so any result citing the snapshot carries its own provenance.
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    soe = eoe = None
    for i, ln in enumerate(lines):
        if "$$SOE" in ln:
            soe = i + 1
        if "$$EOE" in ln:
            eoe = i
            break
    rows = []
    for ln in lines[soe:eoe]:
        s = ln.strip()
        if not s:
            continue
        p = [x.strip() for x in s.split(",")]
        rows.append((float(p[0]), float(p[2]), float(p[3]), float(p[4])))
    arr = np.array(rows)
    return {"t_s": (arr[:, 0] - JD_J2000) * 86400.0, "r_eci_km": arr[:, 1:4],
            "sha256": _sha256(path), "n_points": len(rows)}


def interp_snapshot(t_q: float, snap: dict, apply_precession: bool = True) -> np.ndarray:
    """Correct linear interpolation (regression-pinned after the 2026-09 fix)."""
    ts = snap["t_s"]
    r = snap["r_eci_km"]
    if t_q <= ts[0]:
        rv = r[0]
    elif t_q >= ts[-1]:
        rv = r[-1]
    else:
        idx = int(np.searchsorted(ts, t_q))
        t_lo, t_hi = ts[idx - 1], ts[idx]
        frac = (t_q - t_lo) / (t_hi - t_lo)
        rv = r[idx - 1] + frac * (r[idx] - r[idx - 1])
    if apply_precession:
        return precession_j2000_to_mod(t_q) @ rv
    return rv


def interp_snapshot_buggy(t_q: float, snap: dict, apply_precession: bool = True) -> np.ndarray:
    """Legacy sign-reversed weights, preserved ONLY for bug-impact quantification.

    ``r[i-1] + frac*(r[i-1] - r[i])`` instead of ``r[i-1] + frac*(r[i] - r[i-1])``.
    This defect moved the Mission-2 cross coefficient from -1.49e-4 (SNR 1.79) to
    -7.85e-4 (SNR 6.89), i.e. it manufactured a 3.9-sigma "detection" out of noise.
    Kept so that magnitude stays reproducible and so the regression test can fail
    if anyone reintroduces the weights.
    """
    ts = snap["t_s"]
    r = snap["r_eci_km"]
    if t_q <= ts[0]:
        rv = r[0]
    elif t_q >= ts[-1]:
        rv = r[-1]
    else:
        idx = int(np.searchsorted(ts, t_q))
        t_lo, t_hi = ts[idx - 1], ts[idx]
        frac = (t_q - t_lo) / (t_hi - t_lo)
        rv = r[idx - 1] + frac * (r[idx - 1] - r[idx])
    if apply_precession:
        return precession_j2000_to_mod(t_q) @ rv
    return rv


def third_body_accel(r_eci: np.ndarray, t_s: float, sun: dict, moon: dict, *,
                     lambda_3body: float = 1.0, include_sun: bool = True,
                     include_moon: bool = True, use_buggy_interp: bool = False) -> np.ndarray:
    """Point-mass differential Sun + Moon acceleration (km/s^2).

    ``lambda_3body`` scales the whole contribution and is the perturbation-scaling
    control used by the Mission-2 / Mission-3 lambda sweeps.
    """
    interp = interp_snapshot_buggy if use_buggy_interp else interp_snapshot
    a = np.zeros(3)
    if include_sun:
        rs = interp(t_s, sun, True)
        d = rs - r_eci
        a += SOLAR_GM * (d / np.linalg.norm(d) ** 3 - rs / np.linalg.norm(rs) ** 3)
    if include_moon:
        rm = interp(t_s, moon, True)
        d = rm - r_eci
        a += LUNAR_GM * (d / np.linalg.norm(d) ** 3 - rm / np.linalg.norm(rm) ** 3)
    return lambda_3body * a


def circular_ic(a_km: float, inc_rad: float, t0_s: float, phase: float = 0.0) -> np.ndarray:
    """Circular two-body initial state at inclination ``inc_rad``, RAAN = 0.

    Every mission in the lunisolar chain starts from this constructor, which is
    why it belongs beside the force model rather than inside any one mission.
    """
    n = math.sqrt(MU_EARTH_KM3S2 / a_km ** 3)
    M = phase
    r_pf = np.array([a_km * math.cos(M), a_km * math.sin(M), 0.0])
    v_pf = np.array([-a_km * n * math.sin(M), a_km * n * math.cos(M), 0.0])
    R = _rot3(0.0) @ np.array([[1, 0, 0], [0, math.cos(inc_rad), -math.sin(inc_rad)],
                               [0, math.sin(inc_rad), math.cos(inc_rad)]])
    # Omega0 = 0
    return np.concatenate([R @ r_pf, R @ v_pf])