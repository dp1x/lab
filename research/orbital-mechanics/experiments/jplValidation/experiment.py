"""Experiment 013: JPL Horizons / real-ephemeris validation (ISS primary target).

Research question: how well do progressively enriched deterministic laboratory
force models (M1 two-body, M2 two-body+J2, M3 two-body+J2+drag) reproduce an
authoritative external ephemeris -- NASA JPL Horizons geometric ICRF vector
states of the International Space Station (-125544), TLE/SGP4 provenance --
and how can the discrepancy be decomposed into integration error, reference
sampling/interpolation effects, time-system bookkeeping, initialization,
constants/conventions, and a jointly-attributed remainder?

Deterministic-offline doctrine: all numerics consume ONLY the checksum-pinned
snapshot under ``reference/`` (verified against ``MANIFEST.json`` before any
parsing; hard fail on mismatch). No network access exists anywhere in this
file; online acquisition lives in ``fetch_horizons_snapshot.py`` (one-time).

Pre-registered analysis commitments (fixed before any residual was computed):
  * window/cadence exactly as acquired (2026-Aug-24 .. -27 TDB, 120 s);
  * residuals evaluated ONLY at snapshot epochs (exact-grid alignment; no
    interpolation feeds any headline number);
  * primary metric: RIC along-track RMS (frame built FROM THE REFERENCE state);
  * integration grid: internal substeps landing exactly on snapshot epochs;
    headline dt = 120 s / 8; integration bound = NSUB 16-vs-8 difference;
  * constants frozen to the lab canon (GM 398600.4418, J2 1.082629821e-3,
    Vallado atmosphere verbatim via the Exp 010 donor);
  * ballistic coefficient band {50, 100, 200, 400} kg/m^2 (Exp 010 sweep),
    primary beta = 100 kg/m^2 (implies A_eff = m/(Cd*beta) ~= 1906 m^2 with
    public ISS mass 419,725 kg, Cd = 2.2 -- plausible mid-range projected
    area); NO refitting of any atmospheric or beta parameter permitted;
  * reference acknowledged as SGP4/TLE-provenance whose own uncertainty grows
    ~1-3 km/day from the trajectory revision date (JPL documentation); it is
    an authoritative external reference, NOT metaphysical ground truth;
  * decision rule: an improvement claim requires exceeding BOTH the seeded
    block-bootstrap 95% CI AND the declared reference envelope; otherwise the
    automatic label is "indistinguishable given reference uncertainty";
  * pre-registered rejection: a second-difference jump > 100 m between
    consecutive along-track residuals flags maneuver/TLE-handover
    contamination and rejects the window before any model comparison.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from pathlib import Path

import numpy as np

from lab_utils.integrators import rk4_propagate
from lab_utils.orbits import (
    J2_EARTH,
    MU_EARTH_KM3S2,
    OMEGA_EARTH_RAD_S,
    R_EARTH_KM,
    coe_to_rv_eci,
    j2_rhs,
    mean_motion,
    rv_to_coe_eci,
)

EXP_DIR = Path(__file__).resolve().parent
REFERENCE_DIR = EXP_DIR / "reference"
MANIFEST_PATH = REFERENCE_DIR / "MANIFEST.json"
VECTORS_NAME = "horizons_-125544_iss_vectors_2026-08-24_to_2026-08-27_tdb_2min.txt"
OBJDATA_NAME = "horizons_-125544_iss_objdata.txt"

CADENCE_S = 120.0
EXPECTED_ROWS = 2161
HEADLINE_NSUB = 8
LADDER_NSUB = (1, 2, 4, 8, 16)
BETA_BAND_KG_M2 = (50.0, 100.0, 200.0, 400.0)
BETA_PRIMARY = 100.0
MU_SGP4_WGS72 = 398600.8  # SGP4 heritage GM (Vallado); sensitivity variant only
MU_DE440 = 398600.435507  # JPL DE440 planet-only GM; sensitivity variant only
JUMP_THRESH_M = 100.0  # pre-registered maneuver/TLE-handover flag
ENVELOPE_KM_PER_DAY = 3.0  # upper documented TLE-predict degradation rate
BOOTSTRAP_BLOCKS_DAY = 720  # 1-day blocks (samples)
BOOTSTRAP_B = 200
BOOTSTRAP_SEED = 137
REVISION_OFFSET_DAYS = 1.0  # window starts one day after revision date

# --------------------------------------------------------------------------- #
# Donor machinery (Exp 010 orbitDecay): declared-atmosphere drag propagator.
# --------------------------------------------------------------------------- #
_OD_PATH = EXP_DIR.parent / "orbitDecay" / "experiment.py"
_spec = importlib.util.spec_from_file_location("od010_for_jpl", _OD_PATH)
assert _spec is not None and _spec.loader is not None
_od = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_od)
propagate_3d_rk4_drag = _od.propagate_3d_rk4_drag


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def code_hashes() -> dict:
    """SHA-256 of every source file this result depends on (stale-run guard)."""
    root = Path(__file__).resolve().parents[4]  # repo root
    files = {
        "experiment.py": Path(__file__).resolve(),
        "lab_utils/orbits.py": root / "src" / "lab_utils" / "orbits.py",
        "lab_utils/integrators.py": root / "src" / "lab_utils" / "integrators.py",
        "fetch_horizons_snapshot.py": EXP_DIR / "fetch_horizons_snapshot.py",
    }
    return {name: sha256_bytes(p.read_bytes()) for name, p in files.items()}


# --------------------------------------------------------------------------- #
# Calendar handling: JDTDB column is the sole clock. The calendar string is
# validated against it (column-alignment check) and then discarded. Two
# INDEPENDENT day-ordinal implementations must agree exactly.
# --------------------------------------------------------------------------- #
def _ordinal_formula(y: int, m: int, d: int) -> int:
    """Proleptic Gregorian day ordinal (0001-01-01 = 1), pure arithmetic
    (Hinnant days_from_civil rebased: unix_days - 719468 + 719163)."""
    yy = y - (1 if m <= 2 else 0)
    era = (yy if yy >= 0 else yy - 399) // 400
    yoe = yy - era * 400
    doy = (153 * (m + (-3 if m > 2 else 9)) + 2) // 5 + d - 1
    doe = yoe * 365 + yoe // 4 - yoe // 100 + doy
    return era * 146097 + doe - 305


_MONTH_LEN = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


def _ordinal_table(y: int, m: int, d: int) -> int:
    """Independent ordinal pipeline: cumulative month lengths + leap rules."""
    def leap(year: int) -> bool:
        return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

    ord_ = d
    for mm in range(1, m):
        ord_ += _MONTH_LEN[mm - 1]
        if mm == 2 and leap(y):
            ord_ += 1
    y0 = y - 1
    ord_ += y0 * 365 + y0 // 4 - y0 // 100 + y0 // 400
    return ord_


_CAL_RE = re.compile(
    r"A\.D\.\s+(\d{4})-([A-Za-z]{3})-(\d{2})\s+(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)"
)
_MONTH_NUM = {
    "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
    "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12,
}


def calendar_epoch_s(cal: str, jd_tag: float, jd0_tag: float) -> float:
    """Validate the calendar string against its own row's Julian Date via the
    two independent day-ordinal pipelines (day-level + time-of-day level),
    then return seconds-since-t0 on the SAME JD base as the primary clock.
    Raises on any calendar/tag corruption."""
    m = _CAL_RE.match(cal.strip())
    if m is None:
        raise ValueError(f"unparsable calendar field: {cal!r}")
    y, mon, d, hh, mm, ss = (
        int(m.group(1)), _MONTH_NUM[m.group(2)], int(m.group(3)),
        int(m.group(4)), int(m.group(5)), float(m.group(6)),
    )
    o_f = _ordinal_formula(y, mon, d)
    o_t = _ordinal_table(y, mon, d)
    if o_f != o_t:
        raise RuntimeError(f"independent ordinal pipelines disagree: {o_f} != {o_t}")
    day_noonbased = int(np.floor(jd_tag + 0.5))
    if day_noonbased - 1721425 != o_f:
        raise ValueError(
            f"calendar date {y:04d}-{mon:02d}-{d:02d} inconsistent with JD tag "
            f"{jd_tag!r} (JDN->ordinal {day_noonbased - 1721425})"
        )
    sod_cal = hh * 3600.0 + mm * 60.0 + ss
    frac = jd_tag + 0.5 - np.floor(jd_tag + 0.5)
    if abs(sod_cal - frac * 86400.0) > 1e-3:
        raise ValueError(
            f"time-of-day {sod_cal:.3f} s disagrees with JD fraction {frac * 86400.0:.3f} s"
        )
    return (jd_tag - jd0_tag) * 86400.0


# --------------------------------------------------------------------------- #
# Snapshot loading: hash enforcement BEFORE parsing; structural + identity +
# plausibility gates; dual independent parse pipelines compared cell-by-cell.
# --------------------------------------------------------------------------- #
def load_manifest() -> dict:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if manifest.get("manifest_version") != 1:
        raise RuntimeError("unsupported MANIFEST version")
    return manifest


def verify_snapshot_bytes(manifest: dict) -> tuple[bytes, bytes]:
    """Recompute SHA-256 of both raw artifacts; hard fail on ANY mismatch."""
    paths = {"vectors": REFERENCE_DIR / VECTORS_NAME, "objdata": REFERENCE_DIR / OBJDATA_NAME}
    raw = {}
    for key, p in paths.items():
        data = p.read_bytes()
        rec = manifest["snapshot"]["files"][p.name]
        if len(data) != rec["bytes"]:
            raise RuntimeError(f"{p.name}: byte count {len(data)} != pinned {rec['bytes']}")
        digest = sha256_bytes(data)
        if digest != rec["sha256"] or digest != manifest["acquisition"]["response_sha256"][key]:
            raise RuntimeError(f"{p.name}: SHA-256 mismatch ({digest} != {rec['sha256']})")
        raw[key] = data
    return raw["vectors"], raw["objdata"]


_HEADER_ASSERTIONS = (
    ("target identity", lambda hdr: "International Space Station" in hdr and "-125544" in hdr),
    ("center identity", lambda hdr: "Center body name: Earth (399)" in hdr),
    ("geocentric site", lambda hdr: "BODY CENTER" in hdr or "GEOCENTRIC" in hdr),
    ("ICRF frame", lambda hdr: re.search(r"Reference frame\s*:\s*ICRF", hdr) is not None),
    ("KM-S units", lambda hdr: re.search(r"Output units\s*:\s*KM-S", hdr) is not None),
    ("geometric states", lambda hdr: "GEOMETRIC cartesian states" in hdr),
    ("TDB epochs", lambda hdr: "TDB" in hdr),
)


def parse_header(text: str) -> dict:
    for name, ok in _HEADER_ASSERTIONS:
        if not ok(text):
            raise RuntimeError(f"header assertion failed: {name}")
    meta = {
        "target_body_name": _grab(text, r"Target body name:\s*(.+)"),
        "center_body_name": _grab(text, r"Center body name:\s*(.+)"),
        "reference_frame": _grab(text, r"Reference frame\s*:\s*(\S+)"),
        "output_units": _grab(text, r"Output units\s*:\s*(\S+)"),
        "output_type": _grab(text, r"Output type\s*:\s*(.+)"),
        "eop_file": _grab(text, r"EOP file\s*:\s*(\S+)"),
    }
    return meta


def _grab(text: str, pattern: str) -> str:
    m = re.search(pattern, text)
    return m.group(1).strip() if m else "<missing>"


_FLOAT_RE = re.compile(r"^[+-]?\d+(?:\.\d+)?(?:[Ee][+-]\d\d)?$")


def _parse_rows_split(lines: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Pipeline A: strict split-based CSV parsing (trailing comma fingerprint)."""
    jds, states, cals = [], [], []
    for k, ln in enumerate(lines):
        parts = ln.split(",")
        if len(parts) != 9 or parts[-1] != "":
            raise RuntimeError(f"row {k}: expected 8 fields + trailing-comma artifact, got {len(parts)}")
        jd_str, cal = parts[0].strip(), parts[1].strip()
        if not _FLOAT_RE.match(jd_str):
            raise RuntimeError(f"row {k}: malformed JD tag {jd_str!r}")
        vals = [float(x) for x in parts[2:8]]
        for x in parts[2:8]:
            if not _FLOAT_RE.match(x.strip()):
                raise RuntimeError(f"row {k}: malformed numeric field {x!r}")
        jds.append(float(jd_str))
        cals.append(cal)
        states.append(vals)
    return np.array(jds), np.array(states), cals


def _parse_rows_regex(lines: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Pipeline B: independent whole-row regex capture (no split logic).

    Note: Horizons prints positive mantissas WITHOUT an explicit '+' sign;
    only negatives carry '-'. The pattern accepts either. Returns (jd, states)
    mirroring pipeline A's contract so results can be compared directly."""
    num = r"[+-]?\d+(?:\.\d+)?[Ee][+-]\d\d"
    row_re = re.compile(
        r"\s*(\d+\.\d+),\s*A\.D\.[^,]+," + (r"\s*(%s)," % num) * 6
    )
    jds, out = [], []
    for k, ln in enumerate(lines):
        m = row_re.match(ln)
        if m is None:
            raise RuntimeError(f"row {k}: regex pipeline failed on {ln[:60]!r}")
        jds.append(float(m.group(1)))
        out.append([float(m.group(i)) for i in range(2, 8)])
    return np.array(jds), np.array(out)


def _plausibility_gates(r: np.ndarray, v: np.ndarray, cadence_s: float = CADENCE_S) -> None:
    """Physical self-consistency gates on parsed states (loader + mutants)."""
    # Per-component magnitude gates (catch axis swaps norm gates cannot).
    if not (np.all(np.abs(r) < 8000.0) and np.all(np.abs(v) < 9.0)):
        raise RuntimeError("component-level magnitude gate violated")
    rn = np.linalg.norm(r, axis=1)
    vn = np.linalg.norm(v, axis=1)
    if not (np.all((6500.0 < rn) & (rn < 7500.0)) and np.all((7.0 < vn) & (vn < 8.0))):
        raise RuntimeError(f"norm plausibility gate violated: |r| [{rn.min()}, {rn.max()}], |v| [{vn.min()}, {vn.max()}]")

    # Dynamics continuity: consecutive-position chords must match the
    # published speeds (catches velocity sign/column corruption).
    chord = np.linalg.norm(np.diff(r, axis=0), axis=1)
    speed_mean = 0.5 * (vn[:-1] + vn[1:])
    ratio = chord / (speed_mean * cadence_s)
    if np.any(ratio < 0.95) or np.any(ratio > 1.01):
        raise RuntimeError(f"position-chord/speed continuity violated: [{ratio.min():.4f}, {ratio.max():.4f}]")
    # Angular-momentum direction from consecutive positions must agree with
    # h = r x v (catches velocity sign flips that preserve norms).
    r_mid = 0.5 * (r[:-1] + r[1:])
    h_pub = np.cross(r_mid, 0.5 * (v[:-1] + v[1:]))
    h_imp = np.cross(r[:-1], r[1:]) / cadence_s
    align = np.einsum("ij,ij->i", h_pub, h_imp) / (
        np.linalg.norm(h_pub, axis=1) * np.linalg.norm(h_imp, axis=1)
    )
    if np.min(align) < 0.999:
        raise RuntimeError(f"published h inconsistent with motion-implied h: min cos {np.min(align):.6f}")

    # Frame sanity: osculating inclination must be ISS-like (~51.6 deg); an
    # ecliptic-frame delivery or pole convention error moves this grossly.
    el = rv_to_coe_eci(r, v)
    inc_deg = np.degrees(el["inc"])
    if abs(np.mean(inc_deg) - 51.63) > 0.75:
        raise RuntimeError(f"inclination gate: mean {np.mean(inc_deg):.3f} deg not ISS-like")


def load_reference() -> dict:
    """Verify integrity, parse, and gate the pinned snapshot. Offline-only."""
    manifest = load_manifest()
    vectors_raw, objdata_raw = verify_snapshot_bytes(manifest)
    text = vectors_raw.decode("utf-8", errors="strict")

    header_meta = parse_header(text)
    lines = text.splitlines()
    soe = [i for i, ln in enumerate(lines) if ln.strip() == "$$SOE"]
    eoe = [i for i, ln in enumerate(lines) if ln.strip() == "$$EOE"]
    if len(soe) != 1 or len(eoe) != 1 or eoe[0] - soe[0] - 1 != EXPECTED_ROWS:
        raise RuntimeError(
            f"table structure broken: SOE={soe} EOE={eoe} "
            f"(expected {EXPECTED_ROWS} rows)"
        )
    rows = lines[soe[0] + 1 : eoe[0]]

    jd, states_a, cal = _parse_rows_split(rows)
    jd_b, states_b = _parse_rows_regex(rows)
    if not (np.array_equal(jd, jd_b) and np.array_equal(states_a, states_b)):
        raise RuntimeError("independent parse pipelines disagree")

    # Column layout: JDTDB, Calendar(TDB), X, Y, Z, VX, VY, VZ  (km, km/s).
    r = states_a[:, 0:3]
    v = states_a[:, 3:6]

    # Time: JDTDB is the sole clock; spacing uniform within print quantization.
    t_s = (jd - jd[0]) * 86400.0
    dts = np.diff(t_s)
    if np.any(dts <= 0):
        raise RuntimeError("epoch tags not strictly monotonic")
    if np.max(np.abs(dts - CADENCE_S)) > 2e-4:
        raise RuntimeError(f"cadence deviation beyond print quantization: max {np.max(np.abs(dts - CADENCE_S)):.2e} s")
    for k in (0, EXPECTED_ROWS // 2, EXPECTED_ROWS - 1):
        cal_s = calendar_epoch_s(cal[k], jd[k], jd[0])
        if abs(cal_s - t_s[k]) > 61.0:
            raise RuntimeError(f"calendar/JDTDB disagreement at row {k}: {cal_s - t_s[k]:.3f} s")

    _plausibility_gates(r, v)

    obj_text = objdata_raw.decode("utf-8", errors="strict")
    if "Trajectory is TLE-based" not in obj_text:
        raise RuntimeError("object sheet lacks TLE-provenance disclosure")
    rev = re.search(r"Revised:\s*([A-Za-z]{3}\s+\d{1,2},\s+\d{4})", obj_text)

    return {
        "manifest": manifest,
        "header": header_meta,
        "jd": jd,
        "t_s": t_s,
        "r": r,
        "v": v,
        "cal_first": cal[0],
        "rows": len(rows),
        "revision_date": rev.group(1) if rev else "<missing>",
        "elements_t0": {k: (float(x) if np.ndim(x) == 0 else x) for k, x in rv_to_coe_eci(r[0], v[0]).items()},
    }


# --------------------------------------------------------------------------- #
# Local orbital frame (RTN/RSW), built FROM THE REFERENCE state.
# --------------------------------------------------------------------------- #
def ric_frames(r_ref: np.ndarray, v_ref: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    h = np.cross(r_ref, v_ref)
    hn = np.linalg.norm(h, axis=-1, keepdims=True)
    if np.any(hn < 1e-12):
        raise RuntimeError("rectilinear-degenerate state in RIC construction")
    c_hat = h / hn
    r_hat = r_ref / np.linalg.norm(r_ref, axis=-1, keepdims=True)
    t_hat = np.cross(c_hat, r_hat)
    return r_hat, t_hat, c_hat


def project_ric(dr: np.ndarray, frames: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    r_hat, t_hat, c_hat = frames
    return np.stack([np.einsum("ij,ij->i", dr, r_hat),
                     np.einsum("ij,ij->i", dr, t_hat),
                     np.einsum("ij,ij->i", dr, c_hat)], axis=1)


# --------------------------------------------------------------------------- #
# Model hierarchy on dense substep grids landing exactly on snapshot epochs.
# --------------------------------------------------------------------------- #
def dense_grid(t_s: np.ndarray, nsub: int) -> np.ndarray:
    n_int = len(t_s) - 1
    k = np.arange(n_int * nsub + 1)
    return t_s[0] + k * (CADENCE_S / nsub)


def propagate_models(ref: dict, nsub: int, mu: float = MU_EARTH_KM3S2,
                     betas: tuple[float, ...] = ()) -> dict:
    """Propagate M1/M2 (+mu variants) and optionally M3(beta) from x(t0)."""
    t_dense = dense_grid(ref["t_s"], nsub)
    idx = np.arange(0, len(t_dense), nsub)
    x0 = np.concatenate([ref["r"][0], ref["v"][0]])
    out = {}
    out["M1"] = rk4_propagate(j2_rhs(mu, 0.0), t_dense, x0)[idx]
    out["M2"] = rk4_propagate(j2_rhs(mu, J2_EARTH), t_dense, x0)[idx]
    for b in betas:
        out[f"M3_beta{int(b)}"] = propagate_3d_rk4_drag(
            ref["r"][0], ref["v"][0], mu, t_dense, j2=J2_EARTH, beta=b,
            omega_atm=OMEGA_EARTH_RAD_S,
        )[idx]
    return out


# --------------------------------------------------------------------------- #
# Residual metrics.
# --------------------------------------------------------------------------- #
def rms(x: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


def residual_metrics(name: str, states: np.ndarray, ref: dict) -> dict:
    frames = ric_frames(ref["r"], ref["v"])
    dr = states[:, :3] - ref["r"]
    dv = states[:, 3:] - ref["v"]
    ric = project_ric(dr, frames)
    rad, tra, cro = ric[:, 0], ric[:, 1], ric[:, 2]
    nrm = np.linalg.norm(dr, axis=1)
    t_days = ref["t_s"] / 86400.0

    coef1 = np.polyfit(t_days, tra, 1)
    det1 = tra - np.polyval(coef1, t_days)
    coef2 = np.polyfit(t_days, tra, 2)
    det2 = tra - np.polyval(coef2, t_days)

    return {
        "model": name,
        "rms_ric_km": {"radial": rms(rad), "transverse": rms(tra), "cross": rms(cro)},
        "bias_ric_km": {"radial": float(np.mean(rad)), "transverse": float(np.mean(tra)), "cross": float(np.mean(cro))},
        "rms_3d_km": rms(nrm),
        "rms_dv_kms": rms(np.linalg.norm(dv, axis=1)),
        "max_res_km": float(np.max(nrm)),
        "residual_end_km": float(nrm[-1]),
        "residual_row10_km": float(nrm[min(10, len(nrm) - 1)]),
        "note_row0_zero": "row-0 residual is identically zero by declared initialization",
        "along_track_trend": {
            "c0_km": float(coef1[1]), "c1_km_per_day": float(coef1[0]),
            "c2_km_per_day2_quadratic_fit": float(coef2[0]),
            "detrended_linear_rms_km": rms(det1),
            "secular_ratio_c1T_over_detRMS": float(abs(coef1[0]) * (t_days[-1] - t_days[0]) / max(rms(det1), 1e-15)),
        },
        "equivalent_time_err_s": {
            "end": float(tra[-1] / np.linalg.norm(ref["v"][-1])),
            "max_abs": float(np.max(np.abs(tra)) / np.mean(np.linalg.norm(ref["v"], axis=1))),
        },
        "ric_end_km": {"radial": float(rad[-1]), "transverse": float(tra[-1]), "cross": float(cro[-1])},
    }


def jump_detector(ref: dict, states_m2: np.ndarray) -> dict:
    """Pre-registered contamination flag: second difference of the M2
    along-track residual between consecutive snapshots (> JUMP_THRESH_M)."""
    frames = ric_frames(ref["r"], ref["v"])
    tra = project_ric(states_m2[:, :3] - ref["r"], frames)[:, 1]
    d2 = np.abs(np.diff(tra, 2))
    worst = float(np.max(d2)) if len(d2) else 0.0
    flagged = int(np.sum(d2 > JUMP_THRESH_M / 1000.0))
    return {
        "threshold_m": JUMP_THRESH_M,
        "max_second_difference_m": worst * 1000.0,
        "flagged_samples": flagged,
        "window_rejected": flagged > 0,
    }


def block_bootstrap_skill(a_res: np.ndarray, b_res: np.ndarray) -> dict:
    """Seeded block bootstrap of skill = 1 - RMS(a)/RMS(b) over 1-day blocks."""
    n = len(a_res)
    block = BOOTSTRAP_BLOCKS_DAY
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    nblocks = int(np.ceil(n / block))
    skills = np.empty(BOOTSTRAP_B)
    for i in range(BOOTSTRAP_B):
        idx = np.concatenate(
            [start + np.arange(min(block, n - start))
             for start in rng.integers(0, n - block + 1, size=nblocks)]
        )
        ra = rms(a_res[idx])
        rb = rms(b_res[idx])
        skills[i] = 1.0 - ra / rb
    return {
        "skill_point": float(1.0 - rms(a_res) / rms(b_res)),
        "ci_low": float(np.percentile(skills, 2.5)),
        "ci_high": float(np.percentile(skills, 97.5)),
    }


# --------------------------------------------------------------------------- #
# Error-budget terms (all bounded BEFORE interpreting the remainder).
# --------------------------------------------------------------------------- #
def hermite_eval(t0: np.ndarray, x0: np.ndarray, tq: np.ndarray) -> np.ndarray:
    """Piecewise cubic Hermite using published positions AND velocities."""
    i = np.clip(np.searchsorted(t0, tq) - 1, 0, len(t0) - 2)
    h = (t0[i + 1] - t0[i])[:, None]  # column vector: avoid (Q,1)*(Q,) outer product
    u = ((tq - t0[i]) / h[:, 0])[:, None]
    p0, p1 = x0[i, :3], x0[i + 1, :3]
    m0, m1 = x0[i, 3:], x0[i + 1, 3:]
    return (
        (2 * u**3 - 3 * u**2 + 1) * p0
        + (u**3 - 2 * u**2 + u) * h * m0
        + (-2 * u**3 + 3 * u**2) * p1
        + (u**3 - u**2) * h * m1
    )


def interpolation_study(ref: dict) -> dict:
    """Cadence sensitivity of the REFERENCE itself (never feeds headlines)."""
    out = {}
    rn = np.linalg.norm(ref["r"], axis=1)
    n_mean = mean_motion(float(np.mean(rn)))
    for stride in (5, 10):
        idx = np.arange(0, ref["rows"], stride)
        held = np.array([i for i in range(1, ref["rows"] - 1) if i not in set(idx.tolist())])
        sub = np.concatenate([ref["r"][idx], ref["v"][idx]], axis=1)
        rec = hermite_eval(ref["t_s"][idx], sub, ref["t_s"][held])
        err = np.linalg.norm(rec - ref["r"][held], axis=1)
        h_cad = CADENCE_S * stride
        bound = float(np.mean(rn) * (n_mean * h_cad) ** 4 / 384.0)
        out[f"stride_{stride}"] = {
            "cadence_s": h_cad,
            "max_error_km": float(np.max(err)),
            "analytic_bound_km": bound,
            "within_bound": bool(np.max(err) <= 2.0 * max(bound, 1e-12)),
        }
    return out


# --------------------------------------------------------------------------- #
# Figures (one claim each, deterministic Agg output).
# --------------------------------------------------------------------------- #
def make_figures(payload: dict, series: dict, figdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    t_days = series["t_days"]

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for key, lbl, color in (("M1", "M1 two-body", "tab:red"),
                            ("M2", "M2 +J2", "tab:green"),
                            ("M3_beta100", "M3 +J2+drag (beta=100)", "tab:blue")):
        ax.semilogy(t_days, series["res_norm"][key], label=lbl, color=color, lw=1.2)
    ax.set_xlabel("time since epoch [days]")
    ax.set_ylabel("|r_model - r_ref| [km]")
    ax.set_title("Residual hierarchy vs Horizons ISS reference (3 d)")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "f1_residual_hierarchy.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ric = series["ric"]["M2"]
    ax.plot(t_days, ric[:, 0], label="radial", lw=1.0)
    ax.plot(t_days, ric[:, 1], label="transverse (in-track)", lw=1.0)
    ax.plot(t_days, ric[:, 2], label="cross-track", lw=1.0)
    ax.set_xlabel("time since epoch [days]")
    ax.set_ylabel("M2 residual component [km]")
    ax.set_title("M2 (+J2) residual structure in reference-built RTN frame")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "f2_ric_components.png", dpi=150)
    plt.close(fig)

    conv = payload["error_budget"]["E_integration"]["ladder"]
    dts = [c["dt_s"] for c in conv]
    errs = [c["rms_3d_km"] for c in conv]
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.loglog(dts, errs, "o-", label="M2 RMS3D(dt)")
    ord_ = payload["error_budget"]["E_integration"]["measured_order_selfconvergence"]
    ax.set_xlabel("integration step dt [s]")
    ax.set_ylabel("RMS 3-D residual vs reference [km]")
    ax.set_title(f"dt-convergence ladder (self-convergence order p = {ord_:.2f})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figdir / "f3_dt_convergence.png", dpi=150)
    plt.close(fig)

    interp = payload["error_budget"]["E_interpolation"]["strides"]
    cad = [interp[k]["cadence_s"] for k in sorted(interp)]
    meas = [interp[k]["max_error_km"] for k in sorted(interp)]
    bound = [interp[k]["analytic_bound_km"] for k in sorted(interp)]
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.loglog(cad, meas, "s-", label="Hermite reconstruction error (measured)")
    ax.loglog(cad, bound, "--", label="analytic bound A(w h)^4/384")
    ax.set_xlabel("subsampled cadence h [s]")
    ax.set_ylabel("position error [km]")
    ax.set_title("Reference-snapshot interpolation sensitivity (diagnostic only)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figdir / "f4_interpolation_bound.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    env_days = series["t_days"] + REVISION_OFFSET_DAYS
    ax.fill_between(series["t_days"], 0.0, ENVELOPE_KM_PER_DAY * env_days,
                    alpha=0.18, color="gray", label="declared reference envelope (3 km/day)")
    for key, lbl, color in (("M2", "M2 +J2", "tab:green"),
                            ("M3_beta50", "M3 beta=50", "tab:cyan"),
                            ("M3_beta100", "M3 beta=100", "tab:blue"),
                            ("M3_beta200", "M3 beta=200", "tab:purple"),
                            ("M3_beta400", "M3 beta=400", "tab:orange")):
        ax.plot(series["t_days"], series["res_trans"][key], label=lbl, lw=1.1, color=color)
    ax.plot(series["t_days"], series["res_trans"]["M2_muWGS72"], "--", lw=1.1,
            color="tab:green", label="M2 with SGP4-heritage GM")
    ax.set_xlabel("time since epoch [days]")
    ax.set_ylabel("transverse (in-track) residual [km]")
    ax.set_title("Drag band + GM-convention sensitivity vs declared reference envelope")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "f5_beta_gm_sensitivity.png", dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main experiment.
# --------------------------------------------------------------------------- #
def run() -> dict:
    ref = load_reference()

    # ---- 1. Integration-error ladder FIRST (before residual inspection) -- #
    ladder = []
    prev_states = None
    for nsub in LADDER_NSUB:
        m2 = propagate_models(ref, nsub)["M2"]
        entry = {"nsub": nsub, "dt_s": CADENCE_S / nsub,
                 "rms_3d_km": rms(np.linalg.norm(m2[:, :3] - ref["r"], axis=1))}
        if prev_states is not None:
            entry["max_delta_vs_prev_km"] = float(np.max(np.linalg.norm(m2[:, :3] - prev_states[:, :3], axis=1)))
        ladder.append(entry)
        prev_states = m2
    # Convergence order from SELF-convergence deltas (Exp 009 doctrine: total
    # residuals plateau at the model-mismatch floor and carry no order info).
    deltas = [e["max_delta_vs_prev_km"] for e in ladder[1:]]
    orders = [
        float(np.log2(dp / dc))
        for dp, dc in zip(deltas, deltas[1:])
        if dp > 0 and dc > 0
    ]
    measured_order = float(np.median(orders)) if orders else float("nan")
    richardson_bound = (
        float(deltas[-1] / (2.0**measured_order - 1.0))
        if orders and np.isfinite(measured_order) and measured_order > 1.0
        else float("nan")
    )

    # ---- 2. Headline models at the pre-registered converged grid ---------- #
    states = propagate_models(ref, HEADLINE_NSUB, betas=BETA_BAND_KG_M2)
    states_mu_wgs72 = propagate_models(ref, HEADLINE_NSUB, mu=MU_SGP4_WGS72)
    states_mu_de440 = propagate_models(ref, HEADLINE_NSUB, mu=MU_DE440)

    # Contamination gate BEFORE any comparison is interpreted.
    jump = jump_detector(ref, states["M2"])
    if jump["window_rejected"]:
        raise RuntimeError(
            f"pre-registered rejection triggered: {jump['flagged_samples']} samples "
            f"exceed {JUMP_THRESH_M} m second-difference (maneuver/TLE handover?)"
        )

    metrics = {}
    for name, st in states.items():
        metrics[name] = residual_metrics(name, st, ref)
    metrics["M2_muWGS72"] = residual_metrics("M2_muWGS72", states_mu_wgs72["M2"], ref)
    metrics["M2_muDE440"] = residual_metrics("M2_muDE440", states_mu_de440["M2"], ref)

    # ---- 3. Decision rule: bootstrap CI + declared envelope --------------- #
    frames = ric_frames(ref["r"], ref["v"])
    tra = {k: project_ric(st[:, :3] - ref["r"], frames)[:, 1] for k, st in states.items()}
    boot_m2_m1 = block_bootstrap_skill(tra["M2"], tra["M1"])
    boot_m3_m2 = block_bootstrap_skill(tra["M3_beta100"], tra["M2"])
    envelope_end_km = ENVELOPE_KM_PER_DAY * (ref["t_s"][-1] / 86400.0 + REVISION_OFFSET_DAYS)

    def decide(boot: dict, rms_new: float, rms_old: float) -> str:
        skill = boot["skill_point"]
        improvement_km = rms_old - rms_new
        if skill <= 0.0:
            return "worse-than-previous-model"
        if boot["ci_low"] > 0.0 and improvement_km > envelope_end_km:
            return "improvement-exceeds-CI-and-envelope"
        return "indistinguishable-given-reference-uncertainty"

    decision_m2_vs_m1 = decide(boot_m2_m1, metrics["M2"]["rms_3d_km"], metrics["M1"]["rms_3d_km"])
    decision_m3_vs_m2 = decide(boot_m3_m2, metrics["M3_beta100"]["rms_3d_km"], metrics["M2"]["rms_3d_km"])

    # ---- 4. Remaining budget terms ----------------------------------------- #
    interp = interpolation_study(ref)

    init_seed = coe_to_rv_eci(
        float(ref["elements_t0"]["a"]), float(ref["elements_t0"]["e"]),
        float(ref["elements_t0"]["inc"]), float(ref["elements_t0"]["Omega"]),
        float(ref["elements_t0"]["omega"]), float(ref["elements_t0"]["nu"]),
    )
    x0_direct = np.concatenate([ref["r"][0], ref["v"][0]])
    x0_roundtrip = np.concatenate(init_seed)
    m2_rt = rk4_propagate(
        j2_rhs(MU_EARTH_KM3S2, J2_EARTH), dense_grid(ref["t_s"], HEADLINE_NSUB), x0_roundtrip
    )[::HEADLINE_NSUB]
    e_init = float(np.max(np.linalg.norm(m2_rt[:, :3] - states["M2"][:, :3], axis=1)))

    tau_unc_s = 43.2e-6  # half print quantum of the 9-decimal JD tags
    kappa_kms = float(np.mean(np.linalg.norm(ref["v"], axis=1)))
    e_time_km = kappa_kms * tau_unc_s

    frame_diag = project_ric(states["M2"][:, :3] - ref["r"], ric_frames(states["M2"][:, :3], states["M2"][:, 3:]))
    frame_diag_rms = rms(frame_diag[:, 1] - tra["M2"])

    # ---- 5. Assemble payload ----------------------------------------------- #
    n_mean_t0 = mean_motion(float(ref["elements_t0"]["a"]))
    v_mean = float(np.mean(np.linalg.norm(ref["v"], axis=1)))
    beta_c1 = {f"beta{int(b)}": metrics[f"M3_beta{int(b)}"]["along_track_trend"]["c1_km_per_day"]
               for b in BETA_BAND_KG_M2}
    c1_wgs72 = metrics["M2_muWGS72"]["along_track_trend"]["c1_km_per_day"]
    # Reference-relative mean semi-major-axis decay implied by the GM-corrected
    # M2 transverse drift rate: ds/dt = a*dn  ->  da/day = -(2/3)*c1/(86400*n).
    implied_da_day_km = -(2.0 / 3.0) * (c1_wgs72 / 86400.0) / n_mean_t0

    findings = [
        f"M1-vs-reference: RMS3D {metrics['M1']['rms_3d_km']:.1f} km dominated by a "
        f"{metrics['M1']['along_track_trend']['c1_km_per_day']:.1f} km/day transverse drift "
        f"(constant mean-motion mismatch of the purely-Keplerian model: the mean-vs-osculating "
        f"a offset at the initialization epoch plus absent J2 secular rates); adding J2 removes skill "
        f"{boot_m2_m1['skill_point']:.4f} of the residual (CI [{boot_m2_m1['ci_low']:.4f}, "
        f"{boot_m2_m1['ci_high']:.4f}]) -> {decision_m2_vs_m1}.",
        f"M3(primary beta=100) WORSENS agreement vs M2: skill {boot_m3_m2['skill_point']:.4f} "
        f"(CI [{boot_m3_m2['ci_low']:.4f}, {boot_m3_m2['ci_high']:.4f}], excludes zero on the "
        f"negative side) -> {decision_m3_vs_m2}. Reported verbatim per anti-overfitting doctrine; "
        "no atmospheric or beta parameter was retuned.",
        f"[INFERENCE] The beta band responds monotonically (transverse trend c1: "
        f"{json.dumps({k: round(v, 2) for k, v in beta_c1.items()})} km/day), crossing zero only at the "
        f"band edge (beta=400 -> RMS3D {metrics['M3_beta400']['rms_3d_km']:.2f} km). The drag signature is "
        f"detectable, but the compatible effective beta lies AT OR BEYOND the pre-declared band edge; "
        f"extending the band or picking the nicest member post hoc would be tuning and is deliberately "
        f"not done. Any refined-beta study must be a separately declared experiment.",
        f"[INFERENCE] M2's leftover secular drift ({metrics['M2']['along_track_trend']['c1_km_per_day']:.2f} km/day; "
        f"{c1_wgs72:.2f} km/day under SGP4-heritage GM) implies the reference mean semilatus "
        f"evolves as if decaying at {implied_da_day_km * 1e3 * -1:.1f} m/day relative to our non-decaying model "
        "(consistent with TLE-B*-absorbed drag plus density-weather; not independently separable).",
        f"SGP4-heritage-GM sensitivity shifts M2 RMS3D by "
        f"{metrics['M2_muWGS72']['rms_3d_km'] - metrics['M2']['rms_3d_km']:+.3f} km "
        "(constants/convention term bounded, canon retained).",
        f"Integration bound (NSUB 16 vs 8): {deltas[-1] * 1e3:.1f} m "
        f"(self-convergence order p = {measured_order:.2f}); reference-interpolation diagnostic bounded by "
        f"{max(v['max_error_km'] for v in interp.values()):.2e} km at worst subsampled cadence and never "
        "enters headline numbers (exact-grid alignment). Contamination gate: no sample exceeded the "
        f"{JUMP_THRESH_M:.0f} m second-difference threshold (max {jump['max_second_difference_m']:.1f} m).",
    ]
    limitations = [
        "Reference is TLE/SGP4-provenance (disclosed in the object sheet); its own error grows "
        "~1-3 km/day from the revision date per JPL documentation, so residuals measure "
        "model-vs-reference divergence jointly with reference uncertainty (never separated).",
        "Constant-beta drag with a static exponential atmosphere idealizes ISS attitude-dependent "
        "area and thermospheric weather; no drag-improvement claim is made unless it exceeds both "
        "the bootstrap CI and the declared envelope (decision rule applied automatically).",
        "Unmodelled physics remaining: lunisolar third-body, tesseral/zonal harmonics beyond J2, SRP, "
        "relativity -- jointly attributed with reference uncertainty to the remainder.",
        "TEME->ICRF conversion inside Horizons is undocumented; its effect is inseparable from the "
        "reference's own solution and therefore part of the remainder.",
        "Single 3-day window immediately after one revision date; conclusions are window-scoped.",
    ]

    payload = {
        "provenance": {
            "reference_source": "NASA JPL Horizons API (official REST)",
            "target": "International Space Station (-125544)",
            "reference_provenance": "tle_sgp4_disclosed",
            "trajectory_revision_date": ref["revision_date"],
            "window_tdb": {"start": ref["cal_first"], "rows": ref["rows"],
                           "cadence_s": CADENCE_S, "jdtdb_last": float(ref["jd"][-1])},
            "query_params": ref["manifest"]["acquisition"]["query_params"],
            "header_metadata": ref["header"],
            "reference_snapshot_files": ref["manifest"]["snapshot"]["files"],
            "offline_deterministic": True,
            "network_access_in_experiment": False,
        },
        "constants": {
            "mu_km3_s2": MU_EARTH_KM3S2,
            "j2": J2_EARTH,
            "R_eq_km": R_EARTH_KM,
            "omega_atm_rad_s": OMEGA_EARTH_RAD_S,
            "provenance": "src/lab_utils/orbits.py canon (IAU 2015 nominal GM, WGS-84); "
                          "atmosphere = Exp 010 declared Vallado layers verbatim",
        },
        "models": {
            "M1": "two-body point mass",
            "M2": "two-body + J2",
            "M3": "two-body + J2 + exponential-atmosphere drag (co-rotating), "
                  f"beta band {list(BETA_BAND_KG_M2)} kg/m^2, primary beta = {BETA_PRIMARY}",
            "initialization": "x(t0) taken verbatim from the reference state at the first "
                              "snapshot epoch (declared initialization; no parameter fitted "
                              "against any t > t0 data)",
            "headline_nsub": HEADLINE_NSUB,
        },
        "residual_metrics": metrics,
        "decision_rule": {
            "bootstrap": {"blocks_samples": BOOTSTRAP_BLOCKS_DAY, "draws": BOOTSTRAP_B,
                          "seed": BOOTSTRAP_SEED},
            "M2_vs_M1": {**boot_m2_m1, "verdict": decision_m2_vs_m1},
            "M3_vs_M2": {**boot_m3_m2, "verdict": decision_m3_vs_m2},
            "envelope_km_at_window_end": envelope_end_km,
            "envelope_rate_km_per_day": ENVELOPE_KM_PER_DAY,
        },
        "contamination_gate": jump,
        "error_budget": {
            "E_integration": {
                "ladder": ladder,
                "measured_order_selfconvergence": measured_order,
                "richardson_bound_km": richardson_bound,
                "bound_km_max_ric_diff_nsub16_vs_8": deltas[-1] if deltas else None,
                "headline_nsub": HEADLINE_NSUB,
            },
            "E_interpolation": {"strides": interp,
                                "role": "diagnostic only; headline residuals use exact-grid alignment"},
            "E_timeshift": {
                "kappa_km_per_s": kappa_kms,
                "tau_uncertainty_s": tau_unc_s,
                "charge_km": e_time_km,
                "note": "single TDB scale end-to-end; TT-TDB <= 2 ms never mixed",
            },
            "E_initialization": {
                "element_roundtrip_max_effect_km": e_init,
                "method": "M2 re-seeded via rv_to_coe_eci -> coe_to_rv_eci round trip",
            },
            "E_constants": {
                "mu_sgp4_wgs72_variant_rms3d_km": metrics["M2_muWGS72"]["rms_3d_km"],
                "mu_de440_variant_rms3d_km": metrics["M2_muDE440"]["rms_3d_km"],
                "mu_canon_rms3d_km": metrics["M2"]["rms_3d_km"],
            },
            "E_frame_origin_diagnostic_transverse_rms_km": frame_diag_rms,
        },
        "findings": findings,
        "limitations": limitations,
        "code_sha256": code_hashes(),
    }

    results_dir = EXP_DIR / "results"
    from lab_utils.results import save_json_result

    save_json_result(
        str(results_dir / "results.json"), payload,
        name="jplValidation-013",
        description="ISS (-125544) Horizons ICRF/TDB geometric-state validation of the "
                    "two-body / +J2 / +J2+drag hierarchy with decomposition budget",
    )

    series = {
        "t_days": ref["t_s"] / 86400.0,
        "res_norm": {k: np.linalg.norm(st[:, :3] - ref["r"], axis=1) for k, st in states.items()},
        "res_trans": {
            **{k: tra[k] for k in tra},
            "M2_muWGS72": project_ric(states_mu_wgs72["M2"][:, :3] - ref["r"], frames)[:, 1],
        },
        "ric": {"M2": project_ric(states["M2"][:, :3] - ref["r"], frames)},
    }
    make_figures(payload, series, results_dir / "figures")
    return payload


if __name__ == "__main__":
    p = run()
    print("== Exp 013 complete ==")
    print(f"M1 RMS3D {p['residual_metrics']['M1']['rms_3d_km']:.3f} km | "
          f"M2 {p['residual_metrics']['M2']['rms_3d_km']:.3f} km | "
          f"M3(beta=100) {p['residual_metrics']['M3_beta100']['rms_3d_km']:.3f} km")
    print("decisions:", json.dumps({k: v["verdict"] for k, v in p["decision_rule"].items() if isinstance(v, dict) and "verdict" in v}))
