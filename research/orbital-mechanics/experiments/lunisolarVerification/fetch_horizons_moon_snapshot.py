"""One-time online acquisition of the pinned JPL Horizons Moon reference snapshot.

Exp 017 (lunisolar upper-bound verification) follows the Exp 013/014 acquisition
doctrine: this script runs ONCE while online; every later numerical run consumes
only the checksum-pinned snapshot under ``reference/`` and never touches the
network (deterministic-offline doctrine).

Role of the snapshot: validation-only oracle for the lab's closed-form Lunisolar
RAAN upper bound. It is never a runtime input; byte-reproduction means the
committed files, not re-querying Horizons.

Politeness policy (lab charter): official API, strictly serial requests,
>= 3 s spacing, single-digit request count total, refuse-to-overwrite
idempotence guard. This file follows the eclipseTiming fetch pattern.

Usage:
    python fetch_horizons_moon_snapshot.py            # performs acquisition
    python fetch_horizons_moon_snapshot.py --check    # verify existing artifacts
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
REFERENCE_DIR = Path(__file__).resolve().parent / "reference"
MANIFEST_PATH = REFERENCE_DIR / "MANIFEST.json"

# Moon (301) as seen from Earth body center (399). ICRF, geometric, KM-S, TDB.
# This is the geocentric Moon vector - identical acquisition pattern to the
# existing Sun snapshot (eclipseTiming) and ISS snapshot (jplValidation).
# Cadence matches the Sun snapshot: 1 day, full 2026 calendar year.
MOON_COMMAND = "'301'"
VECTOR_PARAMS = {
    "format": "text",
    "MAKE_EPHEM": "'YES'",
    "OBJ_DATA": "'NO'",
    "EPHEM_TYPE": "'VECTORS'",
    "COMMAND": MOON_COMMAND,
    "CENTER": "'500@399'",
    "REF_PLANE": "'FRAME'",
    "REF_SYSTEM": "'ICRF'",
    "VEC_TABLE": "'2'",
    "VEC_CORR": "'NONE'",
    "OUT_UNITS": "'KM-S'",
    "TIME_TYPE": "'TDB'",
    "CSV_FORMAT": "'YES'",
    "CAL_FORMAT": "'BOTH'",
    # Full year 2026 at daily cadence: covers the lunar synodic month (29.53 d),
    # anomalistic month (27.55 d), draconic month (27.21 d), and nodal
    # regression (18.6 yr - sub-decade). One full year of daily snapshots
    # integrates the dominant Lunisolar perturbations at LEO SSO.
    "START_TIME": "'2026-01-01 00:00'",
    "STOP_TIME": "'2027-01-01 00:00'",
    "STEP_SIZE": "'1d'",
}
EXPECTED_ROWS = 366  # inclusive endpoints over a 365-day span
REQUEST_SPACING_S = 3.0

IDENTITY_TOKENS = {
    "target": "Moon",
    "center": "Earth (399)",
}

# Physical plausibility band for geocentric Moon distance.
# Mean Earth-Moon distance 384400 km; perihelion ~ 356500; aphelion ~ 406700.
# Allow +/- 5000 km margin for measurement rounding and mean-orbit variation.
MOON_DISTANCE_MIN_KM = 350000.0
MOON_DISTANCE_MAX_KM = 412000.0


def _get(params: dict) -> tuple[int, bytes]:
    url = BASE_URL + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=120) as resp:
        return resp.status, resp.read()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _extract(result_text: str, pattern: str) -> str | None:
    m = re.search(pattern, result_text)
    return m.group(1).strip() if m else None


def _header_metadata(text: str) -> dict:
    """Provenance lines extracted verbatim from the response banner."""
    return {
        "target_body_name": _extract(text, r"Target body name:\s*(.+)"),
        "center_body_name": _extract(text, r"Center body name:\s*(.+)"),
        "center_site_name": _extract(text, r"Center-site name:\s*(.+)"),
        "reference_frame": _extract(text, r"Reference frame\s*:\s*(\S+)"),
        "output_units": _extract(text, r"Output units\s*:\s*(\S+)"),
        "output_type": _extract(text, r"Output type\s*:\s*(.+)"),
        "start_time_echo": _extract(text, r"Start time\s*:\s*(.+)"),
        "stop_time_echo": _extract(text, r"Stop\s+time\s*:\s*(.+)"),
        "step_size_echo": _extract(text, r"Step-size\s*:\s*(.+)"),
    }


def _validate_vector_response(text: str) -> dict:
    """Structural validation BEFORE anything is pinned. Raises on failure."""
    lines = text.splitlines()
    soe = [i for i, ln in enumerate(lines) if ln.strip() == "$$SOE"]
    eoe = [i for i, ln in enumerate(lines) if ln.strip() == "$$EOE"]
    if len(soe) != 1 or len(eoe) != 1 or eoe[0] <= soe[0]:
        raise RuntimeError(f"SOE/EOE delimiters malformed: {soe} {eoe}")
    rows = lines[soe[0] + 1 : eoe[0]]
    if len(rows) != EXPECTED_ROWS:
        raise RuntimeError(f"row count {len(rows)} != expected {EXPECTED_ROWS}")
    jds = []
    for k, row in enumerate(rows):
        parts = [p.strip() for p in row.split(",")]
        # CSV columns: JDTDB, Calendar Date (TDB), X, Y, Z, VX, VY, VZ, ''(trailing)
        if len(parts) != 9 or parts[-1] != "":
            raise RuntimeError(f"row {k}: unexpected column structure ({len(parts)} fields)")
        jd = float(parts[0])
        for comp in parts[2:8]:
            float(comp)  # raises ValueError on unparsable scientific notation
        jds.append(jd)
    diffs = {round((b - a) * 86400.0, 6) for a, b in zip(jds, jds[1:])}
    if any(abs(d - 86400.0) > 2e-4 for d in diffs):
        raise RuntimeError(f"non-uniform epoch spacing detected: {sorted(diffs)[:5]}...")

    meta = _header_metadata(text)
    if meta["target_body_name"] is None or IDENTITY_TOKENS["target"] not in meta["target_body_name"]:
        raise RuntimeError(f"target identity gate failed: {meta['target_body_name']!r}")
    if meta["center_body_name"] is None or IDENTITY_TOKENS["center"] not in meta["center_body_name"]:
        raise RuntimeError(f"center identity gate failed: {meta['center_body_name']!r}")
    if meta["reference_frame"] != "ICRF":
        raise RuntimeError(f"reference-frame gate failed: {meta['reference_frame']!r}")
    if meta["output_units"] != "KM-S":
        raise RuntimeError(f"units gate failed: {meta['output_units']!r}")

    # Magnitude plausibility: geocentric Moon distance must sit within the
    # physical perihelion/aphelion band (+/- margin) across the whole year.
    dists = []
    for row in rows:
        parts = [p.strip() for p in row.split(",")]
        x, y, z = (float(c) for c in parts[2:5])
        dists.append((x * x + y * y + z * z) ** 0.5)
    dmin, dmax = min(dists), max(dists)
    if not (MOON_DISTANCE_MIN_KM < dmin < dmax < MOON_DISTANCE_MAX_KM):
        raise RuntimeError(
            f"Moon distance outside physical band "
            f"[{MOON_DISTANCE_MIN_KM}, {MOON_DISTANCE_MAX_KM}]: "
            f"[{dmin:.6e}, {dmax:.6e}]"
        )
    return {"header": meta, "rows": len(rows), "dist_min_km": dmin, "dist_max_km": dmax}


def _write_snapshot(raw: bytes, validation: dict, url: str, status: int) -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    snap_name = "horizons_moon_geocentric_vectors_2026_icrf_tdb_daily.txt"
    snap_path = REFERENCE_DIR / snap_name
    if snap_path.exists():
        raise RuntimeError(f"refuse-to-overwrite guard: {snap_path.name} already exists")
    manifest = {
        "schema": "lab.acquisition.manifest/v1",
        "acquired_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "acquisition": {
            "base_url": BASE_URL,
            "http_method": "GET",
            "http_status": {"vectors": status},
            "query_params": {"vectors": VECTOR_PARAMS},
            "request_order": ["vectors"],
            "request_spacing_s": REQUEST_SPACING_S,
            "response_sha256": {"vectors": _sha256(raw)},
            "response_bytes": {"vectors": len(raw)},
            "validation": validation,
        },
        "snapshot": {
            "files": {
                snap_name: {
                    "bytes": len(raw),
                    "sha256": _sha256(raw),
                    "role": (
                        "validation-only oracle for the lab's closed-form Lunisolar "
                        "RAAN upper-bound formula (Vallado Eq. 9-46 form); never a "
                        "runtime input; the lab's analytical Lunisolar model uses "
                        "mean-orbit constants (LUNAR_DISTANCE_KM=384400, LUNAR_INCLINATION_DEG=5.145) "
                        "for production calculations - this snapshot is the byte-pinned "
                        "JPL DE441 ephemeris against which the closed-form approximation "
                        "is independently verified"
                    ),
                }
            }
        },
    }
    snap_path.write_bytes(raw)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[fetch] pinned {snap_name}: {len(raw)} bytes sha256={_sha256(raw)[:16]}...")


def _verify_existing() -> int:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    ok = True
    for name, rec in manifest["snapshot"]["files"].items():
        p = REFERENCE_DIR / name
        data = p.read_bytes()
        digest = _sha256(data)
        match = digest == rec["sha256"] == manifest["acquisition"]["response_sha256"]["vectors"]
        ok &= match and len(data) == rec["bytes"]
        print(f"[check] {name}: {'OK' if match else 'MISMATCH'} sha256={digest[:16]}...")
    return 0 if ok else 1


def main() -> int:
    if "--check" in sys.argv[1:]:
        if not MANIFEST_PATH.exists():
            print("[check] no manifest present")
            return 1
        return _verify_existing()
    if MANIFEST_PATH.exists():
        print("[fetch] manifest already present; refusing to re-acquire (idempotence guard)")
        return _verify_existing()
    time.sleep(REQUEST_SPACING_S)
    status, raw = _get(VECTOR_PARAMS)
    if status != 200:
        raise RuntimeError(f"HTTP {status} from Horizons API")
    text = raw.decode("utf-8", errors="strict")
    validation = _validate_vector_response(text)
    url = BASE_URL + "?" + urllib.parse.urlencode(VECTOR_PARAMS)
    _write_snapshot(raw, validation, url, status)
    return 0


if __name__ == "__main__":
    sys.exit(main())