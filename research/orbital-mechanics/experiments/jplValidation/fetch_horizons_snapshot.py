"""One-time online acquisition of the pinned JPL Horizons reference snapshot.

Exp 013 (JPL ephemeris validation) separates acquisition from analysis:
this script runs ONCE while online; every later numerical run consumes only
the checksum-pinned snapshot under ``reference/`` and never touches the
network (deterministic-offline doctrine).

Politeness policy (lab charter): official API, strictly serial requests,
>= 3 s spacing, single-digit request count total, no retries beyond the one
documented TIME_DIGITS fallback, refuse-to-overwrite idempotence guard.

Usage:
    python fetch_horizons_snapshot.py            # performs acquisition
    python fetch_horizons_snapshot.py --check    # verify existing artifacts
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

# Target identity is hard-pinned (unique Horizons record, verified against
# the returned object sheet before any vector data is accepted).
ISS_COMMAND = "'-125544'"
ISS_EXPECTED_NAME_TOKEN = "International Space Station"

# Window starts immediately after the trajectory revision date disclosed on
# the object sheet (TLE-based predicts degrade beyond a few days); cadence
# fixed a priori. All epochs TDB (uniform dynamical scale, leap-free).
VECTOR_PARAMS = {
    "format": "text",
    "MAKE_EPHEM": "'YES'",
    "OBJ_DATA": "'NO'",
    "EPHEM_TYPE": "'VECTORS'",
    "COMMAND": ISS_COMMAND,
    "CENTER": "'500@399'",
    "REF_PLANE": "'FRAME'",
    "REF_SYSTEM": "'ICRF'",
    "VEC_TABLE": "'2'",
    "VEC_CORR": "'NONE'",
    "OUT_UNITS": "'KM-S'",
    "TIME_TYPE": "'TDB'",
    "CSV_FORMAT": "'YES'",
    "CAL_FORMAT": "'BOTH'",
    "START_TIME": "'2026-08-24 00:00'",
    "STOP_TIME": "'2026-08-27 00:00'",
    "STEP_SIZE": "'2m'",
}
EXPECTED_ROWS = 2161  # 3 days x (1440/2) intervals, inclusive endpoints
REQUEST_SPACING_S = 3.0


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
        "eop_file": _extract(text, r"EOP file\s*:\s*(\S+)"),
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
    # JD tags print at 9 decimals -> each epoch carries up to +-43 us print
    # quantization; consecutive-difference deviations up to ~1e-4 s are pure
    # rounding, not non-uniform stepping.
    diffs = {round((b - a) * 86400.0, 6) for a, b in zip(jds, jds[1:])}
    if any(abs(d - 120.0) > 2e-4 for d in diffs):
        raise RuntimeError(f"non-uniform epoch spacing detected: {sorted(diffs)[:5]}...")
    meta = _header_metadata(text)
    if meta["reference_frame"] != "ICRF":
        raise RuntimeError(f"frame echo '{meta['reference_frame']}' != ICRF")
    if meta["target_body_name"] is None or "(-125544)" not in meta["target_body_name"]:
        raise RuntimeError(f"unexpected target echo: {meta['target_body_name']}")
    if "GEOMETRIC" not in (meta["output_type"] or ""):
        raise RuntimeError(f"states are not geometric: {meta['output_type']}")
    return {"rows": len(rows), "jd_first": jds[0], "jd_last": jds[-1], "header": meta}


def acquire() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    vectors_path = REFERENCE_DIR / (
        "horizons_-125544_iss_vectors_2026-08-24_to_2026-08-27_tdb_2min.txt"
    )
    objdata_path = REFERENCE_DIR / "horizons_-125544_iss_objdata.txt"
    for p in (vectors_path, objdata_path):
        if p.exists():
            raise SystemExit(
                f"refusing to overwrite existing snapshot artifact: {p.name} "
                "(acquisition is one-time; delete artifacts deliberately to re-acquire)"
            )

    # --- Request 1: object-data provenance sheet -------------------------- #
    t0 = datetime.now(timezone.utc).isoformat()
    status_obj, raw_obj = _get(
        {"format": "text", "COMMAND": ISS_COMMAND, "MAKE_EPHEM": "'NO'", "OBJ_DATA": "'YES'"}
    )
    text_obj = raw_obj.decode("utf-8")
    if status_obj != 200 or ISS_EXPECTED_NAME_TOKEN not in text_obj:
        raise RuntimeError(f"object-sheet request failed (HTTP {status_obj})")
    revised = _extract(text_obj, r"Revised:\s*([A-Za-z]{3}\s+\d{1,2},\s+\d{4})")
    tle_disclosed = "TLE-based" in text_obj
    time.sleep(REQUEST_SPACING_S)

    # --- Request 2: geometric ICRF state vectors --------------------------- #
    try:
        status_vec, raw_vec = _get(VECTOR_PARAMS)
    except urllib.error.HTTPError as e:  # pragma: no cover - network-dependent
        if e.code != 400:
            raise
        # Single documented fallback: FRACSEC -> SECONDS output precision.
        fb = dict(VECTOR_PARAMS)
        fb["TIME_DIGITS"] = "'SECONDS'"
        status_vec, raw_vec = _get(fb)
        VECTOR_PARAMS["TIME_DIGITS"] = "'SECONDS'"
    else:
        VECTOR_PARAMS["TIME_DIGITS"] = "'FRACSEC'"
    if status_vec != 200:
        raise RuntimeError(f"vector request failed (HTTP {status_vec})")
    # format=text carries Horizons-level errors in-band; the structural
    # validation below catches them (no $$SOE table -> RuntimeError).
    text_vec = raw_vec.decode("utf-8")
    validation = _validate_vector_response(text_vec)
    t1 = datetime.now(timezone.utc).isoformat()

    # --- Pin artifacts ------------------------------------------------------ #
    vectors_path.write_bytes(raw_vec)
    objdata_path.write_bytes(raw_obj)

    api_version = _extract(text_vec, r"API VERSION:\s*(\S+)")
    manifest = {
        "manifest_version": 1,
        "experiment": "jplValidation",
        "acquisition": {
            "service": "NASA JPL Horizons",
            "base_url": BASE_URL,
            "http_method": "GET",
            "request_order": ["objdata", "vectors"],
            "retrieved_utc": {"objdata": t0, "vectors_complete": t1},
            "http_status": {"objdata": status_obj, "vectors": status_vec},
            "query_params": {
                "objdata": {
                    "format": "text",
                    "COMMAND": ISS_COMMAND,
                    "MAKE_EPHEM": "'NO'",
                    "OBJ_DATA": "'YES'",
                },
                "vectors": dict(sorted(VECTOR_PARAMS.items())),
            },
            "response_sha256": {
                "objdata": _sha256(raw_obj),
                "vectors": _sha256(raw_vec),
            },
        },
        "snapshot": {
            "files": {
                p.name: {"bytes": n, "sha256": _sha256(p.read_bytes())}
                for p, n in ((vectors_path, len(raw_vec)), (objdata_path, len(raw_obj)))
            }
        },
        "provenance": {
            "horizons_api_version": api_version,
            "trajectory_revised_date": revised,
            "trajectory_tle_based_disclosed": tle_disclosed,
            "header_metadata": validation["header"],
            "attribution": "Source: NASA/JPL Solar System Dynamics Horizons system "
            "(public API, human-use-rate access; see LICENSE-NOTES.md)",
            "acquisition_script": "fetch_horizons_snapshot.py",
            "single_fetch_policy": "one-shot online acquisition; all numerical runs "
            "consume these committed snapshots only and never access the network",
        },
        "structure": {
            "rows": validation["rows"],
            "jdtdb_first": validation["jd_first"],
            "jdtdb_last": validation["jd_last"],
            "cadence_s": 120.0,
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("acquisition OK")
    print("  vectors:", vectors_path.name, f"({len(raw_vec)} bytes)")
    print("  objdata:", objdata_path.name, f"({len(raw_obj)} bytes)")
    print("  rows:", validation["rows"], "revised:", revised, "TLE-based:", tle_disclosed)
    print("  frame:", validation["header"]["reference_frame"],
          "| target:", validation["header"]["target_body_name"])
    print("  manifest:", MANIFEST_PATH.name)


def check() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    ok = True
    for name, rec in manifest["snapshot"]["files"].items():
        p = REFERENCE_DIR / name
        data = p.read_bytes()
        good = len(data) == rec["bytes"] and _sha256(data) == rec["sha256"]
        ok &= good
        print(("OK  " if good else "FAIL"), name, len(data), rec["bytes"])
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    if "--check" in sys.argv:
        check()
    else:
        acquire()
