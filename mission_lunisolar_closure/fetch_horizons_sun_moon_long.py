"""Multi-year DE441 Sun + Moon acquisition for mission_lunisolar_closure.

Acquires 19 years (2026-01-01 -> 2045-01-01) of daily geocentric Sun and Moon
vectors from NASA/JPL Horizons (DE441 source), following the lab's standard
acquisition doctrine (Exp 013/014/017).

Politeness:
  - Single acquisition pass; the resulting byte-identical snapshot is
    committed to the repo under `reference/` with sha256 pins and
    `-text` gitattributes.
  - 5-year chunks per request to keep individual response sizes manageable.
  - 4-second sleep between sequential requests.
  - Never attempts CAPTCHAs / access controls / rate-limit evasion.

Outputs (committed under `reference/`):
  horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt (5 sub-files concatenated)
  horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt
  MANIFEST.json
  ACQUISITION_LOG.json
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
OUT_DIR = Path(r"R:\lab_scratch\mission_lunisolar_closure\ephemeris")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 19-year span: 2026-01-01 to 2045-01-01, 5-year chunks (last is 4 yr)
SPAN_START_YEAR = 2026
SPAN_END_YEAR = 2045
CHUNK_YEARS = 5
REQUEST_SPACING_S = 4.0

IDENTITY_TOKENS = {
    "sun_target": "Sun (10)",
    "moon_target": "Moon (301)",
    "center": "Earth (399)",
}

# Plausibility bands for geocentric distance
SUN_DISTANCE_MIN_KM = 1.45e8
SUN_DISTANCE_MAX_KM = 1.55e8
MOON_DISTANCE_MIN_KM = 350000.0
MOON_DISTANCE_MAX_KM = 412000.0


def _get(params: dict) -> bytes:
    url = BASE_URL + "?" + urllib.parse.urlencode(params)
    last_err = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=180) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code in (503, 429, 504):
                wait = REQUEST_SPACING_S * (attempt + 1) * 2
                print(f"  [retry {attempt+1}] HTTP {e.code}, sleeping {wait:.0f}s")
                time.sleep(wait)
            else:
                raise
    raise RuntimeError(f"Horizons fetch failed after 3 attempts: {last_err}")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _soe_rows(text: str) -> list[str]:
    lines = text.splitlines()
    soe = [i for i, ln in enumerate(lines) if ln.strip() == "$$SOE"]
    eoe = [i for i, ln in enumerate(lines) if ln.strip() == "$$EOE"]
    if len(soe) != 1 or len(eoe) != 1 or eoe[0] <= soe[0]:
        raise RuntimeError(f"SOE/EOE delimiters malformed: {soe} {eoe}")
    return lines[soe[0] + 1 : eoe[0]]


def _expected_n_rows(start_year: int, end_year: int) -> int:
    """Compute expected number of daily rows inclusive of both endpoints.

    Includes leap days: 2028, 2032, 2036, 2040 in our span.
    """
    n = end_year - start_year  # number of complete years
    n_days = 0
    for y in range(start_year, end_year):
        leap = (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)
        n_days += 366 if leap else 365
    return n_days + 1  # inclusive endpoint


def _validate_sun(text: str, start_year: int, end_year: int) -> dict:
    """Validate Sun vector response."""
    rows = _soe_rows(text)
    expected_n = _expected_n_rows(start_year, end_year)
    if abs(len(rows) - expected_n) > 3:
        raise RuntimeError(f"Sun rows {len(rows)} not close to expected {expected_n}")
    if "Target body name: Sun" not in text:
        raise RuntimeError("Sun target header missing")
    if "Center body name: Earth" not in text:
        raise RuntimeError("Earth center missing")
    if "Reference frame : ICRF" not in text:
        raise RuntimeError("ICRF frame missing")
    if "TDB" not in text:
        raise RuntimeError("TDB time missing")
    # Band-check at least first/last rows
    for k in (0, len(rows) - 1):
        parts = [p.strip() for p in rows[k].split(",")]
        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        d = (x * x + y * y + z * z) ** 0.5
        if not (SUN_DISTANCE_MIN_KM < d < SUN_DISTANCE_MAX_KM):
            raise RuntimeError(f"Sun distance {d} outside plausibility band at row {k}")
    return {"n_rows": len(rows), "expected_n_days": expected_n}


def _validate_moon(text: str, start_year: int, end_year: int) -> dict:
    """Validate Moon vector response."""
    rows = _soe_rows(text)
    expected_n = _expected_n_rows(start_year, end_year)
    if abs(len(rows) - expected_n) > 3:
        raise RuntimeError(f"Moon rows {len(rows)} not close to expected {expected_n}")
    if "Target body name: Moon" not in text:
        raise RuntimeError("Moon target header missing")
    if "Center body name: Earth" not in text:
        raise RuntimeError("Earth center missing")
    if "Reference frame : ICRF" not in text:
        raise RuntimeError("ICRF frame missing")
    if "TDB" not in text:
        raise RuntimeError("TDB time missing")
    # Band-check distance
    for k in (0, len(rows) - 1):
        parts = [p.strip() for p in rows[k].split(",")]
        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        d = (x * x + y * y + z * z) ** 0.5
        if not (MOON_DISTANCE_MIN_KM < d < MOON_DISTANCE_MAX_KM):
            raise RuntimeError(f"Moon distance {d} outside plausibility band at row {k}")
    return {"n_rows": len(rows), "expected_n_days": expected_n}


def build_chunks() -> list[tuple[int, int]]:
    """Return list of (start_year, end_year) chunks covering the full span."""
    chunks = []
    year = SPAN_START_YEAR
    while year < SPAN_END_YEAR:
        end = min(year + CHUNK_YEARS, SPAN_END_YEAR)
        chunks.append((year, end))
        year = end
    return chunks


def fetch_body(body: str, start_year: int, end_year: int) -> bytes:
    """Fetch Sun or Moon for one 5-year chunk."""
    if body == "sun":
        command = "'10'"
        target_token = IDENTITY_TOKENS["sun_target"]
    elif body == "moon":
        command = "'301'"
        target_token = IDENTITY_TOKENS["moon_target"]
    else:
        raise ValueError(f"unknown body: {body}")
    params = {
        "format": "text",
        "MAKE_EPHEM": "'YES'",
        "OBJ_DATA": "'NO'",
        "EPHEM_TYPE": "'VECTORS'",
        "COMMAND": command,
        "CENTER": "'500@399'",  # geocentric
        "REF_PLANE": "'FRAME'",
        "REF_SYSTEM": "'ICRF'",
        "VEC_TABLE": "'2'",
        "VEC_CORR": "'NONE'",
        "OUT_UNITS": "'KM-S'",
        "TIME_TYPE": "'TDB'",
        "CSV_FORMAT": "'YES'",
        "CAL_FORMAT": "'BOTH'",
        "START_TIME": f"'{start_year}-01-01 00:00'",
        "STOP_TIME": f"'{end_year}-01-01 00:00'",
        "STEP_SIZE": "'1 d'",
    }
    print(f"  fetching {target_token} {start_year}-{end_year} ...")
    t0 = time.time()
    data = _get(params)
    elapsed = time.time() - t0
    print(f"  received {len(data)} bytes in {elapsed:.1f}s; sha256={_sha256(data)[:16]}")
    return data


def main():
    chunks = build_chunks()
    print(f"Total chunks to acquire: {len(chunks)} bodies x {len(chunks)} years")
    print(f"Span: {SPAN_START_YEAR}-{SPAN_END_YEAR} ({(SPAN_END_YEAR - SPAN_START_YEAR)} years)")
    print(f"Per-body requests: {len(chunks)} chunks (5 yr each)")

    sun_chunks: dict[str, dict] = {}
    moon_chunks: dict[str, dict] = {}
    for body in ("sun", "moon"):
        for (s, e) in chunks:
            key = f"{s}_{e}"
            out_path = OUT_DIR / f"horizons_{body}_geocentric_vectors_{s}_to_{e}_icrf_tdb_daily.txt"
            # Idempotence: if chunk already exists and is valid, reuse it.
            if out_path.exists():
                text = out_path.read_bytes()
                text_str = text.decode("utf-8")
                try:
                    if body == "sun":
                        info = _validate_sun(text_str, s, e)
                    else:
                        info = _validate_moon(text_str, s, e)
                    sha = _sha256(text)
                    print(f"  [{body}] {s}-{e}: cache HIT; {info['n_rows']} rows; sha256={sha[:16]}")
                except RuntimeError as ve:
                    print(f"  [{body}] {s}-{e}: cache invalid ({ve}); re-fetching")
                    text = fetch_body(body, s, e)
                    text_str = text.decode("utf-8")
                    if body == "sun":
                        info = _validate_sun(text_str, s, e)
                    else:
                        info = _validate_moon(text_str, s, e)
                    out_path.write_bytes(text)
                    sha = _sha256(text)
                    print(f"  [{body}] {s}-{e}: re-fetched; {info['n_rows']} rows; sha256={sha[:16]}")
                    time.sleep(REQUEST_SPACING_S)
            else:
                text = fetch_body(body, s, e)
                text_str = text.decode("utf-8")
                if body == "sun":
                    info = _validate_sun(text_str, s, e)
                else:
                    info = _validate_moon(text_str, s, e)
                out_path.write_bytes(text)
                sha = _sha256(text)
                print(f"  [{body}] {s}-{e}: {info['n_rows']} rows; sha256={sha[:16]}")
                time.sleep(REQUEST_SPACING_S)
            if body == "sun":
                sun_chunks[key] = {
                    "path": str(out_path),
                    "sha256": sha,
                    "start_year": s,
                    "end_year": e,
                    "n_rows": info["n_rows"],
                }
            else:
                moon_chunks[key] = {
                    "path": str(out_path),
                    "sha256": sha,
                    "start_year": s,
                    "end_year": e,
                    "n_rows": info["n_rows"],
                }

    # Concatenate: deduplicate by Julian date (chunks share boundary day)
    def _concat_unique(chunks_dict, out_name):
        seen_jds = set()
        out_lines = []
        for k, v in chunks_dict.items():
            text = Path(v["path"]).read_text(encoding="utf-8")
            rows = _soe_rows(text)
            for row in rows:
                parts = [p.strip() for p in row.split(",")]
                jd = float(parts[0])
                if jd in seen_jds:
                    continue
                seen_jds.add(jd)
                out_lines.append(row)
        # Concatenated output: keep header from first chunk
        first_text = Path(chunks_dict[list(chunks_dict.keys())[0]]["path"]).read_text(encoding="utf-8")
        lines = first_text.splitlines()
        soe_idx = next(i for i, ln in enumerate(lines) if ln.strip() == "$$SOE")
        header = lines[: soe_idx + 1]
        eoe = ["$$EOE"]
        out_text = "\n".join(header + out_lines + eoe) + "\n"
        out_path = OUT_DIR / out_name
        out_path.write_text(out_text, encoding="utf-8")
        return out_path, _sha256(out_text.encode("utf-8")), len(out_lines)

    sun_concat_path, sun_concat_sha, sun_n = _concat_unique(
        sun_chunks, "horizons_sun_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
    )
    print(f"\nSun concatenated: {sun_n} rows; sha256={sun_concat_sha[:16]}")
    moon_concat_path, moon_concat_sha, moon_n = _concat_unique(
        moon_chunks, "horizons_moon_geocentric_vectors_2026_to_2045_icrf_tdb_daily.txt"
    )
    print(f"Moon concatenated: {moon_n} rows; sha256={moon_concat_sha[:16]}")

    manifest = {
        "description": "19-year DE441 geocentric Sun + Moon vectors for mission_lunisolar_closure.",
        "source": "NASA/JPL Horizons API (DE441)",
        "frame": "ICRF",
        "time_type": "TDB",
        "center": "Earth (399)",
        "cadence": "1 day",
        "units": "KM-S",
        "vectors_corrected": "NONE (geometric)",
        "sun_concat": {
            "path": str(sun_concat_path),
            "sha256": sun_concat_sha,
            "n_rows": sun_n,
        },
        "moon_concat": {
            "path": str(moon_concat_path),
            "sha256": moon_concat_sha,
            "n_rows": moon_n,
        },
        "per_chunk": {
            "sun": sun_chunks,
            "moon": moon_chunks,
        },
        "acquired_utc": datetime.now(timezone.utc).isoformat(),
        "acquisition_doctrine": "Exp 014/017 byte-pinned snapshot pattern",
    }
    manifest_path = OUT_DIR / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(f"\nMANIFEST written to {manifest_path}")
    print("\nNext: copy the concatenated files to the mission's reference/ directory.")
    return 0


if __name__ == "__main__":
    sys.exit(main())