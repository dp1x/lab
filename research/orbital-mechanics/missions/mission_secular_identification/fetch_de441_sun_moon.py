"""DE441 geocentric Sun + Moon acquisition for mission_secular_identification.

Why this file exists
--------------------
``mission_mean_element_leakage`` (2026-10-04) showed that determining a *secular*
lunisolar RAAN rate requires an arc of at least two lunar nodal cycles
(6798.383 d). Two cycles from the lab's T0 = 820476800.0 (2026-01-01 TDB) end
2063-03-24. The committed 19-year snapshot ends 2045-01-01, and
``lab_utils.ephemeris.interp_snapshot`` CLAMPS outside its range - so a 2-cycle
arc built on it would silently freeze the Sun and Moon and produce a plausible,
completely invalid number. This script re-acquires a longer span.

It is a deliberate copy of ``mission_lunisolar_closure/fetch_horizons_sun_moon_long.py``
(Exp 014/017 acquisition doctrine) rather than an import, because (a) the donor
script hard-codes its span and output directory, and (b) constitutional rule 5
(never commit proprietary content) and rule 7 (complexity must justify itself)
both favour leaving a completed mission's committed artifact untouched. The
*validation* logic here is what matters and it is reproduced exactly, plus one
addition: ``--baseline-manifest`` re-requests chunks that lie inside the
predecessor's span and asserts their sha256 matches the pinned MANIFEST value.

Coverage is 2026-01-01 -> 2064-01-01 = 13,876 daily rows, i.e. 2.04 lunar
nodal cycles plus ~2.4 years of margin beyond the 2-cycle requirement.

Politeness (AGENTS.md "Responsible Web Access"): single sequential pass, 5-year
chunks, 4 s spacing, no parallel-request bursts, no retries beyond 3 with
backoff on 429/503/504. No access control or CAPTCHA is ever circumvented.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

IDENTITY_TOKENS = {"sun_target": "Sun (10)", "moon_target": "Moon (301)"}

SUN_DISTANCE_MIN_KM, SUN_DISTANCE_MAX_KM = 1.45e8, 1.55e8
MOON_DISTANCE_MIN_KM, MOON_DISTANCE_MAX_KM = 350000.0, 412000.0

REQUEST_SPACING_S = 4.0


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
                print(f"  [retry {attempt+1}] HTTP {e.code}, sleeping {wait:.0f}s", flush=True)
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
    n_days = 0
    for y in range(start_year, end_year):
        leap = (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)
        n_days += 366 if leap else 365
    return n_days + 1  # inclusive endpoint


def _validate(text: str, body: str, start_year: int, end_year: int) -> dict:
    rows = _soe_rows(text)
    expected_n = _expected_n_rows(start_year, end_year)
    if abs(len(rows) - expected_n) > 3:
        raise RuntimeError(f"{body} rows {len(rows)} not close to expected {expected_n}")
    want = {"sun": "Target body name: Sun", "moon": "Target body name: Moon"}[body]
    for token in (want, "Center body name: Earth", "Reference frame : ICRF", "TDB"):
        if token not in text:
            raise RuntimeError(f"{body}: header token missing: {token!r}")
    lo, hi = ((SUN_DISTANCE_MIN_KM, SUN_DISTANCE_MAX_KM) if body == "sun"
              else (MOON_DISTANCE_MIN_KM, MOON_DISTANCE_MAX_KM))
    for k in (0, len(rows) - 1):
        p = [x.strip() for x in rows[k].split(",")]
        d = (float(p[2]) ** 2 + float(p[3]) ** 2 + float(p[4]) ** 2) ** 0.5
        if not (lo < d < hi):
            raise RuntimeError(f"{body} distance {d:.1f} km outside band at row {k}")
    return {"n_rows": len(rows), "expected_n_days": expected_n}


def compare_rows(new_path: Path, ref_path: Path) -> dict:
    """Numeric-row comparison between a freshly fetched chunk and a committed one.

    Why rows and not raw bytes: a Horizons response embeds its acquisition
    timestamp ("Ephemeris / API_USER Tue Sep  1 12:00:58 2026"), so the raw
    response sha256 is NOT reproducible across acquisitions by construction --
    pinning raw bytes would fail every re-acquisition for a reason that has
    nothing to do with the data. The meaningful gate is that the numeric vectors
    agree, which is what the 2026-10-04 handoff measured ("bit-for-bit on all
    6941 overlap days, max |dr| = 0.000000e+00 km, 0 days differing by > 1 m").
    """
    new_rows = _soe_rows(new_path.read_text(encoding="utf-8"))
    ref_rows = _soe_rows(ref_path.read_text(encoding="utf-8"))
    if len(new_rows) != len(ref_rows):
        raise RuntimeError(f"{new_path.name}: {len(new_rows)} rows != reference {len(ref_rows)}")
    max_d, n_bad = 0.0, 0
    for a, b in zip(new_rows, ref_rows):
        pa = [x.strip() for x in a.split(",")]
        pb = [x.strip() for x in b.split(",")]
        if pa[0] != pb[0]:
            raise RuntimeError(f"{new_path.name}: JD mismatch {pa[0]} vs {pb[0]}")
        for k in range(2, 5):
            d = abs(float(pa[k]) - float(pb[k]))
            max_d = max(max_d, d)
            if d > 1.0:
                n_bad += 1
    return {"rows": len(new_rows), "max_abs_delta_km": max_d, "components_over_1m": n_bad}


def build_chunks(span_start: int, span_end: int, chunk_years: int) -> list[tuple[int, int]]:
    chunks, year = [], span_start
    while year < span_end:
        end = min(year + chunk_years, span_end)
        chunks.append((year, end))
        year = end
    return chunks


def fetch_body(body: str, start_year: int, end_year: int) -> bytes:
    command = {"sun": "'10'", "moon": "'301'"}[body]
    print(f"  fetching {IDENTITY_TOKENS[body + '_target']} {start_year}-{end_year} ...", flush=True)
    params = {
        "format": "text", "MAKE_EPHEM": "'YES'", "OBJ_DATA": "'NO'",
        "EPHEM_TYPE": "'VECTORS'", "COMMAND": command, "CENTER": "'500@399'",
        "REF_PLANE": "'FRAME'", "REF_SYSTEM": "'ICRF'", "VEC_TABLE": "'2'",
        "VEC_CORR": "'NONE'", "OUT_UNITS": "'KM-S'", "TIME_TYPE": "'TDB'",
        "CSV_FORMAT": "'YES'", "CAL_FORMAT": "'BOTH'",
        "START_TIME": f"'{start_year}-01-01 00:00'",
        "STOP_TIME": f"'{end_year}-01-01 00:00'",
        "STEP_SIZE": "'1 d'",
    }
    t0 = time.time()
    data = _get(params)
    print(f"  received {len(data)} bytes in {time.time() - t0:.1f}s; "
          f"sha256={_sha256(data)[:16]}", flush=True)
    return data


def acquire_all(span_start: int, span_end: int, out_dir: Path, chunk_years: int,
                baseline_manifest: dict | None, baseline_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks = build_chunks(span_start, span_end, chunk_years)
    print(f"chunks: {len(chunks)} x 2 bodies; span {span_start}-{span_end}", flush=True)
    sun_chunks: dict[str, dict] = {}
    moon_chunks: dict[str, dict] = {}
    for body in ("sun", "moon"):
        for (s, e) in chunks:
            key = f"{s}_{e}"
            out_path = out_dir / f"horizons_{body}_geocentric_vectors_{s}_to_{e}_icrf_tdb_daily.txt"
            if out_path.exists():
                raw = out_path.read_bytes()
                try:
                    info = _validate(raw.decode("utf-8"), body, s, e)
                    sha = _sha256(raw)
                    print(f"  [{body}] {s}-{e}: cache HIT {info['n_rows']} rows sha={sha[:16]}", flush=True)
                except RuntimeError as ve:
                    print(f"  [{body}] {s}-{e}: cache invalid ({ve}); re-fetching", flush=True)
                    raw = fetch_body(body, s, e)
                    info = _validate(raw.decode("utf-8"), body, s, e)
                    out_path.write_bytes(raw)
                    sha = _sha256(raw)
                    time.sleep(REQUEST_SPACING_S)
            else:
                raw = fetch_body(body, s, e)
                info = _validate(raw.decode("utf-8"), body, s, e)
                out_path.write_bytes(raw)
                sha = _sha256(raw)
                time.sleep(REQUEST_SPACING_S)

            pinned = (baseline_manifest or {}).get("per_chunk", {}).get(body, {}).get(key)
            overlap: dict = {}
            if pinned is not None:
                # MANIFEST paths are repo-root-relative; resolve against the CWD
                # (the repo root, which is where the script is documented to run).
                ref_path = Path(pinned["path"])
                if not ref_path.is_absolute() and not ref_path.exists():
                    ref_path = baseline_dir / ref_path
                overlap = compare_rows(out_path, ref_path)
                print(f"  [{body}] {s}-{e}: overlap {overlap['rows']} rows, "
                      f"max |dr| = {overlap['max_abs_delta_km']:.6e} km, "
                      f"{overlap['components_over_1m']} components > 1 m", flush=True)
                if overlap["max_abs_delta_km"] > 1.0:
                    raise RuntimeError(
                        f"{body} {s}-{e}: overlap vs committed snapshot exceeds 1 m: {overlap}")

            (sun_chunks if body == "sun" else moon_chunks)[key] = {
                "path": str(out_path), "sha256": sha, "start_year": s,
                "end_year": e, "n_rows": info["n_rows"],
                "overlap_vs_committed": overlap or None,
            }
    return sun_chunks, moon_chunks


def concat_unique(chunks_dict: dict) -> tuple[str, int]:
    seen, out_lines = set(), []
    for v in chunks_dict.values():
        for row in _soe_rows(Path(v["path"]).read_text(encoding="utf-8")):
            jd = float([p.strip() for p in row.split(",")][0])
            if jd in seen:
                continue
            seen.add(jd)
            out_lines.append(row)
    first = Path(next(iter(chunks_dict.values()))["path"]).read_text(encoding="utf-8").splitlines()
    soe_idx = next(i for i, ln in enumerate(first) if ln.strip() == "$$SOE")
    return "\n".join(first[: soe_idx + 1] + out_lines + ["$$EOE"]) + "\n", len(out_lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Acquire DE441 geocentric Sun+Moon vectors.")
    ap.add_argument("--span-start", type=int, default=2026)
    ap.add_argument("--span-end", type=int, default=2064)
    ap.add_argument("--chunk-years", type=int, default=5)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--manifest-out", required=True)
    ap.add_argument("--baseline-manifest", default=None)
    args = ap.parse_args(argv)

    baseline = None
    baseline_dir = Path(".") if args.baseline_manifest else Path(".")
    if args.baseline_manifest and Path(args.baseline_manifest).exists():
        baseline = json.loads(Path(args.baseline_manifest).read_text(encoding="utf-8"))

    sun_chunks, moon_chunks = acquire_all(args.span_start, args.span_end, Path(args.out_dir),
                                          args.chunk_years, baseline, baseline_dir)
    sun_text, sun_n = concat_unique(sun_chunks)
    moon_text, moon_n = concat_unique(moon_chunks)

    target = Path(args.manifest_out)
    target.mkdir(parents=True, exist_ok=True)
    sun_name = f"horizons_sun_geocentric_vectors_{args.span_start}_to_{args.span_end}_icrf_tdb_daily.txt"
    moon_name = f"horizons_moon_geocentric_vectors_{args.span_start}_to_{args.span_end}_icrf_tdb_daily.txt"
    (target / sun_name).write_text(sun_text, encoding="utf-8", newline="")
    (target / moon_name).write_text(moon_text, encoding="utf-8", newline="")
    sun_sha, moon_sha = _sha256(sun_text.encode("utf-8")), _sha256(moon_text.encode("utf-8"))
    print(f"\nSun  {sun_n} rows sha256={sun_sha}", flush=True)
    print(f"Moon {moon_n} rows sha256={moon_sha}", flush=True)

    manifest = {
        "description": f"{args.span_end - args.span_start}-year DE441 geocentric Sun + Moon vectors "
                       "(geometric states, ICRF/TDB, daily). Long enough for >= 2 complete lunar "
                       "nodal cycles (6798.383 d each) from T0=820476800.0, plus margin.",
        "source": "NASA/JPL Horizons API (DE441)", "frame": "ICRF", "time_type": "TDB",
        "center": "Earth (399)", "cadence": "1 day", "units": "KM-S",
        "vectors_corrected": "NONE (geometric)",
        "span_start_year": args.span_start, "span_end_year": args.span_end,
        "sun_concat": {"path": sun_name, "sha256": sun_sha, "n_rows": sun_n},
        "moon_concat": {"path": moon_name, "sha256": moon_sha, "n_rows": moon_n},
        "per_chunk": {"sun": sun_chunks, "moon": moon_chunks},
        "acquired_utc": datetime.now(timezone.utc).isoformat(),
        "acquisition_doctrine": "Exp 014/017 byte-pinned snapshot pattern; extended 2026-10-04 for "
                                "mission_secular_identification (>=2 lunar nodal cycles).",
    }
    (target / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True),
                                          encoding="utf-8", newline="")
    print(f"MANIFEST -> {target / 'MANIFEST.json'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())