# License / Attribution Notes

## Source

- **Service:** NASA/JPL Solar System Dynamics Horizons system
  (`https://ssd.jpl.nasa.gov/api/horizons.api`), official public REST API.
- **Files:** `horizons_-125544_iss_vectors_2026-08-24_to_2026-08-27_tdb_2min.txt`
  (geometric ICRF state vectors of the International Space Station, geocentric,
  TDB epochs, km-s units) and `horizons_-125544_iss_objdata.txt` (object-data
  provenance sheet).
- **Retrieved:** 2026-08-24 (UTC timestamps in `MANIFEST.json`), via one-shot
  serial acquisition (`fetch_horizons_snapshot.py`, human-use-rate access).

## Terms

JPL Horizons is a public government service; responses carry no redistribution
restriction notice. The response footer attributes computations to:

> Computations by ... Solar System Dynamics Group, Horizons On-Line Ephemeris
> System, 4800 Oak Grove Drive, Jet Propulsion Laboratory, Pasadena, CA 91109 USA

These snapshots are retained verbatim (byte-for-byte, SHA-256-pinned in
`MANIFEST.json`) as the immutable external reference ("answer key") for this
experiment, per laboratory policy that authoritative external inputs must be
machine-pinned to keep offline reruns deterministic.

## Reference-quality disclosure (verbatim from the object sheet)

> Trajectory is TLE-based. Predicts run for 4 weeks into future, but are of
> low accuracy for times more than a few days past the revision date above.

Revision date at acquisition: **Aug 23, 2026**. The comparison window
(2026-Aug-24 .. 2026-Aug-27 TDB) starts inside the disclosed validity period;
residuals against this reference therefore measure lab-model-vs-reference
divergence jointly with the reference's own uncertainty (never separated).
