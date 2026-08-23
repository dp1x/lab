"""Exp 013 validation tests: snapshot integrity, parsing, time/frame gates,
model nulls, convergence, interpolation, determinism, adversarial mutants.

Layer map (house style):
  L0 committed-artifact integrity + governance (.gitattributes protection)
  L1 header identity assertions
  L2 parser exactness, dual-pipeline equality, malformed-table mutants
  L3 epoch alignment (independent ordinals, calendar-vs-JD, spacing)
  L4 plausibility/continuity/frame gates on the pinned states
  L5 osculating-element sanity incl. Exp 009 nodal-rate anchor
  L6 model nulls vs donors (J2=0 bit-exact; drag=0 bit-exact; beta<0 rejected)
  L7 RIC frame construction properties
  L8 dt-convergence order
  L9 Hermite interpolation bound (and linear-interpolant mutant detection)
  L10 predictor isolation (anti-peeking: truncation + future-row mutation)
  L11 determinism of the numeric core
  L12 adversarial mutants (frame rotation, UTC shift, unit/column corruption,
      wrong object, tampered snapshot)
  L13 stale-results guard (code hashes in results.json match working tree)
  L14 offline doctrine (no network imports in the analysis path)
"""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pytest

_EXPERIMENTS_DIR = Path(__file__).resolve().parents[2]
_EXP_DIR = _EXPERIMENTS_DIR / "jplValidation"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


exp = _load("jpl013_experiment", _EXP_DIR / "experiment.py")

_006_PATH = _EXPERIMENTS_DIR / "planeChangeManeuvers" / "experiment.py"
pcm = _load("pcm006_for_jpl013_tests", _006_PATH)


@pytest.fixture(scope="session")
def ref():
    return exp.load_reference()


# ------------------------------------------------------------------ L0 ---- #
def test_l0_snapshot_files_match_manifest():
    manifest = exp.load_manifest()
    raw_v, raw_o = exp.verify_snapshot_bytes(manifest)
    assert len(raw_v) == 426243 and len(raw_o) > 0


def test_l0_manifest_schema_and_identity():
    m = exp.load_manifest()
    qv = m["acquisition"]["query_params"]["vectors"]
    assert qv["COMMAND"] == "'-125544'"
    assert qv["CENTER"] == "'500@399'"
    assert qv["REF_PLANE"] == "'FRAME'" and qv["REF_SYSTEM"] == "'ICRF'"
    assert qv["VEC_TABLE"] == "'2'" and qv["VEC_CORR"] == "'NONE'"
    assert qv["OUT_UNITS"] == "'KM-S'" and qv["TIME_TYPE"] == "'TDB'"
    assert qv["CSV_FORMAT"] == "'YES'"
    for rec in m["snapshot"]["files"].values():
        assert re.fullmatch(r"[0-9a-f]{64}", rec["sha256"])
    assert m["provenance"]["trajectory_tle_based_disclosed"] is True


def test_l0_gitattributes_protects_reference_bytes():
    ga = Path(__file__).resolve().parents[5] / ".gitattributes"
    text = ga.read_text(encoding="utf-8")
    assert re.search(r"jplValidation/reference/\*\.txt\s+-text", text), (
        "reference snapshots must be marked -text before commit "
        "(core.autocrlf would corrupt their pinned hashes on fresh clones)"
    )


# ------------------------------------------------------------------ L1 ---- #
def test_l1_header_identity_assertions_pass_on_snapshot(ref):
    hdr = ref["header"]
    assert "-125544" in hdr["target_body_name"]
    assert "Earth (399)" in hdr["center_body_name"]
    assert hdr["reference_frame"] == "ICRF"
    assert hdr["output_units"] == "KM-S"
    assert "GEOMETRIC" in hdr["output_type"]
    assert ref["revision_date"] != "<missing>"


def test_l1_wrong_object_header_rejected():
    fake = "Target body name: Hubble Space Telescope (spacecraft) (-48)\nReference frame : ICRF\nOutput units : KM-S\nGEOMETRIC cartesian states\nTDB\nCenter body name: Earth (399)\nBODY CENTER"
    with pytest.raises(RuntimeError):
        exp.parse_header(fake)


def test_l1_missing_frame_declaration_rejected():
    fake = (
        "Target body name: International Space Station (spacecraft) (-125544)\n"
        "Center body name: Earth (399)\nBODY CENTER\nOutput units : KM-S\n"
        "GEOMETRIC cartesian states\nTDB epochs here\n"
    )  # no 'Reference frame : ICRF' line -> default-plane hazard must fail loud
    with pytest.raises(RuntimeError):
        exp.parse_header(fake)


# ------------------------------------------------------------------ L2 ---- #
def test_l2_dual_parse_pipelines_agree(ref):
    text = (exp.REFERENCE_DIR / exp.VECTORS_NAME).read_text(encoding="utf-8")
    lines = text.splitlines()
    soe = lines.index("$$SOE")
    eoe = lines.index("$$EOE")
    rows = lines[soe + 1 : eoe]
    jd_a, st_a, _ = exp._parse_rows_split(rows)
    jd_b, st_b = exp._parse_rows_regex(rows)
    assert len(jd_a) == exp.EXPECTED_ROWS == len(st_b)
    assert np.array_equal(jd_a, jd_b)
    assert np.array_equal(st_a, st_b)


def test_l2_regex_pipeline_rejects_garbage():
    with pytest.raises(RuntimeError):
        exp._parse_rows_regex(["not a data row at all"])


def test_l2_extra_column_table_fails_loudly():
    bad = "2461276.500000000, A.D. 2026-Aug-24 00:00:00.0000, 1.0E+03, 2.0E+03, 3.0E+03, 1.0E+00, 2.0E+00, 3.0E+00, 9.9E+09,"
    with pytest.raises(RuntimeError):
        exp._parse_rows_split([bad])


def test_l2_truncated_row_fails_loudly():
    bad = "2461276.500000000, A.D. 2026-Aug-24 00:00:00.0000, 1.0E+03, 2.0E+03,"
    with pytest.raises(RuntimeError):
        exp._parse_rows_split([bad])


def test_l2_malformed_float_fails_loudly():
    bad = "2461276.500000000, A.D. 2026-Aug-24 00:00:00.0000, 1,0E+03, 2.0E+03, 3.0E+03, 1.0E+00, 2.0E+00, 3.0E+00,"
    with pytest.raises((RuntimeError, ValueError)):
        exp._parse_rows_split([bad])


# ------------------------------------------------------------------ L3 ---- #
def test_l3_independent_ordinal_pipelines_match_known_dates():
    # Anchors derived by hand: Unix epoch 719163; 2000-01-01 = 730120
    # (JD 2451545 = 2000-01-01T12:00); century spans 36524/36525 days;
    # day-of-year arithmetic for Mar 1 / Aug 24 shown in card provenance.
    cases = [(1970, 1, 1, 719163), (2000, 1, 1, 730120), (2000, 3, 1, 730180),
             (2026, 8, 24, 739852), (1900, 3, 1, 693655), (2100, 3, 1, 766704)]
    for y, m, d, want in cases:
        assert exp._ordinal_formula(y, m, d) == want
        assert exp._ordinal_table(y, m, d) == want


def test_l3_calendar_matches_jd_clock_within_tolerance(ref):
    """Calendar column must be exactly consistent with its row's JDN
    (day level via independent ordinals, time-of-day to < 1 ms)."""
    text = (exp.REFERENCE_DIR / exp.VECTORS_NAME).read_text(encoding="utf-8")
    rows = text.splitlines()[text.splitlines().index("$$SOE") + 1 : text.splitlines().index("$$EOE")]
    parts = [r.split(",") for r in rows]
    for k in (0, 1, 500, 1080, 1440, 2160):
        val = exp.calendar_epoch_s(parts[k][1].strip(), float(parts[k][0]), float(parts[0][0]))
        t_s = (float(parts[k][0]) - float(parts[0][0])) * 86400.0
        assert abs(val - t_s) < 61.0


def test_l3_corrupted_calendar_field_detected():
    # Calendar day shifted by +1 day must fail the JDN<->ordinal consistency.
    with pytest.raises(ValueError):
        exp.calendar_epoch_s("A.D. 2026-Aug-25 00:06:00.0000", 2461276.504166667, 2461276.5)
    # Time-of-day corrupted by >1 ms but same date must also fail.
    with pytest.raises(ValueError):
        exp.calendar_epoch_s("A.D. 2026-Aug-24 00:07:30.0000", 2461276.504166667, 2461276.5)
    # Consistent pair passes and returns the JD-relative seconds.
    ok = exp.calendar_epoch_s("A.D. 2026-Aug-24 00:06:00.0000", 2461276.504166667, 2461276.5)
    assert abs(ok - 360.0) < 1e-3


def test_l3_spacing_monotone_uniform(ref):
    dts = np.diff(ref["t_s"])
    assert np.all(dts > 0)
    assert np.max(np.abs(dts - exp.CADENCE_S)) <= 2e-4
    assert len(ref["t_s"]) == exp.EXPECTED_ROWS


# ------------------------------------------------------------------ L4 ---- #
def test_l4_plausibility_gates_accept_real_states(ref):
    exp._plausibility_gates(ref["r"], ref["v"])  # must not raise


def test_l4_km_to_m_mutant_caught():
    ref = exp.load_reference()
    with pytest.raises(RuntimeError):
        exp._plausibility_gates(ref["r"] * 1000.0, ref["v"] * 1000.0)


def test_l4_state_order_swap_caught():
    ref = exp.load_reference()
    with pytest.raises(RuntimeError):
        exp._plausibility_gates(ref["v"], ref["r"])


def test_l4_velocity_sign_flip_caught():
    ref = exp.load_reference()
    with pytest.raises(RuntimeError):
        exp._plausibility_gates(ref["r"], -ref["v"])


def test_l4_ecliptic_rotation_mutant_caught():
    ref = exp.load_reference()
    th = np.radians(23.4393)  # ecliptic obliquity: silent REF_PLANE=ECLIPTIC trap
    Rx = np.array([[1, 0, 0], [0, np.cos(th), -np.sin(th)], [0, np.sin(th), np.cos(th)]])
    with pytest.raises(RuntimeError):
        exp._plausibility_gates(ref["r"] @ Rx.T, ref["v"] @ Rx.T)


# ------------------------------------------------------------------ L5 ---- #
def test_l5_inclination_is_iss_like(ref):
    el = exp.rv_to_coe_eci(ref["r"], ref["v"])
    inc_deg = np.degrees(el["inc"])
    assert abs(np.mean(inc_deg) - 51.63) < 0.75
    assert np.max(np.abs(inc_deg - np.mean(inc_deg))) < 0.2


def test_l5_nodal_rate_matches_j2_first_order_anchor(ref):
    """Independent physics anchor (Exp 009): reference RAAN regression must
    reproduce the first-order J2 nodal rate for these elements. A wrong J2
    sign/convention or a frame error moves this grossly while along-track
    RMS can look deceptively normal."""
    el = exp.rv_to_coe_eci(ref["r"], ref["v"])
    om = np.unwrap(el["Omega"])
    slope = np.polyfit(ref["t_s"] / 86400.0, om, 1)[0]
    rate_deg_day = np.degrees(slope)
    p = float(el["a"][0]) * (1.0 - float(el["e"][0]) ** 2)
    n = exp.mean_motion(float(el["a"][0]))
    analytic_rad_s = -1.5 * n * exp.J2_EARTH * (exp.R_EARTH_KM / p) ** 2 * np.cos(float(el["inc"][0]))
    analytic_deg_day = np.degrees(analytic_rad_s * 86400.0)
    assert abs(rate_deg_day - analytic_deg_day) < 0.05 * abs(analytic_deg_day)
    assert -5.2 < rate_deg_day < -4.7  # ISS-class anchor band (Exp 009: ~ -4.95)


# ------------------------------------------------------------------ L6 ---- #
def test_l6_j2_off_path_bit_exact_vs_exp006_donor(ref):
    x0 = np.concatenate([ref["r"][0], ref["v"][0]])
    t = ref["t_s"][:61]  # 2 h arc
    ours = exp.rk4_propagate(exp.j2_rhs(exp.MU_EARTH_KM3S2, 0.0), t, x0)
    theirs = pcm.propagate_3d_rk4(ref["r"][0], ref["v"][0], exp.MU_EARTH_KM3S2, t, 120.0)
    assert np.array_equal(ours[:, :3], theirs[:, :3])
    assert np.array_equal(ours[:, 3:], theirs[:, 3:])


def test_l6_drag_zero_bit_exact_vs_j2_canon(ref):
    x0 = np.concatenate([ref["r"][0], ref["v"][0]])
    t_dense = exp.dense_grid(ref["t_s"][:121], 4)
    drag_out = exp.propagate_3d_rk4_drag(
        ref["r"][0], ref["v"][0], exp.MU_EARTH_KM3S2, t_dense,
        j2=exp.J2_EARTH, beta=0.0,
    )
    canon_out = exp.rk4_propagate(
        exp.j2_rhs(exp.MU_EARTH_KM3S2, exp.J2_EARTH), t_dense, x0
    )
    assert np.array_equal(drag_out, canon_out)


def test_l6_negative_beta_rejected_by_donor(ref):
    with pytest.raises(ValueError):
        exp.propagate_3d_rk4_drag(
            ref["r"][0], ref["v"][0], exp.MU_EARTH_KM3S2, ref["t_s"][:10],
            j2=exp.J2_EARTH, beta=-100.0,
        )


# ------------------------------------------------------------------ L7 ---- #
def test_l7_ric_frames_orthonormal_right_handed(ref):
    r_hat, t_hat, c_hat = exp.ric_frames(ref["r"], ref["v"])
    for vec in (r_hat, t_hat, c_hat):
        assert np.allclose(np.einsum("ij,ij->i", vec, vec), 1.0, atol=1e-12)
    cross = np.einsum("ij,ij->i", np.cross(r_hat, t_hat), c_hat)
    assert np.allclose(cross, 1.0, atol=1e-12)


def test_l7_pure_alongtrack_offset_projects_onto_transverse(ref):
    r_hat, t_hat, c_hat = exp.ric_frames(ref["r"], ref["v"])
    dr = t_hat * 3.0  # purely in-track synthetic offset
    ric = exp.project_ric(dr, (r_hat, t_hat, c_hat))
    assert np.allclose(ric[:, 1], 3.0, atol=1e-9)
    assert np.allclose(ric[:, 0], 0.0, atol=1e-9)
    assert np.allclose(ric[:, 2], 0.0, atol=1e-9)


def test_l7_frame_origin_diagnostic_small(ref):
    """Predicted-built triad must nearly reproduce reference-built projections
    (axis-mixing bounded by angular separation ~ |dr|/r)."""
    states = exp.propagate_models(ref, exp.HEADLINE_NSUB)
    frames_ref = exp.ric_frames(ref["r"], ref["v"])
    dr = states["M2"][:, :3] - ref["r"]
    a = exp.project_ric(dr, frames_ref)[:, 1]
    b = exp.project_ric(dr, exp.ric_frames(states["M2"][:, :3], states["M2"][:, 3:]))[:, 1]
    scale = max(float(np.max(np.abs(a))), 1e-9)
    assert float(np.max(np.abs(a - b))) / scale < 0.02


# ------------------------------------------------------------------ L8 ---- #
@pytest.mark.parametrize("nsub_pair", [(1, 4), (2, 8)])
def test_l8_integration_selfconvergence_order(ref, nsub_pair):
    """Pure integration-order check: differences between same-physics runs at
    different dt contain NO model gap, so their decay must show RK4's order.
    (Residuals vs Horizons plateau at the model-mismatch floor and carry no
    order information - Exp 009 doctrine.)"""
    k = 181  # first 6 h of the window
    sub = {"t_s": ref["t_s"][:k], "r": ref["r"][:k], "v": ref["v"][:k]}
    truth = exp.propagate_models(sub, 16)["M2"]

    def err_vs_truth(nsub):
        st = exp.propagate_models(sub, nsub)["M2"]
        return exp.rms(np.linalg.norm(st[:, :3] - truth[:, :3], axis=1))

    e_a, e_b = err_vs_truth(nsub_pair[0]), err_vs_truth(nsub_pair[1])
    order = np.log(e_a / e_b) / np.log(2.0)
    assert order > 3.0, f"RK4 self-convergence order degraded: {order}"


# ------------------------------------------------------------------ L9 ---- #
def test_l9_hermite_within_analytic_bound(ref):
    study = exp.interpolation_study(ref)
    for key, rec in study.items():
        assert rec["within_bound"], f"{key}: {rec}"


def test_l9_linear_interpolant_flagged_as_bound_violation(ref):
    stride = 10
    idx = np.arange(0, ref["rows"], stride)
    held = np.array([i for i in range(1, ref["rows"] - 1) if i not in set(idx.tolist())])
    sub = np.concatenate([ref["r"][idx], ref["v"][idx]], axis=1)
    i = np.clip(np.searchsorted(ref["t_s"][idx], ref["t_s"][held]) - 1, 0, len(idx) - 2)
    u = ((ref["t_s"][held] - ref["t_s"][idx][i]) / (ref["t_s"][idx][i + 1] - ref["t_s"][idx][i]))[:, None]
    lin = (1 - u) * sub[i, :3] + u * sub[i + 1, :3]
    err_lin = float(np.max(np.linalg.norm(lin - ref["r"][held], axis=1)))
    rn = np.linalg.norm(ref["r"], axis=1)
    n_mean = exp.mean_motion(float(np.mean(rn)))
    bound_hermite = float(np.mean(rn) * (n_mean * exp.CADENCE_S * stride) ** 4 / 384.0)
    assert err_lin > 2.0 * bound_hermite  # wrong interpolant is detectable


# ----------------------------------------------------------------- L10 ---- #
def test_l10_predictor_ignores_future_reference_rows(ref):
    """Mutating any FUTURE reference row must leave predictions bit-identical:
    the propagator consumes only x(t0). Catches reference-interpolation peeking."""
    mutated = {**ref, "r": ref["r"].copy(), "v": ref["v"].copy()}
    mutated["r"][900] += 1.0
    k = 301
    sub_full = {"t_s": ref["t_s"][:k], "r": ref["r"][:k], "v": ref["v"][:k]}
    sub_muta = {"t_s": mutated["t_s"][:k], "r": mutated["r"][:k], "v": mutated["v"][:k]}
    a = exp.propagate_models(sub_full, 2)
    b = exp.propagate_models(sub_muta, 2)
    for key in ("M1", "M2"):
        assert np.array_equal(a[key], b[key])


def test_l10_truncation_invariance_no_hidden_consumption(ref):
    full = exp.propagate_models(ref, 2)
    k = 401
    part = exp.propagate_models({"t_s": ref["t_s"][:k], "r": ref["r"][:k],
                                 "v": ref["v"][:k]}, 2)
    for key in ("M1", "M2"):
        assert np.array_equal(full[key][:k], part[key])


# ----------------------------------------------------------------- L11 ---- #
def test_l11_numeric_core_deterministic(ref):
    a = exp.propagate_models(ref, 2, betas=(100.0,))
    b = exp.propagate_models(ref, 2, betas=(100.0,))
    for key in a:
        assert np.array_equal(a[key], b[key])
    m1 = exp.residual_metrics("M2", a["M2"], ref)
    m2 = exp.residual_metrics("M2", b["M2"], ref)
    assert json.dumps(m1, sort_keys=True) == json.dumps(m2, sort_keys=True)


def test_l11_bootstrap_seeded_reproducible(ref):
    states = exp.propagate_models(ref, 2)
    frames = exp.ric_frames(ref["r"], ref["v"])
    tra_m2 = exp.project_ric(states["M2"][:, :3] - ref["r"], frames)[:, 1]
    tra_m1 = exp.project_ric(states["M1"][:, :3] - ref["r"], frames)[:, 1]
    s1 = exp.block_bootstrap_skill(tra_m1, tra_m2)
    s2 = exp.block_bootstrap_skill(tra_m1, tra_m2)
    assert s1 == s2


# ----------------------------------------------------------------- L12 ---- #
def test_l12_utc_shift_time_trap_magnitude_documented(ref):
    """The TT-UTC (69.184 s) labeling catastrophe must be worth ~530 km along-
    track at ISS speed -- quantified so no one can mistake it for noise."""
    v_mean = float(np.mean(np.linalg.norm(ref["v"], axis=1)))
    shift_km = v_mean * 69.184
    assert 450.0 < shift_km < 600.0


def test_l12_seconds_to_days_unit_mutant_caught():
    ref = exp.load_reference()
    with pytest.raises(RuntimeError):
        exp._plausibility_gates(ref["r"], ref["v"] / 86400.0)


def test_l12_axis_swap_yz_caught():
    ref = exp.load_reference()
    r = ref["r"].copy()
    v = ref["v"].copy()
    r[:, [1, 2]] = r[:, [2, 1]]
    v[:, [1, 2]] = v[:, [2, 1]]
    with pytest.raises(RuntimeError):
        exp._plausibility_gates(r, v)


def test_l12_tampered_snapshot_fails_hash(tmp_path, monkeypatch):
    import hashlib
    raw = (exp.REFERENCE_DIR / exp.VECTORS_NAME).read_bytes()
    tampered = bytearray(raw)
    tampered[bytes(tampered).index(b"$$SOE") + 10] ^= 0x01
    (tmp_path / exp.VECTORS_NAME).write_bytes(bytes(tampered))
    (tmp_path / exp.OBJDATA_NAME).write_bytes(
        (exp.REFERENCE_DIR / exp.OBJDATA_NAME).read_bytes()
    )
    manifest = {
        "snapshot": {"files": {
            exp.VECTORS_NAME: {"bytes": len(tampered),
                               "sha256": hashlib.sha256(raw).hexdigest()},
            exp.OBJDATA_NAME: {"bytes": 667, "sha256": "x"},
        }},
        "acquisition": {"response_sha256": {"vectors": hashlib.sha256(raw).hexdigest(),
                                            "objdata": "x"}},
    }
    monkeypatch.setattr(exp, "REFERENCE_DIR", tmp_path)
    with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
        exp.verify_snapshot_bytes(manifest)


def test_l12_stale_byte_count_fails_before_hash():
    manifest = exp.load_manifest()
    name = exp.VECTORS_NAME
    broken = json.loads(json.dumps(manifest))
    broken["snapshot"]["files"][name]["bytes"] += 1
    with pytest.raises(RuntimeError, match="byte count"):
        exp.verify_snapshot_bytes(broken)


def test_l12_jump_detector_thresholds_registered():
    assert exp.JUMP_THRESH_M == 100.0
    assert exp.ENVELOPE_KM_PER_DAY == 3.0
    assert exp.BETA_PRIMARY in exp.BETA_BAND_KG_M2
    assert exp.BOOTSTRAP_SEED == 137  # determinism of the decision rule
    assert exp.HEADLINE_NSUB == 8 and exp.HEADLINE_NSUB in exp.LADDER_NSUB


# ----------------------------------------------------------------- L13 ---- #
def test_l13_results_code_hashes_fresh_when_present():
    results_path = _EXP_DIR / "results" / "results.json"
    if not results_path.exists():
        pytest.skip("results.json not generated yet")
    payload = json.loads(results_path.read_text(encoding="utf-8"))
    recorded = payload["results"]["code_sha256"]
    current = exp.code_hashes()
    assert recorded == current, "results.json was produced by different code - regenerate"


def test_l13_results_snapshot_echo_matches_manifest_when_present():
    results_path = _EXP_DIR / "results" / "results.json"
    if not results_path.exists():
        pytest.skip("results.json not generated yet")
    payload = json.loads(results_path.read_text(encoding="utf-8"))
    echo = payload["results"]["provenance"]["reference_snapshot_files"]
    manifest = exp.load_manifest()["snapshot"]["files"]
    assert echo == manifest


# ----------------------------------------------------------------- L14 ---- #
def test_l14_analysis_path_has_no_network_imports():
    src = (_EXP_DIR / "experiment.py").read_text(encoding="utf-8")
    for token in ("urllib", "requests", "socket", "httpx", "http.client"):
        assert token not in src, f"offline doctrine violated by '{token}' in experiment.py"


def test_l14_acquisition_script_guarded_against_overwrite():
    src = (_EXP_DIR / "fetch_horizons_snapshot.py").read_text(encoding="utf-8")
    assert "refusing to overwrite" in src
    assert "REQUEST_SPACING_S = 3.0" in src
