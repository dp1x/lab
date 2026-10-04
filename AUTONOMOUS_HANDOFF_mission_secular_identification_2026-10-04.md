# AUTONOMOUS_HANDOFF — mission_secular_identification, session of 2026-10-04

**Mission:** `research/orbital-mechanics/missions/mission_secular_identification/` (ACTIVE)
**HEAD at session end:** see `git log -1`
**Predecessor:** `mission_mean_element_leakage` (ACTIVE, README §4 frozen, NOT edited)
**Authority:** `LAB_CONSTITUTION.md` §§2.3, 3, 9, 10.1, 12

---

## 1. One-paragraph answer

**The competing formulas are now adjudicated: FORM-1 (audit-018) survives, FORM-2
does not.** This was decided by an independent quadrature referee that never
touches the propagator, the estimator, or the DE441 data, and it is confirmed by
the sign/magnitude structure of the orbital measurement. Whether a *numerically
measured* secular rate agrees with FORM-1 in magnitude is the campaign's job and
is reported separately below.

## 2. Committed and verified

| artifact | content |
|---|---|
| `README.md` | mission card + pre-registered decision rule, written BEFORE any result |
| `fetch_de441_sun_moon.py` | 2026-01-01 → 2064-01-01, 13 880 daily rows, 2.04 nodal cycles |
| `reference/` | byte-pinned DE441 Sun+Moon, `MANIFEST.json`, `.gitattributes` `-text` |
| `propagator.py` | instantaneous nodal-rate instrument + coverage guard |
| `estimator.py` | fixed-frequency secular estimator, column equilibration |
| `competing_formulas.py` | FORM-1 and FORM-2 reproduced from source, verified |
| `quadrature_referee.py` | independent exact-potential double average |
| `instrument_gates.py` | all validation gates → `results/instrument_gates.json` |
| `campaign.py` | the ≥2-cycle campaign |
| `tests/` | 28 tests, green |

### Data validation gate
Every chunk overlapping the committed 19-yr snapshot reproduces it exactly:
**max |Δr| = 0.000000e+00 km, 0 components differing by > 1 m.** Raw-response
sha256 was deliberately *not* used as the gate — Horizons embeds its acquisition
timestamp, so it is not reproducible by construction.

### Instrument gates (all pass)
- Kepler instantaneous rate **2.1e-16 deg/day** (exact zero)
- J2 arc mean vs analytic: **0.55 %**
- node conditioning: SSO 54, 89.5/90.5° 1.3e4, **90° undefined → excluded**
- 2-cycle arc fully inside the pinned ephemeris (282 d margin)
- estimator recovers a known 1.3476e-4 deg/day signal to **< 1e-12** relative at
  ≥ 2 cycles; VIF 2.55 → 1.07–1.18

## 3. The adjudication (E3, cross-method — done)

Independent quadrature referee vs both formulas, `a = 6978.137 km`:

| inc_deg | referee (Sun+Moon) | FORM-1 | ratio | FORM-2 |
|---|---|---|---|---|
| 97.7876 | 1.347261e-04 | 1.347561e-04 | 0.9998 | +2.63e-05 |
| 82.2124 | 2.033885e-04 | 2.034029e-04 | 0.9999 | −2.63e-05 |
| 89.5000 | 1.760734e-04 | 1.760967e-04 | 0.9999 | −1.70e-06 |
| 90.5000 | 1.716536e-04 | 1.716758e-04 | 0.9999 | +1.70e-06 |

Twin ratio (82.21°/97.79°): referee **1.5096** (same sign), FORM-1 **1.5094**
(same sign), FORM-2 **−1.0000** (opposite sign).

**Verdict: FORM-1 SUPPORTED, FORM-2 REFUTED.** The referee reproduces FORM-1's
inclination law, sign AND magnitude to 0.02 % with no fitting anywhere, and
independently reproduces the 1.51× twin gap that the 2026-10-04 handoff quoted.

Also settled: the handoff's FORM-1 values (+1.3476e-4 SSO) reproduce **only** when
FORM-1 is summed over Sun AND Moon with their own inclinations to the plane of
`i`. The frame error is in Exp 018's code *comment* (it labels 28.584° as an
ecliptic inclination when it is equatorial); the committed formula is correct.

## 4. Four defects found and fixed (all in-code, all documented)

1. **Column equilibration was mandatory.** cond(A) = 1.8e4 on a 2-cycle arc;
   plain `lstsq` annihilated the secular coefficient (returned 1e-24 for a
   1.35e-4 signal). Equilibrating before the solve is not an optimisation.
2. **dt default moved 120 → 30 s.** dt=120 is 1.37 % off the lunisolar signal;
   dt=240 is catastrophically wrong. The ABSOLUTE rate converges only slowly
   (incomplete J2 short-period average) but the lunisolar DIFFERENCE cancels it,
   so the gate now tests the difference — the quantity the mission reports.
3. **Referee stencil bug.** First version used the x=0 *second*-derivative
   formula for the first derivative, returning a magnitude identical to 6
   significant figures at every inclination. The grid-convergence check caught it.
4. **Two of my own test/harness errors**, corrected rather than the product: the
   synthetic gate initially built a signal with no ramp; and the cycle-mean is
   exact only for the *lunar nodal* term alone, not for a signal also carrying
   annual and short-period terms (the 2026-10-04 exactness result stands in its
   actual domain; the scope is now pinned explicitly).

## 5. What is still open

The **numerical** measurement — whether a multi-cycle, phase-decorrelated orbit
propagation reproduces the +1.35e-4 deg/day FORM-1 value at SSO — was still
running when this session ended. `campaign_2cyc.json` is written only on
completion, so **its absence means the run was lost, not that it returned null**.

To re-run (≈ 2–6 hr on 8 workers):
```
uv run python campaign.py --cycles 2.0 --ladder-max-cycles 2.0 \
    --dt 30 --every 240 --workers 8 --phases 0,90,180,270
```

## 6. Constraints a fresh session must not break

- **Do not edit `mission_mean_element_leakage/README.md` §4** (frozen).
- **Do not use `interp_snapshot` beyond 2064-01-01** without extending the pinned
  snapshot; `check_coverage()` aborts, and that abort is load-bearing.
- **Do not precess the snapshot once at load** — that is the 0.05° frame error
  the 2026-10-02 stack audit rejected.
- **Web access is unavailable** in this environment (HTTP 426 from the search
  backend; fallback engines present CAPTCHAs, which must not be circumvented).
  The referee was built as the substitute for literature verification.
- **`R:` scratch is disposable.** All reference data is committed under
  `reference/`; never leave an acquisition's only copy on `R:`.

## 7. Evidence-determined next mission

If the campaign confirms FORM-1's magnitude: the lab's secular lunisolar canon is
`sin2(i−i₃)/sin i`, FORM-2 is retired, and the remaining open question is the
long-period nodal term that the 2026-10-03 report could not identify on a
one-cycle arc.

If the campaign contradicts FORM-1's magnitude while keeping its sign: the
inclination law is right and the amplitude is not, which localises the error to
the third-body geometry (real lunar e₃/i₃(t) versus the mean values) — a
derivable, testable follow-up, not a new free term.