# AUTONOMOUS_HANDOFF — mission_mean_element_leakage, session of 2026-10-04

**Mission:** `research/orbital-mechanics/missions/mission_mean_element_leakage/` (ACTIVE)
**Report:** `localdocs/reports/mission-mean-element-leakage-2026-10-04.md`
**Knowledge:** `localdocs/knowledge/mean-element-leakage.md`
**New tests:** `research/orbital-mechanics/missions/mission_mean_element_leakage/tests/
test_estimator_identifiability.py` (15 tests)

---

## 1. State at end of session

**No measurement was completed.** The session's contribution is entirely
**identifiability analysis + a validated instrument design + two retracted claims**.
The scientific question is exactly where it was on 2026-10-03: **unresolved**, but the
next experiment is now specified to the level of a decision rule.

## 2. The blocking requirement (must be done first)

The pinned DE441 snapshot ends **2045-01-01**. `lab_utils.ephemeris.interp_snapshot`
**clamps** outside that range, freezing the Sun and Moon vectors. One lunar nodal cycle
(6798.383 d) from T0 = 820476800.0 ends 2044-08-12 — **fits**. Two cycles end
2063-03-24 — **does not fit**.

**Action:** re-acquire 2026 → 2064 using the lab's standard pattern
(`mission_lunisolar_closure/fetch_horizons_sun_moon_long.py`): 5-year chunks,
`CENTER='500@399'`, `REF_SYSTEM='ICRF'`, `TIME_TYPE='TDB'`, `VEC_TABLE='2'`,
`STEP_SIZE='1 d'`, 4 s between requests, validate row counts + distance bands.
~16 chunks x 2 bodies, ~3 minutes.

**Acceptance gate (state before running):** the new acquisition must reproduce the
committed 19-year snapshot **bit-for-bit on all 6941 overlap days**
(`max |Δr| = 0.000000e+00 km`, 0 days differing by > 1 m). This was measured on the
extension that was lost to the scratch wipe, so it is known to be achievable. Do not
proceed on a snapshot that fails this gate.

**Trap:** the clamping test must compare the **raw** (`apply_precession=False`)
interpolation. The IAU-1976 precession rotation depends on `t`, so a frozen J2000
vector still rotates and a precessed comparison falsely reads "not clamped". Both
behaviours are pinned in the new test file.

## 3. The instrument (validated, not committed)

Exact instantaneous nodal rate — no angle, no unwrapping, no node detection, no window:

```
h = r × v ;  ḣ = r × a  (exact for a Newtonian field)
Ω  = atan2(−h_x, h_y)
Ω̇  = ( h_y·(−ḣ_x) + h_x·ḣ_y ) / (h_x² + h_y²)
```

Validation gates it must reproduce (measured before the wipe):

| gate | requirement | measured |
|---|---|---|
| Kepler only | Ω̇ = 0 | max \|Ω̇\| = 1.7e-13 deg/day |
| J2 only, 60 d | arc-mean = analytic `−1.5 n J₂(R/a)²cos i` | 0.0797 % |
| finite-difference of Ω | agrees with instantaneous Ω̇ | median 1.7e-4 deg/day |
| dt 240→30 s | p = 4 | ratios 32, 38 |

**Implementation note:** the RK4 must be applied to the full 6-vector
(`f` returns `[v, a(t,x)]`), exactly as `lab_utils.integrators.rk4_step` does. A
hand-rolled version that updates position with stage velocities but velocity with raw
accelerations is **wrong** and drove `a` negative within one sample.

## 4. The estimator decision (settled — do not re-litigate)

- A **bare OLS slope** is wrong by ~1e6 relative to the signal on any arc < 2 cycles.
- A **bare cycle-mean** over a *round* number of years is wrong by ~27x the signal,
  because the offset grows as `|k − 1|` and an 18.6-yr arc is 0.9993 cycles.
- A cycle-mean is **exact** (3e-10 %) at an **integer** number of nodal cycles, using
  the **trapezoid** rule. `np.mean` has an O(1/N) endpoint bias and must not be used.

**Required method:** joint fit of secular + forced terms with frequencies **fixed** to
physical values (lunar nodal 6798.383 d, annual 365.256 d, monthly 29.53 d), over
**≥ 2 nodal cycles**, using an aliasing-aware estimator. No empirical frequency
fitting (`LAB_CONSTITUTION.md` §10.1).

## 5. Pre-registered discrimination rule (fixed before any measurement)

Two competing secular formulas, differing in a **structural** property, hence
decidable by measurement alone:

| | lab audit-018 | Track A (2026-10-04) |
|---|---|---|
| form | `(3/8)n(μ₃/μ)(a/a₃)³ sin2(i−i₃)/sin i` | `−(3/4)n(μ₃/μ)(a/a₃)³ cos i·P₂(cos i₃)` |
| SSO h=600 km | +1.3476e-4 deg/day | +4.0369e-5 deg/day |
| at i = 90° | +1.739e-4 | 0 exactly |
| i=97.79° vs 82.21° | same sign, 1.51x magnitude gap | equal magnitude, **opposite sign** |
| at i = 89.5° | +1.761e-4 | **−2.600e-6** (68x smaller, opposite sign) |

**Rule:**
- `⟨Ω̇_luni⟩(82.2124°) = −⟨Ω̇_luni⟩(97.7876°)` within 15 % → **Track A supported**.
- same sign at 97.79° and 82.21° with `|ratio − 1|` in [1.2, 2.0] → **audit-018 supported**.
- anything else → **neither**.

**Use i = 89.5° and 90.5°, never 90°** (1/cos²i conditioning: 54.5 at SSO, 1.3e4 at
89.5°, undefined at 90°).

## 6. Known-correct, already-done work (do not redo)

- `interp_snapshot` clamping — confirmed, test-pinned.
- Arc arithmetic — the `--years 18.6` arc is 6793.65 d = 0.9993 nodal cycles.
- VIF / non-identifiability on that arc — computed, test-pinned.
- Cycle-mean exactness and the 1/W law — computed, test-pinned.
- Node convention `atan2(−h_x, h_y)` = RAAN (mod 180) — correct, test-pinned.
- `mission_mean_element_leakage/README.md` §4 is **frozen** and was NOT edited.

## 7. What is NOT established

- Whether a secular lunisolar nodal term exists at LEO at the predicted magnitude.
- Which of the two secular formulas is correct.
- Any published SSO (Landsat-family) measured nodal-drift rate — still missing, and
  still the best external validation available. Web search returned HTTP 426 and the
  fallback engines presented CAPTCHAs, which were **not** circumvented.
- Tracks B, C and D of this session were **lost to the `R:` scratch wipe** and must be
  re-run; their questions are listed in the report §6.

## 8. Infrastructure

`R:` scratch is empty and must be recreated before any long arc. Note the wipe: **do
not leave a long campaign's only copy of acquired data on `R:`** — acquire, verify,
and commit reference data to the repo under `reference/` with `-text` gitattributes
before launching a multi-hour campaign.
