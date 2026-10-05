# mission_secular_identification — is there an identifiable secular lunisolar RAAN rate at LEO, and which formula predicts it?

**Session type:** autonomous discrepancy mission (LAB_CONSTITUTION.md §2.5 type 4)
**Date:** 2026-10-04
**Authority:** §§2.3, 3 (evidence doctrine), 9 (FET), 10.1 (hard science gates), 12 (mission selection)
**Predecessor:** `mission_mean_element_leakage` (ACTIVE, README §4 frozen, **not** edited)
**Protocol tag:** **P2** — the decision rule was written and committed before any campaign result existed.
**Evidence tier of the adjudication:** **E3** (cross-method).

---

## 0. One-paragraph answer

**FORM-1 (the lab's audit-018 formula) is correct; FORM-2 is refuted.** The two differ in a
*structural* property rather than a fitted coefficient, so no measurement could have reconciled
them by tuning. An independent quadrature referee — which double-averages the **exact** point-mass
potential and touches neither the propagator, the estimator, nor the DE441 data — reproduces
FORM-1's magnitude, sign *and* inclination law to **0.02 %** at every diagnostic inclination, and
independently reproduces the 1.51× SSO-twin gap. FORM-2 predicts the twin pair at 82.21° and
97.79° to be equal and **opposite** in sign; the exact-potential average says they are equal and
**the same sign**. That is the discriminator, and it is settled on theory grounds alone.

Separately, the estimator question is settled: a secular rate *is* recoverable on a ≥2-cycle arc
with fixed physical frequencies, and the estimator recovers a known 1.3476e-4 deg/day signal to
**<1e-12 relative** where the naive OLS slope is wrong by ~100 % and the bare cycle-mean by
~100–500 %. Whether a real orbit reproduces FORM-1's amplitude is the question the campaign was
meant to answer; **it ran out of budget (§5) and is unresolved.** That gap is a planning failure
on the lead agent's part, not a scientific one.

## 1. What was blocking, and what was done about it

`mission_mean_element_leakage` (2026-10-04) identified the blocker precisely: the committed
18.6-yr arc is **0.9993 of one lunar nodal cycle**, the single worst window for separating a
secular rate from a forced term, and the secant's slope there is not identifiable (VIF 2.6–365;
periods from 6798 d to 20 000 d all explain ≥99.96 % of the variance while moving the fitted rate
by −937 % to +7671 %, including a sign flip). It also named the hard prerequisite:

> the pinned DE441 snapshot ends 2045-01-01; `interp_snapshot` **clamps** outside its range, so a
> 2-cycle arc (needing 2063-03-24) would silently freeze the Sun and Moon.

### 1.1 Data (re-acquired durably, in-repo)

| item | value |
|---|---|
| span | 2026-01-01 → **2064-01-01** |
| rows | **13 880** daily, geocentric Sun + Moon, ICRF/TDB, DE441 |
| coverage | **2.0415** lunar nodal cycles |
| committed at | `research/orbital-mechanics/missions/mission_secular_identification/reference/` |
| protection | `.gitattributes` `-text`, both dir-relative and repo-root-relative |

**Validation gate.** Every chunk overlapping the committed 19-yr snapshot reproduces it exactly:
**max |Δr| = 0.000000e+00 km, 0 components differing by > 1 m.**

**The gate was wrong on the first attempt, and the correction is instructive.** I initially pinned
the *raw response sha256* and the very first chunk failed. The cause was not data corruption: a
Horizons response embeds its acquisition timestamp (`Ephemeris / API_USER Tue Sep 1 ...`), so the
raw sha256 is **not reproducible across acquisitions by construction**. Pinning it would fail every
re-acquisition for a reason unrelated to the data. The gate now compares **numeric rows**, which is
what the prior session actually measured. This is the same defect class the lab has hit twice
before (the 86400× units bug, the `rolling_mean` pad bug): an apparent discrepancy that is purely a
measurement artefact.

### 1.2 Silent clamping eliminated

`check_coverage()` is called before every propagation and **aborts** an arc exceeding the pinned
span. It reports the margin explicitly. Tests assert both that a 2-cycle arc passes and that a
3-cycle arc is refused. A separate test asserts that clamping *is* real (the raw J2000 vector at
end+1 s is bit-identical to end+400 d), because a guard against a defect that does not exist is
worthless.

One honest wrinkle: the lab's `T0_S = 820476800.0` sits **0.24 d before** the first snapshot row.
That is pre-existing in every mission in the chain and harmless (it falls inside the first
interpolation segment), but it must not make a healthy arc read as uncovered. The guard therefore
reports `start_slack_days_used` explicitly rather than assuming it.

## 2. The instrument

Every previously committed lunisolar number in the Lab is an OLS slope of an unwrapped node angle
sampled at equator crossings. That estimator needs node detection, needs unwrapping, and has the
whole arc as its window. This mission replaces it with the exact instantaneous rate:

```
Ω̇ = ( h_y·(−ḣ_x) + h_x·ḣ_y ) / (h_x² + h_y²),    h = r × v,   ḣ = r × a
```

No angle, no unwrapping, no node detection, no window. `ḣ = r × a` is exact for a Newtonian field.

### Gates (all pass)

| gate | requirement | measured |
|---|---|---|
| Kepler | exactly zero | **2.5e-13 deg/day** |
| J2 arc mean vs analytic | < 1 % | **0.49 %** (0.990487 vs 0.985641) |
| J2 wobble structure | bounded, structured | rms 0.71× the mean — the known short-period term |
| node conditioning | all diagnostics measurable | SSO 54, 89.5/90.5° 1.3e4 |
| 2-cycle coverage | inside pinned span | yes, **282 d** margin |

**i = 90° is excluded by construction.** The node longitude comes from `h_x : h_y`, whose magnitude
is `|cos i|`. At i = 90° exactly, `h_x = h_y = 0` and the node line does not exist — the one
inclination at which the instrument cannot resolve the quantity being discriminated. This is why the
2026-10-03 report's "correction" of 108.45 % at 90° was a sign reversal (conditioning failure), not
physics. The diagnostic set is **97.7876°, 82.2124°, 89.5°, 90.5°**, which preserves the
`sin i` vs `cos i` discrimination at 6× better conditioning than SSO.

## 3. The estimator, and two ways it would have lied

Fits `Ω̇ = a₀ + a₁t + Σ_p [b_p cos(ω_p t) + c_p sin(ω_p t)]` with every ω_p **fixed** to a physical
value (lunar nodal 6798.383 d, annual 365.256363 d, anomalistic month 29.530589 d, evection
14.765 d). No empirical frequency fitting — LAB_CONSTITUTION.md §10.1 forbids it.

| arc (nodal cycles) | estimator rel. error | VIF | naive OLS | bare cycle-mean |
|---|---|---|---|---|
| 1.0 | 3.4e-13 | 2.550 | 4.08 | 3.07 |
| 1.5 | 9.9e-13 | 1.012 | 0.99 | 9.24 |
| 2.0 | 6.4e-15 | 1.179 | 1.03 | 1.42 |
| 2.5 | 2.0e-13 | 1.002 | 0.22 | 5.63 |
| 3.0 | 4.3e-13 | 1.072 | 0.45 | 0.89 |

*(errors are relative; the contaminants are 5e-3 deg/day, ~37× the signal)*

### 3.1 Column equilibration is mandatory, not an optimisation

The design matrix is badly scaled: over a multi-year arc the `t` column has norm ~2e4 against ~140
for every harmonic column, and the secular direction is nearly parallel to the constant column.
Measured **cond(A) = 1.8e4**. Plain `lstsq` applied to a *purely secular* signal returned **1e-24**
instead of 1.35e-4 — the secular coefficient was annihilated, not merely made imprecise. Dividing
each column by its norm before the solve and rescaling afterwards is algebraically identical and is
what makes the estimator able to recover a known secular rate at all.

### 3.2 The cycle-mean is exact only in its actual domain

The 2026-10-04 session concluded a trapezoid cycle-mean over an integer number of nodal cycles is
exact (3e-10 %). **That result is correct, but its domain is narrower than stated.** A cycle-mean
cancels only harmonics whose period *divides the arc*. Two lunar nodal cycles is not an integer
number of years or of months, so on a signal that also carries an annual term the bare cycle-mean
is wrong by ~140 % of the signal while the fixed-frequency fit is exact to 1e-6. The result stands
for the nodal term alone; the scope is now pinned by an explicit test.

### 3.3 One harness error of my own

The synthetic gate first built `rate = true + amplitude·sin(ωt)` and asked the fit to recover a
*slope* from it. That signal contains no ramp: a constant plus sinusoids is exactly what the
intercept plus harmonic columns already represent, so slope = 0 was the correct answer and the
estimator returned it. The gate was testing the harness, not the estimator. Fixed by integrating
the rate into the node angle the mission actually measures.

## 4. The adjudication (E3, independent of the orbital data)

`quadrature_referee.py` double-averages the **exact** potential `−μ₃/|r_b − r|` over the
satellite's mean anomaly and the third body's mean anomaly, with a 4th-order 5-point stencil and an
explicit grid-convergence check (converged to 8 significant figures from 360×90 upward).

| inc_deg | referee (Sun+Moon) | FORM-1 | ratio | FORM-2 |
|---|---|---|---|---|
| 97.7876 | 1.347261e-04 | 1.347561e-04 | **0.9998** | +2.63e-05 |
| 82.2124 | 2.033885e-04 | 2.034029e-04 | **0.9999** | −2.63e-05 |
| 89.5000 | 1.760734e-04 | 1.760967e-04 | **0.9999** | −1.70e-06 |
| 90.5000 | 1.716536e-04 | 1.716758e-04 | **0.9999** | +1.70e-06 |

**Structural discriminator** — the pre-registered decisive test:

| | twin ratio 82.21°/97.79° | sign |
|---|---|---|
| referee | **1.5096** | same |
| FORM-1 | **1.5094** | same |
| FORM-2 | **−1.0000** | **opposite** |

**Verdict: FORM-1 SUPPORTED. FORM-2 REFUTED.** The referee fits nothing; it agrees with FORM-1's
1.51× twin gap to 0.01 % because both follow from the same geometry.

### 4.1 What the handoff numbers actually were

The 2026-10-04 handoff quotes FORM-1 at SSO as +1.3476e-4 deg/day. A naive reading — FORM-1 applied
to the Moon alone — gives +9.91e-5 (i₃ = 28.584°) or +5.35e-5 (i₃ = 18.294°); neither matches.
The values reproduce **only** when FORM-1 is summed over **Sun and Moon with their own
inclinations to the plane of i**:

```
FORM-1 = Σ₃ (3/8)·n·(μ₃/μ)·(a/a₃)³·sin2(i − i₃)/sin i
         Sun:  i₃ = 23.4393°        Moon: i₃ = 28.5843°
```

Verified to 4 significant figures at all three quoted inclinations (+1.3476e-4, +1.761e-4 at
89.5°, +1.739e-4 at 90°).

**A comment error, not a formula error.** Exp 018's `corrected_secular_lunisolar_raan_rate_rad_s`
labels 28.584° as the lunar "secular mean inclination **to the ecliptic**". It is not: the lunar
orbit is inclined to the ecliptic by 5.145°, so 28.584° is its inclination to the **equator** — the
same frame as the equatorial `i` it is subtracted from. The committed formula is correct; its
comment mislabels the frame. Worth fixing in a remediation commit, but it does not change any number.

### 4.2 Literature status

**Not verified against primary literature.** Web search was unavailable in this environment
(HTTP 426 from the search backend; fallback engines presented CAPTCHAs, which were **not**
circumvented, per `AGENTS.md` Responsible Web Access). The referee was built as the substitute, and
it is a stronger check for this specific question than a citation would be: it tests the formula
against the exact potential it claims to approximate, with no fitted constants anywhere. The
external check that remains genuinely missing is a **published SSO (Landsat-family) measured nodal
drift rate** — that is still un-obtained.

## 5. The multi-phase ≥2-cycle orbital campaign — ABANDONED (resource overrun)

**No measurement was completed. This is a lost run, not a negative result.**

The campaign ran ~10 hr wall on 8 workers (~64 CPU-hours) and was terminated while on its third
of four waves. Measured throughput: **~1.9 hr per full-mode case per worker** at dt = 30 s. The
requested matrix (4 inclinations × 4 phases × 2 modes × 13 597 d = 32 propagations) needed ~30 hr
wall on 8 workers — **three times the mission's own declared budget** and well past
`LAB_CONSTITUTION.md` §4.3's "> 10 hr single-core" stop condition.

This is the lead agent's planning error: the dt = 30 s gate requirement was correctly identified,
but its per-step cost was not multiplied across the full case matrix before committing. The
governance rule that should have caught it — declare a budget, then check the estimate against it
— was satisfied on paper (the card says ≤ 10 hr) and violated in practice.

`campaign.py` writes its JSON only after `pool.map` returns, so **no partial artifact exists and
nothing from the run can be misread**. The absence of `results/campaign_2cyc.json` means lost.

**What this cost the mission:** the independent *amplitude* check. The mission card's §4.2
per-inclination scoring (does `|measured/formula| ∈ [0.5, 2.0]?`) was therefore never evaluated.

**What it did not cost:** the §4.3 structural verdict. That rests on exact-potential quadrature,
which never touches the propagator, the estimator, or the ephemeris data. FORM-1's inclination law,
sign and magnitude to 0.02 %, and FORM-2's refutation by its opposite-sign twin prediction, stand
unchanged and are E3 on their own.

Correctly-sized variant if resumed: 2 inclinations (the 97.7876°/82.2124° twin pair the
discriminator actually needs) × 1 phase × 2 modes at dt = 30 s ≈ **7.6 hr** on 8 workers, or
dt = 45 s (still inside the dt-stability gate) ≈ 5 hr. Alternatively, dt = 120 s would cut it to
~1.9 hr but sits outside the 0.5 % gate and must not be used.

## 6. Defects found and fixed this session

| # | defect | consequence if unfixed |
|---|---|---|
| 1 | raw-response sha256 as the data gate | every re-acquisition fails for a non-data reason |
| 2 | design matrix not column-equilibrated | secular coefficient annihilated (1e-24 for 1.35e-4) |
| 3 | campaign dt default 120 s | 1.37 % error on a 1.35e-4 signal; dt=240 s catastrophically wrong |
| 4 | dt gate tested the ABSOLUTE rate | asserted a meaningless order (2.0 and 7.2) |
| 5 | referee used the 2nd-derivative formula for the 1st | magnitude identical to 6 sig figs at every i |
| 6 | window ladder re-propagated every rung | 5× cost; avoidable — sub-windows are bit-identical truncations |
| 7 | cycle-mean exactness stated too broadly | a wrong claim would have been inherited |
| 8 | synthetic gate built a ramp-free signal | gate tested itself, not the estimator |

Defect 6 deserves a note: fixed-step RK4 advances a state independently of the horizon, so a
2-cycle run is the **bit-identical prefix** of a 3-cycle run (verified:
`np.array_equal(short, long[:n]) → True`). One propagation therefore serves the entire window
ladder — turning a 5× cost into 1×.

## 7. Supersession and what is NOT reopened

- `mission_mean_element_leakage` README §4 remains **frozen and unedited**; this mission supplies the
  instrument and estimator its Phase 1/2 specify.
- `mission_j2_lunisolar_coupling` and `mission_j2_precession_modulated_lunisolar_term`: mechanism
  claims already retracted; **not resumed**. Mechanism work remains unfalsifiable until a secular
  rate is identifiable.
- The 2026-10-03 "99.4 % explained by a ~9175-day sinusoid" identification stays **RETRACTED**.

## 8. Limitations

- The measured rate is the rate *this force model* produces: point-mass Sun and Moon at true DE441
  positions, first-order J2, no tesseral gravity, no SRP, no drag, no finite spacecraft size.
- Osculating nodal rate, not a Brouwer mean element. The Lab's prior position (audit-019 Track F,
  audit-020 §7.3) that a full Brouwer transform is the wrong tool here is untouched.
- FORM-1 and the referee both assume a **circular** third-body orbit. Real lunar e₃/i₃(t) is a
  time-dependent modulation; its effect on the secular structure is a follow-up question, not a
  resolved one.
- The referee's agreement with FORM-1 is at the **quadrupole order in (a/a₃)⁴**; higher orders are
  ~1e-4 relative and are below the comparison's precision.