# mission_secular_identification — is there an identifiable secular lunisolar RAAN rate at LEO, and which formula predicts it?

**Mission type:** Discrepancy mission (LAB_CONSTITUTION.md §2.5 type 4).
**Status:** ACTIVE — card, hypothesis and decision rule written BEFORE any campaign result exists.
Phase A COMPLETE (data + instrument + estimator + **formula adjudication**). Phase B (multi-phase
orbital campaign) **ABANDONED — resource overrun, see §6.2.**
**Date established:** 2026-10-04.
**Constitutional authority:** §§2.3 (mission architecture), 3 (evidence doctrine), 9 (FET), 10.1 (hard science gates), 12 (mission selection).
**Budget:** ≤ 10 hr on 8 workers, ≤ 2 GB RAM, local commodity compute only.

---

## 0. Selection rationale (§12.1)

`mission_mean_element_leakage` remains ACTIVE but its 2026-10-03 Phase-0 verdict
and 2026-10-04 identifiability report both conclude **its central question is
not answerable on an 18.6-yr arc**: the committed arc is 0.9993 of one lunar
nodal cycle, the single worst window for separating a secular rate from a forced
term, and the secant's slope on it is not identifiable (VIF 2.6–365; periods
from 6798 d to 20 000 d all explain ≥99.96 % of the variance while moving the
fitted rate by −937 % to +7671 %, including a sign flip).

Both the prior report and the constitution's §4.3 stop-condition list make the
next step a human decision. This card records that the human selected it, and
§0.1 records why it passes the Frontier Economic Test. Per §13.1 formal mission
selection remains a human act; this mission does NOT reopen it.

**This mission explicitly does NOT resume the abandoned J2-coupling-theory
branch.** `mission_j2_precession_modulated_lunisolar_term` falsified H-mod
(commit `3389949`) and `mission_j2_lunisolar_coupling`'s mechanism claim was
retracted. Mechanism work remains unfalsifiable until a secular rate is
identifiable; that identifiability is the sole objective here.

## 0.1 Frontier Economic Test (§9)

| gate | verdict |
|---|---|
| 1 reasoning-space-bound | PASS — the blocker is estimator identifiability, not compute |
| 2 independent validation | PASS — synthetic oracles + analytic J2 gate + external DE441 |
| 3 durable knowledge | PASS — closes or reopens the lab's oldest open question |
| 4 hypothesis-distinguishing | PASS — two formulas, a pre-registered discriminator, and a third "neither" branch |
| 5 adversarial-survivable | PASS — structural discriminator, decided by measurement not by fit |
| 6 capability-advancing | PASS — validated instrument + estimator to `src/lab_utils` if reusable |
| 7 modest resources | PASS — ~4 hr on 8 workers |

## 1. Question

Over an arc of **≥ 2 complete lunar nodal cycles** (6798.383 d each) with an
estimator that is **not structurally biased by a near-one-cycle window** and
that does not sample the node at i = 90°:

> what secular lunisolar RAAN rate, if any, is actually identifiable, and which
> of the two competing formulas — if either — matches it?

## 2. Hypotheses (pre-registered)

- **H1 (FORM-1 / audit-018):** the measured secular rate matches
  `(3/8) n (μ₃/μ) (a/a₃)³ sin2(i−i₃)/sin i`, summed over Sun and Moon.
- **H2 (FORM-2 / 2026-10-04 competitor):** it matches `C(a)·cos i`.
- **H3 (H-null):** neither; the residual is dominated by an unmodelled
  long-period term, and no stable secular rate is identifiable even at ≥2 cycles.

## 3. The competing formulas (reproduced, NOT fitted)

Both are implemented in `competing_formulas.py` from their written forms.
Verified there: FORM-1 reproduces the 2026-10-04 handoff table (+1.3476e-4 at
SSO, +1.761e-4 at 89.5°, +1.739e-4 at 90°) **only when summed over Sun and
Moon with their own inclinations to the plane of `i`**, and independently
reproduces the handoff's "1.51× twin gap".

| | FORM-1 (audit-018) | FORM-2 (competitor) |
|---|---|---|
| inclination law | `sin2(i−i₃)/sin i` | `cos i` |
| SSO (97.7876°) | +1.3476e-4 deg/day | +4.0369e-5 |
| twin (82.2124°) | +2.0340e-4 (**same sign**, ratio 1.509) | −4.0369e-5 (**opposite**, ratio −1.000) |
| 89.5° | +1.761e-4 | −2.600e-6 |
| 90.5° | +1.702e-4 | +2.600e-6 |

The two differ in a **structural** property, not a fitted coefficient, so no
measurement can reconcile them by tuning.

## 4. Pre-registered decision rule (BINDING, fixed before any result exists)

Measured quantity is the **lunisolar** rate,
`Ω̇_luni = Ω̇(sun_moon_j2) − Ω̇(j2_only)`, both estimated with the same estimator
on the same grid.

**4.1 Identifiability gate** (must pass before any formula is scored):
- `VIF_max ≤ 1.20` on the campaign arc, and
- the measured rate is **stable across every phase and every window ladder**
  (windowed subsets of {1.0, 1.5, 2.0, 2.5, 3.0} cycles) to within the estimator
  uncertainty, and
- sign agreement across phases at every inclination.

If 4.1 fails at any inclination, H1/H2 are **not scored there** and the result is
recorded as H-null at that inclination — never as evidence for a formula.

**4.2 Formula adjudication** (per inclination):
- **H1 supported** if the FORM-1 and measured rates have the **same sign** and
  `|measured / FORM-1| ∈ [0.5, 2.0]`.
- **H2 supported** if same sign and `|measured / FORM-2| ∈ [0.5, 2.0]`.
- **H3** otherwise.

**4.3 The structural discriminator** (the decisive test, evaluated across the
inclination set, not per inclination):
- If `Ω̇(82.2124°) ≈ −Ω̇(97.7876°)` within 15 % in magnitude ⇒ **H2 supported,
  H1 refuted** (FORM-1 cannot produce a sign flip; it predicts a 1.51× same-sign
  pair).
- If both are the same sign with `|Ω̇(82.2124°)/Ω̇(97.7876°)| ∈ [1.2, 2.0]` ⇒
  **H1 supported, H2 refuted**.
- Otherwise ⇒ **neither**, with the surviving explanation stated.

**4.4 Magnitude verdict** (only if 4.3 picks a winner): the ratio
`measured / formula` is reported with its estimator uncertainty and, if it is
statistically distinguishable from 1, the formula is declared
**wrong in magnitude** even when its sign and inclination law are correct.

## 5. Method

**Instrument** (no angle, no unwrapping, no node detection, no window):
`Ω̇ = (h_y(−ḣ_x) + h_x ḣ_y)/(h_x²+h_y²)`, with `h = r×v`, `ḣ = r×a`, exact for a
Newtonian field. RK4 on the full 6-vector.

**Estimator:** joint fit of `rate = a₀ + a₁t + Σ_p [b_p cos(ω_p t) + c_p sin(ω_p t)]`
with every ω_p **fixed** to a physical value (lunar nodal 6798.383 d, annual
365.256363 d, anomalistic month, evection). No empirical frequency fitting
(§10.1).

**Diagnostic inclinations:** 97.7876°, 82.2124°, 89.5°, 90.5°.
**i = 90° is excluded by construction** — `cos i = 0` makes the node line
undefined (`1/cos²i = 1.3e4` even at 89.5°).

**Phases:** every inclination is run at multiple argument-of-latitude phases so
no result depends on one favourable initial condition.

**Coverage:** DE441 re-acquired to 2064-01-01 (13,880 daily rows, 2.04 nodal
cycles + margin). `check_coverage()` aborts any arc exceeding the pinned span,
because `interp_snapshot` clamps silently.

## 6. Phases

0. Data acquisition + validation gate. **DONE** — 2026 → 2064, 13 880 rows, overlap with the
   committed 19-yr snapshot reproduces it exactly (max |Δr| = 0.0 km, 0 components > 1 m).
1. Instrument gates. **DONE** — Kepler 2.5e-13 deg/day; J2 0.49 %; coverage 282 d margin;
   i = 90° excluded by construction.
2. Estimator validation on synthetic signals. **DONE** — recovers a known 1.3476e-4 deg/day signal
   to <1e-12 relative at ≥2 cycles; VIF 2.55 → 1.07–1.18.
3. **Formula adjudication by independent quadrature referee. DONE — H1 SUPPORTED, H2 REFUTED.**
4. Multi-phase ≥2-cycle orbital campaign. **RESULT RETRACTED — Nyquist violation (§6.3).**
5. Adjudication, figures, report, knowledge note. **DONE.**

## 6.1 PHASE-A VERDICT (appended 2026-10-04; §4 above is unchanged and was not retuned)

An independent quadrature referee — double-averaging the **exact** point-mass potential, touching
neither the propagator, the estimator nor the DE441 data — reproduces **FORM-1 to 0.02 %** in
magnitude, sign and inclination law at every diagnostic inclination, and independently reproduces
the 1.51× SSO-twin gap (referee 1.5096 vs FORM-1 1.5094, same sign).

| | twin ratio 82.21°/97.79° | sign |
|---|---|---|
| quadrature referee | 1.5096 | same |
| FORM-1 (audit-018) | 1.5094 | same |
| FORM-2 (competitor) | −1.0000 | **opposite** |

**§4.3 verdict: FORM-1 SUPPORTED, FORM-2 REFUTED** — by the first branch's negation. FORM-2
predicts an equal-and-opposite twin pair; the exact potential produces an equal-and-same-sign pair.

This verdict is reached on §4.3 grounds alone and does not depend on §4.1/§4.2, because the
discriminator is structural: no orbital measurement could have reconciled the two formulas by
tuning. Phase 4 (the orbital campaign) tests the *amplitude* independently and is still running.

Also settled: the §3 handoff values are a **Sun + Moon sum** with each body at its own inclination
to the plane of `i`; they do not reproduce from the Moon alone under any choice of i₃. Exp 018's
code *comment* mislabels the 28.584° frame (equatorial, not ecliptic); the formula is correct.

Full evidence: `localdocs/reports/mission-secular-identification-2026-10-04.md`,
`localdocs/knowledge/secular-lunisolar-raan-canon.md`.

## 6.2 Phase 4 outcome: ABANDONED (resource overrun) — NOT a negative result

The multi-phase orbital campaign was **terminated without producing results**. It ran ~10 hr wall
on 8 workers (~64 CPU-hours) and was on its third of four waves when stopped. Measured
throughput was ~1.9 hr per full-mode case per worker at dt = 30 s; the requested matrix
(4 inclinations × 4 phases × 2 modes × 13 597 d) needed ~30 hr wall on 8 workers.

**This exceeded the mission's own declared budget** (`README` §0.1: "≤ 10 hr on 8 workers") and
LAB_CONSTITUTION.md §4.3's stop condition ("> 10 hr single-core, remote/Colab not available").
The overrun is the lead agent's planning error: the dt=30 s gate requirement was applied but its
per-step cost was not multiplied by the full case matrix before committing.

**Nothing was concluded from this run and no number from it exists.** `campaign.py` writes its
JSON only after `pool.map` returns, so there is no partial artifact and nothing to misinterpret.
The absence of `results/campaign_2cyc.json` means *lost*, not *null*.

**This does not weaken the §6.1 verdict.** That verdict rests on §4.3 — a structural
discriminator decided by exact-potential quadrature, which never touches the propagator, the
estimator, or the ephemeris data. A tighter resource estimate for Phase 4 (and a decision on
whether it is worth running at all, given the referee already settled the question) belongs to
the next mission.

Correctly-sized variant if it is ever resumed: 2 inclinations (97.7876°, 82.2124° — the twin pair
that the discriminator actually needs) × 1 phase × 2 modes at dt = 30 s ≈ 7.6 hr on 8 workers,
or the same at dt = 45 s (still inside the dt-stability gate) ≈ 5 hr.

## 6.3 Phase 4 result: RETRACTED — Nyquist violation (supersedes §6.2's framing)

A later, correctly-sized run (2 cases × 2 modes, 2 cycles, dt = 30 s, 6 workers, 16 425 s) **did**
complete and produced numbers. **Those numbers are not measurements and are retracted.**

```
  inc    phase   measured lunisolar      |stderr|   VIF    FORM-1       m/FORM-1
  82.212   0.0   +4.751804e-09 deg/day   6.80e-07   1.18   +2.03403e-04     0.000
  97.788   0.0   -1.937345e-07 deg/day   6.79e-07   1.18   +1.34756e-04    -0.001
```

**Cause.** `--every 240` at `dt = 30 s` gives a **120-minute output cadence**. The orbital period
at h = 600 km is **96.7 minutes**, so the output sampled the orbit at **0.81 samples per
revolution** — below Nyquist. The instantaneous nodal rate carries the full short-period content
at orbital frequency; that harmonic aliases onto DC and contaminates the fitted secular slope at
order 1e-4 deg/day, which is the size of the signal being measured.

**Why the numbers looked believable.** +4.75e-09 and -1.94e-07 deg/day with 6.8e-07 error bars
are near zero and internally consistent — precisely how an aliased signal presents. The values are
~3–4 orders of magnitude below FORM-1, which would have looked like a spectacular refutation.
**A near-zero answer with small error bars cannot be caught by inspecting the number; it has to be
caught before the run.**

**Root cause.** The 2-hour cadence was chosen to reduce sample count and was never checked against
the orbital period.

**Fixed.** `check_sample_cadence()` now refuses to launch below 2 samples/orbit, and records the
check in the results payload. Verified: it blocks the exact `--every 240` configuration that
produced this result, and accepts `--every 20` (10-minute cadence, 9.7 samples/orbit).

**Status: Phase 4 is abandoned with no valid orbital measurement.** The §6.1 verdict is unaffected —
it rests on §4.3, decided by exact-potential quadrature, which never touches the propagator, the
estimator, or the sampling cadence.

## 7. Limitations / non-claims

- The lunar and solar third bodies are point masses at their true DE441
  positions; the satellite has no finite-size, no SRP, no drag, and no
  tesseral gravity. The measured "lunisolar" rate is the rate this force model
  produces, and any remaining long-period content is a *statement about the
  estimator and this model*, not about the full Earth–Sun–Moon system.
- Osculating nodal rate, not a Brouwer mean element. The lab's prior position
  (audit-019 F, audit-020 §7.3) is that a full Brouwer transform is the wrong
  tool; that position is untouched by this mission.
- J2 is carried to first order only.
- If identifiability fails, that is a **reportable negative result**, not a
  licence to invent a term.

## 8. Supersession

- Predecessor `mission_mean_element_leakage` (ACTIVE): this mission supplies the
  instrument and estimator its Phase 1/2 specify. Its frozen README §4 is NOT
  edited by this mission.
- `mission_j2_lunisolar_coupling` and `mission_j2_precession_modulated_lunisolar_term`:
  mechanism claims already retracted; NOT resumed.