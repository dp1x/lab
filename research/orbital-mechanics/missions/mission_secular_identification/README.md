# mission_secular_identification — is there an identifiable secular lunisolar RAAN rate at LEO, and which formula predicts it?

**Mission type:** Discrepancy mission (LAB_CONSTITUTION.md §2.5 type 4).
**Status:** ACTIVE — card, hypothesis and decision rule written BEFORE any campaign result exists.
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

0. Data acquisition + validation gate (done, pre-result).
1. Instrument gates: Kepler Ω̇=0; J2 arc-mean vs analytic; dt convergence p=4.
2. Estimator validation on synthetic signals with a known answer.
3. ≥2-cycle campaign, parallel over inclinations × phases × modes.
4. Adjudication, figures, results.json, report, knowledge note.

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