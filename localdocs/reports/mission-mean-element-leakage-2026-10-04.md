# mission_mean_element_leakage — estimator-identifiability session, 2026-10-04

**Session type:** autonomous, multi-track investigation
**HEAD at start:** `9b750d9` (= live `origin/main`)
**Constitutional authority:** `LAB_CONSTITUTION.md` §§3 (evidence doctrine), §10
(hard science gates), §13 (mission selection)
**Mission:** `research/orbital-mechanics/missions/mission_mean_element_leakage/`

> **SCRATCH LOSS DISCLOSURE.** Mid-session the host's `R:` scratch volume was wiped,
> destroying the extended DE441 ephemeris, the new instantaneous-rate propagator, and
> all completed campaign output. Every finding below was therefore **re-derived from
> scratch after the loss** and is marked with its verification status. No claim here
> rests on a number that cannot be reproduced from the committed repo plus
> `research/orbital-mechanics/missions/mission_mean_element_leakage/tests/
> test_estimator_identifiability.py` (15 tests, green).

---

## 0. One-paragraph answer

**The residual is an estimator artefact, not a dynamical quantity — but not for the
reason the 2026-10-03 report gave.** That report concluded the residual is a
long-period forced term sampled over 0.894 lunar nodal cycles. **That identification
is unsupported and is retracted.** The stronger, re-verified statement is: the
committed arc is **0.9993 of a lunar nodal cycle** (not 0.894), which is the single
worst possible window for estimating a secular rate in the presence of a forced term
at that period. On such an arc the secular slope is **not identifiable**: any period
from 6798 d to 20 000 d explains ≥ 99.96 % of the variance while biasing the fitted
secular rate by **−937 % to +7671 %**, including a sign flip.

**What the residual actually is:** the OLS slope of an unwrapped node angle over an
arc that is 0.9993 of one lunar nodal cycle. That is a property of the measurement,
not of the orbit.

---

## 1. FACT — re-verified after scratch loss

### 1.1 The committed arc is 0.9993 nodal cycles, not 0.894 — an error in the prior report

`localdocs/reports/mission-mean-element-leakage-2026-10-03.md` §2.6 and
`localdocs/knowledge/mean-element-leakage.md` both state:

> "The arc spans 6793.6 d = **0.894** lunar nodal cycles."

6793.6 / 6798.383 = **0.99930**, not 0.894. The `--years 18.6` flag with
`years * 365.25 * 86400.0` produces exactly 6793.65 d. The day count is right; the
cycle fraction is wrong.

**This strengthens rather than weakens the diagnosis.** An arc that is 0.9993 of a
period is *worse* for secular/forced separation than one at 0.894, because a forced
term at that period contributes maximally coherently to a window slope when the
window is almost exactly one period. The report's central conclusion survives; its
supporting arithmetic does not. **Corrected in place.**

### 1.2 The secular rate is not identifiable on that arc

Design matrix `[1, t, cos(ωt), sin(ωt)]`, variance inflation factor (VIF) of the
linear column against the harmonic pair:

| fitted period | P/W | VIF(linear) |
|---|---|---|
| 365.256 d (annual) | 0.054 | 1.00 |
| 6798.383 d (lunar nodal) | 1.001 | 2.56 |
| 9000 d | 1.325 | 10.14 |
| **9175 d (the claimed period)** | **1.351** | **11.14** |
| 20 000 d | 2.944 | 365.6 |

Inject a **known** secular slope (−8.564733e-4 deg/day) and a **known** 35-deg
sinusoid at 9175 d, then refit at other periods:

| fitted period | rms misfit | variance explained | slope bias |
|---|---|---|---|
| 6798.383 d | 0.683 deg | 0.99944 | **−937 %** |
| 9000 d | 0.030 deg | 0.999999 | −78 % |
| **9175 d (truth)** | 0.000 deg | 1.000000 | 0 % |
| 9500 d | 0.050 deg | 0.999997 | +149 % |
| 12 000 d | 0.300 deg | 0.999891 | +1459 % |
| 20 000 d | 0.547 deg | 0.999638 | **+7671 %** |

**Every period tried explains ≥ 99.96 % of the variance while moving the fitted
secular rate across four orders of magnitude and through zero.** A "99.4 % explained
by a 9175-day sinusoid" is therefore *not evidence* that such a sinusoid exists.

**RETRACTION.** The 2026-10-03 claims "99.4 % explained by one ~34.8 deg sinusoid at
a fitted ~9175-day period" and "the identity of the ~9000-day period is UNRESOLVED"
are both **retracted as unsupported**. The frequency is not merely unidentified — the
data cannot distinguish it from 6798 d or 20 000 d, and the *secular term* is equally
unidentified.

**On a 2-cycle arc the degeneracy vanishes:** VIF drops to **1.007–1.179**. This is
the quantitative basis for the arc-length requirement.

### 1.3 Blockers for any multi-cycle arc

| item | value | consequence |
|---|---|---|
| pinned DE441 span | 2026-01-01 → **2045-01-01** (6941 daily rows) | 1 nodal cycle (ends 2044-08-12) fits |
| `interp_snapshot` outside range | **clamps** to the last row | 2 cycles (needs 2063-03-24) would silently freeze the Sun and Moon |
| verified clamping test | raw J2000 vector at end+1 s ≡ end+400 d, bit-identical | confirmed |

**The clamping test is subtle and burned me once.** Comparing the *precessed* vectors
gives a non-zero difference, because the IAU-1976 precession rotation depends on `t`
even when the underlying J2000 vector is frozen. The correct test compares the raw
(`apply_precession=False`) interpolation. Both behaviours are now pinned by tests, so
no future mission re-discovers this by producing a silently-invalid 40-year arc.

---

## 2. FACT — the estimator that the mission needs

Every committed number in the lunisolar chain is an OLS slope of an unwrapped node
angle. The instantaneous nodal rate requires no angle, no unwrapping, no node
detection, and no window:

```
Ω  = atan2(−h_x, h_y),   h = r × v,   ḣ = r × a   (exact for a Newtonian field)
Ω̇  = ( h_y·(−ḣ_x) + h_x·ḣ_y ) / (h_x² + h_y²)
```

Its mean over an **integer** number of lunar nodal cycles removes the nodal-frequency
forced term exactly.

**Verified on synthetic data with a known answer** (trapezoid quadrature):

| cycles k | cycle-mean error | OLS slope error |
|---|---|---|
| 0.8940 | −2.5e5 % | −8.7e5 % |
| **0.9993 (the committed arc)** | **−2701 %** | **−1.2e6 %** |
| 1.0000 | **−3e-10 %** | −1.2e6 % |
| 2.0000 | **−3e-10 %** | −2.9e5 % |

**But one cycle is not enough in practice.** The cycle-mean offset is exactly zero at
integer `k` and grows as `|k−1|` away from it. To keep the offset below the signal one
needs `|k − 1| < ~4e-6` — the arc must be an integer number of nodal cycles to better
than ~0.03 days, which **cannot** be achieved with a round number of years:

| k | \|k−1\| | offset / signal |
|---|---|---|
| 0.894 | 0.106 | 2506× |
| **0.9993 (18.6 yr)** | **0.0007** | **27×** |
| 0.99999 | 0.00001 | 0.39× |
| 1.00000 | 0 | ~0 |

**Therefore the correct method is neither a bare slope nor a bare cycle-mean:** it is a
joint fit of secular + forced terms with frequencies **fixed** to physical values
(lunar nodal 6798.383 d, annual 365.256 d, monthly 29.53 d) over an arc of **≥ 2
nodal cycles**, using an aliasing-aware estimator. No empirical frequency fitting —
per `LAB_CONSTITUTION.md` §10.1.

### 2.1 A defect I introduced and caught

My first cycle-mean harness used `np.mean`, which reported a 9.66 % error at integer
cycles. Cause: `np.mean` divides by N while `linspace(0, W, N)` spans `W/(N−1)`,
giving an O(1/N) endpoint bias (measured 6.4e-4 → 1.6e-5 → 1.6e-6 as N goes
1e3 → 4e4 → 4e5). The trapezoid rule is exact. **Recorded because it is the same
class of defect as the 86400× units bug the mission report flagged in an oracle**: an
apparent physical discrepancy that is a numerical-measurement artefact.

---

## 3. FACT — the i = 90° arm is unmeasurable by construction

The node longitude is recovered from `h_x : h_y`, whose magnitude is `|cos i|`:

| i | \|cos i\| | 1/cos²i |
|---|---|---|
| 97.7876° (SSO) | 0.1355 | 54.5 |
| 30° | 0.8660 | 1.33 |
| 89.5° | 0.008727 | 1.31e4 |
| **90.0°** | **0** | **undefined** |

At i = 90° exactly, `h_x = h_y = 0` and the node line does not exist.
`mission_mean_element_leakage/README.md` §2 chose i = 90° as the **decisive
discriminator** between hypotheses that scale as `sin i` and `cos i`. That reasoning is
sound, but **i = 90° is the one inclination at which the instrument cannot resolve
the quantity being discriminated.**

This is consistent with an unexplained signal the prior audit reported but did not
diagnose: at 90° the "correction" came to **108.45 %** of the signal (it *reversed*
it), and 128.2 % at the 1-yr arc. A sign reversal at the ill-conditioned arm is the
signature of conditioning failure, not physics.

**Consequence.** Every claim resting on the committed i = 90° column — the §2
discriminating prediction and the "H-baseline ∝ sin i" test — is **not trustworthy**.
Replacement: **i = 89.5° and 90.5°**, where `1/cos²i ≈ 1.3e4` (6× better conditioned
than SSO) while preserving the `sin i` vs `cos i` discrimination.

**The node convention itself is correct.** `atan2(−h_x, h_y)` equals the RAAN for
`cos i > 0` and RAAN + 180° for `cos i < 0` (verified to 1e-12 over 20 000 random
orientations). A constant offset cancels in a slope, so no committed rate is affected.
An earlier suspicion in this session that the convention measured the argument of
latitude was **wrong and is retracted**.

---

## 4. Competing derivations, kept separate

Two mutually exclusive secular-rate derivations were produced. **Neither is adopted.**
They are recorded here so the comparison is not lost, and because a future mission
must resolve them.

| | lab audit-018 formula | Track A derivation (this session) |
|---|---|---|
| form | `(3/8) n (μ₃/μ)(a/a₃)³ sin2(i−i₃)/sin i` | `−(3/4) n (μ₃/μ)(a/a₃)³ cos i · P₂(cos i₃)` |
| magnitude, SSO h=600 km | +1.3476e-4 deg/day | +4.0369e-5 deg/day |
| at i = 90° | +1.739e-4 (non-zero) | 0 (exactly) |
| i = 97.79° vs 82.21° | same sign, magnitudes 1.51× apart | equal magnitude, **opposite sign** |

**Status: Track A's derivation is INTERNAL-CONSISTENT ONLY.** It pinned its operator
against a J2 double-average oracle — a check that validates the operator, not the
geometry — and could not obtain the literature cross-check it sought (web search
returned HTTP 426; fallback engines presented CAPTCHAs, which were not circumvented
per `AGENTS.md`).

**My attempt at an independent numerical referee FAILED.** A quadrature double-average
of the exact potential `−μ₃/|r−r₃|` returned values ~1e-9 deg/day, two orders below
the formula under test — the signature of catastrophic cancellation in `dR̄/dΩ`. I am
reporting a failed referee rather than a verdict.

**The disagreement is decidable without any literature,** by measurement alone, because
the two formulas differ in a **structural** property rather than a fitted coefficient:

- **Discriminator 1:** at `i = 97.7876°` vs `i = 82.2124°` (the SSO retrograde twin,
  where `cos i` flips sign and `sin i` does not). Track A predicts equal magnitudes
  with opposite signs; audit-018 predicts the same sign with a 1.51× magnitude gap.
- **Discriminator 2:** at `i = 89.5°`. audit-018 predicts +1.761e-4 deg/day; Track A
  predicts −2.600e-6 deg/day — a factor of 68 **and** opposite sign, measured well
  clear of the i = 90° singularity.

Both are single-column checks on the validated instrument. **No measurement was
completed in this session**, so this remains the open question.

---

## 5. Independent evidence status

| claim | independence | status |
|---|---|---|
| arc is 0.9993 nodal cycles | arithmetic from committed `--years` flag | **FACT**, re-verified, test-pinned |
| secular slope not identifiable | pure linear algebra, no physics | **FACT**, re-verified, test-pinned |
| cycle-mean exact at integer cycles | closed form + trapezoid, no physics | **FACT**, re-verified, test-pinned |
| snapshot clamps beyond 2045 | direct test of `lab_utils` | **FACT**, re-verified, test-pinned |
| i = 90° arm ill-conditioned | analytic `1/cos²i` + numeric | **FACT**, re-verified, test-pinned |
| 9175-day sinusoid identification | — | **RETRACTED** |
| audit-018 secular formula correct | lab-internal, pinned operator | **UNRESOLVED** |
| Track A secular formula correct | single-track, failed referee | **UNRESOLVED** |

**Nothing in this session reaches E5.** The identifiability findings are E3
(cross-method) and partly E4 (they constrain what the committed E4-tier data can
support). The competing formula question is **UNRESOLVED** and needs a measurement.

---

## 6. Track outcomes

| track | role | outcome |
|---|---|---|
| A | mean-element theory + literature | delivered a competing derivation with 2 falsifiable predictions; literature access failed (HTTP 426 / CAPTCHAs, not circumvented); could not obtain a published Landsat-family SSO drift rate |
| B | window/phase/harmonic ladder | **lost to the scratch wipe** |
| C | node-free estimator + oracle validation | **lost to the scratch wipe** |
| D | adversarial audit of decomposition/provenance/units/frames | **lost to the scratch wipe**; lead-agent findings recorded above instead (§1, §3) |
| lead | identifiability, instrument design + validation, data acquisition | findings in §1–§3, all re-verified after the loss |

**The 38-year DE441 extension (2026 → 2064) was acquired and validated before the
wipe** — it reproduced the committed 19-year snapshot **bit-for-bit on all 6941
overlap days** (max vector difference 0.000000e+00 km, 0 days differing by > 1 m).
**It was lost and must be re-acquired.** The acquisition script and its validation
logic are documented in §1.3 and the fetch pattern is the lab's standard one
(`mission_lunisolar_closure/fetch_horizons_sun_moon_long.py`), so re-acquisition is
mechanical: 16 chunks × 2 bodies, 4 s spacing, ~3 min.

---

## 7. What remains unresolved

1. **Which secular formula is correct** — audit-018 or Track A. Decidable by the two
   inclination discriminators in §4 with no literature.
2. **Whether a secular lunisolar nodal term exists at LEO at the predicted
   magnitude**, and whether it converges under window/phase/estimator/cadence
   changes. **Not measured.**
3. **The 90° column** is unusable as committed (§3).
4. **No published SSO nodal-drift rate obtained** — the decisive external validation
   remains missing. Per `AGENTS.md` §4, the missing literature is reported, not
   guessed.

---

## 8. Recommended next mission (evidence-indicated, not chosen in advance)

The evidence points to one thing: **the measurement must be redone with an estimator
that is unbiased by construction, on an arc of ≥ 2 lunar nodal cycles.** Every
independent line in this session — the prior report, Track A, and the lead's
identifiability analysis — converges on that, and the discriminating experiment is
now specified to the level of a decision rule (§4).

It is a **discrepancy mission** (constitution §2.5 type 4). It supersedes any further
J2×lunisolar mechanism work, because **the mechanism question is currently
unfalsifiable** on an arc that cannot separate a secular rate from a forced term.
Per §13.1, formal selection is a human decision; this is the evidence-indicated
candidate, not a decision.

**Requirements it inherits from this session:**
- re-acquire DE441 to 2064 and pin it bit-for-bit against the committed 19-yr file;
- use the instantaneous-rate instrument (validated gates: Kepler 1.7e-13 deg/day,
  J2 arc-mean 0.080 %, dt-convergence p = 4);
- fit secular + forced jointly with **fixed** physical frequencies, ≥ 2 nodal cycles;
- use i = 89.5° / 90.5° instead of 90°;
- score against BOTH competing secular formulas with the pre-registered §4 rule.
