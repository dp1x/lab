# Post-Roadmap Strategic Probe — Research Lab Operating Model

> **Probe date**: 2026-08-31
> **Author**: Autonomous strategic-probe session (read-only against repo; no
> experiments implemented, no scientific code modified).
> **Status**: PLANNING ONLY. Not a directive. The probe asked "what should this
> laboratory become now that it has outgrown the assumptions under which it was
> originally designed?"; this document is the canonical concise answer.
> **FACT / INFERENCE / RECOMMENDATION** separation is enforced throughout.

---

## 0. Executive summary (FACT)

The Computational Research Laboratory has reached a structural inflection
point. Twenty numbered experiments are complete; the lab's original
finite roadmap (orbital-mechanics flagship Phase 2) has effectively been
consumed and surpassed; the most recent work (017→018→019→020) is a
four-experiment investigation into a single scientific question (the
lunisolar RAAN secular rate) and that investigation is still **unresolved**
at the 1-year numerical arc.

The laboratory is now operating at a scale and depth that the original
"finite roadmap" was never designed to manage. The next experiment, the
next delegation pattern, the next governance rule, and the next
organizational model all need explicit design. Continuing with the
existing fixed-roadmap model would consume the lunisolar thread indefinitely
without yielding durable improvement and would block the lab from
exercising the capabilities it has accumulated.

**Recommended operating model after Experiment 020** (full reasoning in §11):

- Replace "Experiment NNN" as the unit of work with a **research mission**
  (`mission_<topic>`) selected dynamically from a maintained
  capability/uncertainty frontier.
- Each mission contains one or more **experiments** (the existing artifact
  format stays unchanged); remediation, capability-graduation, and
  exploratory work are first-class mission types, not numbered experiments.
- Numbers persist only as historical ordinals for completed work; new
  work uses topic-named mission IDs.

**Recommended next mission** (full ranking in §13):

The single highest-leverage move is **NOT** continuing the lunisolar
thread with a 5-yr DE441 arc (the 020 recommendation), but rather a
**Capability Audit + Frontier Expansion Mission** that pauses the
lunisolar investigation at its current empirical plateau, runs one
final multi-year harmonic-regression cross-check (delegated, cheap,
deterministic), and then pivots to a portfolio of new research programs
in: (a) **end-to-end SSO station-keeping as a closed validation chain**,
(b) **a new domain — either celestial-mechanics three-body or
energy-systems battery degradation — chosen for high
information-value / low marginal compute**, and (c) **a
verification-architecture experiment that turns the 020 lessons
(thin-sample estimator fragility, frame-mismatch bugs, finite-window
bias) into a permanent reusable estimation library**.

---

## 1. What the laboratory was supposed to be (FACT)

The charter (`localdocs/charter.md`) defines the laboratory as a
research environment producing reproducible experiments, validated
models, simulations, documentation, reusable software, a knowledge
base, and open-source-quality artifacts. The phrase is exact:

> "The goal is not code generation. The goal is a continuously improving
> research system."

The original roadmap (`localdocs/roadmap.md`) explicitly defined a
**finite Phase 2 of orbital mechanics with experiments 002–014+**:

| # | Experiment | Question |
|---|-----------|----------|
| 002 | Kepler orbit validation | Newtonian gravity → ellipses |
| 003 | Kepler equation solvers | Newton vs bisection vs series |
| 004 | Hohmann transfer | Least-fuel two-burn transfer |
| 005 | Bi-elliptic vs Hohmann | Crossover radius law |
| 006 | Plane-change maneuvers | Inclination change cost |
| 007 | Gravity assist | Flyby velocity boost |
| 008 | Ground tracks | Path over Earth |
| 009 | J2 precession | Nodal drift |
| 010 | Orbit decay | Atmospheric drag |
| 011 | Lagrange points | Three-body stability |
| 012 | Orbit classes | SSO/Molniya/GTO families |
| 013 | JPL ephemeris validation | Full propagator vs Horizons |
| 014 | Eclipse timing / launch windows | Event-driven geometry |
| 015+ | "Eclipse-aware station-keeping, ground-track targeting…" | "each seeds the next" |

Phase 1 (numerics, Exp 001) was explicitly labelled "foundation".
Phases 3, 4, 5 (energy systems, computer architecture, cybersecurity)
were sketched as future pillars but never detailed.

**The original intended scope was therefore a finite orbital-mechanics
flagship of ~14–16 experiments**, with the explicit understanding that
numerics foundation work would be limited to what directly serves the
flagship.

## 2. How the laboratory expanded from the finite roadmap (FACT + INFERENCE)

Twenty experiments are complete (001–020). The path beyond 014 was
**not** "roadmap inertia" — every post-014 experiment was triggered by a
specific defect or knowledge gap identified in the prior experiment.
The expansion was driven by five distinct forces, each verifiable in
the git history.

### 2.1 Audit-driven expansion (Exp 015 remediation)

**FACT**: Exp 015 (dawn-dusk SSO launch-window targeting, 2026-08-29)
made a public scientific claim — "LST at the ascending node drifts
through 24 h/year at the sidereal-solar differential (4 min/day)" — that
was **wrong** (frame/convention error confusing sidereal rotation rate
with SSO nodal rate). The 8-track audit (audit-015) retracted the
claim.

**INFERENCE**: The remediation contract produced Exp 016 (SSO LST-drift
correction, 2026-08-30) to provide a first-principles derivation of
the actual drift rate. Exp 016 then made a derived claim about a
lunisolar "upper bound" that itself turned out to be wrong. This is
the genesis of the lunisolar chain.

### 2.2 Discrepancy-driven expansion (Exp 017 → 020)

**FACT**: The Exp 016 Lunisolar upper-bound formula was tested in
Exp 017 (lunisolar upper-bound verification, 2026-08-30) and found to
over-estimate the numerical by a signed ratio of ~170× at h=600 km
(opposite sign). The 8-track audit (audit-018) showed the formula was
mathematically wrong in three compounded ways (wrong radial scale,
wrong geometric factor, wrong sign at SSO retrograde). Exp 018 fixed
the formula and reported a 9.78× residual; Exp 019 attributed the
residual to finite-window bias and produced a 27× extrapolation; Exp
020 found the 019 extrapolation is not theoretically justified.

**INFERENCE**: This is a textbook discrepancy-driven scientific
investigation: hypothesis → test → fail → revise → retest → revised
hypothesis → retest. The chain is legitimate science, not
gold-plating.

### 2.3 Capability-driven expansion (Exp 011 → 013 → 014)

**FACT**: Exp 011 (Lagrange points, 2026-08-22) graduated the
first shared infrastructure (`src/lab_utils/integrators.py`,
`src/lab_utils/orbits.py`) because a second consumer (the CR3BP
rotating-frame dynamics) emerged. Exp 012 (orbit classes, 2026-08-23)
added `j2_rhs` to `src/lab_utils/orbits.py` (second consumer). Exp 013
(JPL ephemeris validation, 2026-08-24) introduced the
"byte-pinned Horizons snapshot + `-text` gitattributes + offline
analysis" pattern, which was reused by Exp 014 (Sun snapshot), Exp 017
(Moon snapshot), and Exp 020 (multi-year acquisition plan).

**INFERENCE**: The capability is real, the reuse is real, and it
materially changed the laboratory's external-validation capacity
(byte-pinned JPL data is what makes the lab a verification-quality
research organization rather than a numerical-experiment-only one).

### 2.4 Mission-composition expansion (Exp 014 → 015)

**FACT**: Exp 015 (dawn-dusk SSO launch-window targeting) is
described in AGENTS.md as "the first end-to-end multi-constraint
mission analysis" — it composes four prior capabilities (SSO
inclination lock, LST-at-node, J2 nodal drift, eclipse event-finder).

**INFERENCE**: Exp 015 is qualitatively different from Exp 002–014 in
that it is no longer "validate one physics formula against one
closed-form answer"; it is "compose several validated physics modules
into a multi-constraint mission feasibility analysis". The
composition succeeded enough to make a (wrong) claim, and the
subsequent audit chain.

### 2.5 Audit-culture expansion (audit-015, audit-017, audit-018, audit-019, audit-020)

**FACT**: Five separate 8-track independent audits have been run
(2026-08-29 through 2026-08-31), each consuming 8 parallel
sub-agents and producing ~30–60 KB of audit reports. The audits are
the primary reason three published findings (Exp 015 LST claim, Exp
016/017 Lunisolar closed form, Exp 019 extrapolation) have been
remediated rather than left as published errors.

**INFERENCE**: The 8-track audit pattern is itself a major capability.
It is the lab's mechanism for catching its own errors before they
become durable. This capability is **not** documented in the
original charter or roadmap; it emerged organically.

## 3. Capabilities now established (FACT + status classification)

For each capability, status is one of:

- **VERIFIED** — confirmed by ≥1 experiment against external reality,
  with peer-review-quality audit.
- **VERIFIED-WITH-LIMITATION** — confirmed within an explicit
  pre-declared envelope; the limitation is documented.
- **PARTIALLY-VERIFIED** — established but with known gaps.
- **UNRESOLVED** — open scientific question.
- **DEPRECATED** — formerly VERIFIED, now known to be wrong;
  preserved with deprecation warnings.
- **REUSABLE-INFRASTRUCTURE** — graduated to `src/lab_utils/` and
  exercised by ≥2 consumers.

### 3.1 Numerical integration (VERIFIED)

- **Fixed-step RK4** at design order p≈4.5 (verified by self-convergence
  ladders in Exp 001, 009, 010, 017, 018, 020 across multiple
  dynamical regimes).
- **Symplectic Euler + velocity Verlet** verified for conservative
  dynamics (Exp 001).
- **Reusable**: `src/lab_utils/integrators.py` (`rk4_step`,
  `rk4_propagate`), 2nd consumer Exp 011 graduated it.

### 3.2 Orbital element canon (VERIFIED)

- **Keplerian 6-element state ↔ Cartesian (r, v)** is machine-precision
  exact (Exp 002, 003) — `state_to_elements` / `elements_to_state`
  in `src/lab_utils/orbits.py`.
- **J2 secular nodal/apsidal rates** (Exp 009) match closed-form to
  ~0.5% at LEO; the model-order residual is documented as
  mean-vs-osculating + small-divisor near `i_crit`.
- **REUSABLE**: 5+ consumers across Exp 002–020.

### 3.3 Conservative perturbed dynamics (VERIFIED)

- **Hohmann, bi-elliptic, plane-change** closed forms vs RK4 to
  machine precision (Exp 004–006).
- **Patched-conic gravity assist** analytical ceilings reproduced
  (Exp 007; Voyager 1/2 Jupiter + Saturn anchors).

### 3.4 Dissipative dynamics (VERIFIED-WITH-LIMITATION)

- **Atmospheric drag (Exp 010)**: decay law rediscovered vs
  erfi/quadrature oracles; 3.6 m max residual over 500 revs;
  benchmark PASS decade band.
- **Limitation**: the drag model is exponential-atmosphere + spherical
  Earth; actual atmosphere is Jacchia/Bowman, terrain-resolved. The
  drag regime in Exp 010 is realistic for first-order mission design
  but not for precision orbit determination.

### 3.5 Rotating-frame / three-body dynamics (VERIFIED)

- **CR3BP equilibria + Jacobi + Routh** (Exp 011) reproduced to
  mpmath precision (≤7.3e-16); EML1/2/3 + SEL1/2 mission anchors
  within 0.03%.

### 3.6 External ephemeris validation (VERIFIED-WITH-LIMITATION)

- **Byte-pinned JPL Horizons DE441** snapshots for ISS, Sun (1 yr,
  Exp 014), Moon (1 yr, Exp 017), with `-text` gitattributes and
  offline-deterministic analysis (Exp 013).
- **Limitation**: snapshots are 1 year; multi-year acquisition is
  feasible (Track 6 design in audit-020) but not yet executed.

### 3.7 Event detection (VERIFIED)

- **Closed-form shadow geometry event finder** on analytic Kepler
  states (Exp 014) decouples event error from integration step;
  agreement with pinned ISS trajectory to 5.5–13.5 s for first 4
  events.

### 3.8 Multi-constraint mission analysis (PARTIALLY-VERIFIED)

- **SSO launch-window feasibility** (Exp 015) was structurally
  correct (cardinality, equinox dominance, sensitivity matrix,
  i_SSO anchors) but the headline physical claim was wrong (audit
  retraction). The Exp 016/017/018 chain produced a defensible
  station-keeping budget via corrected theory.
- **Limitation**: only one mission type has been validated; the
  general capability of multi-constraint composition is plausible
  but not yet proven across multiple mission types.

### 3.9 Lunisolar secular theory (PARTIALLY-VERIFIED, with one UNRESOLVED)

- **Corrected doubly-averaged quadrupole secular formula** (Exp 018):
  matches 1-year numerical in sign and within 9.78× at LEO SSO;
  within 2.81× at i=90°. **VERIFIED-WITH-LIMITATION** for
  short-arc asymptotic behavior.
- **The Lunisolar RAAN secular limit at W → ∞** (the actual long-term
  asymptotic rate): **UNRESOLVED**. The 019 polynomial-in-1/W
  extrapolation (+0.0036 deg/day, 27× the corrected formula) is NOT
  validated; the 020 1-year multi-phase ensemble reproduces the 018
  finding (~9.3× ratio). The 018 corrected formula gives the correct
  SIGN but under-estimates the 1-year numerical rate by ~10×. The
  secular limit at W → ∞ cannot be established from 1-year data
  regardless of estimator sophistication.

### 3.10 Long-period estimator design (PARTIALLY-VERIFIED)

- **Theory-driven harmonic regression** (Exp 020 Track 3) is shown
  by synthetic oracle to recover known secular to machine precision
  (7e-16 deg/day bias).
- **On real data, the same estimator is fragile** (8.89e-1 to
  1.11e+0 deg/day swings); unmodelled short-period content aliases
  into long-period harmonics.
- **VERIFIED in principle, UNVERIFIED in practice** until a multi-year
  arc is run.

### 3.11 Frame-convention and time-system discipline (VERIFIED, after remediation)

- **ECI ↔ ECEF ↔ lat/lon** (Exp 008, 015, 020) consistent to < 1 deg
  after the audit-015 / audit-019 remediation commits.
- **IAU-1976 precession** correctly applied after the audit-019
  `_rot3` transpose fix.
- **Sun/Moon ICRF/TDB geocentric** consistent with byte-pinned
  Horizons.

### 3.12 8-track independent audit pattern (REUSABLE-INFRASTRUCTURE)

- **8 parallel sub-agents, each on a distinct epistemic track**,
  producing per-track outputs that are integrated by the lead agent.
- **Used** in audit-006 (closure), audit-015 (LST drift), audit-017
  (lunisolar verification), audit-018 (lunisolar discrepancy),
  audit-019 (long-period), audit-020 (secular limit).
- **Limitation**: cost is high (~10–30 min wall-clock per audit);
  should not be invoked for trivial corrections.

### 3.13 Byte-pinned external data acquisition (REUSABLE-INFRASTRUCTURE)

- **Pattern**: Horizons API → sha256-pinned JSON with `-text`
  gitattributes → MANIFEST.json under `reference/` → offline
  deterministic analysis.
- **Used** in Exp 013 (ISS), Exp 014 (Sun), Exp 017 (Moon).
- **NOT YET EXERCISED** for multi-year (≥5-yr) data.

## 4. Scientific questions genuinely unresolved (FACT + INFERENCE)

The following questions are open as of 2026-08-31. Each is sized by
estimated information value vs required resources.

### 4.1 Lunisolar RAAN secular limit at W → ∞

**STATUS**: UNRESOLVED at 1-yr arc. The 018 corrected formula under-estimates
the 1-yr numerical by ~9.3× at LEO SSO; the cause is documented as
mean-vs-osculating finite-window bias, but the exact asymptotic
correction is not known.

**WHY IT MATTERS**: This is the third-largest secular perturbation on
LEO satellites; an order-of-magnitude error in its prediction directly
affects station-keeping Δv budgets and constellation design.

**WHAT WOULD RESOLVE IT**: A 5-year DE441 arc at i_sso + i=90° + i=30° +
the 2-window phase-locked estimator (audit-020 Track 6) would cancel
the lunar nodal modulation and bound the secular to ±10% of the corrected
formula at ~3.5 hr single-core compute.

**EXPECTED INFORMATION VALUE**: High — closes the open 020 question;
produces a publication-quality lunisolar secular benchmark.

**EXPECTED RESOURCE COST**: Low — 5–10 hr single-core; ~2 MB multi-year
ephemeris; ~40 new tests.

**VERDICT**: SHOULD be addressed but **bundled** with broader Lunisolar
Capability Graduation (see mission ranking in §13).

### 4.2 Lunisolar coupling with J2 in the secular regime

**STATUS**: UNVERIFIED. The corrected formula assumes a ≪ a₃ and
quadrupole-only; J2 × Lunisolar cross-coupling is documented in Exp 020
Track 3 as the dominant unmodelled content in the i_sso 1-yr fit
(J2 cos i ≈ 0 at i_sso, but the derivative is non-zero).

**EXPECTED INFORMATION VALUE**: Medium — explains ~10–30% of the
019/020 residual.

**VERDICT**: Same mission as 4.1; not worth a standalone experiment.

### 4.3 Operational Δv budgets under lunisolar + J2 + drag

**STATUS**: PARTIALLY VERIFIED (Exp 016 with corrected formula).
The operational envelope (Sentinel-1 ~15 m/s/yr, Landsat ~5–15 m/s/yr)
is consistent with the corrected formula but **not** validated by
end-to-end numerical propagation of a realistic spacecraft.

**VERDICT**: Could be folded into a station-keeping benchmark mission.

### 4.4 SSO ground-track repeat-cycle targeting

**STATUS**: UNVERIFIED. Exp 015 explored dawn-dusk SSO feasibility but
did not solve the **sub-cycle repeat-ground-track** problem (e.g., the
Landsat-7 16-day repeat at 705 km). The lab has the machinery (Exp
008 ground tracks, Exp 012 SSO lock, Exp 015 launch windows) but has
not composed them.

**EXPECTED INFORMATION VALUE**: Medium–High. Repeat-ground-track
targeting is a real mission-design problem with published closed-form
theory (Lansdale, 1991) that could be reproduced and audited.

**EXPECTED RESOURCE COST**: Low — 1–2 days; ~30 tests; no new
infrastructure.

**VERDICT**: Strong candidate for next mission; recomposes existing
capabilities rather than adding new physics.

### 4.5 Three-body Earth–Moon mission trajectories (beyond libration points)

**STATUS**: PARTIALLY VERIFIED (Exp 011 establishes L1–L5 + Jacobi).
Real-world mission design (ARTEMIS, GRAIL, CAPSTONE) requires
**low-energy transfers**, **ballistic lunar capture**, **NRHO
maintenance**, and **B-plane targeting from Earth-escape**. None of
these have been attempted.

**EXPECTED INFORMATION VALUE**: High — opens an entirely new
mission-design domain with mature external validation (GRAIL, CAPSTONE
tracking data).

**EXPECTED RESOURCE COST**: Medium — needs JPL ephemeris + lunar
gravity model + RTN/CR3BP hybrid propagator.

**VERDICT**: Strong long-term candidate; not a fast follow-up because
it requires new infrastructure (lunar gravity harmonics).

### 4.6 Beyond LEO: cislunar + interplanetary

**STATUS**: NOT YET DEVELOPED. The lab has been 100% LEO-focused.
No MEO, GEO station-keeping, cislunar, or interplanetary work.

**VERDICT**: Major new domain; would require deep investigation of
mission priorities before commitment.

### 4.7 Non-orbital-mechanics domains (energy, computer architecture, cybersecurity)

**STATUS**: PLANNED ONLY per the original charter Phases 3–5.
Nothing has been built; the lab has no track record in these areas.

**VERDICT**: Out-of-scope for current capabilities; would require
significant new infrastructure. The "second pillar" energy promise
from the original charter has not materialized in 31 days of
operation; the deprioritization is implicit.

## 5. Capability / uncertainty matrix (FACT)

| Capability | Status | Cost to graduate further |
|-----------|--------|--------------------------|
| Numerical integration | VERIFIED | none needed |
| Orbital element canon | VERIFIED | none needed |
| Conservative perturbed dynamics | VERIFIED | none needed |
| Dissipative dynamics | VERIFIED-WITH-LIMITATION | higher-fidelity drag = 2–4 wk |
| Rotating-frame / 3-body | VERIFIED | NRHO + lunar gravity = 2–4 wk |
| External ephemeris validation | VERIFIED-WITH-LIMITATION | multi-year = 1–2 days |
| Event detection | VERIFIED | none needed |
| Multi-constraint mission analysis | PARTIALLY-VERIFIED | one more mission type = 1 wk |
| Lunisolar secular theory | PARTIALLY-VERIFIED + UNRESOLVED W→∞ | 5-yr arc = 1–2 wk |
| Long-period estimator design | PARTIALLY-VERIFIED | multi-year validation = 1–2 wk |
| Frame-convention discipline | VERIFIED (post-remediation) | ongoing vigilance |
| 8-track independent audit | REUSABLE-INFRASTRUCTURE | mature; document the pattern |
| Byte-pinned external data | REUSABLE-INFRASTRUCTURE | graduate to `lab_utils/`? |

The highest-leverage graduation is **multi-year DE441 acquisition**,
which unlocks both §4.1 (lunisolar W→∞) and §4.2 (L+S coupling). The
second-highest is **repeat-ground-track mission analysis** (§4.4),
which exercises the lab's multi-constraint composition capability on
a real-world mission class.

## 6. The 017→018→019→020 chain: what is established vs unresolved

(FACT, drawn from `audit-019-synthesis-2026-08-30.md`,
`audit-020-track-{1..8}-*.md`, and `localdocs/knowledge/lunisolar-secular-limit-020.md`)

### 6.1 Established (VERIFIED)

- **Correct sign** of Lunisolar RAAN perturbation at SSO retrograde
  (both corrected formula and 1-yr numerical give prograde).
- **Order of magnitude** of the rate at LEO SSO (~10⁻³ deg/day;
  ~10⁻⁴ deg/day for the corrected doubly-averaged quadrupole).
- **The 016/017 "Vallado Eq. 9-46" closed-form is mathematically wrong**;
  it is preserved with DeprecationWarning.
- **Estimator (f) (theory-driven harmonic regression) on synthetic
  data recovers known secular to machine precision** (7e-16 deg/day
  bias).
- **The 019 polynomial-in-1/W extrapolation is NOT theoretically
  justified** (audit-020 Track 3) — it is reported only as a
  diagnostic.
- **OLS bias scaling** for fast harmonics is O(1/W²), not O(1/W);
  for slow harmonics (e.g., 18.6-yr lunar nodal) the bias
  asymptotes to a constant O(A ω sin φ).
- **Estimator (g) (secant) has structural bias for slow harmonics**;
  it cannot be the headline at any W.

### 6.2 UNRESOLVED (open scientific question)

- **The Lunisolar RAAN secular limit at W → ∞** at LEO SSO. The 018
  corrected formula + harmonic regression on a 5-yr arc would bound
  it, but it has not been done.
- **The exact evection / variation amplitudes** in the osculating Ω
  at i_sso (Track B estimates O(0.05 deg); 019 FFT threshold O(0.005
  deg)). Cannot be resolved without re-extracting from the 019
  propagation cache.
- **Whether the corrected formula is the right asymptotic prediction
  at all**. The sign is right; the magnitude at finite W has a
  ~10× positive bias; whether that bias fully decays to 0 or
  asymptotes to a constant O(A_horizon ω_horizon) is unknown without
  a multi-year arc.

### 6.3 DEPRECATED

- `closed_form_lunisolar_raan_rate_rad_s` (017) — wrong.
- `luni_solar_raan_rate_rad_s` (016) — wrong.
- The 019 Ω̇_fit(W → ∞) extrapolation — theoretically unjustified.

### 6.4 INFERENCE (working hypothesis)

The evidence supports: **the 018 corrected formula gives the correct
asymptotic MEAN drift; the 1-yr numerical OSCULATING fit is biased
high by O(1/W²) for fast harmonics + O(A ω sin φ) for slow harmonics;
the secular limit is the 018 corrected value with an as-yet-unquantified
correction term from J2 × Lunisolar coupling + lunar-evection + variation
+ annual solar forcing aliases.** The multi-year arc + harmonic regression
+ 2-window phase-locked estimator would confirm or refute this hypothesis
quantitatively.

## 7. Frontier economic test (RECOMMENDATION, drafted here)

A new mission should be evaluated against the following test before
being approved. This replaces the implicit "next number in the
roadmap" rule.

### 7.1 Frontier economic test (FET)

A mission passes the FET iff it satisfies **all** of the following:

1. **Reasoning-space / compute-space ratio is high.** The mission is
   more constrained by analytical reasoning and adversarial review than
   by CPU-hours. A 30-minute problem requiring 30 hours of code is
   not Frontier; a 30-hour problem requiring 30 minutes of compute is.
2. **Independent validation is possible.** The result can be checked
   against external reality (published data, byte-pinned ephemeris,
   closed-form theory) OR against an internal oracle of equivalent
   strength.
3. **Output becomes durable scientific knowledge.** The result will be
   cited, reused, or referenced in future missions, not just filed.
4. **Hypothesis-distinguishing.** The mission distinguishes between
   competing hypotheses rather than merely producing numbers that
   confirm an existing hypothesis.
5. **Adversarial-survivable.** The result is robust against an
   8-track audit. The mission declares its pre-registered decision
   rules; an audit could falsify the conclusion if the rules were
   not satisfied.
6. **Capability-advancing.** The mission either unlocks a new
   capability or refines an existing one. Pure demonstration without
   reuse is penalized.
7. **Deterministic on modest resources.** The mission runs to
   completion in ≤ 10 hr single-core on commodity hardware and
   ≤ 1 GB RAM. Anything larger must justify remote/Colab use
   explicitly.

### 7.2 Goodhart avoidance

The FET is **NOT** a numerical score. It is a checklist of gates. A
mission that satisfies 6 of 7 is REJECTED if it fails gate 1. This
prevents optimizing for high scores on the wrong dimensions.

A mission is NOT evaluated by its expected information gain alone,
because that metric is hard to estimate honestly and is easily
gamed. The lab should rely on the audit culture rather than
self-reported metrics.

### 7.3 Anti-bureaucracy rules

- The FET is a checklist, not a scorecard. There is no "FET pass
  percentage"; a mission either passes all 7 gates or it doesn't.
- The FET is documented in this single file; no other file is
  permitted to redefine it.
- The FET does not require a written mission proposal. The lead
  agent may decide a mission passes the FET based on a short
  argument; the audit (if requested by a second agent) is the
  check, not a separate FET review.

## 8. Delegation and parallelism design (RECOMMENDATION)

### 8.1 The right unit of parallelism

The 12-track exploration probe that ran at the start of this session
**FAILED** — all 12 sub-agents were cancelled before producing useful
output. This is informative: parallel delegation is not free, and
12 sub-agents on loosely-related topics creates coordination overhead
that exceeds the benefit. The lab should adopt explicit parallelism
rules.

### 8.2 Delegation sizing rules

| Mission phase | Recommended delegation size | Rationale |
|--------------|-----------------------------|-----------|
| **Pre-design exploration** (1–2 parallel tracks) | 2–4 | Distinct epistemic tracks (e.g., literature scan + numerical design + adversarial pre-mortem); leads to genuine independence. |
| **Implementation** | 1 (the lead agent) | Serial; sub-agent parallelism introduces merge conflicts in shared files. |
| **Test design** | 1–2 | The lead agent designs; a second agent reviews with fresh eyes. |
| **8-track audit** | **8** (this number is justified) | Each track has a distinct epistemic role (derivation, bias theory, numerical experiment, implementation audit, independent estimator, arc design, hostile review, compute feasibility). 8 is the minimum that gives meaningful track independence; 4 is too few (tracks begin to overlap); 12 is too many (coordination overhead dominates). |
| **Cross-domain discovery** | 4–6 (only when looking outside orbital mechanics) | Distinct domains (numerical methods frontier, scientific methodology frontier, adjacent computational physics, formal verification opportunities, mission-design/optimization, cross-domain opportunities). 6 is the sweet spot; >8 produces overlapping outputs. |

### 8.3 Synchronous vs asynchronous

- **Synchronous (lead waits for all)** is correct for the 8-track
  audit pattern. The integration step requires all tracks.
- **Asynchronous (lead integrates incrementally)** is correct for
  pre-design exploration. A short feedback loop between leads and
  tracks accelerates discovery.
- **Fire-and-forget** is NEVER correct. Every sub-agent must produce
  a structured deliverable that the lead agent can integrate. Sub-agents
  that produce no output (as in the failed 12-track probe at the
  start of this session) are wasted resources.

### 8.4 Anti-Goodhart rules

- No mission should be approved solely on the basis of "it would
  require many parallel agents to do well". A mission's value is in
  its epistemic content, not in its parallelism requirement.
- Sub-agents MUST have distinct epistemic roles; they MAY NOT be
  asked to "independently verify" the same hypothesis in different
  ways (this is repetition, not independence).
- The lead agent integrates. Sub-agents do not edit the same file
  in parallel.

## 9. Autonomy model and hard gates (RECOMMENDATION)

### 9.1 Decision hierarchy

The autonomous agent MAY make the following decisions without
human approval:

- Choose the next mission within an already-approved portfolio.
- Spawn parallel sub-agents for an audit or pre-design exploration.
- Acquire public data (byte-pinned, with `-text` gitattributes).
- Modify `src/lab_utils/` when a third consumer justifies graduation
  (and the change is equivalence-pinned against the existing
  donors).
- Modify `localdocs/knowledge/` notes.
- Commit with the canonical Git identity (per AGENTS.md
  Remote-State Safety).
- Update `AGENTS.md`, `localdocs/charter.md`, `localdocs/roadmap.md`
  with corrections that are clearly necessary to prevent unsafe
  autonomous behavior.

The autonomous agent MUST STOP and report to the human when:

- A `git push` is needed and the live remote tip differs from the
  expected tip (per AGENTS.md Remote-State Safety).
- A site presents an antibot / CAPTCHA / login / verification gate
  (per AGENTS.md Responsible Web Access).
- A test fails after remediation AND the failure cannot be
  attributed to a known-bug-with-deprecation-warning.
- A scientific claim contradicts a published peer-reviewed result
  by > 3σ in a regime that has been previously validated.
- An experiment requires resources exceeding the local envelope
  (≤ 10 hr single-core, ≤ 1 GB RAM) and remote/Colab is not
  available.
- A destructive Git operation (history rewrite, force push) is
  needed for any reason other than remediation of a recently
  committed error.

### 9.2 Distinguishing uncertainty from blocker

Normal scientific uncertainty (a test that exposes a small
discrepancy, a result that requires more data to confirm) is NOT a
stop condition. It is a signal to spawn more research, not to ask
the human.

A genuine blocker is one of: missing credential, external access
restriction, destructive Git conflict, unsafe security boundary, or
a scientifically impossible task (e.g., resolving a measurement that
does not exist).

### 9.3 Minimum interaction principle

The autonomous agent should aim for one human interaction per
session: a final report. Internal scientific uncertainties are
handled by spawning audits or follow-up missions, not by asking the
human.

## 10. Resource routing (RECOMMENDATION)

### 10.1 R: drive usage policy

The R: drive has 2.27 GB free / 2.02 GB used. It is **valuable
disposable workspace**. The current usage pattern (commit only
compact summaries, use R: for large sweeps) is correct in
principle but undocumented as a policy.

The policy SHOULD be formalized:

| Resource class | Default location | Commit to repo? | R: as alternative? |
|----------------|-----------------|-----------------|--------------------|
| Source code | repo (`research/`, `src/`) | yes | no |
| Tests | repo (`research/.../tests/`) | yes | no |
| Knowledge notes | repo (`localdocs/knowledge/`) | yes | no |
| Experiment cards | repo (`research/.../README.md`) | yes | no |
| Results JSON (small) | repo (`research/.../results/results.json`) | yes | no |
| Figures | repo (`research/.../results/figures/`) | yes | no |
| Byte-pinned reference data | repo (`research/.../reference/`, `-text`) | yes (small) | yes (large multi-year) |
| Large sweeps (>10 MB raw CSV) | R: scratch by default | commit only summary | required |
| Intermediate numerics (RK4 trajectories > 100 MB) | R: scratch | commit only ascending-node crossings | required |
| Profiles, logs, caches | R: scratch | no | required |
| Virtual environments | R: scratch | no | required |
| Temporary clones / worktrees | R: scratch | no | required |

### 10.2 Local vs remote compute routing

The default is **local**. Remote (Colab) is invoked only when:

- The mission would consume > 1 GB of C: drive SSD writes, OR
- The mission would take > 10 hr single-core on commodity hardware, OR
- The mission needs a hardware configuration (GPU, large RAM) that
  local does not have.

The default is **commodity computation + abundant inference**, not
"high-performance compute". The laboratory has been successful at
sub-1-hour single-core computations; nothing in the current portfolio
has required remote resources.

The 5-yr DE441 acquisition for the proposed Lunisolar Mission Closure
(§13) is ~3.5 hr single-core — well within the local envelope.

### 10.3 Deterministic across environments

All experiments MUST remain runnable on commodity Windows
(PowerShell) + Python + uv + numpy + matplotlib. This is documented
in AGENTS.md §Environment; it is also the implicit constraint
behind every experiment's design.

External environment dependencies (specific NumPy BLAS, BLAS thread
count, GPU presence) MUST NOT be relied on for numerical results.
The lab's RK4 self-convergence tests + force-level identity tests
catch environment drift at the experiment level.

## 11. Recommended operating model after Experiment 020 (RECOMMENDATION)

### 11.1 Replace "Experiment NNN" with a mission queue

The current `Experiment NNN` model is becoming a bottleneck because:

1. **Experiment 020 is already a research program** (8 tracks, 5
   figures, multiple estimators, multi-phase ensemble, three
   remediation rounds). Calling it "Experiment 020" understates
   its scope.
2. **Remediation is a first-class activity** that the numbered
   model treats as an afterthought. Exp 015 → 016, 016/017 → 018,
   018 → 019, 019 → 020 are all "the same scientific question,
   revisited because the prior answer was wrong". This deserves its
   own mission structure.
3. **Capability graduation is invisible** in the numbered model.
   `j2_rhs` moved into `src/lab_utils/orbits.py` after Exp 009 → 012,
   but this was a side effect, not a tracked outcome.
4. **Portfolio diversity is invisible**. Twenty orbital-mechanics
   experiments with no energy-systems / computer-architecture /
   cybersecurity work is a portfolio concentration that the
   numbered model does not surface.

The recommended model:

```
research/<domain>/missions/<topic>/
├── README.md              # mission card: question, hypothesis, FET, scope
├── plan.md                # plan-of-record: phases, sub-experiments, criteria
├── experiments/<exp>/     # one or more experiments (existing format)
├── results/               # mission-level synthesis
└── knowledge.md           # mission-level knowledge note (links to localdocs/knowledge/)
```

### 11.2 Mission types

Five first-class mission types:

1. **Validation mission** — re-runs prior knowledge against new data
   or new analytical insight (e.g., the proposed 5-yr DE441 closure).
2. **Capability mission** — builds new shared infrastructure
   (e.g., the proposed multi-year ephemeris acquisition module).
3. **Composition mission** — composes existing capabilities into a
   new mission class (e.g., the proposed repeat-ground-track
   targeting, or a NRHO + low-energy transfer study).
4. **Discrepancy mission** — investigates a published result that
   contradicts expectation (e.g., the Lunisolar chain).
5. **Frontier exploration mission** — opens a new domain
   (e.g., a first foray into battery degradation).

The mission type is declared in the mission card and informs the
audit's framing.

### 11.3 Backward compatibility

Existing experiments 001–020 are NOT renumbered. Their directory
names (`research/orbital-mechanics/experiments/keplerOrbitValidation/`)
remain unchanged. The mission model is added ABOVE the experiment
level; experiments are sub-artifacts of missions.

A retroactive "Mission 001: Orbital Mechanics Flagship Phase 2"
could be created to group Exp 002–020 under a single retrospective
mission card, but this is optional and should not block the new
model.

### 11.4 The numbering controversy

**RECOMMENDATION**: Do NOT continue "Experiment NNN" numbering.
**RECOMMENDATION**: Do NOT introduce "Mission NNN" numbering.
**RECOMMENDATION**: Use topic-named mission IDs (`mission_<topic>`).
This avoids the "we're up to Mission 47, are we done yet?" problem
and makes each mission's scope visible in its name.

If historical ordinals are needed for citation, use the canonical
form `Exp NNN / mission_<topic>`.

## 12. What should NOT change (RECOMMENDATION, derived from §3 success)

These principles have demonstrably worked and should remain stable:

1. **Deterministic execution.** Fixed seeds, fixed inputs, byte-pinned
   external data, deterministic figures, reproducible JSON payloads.
2. **Reality is the verification layer.** Byte-pinned JPL / Horizons
   data + closed-form theory + conservation laws.
3. **Hostile review.** The 8-track audit pattern has caught three
   published errors and is the lab's primary defense against
   self-deception.
4. **Independent validation.** Synthetic oracles, mutual
   cross-checks, adversarial mutants, pre-registered decision rules.
5. **Durable source of truth.** All scientifically valuable material
   under the repo root; no dependence on R: for reproducibility.
6. **Scratch separation.** R: for large / ephemeral work; repo for
   durable records.
7. **Signed Git history + live remote checks.** Per AGENTS.md
   Remote-State Safety; the 2026-08-22 incident motivated this rule.
8. **Transparent remediation.** Errors are retracted, not hidden;
   deprecation warnings preserve prior artifacts.
9. **Minimal dependency growth.** numpy + matplotlib + pytest is
   the lab's effective stack; everything else is opt-in per mission.
10. **Reuse of verified infrastructure.** `src/lab_utils/` is
    preferred over per-experiment re-implementations.
11. **Evidence-driven stopping.** A mission stops when it has answered
    its question to the pre-registered criterion, not when it has
    consumed its budget.

The recommendation is NOT "more machinery"; it is "preserve what
works, replace what is becoming a bottleneck".

## 13. Ranked shortlist of next missions (RECOMMENDATION)

These are ranked by expected epistemic value / marginal compute ratio
and capability-advancement, with the FET applied to each. They are
not executed in this probe.

### 13.1 Mission 1: Lunisolar Capability Closure + Estimation Library Graduation

**Type**: Validation + Capability hybrid.
**Question**: Does the 018 corrected secular formula match the
asymptotic limit at LEO SSO when estimated by theory-driven harmonic
regression on a 5-yr DE441 arc with the 2-window phase-locked estimator?

**FET gates**:
1. Reasoning/compute ratio: High — analytical estimation theory +
   3.5 hr compute.
2. Independent validation: Byte-pinned 5-yr DE441 + corrected
   formula + harmonic regression oracle.
3. Durable knowledge: Yes — closes the open 020 question.
4. Hypothesis-distinguishing: Yes — sign, magnitude, and bias-
   correction form are all testable.
5. Adversarial-survivable: Yes — pre-registered estimators and
   decision rules.
6. Capability-advancing: High — graduates a reusable
   multi-year-ephemeris acquisition + harmonic-regression estimator
   library to `src/lab_utils/`.
7. Deterministic on modest resources: Yes — 3.5 hr single-core, ~2 MB
   data.

**Expected outcome**: Either the corrected formula is confirmed
(within 10–20%) and the 019/020 extrapolation is finally falsified or
confirmed, OR a missing physics term is identified (likely J2 ×
Lunisolar coupling). In either case, the open question is closed.

**Why rank #1**: This is the lab's live open scientific question, and
it is the cheapest credible closure. The 5-yr arc + phase-locked
estimator was the explicit 020 Track 6 recommendation. NOT doing it
leaves the chain unresolved indefinitely.

### 13.2 Mission 2: Repeat-Ground-Track Targeting (Landsat 16-day at SSO)

**Type**: Composition mission.
**Question**: Can the lab compose the Exp 008 ground-track
machinery + Exp 012 SSO lock + Exp 015 launch-window feasibility into
a defensible end-to-end repeat-ground-track targeting analysis for
the Landsat-7 16-day / 705 km reference orbit?

**FET gates**:
1. Reasoning/compute ratio: High — closed-form repeat-cycle theory
   exists (Lansdale 1991); verification is 1–2 hr compute.
2. Independent validation: Published Landsat-7 reference orbit
   + repeat-cycle ground-truth at 16-day repeat.
3. Durable knowledge: Yes — new mission class opened.
4. Hypothesis-distinguishing: Yes — distinguishes
   "machinery-composes-correctly" from "machinery-has-subtle-bug".
5. Adversarial-survivable: Yes — pre-declared repeat-cycle accuracy
   ±1 km.
6. Capability-advancing: Medium — exercises multi-constraint
   composition on a real-world mission.
7. Deterministic on modest resources: Yes — 1–2 hr single-core.

**Expected outcome**: A verified repeat-ground-track design
methodology, with a specific Landsat-7 reproduction as the test
case. Independently audits Exp 015's composition capability.

### 13.3 Mission 3: Estimation Doctrine Graduation (library)

**Type**: Capability mission (infrastructure).
**Question**: Can the lab extract the lessons from audit-020 (exact
OLS bias formula, regime-based convergence, harmonic-regression
estimator, phase-locked windowing, frame-convention checks) into a
reusable `src/lab_utils/estimation.py` library with documented
applicability regimes?

**FET gates**:
1. Reasoning/compute ratio: High — pure design + documentation.
2. Independent validation: Each estimator validated against the
   audit-020 synthetic oracle and against a synthetic-data
   benchmark.
3. Durable knowledge: Yes — codifies the lab's estimation doctrine.
4. Hypothesis-distinguishing: N/A — capability mission.
5. Adversarial-survivable: Yes — pure code + tests; deterministic.
6. Capability-advancing: High — every future mission that estimates
   a secular rate reuses this library.
7. Deterministic on modest resources: Yes — pure Python + numpy.

**Expected outcome**: A `lab_utils/estimation.py` library with
documented estimators (direct OLS, secant, cycle-averaged, harmonic
regression, phase-locked), each with bias formulas and convergence
regimes. Future estimation work references this library rather than
re-deriving from first principles.

**Why this is rank #3**: It is the most efficient way to convert
the lunisolar investigation's hard-won lessons into durable
infrastructure.

### 13.4 Mission 4: NRHO / Low-Energy Cislunar Trajectory Study

**Type**: Frontier exploration mission (within orbital mechanics).
**Question**: Can the lab extend the Exp 011 CR3BP framework to
near-rectilinear halo orbit (NRHO) maintenance + low-energy transfer
design, with byte-pinned JPL DE441 + GRAIL-derived lunar gravity
field as the validation source?

**FET gates**:
1. Reasoning/compute ratio: High — analytical theory + small
   compute budget (days, not weeks).
2. Independent validation: GRAIL + CAPSTONE tracking data +
   ARTEMIS published reconstruction.
3. Durable knowledge: Yes — opens cislunar mission design as a
   lab capability.
4. Hypothesis-distinguishing: Yes — Lyapounov orbit family structure
   + low-energy transfer cost vs direct Hohmann.
5. Adversarial-survivable: Yes — independent of LEO machinery;
   adversarial review would re-derive the NRHO existence condition.
6. Capability-advancing: High — new regime (3-body + lunar gravity).
7. Deterministic on modest resources: Mostly yes; needs lunar
   gravity harmonics (small data, available publicly).

**Expected outcome**: NRHO existence conditions + a CAPSTONE
trajectory reproduction + a station-keeping Δv budget for an
NRHO-resident asset. Opens a new mission-design class.

**Why this is rank #4**: It is the most strategically important
new domain within orbital mechanics, but it requires more new
infrastructure than Missions 1–3.

### 13.5 Mission 5: Adversarial Open-Data Audit (cross-domain exploration)

**Type**: Frontier exploration (outside orbital mechanics).
**Question**: What is the cheapest credible scientific-computing
problem outside orbital mechanics where the lab's hostile-review
infrastructure could produce a publication-quality result?

**FET gates**:
1. Reasoning/compute ratio: TBD per candidate.
2. Independent validation: TBD per candidate.
3. Durable knowledge: TBD per candidate.
4. Hypothesis-distinguishing: TBD per candidate.
5. Adversarial-survivable: TBD per candidate.
6. Capability-advancing: TBD per candidate.
7. Deterministic on modest resources: TBD per candidate.

**Recommendation**: NOT to be approved until Missions 1–3 have
landed. Cross-domain exploration without a working estimation
library + multi-year acquisition capability is premature.

## 14. Anti-bureaucracy self-check (RECOMMENDATION)

This probe has produced the following artifacts:

- **This document** (`POST_ROADMAP_PROBE.md`).
- A modified `AGENTS.md` is NOT being recommended (the probe was
  read-only and the rules are stable).
- A modified `localdocs/roadmap.md` is NOT being recommended; the
  recommended Mission model in §11 supersedes the numbered roadmap
  implicitly and should be proposed separately if adopted.

The artifact set is intentionally minimal. The next-session
execution should produce only the mission card and one or two
experiment directories, not a new governance framework.

## 15. Limits of this probe (FACT + INFERENCE)

### 15.1 What this probe did NOT do

- Did NOT run any experiment.
- Did NOT modify any experiment artifact, scientific code, test,
  result, or knowledge note.
- Did NOT modify `AGENTS.md`, `localdocs/charter.md`, or
  `localdocs/roadmap.md`.
- Did NOT exercise the parallel-delegation pattern (the 12-track
  exploration probe failed before producing outputs; this probe
  was conducted serially by the lead agent).
- Did NOT acquire new external data.
- Did NOT commit anything to Git.

### 15.2 Open questions for human review

The following decisions are explicitly NOT made by this probe and
require human input if they are to be adopted:

1. **Adoption of the Mission model (§11)** vs continuation of the
   numbered-experiment model.
2. **Adoption of the Frontier economic test (§7)** as a binding
   gate.
3. **Adoption of the resource routing policy (§10)** as a formal
   rule (vs the current implicit policy).
4. **Whether Mission 1 (lunisolar closure) is the right next
   mission** vs continuing the chain with a longer arc (10–18 yr).
5. **Whether Missions 2–5 are the right portfolio** vs an
   alternative.

### 15.3 Known unknowns

- The exact cost of multi-year DE441 acquisition (Track 8 estimated
  ~120 s HTTP; not yet executed for 5+ years).
- The exact RAM/SSD envelope of multi-year integration (a 10-yr arc
  at dt=60s with full state logging is ~250 GB; subsampled at
  ascending-node crossings is <1 MB).
- Whether the corrected 018 formula is the right asymptotic
  prediction. This is what Mission 1 is designed to test.
- Whether the lab can sustain a 4–6 sub-agent swarm across mission
  boundaries without coordination overhead dominating.

---

## Appendix A: Open-vs-closed status of every claimed result

This is the canonical "what does the lab actually know" table,
derived from §3 + §6 + the experiment READMEs.

| Claim | Status | Evidence | Audit |
|-------|--------|----------|-------|
| Newtonian gravity → elliptical orbits | VERIFIED | Exp 002 | synthesis 001–006 |
| Kepler equation solvers agree to ≤1.8e-13 | VERIFIED | Exp 003 | — |
| Hohmann Δv closed forms vs RK4 | VERIFIED | Exp 004 | — |
| Bi-elliptic crossover R* = 15.5817 | VERIFIED | Exp 005 | — |
| Plane-change s→∞ limit = bi-parabolic | VERIFIED | Exp 006 | audit-006 |
| Gravity assist Δε_max ceiling | VERIFIED | Exp 007 | — |
| Ground-track invariants | VERIFIED | Exp 008 | — |
| J2 secular nodal/apsidal rates | VERIFIED | Exp 009 | — |
| Drag decay law | VERIFIED-WITH-LIMITATION | Exp 010 | — |
| L1–L5 equilibria + Jacobi | VERIFIED | Exp 011 | 6 mutants |
| SSO inclination lock | VERIFIED | Exp 012 | adversarial battery |
| JPL DE441 vs ISS J2-only | VERIFIED-WITH-LIMITATION | Exp 013 | 6-track panel |
| Eclipse event finder | VERIFIED | Exp 014 | dual-route agreement |
| Dawn-dusk SSO feasible launch windows | PARTIALLY-VERIFIED | Exp 015 | 8-track retraction |
| LST at ascending node ~constant (EoT envelope) | VERIFIED-WITH-LIMITATION | Exp 016 | first-principles |
| Lunisolar closed-form over-estimates by 170× | DEPRECATED | Exp 017 | audit-018 |
| Corrected doubly-averaged quadrupole formula | VERIFIED-WITH-LIMITATION (sign + order of magnitude) | Exp 018 | — |
| 1-yr linear fit bias = O(1/W²) | VERIFIED | audit-019 Track F | — |
| 019 polynomial-in-1/W extrapolation = secular limit | REFUTED | audit-020 Track 3 | — |
| Lunisolar RAAN secular limit at W → ∞ | **UNRESOLVED** | audit-020 | — |
| Harmonic regression recovers known secular | VERIFIED (synthetic); UNVERIFIED (real data) | audit-020 Track 3 | — |
| 8-track audit pattern catches published errors | VERIFIED | audits 015, 018, 019, 020 | — |

## Appendix B: References

- `AGENTS.md`, `localdocs/charter.md`, `localdocs/roadmap.md`
- All experiment READMEs in `research/orbital-mechanics/experiments/`
- All knowledge notes in `localdocs/knowledge/`
- All audit reports in `localdocs/reports/`
- Synthesis report: `localdocs/reports/orbital_mechanics_001_006_synthesis.md`
- Lab_utils: `src/lab_utils/{integrators,orbits,earth_frames,metrics,results}.py`

---

**End of probe.**