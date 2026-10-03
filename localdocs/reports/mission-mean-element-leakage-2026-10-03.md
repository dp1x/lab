# mission_mean_element_leakage — Phase-0 verdict and window-closure audit

**Date:** 2026-10-03 · **Session type:** autonomous, multi-track adversarial audit
**Mission:** `research/orbital-mechanics/missions/mission_mean_element_leakage/` (ACTIVE)
**Authority:** `LAB_CONSTITUTION.md` §§0.1, 3, 9, 10; mission `README.md` §4 (frozen)
**HEAD at audit:** `0d8f7e4` (= live `origin/main` at audit time)

---

## 0. One-paragraph answer

The surviving **−8.564733e-4 deg/day** SSO residual is **NOT** a new dynamical
interaction, and it is **not** explained by the mean-element treatment, by node-crossing
sampling, by a step-size artifact, or by the OLS conditioning. It is **not a
window-independent quantity at all**: measured over the same propagation with the same
code, merely by sliding the analysis window, the "residual" ranges from **+1.69e-3** to
**−1.90e-3 deg/day** and changes sign. Its magnitude at the committed window
(8.53e-4) is **smaller than its own window-to-window standard deviation**
(4.10e-4 → the value sits 2.1σ out; at the half-cycle window the σ is 1.17e-3 and the
committed value is only 0.73σ). The detrended node-difference series is 99.4% explained
by a single sinusoid of ~35 deg amplitude at a ~9000-day period, i.e. **long-period
forced content sampled over ~0.89 of a lunar nodal cycle**. Over less than one cycle a
sinusoid's OLS slope is a phase-dependent quantity of order `A·ω`, which is exactly the
size of the residual being discussed.

**The residual therefore does not currently constitute evidence for any physical
effect.** It is a single-phase sample of a long-period forced term.

---

## 1. FACT — what the committed evidence proves

### 1.1 The pre-registered rule is *not* satisfied as written, and is *partly undefined*

Evaluated clause by clause against both committed arms:

| §4 clause | Requires | Verdict |
|---|---|---|
| 1a. ≥70% of `R` at **all three** inclinations, 1 yr | 1-yr 5-mode `R` | **FAIL** (30°: 9.8% smoothed / 27.5% raw) |
| 1b. 18.6-yr SSO residual < 5.72e-3 | 18.6-yr "residual" | **PASS but definition-dependent** (§4 never defines "residual") |
| H-baseline CONFIRMED | 1a ∧ 1b | **FAIL** |
| H-sampling CONFIRMED | estimators A and C | **UNDEFINED — data does not exist** |
| PARTIAL (both) | both quantified, per incl. **and altitude** | **UNDEFINED — unreachable** |
| H-phys required | `R` at the arc | **UNDEFINED** (18.6-yr ran 2 modes; `R_measured` correctly omitted, `mode_set_complete: false`) |
| FALSIFIED (all artifact arms) | ≥50% at SSO **and** 30° | **UNDEFINED** |
| Oracle gate <1% | synthetic oracle | **NOT SATISFIED at the time Phase 0 was written** |

**The oracle gate is binding and gates every real-data claim.** It was scheduled for
Phase 1, which meant Phase 0 could not be reported as evidence *by construction*. That is
now discharged: the estimator ladder passed its oracle this session (§2.4), and the
oracle gate is satisfied for estimators EST-1/EST-2/EST-1b.

### 1.2 Arithmetic identities hold; two reported columns are definition-inconsistent

All identities reproduce exactly from the committed scalars. Two defects in reporting:

1. **"% removed" hides sign-reversing overshoots.** At 90° the correction is **108.45%**
   of the signal (it *reverses* it); the handoff's "91.6%" silently discards that. At
   1-yr 90° it is 128.2%. The column should report the signed ratio.
2. **The 1-yr and 18.6-yr artifacts came from different code** (`provenance.code`
   `49d20cf9…` vs `646d68f7…`), and the current code does not reproduce the 1-yr artifact's
   algebra (identity mismatch up to 16.6% relative at 30°). The 18.6-yr artifact *does*
   satisfy the current algebra.

### 1.3 A contaminating nuisance term is comparable to the residual at SSO

`D = ⟨Ω̇_J2⟩_full − ⟨Ω̇_J2⟩_j2` differences two *mean analytic* rates, but the chain
subtracts `j2`'s *OLS slope*. The gap is a systematic offset:

| i | contaminant `⟨Ω̇_J2⟩_j2 − rate_j2` | committed residual | ratio |
|---|---|---|---|
| 97.7876° | **+1.059212e-3** | −8.564733e-4 | **1.24×** |
| 90.0° | 1.6e-18 | −4.243816e-4 | ~0 (cos i = 0) |
| 30.0° | **−1.237715e-2** | −1.951493e-4 | **63×** |

**At 30° the estimand sits below its own nuisance floor by a factor of 63**, so the 30°
column carries no usable information about H-baseline at either window.

---

## 2. What was falsified this session

### 2.1 FALSIFIED — "the semi-major axis decays 165 km over 18.6 years"

The committed `a_mean ≈ 6891` vs `a_first = 6978.137` is a **step-size artifact**, not
physics. A conservative model cannot lose `a`.

- `elements()` is correct vis-viva (bit-identical to an independent `−μ/2E`).
- The drift is **linear in time and 5th order in dt** (measured p = 5.001 / 5.000 / 5.000
  at dt = 60/30/15 s), the signature of RK4 truncation whose per-orbit O(dt⁵) error
  accumulates as O(t).
- Extrapolated `kepler_only` rate −0.021055 km/day × 18.6 yr = **−143.0 km**, within
  1.14× of the committed −164.8 km.
- A **second, expected** term is real: the J2 short-period wobble of `a` (±9.2 km at
  SSO, ±2.1 km at 30°) sampled only at crossings aliases to +9.14 km / +2.29 km.
- **Fix:** dt ≤ 30 s bounds the arc drift to ≈ −4.5 km; dt = 15 s gives ≈ −0.14 km.

### 2.2 FALSIFIED — "the `a` drift contaminates the residual"

It is **common-mode and cancels**. `D`'s pure-`a` channel is
−7.50e-6 deg/day at SSO (**0.034%** of `D`) and +1.58e-8 at 90° (**0.000%**) where `a`
still dropped 86 km but `cos 90° = 0`. Confirmed three ways: closed-form decomposition;
dt-ladder (residual shifts 0.5% from dt 60→30); and my own re-propagation, where the
crossings-sampled `⟨a⟩` deficit is −13.74 km at dt=60 but −2.58 km at dt=30 while the
residual is unchanged to 1.1%.

**It does contaminate the ABSOLUTE `Ω̇_J2` level (by 0.50%),** which matters for any
absolute-rate claim but not for the difference-based residual.

### 2.3 FALSIFIED — "the mean-element treatment is inadequate" (Track B, literature + derivation)

For a J2-perturbed near-circular orbit, the short-period terms are

```
Ω_sp = +3 J2 (R/a)² cos i · sin u
i_sp = −(3/2) J2 (R/a)² sin 2i · cos u
a_sp = +2 J2 (R/a) e · sin u
```

At the ascending node `u ≡ 0`, so **`Ω_sp` is exactly zero there** — the sampling point is
the node of the short-period node term's own argument. Further, a per-revolution harmonic
sampled at one point per orbit aliases to **exactly zero** contribution to an OLS slope.
**There is no small-divisor channel between short-period content and the nodal drift in
this sampling scheme.** Directly confirmed numerically at dt = 7.5 s: the detrended
short-period `Ω` is ≤ 3.2e-4 deg, four orders of magnitude too small to bias an
8.6e-4 deg/day rate.

Independently, **Vallado (2013) §9.6 p.654**: *"for a first-order theory, it's immaterial
whether we use mean or osculating elements on the right-hand sides … the resulting
errors will be on the order of J²."* The lab's canon is first-order, so the
mean-vs-osculating substitution error is **O(J2²) by construction** — below the
truncation level of the canon itself.

The boxcar's one real defect is the **Jensen gap** (averaging `a,i` then evaluating the
nonlinear `Ω̇_J2`), measured from the committed scalars at **0.17% of the residual** at
SSO. Three orders of magnitude too small.

> **Retraction.** `audit-019-track-F` §4's "~0.86 deg of short-period scatter at each
> crossing" is **not reproducible** and should be treated as retracted: it multiplies the
> short-period amplitude by `T_snap/T_orb ≈ 15`, which is invalid. The lab's earlier
> position that a full Brouwer transform is the wrong tool is **confirmed and strengthened**.
>
> **Citation correction.** "Brouwer (1958) *The artificial satellite orbits of Earth*" does
> not exist. The correct source is **Brouwer, D. (1959), "Solution of the Problem of
> Artificial Satellite Theory Without Drag," *Astron. J.* 64(1274), 378–397**, DOI
> 10.1086/107958 (verified via Crossref and Vallado's reference list).

### 2.4 EXONERATED — node-crossing sampling (independent estimator, oracle-gated)

A node-detection-free estimator (fixed-cadence `Ω(t)` from the angular-momentum vector, no
equator-crossing search) reproduces the committed node-crossing `combined` to
**5.8e-9 deg/day at 1 yr across all three inclinations** — agreement to 9 significant
figures. Oracle gate: **0.083–0.257%** recovery of known injected drift over 1-yr and 4-yr
synthetic J2 arcs; ≤ 2e-15 on pure 2-body. **Node-crossing sampling geometry is exonerated.**

### 2.5 NOT REPRODUCED — the claimed OLS conditioning artifact

A track reported that `ols_slope`'s uncentred `[1, t_s]` fit yields a spurious
**+2.17e-3 deg/day** at the 18.6-yr grid, 2.5× the residual. **This does not reproduce.**
On the real committed grid (102 940 crossings, `t₀ = 8.205e8`, cond 7.5e9) the constant-
series slope is **+4.7e-18 deg/day**, and re-scoring the actual 18.6-yr series with a
centred OLS moves the residual by **8.9e-16 deg/day** — zero. Reproduced on four grid
variants; all are at 1e-17 or below. **The committed slope computation is numerically
sound; the claimed artifact is retracted.**

### 2.6 CONFIRMED — the residual is window-dependent (the central finding)

Same propagation, same code, same data (102 940 crossings), only the analysis window
varies.

**Window ladder** (start = arc start):

| window | `combined` | `D_smoothed` | residual |
|---|---|---|---|
| 1.00 yr | +1.1653e-3 | +1.0260e-3 | **+1.3929e-4** |
| 4.65 yr | +2.3365e-3 | +1.8533e-3 | **+4.8318e-4** |
| 9.30 yr | −3.1942e-3 | −4.8876e-3 | **+1.6934e-3** |
| 12.00 yr | −9.0582e-3 | −1.0975e-2 | **+1.9164e-3** |
| 18.60 yr | −2.2863e-2 | −2.2005e-2 | **−8.5792e-4** |

**Offset ladders** (sliding start):

| ladder | min | max | mean ± std | peak-to-peak | sign changes |
|---|---|---|---|---|---|
| 4.65 yr (28 windows) | −6.574e-4 | +5.615e-4 | +1.91e-5 ± 4.10e-4 | 1.219e-3 | 1 |
| 9.30 yr (19 windows) | −1.899e-3 | +1.693e-3 | −9.19e-5 ± 1.175e-3 | 3.593e-3 | 1 |

**The residual changes sign and its peak-to-peak swing is 1.4× (4.65-yr) to 4.2×
(9.30-yr) its own central value.** A secular rate cannot do this.

**Spectral identification:** a single sinusoid at a fitted ~9175-day period with ~34.8 deg
amplitude reduces the residual RMS of a linear fit to the node-difference series from
18.70 deg to **0.106 deg (99.4%)**. Every other period is ≤ 0.26 deg (annual 0.256,
semiannual 0.093, lunar monthly ≤ 0.024). The arc spans 6793.6 d = **0.894 lunar nodal
cycles**. Fitting the linear and 6798-day terms jointly is **degenerate** over 0.894
cycles (the fitted period ran to the search bound), which is itself the diagnosis: over
this arc the secular term and the lunar-nodal term are not separable, and the OLS slope
of one is absorbing the other.

For a sinusoid `A·sin(ωt+φ)` fitted over a window shorter than a cycle, the OLS slope is
`≈ A·ω·g(phase)/2π` — a **constant** bias independent of window length (this is
audit-020 Track 3's slow-harmonic regime). Here `A·ω/2π ≈ 3.3e-3 deg/day`, i.e. **4× the
residual**. The residual is not separable from one-cycle windowing.

---

## 3. What remains unresolved

1. **The identity of the ~9000-day period.** The fit prefers 9175 d (+35% vs the 6798.4-d
   lunar nodal regression) and is degenerate with the linear term. Separating them needs
   **more than one lunar nodal cycle** (> 37 yr of DE441, or two phase-locked windows),
   which exceeds the mission's ≤4 hr budget and the byte-pinned snapshot's span.
2. **Whether any secular term survives at all.** Not measured. Requires a ≥2-cycle arc.
3. **The 30° column.** Uninformative (§1.3) until the nuisance term is removed by
   regressing the correction against the analytic single-baseline form
   `full.rate − Ω̇_J2(⟨a⟩,⟨i⟩)` instead of differencing two baselines.
4. **Absolute rates.** Contaminated at the 0.50% level by the dt-60 `a` bias (§2.2).
5. **Track F's exact-`D` recomputation** moves the SSO residual by 46.6%
   (−8.56e-4 → −1.26e-3), so even the central value is estimator-dependent at that level.

---

## 4. Evidence tier and protocol tag

Per `LAB_CONSTITUTION.md` §3.2 / §3.4 (note §3.2 removed tier E6):

| Claim | Tier | Tag |
|---|---|---|
| Committed rates reproduce an independent propagation (0.085–0.088% cross-mission) | E3 | P3 |
| Node-crossing sampling is not a bias source | **E3** (cross-method, oracle-gated) | P1 |
| Mean-element/Jensen treatment cannot explain the residual | E3 (cross-method + literature) | P1 |
| `a` deficit is RK4 step-size artifact (p=5 ladder + conservation) | **E4** (conservation law) | P1 |
| Residual is window-dependent (model-free ladder) | **E3** | P0 |
| Exact attribution of the ~9000-day period | — | P0, UNRESOLVED |

**Nothing here is VERIFIED (E5).** The strongest claims are E3–E4.

---

## 5. Governance items

- `LAB_CONSTITUTION.md` §11.1 item 5 and §10.2 still reference "E0–E6"; §3.2 removed
  E6. **§12.2 still names `mission_lunisolar_closure` as the active candidate**; that
  mission completed 2026-09-03. Both are **human-gated** (§13.1 forbids autonomous
  amendment of the evidence hierarchy and the hard-stop list).
- `mission_mean_element_leakage/README.md` §4 is **frozen and was not edited**. Its clause
  1a remains unsatisfied; no threshold was retuned. If a dated amendment is warranted it
  must *disambiguate*, not relax — and should wait for §3 item 3.

---

## 6. Provenance remediation performed this session

1. **`*.py`/`*.json`/`*.md`/`*.toml`/`*.txt` pinned to `eol=lf`** in the root
   `.gitattributes`. Verified: does not change existing blobs; makes worktree == blob.
2. **Fixed an inert `.gitattributes`.** `mission_lunisolar_closure/.gitattributes` used
   repo-root-relative patterns; a nested `.gitattributes` resolves patterns **relative to
   its own directory**, so all four `-text` rules matched nothing. `git check-attr text`
   reported `unspecified` for all 12 DE441 files — **2.7 MB of byte-pinned data was
   silently unprotected**. Proven in an isolated scratch repo (root-relative → `unspecified`;
   dir-relative → `unset`), fixed at both root and nested levels, now `unset`.
   `MANIFEST.json` files for the three experiments were additionally unprotected and are
   now pinned `-text`.
3. **New regression test** `src/lab_utils/tests/test_reference_data_provenance.py`
   (7 tests) — it caught a real regression during authoring (my `*.json text eol=lf`
   rule overrode the `MANIFEST.json` protections, since gitattributes is last-match-wins).
4. **New fingerprint scheme** `lf-normalized-v1`: hashes `bytes.replace(b"\r\n", b"\n")`
   so new artifacts are clone-verifiable. **Pre-2026-10-03 pins are raw worktree-byte
   hashes and remain valid under either convention** (line endings are semantically inert
   here: `load_snapshot` reads in text mode with universal newlines). **No committed
   artifact was regenerated** — that would rewrite pre-registered evidence to satisfy a
   cosmetic check.
5. **Fixed `rolling_mean`** (see §7).

---

## 7. Code changes made

| File | Change | Justification |
|---|---|---|
| `.gitattributes` | `eol=lf` for text types; `-text` for all pinned snapshots | §6.1–6.2 |
| `mission_lunisolar_closure/.gitattributes` | dir-relative patterns + explanation | §6.2 |
| `.gitignore` | ignore `.kilo/` agent scratch | tooling hygiene |
| `phase0_baseline_breathing.py` | `lf-normalized-v1` fingerprints + `code_hash_scheme` key | §6.4 |
| `phase0_baseline_breathing.py` | `rolling_mean` trailing pad `win-1-pad` → `pad` | the window was effectively 58 d, not 30 d, and element 0 was inflated ~5 km |
| `src/lab_utils/tests/test_reference_data_provenance.py` | new (7 tests) | §6.3 |

**Not changed:** `dt_s = 60.0` default. The dt defect (§2.1) is real and the fix is dt ≤ 30 s
(15 s preferred), but changing the default would silently alter the meaning of every
committed Phase-0 artifact. It is left as a **pre-registered change for the next mission**,
to be made together with a documented re-run.

**Not changed:** `mission_mean_element_leakage/phase0_baseline_breathing.py` still imports
the predecessor mission's `mission_experiment.py`. Rewiring to `lab_utils.ephemeris` is
low-risk and equivalence-pinned, but it changes the provenance fingerprint of the code
that produced the committed artifacts — deferred to the same commit as the dt change, so
the artifacts are regenerated once, for both reasons, under a stated rule.

---

## 8. Reproduction

All figures in this report are reproducible from `.kilo/` scratch (ephemeral, not
committed). The load-bearing ones, in order of decisiveness:

1. **Window/offset ladders** — 18.6-yr propagation, 102 940 crossings, `dt = 60 s`,
   sliding analysis windows. Model-free; no fit assumptions.
2. **Spectral fit** — `dΩ = Ω_full − Ω_j2`, detrended, Lomb-Scargle + period scan.
3. **`a` conservation ladder** — `kepler_only` (where `a` is exactly conserved) at
   dt = 60/30/15/7.5 s; measured p = 5.001/5.000/5.000.
4. **Crossings-sampled vs continuous-mean `a`** at dt = 60 → 7.5 s (200 orbits).
5. **Oracle-gated node-free estimator** — reproduces committed `combined` to 5.8e-9.

---

## 9. Recommended next mission (selected by evidence, not in advance)

The evidence points to one thing: **the arc is too short to contain the signal.** Every
independent line converges on that. The single highest-information question is therefore
**whether a secular lunisolar nodal rate exists at all at LEO**, and the experiment that
answers it is a **multi-cycle, phase-decorrelated determination of the long-period
nodal term**, requiring:

- a DE441 span of ≥ 2 lunar nodal cycles (> 37 yr) or an explicitly stated
  two-phase-locked window design;
- an estimator that fits the secular term **jointly with** the long-period forced terms
  rather than letting the OLS slope absorb them (audit-020 Track 3 already showed the
  polynomial-in-1/W extrapolation has no theoretical basis and must not be reused);
- the contaminant `⟨Ω̇_J2⟩_run − rate_run` reported as an uncertainty budget on every
  slope, and the single-analytic-baseline form of the correction (§3 item 3);
- `dt ≤ 15 s`, so the absolute rate is meaningful;
- the estimator ladder's oracle gate re-run and recorded **before** the real arc.

This is a **discrepancy mission** (constitution §2.5 type 4), and it supersedes any further
J2×lunisolar-mechanism work — because the mechanism question is currently **unfalsifiable**
on an arc that cannot separate the terms. Per the constitution's own rule
(`§12.1` step 6), if no candidate passes the Frontier Economic Test the Lab stops; this
candidate is recorded as *the evidence-indicated* next mission, and selecting it formally
remains a human decision per §4.3.