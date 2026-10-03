# The Skeptic's Stack — how to build trustworthy computational science without a cluster

**Date:** 2026-10-02. **Status:** METHODOLOGY (not a result; no scientific claim).
**Type:** Capability note. Companion to `localdocs/reports/scientific-stack-audit-2026-10-02.md`.
**Why it exists:** the Lab reached non-trivial, adversarially-audited science on **four
dependencies, 8 CPU cores, no GPU**. This note records *why that worked* and *what would
have to transfer* if the same method is applied to electricity, nuclear, or aerospace —
domains that conventionally start from a datacentre.

---

## The one-paragraph version

Real science (weather, nuclear, aerospace) runs on datacentres. But the results that
actually change a field are overwhelmingly produced by small teams running
well-understood open-source libraries against well-posed questions. The scarce resource
is **not compute — it is doubt, directed correctly**. The Lab's thesis is that you can
buy almost no compute and still do real science, *provided* the doubt is aimed at the
physics and never at the floor. This note is the specification of that stack.

---

## 1. Three admissible lists

### 1.1 An admissible library

- **Open source**, inspectable end to end. If the tool cannot be read, its error bar is fictional.
- **Reference or peer-reviewed implementation** — a published algorithm with a citable
  provenance, not a blog post or a wrapper someone wrote on a Tuesday.
- **Version-pinned**, and the pin recorded in the same commit as the result that used it.
- **Deterministic** on commodity hardware.
- **Agent-operable**: installable, importable, and usable by an autonomous agent without
  credentials, licences, or interactive setup. A tool that needs a login cannot be used
  in an unattended run, so it cannot be a floor.

Explicitly **not** sufficient on its own: a large multi-install dependency tree with no
domain-specific validation value (see `astropy`, rejected in the stack audit §6 — the
Lab works in ICRF/TDB against byte-pinned vectors and has no time-scale ambiguity to fix).

### 1.2 An admissible claim

The evidence hierarchy (`LAB_CONSTITUTION.md` §3.2) is the operative list; restated here
because it is the transferable part:

| Level | Test |
|---|---|
| E0 | tests pass on synthetic input |
| E1 | re-running produces identical bytes |
| E2 | an independent re-implementation agrees |
| E3 | two different methods agree |
| E4 | conforms to byte-pinned external data |
| **E5** | **survives an adversarial audit that was trying to kill it** |

**Only E5 is "VERIFIED."** E0–E4 are necessary and insufficient.

Two clauses that are routinely violated and are the highest-value part of the doctrine:

- **Pre-registration is protocol quality, not evidence.** Writing the decision rule before
  seeing the data reduces selection bias; it does not make a number true. These are graded
  separately and must not be conflated.
- **Negative results are first-class artifacts.** A falsified hypothesis with a recorded
  reason is a *product*, not a failure. The Lab's repository contains more value in its
  dead formulas than in most projects' best results.

### 1.3 An admissible doubt

Doubt must be aimed at the layer that can actually be wrong. The Lab's own history makes
the priority ordering empirical, not stylistic — see §3.

---

## 2. The floor is not the product

**The deliverable of an autonomous research lab is the trail of killed ideas, not the
surviving number.**

Practically, that means:

- Wrong formulas are **renamed and kept with a deprecation marker**, not deleted.
  `closed_form_lunisolar_raan_rate_rad_s` (Exp 017) still exists, with a
  `DeprecationWarning`, beside its corrected replacement.
- Superseded knowledge notes **link forward** to the current claim, never backward to
  their own superseded self.
- A retraction is a signed commit with its own rationale, exactly like a result.
- The audit trail records which claim died, on what date, and against which pre-registered
  rule.

This is not fastidiousness. A repository that deletes its mistakes is a repository whose
code is a lie about its own past, and whose next agent will confidently re-derive the same
error.

---

## 3. The empirical finding that shapes everything else

The Lab's stack audit (`scientific-stack-audit-2026-10-02.md` §5) checked whether better
libraries would have prevented the errors that actually occurred. **None of them would
have.** Every serious failure in the Lab's history was a *physics, convention, or
modelling* error:

| Error | Would a stronger numeric library catch it? |
|---|---|
| Force-mode contamination (`use_j2 = mode != 'kepler_only'`) | No — logic |
| DE441 interpolation weights (moved a coefficient 5×) | No — data pipeline |
| IAU-1976 precession matrix transpose | No — convention |
| "SNR 6.89" from a mis-specified 2nd-order polynomial | No — modelling |
| Sidereal-vs-solar rotation-rate confusion | No — concept |
| Three-compound closed-form derivation error | No — derivation |

**Consequence for planning:** library quality is *hygiene*, not *risk reduction*. Spending
effort there is legitimate and cheap, but it is not where the Lab's ability to catch
errors comes from.

That ability came from four things, all of which are cheap and all of which transfer:

1. **Sanity-check gates inside the analysis itself.** The force-mode bug was caught by an
   in-script check that printed ~1 deg/day residuals — a *wrong answer that looked wrong*.
2. **Pre-registered decision rules** written before the result existed, with numeric
   thresholds, so a later session cannot rationalize a failure.
3. **The 8-track adversarial audit**, where tracks are separated by *epistemic role*
   (derivation, bias theory, implementation audit, independent estimator, hostile review),
   not by author. Four audits have been run this way.
4. **Kill criteria.** The frozen-plane control had a pre-registered "<25% survives" rule.
   When 52% survived, the hypothesis died. It could not be argued into survival.

---

## 4. The independence trap (documented, not solved)

`LAB_CONSTITUTION.md` Appendix C **Q8** records an unresolved gap: the independence test
requires two pieces of evidence to be able to disagree without coordination, but it
**cannot detect shared lineage**.

This is not hypothetical in the Lab. A DE441 *interpolation* defect would have silently
corrupted every mission that used those snapshots simultaneously — four missions,
perfectly consistent with each other, all wrong together. The mission that found it did
so by re-deriving with a corrected loader, not by being independent.

**Rule of thumb this yields:** agreement among N estimators is only evidence if the
estimators differ in their *inputs* as well as their *methods*. Estimators sharing a
propagator, a snapshot loader, or a covariance convention measure **estimator-family
spread**, not epistemic independence. When reporting "four estimators agree within 4%",
state what they share.

**Practical mitigation available today:** require at least one discriminating track or
estimator to consume a genuinely independent input (different ephemeris, different
propagator class, or a synthetic oracle with known ground truth). This is currently
practice, not doctrine; formalizing it is an open amendment candidate.

---

## 5. What transfers, and what does not

### Transfers intact

- Evidence hierarchy E0–E5.
- Pre-registration as protocol quality.
- Negative results as products.
- Role-separated adversarial audit.
- Sanity-check gates and kill criteria.
- Byte-pinned inputs with checksums recorded beside results.
- The deprecation discipline.

### Transfers with modification

- **The estimator ladder.** For orbital mechanics the ladder was: node-crossing OLS →
  secant → node-vector → mean-element integration → angular-momentum vector. Its power
  came from including estimators that share no failure mode. Any new domain needs the
  same property, and the concrete ladder differs.
- **Independence by role.** The roles ("bias theory", "implementation audit") are
  domain-agnostic; the *contents* are not.

### Does not transfer

- **The 1% synthetic-oracle gate as a universal number.** That threshold is justified
  here because the estimator is one-dimensional and smooth. Another domain needs its own
  oracle design and its own threshold, justified the same way.
- **Byte-pinned Horizons/DE441 snapshots.** Their role is filled by whatever the
  domain's authoritative public dataset is, under the same discipline: pinned, hashed,
  committed, offline-analysable.

### Does the domain need more compute?

Per-domain, honestly:

| Domain | Extra stack likely required | Why |
|---|---|---|
| **Orbital mechanics (current)** | nothing | point-mass + J2 + third-body is a 6-D ODE; fixed-step RK4 at dt=60 s is ample |
| **Electricity / power systems** | `xarray` + netCDF/HDF5, `dask`, a solver (HiTAMAT/POWERFLOW), real network topology | starts from large gridded/time-series data, not a handful of scalars |
| **Nuclear** | lattice-QCD or transport codes, HPC | genuinely compute-bound; FET gate 7 would fail |
| **Aerospace / CFD** | HPC, mesh generators | same — the roadmap explicitly ruled aerospace out as "supercomputer-only and NOT our domain" |
| **Materials / chemistry** | `pymatgen`, `ase`, ML potentials | compute scales with structure count |

The **method** transfers to all of them. The **compute envelope** does not, and FET gate
7 (`≤10 hr single-core, ≤1 GB RAM`) would legitimately exclude several. That is the gate
doing its job.

---

## 6. How to tell whether a skeptic is aimed correctly

Five questions. If any answer is "no", the doubt is being spent on the floor instead of
the claim:

1. Could this number be wrong because of a **convention** (frame, time scale, sign,
   unit, epoch)? — the Lab lost three findings to exactly this.
2. Could it be wrong because two "independent" results **share an input**?
3. Could it be wrong because a **threshold was chosen after** seeing the data?
4. Could it be wrong because a **sane-looking bug** passed the tests (the 017 closed-form
   had passing tests against a formula that was mathematically wrong)?
5. Could it be wrong because the **estimator's error bar assumes noise that does not
   exist**? — the `a11` SNR 6.89 case.

None of those five is fixed by installing a better library. All five are fixed by method.

---

## 7. One-line statement of the transferable thesis

> The scarce resource in computational science is not compute; it is **doubt aimed at
> the right layer**. A lab with four dependencies, a byte-pinned dataset, a
> pre-registered decision rule, and a role-separated adversarial audit can out-produce a
> datacentre that has none of those — provided it is honest about which layer its doubt
> is aimed at, and willing to kill its own favourite result in writing.

---

## What to cite

- `localdocs/reports/scientific-stack-audit-2026-10-02.md` — the dependency audit, the
  three-bucket classification, the astropy decision, and the §5 finding that library
  upgrades do not prevent this repository's error classes.
- `LAB_CONSTITUTION.md` §3 (evidence doctrine), §9 (FET), §13 (mission selection),
  Appendix C Q8 (the unresolved shared-lineage gap).
- `localdocs/knowledge/j2-precession-modulated-lunisolar-term.md` — a worked example of
  a mission that killed its own mechanism claim against pre-registered thresholds.