# Experiment Card: Orbit Decay — Atmospheric Drag (Experiment 010)

> Status: complete | Date: 2026-08-22 | Domain: orbital-mechanics |
> Experiment dir: `research/orbital-mechanics/experiments/orbitDecay`

## Research Question

Does deterministic Cowell propagation with atmospheric drag rediscover pointwise
orbital-energy dissipation, the circular-orbit decay law, structural drag scaling
laws, eccentric-orbit decay behavior, and deterministic re-entry timing — while
separating numerical integration error from atmospheric-model uncertainty?

## Background Theory

### Frames and units

ECI frame as in Experiments 006–009 (km, km/s, s). Altitude is spherical-Earth
geocentric, `h = |r| − R_E`. Density and ballistic-coefficient quantities are SI
(kg/m³, kg/m²) inside the drag term; a single explicit factor
(`DRAG_SI_TO_KKM = 1e3`, derived in-source and unit-fired against a pure-SI hand
calculation) converts to km/s². The frozen ballistic-coefficient convention is

```
beta = m / (C_D A)  [kg/m^2],   kappa = C_D A / m = 1/beta,
a_drag = -(1/2) * kappa * rho(h) * |v_rel| * v_rel ,   v_rel = v - omega_atm x r .
```

(The handoff's preliminary contract wrote `B = C_D A/m` while sweeping magnitudes
that only make sense for `m/(C_D A)`; the kg/m² form is the entry-community
standard and is adopted here. NORAD's TLE `B*` is the inverse-flavored term.)

### Drag dynamics and dissipation accounting

Drag is the lab's first non-conservative force: energy conservation is replaced by
energy MONOTONICITY plus pointwise dissipation accounting. The exact identities are

```
d(eps_total)/dt = a_drag . v = -(1/2) * kappa * rho(h) * |v_rel| (v_rel . v)  <= 0 ,
eps_total = v^2/2 - mu/r (+ static J2 potential when J2 is active).
```

The popular `-kappa rho |v|^3 / 2` form is the inertial-atmosphere special case
(`v_rel` parallel `v`); using it under a co-rotating atmosphere is itself the
missing-rotation bug. In LEO `omega h_z/v^2 <= 0.0645`, so `v_rel.v > 0` always:
dissipation is strictly negative whenever `rho > 0`.

Circular decay law (`e=0`, inertial atmosphere, exponential layer):
`da/dt = -kappa rho(h(a)) sqrt(mu a)`; exact closed-form lifetime over one layer,

```
t(a0 -> af) = sqrt(pi H)/K [erfi(sqrt(a0/H)) - erfi(sqrt(af/H))] ,
K = kappa rho0 sqrt(mu) exp((R_E + h0)/H).
```

(The research round's contrary gamma-function claim was rejected after
re-derivation: the integrand is `a^(-1/2) e^(+a/H)` because `1/rho` grows with
altitude.) Constant-density short-window form: `sqrt(a(t)) = sqrt(a0) - kappa rho
sqrt(mu) t / 2`.

### Decay-law oracles and their limits

Five independent oracle paths (frozen contract §4): O1 Gauss–Legendre quadrature of
the separable lifetime integral with panels aligned to density-layer joints;
O2 the erfi closed form (single layers with h0 ≥ 120 km — lower layers overflow
`exp((R_E+h0)/H)` in float64), evaluated with mpmath; O3 the constant-density form;
O4 direct quadrature of the osculating identities `da/dt = (2a²/μ)(a_drag·v)` and
`d(e²)/dt = 2h²(dε/dt)/μ² + 4εh(dh/dt)/μ³·h` along one conic revolution (independent
of Cowell integration, valid at all eccentricities); O5 King-Hele first-order
per-revolution deltas with modified Bessel functions, gated to `e ≤ 0.1`.

## References

- Vallado, *Fundamentals of Astrodynamics and Applications* — exponential-atmosphere
  data file `ATMOSEXP.DAT`, official CelesTrak repository (retrieved 2026-08-22;
  verified byte-for-byte). Declared atmosphere, used verbatim and untuned.
- U.S. Standard Atmosphere 1976, NASA-TM-X-74335 (NTRS 19770009539); numeric tables
  via PDAS (`pdas.com/bigtables.html`). Used as plausibility spots only.
- Ray & Scheeres, "King-Hele orbit theory for periodic orbit and attitude
  variations", MNRAS 501, 1168 (2021), arXiv:2008.10644 — King-Hele series forms.
- Peet, "Lecture 12: Orbital Perturbations", ASU MAE462 — Gauss forms, drag work,
  `B = m/(C_D A)` convention.
- *Project Space Track Report #3* (NORAD, 1980) — `B*` definition.
- Ann. Geophysicae 39, 397–410 (2021) — canonical `C_D = 2.2` for LEO spheres
  (quoting Cooke 1965).
- ESA, "ISS reboost" (2016) — binding quiet-time decay band (~2 km/month at ~400 km).
- Oliveira, Zesta & Garcia-Sage, arXiv:2505.13752 (2025) — Starlink storm-time
  decay; CONTEXT ONLY (geomagnetic augmentation outside model scope).
- Cook, King-Hele & Walker, Phil. Trans./Proc. R. Soc. A (1960–1963) — primary
  King-Hele sources sit behind a bot challenge; NOT accessed per lab rules, cited
  here as bibliographic context only (formulas taken from the open sources above).

## Assumptions

- Declared atmosphere = full Vallado exponential stack (27 layers, 0.01–900 km),
  half-open `[h_i, h_{i+1})` selection, top row extended above 900 km. Known bias vs
  US76 spot values: +26 % @300 km, +33 % @400 km, +34 % @500 km — documented, never
  tuned away.
- Spherical non-rotating Earth unless `omega_atm` is set; isotropic density depends
  on `|r|` only.
- `omega_atm = 0` for every closed-form/scaling anchor (declared inertial
  atmosphere); `omega_atm = omega_E` only in the rotation-asymmetry cases.
- Canonical `C_D = 2.2`; `beta` sweep {50, 100, 200, 400} kg/m².
- Headline decay numbers run with J2 OFF (no secular J2-decay coupling under a
  spherical altitude-only atmosphere); the J2 interaction is validated separately.

## Methodology

1. Clone Exp 009's fixed-step RK4 loop verbatim (`_rk4_core`); add a gated drag
   branch (`propagate_3d_rk4_drag`) following the copy-with-one-branch doctrine:
   disabled branches are skipped entirely, so `(beta=0, j2=0)` is bit-exact with
   Exp 006 and `(beta=0, j2=J2)` bit-exact with Exp 009 (`np.array_equal`).
2. Resolution rule inherited from Exp 009: `steps_per_orbit(e) =
   max(512, ceil(720/(1-e)^{3/2}))`; steep-gradient case additionally gated on an
   empirical dt-halving invariance (<0.1%) before any number from it is recorded.
3. Oracles O1–O5 as above; same-law-only equality assertions, published numbers used
   exclusively as decade-wide bands.
4. Re-entry events: bracket crossing on the fixed grid, then reintegrate ONLY the
   bracketed interval at dt/2^j (j ≤ 6); thresholds 120 km and 100 km both reported.
5. Batteries: convergence (time-to-fall observable), plateau separation (law swap at
   fixed dt across a dt grid), structural scalings, rotation twins, benchmarks,
   pathological sentinels, adversarial mutants (sign flip, kg→km³ unit error,
   B inversion).

## Implementation

`experiment.py` (~1450 lines): declared atmosphere table + density helper; drag-gated
propagator; element-series and per-step relative dissipation-residual machinery;
oracles O1–O5; window/reentry drivers; convergence/order/plateau/scaling/rotation/
benchmark/eccentric/lifetime batteries; pathological + mutant harnesses; deterministic
figure generation. No new dependencies (numpy + matplotlib + mpmath already pinned).

## Validation Method

39 focused pytest tests (`tests/test_orbit_decay.py`, banners L1–L7 mirroring the
pre-registered failure catalog). Theory constants are duplicated inline in the tests
and oracles reimplemented independently where shared algebra would fake a pass.
Full repository suite green before commit.

## Results

All numbers from `results/results.json` (deterministic regeneration; figures
regenerated from recorded data only, dpi=150).

**Decay law (pillar 2).** leoRef (420 km circular, beta=100, 500 revs, single
density layer): numerical a(t) matches the quadrature oracle to max 3.6 m / RMS
2.0 m over a 3.98 km decay (0.09%); window-time agreement −8.7e-4 rel. Starlink-like
(550 km) and SSO-like (600 km) 100-rev windows agree with the quadrature oracle to
−8.3e-4 and −2.7e-3 rel respectively. erfi closed form vs quadrature: −8.2e-15 rel.
Quadrature node-doubling self-convergence <3e-15.

**Pointwise dissipation (pillar 1).** Relative per-step identity residual (finite-
difference d eps/dt vs inline power): median 1.07e-3 at default resolution, dropping
to 6.7e-5 at 4x resolution (order-4 shrink, factor 15.9); strictly dissipative
trajectories (zero monotonicity violations) with J2 off AND total-energy monotonicity
with J2 on. All aggregates recorded in `results.json`.

**Structural laws (pillar 3).** beta ratios {50/100, 200/100} within 1.7e-2 of exact
at the frozen 250-rev windows (finite-window drift nonlinearity O(kappa*Delta_a/H),
documented); rho0 uniform scalings identical to beta ratios (same-trajectory
consistency); acceleration-level linearity machine-exact (|ratio − 0.25| < 1e-12).
Omega-twin frame symmetry: 3.0e-14. Co-rotation asymmetry (equatorial i=0 vs 180 deg):
measured rate ratio 1.3032 vs exact theory 1.2946 (+0.66%); inclined pair
(63.6/116.4 deg): 1.1247 vs first-order 1.1215 (+0.29%).

**J2 interaction (finding).** Switching J2 on does NOT leave the decay rate
unchanged: seeding a Kepler state relaxes the mean elements by
`(2a^2/mu)<U_J2>` ≈ −5.7 km at i=51.6° (nonzero inclination-dependent mean of the
J2 potential), so decay proceeds in denser air. Settled-tail energy-bookkeeping
rates: measured ratio 1.1225 vs altitude-shift prediction 1.1229 — residual −0.04%.
Raw osculating-a slopes under J2 are invalid (offset transient + ripple aliasing);
documented as the estimator lesson of this experiment.

**Eccentric behavior.** Steep case (perigee 250 km, e=0.3, 300 revs): apogee drops
22.31 km while perigee moves 0.05 km (>400:1); e strictly decreasing (−8.2e-4);
measured per-rev delta-a vs Gauss-conic oracle O4: +0.05%; per-rev apsis-sampling
estimator agrees with O4 to 3.8e-4 (steep) / 1.5e-4 (Molniya-like); dt-halving
invariance 5e-6 (gate passed). Molniya-like (perigee 500 km, e=0.74, 100 revs):
measured vs O4 +0.006%. Small-e King-Hele check (e=0.05, perigee 350 km):
delta-a +2.27%, signed delta-e −6.51e-7/rev vs King-Hele −6.45e-7/rev (+0.9%).

**Re-entry timing.** 200 km circular start (beta=100): crosses 120 km at
2.01551 d and 100 km at 2.03016 d (both refined to µs-level bracket stability);
full decay from 300 km takes 38.76 d to the 120 km threshold. Numerical-vs-oracle
spot agreement <5e-3 rel (280 km start verified against quadrature).

**Lifetimes (oracle, beta=100, to 120 km).** 250 km: 10.56 d; 280 km: 23.73 d;
300 km: 38.76 d; full curve 250–800 km in `results.json`/figure F3.

**Convergence & separation.** Time-to-fall observable: rates [3.78, 5.00, 5.03]
across spp 64→512 vs 1024 reference — pre-asymptotic SUPERCONVERGENCE of this
observable (documented; raw Kepler position-error order of the clone stays ≥4.42,
i.e. no degradation below design order 4). Plateau separation: law swap moves the
transit time 103,000+ s (flat in dt to 0.73%) while dt refinement contributes only
2781 s → separation ratio 37.2×.

**Benchmarks (pillar 4).** BINDING: model quiet-time decay at 420 km = 3.64
(beta=100) / 1.91 (beta=190) km/month vs published ~2 km/month — PASS within the
pre-frozen decade band [0.2, 20]. CONTEXT-ONLY (recorded, not asserted): Starlink
storm-time 95–176 km/day.

**Adversarial mutants.** Sign flip → energy grows → monotonicity detector fires.
kg/m³-fed-as-kg/km³ → drag vanishes → pinned-probe detector fires (measured delta
7.8e-6 of the hand value). B inversion (kappa/beta swapped) → rate ×100 → detector
fires. All three plausible-but-wrong results would have passed a naive "it decays"
check.

**Hand-computed holdout.** Delta-a/rev at 420 km, beta=100:
`rho(420) = 3.725e-12 · exp(-20/58.515) = 2.6465956e-12 kg/m³`;
`Delta_a = 2 pi kappa rho a² = 2 pi · 0.01 · 2.6465956e-12 · (6.798137e6)² = 7.68506 m`.
Measured over a 30-rev window: 7.707 m (+0.28%, window drift as predicted).
(This replaces the handoff's unverifiable "~8.66 m" probe.)

Test count: 290 → **329 passed** (39 new). Runtime ~7 min (experiment) / ~6 min (new tests).

## Limitations

- Geodetic-vs-geocentric altitude differs by up to 21.4 km (pole) / 13.1 km (ISS
  inclination) → external density comparisons carry up to ~1.4x bias; folded into
  decade-wide bands only.
- Quiet-time exponential atmosphere; geomagnetic storms, diurnal/seasonal variation
  excluded (Starlink storm benchmarks recorded as context only, with reason).
- Lifetimes above ~350 km come from the quadrature oracle validated by numerical
  spot decays at 250–300 km; direct end-to-end propagation there exceeds the step budget.
- erfi closed form restricted to single layers h0 ≥ 120 km; multi-layer windows use
  quadrature.
- King-Hele few-% checks gated to e ≤ 0.1; large-e anchors use the Gauss-conic oracle.
- kappa ≥ ~1e6 m²/kg removes more energy per step than low-altitude orbits hold —
  outside the propagator's valid regime (documented; not run).
- `save_json_result` rounds |x| ≥ 1e-10 to 12 decimals; rates stored in m/day.

## Future Improvements

- Graduate shared orbital machinery (propagators, element conversion, constants) to
  `src/lab_utils/orbits.py` (now justified — five experiments share it).
- Layered-atmosphere uncertainty quantification: propagate the US76-vs-fit bias into
  lifetime confidence bounds instead of decade bands.
- Wind models (non-co-rotating atmosphere) and diurnal density bulges.
- Feeds roadmap 013 (JPL Horizons validation): unmodeled decay produces
  quadratically growing along-track error at LEO — now quantified.

---

### Reproducibility Notes

Run: `.venv/Scripts/python.exe research/orbital-mechanics/experiments/orbitDecay/experiment.py`
(~7 min, deterministic; results.json carries commit hash + timestamp + python version).
Tests: `.venv/Scripts/python.exe -m pytest research/orbital-mechanics/experiments/orbitDecay`.
Figures regenerate deterministically from `results.json` data (Agg, dpi=150):
`f1_leoref_decay_vs_oracle.png`, `f2_circularization.png`,
`f3_lifetime_vs_altitude.png`, `f4_convergence_dissipation.png`.
Declared-atmosphere provenance incl. retrieval date and blocked-source notes live in
`results.json` (`atmosphere`, `limitations`).
