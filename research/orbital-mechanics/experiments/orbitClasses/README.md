# Experiment Card: Orbit Classes — Constraint-Defined Families (SSO / Molniya / GTO / GEO anchor)

> Status: complete
> Date: 2026-08-23
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/orbitClasses/`

## Research Question

Can the operationally important Earth-orbit classes be **recovered and classified as
solution sets of coupled dynamical constraints** under a declared model (two-body
Kepler + first-order secular J2 + spherical Earth rotating uniformly at omega_E),
with each class-defining quantity independently reproduced by closed-form algebra,
by full-force numerical propagation, and by structural identities — rather than
asserted from mission folklore?

Class constraint equations:

| Class | Defining constraints | Nature |
|---|---|---|
| SSO | Omega_dot(a,e,i) = +360/365.2422 deg/day, retrograde branch | analytic lock on i; finite existence boundary |
| Molniya | omega_dot = 0 <=> cos^2 i = 1/5; n = 2*omega_E; high e; omega = 270 deg | analytic lock + resonance + engineering family |
| GEO | n = omega_E, e = 0, i = 0 | 1:1 spin-orbit fixed point |
| GTO | r_p = LEO injection, r_a = a_GEO | two-body boundary conditions (connector class) |

## Background Theory

Secular J2 rates (Vallado ch. 9, chapter-level; independently re-derived):

```
Omega_dot = -(3/2) n J2 (R/p)^2 cos i
omega_dot = +(3/4) n J2 (R/p)^2 (5 cos^2 i - 1)
M_dot     =   n + (3/4) n J2 (R/p)^2 sqrt(1-e^2) (3 cos^2 i - 1)
n = sqrt(mu/a^3), p = a(1-e^2)
```

SSO: setting Omega_dot equal to the mean-sun apparent rate gives the closed form
`cos i_SSO = -(a/a_max(e))^(7/2)` (e = 0), retrograde branch only. The existence
limit follows from `|cos i| <= 1`: `a_max^7/2 = 1.5 J2 sqrt(mu) R^2 / (lambda(1-e^2)^2)`
— eccentricity EXTENDS the limit ((1-e^2)^(-4/7)). Beyond a_max no sun-synchronous
orbit exists at any inclination.

Molniya: apsidal freeze at `cos^2 i = 1/5` (i_crit = arccos(1/sqrt5) = 63.43494882 deg,
supplement 116.56505118 deg); semi-synchronous resonance `T = P_sidereal/2 <=> n = 2 omega_E`
gives a = 26561.762 km. Dwell fraction near apogee is exactly closed-form:
for window +/-Delta about apogee, `f = (pi - E_1 + e sin E_1)/pi` with E_1 = E(nu = pi - Delta)
— slow motion at apogee concentrates >92% of the period inside ±90 deg at e = 0.74.

GEO: `a_GEO = (mu/omega_E^2)^(1/3) = 42164.169 km`. Individual Omega_dot/omega_dot are
NONZERO there (-0.013414/+0.026828 deg/day, tied by omega_dot = -2 Omega_dot at i = 0)
but act only on unobservable degenerate elements: in this model the fixed point is exact.
The apparent-longitude stationarity residual `M_dot + Omega_dot + omega_dot - omega_E`
at the Keplerian radius is +0.02683 deg/day — recorded as a NEGATIVE CONTROL against
"everything vanishes at GEO" mutants.

GTO: pure vis-viva budgets between LEO injection and the synchronous radius;
transfer time = half the Kepler period of the transfer ellipse.

## References

- D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed.,
  Microcosm, 2013 — Ch.9 general perturbations (secular rates incl. M_dot);
  Ch.3 time frames/constants (omega_E).
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier,
  2021 — Ch.10 perturbations (SSO/Molniya design conditions).
- R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*, Dover,
  1971 — Ch.9 perturbations (critical inclination).
- NIMA, WGS-84, TR8350.2 — R_E = 6378.137 km; J2 = sqrt(5)|C20_bar|
  = 1.082629821e-3; omega_E = 7.2921159e-5 rad/s.
- IAU 2015 Resolution B3 (arXiv:1510.07674) — nominal GM_E = 398600.4418 km^3/s^2.
- Real-mission anchors (Sentinel/Landsat-class SSO inclinations, classic Molniya
  geometry ~500 x ~39,900 km, GEO altitude 35786 km) are CONTEXTUAL screens only.

## Assumptions

- Two-body + first-order secular J2 + spherical-Earth uniform sidereal rotation
  (verified model boundary: tesseral/luni-solar/SRP effects are OUT of scope).
- Mean-solar-year target 365.2422 d for the sun rate (Exp 009 continuity). The
  tropical-year variant differs by 2.8e-8 deg/day and is behaviorally
  indistinguishable (documented blindness, pinned by constant literal).
- omega_E is the MASTER sidereal constant; P_sidereal always derived as 2*pi/omega_E.
- GTO injection altitude h_p = 300 km canonical (200 km kept as the Exp 004
  continuity case); coplanar budgets unless stated.
- Near-circular seeds under J2 develop induced eccentricity; omega-dot trends are
  claimed only for e >= 0.01 seeds (Exp 009 claims policy).

## Methodology

Deterministic pipeline (`experiment.py`, ~140 s single core):

1. **Closed-form layer**: domain-safe SSO solver (returns margin/status, NEVER clips
   arccos), existence-limit formula + bracketing grid, convention probes
   (sidereal/Julian/tropical years, Earth-rotation confusion), repeat-lattice solver
   m:T = k:P_sidereal, GTO budgets via borrowed Exp 004 vis-viva forms.
2. **Full-force numerical layer** (graduated `lab_utils.j2_rhs` + `rk4_propagate`):
   - SSO numeric closure at solved inclinations (Paths A/B estimators);
   - critical-inclination staged sweep (stage A 61.5..65.5 step 0.5 deg; stage B
     62.95..63.95 step 0.05 deg; dedicated offset set at icrit+{-2,-1,-0.5,+0.5,+1,+2}
     deg; plateau proof at 1024/2048/4096 spp);
   - Molniya family sweep e in {0.60..0.75} at a_semisync (spp capped at 2048,
     documented deviation);
   - period-clock runs at/near the lock (32-orbit windows, full resolution law);
   - 12.5-day apogee-drift arc; inclined-GEO nodal-period measurement (4096 spp x 12
     orbits — coarse grids quantize refined crossings to the sample step);
   - GTO RK4 flight anchor (Kepler-only + J2-on arrival).
3. **Adversarial battery** computed live and recorded: sign flips, wrong branches,
   year/day conventions, p:=a substitution, circular-v-on-ellipse, dwell-linear,
   altitude/radius bugs — each with its named catch layer.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib (Agg)
- Runtime: `uv run python experiment.py` (~140 s single core)
- Determinism: pure float64; no RNG; no wall-clock in experiment code; fixed grids;
  figures regenerate deterministically from recorded results (Agg dpi=150).
- Dependencies: numpy, matplotlib + `lab_utils` (orbits/integrators/metrics/results);
  single-hop importlib borrows: `ols_fit`, `measure_secular_rates`,
  `node_crossing_raan_rate`, `analytic_rates` (Exp 009 estimator plumbing, Exp 011
  precedent), `hohmann_dv1/dv2/transfer_elements/hohmann_transfer_time` (Exp 004
  verified closed forms). New shared machinery: `j2_rhs` graduated into
  `src/lab_utils/orbits.py` at this experiment (second consumer after Exp 009),
  pinned against the donor propagator (J2-on < 1e-12 rel; J2-off bit-exact).

## Validation Method

Five layers (`tests/test_orbit_classes.py`, 43 tests):

- **L1 closed-form identities & convention firewalls** — inline-duplicated oracles
  (anti-shared-algebra doctrine): SSO anchors to 5e-5 deg (binding solar-vs-Julian
  separation 1.67e-4 -> tolerance below half of it), wrong-branch mirror producing
  EXACTLY the negated rate, cos identity `-(a/a_max)^3.5`, monotonicity, eccentric-SSO
  both-directions trade, existence sentinel (no silent clip), year-convention
  discriminators, Earth-rotation confusion guard, critical-i identity + supplement,
  semi-sync/GEO radii, signed per-burn GTO asserts + speed-ordering chain +
  h-conservation dual form, swap-trap demonstration, bug quantifications, dwell
  values/limits/monotonicity, M_dot bracket distinctness.
- **L2 numerical recovery** — SSO closure residual <= 1% model-order band with
  positive-sign and dual-path (element-regression vs node-crossing) checks; freeze on
  the element-regression path at the lock; Kepler machinery check bit-tight (6.8e-16);
  measured Kepler excess with plateau; zero-crossing localization +/-0.15 deg with
  slope/antisymmetry/plateau proofs; dwell numeric-vs-closed; apogee event-rate
  identity; repeat lattice closures; corrected-radius disclosure; GEO negative control
  + inclined nodal shift; GTO RK4 arrival.
- **L3 convergence/invariants/determinism** — state-space order 4.09 in [3.6, 4.4]
  vs closed-form truth; rate-ladder orders all in [3.6, 5.0]; period ladder
  REPORT-ONLY near the lock (short-period jitter floor ~0.45 s documented);
  pathological grid (20 rows) all-ok with NaN sentinels; bitwise repeat-call checks.
- **L4 adversarial battery** — pre-registered survivors (omega-tests blind to branch
  flips; total-dv blind to swaps; null tests blind to p:=a; no threshold catches the
  tropical year) each covered by its named compensating discriminator; km/m and
  deg/rad firewalls; unwrap aliasing guard.
- **L5 committed-artifact integrity** — headline pins against results.json, figure
  registry, contract-block disclosure completeness.

## Results

Headline (full detail in `results/results.json`; figures carry one claim each):

| Quantity | Value | Layer |
|---|---|---|
| Sun rate (mean-solar target) | 0.985647332099 deg/day | analytic pin |
| i_SSO @ 500/600/800 km, e=0 | 97.401786 / 97.787647 / 98.603085 deg | closed form + numeric |
| SSO existence limit | a_max = 12352.505076 km (h_max = 5974.368 km) | closed form + sentinel bracket |
| SSO numeric closure residual | worst 6.14e-03 rel (model-order band <= 1%) | full-force RK4, Paths A+B |
| Critical inclination | 63.43494882 deg; supplement 116.56505118 deg | analytic lock |
| Freeze at lock (element regression) | omega_dot = -5.97e-05 deg/day (bound +/-5e-3) | 12-orbit propagation |
| Antisymmetry ratio at icrit+/-0.5 deg | 6.75e-03 (<= 2%) | staged sweep |
| d(omega_dot)/di at lock | measured -0.011412 vs theory -0.011531 deg/day/deg | sweep + theory |
| Zero crossing localized | 63.42989 deg (-0.0051 deg from exact) | stage-B bracket |
| Molniya semi-sync radius | a = 26561.762328 km; h_p/h_a = 527.92/39839.33 km @ e=0.74 | resonance lock |
| Apogee dwell (+/-90 deg, e=0.74) | f = 0.923607 (numeric err 9.8e-04) | closed form + numeric |
| J2-on Kepler excess near lock | +323.0 s/orbit (plateau ratio 1.0) | FINDING, converged |
| First-order draconitic split | +24.06 ms disclosed, NOT claimed detected | honest bound |
| Apogee event-rate identity | meas 355.46024 vs pred 355.46527 deg/day | machine-verified identity |
| Repeat-corrected radius (first-order) | 26553.420405 km (naive - 8.34 km) | design disclosure |
| GEO radius / period identity | 42164.169462 km; match 8.4e-16 rel | construction |
| GEO stationarity residual (negative control) | +0.02683 deg/day, nonzero | kills zero-mutants |
| Inclined-GEO (i=5 deg) nodal shift | pred -9507 ms, num -10199 ms (band +/-1.5 s) | first-order + measurement |
| GTO budget (300 km -> GEO) | dv1 = 2.42573, dv2 = 1.46682, total 3.89256 km/s | vis-viva + RK4 flight |
| GTO RK4 arrival | kepler 2.0e-11 rel; j2-on 4.3e-3 (disclosed) | flight anchor |
| Convergence | state order 4.09; rate orders >= 4.47; pathological 20/20 ok | numerics |

Figures (each carries one claim): `f1_sso_existence.png` (family rises to a finite
existence boundary), `f2_omegadot_vs_inclination.png` (apsidal lock made visible),
`f3_molniya_dwell.png` (dwell geometry over one orbit), `f4_gto_budget.png` (budget
decomposition), `f5_convergence.png` (order proofs), `f6_repeat_tracks.png` (apogee
event-rate identity bar comparison).

## Limitations

- First-order secular J2 only: near i_crit the FULL problem carries
  small-divisor-amplified short-period dynamics — MEASURED here as osculating-a
  excursions ~ +160 km and apsis/node event-period excess ~ +325 s/orbit beyond any
  first-order clock (energy-conserving to the integrator floor, plateaued under step
  refinement). Consequences: (a) millisecond period-split detection is infeasible by
  event timing (disclosed, not claimed); (b) event-based "periods" are window-length
  sensitive (~0.45 s floor at 8-orbit ladders); (c) the repeat-corrected radius is a
  FIRST-ORDER DESIGN DISCLOSURE, not a propagated closure.
- Real GEO fights tesseral + luni-solar + SRP drift (outside model); real SSO needs
  J3/J4/luni-solar refinements and true-Sun seasonality; real Molniya adds lunisolar
  perturbations and station-keeping deadbands.
- Mission anchors used contextually; none tuned to.
- The tropical-year constant variant (2.1e-7 deg in i_SSO) is below behavioral
  discriminability — pinned by literal instead (documented blindness).

## Future Improvements

- Second-order secular (J2^2) or mean-element propagation near the critical
  inclination to explain the +325 s/orbit excess quantitatively.
- Frozen-orbit families via J3 (perigee libration, non-critical frozen inclinations).
- Ground-track lattice explorer over (m, k) with drag/J2-coupled repeat decay.

---

### Reproducibility Notes

- `uv.lock` pins exact dependency versions.
- Command: `uv sync && uv run pytest && uv run python experiment.py`
- Two independent runs produce identical `results.results` payloads (only
  `meta.timestamp_utc/git_commit` differ); figures regenerate deterministically
  from committed data given the pinned matplotlib version.
