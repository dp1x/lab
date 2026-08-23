# Experiment Card: J2 Precession — Secular Nodal & Apsidal Drift from Earth's Oblateness

> Status: complete
> Date: 2026-08-22
> Domain: orbital-mechanics
> Experiment dir: `research/orbital-mechanics/experiments/j2Precession/`

## Research Question

For Earth orbits (a, e, i) perturbed only by the J2 zonal term, does **full-force numerical
Cowell RK4 propagation** (`a = a_Kepler + a_J2`) independently rediscover the **first-order
analytical secular rates**

```text
Omega_dot = -(3/2) n J2 (R_E/p)^2 cos i
omega_dot =  (3/4) n J2 (R_E/p)^2 (5 cos^2 i - 1)
n = sqrt(mu/a^3),  p = a(1-e^2)
```

when rates are measured from propagated Cartesian states via an independent
state->element->trend estimator — with (A) numerical convergence to a high-accuracy
numerical reference at RK4 order ~4 cleanly separated from (B) the legitimate
first-order model-order residual? Sharpened from `localdocs/roadmap.md` ("Node drift
from Earth's bulge | Analytic secular rates"). No M_dot anywhere in this experiment
(explicitly out of contract).

## Background Theory

### Frames and units

ECI J2000 pseudo-inertial, Z = Earth spin axis; `a_J2` is axisymmetric so any Z-rotation
leaves the dynamics invariant. Units km, km^3/s^2, s; angles rad internal / deg I/O.
ECEF/GMST not needed for pure J2 dynamics.

### J2 perturbed dynamics

```text
a_J2 = -(3 mu J2 R_E^2)/(2 r^5) * [ x(1 - 5z^2/r^2),  y(1 - 5z^2/r^2),  z(3 - 5z^2/r^2) ]
```

Verified as the exact gradient of the static potential
`U_J2 = +mu J2 R_E^2 P2(z/r)/r^3` (`P2(u) = (3u^2-1)/2`): componentwise differentiation
reproduces all three acceleration signs. Consequences used as invariants: **total specific
energy** `E = v^2/2 - mu/r + U_J2` is exactly conserved by the ODE (integrator drift only),
and **h_z = (r x v)_z is exactly conserved** (axisymmetry) while |h| physically oscillates at
O(J2) — an axisymmetry signature the tests assert directly (h_z range / |h| range ~ 1e-6..1e-8).
|a_J2|/|a_kepler| ~ 1.43e-3 at LEO (r = 6798 km).

### Secular oracle and its limits

The double-averaged (mean-element) first-order rates are the two boxed formulas above.
`omega_dot` vanishes at the critical inclinations `i = arccos(1/sqrt(5)) = 63.43494882 deg`
and `116.56505118 deg`; `(3/4)(4-5sin^2 i) == (3/4)(5cos^2 i - 1)` identically.
SSO condition: `cos i_SSO = -Omega_dot_target p^2 / (1.5 n J2 R_E^2)` with target
`360/365.2422 = 0.9856473321 deg/day` (exact mean-solar-year quotient; the commonly printed
0.98564736 corresponds to the tropical-year variant 365.24219 d — the two differ by
2.8e-8 deg/day, negligible vs every tolerance here; sidereal-year rate kept separate).

**Mean vs osculating (the core distinction).** The oracle consumes MEAN elements; the
estimator recovers OSCULATING elements. Even at zero integration error they differ at
O(J2): short-period (T, T/2) oscillations average out over integer-orbit windows, but the
constant mean-minus-osculating element offset remains (implied delta-a/a ~ -1.2e-3 gives
delta-Omega_dot/Omega_dot = -3.5 delta_a/a ~ +0.4%), plus second-order-in-J2 secular terms
with small divisors `(1 - 5cos^2 i)` near the critical inclination. This residual PLATEAUS
with step refinement (verified: it does NOT halve when dt halves) and is reported as
model-order difference, never as integration error.

## References

(title/edition/chapter-level only — equation numbers omitted where they could not be
verified against a specific printing; citation != truth)

- D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed., Microcosm,
  2013 — Ch.9 general perturbations, central-body J2 secular analysis.
- H. D. Curtis, *Orbital Mechanics for Engineering Students*, 4th ed., Elsevier, 2021 —
  Ch.10 introduction to orbital perturbations (J2 nodal/apsidal drift).
- R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*, Dover, 1971 —
  Ch.9 perturbations (nodal regression).
- C. D. Murray & S. F. Dermott, *Solar System Dynamics*, Cambridge UP, 1999 — Ch.6
  disturbing function, planetary oblateness effects.
- NIMA, *WGS-84*, TR8350.2 — R_E = 6378.137 km; J2 = sqrt(5)|C20_bar| =
  1.082629821e-3 (C20_bar = -0.484166774985e-3). EGM2008's 1.08262668e-3 documented, NOT used.
- IAU 2015 Resolution B3 (arXiv:1510.07674) — nominal GM_E = 398600.4418 km^3/s^2.
- Exp 002 RK4 machinery; Exp 006 verified 3-D Cowell; Exp 008 constants/element helpers — reused.

## Assumptions

* Pure first-order J2 acceleration; no J3/J4/tesserals/drag/SRP/lunisolar — **idealization** (Exp 010+).
* WGS-84 constant set (mu IAU nominal, R_E, J2 via C20_bar) used coherently — **verified** (provenance in results.json; km/m firewall test).
* First-order secular formulas evaluated at INITIAL OSCULATING elements — **convention** (the resulting O(J2) offset is measured and reported, not absorbed).
* omega_dot claimed only for seed e >= 0.01 — **policy**: J2 induces real eccentricity (~1e-3 at LEO) that dominates smaller seeds, making recovered omega sweep once per orbit with no secular content.
* Classical elements only: i = 0/180 -> RAAN structurally undefined (no Omega_dot claim); omega, measured FROM the node, is undefined there too (longitude of periapsis would be the alternative — out of scope).
* Fixed-step RK4, uniform dt per case — **idealization** (non-symplectic; energy drift <= 4.9e-9 rel over 100 orbits).

## Methodology

1. **Constants pinned** with provenance strings (results.json:constants); one coherent set everywhere.
2. **Analytical oracle** `analytic_rates(a,e,i)` consumes ONLY (a,e,i); independent inline reimplementation in tests.
3. **Propagator** `propagate_3d_rk4_j2`: loop structure and j2==0 expression operation-for-operation identical to verified Exp 006 `propagate_3d_rk4` — regression-tested BIT-EXACT (`np.array_equal`) against it.
4. **Resolution rule** (documented, per case): `steps_per_orbit = max(512, ceil(720/(1-e)^{3/2}))`; one consistent timestep per case (Molniya e=0.74 -> 5432/orbit).
5. **Estimator** (independent of oracle algebra): per-sample osculating elements from h = rxv, node = zhat x h, e_vec = (v x h)/mu - r/r; np.unwrap (safe by ~7 decades at documented resolution); closed-form OLS over integer-orbit windows; windows 20/50/100 orbits ALL reported, primary = longest declared a priori; stabilization criterion max pairwise rel diff < 1e-3.
6. **Third-path estimator**: ascending-node-crossing inertial longitudes vs parabolically refined crossing times — first-order Omega short-period terms vanish at nodes, giving clean nulls (polar case measures 1e-16 deg/day) and independence from element algebra (agrees with primary estimator to 7.7e-4 rel at 100 orbits).
7. **Convergence protocol**: grids 128..1024 steps/orbit vs 2048/orbit NUMERICAL reference, 20 orbits, identical sample phases across grids (stride subsampling) so short-period sampling cancels; rate-metric orders AND raw integrator order via `kepler_order_check` (max full-vector position error vs closed-form Kepler truth — phase-sensitive, unlike final-|r| which hides along-track error and decays at order ~5).
8. **Adversarial layer**: J2=0 nulls; J2->-J2 sign-flip sensitivity; wrong-p discriminator (Molniya, bug signature (1-e^2)^{-2} = 4.89x); pathological i x e grid {0, 63.435, 90, 116.565, 180} x {0, 0.05, 0.2, 0.74}.

## Implementation

- Script: `experiment.py`
- Language/runtime: Python 3.12, numpy, matplotlib (Agg before pyplot import)
- Runtime: `$REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/j2Precession/experiment.py` (~53 s single core)
- Determinism: pure float64, no RNG, fixed grids; two full runs byte-identical apart from `meta.timestamp_utc`; figure MD5s identical across runs
- Reuse: Exp 006 `propagate_3d_rk4` (importlib `pcm_006_for_j2`, J2=0 bit-exact oracle); Exp 008 Kepler/element helpers + constants (importlib `groundtracks_for_j2`); `src/lab_utils/results.py` save_json_result; `src/lab_utils/metrics.py` convergence_rate (eps-guarded). No scaffolding rebuilt; no previous experiment modified.
- Dependencies: numpy, matplotlib only (already pinned in uv.lock).

## Validation Method

32 focused pytest tests (`tests/test_j2_precession.py`, importlib-loaded, independence doctrine:
expected values derived inline from theory or separate code paths, never by calling the same helper):
L1 analytical identities (J2=0 nulls, polar zero, exact-critical zero, sign conventions,
form identity, a^{-7/2} scaling, e-only-through-p identity (1-e^2)^{-2}, SSO table vs anchors);
L2 element/state consistency (round-trip through imported coe_to_rv_eci across quadrants,
M0->E0->nu0 chain, energy + h_z axisymmetry signature, km/m firewall);
L3 numerical-vs-reference (bit-exact vs 006, anchor rates, node-crossing agreement,
window stabilization, residual-plateau proof that model-order error is not integration error);
L4 convergence (Kepler-truth order band [3.6,4.4]; rate-metric >= 3.6 monotone with floor check);
L5 anchors (Starlink, SSO target rate, polar null < 1e-8 deg/day via node crossings,
critical-i apsidal freeze, Molniya incl. wrong-p discriminator);
L6 pathological grid sentinels; L7 adversarial regressions (J2=0 drift bound, sign-flip ratios,
deg/day twin consistency, determinism).

Full repository suite: **290 passed** (258 pre-existing + 32 new).

## Results

Headline (full tables in `results/results.json`; rates deg/day):

| case | a [km] | e | i [deg] | analytic Omega_dot | numeric Omega_dot | resid rel |
|------|--------|---|---------|--------------------|-------------------|-----------|
| ISS | 6798.137 | 0.0003 | 51.6 | -4.951018 | -4.972394 | 4.3e-3 |
| STARLINK | 6928.137 | 0.0003 | 53.0 | -4.489207 | -4.507991 | 4.2e-3 |
| SSO600 (solved i_SSO=97.787647) | 6978.137 | 0 | 97.787647 | +0.9856473 | +0.990102 | 4.5e-3 |
| POLAR | 6878.137 | 0 | 90 | 0 (exact) | ~0 (abs 2.2e-18) | null |
| MOLNIYA | 26560.0 | 0.74 | 63.4 | -0.147933 | -0.146589 | 9.1e-3 |
| CRITICAL | 6878.137 | 0.2 | 63.434949 | -3.712691 | -3.737927 | 6.8e-3 |

Apsidal: ECC_REF (i=30, e=0.2) omega_dot -11.41 analytic vs -11.49 numeric (+0.67%);
CRITICAL omega_dot = 0 analytic vs +6.36e-3 numeric (estimator floor, suppression >185x
vs typical LEO signal); MOLNIYA omega_dot +4.03e-4 vs +3.11e-4 (near-critical small-divisor regime).

* **Physics residual structure**: systematic +0.42..0.68% relative offset, SAME SIGN across
  all seven cases (|numeric| > |analytic|), stable across 20/50/100-orbit windows
  (stabilization 1.8e-4 for LEO) — identified as the O(J2(R_p)^2) mean-vs-osculating offset
  (implied delta-a/a ~ -1.2e-3) plus second-order small-divisor amplification near critical
  inclination (MOLNIYA 0.91%, CRITICAL 0.68%). Plateau proof: halving dt changes the residual
  by <50% while integration-error metrics shrink 16x+. NOT integration error.
* **Convergence**: Kepler-truth order 4.092 (intervals 4.163/4.077/4.038, band [3.6,4.4]);
  rate-metric orders 4.72/4.56/4.47 (super-fourth: orbit-averaging cancels RK4's leading
  phase-error mode; rule: every interval >= 3.6, floor check vs round-off).
* **Nulls/adversarial**: J2=0 slopes -4.2e-14 (Omega) and 1.9e-6 deg/day (omega, documented
  1/e noise floor) vs frame-bug artifact scale O(0.98); sign-flip ratios -1.0088 (Omega),
  -1.0125 (omega) — deviation from exactly -1 is the even-in-J2 second-order secular term,
  not estimator asymmetry (identical at 10 and 25 orbits).
* **Pathological grid**: 20 rows all finite/bounded/sentinel-correct; equatorial rows have NO
  RAAN claim (node undefined => omega undefined too); circular seeds develop measured induced
  eccentricity +/-9.60 km = a*(3/2)J2(R_E/p)^2 (band width 5e-3 justified from this).
* **SSO table (from pinned constants)**: 500 km -> 97.401786, 600 km -> 97.787647,
  800 km -> 98.603085 deg; integrated SSO600 reproduces the solar-rate target to 0.45%.

Figures (`results/figures/`, regenerated deterministically from results.json data, dpi=150):
`raan_vs_time_fit.png`, `omega_vs_time_fit.png`, `convergence_order.png`,
`analytic_vs_numeric.png`.

## Limitations

* First-order J2 only: the measured residuals are dominated by model order (mean-vs-osculating
  offset, second-order secular small divisors near i_crit), quantified but not modeled away —
  a Brouwer-Lyddane mean-element layer would close most of the gap and is deliberately out of scope.
* Critical-inclination apsidal freeze demonstrated to |omega_dot| < 0.02 deg/day (estimator
  leakage floor 1.91*A/N with short-period amplitude ~0.25 deg at e=0.2); the EXACT zero is
  established analytically, not numerically.
* Molniya 12-orbit window still converging toward the 48-orbit primary (~4e-3 pairwise);
  high-e windows are expensive, so the residual is reported with this caveat instead of
  tightened by cherry-picking.
* omega trend claims restricted to seed e >= 0.01; ISS/Starlink omega fields record the
  induced-eccentricity sweep (slope ~= mean motion) flagged with an explicit note.
* Non-symplectic fixed-step RK4; energy drift <= 4.9e-9 rel over 100 orbits — adequate for
  secular-rate extraction over <= 100 orbits, not for year-long ephemerides (Exp 013).
* Anchors are MODEL predictions from published orbital elements, not TLE observations;
  comparison against real catalogs is Exp 013 scope.

## Future Improvements

* Brouwer-Lyddane mean-element conversion to separate "secular theory at mean elements"
  from higher-order terms (kills the systematic 0.4% plateau).
* J3 (and Kozai-style long-period dynamics) — breaks axisymmetry, gives h_z something to do.
* Nodal-period/anomalistic-period corrections and repeat-ground-track coupling back to Exp 008.
* Drag (Exp 010) reuses the same generalized-propagator pattern with a non-conservative force
  (energy invariant then becomes a diagnostic, not an invariant).

---

### Reproducibility Notes

* `uv.lock` pins exact dependency versions.
* Command: `$REPO_ROOT\.venv\Scripts\python.exe -m pytest && $REPO_ROOT\.venv\Scripts\python.exe research/orbital-mechanics/experiments/j2Precession/experiment.py`
  (generic `uv sync && uv run pytest && uv run python ...` when uv is on PATH).
* Figures regenerate deterministically from recorded data in results.json (Agg, fixed grids, no RNG).
* Provenance: `results/results.json:meta` stores name, description, timestamp_utc, git_commit,
  python_version via `src/lab_utils/results.py`.
