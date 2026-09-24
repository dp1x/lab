# DERIVATION — J2-precession-modulated lunisolar node response (lead track)

**Date:** 2026-09-05. **Status:** frozen theory before final 18.6-yr validation.
**Companion independent track:** subagent Kaula derivation (reconciled in final report).
**No empirical coefficient is fitted to RK4 output.** All constants come from
`src/lab_utils` (MU, J2, R_EARTH), mission constants (MU_SUN/MOON, AU), and
ephemeris frequencies. Amplitudes come from the quadrupole disturbing function
via quadrature/Fourier analysis, not from regressing measured `R`.

---

## 1. Setup, conventions, observable

- Inertial ECI mean-of-date, TDB. Satellite circular LEO (`e≈0`), `a=R_EARTH+h`.
- Mean motion `n=sqrt(mu/a³)`. J2 secular (first-order, circular):
  `Omega_dot_J2 = -(3/2)·n·J2·(R/p)²·cos i`, `p=a(1−e²)=a`.
- Third-body quadrupole disturbing function (geocentric, direct+indirect;
  indirect averages to zero for secular/long-period at quadrupole order):
  `R_3b = mu3·[1/|r−r3| − r·r3/r3³]`, quadrupole part
  `R_quad = (mu3/a3)·(a/a3)²·P2(cos S)`, `P2(x)=(3x²−1)/2`,
  `cos S = r_hat·r3_hat`.
- Lagrange (circular): `dOmega/dt = [n·a²·sin i]⁻¹·∂R/∂i`
  (exact `1/(n·a²·sqrt(1−e²)·sin i)`; `sqrt→1` for `e=0`).
- Numerical observable (osculating, ascending-node OLS):
  `s_hat(W) = OLS slope of unwrapped Omega_node(t) over window W`.
  Theory predicts the **mean-element** rate; the map to osculating OLS is
  through the finite-window bias formula below plus short-period averaging
  (orbital ~96 min averages to ~0 at node crossings; evection/variation/
  annual/monthly/nodal remain as `A_k` terms).

---

## 2. Singly-averaged rate retaining Omega

Write inertial directions with `R_z(Omega)·R_x(i)·R_z(u)` (`u` argument of
latitude; `u=M` for circular `omega=0`):

```
r_hat(Omega,i,u),  r3_hat(Omega3,i3,u3).
cos S = r_hat·r3_hat
      = c1(i,i3,Omega−Omega3)·cos u·cos u3 + s-terms in sin u·sin u3, etc.
```

`P2(cos S)` contains harmonics `0, ±u±u3, ±2u±2u3` crossed with
`0, ±(Omega−Omega3), ±2(Omega−Omega3)`. Averaging over satellite `u∈[0,2π)`
(single average, holding `Omega,i,u3,Omega3,i3` fixed) kills all `u`-dependent
terms. What remains is

```
R_bar(Omega,i; u3,Omega3,i3) = (mu3·a²/a3³)·Σ_{m=0..2} Σ_{j} C_{mj}(i,i3)·cos(m(Omega−Omega3)+j·u3+phi).
```

Differentiating w.r.t `i` and dividing by `n·a²·sin i` gives the
**singly-averaged node rate**

```
Omega_dot_3b(t) = S0(i,i3) + Σ_k A_k(i,i3)·cos(theta_k(t)),   (1)
theta_k = m_k·(Omega(t)−Omega3(t)) + j_k·u3(t) + phi_k,
```

with `S0` the `m=0,j=0` (doubly-averaged secular) piece. `S0` is exactly the
lab Convention-B quadrupole

```
S0 = (3/8)·n·(mu3/mu)·(a/a3)³·sin2(i−i3)/sin i,   (2)
```

summed over Sun (`i3=obliquity`, `a3=AU`) and Moon (`i3≈mean lunar inclination
to equator`, `a3≈384400 km`). Prograde (`>0`) at all three headline
inclinations for 2026 geometry. Dimensions: `n[rad/s]→deg/day` conversion
explicit in code; `∂R/∂i / (n·a²·sin i)` is `1/time`. Limit checks:
`S0→0` as `a→0` (`(a/a3)³`), as `mu3→0`, as `sin2(i−i3)→0` (`i=i3` or
`i=i3+90°` gives zero secular but nonzero long-period — important
discriminator at `i=90°`).

`A_k` are the same order as `S0`: each `A_k = n·(mu3/mu)·(a/a3)³·F_k(i,i3)/sin i`
with dimensionless `F_k=O(0.1–1)` from `∂C/∂i`. They are computed in
`theory_modulation.py` by **quadrature over u and Fourier analysis in
(Omega−Omega3, u3)**, not by fitting RK4 `R`. Eccentricities (`e3≈0.05` lunar,
`0.0167` solar) and `i3(t)` nodal variation enter as `O(e3)`/`O(delta i3)`
corrections to `A_k` (bounded in report; synthetic-circular control isolates them).

J2 enters (1) **only** through `Omega(t)`:

```
Omega(t) = Omega0 + Omega_dot_J2·t + (lunisolar-driven Omega drift, O(3b), neglected at leading modulation order).
```

This is an **adiabatic/resummation** treatment, not a conventional
second-order Lie cross term. Its perturbation order is: first-order in 3body
for `A_k,S0`, with `Omega(t)` resumming all-orders J2 secular phase. Taylor
expansion for small `Omega_dot_J2·W` yields apparent mixed scaling
`∝ lambda_J2·lambda_3body` (see §4) without any true `J2·3b` secular term.
The direct Lie cross term `∝ J2·(n3/n)² ~ 1e-9` relative is separately bounded
and negligible — H-x is rejected on magnitude grounds before numerics.

---

## 3. Finite-window OLS bias (what the estimator actually measures)

Model the mean node angle as

```
Omega(t) = s·t + c + Σ_k (A_k/omega_k)·sin(omega_k·t+phi_k),   (3)
```

where `s = Omega_dot_J2 + S0` (secular), `omega_k = dtheta_k/dt =
m_k·(Omega_dot_J2 − Omega_dot_3) + j_k·n3 + …` (beats of J2 precession with
third-body node/mean motions), `A_k` from §2. For `omega_k→0` (resonance) the
term is `(A_k·cos phi_k)·t` (secular-like over any finite W).

For dense sampling on `[0,W]`, the OLS slope bias from harmonic `k` is exactly

```
B_k = A_k · g(omega_k·W, phi_k),   (4)
g(x,phi) = [12·( (x·cos(phi+x) − sin(phi+x) + sin phi)/x² ) − 6·(cos(phi+x)−cos phi)/x ] / … (closed form in code),
```

with transparent limits (derived by exact OLS projection of `sin` onto `t`):

- Short window `|x|=|omega_k|W << 1`: `B_k = A_k·cos(phi_k+ x/2)·[1 − x²/…] + O(x)`,
  i.e. `|B_k| ≤ |A_k|`, **scales as amplitude (∝3body) with phase set by
  epoch**, and to first order in `Omega_dot_J2` varies as
  `dB_k/dOmega_dot_J2 ∝ A_k·m_k·W·sin(…)`, giving **linear precession-rate
  dependence** and `lambda_J2×lambda_3body` scaling when `omega_k` contains
  `m_k·Omega_dot_J2`.
- Long window `|x| >> 1` (fast, non-resonant): `B_k = O(A_k/(omega_k·W))` with
  `1/W` envelope and `1/W²` for integer-cycle harmonics; averages out.
- Resonant `omega_k→0`: `B_k → A_k·cos phi_k` (quasi-secular on any
  decadal W; still long-period asymptotically, period `2π/|omega_k|→∞`).

`theory_modulation.py:ols_bias()` implements the exact `g()` and is verified
against synthetic oracles to <1% (test). No fitted constant enters.

**Three mechanisms distinguished:**
- (i) Direct cross-order secular: would contribute to `s` itself, survive
  `W→∞` and frozen-`Omega`. Bounded ~1e-9 relative → negligible.
- (ii) Precession modulation (H-mod): `B_k` with `omega_k` containing
  `Omega_dot_J2`. Vanishes for frozen `Omega` (`omega_k` loses J2 part;
  `B_k` changes as predicted), scales with prescribed `Omega_dot_p`,
  vanishes for `J2→0` at fixed W up to non-J2 beats.
- (iii) Pure aliasing (H-alias): `B_k` with `omega_k` independent of J2
  (evection `~2π/31.8d`, variation `~2π/14.8d`, annual, monthly). Survives
  frozen-plane; independent of `Omega_dot_p`. Discriminated by controls.

---

## 4. Predicted scalings and resonance

- `S0 ∝ n·(a/a3)³` (extra `1/sin i·sin2(i−i3)`).
- `A_k ∝ n·(mu3/mu)·(a/a3)³·F_k/sin i`.
- `omega_k = m_k·Omega_dot_J2 + rest_k`, `Omega_dot_J2 ∝ n·J2·(R/p)²·cos i`.
- Short-window extra slope `B_k ∝ A_k ∝ n·(a/a3)³` with J2 entering via
  `cos(phi_k)` phase and via `dB/dOmega_dot_J2 ∝ A_k·W`; the **perturbative
  cross coefficient** `a11 = ∂²s_hat/∂lambda_J2∂lambda_3body|_0 ∝
  Σ_k (∂A_k/∂lambda_3body)·(∂g/∂omega_k)·m_k·(∂Omega_dot_J2/∂lambda_J2)`
  is nonzero generically (explains why Phase-B `a11≠0` does not imply a true
  secular cross term).
- Altitude: `B_k/S0 ∝ (extra n·J2·(R/p)²·W or phase factor)` — **steeper**
  than pure `(a/a3)³`. `R(500)/R(800)` predicted ≠ `S0(500)/S0(800)`;
  tested within factor 2.
- Inclination: `F_k(i,i3)·(cos i in omega_k)` structure; secular zero-crossings
  (`sin2(i−i3)=0`) do not zero `A_k`. Extra discriminator at `i=60°,0°`.
- **Sun-sync resonance:** `omega = Omega_dot_J2 − n_sun ≈ 0` near SSO
  (`Omega_dot_J2≈+0.986 deg/day`, `n_sun≈0.9856 deg/day`). Solar `m=1`
  long-period bias persists on decadal W (quasi-secular, still long-period).
  At `i=90°` (`Omega_dot_J2=0`) solar beat is `-n_sun` (1-yr period, averages
  over 18.6-yr); lunar nodal beat dominates. Resonance predicts: frozen-plane
  kills SSO extra; detuning `Omega_dot_p` away from `n_sun` restores `1/W`
  decay. Tested explicitly.

---

## 5. Mean vs osculating

Theory (1)–(4) is **mean-element** (averaged over satellite M). Numerics sample
**osculating** `Omega` at ascending nodes (near `u≈0`, not full-`u` average).
Short-period (`~96 min`) content averages to <1e-6 deg/day at node sampling
(verified by dt ladder + cadence decimation). Remaining gap is precisely the
`A_k` long-period family above (evection/variation included as `j_k·n3`
sidebands with `A_k` from quadrature bounds). No separate "mean-to-osculating
artifact" beyond (4) is invoked; if frozen-plane + prescribed-rate tests fail,
the residual is assigned to H-alias/H-other, not to H-mod.

---

## 6. Falsifiers (binding)

H-mod is falsified if: (i) `R` does not vanish for `J2→0`+`Omega_dot_p=0`
within detection limit; (ii) `dR/dOmega_dot_p|_0 = 0` within 3σ while `|R|>0`;
(iii) frozen-plane does not reduce `|R|` by ≥4×; (iv) sign wrong at ≥2
headline inclinations; (v) magnitude off by >5× at ≥2 inclinations with no
fitted scale; (vi) an `Omega_dot_p`-independent alternative predicts `R`
within factor 2 while H-mod is >5× off. Partial support = correct
frequency-dependence but factor 3–5 magnitude or 1-inclination sign miss.

Dimensional check: all rates `deg/day`; `A_k·g` is `deg/day`; `a11` is
`deg/day` (per unit lambda²); `dR/dOmega_dot_p` dimensionless. Limits:
`R→0` as `mu3→0`, as `(a/a3)→0`, as `W→∞` non-resonant; `B_k→A_k·cos phi`
resonant. Code asserts each.

