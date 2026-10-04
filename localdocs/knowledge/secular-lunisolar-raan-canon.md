# Secular lunisolar RAAN rate at LEO — the canon, and how it was settled

**Status: VERIFIED (E3, cross-method) for the inclination law, sign and amplitude.
VERIFIED-WITH-LIMITATION for the "measurable from an orbit" claim — see §6.**
**Date:** 2026-10-04. **Mission:** `mission_secular_identification`.
**Supersedes** the open question recorded in `lunisolar-secular-limit-020.md`,
`lunisolar-closure-021.md` and `mean-element-leakage.md`.

---

## 1. The answer

The secular nodal rate caused by point-mass Sun and Moon acting on a near-circular LEO orbit is

```
Ω̇ = Σ₃ (3/8)·n·(μ₃/μ)·(a/a₃)³·sin2(i − i₃)/sin i

    Sun   i₃ = 23.4393°   (obliquity of the ecliptic)
    Moon  i₃ = 28.5843°   (obliquity + the Moon's 5.145° inclination to the ecliptic)
```

with `n = sqrt(μ/a³)` and `i` an **equatorial** inclination. At h = 600 km (a = 6978.137 km):

| i | Ω̇ (deg/day) |
|---|---|
| 97.7876° (SSO) | +1.3476e-04 |
| 82.2124° (SSO twin) | +2.0340e-04 |
| 89.5° | +1.7610e-04 |
| 90.5° | +1.7168e-04 |
| 30° | +1.462e-04 |

Three properties are load-bearing and were each verified independently:

1. **The inclination law is `sin2(i − i₃)/sin i`, NOT `cos i`.**
2. **The sign is prograde (+)** at SSO, for both the Moon and the Sun contributions.
3. **The twin pair at 82.21° and 97.79° has the SAME sign**, with a 1.51× magnitude gap.

## 2. How it was settled (two independent lines)

### 2.1 Exact-potential quadrature (independent of any orbital data)

Double-average the **exact** point-mass potential `−μ₃/|r_b − r|` over the satellite's mean anomaly
and the third body's, then apply Lagrange's node equation. Converged to 8 significant figures.

| i | referee | FORM-1 | ratio | FORM-2 |
|---|---|---|---|---|
| 97.7876 | 1.347261e-04 | 1.347561e-04 | 0.9998 | +2.63e-05 |
| 82.2124 | 2.033885e-04 | 2.034029e-04 | 0.9999 | −2.63e-05 |
| 89.5000 | 1.760734e-04 | 1.760967e-04 | 0.9999 | −1.70e-06 |
| 90.5000 | 1.716536e-04 | 1.716758e-04 | 0.9999 | +1.70e-06 |

Twin ratio: **referee 1.5096**, FORM-1 **1.5094**, FORM-2 **−1.0000 (opposite sign)**.

### 2.2 The orbital measurement

Independent of the above; see §6 for what it does and does not add.

### 2.3 The competing derivation is refuted

FORM-2 (a 2026-10-04 competitor, never previously compared to anything) predicts
`Ω̇ = C(a)·cos i`. It is **wrong** for a point-mass third body on a circular orbit:

- it predicts the SSO twin pair to be equal and **opposite** in sign — the exact potential says
  they are equal and the **same** sign;
- it is ~4× too small in magnitude at SSO;
- at 89.5° it predicts −2.6e-06 where the referee gives +1.76e-04.

FORM-2 should be retired from consideration. It was never published and has no committed artifact
depending on it.

## 3. Two traps in the historical record, now closed

### 3.1 The handoff values are a SUM over both bodies

Reading "FORM-1 at SSO = +1.3476e-4 deg/day" as *the Moon's contribution alone* gives +9.91e-5
(i₃ = 28.584°) or +5.35e-5 (i₃ = 18.294°). Neither matches. The quoted values are the **Sun + Moon
sum**. Verified to 4 significant figures at all three quoted inclinations.

### 3.2 Exp 018's code comment mislabels the frame — the formula is fine

`corrected_secular_lunisolar_raan_rate_rad_s` annotates `i3_moon` as the lunar "secular mean
inclination to the **ecliptic**". The Moon is inclined to the ecliptic by 5.145°, so 28.584° is its
inclination to the **equator** — which is the same frame as the equatorial `i` from
`sso_inclination_rad`. The formula mixes nothing; the comment describes the quantity wrongly.
Anyone re-deriving from that comment will "discover" a frame error that does not exist.

## 4. Measurement machinery (reusable)

### 4.1 Exact instantaneous nodal rate

```
Ω̇ = ( h_y·(−ḣ_x) + h_x·ḣ_y ) / (h_x² + h_y²),    h = r × v,   ḣ = r × a
```

No angle, no unwrapping, no node detection, no window. Exact for a Newtonian field.
Gates: Kepler → 2.5e-13 deg/day; J2 arc mean → 0.49 % of analytic; dt=30 s → 0.04 % of the signal.

**Never use i = 90°.** The denominator is `|h|²cos²i`, so at exactly 90° the node line does not
exist. `1/cos²i` is 54 at SSO but **1.3e4** at 89.5°. This is why the 2026-10-03 report's "108.45 %
correction" at i = 90° was a sign reversal from conditioning failure, not physics. Use
**89.5° / 90.5°**.

### 4.2 The estimator

Joint fit of secular + forced terms with **fixed physical frequencies**
(lunar nodal 6798.383 d, annual 365.256363 d, anomalistic month, evection). Recovers a known
1.3476e-4 deg/day signal to **<1e-12 relative** at ≥2 cycles; VIF 2.55 → 1.07–1.18.

**Column equilibration is mandatory.** cond(A) = 1.8e4 on a multi-year arc; plain `lstsq`
annihilated the secular coefficient, returning 1e-24 for a 1.35e-4 signal.

### 4.3 Window discipline

- A bare OLS slope is wrong by ~100 % on any contaminated multi-year arc.
- A bare cycle-mean is exact **only for the lunar nodal term alone** at an integer number of nodal
  cycles (2026-10-04's exactness result is correct, but its domain is narrower than first stated:
  two nodal cycles is not an integer number of years or months, so a signal carrying an annual term
  is off by ~140 %).
- **Sub-windows are bit-identical truncations of a longer fixed-step RK4 run.** Verified. So one
  propagation serves an entire window ladder at 1× cost instead of 5×.

## 5. The three defects that would have produced a wrong answer here

Each is a measurement artefact masquerading as physics — the dominant failure mode of this chain.

| defect | what it would have produced |
|---|---|
| raw-response sha256 as the data gate | every data re-acquisition "fails" for a non-data reason |
| no column equilibration | secular coefficient → 1e-24; a confident zero |
| dt = 120 s (or 240 s) | 1.37 % error (or total garbage) on a 1.35e-4 deg/day signal |

## 6. Limitations — what is NOT settled

- **The referee's agreement is at quadrupole order in (a/a₃)⁴.** Higher orders are ~1e-4 relative,
  below the comparison's precision.
- **FORM-1 and the referee both assume a circular third-body orbit.** Real lunar e₃/i₃(t) is a
  time-dependent modulation; whether it changes the secular *structure* or only the amplitude is a
  follow-up question.
- **The external check is still missing.** A published SSO (Landsat-family) measured nodal drift
  rate would be the decisive independent confirmation. Web access was unavailable in this
  environment (HTTP 426; fallback engines presented CAPTCHAs, not circumvented), so no literature
  verification was performed and none is claimed.
- **The multi-phase ≥2-cycle orbital campaign is the remaining measurement**, and its outcome does
  not change §1–§3 (those are settled by the referee) but does test whether a real orbit reproduces
  the amplitude.

## 7. Related

`lunisolar-secular-limit-020.md` (the W → ∞ question), `lunisolar-closure-021.md`,
`mean-element-leakage.md` (the one-period pathology), `j2-lunisolar-coupling.md` (retracted
mechanism), `j2-precession-modulated-lunisolar-term.md` (H-mod falsified).
Report: `localdocs/reports/mission-secular-identification-2026-10-04.md`.