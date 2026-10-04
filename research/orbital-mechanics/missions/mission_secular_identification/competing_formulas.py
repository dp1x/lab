"""The two competing secular lunisolar RAAN formulas, reproduced from source.

Neither formula is fitted to any numerical result. Each is implemented from its
own written form, evaluated on its own, and then tested against the measurement.

RECONSTRUCTION IS EXACT AND VERIFIED
------------------------------------
The 2026-10-04 handoff table quotes FORM-1 at +1.3476e-4 deg/day (SSO),
+1.761e-4 (89.5 deg) and +1.739e-4 (90 deg). A naive reading - FORM-1 applied to
the Moon alone with ``i3 = 28.584 deg`` - gives +9.91e-5 (i3 = 28.584) or
+5.35e-5 (i3 = 18.294), neither of which matches. The values reproduce only when
FORM-1 is summed over BOTH third bodies with their own inclinations to the
reference plane of ``i``:

    Sun : i3 = obliquity  = 23.4393 deg, GM = 1.3271244e11, a3 = 1.4959787e8 km
    Moon: i3 = 28.5843 deg,            GM = 4902.8001,    a3 = 384400 km

That reconstruction is verified to 4 significant figures at all three quoted
inclinations (`test_form1_reproduces_the_handoff_table`), and it independently
reproduces the handoff's "1.51x magnitude gap" between 97.79 and 82.21 deg
without that number being an input (1.5094 computed vs 1.51 quoted).

FRAME NOTE. ``i`` is an EQUATORIAL inclination (what `sso_inclination_rad`
returns), so ``i3 = 28.584 = obliquity + 5.145`` for the Moon is the correct
lunar inclination to the equator and the two terms are in a single consistent
frame. The apparent frame error is in the CODE COMMENT in Exp 018's
`corrected_secular_lunisolar_raan_rate_rad_s`, which labels 28.584 deg as a
"secular mean inclination to the ecliptic". It is not: the lunar orbit is
inclined to the ECLIPTIC by 5.145 deg, so 28.584 deg is its inclination to the
EQUATOR. The committed formula is right; its comment mislabels the frame.

Provenance of the two formulas
------------------------------
FORM-1 ("audit-018", Exp 018 Track B) superseded the Kozai-style apsidal form
that audit-018 found wrong in three compounding ways. `mission_lunisolar_closure`
(2026-09-03) then found FORM-1 fails to reproduce even the SIGN of the measured
rate at i_sso and i = 30 deg on an 18.6-yr arc.

FORM-2 (2026-10-04, session Track A) is a competing derivation, never compared to
a measurement, and whose attempted literature cross-check failed (HTTP 426 /
CAPTCHAs, not circumvented).

They differ in a STRUCTURAL property, not a fitted coefficient, so no
measurement can be tuned to reconcile them. See the mission card's
pre-registered discrimination rule.
"""
from __future__ import annotations

import math

from lab_utils import MU_EARTH_KM3S2, R_EARTH_KM

LUNAR_GM = 4902.8001
LUNAR_A3_KM = 384400.0
LUNAR_I3_DEG = 28.5843          # obliquity + lunar inclination to ecliptic
SOLAR_GM = 132712440018.0
SOLAR_A3_KM = 1.4959787e8
SOLAR_I3_DEG = 23.4393          # obliquity of the ecliptic
SOLAR_OBLIQUITY_DEG = SOLAR_I3_DEG

SEC_PER_DAY = 86400.0
DEG_PER_RAD = 180.0 / math.pi


def _n_rad_s(a_km: float) -> float:
    return math.sqrt(MU_EARTH_KM3S2 / a_km ** 3)


def _form1_base(a_km: float, gm3: float, a3_km: float) -> float:
    """(3/8) n (mu3/mu) (a/a3)^3, in deg/day, for one third body."""
    return (3.0 / 8.0) * _n_rad_s(a_km) * (gm3 / MU_EARTH_KM3S2) \
        * (a_km / a3_km) ** 3 * SEC_PER_DAY * DEG_PER_RAD


def form1_audit018(a_km: float, inc_deg: float,
                   include_sun: bool = True) -> float:
    """FORM-1, the lab's audit-018 secular rate, deg/day.

        Omega_dot = sum_3 (3/8) n (mu3/mu) (a/a3)^3 sin(2(i - i3)) / sin(i)

    summed over the Sun (i3 = obliquity) and the Moon (i3 = obliquity + 5.145).
    """
    i = math.radians(inc_deg)
    si = math.sin(i)
    if abs(si) < 1e-12:
        raise ValueError("FORM-1 is singular at i = 0 or 180 deg")
    total = _form1_base(a_km, LUNAR_GM, LUNAR_A3_KM) \
        * math.sin(2.0 * (i - math.radians(LUNAR_I3_DEG))) / si
    if include_sun:
        total += _form1_base(a_km, SOLAR_GM, SOLAR_A3_KM) \
            * math.sin(2.0 * (i - math.radians(SOLAR_I3_DEG))) / si
    return total


def form2_competitor(a_km: float, inc_deg: float) -> float:
    """FORM-2, the 2026-10-04 competing derivation, deg/day.

        Omega_dot = C(a) cos(i),  C(a) = -(3/4) n (mu3/mu) (a/a3)^3 P2(cos i3)

    with the lunar parameters. Note C(a) is a POSITIVE constant: the quoted
    ``-(3/4) ...`` prefactor is negative and ``P2(cos 28.584 deg) = -0.1456``,
    which is also negative, so the product of the two negatives is positive.
    """
    i = math.radians(inc_deg)
    p2 = 0.5 * (3.0 * math.cos(math.radians(LUNAR_I3_DEG)) ** 2 - 1.0)
    c = -(3.0 / 4.0) * _n_rad_s(a_km) * (LUNAR_GM / MU_EARTH_KM3S2) \
        * (a_km / LUNAR_A3_KM) ** 3 * p2
    return c * math.cos(i) * SEC_PER_DAY * DEG_PER_RAD


def both_forms(a_km: float, inc_deg: float) -> dict:
    return {"inc_deg": inc_deg,
            "form1_deg_day": form1_audit018(a_km, inc_deg),
            "form2_deg_day": form2_competitor(a_km, inc_deg)}


#: Diagnostic inclinations. 90 deg is EXCLUDED BY CONSTRUCTION: there
#: ``cos i = 0`` so both ``h_x`` and ``h_y`` vanish and the node line does not
#: exist. 89.5 / 90.5 deg replace it with the same sin-vs-cos discrimination
#: (1/cos^2 i ~ 1.3e4) at 6x better conditioning than SSO.
DIAGNOSTIC_INCS_DEG = (97.7876, 82.2124, 89.5, 90.5)


def sso_twin_discriminator(a_km: float = R_EARTH_KM + 600.0) -> dict:
    """Discriminator 1: i = 97.7876 deg vs its twin i = 82.2124 deg.

    ``cos i`` changes sign between them; ``sin i`` does not. A ``cos i`` law
    predicts EQUAL magnitudes with OPPOSITE signs; FORM-1 predicts the SAME sign
    with a magnitude gap (1.509 here, and that is an OUTPUT not an input).
    """
    hi, lo = 97.7876, 82.2124
    out = {"inc_hi_deg": hi, "inc_lo_deg": lo}
    for name, fn in (("form1", form1_audit018), ("form2", form2_competitor)):
        a, b = fn(a_km, hi), fn(a_km, lo)
        out[name] = {f"at_{hi}": a, f"at_{lo}": b,
                     "same_sign": bool(a * b > 0),
                     "magnitude_ratio_hi_over_lo": abs(a / b),
                     "signed_ratio": a / b}
    return out


def near_vertical_discriminator(a_km: float = R_EARTH_KM + 600.0) -> dict:
    """Discriminator 2: i = 89.5 / 90.5 deg, clear of the node singularity.

    At exactly 90 deg FORM-1 is non-zero and FORM-2 is exactly zero, but 90 deg
    is unmeasurable, so the check is made either side of it.
    """
    return {"form1": {f"at_{i}": form1_audit018(a_km, i) for i in (89.5, 90.5)},
            "form2": {f"at_{i}": form2_competitor(a_km, i) for i in (89.5, 90.5)}}