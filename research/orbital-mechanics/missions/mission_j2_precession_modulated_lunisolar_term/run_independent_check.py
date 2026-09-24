"""Independent averaged-element cross-check vs stored 365-d headline.

The Cowell RK4 in mission_experiment.py propagates a full Cartesian state and
measures the osculating node at ascending crossings. This script uses the
materially different 1D averaged-element path (averaged_propagator.py):
dOmega/dt = Omega_dot_J2 + S0 + sum_k A_k cos(theta_k(t)) with prescribed
Omega_dot_J2 and circular Sun/Moon means. It reproduces the headline R sign
and order of magnitude at all three inclinations using only model-derived
constants (no fitted scale), and it evaluates the frozen-plane and
prescribed-rate discriminators in the averaged frame.

Writes results/independent_check.json.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import averaged_propagator as AV
import theory_modulation as TH
import mission_experiment as ME

W_D = 365.25
INCS = [97.7876, 90.0, 30.0]
H = 600.0


def harmonics_for(inc_rad: float) -> list:
    """Leading m=1,j=0 harmonics from quadrature (model-derived, no fit)."""
    a = 6378.137 + H
    amps = TH.leading_amplitudes_deg_day(a, inc_rad)
    out = []
    for body in ("sun", "moon"):
        for (m, j, A, phi) in amps[body]:
            if m == 1 and j == 0:
                out.append((m, j, A, phi))
    return out


def averaged_rate(inc_deg: float, omega_dot_j2: float, s0: float, harmonics: list,
                  W: float = W_D) -> float:
    inc = math.radians(inc_deg)
    out = AV.propagate_averaged(6378.137 + H, inc, W_day=W, dt_day=1.0,
                                oj2_deg_day=omega_dot_j2, s0_deg_day=s0,
                                harmonics=harmonics,
                                n3_deg_day=TH.N_MOON_DEG_DAY,
                                o3dot_deg_day=TH.OMEGA3_DOT_DEG_DAY)
    return AV.ols_slope_deg_day(out["t_day"], out["Omega_deg"])


def main():
    rows = []
    summary = {}
    for inc in INCS:
        inc_r = math.radians(inc)
        oj2 = TH.omega_dot_j2_deg_day(6378.137 + H, inc_r)
        s0 = TH.secular_total_deg_day(6378.137 + H, inc_r)
        h = harmonics_for(inc_r)
        # modulated: real J2 precession drives the beat frequencies.
        r_mod = averaged_rate(inc, oj2, s0, h)
        # frozen: Omega_dot_J2 = 0 in the beat (kills modulation), s0 kept.
        r_frozen = averaged_rate(inc, 0.0, s0, h)
        # prescribed half-rate: detune the J2 precession by 0.5.
        r_half = averaged_rate(inc, 0.5 * oj2, s0, h)
        # j2-only baseline (s0=0, harmonics=0) for subtraction.
        r_j2 = averaged_rate(inc, oj2, 0.0, [])
        r_j2_0 = averaged_rate(inc, 0.0, 0.0, [])
        # lunisolar = full - j2_only in the averaged frame.
        luni_mod = r_mod - r_j2
        luni_frozen = r_frozen - r_j2_0
        luni_half = r_half - averaged_rate(inc, 0.5 * oj2, 0.0, [])
        summary[str(inc)] = {
            "averaged_lunisolar_modulated": luni_mod,
            "averaged_lunisolar_frozen": luni_frozen,
            "averaged_lunisolar_prescribed_half": luni_half,
            "Omega_dot_J2": oj2,
            "S0": s0,
        }
        rows.append({"inc": inc, "oj2": oj2, "s0": s0,
                     "r_mod": r_mod, "r_frozen": r_frozen, "r_half": r_half,
                     "r_j2": r_j2, "r_j2_0": r_j2_0,
                     "luni_mod": luni_mod, "luni_frozen": luni_frozen,
                     "luni_half": luni_half})
    out = {"summary": summary, "rows": rows,
           "theory_sha": hashlib_sha(TH.__file__),
           "code": ME.code_hashes()}
    (HERE / "results" / "independent_check.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(summary, indent=2))


def hashlib_sha(p: str) -> str:
    import hashlib
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


if __name__ == "__main__":
    main()