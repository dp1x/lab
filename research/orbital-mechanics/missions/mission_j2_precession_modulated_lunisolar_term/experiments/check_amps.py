import math, sys, time
sys.path.insert(0, 'src')
sys.path.insert(0, 'research/orbital-mechanics/missions/mission_j2_precession_modulated_lunisolar_term')
import theory_modulation as T
a = 6378.137 + 600
for inc in [97.7876, 90.0, 30.0]:
    i = math.radians(inc)
    t0 = time.monotonic()
    amps = T.leading_amplitudes_deg_day(a, i)
    dt = time.monotonic() - t0
    print(f"i={inc} ({dt:.1f}s)")
    for body in ("sun", "moon"):
        print(f" {body}:")
        for (m, j, A, phi) in amps[body]:
            print(f"  m={m} j={j} A={A:.3e} phi={phi:.3f}")
    pr = T.predicted_R_deg_day(a, i, 365.25, phi_sun=0.0, phi_moon=0.0)
    print(f" predicted B_total(365d,phi=0)={pr['B_total']:.3e} (S0={pr['S0']:.3e} A_sun={pr['A_sun']:.3e} A_moon={pr['A_moon']:.3e})")
