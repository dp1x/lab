import math, sys
sys.path.insert(0, 'src')
sys.path.insert(0, 'research/orbital-mechanics/missions/mission_j2_precession_modulated_lunisolar_term')
import theory_modulation as T
a = 6378.137 + 600
for inc in [97.7876, 90.0, 30.0, 60.0]:
    i = math.radians(inc)
    s0 = T.secular_total_deg_day(a, i)
    oj2 = T.omega_dot_j2_deg_day(a, i)
    print(f"i={inc} S0={s0:.3e} oj2={oj2:.4f}")
print("oct90", T.octupole_bound_deg_day(a, math.radians(90.0)))
