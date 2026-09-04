__test__ = False  # not a pytest test module
"""Quick analytical diagnostic for the mission."""
import sys
sys.path.insert(0, '.')
from mission_experiment import forced_secular_lunar_nodal_node_rate_deg_day, octupole_lunisolar_raan_rate_rad_s

print('=== Forced-secular lunar nodal mode diagnostic ===')
for i_deg in [97.7876, 90.0, 30.0]:
    res = forced_secular_lunar_nodal_node_rate_deg_day(600.0, i_deg)
    print(f'i={i_deg:6.2f} deg: standard secular = {res["standard_secular_lunar_deg_day"]:+.4e} deg/day, '
          f'forced-sec amplitude bound = {res["forced_secular_amplitude_bound_deg_day"]:+.4e}, '
          f'ratio = {res["ratio_bound_to_standard"]:.4f}')

print()
print('=== Octupole diagnostic ===')
for i_deg in [97.7876, 90.0, 30.0]:
    res = octupole_lunisolar_raan_rate_rad_s(600.0, i_deg)
    print(f'i={i_deg:6.2f} deg: octupole lunar = {res["octupole_lunar_deg_day"]:+.4e} deg/day')
