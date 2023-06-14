from wing_aero import max_moment
import structural_sizing as ss
import numpy as np


print(np.max(abs(np.array(ss.wing.normal_stress_bending))))

print(np.max(abs(np.array(ss.wing.shear_stress_torsion(max_moment)))))
print(np.max(abs(np.array(ss.wing.calculate_shear_stress(ss.wing.second_moment_of_area, ss.wing.first_moment_of_area)))))


