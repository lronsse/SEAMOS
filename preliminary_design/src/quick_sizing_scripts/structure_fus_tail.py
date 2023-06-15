import numpy as np
import matplotlib.pyplot as plt

class Material:
    def __init__(self, rho, E, sigma_yield, tau_yield):
        self.rho = rho
        self.E = E
        self.sigma_yield =sigma_yield
        self.tau_yield = tau_yield



def get_points_on_circle(number_of_points, radius):
    points = []
    for i in range(number_of_points):
        angle = np.pi * 2 * i / number_of_points
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        points.append((x, y))
    return points

def calc_first_moment_of_area(locations, area):
    Q = 0
    for i in range(n_booms):
        Q += abs(locations[i][1]) * area
    return Q

def calc_second_moment_of_area(locations, area, second_moment_boom):
    I = n_booms * second_moment_boom
    for i in range(n_booms):
        I += locations[i][1] ** 2 * area
    return I


n_booms = 3
radius_boom = 4 * 10 ** -3

area_booms = np.pi * radius_boom ** 2
mass_puffin = 16
cfrp = Material(1800, 180*10**9, 100*10**6, 50000000)

radius_fuselage = 0.03
thickness_fuselage = 1 * 10 ** -3
length_fuselage = 1
boom_locations = get_points_on_circle(n_booms, radius_fuselage)
first_moment_of_area = calc_first_moment_of_area(boom_locations, area_booms)
second_moment_of_area_boom = 0.25 * np.pi * radius_boom ** 4

second_moment_of_area = n_booms * second_moment_of_area_boom
for i in range(n_booms):
    second_moment_of_area += boom_locations[i][1] ** 2 * area_booms


alpha = np.radians(23.4) * 2 / (1.5 ** 2)
test_I = mass_puffin * (radius_fuselage * 10) ** 2

normal_force = 50
bending_moment = mass_puffin * 9.81
torque = test_I * alpha
shear_force = mass_puffin * 9.81

area = n_booms * area_booms

polar_moment_of_inertia = n_booms * radius_fuselage * area_booms


normal_stress_normal_force = normal_force / area
normal_stress_bending = bending_moment * radius_fuselage / second_moment_of_area
shear_stress_torsion = torque * radius_fuselage / polar_moment_of_inertia
shear_stress_shear_force = shear_force * first_moment_of_area / (second_moment_of_area * thickness_fuselage)

max_force_buckling = (np.pi) ** 2 * cfrp.E * second_moment_of_area / (length_fuselage ** 2)
mass_structure = area * length_fuselage * cfrp.rho


print(f'Maximal normal stress = {(abs(normal_stress_normal_force) + abs(normal_stress_bending))*10**-6 * 1.25} MPa')
print(f'Maximal shear stress = {(abs(shear_stress_shear_force) + abs(shear_stress_torsion)) * 10 **-6 * 1.25} MPa')
print(f'Max force for buckling = {max_force_buckling * 0.75} N')
print(f'Mass of the {n_booms} sticks = {mass_structure} kg')
if max_force_buckling <= normal_force:
    print('normal force is too high')




