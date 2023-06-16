import numpy as np
from scipy.optimize import root

force_targets = np.array([50, 0, 0, 0, 5, 0])
fin_angle = np.deg2rad(30)
neutral_angle = np.deg2rad(23.1985905136)

x_tm = 0.7
w_tm = 0.3


def force_function(inputs, targets):
    thrust_1 = inputs[0]
    thrust_2 = inputs[1]
    thrust_3 = inputs[2]
    psi_1 = inputs[3]
    psi_2 = inputs[4]
    psi_3 = inputs[5]
    thrust_x = thrust_1 * np.cos(psi_1) + thrust_2 * np.cos(psi_2) + thrust_3 * np.cos(psi_3)
    thrust_y = - thrust_2 * np.cos(fin_angle) * np.sin(psi_2) + thrust_3 * np.cos(fin_angle) * np.sin(psi_3)
    thrust_z = - thrust_1 * np.sin(psi_1) + thrust_2 * np.sin(fin_angle) * np.sin(psi_2) + thrust_3 * np.sin(fin_angle) * np.sin(psi_3)
    moment_x = 0
    moment_y = thrust_1 * (-w_tm * np.cos(psi_1) + x_tm * np.sin(psi_1)) + \
        thrust_2 * (w_tm * np.cos(psi_2) * np.sin(fin_angle) + x_tm * np.sin(psi_2) * np.sin(fin_angle)) + \
        thrust_3 * (w_tm * np.cos(psi_3) * np.sin(fin_angle) + x_tm * np.sin(psi_3) * np.sin(fin_angle))
    moment_z = thrust_2 * (w_tm * np.cos(psi_2) * np.cos(fin_angle) + x_tm * np.sin(psi_2) * np.cos(fin_angle)) + \
        thrust_3 * (w_tm * np.cos(psi_3) * np.cos(fin_angle) + x_tm * np.sin(psi_3) * np.cos(fin_angle))
    return np.array([thrust_x, thrust_y, thrust_z, moment_x, moment_y, moment_z]) - targets


solution = root(force_function, np.array([0, 0, 0, neutral_angle, neutral_angle, neutral_angle]), force_targets, method='lm').x

print(f"Targets: {force_targets}")
print(f"\nThrust XYZ: {solution[0]:3.2f}, {solution[1]:3.2f}, {solution[2]:3.2f} \nAngle XYZ: {np.rad2deg(solution[3]):3.2f}, {np.rad2deg(solution[4]):3.2f}, {np.rad2deg(solution[5]):3.2f}")
forces = force_function(solution, np.zeros(6))
print(f"\nResults: {forces[0]:3.2f}, {forces[1]:3.2f}, {forces[2]:3.2f} \nMoment XYZ: {forces[3]:3.2f}, {forces[4]:3.2f}, {forces[5]:3.2f}")
