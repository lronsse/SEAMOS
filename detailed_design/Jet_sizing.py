import numpy as np
import matplotlib.pyplot as plt
import math


# UAV data
length_uav = 2  # [m]
diameter_uav = 0.2  # [m]
mass_uav = 16  # [kg]
Cd = 0.3

# path calculations
min_depth = -10  # [m]
jet_depth = -2  # [m]
buoyancy_accel = 0.3  # [m/s^2]
buoyancy_thrust = buoyancy_accel * mass_uav

## Water parameters
nu = 1.3 * 10 ** (-6)


def drag(Cd, V, Diam_uav, len_uav):
    Re = (V * length_uav) / nu
    rho = 998
    radius = Diam_uav / 2
    Area = np.pi * radius ** 2
    Drag = Cd * (1 / 2) * rho * V ** 2 * Area
    return Drag


V = [0]
h = [min_depth]
t = [0]
dt = 1
i = 0
if __name__ == '__main__':
    while h[i] <= 5:  # maneuver up to 5m in height
        if h[i] <= jet_depth:
            a = (buoyancy_thrust - drag(Cd, V[i], diameter_uav, length_uav)) / mass_uav
            dV = a * dt
        V.append(V[i] + dV)
        h.append(h[i] + V[i+1] * dt)
        t.append(t[i] + dt)
        i += 1
        print(h[i])

