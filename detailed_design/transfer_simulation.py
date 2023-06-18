import matplotlib.pyplot as plt
import numpy as np

# Material properties
rho = 2710  # [kg/m3]  Aluminum
sigma_yield = 240E6  # [Pa]  Aluminum

# UAV characteristics
m_tot = 16  # [kg]
length_uav = 2
nu = 1.3 * 10 ** (-6)

# Mission characteristics
rho_w = 1000  # [kg/m3]
target_v = 12  # [m/s]
target_h = 15  # [m]
x_horizontal = 100  # [m]
x_vertical = 100  # [m]
x_direct = np.sqrt(x_horizontal ** 2 + x_vertical ** 2)  # hypothenuse distance of mission path
v_current = 1.3  # [m/s]



