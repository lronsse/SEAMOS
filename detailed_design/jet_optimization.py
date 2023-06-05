import scipy.optimize as optimize
import math
import matplotlib.pyplot as plt

# Material properties
rho = 2710  # [kg/m3]  Aluminum
sigma_yield = 240E6  # [Pa]  Aluminum

# UAV characteristics
m_tot = 16  # [kg]

# Mission characteristics
rho_w = 1000  # [kg/m3]
dV = 1
pressure_x0 = 10E5  # [Pa]  Initial guess for optimization


def m_jet(P):
    vol_w = (math.exp(dV / math.sqrt(2 * P / rho_w)) - 1) * m_tot / rho_w
    r = math.pow(3/4 * vol_w / math.pi, (1/3))
    t = P * r / (2 * sigma_yield)
    return 4 * math.pi * r ** 2 * t * rho

results = []
V = []
mass = []
for i in range(20):
    res = optimize.minimize_scalar(m_jet, pressure_x0, bounds=(0., 200E5), method='bounded')
    mass.append(m_jet(res.x))
    results.append(res.x)
    V.append(dV)
    dV += 1

plt.plot(V, results)
plt.show()
plt.plot(V, mass)
plt.show()