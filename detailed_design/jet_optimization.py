import scipy.optimize as optimize
import math
import matplotlib.pyplot as plt
import numpy as np

# Material properties
rho = 2710  # [kg/m3]  Aluminum
sigma_yield = 200E6  # [Pa]  Aluminum

# UAV characteristics
m_tot = 16  # [kg]

# Mission characteristics
rho_w = 1000  # [kg/m3]
dV = 1 # start!
pressure_x0 = 10E5  # [Pa]  Initial guess for optimization

# extra bonus
t_min = 1E-3 # [m]

# def vol_w_implicit(v):
#     return np.abs((math.exp(dV / math.sqrt(2 * P / rho_w)) - 1) * (m_tot + 6 * v * P * rho / sigma_yield) / rho_w - v)

def m_jet(P, verbose=False):
    # def vol_w_implicit(v):
    #     return np.abs(
    #         (math.exp(dV / math.sqrt(2 * P / rho_w)) - 1) * (m_tot + 6 * v * P * rho / sigma_yield) / rho_w - v)
    # vol_w = optimize.minimize_scalar(vol_w_implicit, bounds=[0., 100], method='bounded').x
    vol_w = (math.exp(dV / math.sqrt(2 * P / rho_w)) - 1) * (m_tot+6*P*rho/sigma_yield) / rho_w
    print(P/1e5, vol_w*1000)

    print(dV)
    r = (3/4*1/math.pi*vol_w)**(1/3)
    t = max((P * r / (2 * sigma_yield)), t_min)
    return 4 * math.pi * r ** 2 * t * rho
    # return 6*vol_w*P*rho/sigma_yield

def dv(P, v):
    return np.sqrt(2*P/rho_w)*np.log(1+(rho_w*v)/(m_tot+6*v*P*rho/sigma_yield))

print(dv(100E5, 5E-3))

# res = []
# pr = 1e5*np.logspace(1, 5, 200)
# for p in pr:
#     res.append(dv(p, .150/1000))
#
# plt.plot(pr, res)
# plt.show()


# results = []
# V = []
# mass = []
# for i in range(20):
#     res = optimize.minimize_scalar(m_jet, pressure_x0, bounds=(0., 200E5), method='bounded')
#     mass.append(m_jet(res.x, verbose=True))
#     results.append(res.x)
#     V.append(dV)
#     dV += 1
#
# results = list(entry/1e5 for entry in results)
#
# plt.plot(V, results)
# plt.show()
# plt.plot(V, mass)
# plt.show()