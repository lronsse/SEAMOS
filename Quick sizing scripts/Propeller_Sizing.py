import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

Treq = 30
rho = 0.8
v0 = 20
#omega_rpm = 1

def func(r):
    return rho * np.pi * r ** 2 * (v0 + 0.5 * np.sqrt(v0 ** 2 + (omega * r) ** 2)) * np.sqrt(v0 ** 2 + (omega * r) ** 2) - Treq


#res = op.root(func, 1)

#print(res)
sol = []
rpm = []
for i in range(10000):
    rpm.append(i)
    omega_rpm = i
    omega = omega_rpm * 2 * np.pi / 60
    sol.append(op.root(func, 1).x)


plt.plot(sol, rpm)
plt.show()

