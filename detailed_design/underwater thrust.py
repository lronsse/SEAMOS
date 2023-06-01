import numpy as np

V_current = 3.5 # [ m / s ]
l = 2 # [ m ]
nu = 1.3 * 10 ** ( -6 )
Re = ( V_current * l ) / nu
Cd = 0.3


def thrust_required(Cd, V, D):
    rho = 998
    r = D / 2
    S = np.pi * r ** 2
    T = Cd * ( 1/ 2 ) * rho * V ** 2 * S
    return T

def power_required():
    P_req = T * V
    return P_req

print(thrust_required(Cd, V_current, 0.5))

