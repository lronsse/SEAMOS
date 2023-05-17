import numpy as np
#from sympy import symbols, Eq, solve
import sympy as sm

###
D_ts = 0.2064
D_fl = 0.2184
L_ts = 5.5 * D_ts
L_fl = 5.5 * D_fl
Cd_fl = 0.85
Cd_ts = 1.2
V_current = 0.65
rho=1027
V=0.5
A_ts = 0.6 * D_ts * L_ts
A_fl = np.pi * ( D_fl / 2 ) ** 2 * 0.3

def drag(rho, V, A , Cd):
    D = ( 1 / 2 ) * rho * V ** 2 * A * Cd
    return D

def Vrel(D, rho, S, Cd ):
    V = np.sqrt ( D / ( ( 1 / 2 ) * rho * S * Cd ) )
    return V

### FOR TAILSITTER

#Cd for wing perpendicular to wind = 1.28
#relative velocity wind - current

Drag_ts_wind = drag(1.225, 3, 10, 1.28)
print('wind drag:', Drag_ts_wind)
A_under = 1.5 * 0.5
V_relative = Vrel(Drag_ts_wind, 1028, A_under, 1.2 )
print('Drifting velocity for tailsitter, wind & stream aligned:', V_relative + V_current )
print('Drifting velocity for tailsitter, wind & stream opposite', )


Cd_a = 1.28
Cd_w = 1.2
V_a =  5.6
V_w = 0.65
r_a = 1.225
r_w = 1028
S_a = 6
S_w = 0.5 * 1.5


p =[- r_a * Cd_a * S_a + S_w *  r_w * Cd_w , \
    - 2 * V_w * r_w * Cd_w * S_w + 2 * S_a * V_a * r_a * Cd_a, \
    - S_a * r_a * Cd_a * V_a**2 + S_w*  r_w * Cd_w * V_w ** 2 ]

print(np.roots(p))
