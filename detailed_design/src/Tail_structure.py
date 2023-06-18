import numpy as np

### tail planform
C_root = 0.2581
C_tip = 0.5 * C_root
b = 0.3
S = ( ( C_root + C_tip ) / 2 ) * b
t_tail = 0.09 * C_tip
Cl = 0.7 # naca0009 coefficient at cruise AoA of 6^
t_skin = 0.001

### tail cross section
# Ixx = Ixx + A * d^2
# Ixx rectangle = ( b * h^3 ) / 12

# ==> to be conservative: cross section at tip
Ixx = ( ( 1 / 12 ) * C_tip * t_skin ** 3 + ( C_tip * t_skin ) * ( t_tail / 2 ) ** 2 ) \
    + ( 1 / 12 ) * t_skin * t_tail ** 3 * 2

Izz = 2 * ( 1 / 12 ) * C_tip ** 3 * t_skin + 2 * ( ( 1 / 12 ) * t_tail * t_skin **  3 + t_skin * t_tail * ( C_tip / 2 ) ** 2)


### aerodynamic loads

def lift(V, S, Cl):
    L = ( 1 / 2 ) * 1.225 * V ** 2 * S * Cl * np.cos((np.pi / 6)) ** 2
    return L

def stress_lift(L, b, I, y):
    M = L * ( b / 2)
    sigma = ( M * y ) / I
    sigma = sigma / 1000000
    return sigma

def shear_lift(V, t, h, Ixx):
    q = ( V * t * h ** 2 ) / ( 8 * Ixx )
    tau = q / t
    tau = tau / 1000000
    return tau

L = lift(20, S, Cl)
sigma_ad = stress_lift(L, b, Ixx, (t_tail/2))
tau_ad = shear_lift(L, t_skin, t_tail, Ixx)
print('lift:', L)
print('highest normal stress in MPa:', sigma_ad)
print('highest shear stress in MPa:',tau_ad )
### launch loads --> for fuselage

### impact forces when diving
# --> at tail level around 800 [N]
# LEADING EDGE SURFACE
F = 800 / 3
t_tip = 0.09 * C_tip
t_root = 0.09 * C_root
S_le = ( t_tip + t_root ) / 2 * b
w = F / b

sigma_impact = stress_lift(F, b, Izz, (C_tip/2))
tau_impact = shear_lift(F, t_skin, C_tip, Izz)

print('normal stress due to impact:', sigma_impact)
print('shear stress due to impact:',tau_impact)

