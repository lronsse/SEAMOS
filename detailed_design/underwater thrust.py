import numpy as np

# required speed is 2.5km/h which is 0.7 m/s
V_current = 1.3 # [ m / s ] --> this gives room for moving up current for current speeds up to 1.3 m/s
V_move = 0.7
V_flow = V_current + V_move
l = 2 # [ m ] --> bit random estimated diameter of the ball ( revise )
nu = 1.3 * 10 ** ( -6 ) # kinematic viscosity of water at a tempersture om 10^C
Re = ( V_current * l ) / nu # calculate reynoldsnumber, used to pick a Cd
Cd = 0.3 # The chosen Cd
ozin_to_Nm = 0.0070615 # convert random internet unit
Q = 180 * ozin_to_Nm # the conversion
m = 16 # [ kg ] estimated mass of the drone
a = 0.1 # acceleration required to reach speed of 1 [ m/s ] in 10 seconds

def thrust_required(Cd, V, D):
    rho = 1023
    r = D / 2
    S = np.pi * r ** 2
    T = Cd * ( 1/ 2 ) * rho * V ** 2 * S + m * a # just F = m*a with thrust and drag
    return T

def power_required(Cd, V, D):
    T = thrust_required(Cd, V, D) # power is T * V
    P_req = T * V
    return P_req

def rpm_from_thrust(T, V, Q, e_prop):
    rpm = ( T * V ) / ( Q * e_prop * 2 * np.pi * 60 )
    return rpm

def eta_prop(T, V, rpm, Q):
    e_p = ( T * V ) / ( 2 * np.pi * rpm / 60 * Q )
    return e_p

def eta_total(T, V, volt, I):
    e_tot = ( T * V ) / ( volt * I )
    return e_tot



def thrust_from_motor(eta_p, eta_m, volt, I, V):
    thrust_delivered = ( eta_p * eta_m * volt * I ) / V
    return thrust_delivered

def diameter_calc(Q, omega, K_q):
    rho = 1023
    D  = ( ( 4 * np.pi ** 2 * Q ) / ( rho * omega ** 2 * K_q) ) ** ( 1 / 5 )
    return D

def advance_coef(V_a, omega, D):
    J = ( 2 * np.pi * V_a ) / ( omega * D )
    return J

def thrust_coef():
    v
    return

# now using: Q = 0.19, rpm = 16,000

T_req = thrust_required(Cd, V_flow, 0.5)
P_req = power_required(Cd, V_flow, 0.5)
rpm_req = rpm_from_thrust(T_req, V_flow, 2.39, 0.7 )
T_del = thrust_from_motor(0.7, 0.7, 220, 4.29, V_flow)
P_i = eta_total(T_req, V_flow, 24, 5.1) * P_req
e_p = eta_prop(T_req, V_flow, 16000, 0.19 )



print('thrust required:', T_req)
print('Power required:', P_req)
print('thrust delivered:', T_del)
print('input power required:', P_i)
print('propeller efficiency:', e_p)
print()


