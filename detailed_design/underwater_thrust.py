import matplotlib.pyplot as plt
import numpy as np

# required speed is 2.5km/h which is 0.7 m/s
V_current = 1 # [ m / s ] --> this gives room for moving up current for current speeds up to 1.3 m/s
V_move = 0.5
V_flow = 3 #V_current + V_move
l = 2 # [ m ] --> bit random estimated diameter of the ball ( revise )
nu = 1.3 * 10 ** ( -6 ) # kinematic viscosity of water at a tempersture om 10^C
Re = ( V_current * l ) / nu # calculate reynoldsnumber, used to pick a Cd
Cd = 0.3 # The chosen Cd
ozin_to_Nm = 0.0070615 # convert random internet unit
Q = 180 * ozin_to_Nm # the conversion
m = 16 # [ kg ] estimated mass of the drone
accel = 0.1 # acceleration required to reach speed of 1 [ m/s ] in 10 seconds
diam = 0.4
n_blades=1
# diam = np.linspace(0.05, 0.4, 100)
a=1.5
rho=1023
print(Re)


# diam = 0.2 #np.linspace(0.05, 0.4, 100)

def thrust_required(Cd, V, D):
    rho = 1023
    r = D / 2
    S = np.pi * r ** 2
    #T = Cd * ( 1/ 2 ) * rho * V ** 2 * S + m * accel # just F = m*a with thrust and drag
    T = 70 + m * accel
    return T

def power_required(Cd, V, D):
    T = thrust_required(Cd, V, D) # power is T * V
    P_req = T * V
    return P_req

def rpm_from_thrust(T, V, Q, e_prop):
    rpm = (( T * V ) / ( Q * e_prop * 2 * np.pi * 60 ) ) / n_blades
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

def rpm_from_prop(V, D, phi, alpha,n_blades):
    phi = (phi / 360 * (2 * np.pi))
    alpha = (alpha / 360 * (2 * np.pi))
    rpm =( ( 2 * V ) / ( D * np.tan( phi - alpha ) )) / n_blades
    return rpm

def find_torque(phi, alpha, T, D,a):
    phi = (phi / 360 * (2 * np.pi))
    alpha = (alpha / 360 * (2 * np.pi))
    a = 1.5
    Q = ( a * np.tan(phi - alpha) / 4 ) * T * D
    return Q

# print("ahhh",find_torque(phi,alpha,T,D,a))

def torque_from_rpm(rpm, Kq, D):
    omega = ( rpm * 2 * np.pi ) / 360
    rho = 1023
    Q = ( Kq * rho * omega ** 2 * D ** 5 ) / ( 4 * np.pi ** 2 )
    return Q


def driving_pwr(T, D,a,rho):
    rho = 1023
    a = 1.5
    P = ( a / ( rho * np.pi )) ** ( 1 / 2 ) * ( T ** (3 / 2) ) / D
    return P

def battery_weight(P, V_ground):
    P = P / 0.57 / 0.8
    t = 100 / V_ground
    n_hours = t / 3600
    Wh = n_hours * P
    m = Wh * ( 1 / 200 )
    return m, Wh



# now using: Q = 0.19, rpm = 16,000
# phi = 10degr, alpha = 5 degr

T_req = thrust_required(Cd, V_flow, 0.4) # / 2 # uncomment when thinking of 2 props
P_req = power_required(Cd, V_flow, 0.4)
rpm_req = rpm_from_thrust(T_req, V_flow, 2.39, 0.7 )
T_del = thrust_from_motor(0.7, 0.7, 220, 4.29, V_flow)
P_i = eta_total(T_req, V_flow, 24, 5.1) * P_req
e_p = eta_prop(T_req, V_flow, 16000, 0.19 )

#print('thrust delivered:', T_del)
#print('input power required:', P_i)
#print('propeller efficiency:', e_p)
#print(P_req / e_p)

print('thrust required:', T_req)
print('Power required:', P_req)
print('rpm:',rpm_from_prop(V_flow, diam, 10,  5 ,n_blades) * 60 / ( 2 * np.pi))
print('applied torque:', find_torque(10, 5, T_req , diam,a))
print('driving power:', driving_pwr(T_req , diam,a,rho))

# plt.plot(diam, rpm_from_prop(V_flow, diam, 10,  5 ,n_blades) * 60 / ( 2 * np.pi))
# plt.title('Rotations per minute vs. Propeller diameter')
# plt.xlabel('Propeller diameter [m]')
# plt.ylabel('Angular speed [rpm]')
# #plt.show()
# plt.plot(diam,find_torque(10, 5, T_req , diam,a) )
# plt.title('Shaft torque vs. Propeller diameter')
# plt.xlabel('Propeller diameter [m]')
# plt.ylabel('Shaft torque [Nm] ')
# #plt.show()
# plt.plot(diam, driving_pwr(T_req , diam,a,rho) )
# plt.title('Propeller driving power vs. Propeller diameter')
# plt.xlabel('Propeller diameter [m]')
# plt.ylabel('Propeller driving power required [W]')
# #plt.show()

print('torque:', torque_from_rpm(800, 0.0877, 0.2))


######### Verification ########

def test_thrust_required():
    Cd=1
    V=2
    D=3
    r= thrust_required(Cd, V, D)
    S= thrust_required(Cd, V, D)
    T= thrust_required(Cd, V, D)
    r_test = 1.5
    S_test = 7.06858347058
    T_test = 14463.9217808
    assert r==r_test
    assert S==S_test
    assert T==T_test

def test_power_required():
    Cd=1
    V=2
    D=3
    T = thrust_required(Cd, V, D) # power is T * V
    P_req_test = 28927.8435616
    P_req=power_required(Cd,V,D)
    assert abs(P_req_test-P_req) <1e-6

def test_rpm_from_prop():
    V=1
    D=2
    phi=3
    alpha=4
    n_blades=1
    phi = rpm_from_prop(V,D,phi,alpha,n_blades)
    alpha = rpm_from_prop(V,D,phi,alpha,n_blades)
    rpm=rpm_from_prop(V,D,phi,alpha,n_blades)
    phi_test = 0.0523598775598
    alpha_test = 0.0698131700798
    rpm_test = -0.642092615934
    assert abs(phi_test-phi) < 1e-6
    assert abs(alpha_test - alpha) < 1e-6
    assert abs(rpm_test - rpm) < 1e-6

def test_find_torque():
    phi_t=1
    alpha_t=2
    T_t=3
    D_t=4
    a_T = 5
    Q=find_torque(phi_t,alpha_t,T_t,D_t,a_T)
    phi=find_torque(phi_t,alpha_t,T_t,D_t,a_T)
    alpha =find_torque(phi_t,alpha_t,T_t,D_t,a_T)
    test_phi = 0.0174532925199
    test_alpha = 0.0349065850399
    test_Q = -23.3611158698
    assert abs(phi-test_phi) < 1e-6
    assert abs(alpha-test_alpha) < 1e-6
    assert abs(Q-test_Q) < 1e-6

def test_driving_pwr():
    rho = 1023
    a = 1.5
    T=2
    D=3
    P=driving_pwr(T,D,a,rho)
    test_P = 0.02036838592618265
    assert abs(P-test_P) < 1e-6