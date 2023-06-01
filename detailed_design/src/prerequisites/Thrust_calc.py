import numpy as np
import math
#from isa import


def Advance_ratio(V_airspeed,omega_rotor,D_rotor):
    J = 2*math.pi*V_airspeed/(omega_rotor*D_rotor)
    return J


def Torque_prop(rho,D_rotor,C_Q_0,C_Q_1,J,omega_rotor):
    torque = rho*D_rotor**5*(C_Q_0+C_Q_1*J)*omega_rotor**2
    return torque


def Thrust_prop(rho,D_rotor,C_T_0,C_T_1,J,omega_rotor):
    thrust = rho * D_rotor ** 4 * (C_T_0 + C_T_1 * J) * omega_rotor ** 2
    return thrust


def Thrust_BeardMclain(rho,S_swept_prop,eta_rotor,k_motor,delta_t,V_airspeed):
    thrust = 0.5*rho*S_swept_prop*eta_rotor*((k_motor*delta_t)**2-V_airspeed**2)
    return thrust


def Thrust_Fitzpatrick(rho,S_swept_prop,eta_rotor,k_motor,delta_t,V_airspeed):
    thrust = rho * S_swept_prop*eta_rotor*(V_airspeed+delta_t*(k_motor- V_airspeed))*delta_t*(k_motor- V_airspeed)
    return thrust


def Thrust():
    V_airspeed = 70  # m/s
    rho = 1.25  # kg/m3
    D_rotor = 0.2  # m
    omega_rotor = 1  #
    J = Advance_ratio(V_airspeed, omega_rotor, D_rotor)
    C_Q_0 = 1
    C_Q_1 = 1
    C_T_0 = 0.126
    C_T_1 = -0.01378
    eta_rotor = 1
    k_motor = 1
    delta_t = 1
    S_swept_prop=np.pi*D_rotor**2/4
    V_exit = k_motor*delta_t
    thrust_prop=Thrust_prop(rho,D_rotor,C_T_0,C_T_1,J,omega_rotor)
    thrust_fritz=Thrust_Fitzpatrick(rho,S_swept_prop,eta_rotor,k_motor,delta_t,V_airspeed)
    thrust_BM = Thrust_BeardMclain(rho, S_swept_prop , eta_rotor, k_motor, delta_t, V_airspeed)
    return thrust_prop,thrust_fritz,thrust_BM

## Assumed Thrust