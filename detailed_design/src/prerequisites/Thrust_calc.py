import numpy as np
import math
import matplotlib.pyplot as plt
#from isa import


def Advance_ratio(V_airspeed,omega_rotor,D_prop):
    J = 2*math.pi*V_airspeed/(omega_rotor*D_prop)
    return J


def Torque_prop(rho,D_rotor,C_Q_0,C_Q_1,J,omega_rotor):
    torque = rho*D_rotor**5*(C_Q_0+C_Q_1*J)*omega_rotor**2
    return torque


def Thrust_prop(rho):
    C_T_0 = 0.126
    C_T_1 = -0.1378
    omega_rotor = 800
    D_rotor = 0.356
    J=Advance_ratio(18.5,omega_rotor,D_rotor)
    thrust = rho * D_rotor ** 4 * (C_T_0 + C_T_1 * J) * omega_rotor ** 2/(4*np.pi)
    return thrust


def Thrust_BeardMclain(rho,S_swept_prop,eta_rotor,v_exit,V_airspeed):
    thrust = 0.5*rho*S_swept_prop*eta_rotor*((v_exit)**2-V_airspeed**2)
    return thrust


def Thrust_Fitzpatrick(rho,S_swept_prop,eta_rotor,k_motor,delta_t,V_airspeed):
    thrust = rho * S_swept_prop*eta_rotor*(V_airspeed+delta_t*(k_motor- V_airspeed))*delta_t*(k_motor- V_airspeed)
    return thrust


def Thrust():
    V_airspeed = 70  # m/s
    rho = 0.8194  # kg/m3
    D_rotor = 0.4  # m
    omega_rotor = 1  #
    J = Advance_ratio(V_airspeed, omega_rotor, D_rotor)
    C_Q_0 = 1
    C_Q_1 = 1
    C_T_0 = 0.126
    C_T_1 = -0.01378
    eta_rotor = 0.8
    k_motor = 1
    delta_t = 1
    S_swept_prop=np.pi*D_rotor**2/4
    print(rho*S_swept_prop*eta_rotor*V_airspeed)
    V_exit = k_motor*delta_t
    thrust_prop=Thrust_prop(rho,D_rotor,C_T_0,C_T_1,J,omega_rotor)
    thrust_fritz=Thrust_Fitzpatrick(rho,S_swept_prop,eta_rotor,k_motor,delta_t,V_airspeed)
    thrust_BM = Thrust_BeardMclain(rho, S_swept_prop , eta_rotor, k_motor, delta_t, V_airspeed)
    return thrust_prop,thrust_fritz,thrust_BM



## Assumed Thrust 8 N
V_airspeed = 20  # m/s
rho = 0.8194  # kg/m3
D_rotor = 0.4  # m
omega_rotor = 1  #
J = Advance_ratio(V_airspeed, omega_rotor, D_rotor)
C_Q_0 = 1
C_Q_1 = 1
C_T_0 = 0.126
C_T_1 = -0.01378
eta_rotor = 0.8
k_motor = 1
delta_t = 1
S_swept_prop=np.pi*D_rotor**2/4
thrust=20

def Propeller_from_Thrust(thrust,v_airspeed,n_propellers,eta_rotor,v_prop_out):
    #omega_rotor=2*np.pi*rpm/60
    thrust=thrust/n_propellers
    #S_swept_prop = np.pi * D_rotor ** 2 / 4
    s_swept_prop=thrust/(0.5*rho*eta_rotor*((v_prop_out)**2-v_airspeed**2))
    d_rotor=np.sqrt(s_swept_prop*4/np.pi)
    # thrust = 0.5*rho*S_swept_prop*eta_rotor*((v_prop_out)**2-V_airspeed**2)


    return d_rotor

def Thrust_in_climb(climb_angle,drag,mass):
    thrust=drag+mass*9.81*np.sin(climb_angle)
    return thrust

def Pitch_calc(rpm,v):
    return 60/rpm*v
roc=1
v_airspeed=20
climb_angle=np.arcsin(float(roc/v_airspeed)) ##degree
print('V airspeed [m/s] =',v_airspeed,'Rate of climb[m/s] =',roc)
print('Thrust in climb [N]',Thrust_in_climb(climb_angle,10,15))
print('ROC [m/s]',v_airspeed*np.sin(climb_angle))
print('Climb angle radians:', climb_angle)

#print(Thrust_in_climb(climb_angle,10,15)/9.81)
#print('Propeller diameter from thrust[m]', Propeller_from_Thrust(Thrust_in_climb(climb_angle,10,15),v_airspeed,1,0.9,v_exit))
print('thrust using their coefficients ', Thrust_prop(1.225))
pitch=Pitch_calc(7000,20)
print('Pitch',pitch)
D_tab=[]
v_exit_tab=[]
rpm_tab=np.linspace(5000,8000,300)
for rpm in rpm_tab:
    D = 0.1
    for i in range(100):
        v_exit=np.sqrt(v_airspeed**2+rpm*2*np.pi/60*D)
        D =Propeller_from_Thrust(Thrust_in_climb(climb_angle,10,15),v_airspeed,1,0.62,v_exit)
    D_tab=np.append(D_tab,D)
plt.plot(rpm_tab,D_tab)
plt.show()
print(D)
print('pitch to diameter',pitch/D)
print(Thrust_prop(1.225)*20*60000/20/3600)
print(thrust*20*60000/20/3600)
#Thrust = rho*eta*v*S
#S=D2/4*pi
#print(np.sqrt(5*4/(np.pi*rho*eta_rotor*V_airspeed)))
print(60000/20/3600)
print(D)
def Thrust_in_climb_average(drag,mass):
    thrust=drag+mass*9.81*np.sin(np.linspace(np.pi/2,0,1000))
    thrust=np.average(thrust)
    return thrust
def Energy_hop_flight(v_cruise,cruise_alt,drag_climb):
    energy_climb=[]
    energy_cruise=[]

    roc_tab=np.linspace(0.01,10,1000)
    roc_tab_2=[]
    climb_distance_tab=[]
    for roc in roc_tab:
        climb_time=cruise_alt/roc
        climb_distance=(climb_time*np.sqrt(v_cruise**2-roc**2))

        hop_distance = np.sqrt(100 ** 2 + 100 ** 2)
        if climb_distance<hop_distance:

        #print(climb_time,climb_distance)
            energy_climb=np.append(energy_climb,climb_time*drag_climb*v_cruise)

            energy_cruise =np.append(energy_cruise, (hop_distance-climb_distance)*10)
            roc_tab_2 = np.append(roc_tab_2, roc)
            climb_distance_tab = np.append(climb_distance_tab, climb_distance)

    print('energy_climb[J]',energy_climb)
    print('energy_cruise[J]', energy_cruise)
    print(len(energy_climb),len(energy_cruise))
    plt.plot(roc_tab_2,energy_cruise,label='energy_cruise')
    plt.plot(roc_tab_2, energy_climb,label='energy_climb')
    energy_total_hop=energy_cruise + energy_climb
    print(print('energy_total_hop[J]', energy_total_hop[-1]))
    plt.plot(roc_tab_2, energy_total_hop,label='total energy per hop')
    plt.legend(loc='upper right')

    plt.show()
    #plt.plot(roc_tab_2,climb_distance_tab)
    #plt.show()
    return energy_climb,energy_cruise
Energy_hop_flight(20,10,Thrust_in_climb(climb_angle,10,15))



