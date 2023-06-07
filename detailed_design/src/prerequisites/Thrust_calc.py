import numpy as np
import math
import matplotlib.pyplot as plt
#from isa import

"""
(written by Mikolaj)
Propeller sizing based on thrust requirement
Energy consumption for the hop(initial stall velocity, no initial velocity)
Assumed Thrust 11 N
"""

def Advance_ratio(v_airspeed,omega_rotor,D_prop):
    J = 2*math.pi*v_airspeed/(omega_rotor*D_prop)
    return J


def Propeller_from_Thrust(thrust,v_airspeed,n_propellers,eta_rotor,v_prop_out):
    """
    Propeller diameter based on thrust needed, speed and efficiency, rho should be revisited
    """
    #omega_rotor=2*np.pi*rpm/60
    rho=1.225
    thrust=thrust/n_propellers
    #S_swept_prop = np.pi * D_rotor ** 2 / 4
    s_swept_prop=thrust/(0.5*rho*eta_rotor*((v_prop_out)**2-v_airspeed**2))
    d_rotor=np.sqrt(s_swept_prop*4/np.pi)
    # thrust = 0.5*rho*S_swept_prop*eta_rotor*((v_prop_out)**2-V_airspeed**2)
    return d_rotor
#print(Propeller_from_Thrust(15*9.81,12,1,0.9,24))


def Thrust_in_climb(climb_angle,drag,mass):
    thrust=drag+mass*9.81*np.sin(climb_angle)
    return thrust


def Pitch_calc(rpm,v):
    return 60/rpm*v


def Prop_calc():
    """
    Just run this for the propeller sizing
    """
    mass=15
    print('Puffin mass')
    drag_cruise=11
    print('Drag in cruise',drag_cruise,'[N]')
    roc=1
    print('Rate of climb', roc, '[m/s]')
    v_airspeed=20
    print('Cruise velocity', v_airspeed, '[m/s]')
    climb_angle = np.arcsin(float(roc / v_airspeed))  ##degree
    print('Thrust in climb [N]', Thrust_in_climb(climb_angle, 10, 15))
    print('Climb angle radians:', climb_angle)
    number_of_prop=1
    print('number of propellers',number_of_prop)
    D_tab=[]
    v_exit_tab=[]
    rpm_tab=np.linspace(400,8000,1000)
    for rpm in rpm_tab:
        D = 0.1
        for i in range(100):
            v_exit=np.sqrt(v_airspeed**2+rpm*2*np.pi/60*D/2)
            D =Propeller_from_Thrust(Thrust_in_climb(climb_angle,drag_cruise,mass),v_airspeed,1,0.9,v_exit)
        v_exit_tab=np.append(v_exit_tab,v_exit)
        D_tab=np.append(D_tab,D)
    plt.plot(rpm_tab,D_tab)
    #plt.plot(rpm_tab, v_exit_tab)
    plt.xlabel('rpm')
    plt.ylabel('diameter [m]')
    plt.title('Diameter vs rpm')
    plt.show()
    print('Propeller diameter',D,'[m]')
    return climb_angle


def Thrust_in_climb_average(drag,mass):
    thrust=drag+mass*9.81*np.sin(np.linspace(np.pi/2,0,1000))
    thrust=np.average(thrust)
    return thrust


def Energy_hop_VTOL(v_cruise,cruise_alt,mass):
    energy_climb=[]
    energy_cruise=[]

    roc_tab=np.linspace(0.01,10,1000)
    roc_tab_2=[]
    climb_distance_tab=[]
    for roc in roc_tab:
        climb_time=cruise_alt/roc
        climb_distance=(climb_time*np.sqrt(v_cruise**2-roc**2))
        drag_climb_tab=[]
        climb_angle=np.pi/2
        velocity=0
        for i in range(100):

            drag_climb_tab=np.append(drag_climb_tab,Thrust_in_climb(climb_angle,10*(velocity/20)**2,15))
            climb_angle = climb_angle - np.pi / 100
            velocity=velocity+.12
        drag_climb=np.average(drag_climb_tab)
        hop_distance = np.sqrt(100 ** 2 + 100 ** 2)
        if climb_distance<hop_distance:

        #print(climb_time,climb_distance)
            energy_climb=np.append(energy_climb,climb_time*(drag_climb+mass*v_cruise/climb_time)*v_cruise)

            energy_cruise =np.append(energy_cruise, (hop_distance-climb_distance)*10)
            roc_tab_2 = np.append(roc_tab_2, roc)
            climb_distance_tab = np.append(climb_distance_tab, climb_distance)

    #print('energy_climb[J]',energy_climb)
    #print('energy_cruise[J]', energy_cruise)
    print(len(energy_climb),len(energy_cruise))
    plt.plot(roc_tab_2,energy_cruise,label='energy_cruise')
    plt.plot(roc_tab_2, energy_climb,label='energy_climb')
    energy_total_hop=energy_cruise + energy_climb
    print(print('energy_total_hop VTOL[J]', energy_total_hop[-1]))
    plt.plot(roc_tab_2, energy_total_hop,label='total energy per hop')
    plt.legend(loc='lower right')
    plt.xlabel('ROC m/s')
    plt.ylabel('Energy J')
    plt.title('Energy vs roc VTOL')
    plt.show()
    plt.plot(climb_distance_tab,energy_total_hop)
    plt.xlabel('climb distance m')
    plt.ylabel('Energy J')
    plt.title('Energy vs climb distance VTOL')
    plt.show()
    return energy_climb,energy_cruise


def Energy_hop_NOT_VTOL(v_cruise,cruise_alt,mass):
    energy_climb=[]
    energy_cruise=[]

    roc_tab=np.linspace(0.01,10,1000)
    roc_tab_2=[]
    climb_distance_tab=[]
    for roc in roc_tab:
        climb_time=cruise_alt/roc
        climb_distance=(climb_time*np.sqrt(v_cruise**2-roc**2))
        drag_climb_tab = []
        climb_angle = np.pi / 2
        velocity = 0
        for i in range(100):
            drag_climb_tab = np.append(drag_climb_tab, Thrust_in_climb(climb_angle, 10 , 15))
            climb_angle = climb_angle - np.pi / 100
        drag_climb = np.average(drag_climb_tab)
        hop_distance = np.sqrt(100 ** 2 + 100 ** 2)
        if climb_distance<hop_distance:

        #print(climb_time,climb_distance)
            energy_climb=np.append(energy_climb,climb_time*(drag_climb*(12/20)**2)*v_cruise)

            energy_cruise =np.append(energy_cruise, (hop_distance-climb_distance)*10)
            roc_tab_2 = np.append(roc_tab_2, roc)
            climb_distance_tab = np.append(climb_distance_tab, climb_distance)

    #print('energy_climb[J]',energy_climb)
    #print('energy_cruise[J]', energy_cruise)
    print(len(energy_climb),len(energy_cruise))
    plt.plot(roc_tab_2,energy_cruise,label='energy_cruise')
    plt.plot(roc_tab_2, energy_climb,label='energy_climb')
    energy_total_hop=energy_cruise + energy_climb
    print(print('energy_total_hop not vtol[J]', energy_total_hop[-1]))
    plt.plot(roc_tab_2, energy_total_hop,label='total energy per hop')
    plt.legend(loc='lower right')
    plt.xlabel('ROC m/s')
    plt.ylabel('Energy J')
    plt.title('Energy vs roc not VTOL')
    plt.show()
    plt.plot(climb_distance_tab,energy_total_hop)
    plt.xlabel('climb distance m')
    plt.ylabel('Energy J')
    plt.title('Energy vs climb distance not VTOL')
    plt.show()

    return energy_climb,energy_cruise


def Energy_hop_cruise(v_stall,drag):
    energy_cruise=np.sqrt(100**2+100**2)*drag*(v_stall/20)**2
    print('energy_cruise [J]',energy_cruise)
    return energy_cruise



rho=1.225
Prop_calc()
Energy_hop_VTOL(12,15,15)
Energy_hop_NOT_VTOL(12,15,15)







