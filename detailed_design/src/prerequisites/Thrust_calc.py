import numpy as np
import math
import matplotlib.pyplot as plt
#from isa import
from ambiance import Atmosphere

"""
(written by Mikolaj)
Propeller sizing based on thrust requirement
Energy consumption for the hop(initial stall velocity, no initial velocity)
Assumed Thrust no longer 11 N yay
"""

def Advance_ratio(v_airspeed,omega_rotor,D_prop):
    J = 2*math.pi*v_airspeed/(omega_rotor*D_prop)
    return J


def Propeller_from_Thrust(thrust,v_airspeed,n_propellers,eta_rotor,v_prop_out):
    """
    Propeller diameter based on thrust needed, speed and efficiency, rho should be revisited
    """
    #omega_rotor=2*np.pi*rpm/60
    rho=Atmosphere(400).density
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


def Prop_calc(roc):
    """
    Just run this for the propeller sizing
    """
    mass=15
    print('Puffin mass',mass)
    drag_cruise = 4
    print('Drag in cruise',drag_cruise,'[N]')

    print('Rate of climb', roc, '[m/s]')
    v_airspeed=20
    print('Cruise velocity', v_airspeed, '[m/s]')
    climb_angle = np.arcsin(float(roc / v_airspeed))  ##degree
    print('Thrust in climb [N]', Thrust_in_climb(climb_angle, drag_cruise, 15))
    print('Climb angle radians:', climb_angle)
    number_of_prop=1
    print('number of propellers',number_of_prop)
    D_tab=[]
    v_exit_tab=[]
    rpm_tab=np.linspace(100,8100,100)
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
    rpm1 = float(input("Enter RPM:"))
    D = 0.1
    for i in range(100):
        v_exit = np.sqrt(v_airspeed ** 2 + rpm1 * 2 * np.pi / 60 * D / 2)
        D = Propeller_from_Thrust(Thrust_in_climb(climb_angle, drag_cruise, mass), v_airspeed, 1, 0.9, v_exit)
    print('Propeller diameter',D,'[m]')
    return D


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
    #print(len(energy_climb),len(energy_cruise))
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
    #print(len(energy_climb),len(energy_cruise))
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


def Glide_range(cruise_altitude, lift_to_drag):
    glide_range=cruise_altitude*lift_to_drag
    return glide_range


def Rho_average(altitude):
    altitudes=np.linspace(0,altitude,1000)
    atmosphere=Atmosphere(altitudes)
    return np.average(atmosphere.density)


def Optimal_flight_energy(roc):
    #roc=1
    range=30000
    v_airspeed=20
    eta_motor=.8
    drag_in_cruise=4
    lift_to_drag=25
    c_l=1.4
    c_d= 0.013#cruise
    c_d_0=0.05
    c_d_max=0.1
    #cl/cd 25
    mass=15
    wing_area=1.25
    eta_prop=0.8
    time_tab=[[],[],[],[]]
    distance_tab=[[],[],[]]
    energy_tab=[[],[],[]]
    altitude_tab=np.linspace(0,1500,1000)
    for altitude in altitude_tab:
        climb_time= altitude/roc
        climb_distance = (climb_time*np.sqrt(v_airspeed**2-roc**2))
        climb_angle = np.arcsin(float(roc / v_airspeed))
        energy_climb = Thrust_in_climb(climb_angle,drag_in_cruise,mass)*climb_distance/eta_motor/eta_prop

        glide_range=Glide_range(altitude,lift_to_drag)
        glide_time=altitude/(c_d/c_l*np.sqrt(mass*9.81/(0.5*Rho_average(altitude)*wing_area*c_l)))

        if glide_range+climb_distance>=range:
            cruise_range=0
            cruise_time=0
            energy_cruise=0
            time_tab[0] = np.append(time_tab[0], climb_time)
            time_tab[1] = np.append(time_tab[1], cruise_time)
            time_tab[2] = np.append(time_tab[2], glide_time)
            time_tab[3] = np.append(time_tab[3], climb_time + cruise_time + glide_time)
            distance_tab[0] = np.append(distance_tab[0], climb_distance)
            distance_tab[1] = np.append(distance_tab[1], cruise_range)
            distance_tab[2] = np.append(distance_tab[2], glide_range)
            energy_tab[0] = np.append(energy_tab[0], energy_climb)
            energy_tab[1] = np.append(energy_tab[1], energy_cruise)
            energy_tab[2] = np.append(energy_tab[2], energy_climb + energy_cruise)
            altitude_tab=altitude_tab[0:len(time_tab[0])]
            break
        else:
            cruise_range = range - climb_distance - glide_range
            atmosphere=Atmosphere(altitude)

            cruise_drag=v_airspeed**2*float(atmosphere.density[0])*0.5*c_d*wing_area
            energy_cruise = cruise_drag/eta_motor*cruise_range/eta_prop
            cruise_time=cruise_range/v_airspeed

        time_tab[0] = np.append(time_tab[0],climb_time)
        time_tab[1] = np.append(time_tab[1],cruise_time )
        time_tab[2] = np.append(time_tab[2], glide_time)
        time_tab[3] = np.append(time_tab[3], climb_time+cruise_time+glide_time)
        distance_tab[0] = np.append(distance_tab[0], climb_distance)
        distance_tab[1] = np.append(distance_tab[1], cruise_range )
        distance_tab[2] = np.append(distance_tab[2], glide_range )
        energy_tab[0] = np.append(energy_tab[0], energy_climb)
        energy_tab[1] = np.append(energy_tab[1],energy_cruise )
        energy_tab[2] = np.append(energy_tab[2], energy_climb+energy_cruise)
    return time_tab,distance_tab,energy_tab,altitude_tab


def Plot_optimal(time_tab,distance_tab,energy_tab,altitude_tab,roc):

    plt.plot(altitude_tab, energy_tab[2],label='total energy')
    plt.plot(altitude_tab, energy_tab[0],label='energy climb')
    plt.plot(altitude_tab, energy_tab[1],label='energy cruise')
    plt.xlabel('altitude [m]')
    plt.ylabel('energy[J]')
    plt.legend(loc='upper left')
    plt.title('Energy vs cruise altitude'+' for ROC '+str(roc)+'m/s')
    plt.show()
    plt.plot(altitude_tab, time_tab[3], label='total time')
    plt.plot(altitude_tab, time_tab[0], label='time climb')
    plt.plot(altitude_tab, time_tab[1], label='time cruise')
    plt.plot(altitude_tab, time_tab[2], label='glide time')
    plt.xlabel('altitude [m]')
    plt.ylabel('time[s]')
    plt.legend(loc='upper left')
    plt.title('Time vs cruise altitude'+' for ROC '+str(roc)+'m/s')
    plt.show()

    plt.plot(altitude_tab, distance_tab[0], label='distance climb')
    plt.plot(altitude_tab, distance_tab[1], label='range cruise')
    plt.plot(altitude_tab, distance_tab[2], label='glide range')
    plt.xlabel('altitude [m]')
    plt.ylabel('distance[m]')
    plt.legend(loc='upper left')
    plt.title('distance vs cruise altitude'+' for ROC '+str(roc)+'m/s')
    plt.show()


def Plot_for_different_roc():
    time_tab1, distance_tab1, energy_tab1, altitude_tab1 = Optimal_flight_energy(0.5)
    time_tab2, distance_tab2, energy_tab2, altitude_tab2 = Optimal_flight_energy(1)
    time_tab3, distance_tab3, energy_tab3, altitude_tab3 = Optimal_flight_energy(2)
    time_tab4, distance_tab4, energy_tab4, altitude_tab4 = Optimal_flight_energy(4)
    time_tab5, distance_tab5, energy_tab5, altitude_tab5 = Optimal_flight_energy(6)
    plt.plot(altitude_tab1[0:400], energy_tab1[2][0:400], label='total energy 0.5 [m/s]')
    plt.plot(altitude_tab2[0:400], energy_tab2[2][0:400], label='total energy 1 [m/s]')
    plt.plot(altitude_tab3[0:400], energy_tab3[2][0:400], label='total energy 2 [m/s]')
    plt.plot(altitude_tab4[0:400], energy_tab4[2][0:400], label='total energy 4 [m/s]')
    plt.plot(altitude_tab5[0:400], energy_tab5[2][0:400], label='total energy 6 [m/s]')

    plt.xlabel('altitude [m]')
    plt.ylabel('energy[J]')
    plt.legend(loc='upper left')
    plt.title('Energy vs cruise altitude' + ' for different ROC ')
    plt.show()

def Prop_calc_lite(roc,rpm1):
    mass = 15
    drag_cruise = 4
    v_airspeed = 20
    climb_angle = np.arcsin(float(roc / v_airspeed))  ##degree
    thrust_max= Thrust_in_climb(climb_angle, drag_cruise, mass)
    number_of_prop = 1
    D_tab = []
    v_exit_tab = []

    D = 0.1
    for i in range(100):
        v_exit = np.sqrt(v_airspeed ** 2 + rpm1 * 2 * np.pi / 60 * D / 2)
        D = Propeller_from_Thrust(Thrust_in_climb(climb_angle, drag_cruise, mass), v_airspeed, 1, 0.9, v_exit)
        v_exit_tab = np.append(v_exit_tab, v_exit)
        D_tab = np.append(D_tab, D)
    return D,thrust_max


def Plot_diameter_vs_roc():
    n=4
    rpm=5000
    roc=np.linspace(1,6,n)
    print(roc)
    D_tab=np.zeros(n)
    for i in range(len(roc)):
        D_tab[i]=Prop_calc_lite(roc[i],rpm)[0]
    plt.plot(roc,D_tab)
    plt.ylabel('diameter [m]')
    plt.xlabel('Rate of Climb [m/s]')

    plt.title('Propeller diameter' + ' for different ROC ')
    plt.show()

print(0.5*Rho_average(200)*20**2*1.25*0.013)
#rho=1.225

#Energy_hop_VTOL(12,15,15)
#Energy_hop_NOT_VTOL(12,15,15)
print(20*np.sin(20/180*np.pi))

roc=1
Prop_calc(roc)
time_tab,distance_tab,energy_tab,altitude_tab=Optimal_flight_energy(roc)
Plot_optimal(time_tab,distance_tab,energy_tab,altitude_tab,roc)
Plot_for_different_roc()
Plot_diameter_vs_roc()





