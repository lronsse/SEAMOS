import numpy as np
import math
import matplotlib.pyplot as plt
# from isa import
from ambiance import Atmosphere

"""
(written by Mikolaj)
Propeller sizing based on thrust requirement
Energy consumption for the hop(initial stall velocity, no initial velocity)
Assumed Thrust no longer 11 N yay now its more 4
"""


def Advance_ratio(v_airspeed, omega_rotor, D_prop):
    J = 2 * math.pi * v_airspeed / (omega_rotor * D_prop)
    return J


def Propeller_from_Thrust(thrust, v_airspeed, n_propellers, eta_rotor, v_prop_out):
    """
    Propeller diameter based on thrust needed, speed and efficiency, rho should be revisited
    """
    # omega_rotor=2*np.pi*rpm/60
    rho = Atmosphere(120).density
    thrust = thrust / n_propellers
    # S_swept_prop = np.pi * D_rotor ** 2 / 4
    s_swept_prop = thrust / (0.5 * rho * eta_rotor * ((v_prop_out) ** 2 - v_airspeed ** 2))
    d_rotor = np.sqrt(s_swept_prop * 4 / np.pi)
    # thrust = 0.5*rho*S_swept_prop*eta_rotor*((v_prop_out)**2-V_airspeed**2)
    return d_rotor


# print(Propeller_from_Thrust(15*9.81,12,1,0.9,24))


def Thrust_in_climb(climb_angle, drag, mass):
    thrust = mass * 9.81 * np.sin(climb_angle) + drag  # *0.08/0.05
    return thrust


def Pitch_calc(rpm, v):
    return 60 / rpm * v


def Prop_calc(roc, drag_cruise):
    """
    Just run this for the propeller sizing
    """
    # mass=15
    print('Puffin mass [kg]', mass)
    # drag_cruise = 24.5
    print('Drag in cruise', drag_cruise, '[N]')

    print('Rate of climb', roc, '[m/s]')
    v_airspeed = 20
    print('Cruise velocity', v_airspeed, '[m/s]')
    climb_angle = np.arcsin(float(roc / v_airspeed))  ##degree
    maxthrust = Thrust_in_climb(climb_angle, drag_cruise, mass)
    print('Thrust in climb [N]', maxthrust)
    print('Climb angle radians:', climb_angle / np.pi * 180)
    number_of_prop = 1
    print('number of propellers', number_of_prop)
    power = maxthrust * 20 / 0.8
    print(power, 'power')
    D_tab = []
    v_exit_tab = []
    torque_tab = []
    rpm_tab = np.linspace(100, 8100, 80)
    for rpm in rpm_tab:
        D = 0.1
        for i in range(100):
            v_exit = np.sqrt(v_airspeed ** 2 + rpm * 2 * np.pi / 60 * D / 2)
            D = Propeller_from_Thrust(Thrust_in_climb(climb_angle, drag_cruise, mass), v_airspeed, 1, 0.9, v_exit)

        v_exit_tab = np.append(v_exit_tab, v_exit)
        D_tab = np.append(D_tab, D)
        torque_tab = np.append(torque_tab, power / (rpm * 2 * np.pi / 60))
    plt.plot(rpm_tab, D_tab, label='Diameter vs RPM')
    # plt.plot(rpm_tab, v_exit_tab)
    plt.xlabel('RPM')
    plt.ylabel('Diameter [m]')
    plt.axvline(x=5000, color='r', label='The RPM chosen')
    plt.title('Diameter vs RPM')
    plt.legend(loc='upper right')
    plt.show()
    plt.plot(rpm_tab, torque_tab, label='Torque vs RPM')
    # plt.plot(rpm_tab, v_exit_tab)
    plt.xlabel('RPM')
    plt.ylabel('Torque [Nm]')
    plt.axvline(x=5000, color='r', label='The RPM chosen')
    plt.title('Torque vs RPM')
    plt.legend(loc='upper right')
    plt.show()
    print('Torque Nm', power / (5000 * 2 * np.pi / 60))

    rpm1 = float(input("Enter RPM:"))
    D = 0.1
    for i in range(100):
        v_exit = np.sqrt(v_airspeed ** 2 + rpm1 * 2 * np.pi / 60 * D / 2)
        D = Propeller_from_Thrust(Thrust_in_climb(climb_angle, drag_cruise, mass), v_airspeed, 1, 0.9, v_exit)
    print('Propeller diameter', D, '[m]')
    print('Propeller diameter', D * 39.37, '[inches]')
    #D = 0.47
    return D


def Thrust_in_climb_average(drag, mass):
    thrust = drag + mass * 9.81 * np.sin(np.linspace(np.pi / 2, 0, 1000))
    thrust = np.average(thrust)
    return thrust


def Energy_hop_VTOL(v_cruise, cruise_alt, mass):
    energy_climb = []
    energy_cruise = []

    roc_tab = np.linspace(0.01, 6.01, 6)
    roc_tab_2 = []
    climb_distance_tab = []
    for roc in roc_tab:
        climb_time = cruise_alt / roc
        climb_distance = (climb_time * np.sqrt(v_cruise ** 2 - roc ** 2))
        drag_climb_tab = []
        climb_angle = np.pi / 2
        velocity = 0
        for i in range(100):
            drag_climb_tab = np.append(drag_climb_tab, Thrust_in_climb(climb_angle, 10 * (velocity / 20) ** 2, 17))
            climb_angle = climb_angle - np.pi / 100
            velocity = velocity + 0.15
        drag_climb = np.average(drag_climb_tab)
        hop_distance = np.sqrt(100 ** 2 + 100 ** 2)
        if climb_distance < hop_distance:
            # print(climb_time,climb_distance)
            energy_climb = np.append(energy_climb, climb_time * (drag_climb + mass * v_cruise / climb_time) * v_cruise)

            energy_cruise = np.append(energy_cruise, (hop_distance - climb_distance) * 10)
            roc_tab_2 = np.append(roc_tab_2, roc)
            climb_distance_tab = np.append(climb_distance_tab, climb_distance)
        # print('D Vtol',Prop_calc_vtol(roc,drag_climb))
    # print('energy_climb[J]',energy_climb)
    # print('energy_cruise[J]', energy_cruise)
    # print(len(energy_climb),len(energy_cruise))
    plt.plot(roc_tab_2, energy_cruise, label='energy_cruise')
    plt.plot(roc_tab_2, energy_climb, label='energy_climb')
    energy_total_hop = energy_cruise + energy_climb
    print(print('energy_total_hop VTOL[J]', energy_total_hop[-1]))
    plt.plot(roc_tab_2, energy_total_hop, label='total energy per hop')
    plt.legend(loc='lower right')
    plt.xlabel('ROC m/s')
    plt.ylabel('Energy J')
    plt.title('Energy vs roc VTOL')
    plt.show()
    plt.plot(climb_distance_tab, energy_total_hop)
    plt.xlabel('climb distance m')
    plt.ylabel('Energy J')
    plt.title('Energy vs climb distance VTOL')
    plt.show()
    return energy_climb, energy_cruise


def Energy_hop_NOT_VTOL(v_cruise, cruise_alt, mass):
    energy_climb = []
    energy_cruise = []

    roc_tab = np.linspace(0.01, 10, 1000)
    roc_tab_2 = []
    climb_distance_tab = []
    for roc in roc_tab:
        climb_time = cruise_alt / roc
        climb_distance = (climb_time * np.sqrt(v_cruise ** 2 - roc ** 2))
        drag_climb_tab = []
        climb_angle = np.pi / 2
        velocity = 0
        for i in range(100):
            drag_climb_tab = np.append(drag_climb_tab, Thrust_in_climb(climb_angle, 10, 17))
            climb_angle = climb_angle - np.pi / 100
        drag_climb = np.average(drag_climb_tab)
        hop_distance = np.sqrt(100 ** 2 + 100 ** 2)
        if climb_distance < hop_distance:
            # print(climb_time,climb_distance)
            energy_climb = np.append(energy_climb, climb_time * (drag_climb * (15 / 20) ** 2) * v_cruise)

            energy_cruise = np.append(energy_cruise, (hop_distance - climb_distance) * 10)
            roc_tab_2 = np.append(roc_tab_2, roc)
            climb_distance_tab = np.append(climb_distance_tab, climb_distance)

    # print('energy_climb[J]',energy_climb)
    # print('energy_cruise[J]', energy_cruise)
    # print(len(energy_climb),len(energy_cruise))
    plt.plot(roc_tab_2, energy_cruise, label='energy_cruise')
    plt.plot(roc_tab_2, energy_climb, label='energy_climb')
    energy_total_hop = energy_cruise + energy_climb
    print(print('energy_total_hop not vtol[J]', energy_total_hop[-1]))
    plt.plot(roc_tab_2, energy_total_hop, label='total energy per hop')
    plt.legend(loc='lower right')
    plt.xlabel('ROC m/s')
    plt.ylabel('Energy J')
    plt.title('Energy vs roc not VTOL')
    plt.show()
    plt.plot(climb_distance_tab, energy_total_hop)
    plt.xlabel('climb distance m')
    plt.ylabel('Energy J')
    plt.title('Energy vs climb distance not VTOL')
    plt.show()

    return energy_climb, energy_cruise


def Energy_hop_cruise(v_stall, drag):
    energy_cruise = np.sqrt(100 ** 2 + 100 ** 2) * drag * (v_stall / 20) ** 2
    print('energy_cruise [J]', energy_cruise)
    return energy_cruise


def Glide_range(cruise_altitude, lift_to_drag):
    glide_range = cruise_altitude * lift_to_drag
    return glide_range


def Rho_average(altitude):
    altitudes = np.linspace(0, altitude, 1000)
    atmosphere = Atmosphere(altitudes)
    return np.average(atmosphere.density)


def Optimal_flight_energy(roc):
    # roc=1
    range = 30000
    v_airspeed = 20
    eta_motor = .8

    lift_to_drag = 25
    c_l = 1.365
    c_d = 0.05  # cruise
    c_d_0 = 0.05
    c_d_max = 0.1
    # cl/cd 25
    mass = 17
    wing_area = 0.75
    eta_prop = 0.8
    drag_in_cruise = 0.5 * wing_area * v_airspeed ** 2 * Atmosphere(120).density * c_d
    time_tab = [[], [], [], []]
    distance_tab = [[], [], []]
    energy_tab = [[], [], []]
    altitude_tab = np.linspace(0, 1000, 1000)
    for altitude in altitude_tab:
        climb_time = altitude / roc
        climb_distance = (climb_time * np.sqrt(v_airspeed ** 2 - roc ** 2))
        climb_angle = np.arcsin(float(roc / v_airspeed))
        energy_climb = Thrust_in_climb(climb_angle, drag_in_cruise, mass) * climb_distance / eta_motor / eta_prop

        glide_range = Glide_range(altitude, lift_to_drag)
        glide_time = altitude / (c_d / c_l * np.sqrt(mass * 9.81 / (0.5 * Rho_average(altitude) * wing_area * c_l)))

        if glide_range + climb_distance >= range:
            cruise_range = 0
            cruise_time = 0
            energy_cruise = 0
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
            altitude_tab = altitude_tab[0:len(time_tab[0])]
            break
        else:
            cruise_range = range_meters - climb_distance - glide_range
            atmosphere = Atmosphere(altitude)

            cruise_drag = v_airspeed ** 2 * float(atmosphere.density[0]) * 0.5 * c_d * wing_area
            energy_cruise = cruise_drag / eta_motor * cruise_range / eta_prop
            cruise_time = cruise_range / v_airspeed

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

    return time_tab, distance_tab, energy_tab, altitude_tab


def Plot_optimal(time_tab, distance_tab, energy_tab, altitude_tab, roc):
    altitude1 = int(input("Enter altitude:"))
    for altitude in altitude_tab:
        if altitude > altitude1:
            index1 = list(altitude_tab).index(altitude)
            break
    print('climb time [min]', time_tab[0][index1] / 60)
    print('cruise time [min]', time_tab[1][index1] / 60)
    print('glide time [min]', time_tab[2][index1] / 60)
    print('climb distance [m]', distance_tab[0][index1])
    print('cruise distance [m]', distance_tab[1][index1])
    print('glide distance [m]', distance_tab[2][index1])
    print('climb energy [J]', energy_tab[0][index1])
    print('cruise energy [J]', energy_tab[1][index1])
    print('Total energy [J]', energy_tab[2][index1])
    plt.plot(altitude_tab, energy_tab[2], label='total energy')
    plt.plot(altitude_tab, energy_tab[0], label='energy climb')
    plt.plot(altitude_tab, energy_tab[1], label='energy cruise')
    plt.xlabel('altitude [m]')
    plt.ylabel('energy[J]')
    plt.legend(loc='upper left')
    plt.title('Energy vs cruise altitude' + ' for ROC ' + str(roc) + 'm/s')
    plt.show()
    plt.plot(altitude_tab, time_tab[3], label='total time')
    plt.plot(altitude_tab, time_tab[0], label='time climb')
    plt.plot(altitude_tab, time_tab[1], label='time cruise')
    plt.plot(altitude_tab, time_tab[2], label='glide time')
    plt.xlabel('altitude [m]')
    plt.ylabel('time[s]')
    plt.legend(loc='upper left')
    plt.title('Time vs cruise altitude' + ' for ROC ' + str(roc) + 'm/s')
    plt.show()

    plt.plot(altitude_tab, distance_tab[0], label='distance climb')
    plt.plot(altitude_tab, distance_tab[1], label='range cruise')
    plt.plot(altitude_tab, distance_tab[2], label='glide range')
    plt.xlabel('altitude [m]')
    plt.ylabel('distance[m]')
    plt.legend(loc='upper left')
    plt.title('distance vs cruise altitude' + ' for ROC ' + str(roc) + 'm/s')
    plt.show()


def Plot_for_different_roc(drag_in_cruise):
    time_tab1, distance_tab1, energy_tab1, altitude_tab1 = Optimal_flight_energy2(0.5, drag_in_cruise)
    time_tab2, distance_tab2, energy_tab2, altitude_tab2 = Optimal_flight_energy2(1, drag_in_cruise)
    time_tab3, distance_tab3, energy_tab3, altitude_tab3 = Optimal_flight_energy2(2, drag_in_cruise)
    time_tab4, distance_tab4, energy_tab4, altitude_tab4 = Optimal_flight_energy2(4, drag_in_cruise)
    time_tab5, distance_tab5, energy_tab5, altitude_tab5 = Optimal_flight_energy2(6, drag_in_cruise)
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


def Prop_calc_lite(roc, rpm1):
    mass = 17
    drag_cruise = drag_in_cruise
    v_airspeed = 20
    climb_angle = np.arcsin(float(roc / v_airspeed))  ##degree
    thrust_max = Thrust_in_climb(climb_angle, drag_cruise, mass)
    number_of_prop = 1
    D_tab = []
    v_exit_tab = []

    D = 0.1
    for i in range(100):
        v_exit = np.sqrt(v_airspeed ** 2 + rpm1 * 2 * np.pi / 60 * D / 2)
        D = Propeller_from_Thrust(Thrust_in_climb(climb_angle, drag_cruise, mass), v_airspeed, 1, 0.9, v_exit)
        v_exit_tab = np.append(v_exit_tab, v_exit)
        D_tab = np.append(D_tab, D)
    return D, thrust_max


def Plot_diameter_vs_roc(rpm):
    n = 4
    v_airspeed = 20
    roc = np.linspace(1, 6, n)
    print(roc)
    D_tab = np.zeros(n)
    thrusts_tab = np.zeros(n)
    for i in range(len(roc)):
        D_tab[i] = Prop_calc_lite(roc[i], rpm)[0]
        climb_angle = np.arcsin(float(roc[i] / v_airspeed))
        thrusts_tab[i] = Thrust_in_climb(climb_angle, drag_in_cruise, 17)
    plt.plot(roc, D_tab)
    plt.ylabel('diameter [m]')
    plt.xlabel('Rate of Climb [m/s]')

    plt.title('Propeller diameter' + ' for different ROC ')
    plt.show()
    plt.plot(roc, thrusts_tab)
    plt.ylabel('Thrust [N]')
    plt.xlabel('Rate of Climb [m/s]')

    plt.title('Thrust required' + ' for different ROC ')
    plt.show()


def Prop_calc_vtol(roc, drag_cruise):
    mass = 17
    print('Puffin mass', mass)
    # drag_cruise = 4
    print('Drag in cruise', drag_cruise, '[N]')

    print('Rate of climb', roc, '[m/s]')
    v_airspeed = 20
    print('Cruise velocity', v_airspeed, '[m/s]')
    climb_angle = np.arcsin(float(roc / v_airspeed))  ##degree
    maxthrust = Thrust_in_climb(climb_angle, drag_cruise, mass)
    print('Thrust in climb [N]', maxthrust)
    print('Climb angle radians:', climb_angle)
    number_of_prop = 1
    print('number of propellers', number_of_prop)
    power = maxthrust * 20 / 0.8
    print(power, 'power')
    D_tab = []
    v_exit_tab = []
    torque_tab = []
    rpm_tab = np.linspace(100, 8100, 100)
    for rpm in rpm_tab:
        D = 0.1
        for i in range(100):
            v_exit = np.sqrt(v_airspeed ** 2 + rpm * 2 * np.pi / 60 * D / 2)
            D = Propeller_from_Thrust(Thrust_in_climb(climb_angle, drag_cruise, mass), v_airspeed, 1, 0.9, v_exit)

        v_exit_tab = np.append(v_exit_tab, v_exit)
        D_tab = np.append(D_tab, D)
        torque_tab = np.append(torque_tab, power / (rpm * 2 * np.pi / 60))

    print('Torque Nm', power / (5000 * 2 * np.pi / 60))
    rpm1 = float(input("Enter RPM:"))
    D = 0.1
    for i in range(100):
        v_exit = np.sqrt(v_airspeed ** 2 + rpm1 * 2 * np.pi / 60 * D / 2)
        D = Propeller_from_Thrust(Thrust_in_climb(climb_angle, drag_cruise, mass), v_airspeed, 1, 0.9, v_exit)
    print('Propeller diameter', D, '[m]')
    print('Propeller diameter', D * 39.37, '[inches]')

    return D


def Optimal_flight_energy2(roc, drag_in_cruise):
    # roc=1
    range_meters = 30000
    time_tab = [[], [], [], []]
    distance_tab = [[], [], []]
    energy_tab = [[], [], []]
    altitude_tab = np.linspace(0, 1000, 1000)
    for altitude in altitude_tab:
        climb_time = altitude / roc
        climb_distance = (climb_time * np.sqrt(v_airspeed ** 2 - roc ** 2))
        climb_angle = np.arcsin(float(roc / v_airspeed))
        energy_climb = 278 * climb_time

        glide_range = Glide_range(altitude, lift_to_drag)
        glide_time = altitude / (c_d / c_l * np.sqrt(mass * 9.81 / (0.5 * Rho_average(altitude) * wing_area * c_l)))

        if glide_range + climb_distance >= range_meters:
            cruise_range = 0
            cruise_time = 0
            energy_cruise = 0
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
            altitude_tab = altitude_tab[0:len(time_tab[0])]
            break
        else:
            cruise_range = range_meters - climb_distance - glide_range
            atmosphere = Atmosphere(altitude)

            cruise_drag = v_airspeed ** 2 * float(atmosphere.density[0]) * 0.5 * c_d * wing_area

            cruise_time = cruise_range / v_airspeed
            energy_cruise = 123 * cruise_time

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

    return time_tab, distance_tab, energy_tab, altitude_tab


def Sizing_for_VTOL(mass, climb_altitude, stall_speed):
    thrust_tab = []
    eta_rotor = 0.9
    rpm1 = 5000
    for acceleration in range(5, 44, 1):
        thrust_max = 0
        speed = 0
        time = 0
        energy = 0
        while speed < stall_speed:
            thrust_append = mass * 9.81 + 0.5 * 0.75 * speed ** 2 * 0.05 * Atmosphere(0).density + mass * acceleration
            speed = speed + acceleration * 0.01
            time = time + 0.01
            energy = energy + 0.01 * speed * thrust_append
            if thrust_append > thrust_max: thrust_max = thrust_append
        print('acceleration', acceleration, 'thrust max[N]', thrust_max, 'energy', energy)
        D = 0.1
        for i in range(100):
            v_exit = np.sqrt(v_airspeed ** 2 + rpm1 * 2 * np.pi / 60 * D / 2)
            s_swept_prop = thrust_max / (0.5 * Atmosphere(0).density * eta_rotor * ((v_exit) ** 2 - v_airspeed ** 2))
            D = np.sqrt(s_swept_prop * 4 / np.pi)

        height = acceleration * time ** 2 / 2
        print(height)
        print('diameter', D)
    acceleration = float(input('Pick the acceleration buckaroo'))
    thrust_max = 0
    speed = 0
    time = 0
    energy = 0
    while speed < stall_speed:
        thrust_append = mass * 9.81 + 0.5 * 0.75* speed ** 2 * 0.05 * Atmosphere(0).density + mass * acceleration
        speed = speed + acceleration * 0.01
        time = time + 0.01
        energy = energy + (0.01 * speed * thrust_append) / 0.8
        if thrust_append > thrust_max: thrust_max = thrust_append
    height = acceleration * time ** 2 / 2
    print('acceleration', acceleration, 'thrust max[N]', thrust_max, 'energy', energy, 'height', height)
    D = 0.1
    for i in range(100):
        v_exit = np.sqrt(v_airspeed ** 2 + rpm1 * 2 * np.pi / 60 * D / 2)
        s_swept_prop = thrust_max / (0.5 * Atmosphere(0).density * eta_rotor * ((v_exit) ** 2 - v_airspeed ** 2))
        D = np.sqrt(s_swept_prop * 4 / np.pi)
    roc = 1
    climb_angle = np.arcsin(float(roc / stall_speed))
    thrust_climb = mass * 9.81 * np.sin(climb_angle) + 0.5 * 0.75 * speed ** 2 * 0.05 * Atmosphere(
        0).density  # *0.08/0.05
    climb_time = (climb_altitude - height) / roc
    energy_climb = thrust_climb * stall_speed * climb_time
    climb_distance = (climb_time * np.sqrt(stall_speed ** 2 - roc ** 2))
    energy = energy + energy_climb
    print(energy)
    glide_range = climb_altitude * lift_to_drag
    hop_distance = np.sqrt(100 ** 2 + 100 ** 2)
    cruise_distance = hop_distance - glide_range - climb_distance
    if cruise_distance > 0:
        energy_cruise = 0.5 * 0.75 * speed ** 2 * 0.05 * Atmosphere(0).density * cruise_distance
    else:
        energy_cruise = 0
    energy = energy + energy_cruise
    print(energy)


# print(0.5*Rho_average(100)*20**2*1.25*0.05)
# Energy_hop_VTOL(12,15,15)
# Energy_hop_NOT_VTOL(12,15,15)
# print(20*np.sin(20/180*np.pi))
#stall speed 15
roc = 1
range_meters = 30000
v_airspeed = 20
eta_motor = .8
lift_to_drag = 20
c_l = 1.356
c_d = 0.05  # cruise
mass = 17
wing_area = 0.75
eta_prop = 0.8
drag_in_cruise = 0.5 * wing_area * v_airspeed ** 2 * Atmosphere(120).density * c_d
Sizing_for_VTOL(17,15,15)
print(drag_in_cruise,'drag in cruise [N]')
Prop_calc(roc,drag_in_cruise)
time_tab,distance_tab,energy_tab,altitude_tab=Optimal_flight_energy(roc)
Plot_optimal(time_tab,distance_tab,energy_tab,altitude_tab,roc)
time_tab,distance_tab,energy_tab,altitude_tab=Optimal_flight_energy2(roc,drag_in_cruise)
Plot_optimal(time_tab,distance_tab,energy_tab,altitude_tab,roc)
Plot_for_different_roc(drag_in_cruise)

Plot_diameter_vs_roc(float(input('RPM again pls')))
