import numpy as np
import matplotlib.pyplot as plt
import copy

# Chemical parameters
M_cac2 = 64.099  # [g/mol]
M_c2h2 = 26.04  # [g/mol]
M_co2 = 44.01  # [g/mol]
M_h20 = 18.01528  # [g/mol]
R = 8.314  # [J/(mol*K)]
T = 3000  # [K] # In air combustion
impulse = 16 * 18  # [Mass x dV]
rho_w = 1000  # [kg/m^3]

colors = ['tab:red', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:pink', 'tab:gray',
          'tab:olive', 'tab:cyan']

def n_to_mass(n, M):
    return n * M


def mass_to_n(mass, M):
    return mass / M


def mass():
    target_mass = 10  # g
    mass_step = 0.1  # [g]
    target_V = 2  # [L]
    V_step = 0.1  # [L]
    V = [V_step]
    for i in range(int(target_V / V_step)):
        mass_cac2 = [0]  # [g]
        P_gas = []
        for j in range(int(target_mass / mass_step)):
            n_cac2 = mass_to_n(mass_cac2[j], M_cac2)
            n_gas = 2 * n_cac2
            water_volume = n_cac2 * 18E-5
            n_oxygen = 5/2 * n_cac2
            V_oxygen = n_oxygen * 22.4E-3  # [m^3] # Required oxygen volume for combustion
            if V_oxygen > V[i] / 1000:
                print('oxygen volume required too high')
                break

            P_gas.append(n_gas * R * T / (V[i] / 1000))
            mass_cac2.append(round(mass_cac2[j] + mass_step, 2))
        plt.plot(mass_cac2[:-1], [pressure / 10000 for pressure in P_gas], label=f'Volume = {V[i]} [L]')
        V.append(round(V[i] + V_step, 5))

    plt.xlabel('Mass calcium carbide [g]')
    plt.ylabel('Pressure [bar]')


def volume():
    target_mass = 1.2  # g
    mass_step = 0.3  # [g]
    target_V = 10  # [L]
    V_step = 1  # [L]

    mass = [mass_step]
    for i in range(int(target_mass / mass_step)):
        V = [V_step]  # [L]
        P_gas = []
        for j in range(int(target_V / V_step)):
            n_cac2 = mass_to_n(mass[i], M_cac2)
            n_gas = 2 * n_cac2
            water_volume = n_cac2 * 18E-5
            P_gas.append(n_gas * R * T / (V[j] / 1000))
            V.append(round(V[j] + V_step, 5))

        n_oxygen = 5 / 2 * n_cac2
        V_oxygen = n_oxygen * 22.4E-3  # [m^3] # Required oxygen volume for combustion
        plt.axvline(V_oxygen * 10000, color=colors[i], linestyle='--')
        plt.plot(V[:-1], [pressure / 10000 for pressure in P_gas], label=f'mass = {mass[i]} [L]', color=colors[i])
        plt.xlabel('Chamber Volume [L]')
        plt.ylabel('Chamber Pressure [bar]')
        mass.append(round(mass[i] + mass_step, 5))


def volume_opti():
    target_mass = 1  # g
    mass_step = 1  # [g]
    target_V = 10  # [L]
    V_step = 0.1  # [L]

    mass = [mass_step]
    V_oxygen = []
    P_gas = []
    for i in range(int(target_mass / mass_step)):
        n_cac2 = mass_to_n(mass[i], M_cac2)
        n_gas = 3 * n_cac2
        n_oxygen = 5 / 2 * n_cac2
        oxy_volume = n_oxygen * 22.4E-3 * 1000 * 5  # [m^3] # Required oxygen volume for combustion
        P_gas.append(n_gas * R * T / (oxy_volume / 1000))
        V_oxygen.append(oxy_volume)
        mass.append(round(mass[i] + mass_step, 5))
        print('something is happening')
    print('out of loop', P_gas[-1])
    plt.plot(V_oxygen, [pressure / 10000 for pressure in P_gas], label=f'mass = {mass[i]} [L]', color=colors[i])
    print('plotted')


def solid():
    water_volume = 5  # [L]
    P_cc = 25E5  # [bar]
    a_target = 4  # [m/s^2]
    a_step = 0.5  # [m/s^2]
    m_uav = 16  # [kg]
    g = 9.81  # [m/s^2]
    A_cc = np.pi * (0.25/2) ** 2  # [m^2]
    V0_res = 0.1E-3  # [m^3]
    T_cc = 1400  # [K]
    r_rate = 4E-3  # [m/s]
    M_adn = 0.12406  # [kg/mol]
    rho_adn = 1.81E-3  # [kg/m^3]

    a = [0]
    for i in range(int(a_target / a_step)):
        a.append(a[i] + a_step)

        Ve_w = np.sqrt((P_cc - 1) * 2 / rho_w)
        Thrust = m_uav * a[i]
        m_dot = Thrust / Ve_w

        nozzle_diam = [0.05]
        A_solid, mass_solid, len_solid, t_water, Vol_gas = [], [], [], [], []
        j = 0
        while nozzle_diam[j] <= 0.15:
            t_water.append(water_volume / ((np.pi * (nozzle_diam[j] / 2) ** 2) * Ve_w))
            V_cc = (np.pi * (nozzle_diam[j] / 2) ** 2 * Ve_w) / A_cc
            Vol_gas.append(V0_res)  # + (A_cc * V_cc)
            n_gas = P_cc * Vol_gas[j] / (R * T_cc)
            mass_solid.append(n_gas * M_adn)  # in kg
            len_solid.append(r_rate * t_water[j])  # in m
            A_solid.append(mass_solid[j] / (rho_adn * len_solid[j]))
            nozzle_diam.append(round(nozzle_diam[j] + 0.02, 3))
            print('working')
            j += 1
        fig = plt.figure()
        plt.subplot(2, 2, 1)
        plt.plot(nozzle_diam[:-1], t_water)
        plt.xlabel('nozzle diameter [m]')
        plt.ylabel('water ejection time')
        plt.subplot(2, 2, 2)
        plt.plot(nozzle_diam[:-1], mass_solid)
        plt.xlabel('nozzle diameter [m]')
        plt.ylabel('mass of solid propellant')
        plt.subplot(2, 2, 3)
        plt.plot(nozzle_diam[:-1], len_solid)
        plt.xlabel('nozzle diameter [m]')
        plt.ylabel('length of solid propellant')
        plt.subplot(2, 2, 4)
        plt.plot(nozzle_diam[:-1], A_solid)
        plt.xlabel('nozzle diameter [m]')
        plt.ylabel('Area of solid propellant')
    plt.show()

    return a, diam_nozzle, V_gas, t_water, len_solid, mass_solid, A_solid

#mass()
#plt.legend()
#plt.show()

volume()
plt.legend()
plt.savefig('CaC2_sizing.png')

