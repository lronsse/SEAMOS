"""
Water chamber sizing class
- Jan
"""
import numpy as np
from math import frexp
import matplotlib.pyplot as plt
# <editor-fold desc="CONSTANTS">
DENSITY_WATER = 998.
MASS_SYSTEM = 16.

RESERVOIR_PRESSURE = 300e5
RESERVOIR_VOLUME = 1.5/1e3
RESERVOIR_GAMMA = 1.4

CHAMBER_MAT_DENSITY = 2710.
CHAMBER_MAT_YIELD = 200E6
SAFETY_FACTOR = 2

SIGMA_ALLOW = CHAMBER_MAT_YIELD/SAFETY_FACTOR
C1 = CHAMBER_MAT_DENSITY / SIGMA_ALLOW * 1.5
# </editor-fold>


class WaterChamber:
    def __init__(self, pressure_reservoir=RESERVOIR_PRESSURE, volume_reservoir=RESERVOIR_VOLUME,
                 reservoir_gamma=RESERVOIR_GAMMA, density_water=DENSITY_WATER, mass_system=MASS_SYSTEM, c1=C1):
        self.pressure_reservoir = pressure_reservoir  # Pa
        self.volume_reservoir = volume_reservoir  # m3
        self.reservoir_gamma = reservoir_gamma  # -

        mantissa, exp = frexp(RESERVOIR_PRESSURE)
        print(mantissa, exp)
        self.pressure_plot_range = mantissa*np.logspace(-2, exp, 300, base=2)

        self.density_water = density_water  # kg/m3
        self.mass_system = mass_system  # kg
        self.c1 = c1

    def dV_from_PV(self, chamber_pressure, chamber_volume):
        # SI units
        if chamber_pressure < 1e5:
            print(f'WARNING: chamber pressure low {chamber_pressure/1e5} bar')
        if chamber_volume > 1.:
            print(f'WARNING: chamber volume high {chamber_pressure * 1e3} L')

        return np.sqrt(2*chamber_pressure/self.density_water) \
               * np.log(1 + self.density_water * chamber_volume/self.mass_system)

    def V_from_PdV(self, delta_velocity, chamber_pressure):
        # SI units
        if chamber_pressure < 1e5:
            print(f'WARNING: chamber pressure low {chamber_pressure/1e5} bar')
        if chamber_volume > 1:
            print(f'WARNING: chamber volume high {chamber_pressure * 1e3} L')

        return self.mass_system / self.density_water \
               * (np.exp(delta_velocity / np.sqrt(2 * chamber_pressure / self.density_water)) - 1)

    def plot_isentropic(self, discharges_per_reservoir=[1, 2, 3]):
        if type(discharges_per_reservoir) is not list:  # to list
            discharges_per_reservoir = [discharges_per_reservoir]
        capacity = self.pressure_reservoir * (self.volume_reservoir**self.reservoir_gamma)

        for n_discharge in discharges_per_reservoir:
            dv_list = []
            for chamber_pressure in self.pressure_plot_range:
                chamber_volume = (capacity/chamber_pressure/n_discharge)**(1/self.reservoir_gamma)
                dv_list.append(self.dV_from_PV(chamber_pressure, chamber_volume))
            plt.plot(self.pressure_plot_range/1e5, dv_list, label=f'isentropic, n = {n_discharge}')

    def plot_isoenthalpic(self, chamber_scaling=[0.5, 1, 2]):
        if type(chamber_scaling) is not list:
            chamber_scaling = [chamber_scaling]
        PV = self.pressure_reservoir*self.volume_reservoir

        for plot_counter, mass_factor in enumerate(chamber_scaling):
            PV_plot = PV*mass_factor
            mass = PV_plot * self.c1
            dv_list = []
            for chamber_pressure in self.pressure_plot_range:
                chamber_volume_plot = PV_plot/chamber_pressure
                dv_list.append(self.dV_from_PV(chamber_pressure, chamber_volume_plot))
            plt.plot(self.pressure_plot_range/1e5, dv_list,
                     label=f'const. v_c*p_c: m_c =  {mass*1E3} g')

    def plot_isochoric(self, chamber_volume, bounding_lines=[.8, 1, 1.2]):
        line_colors = ['grey', 'black', 'grey']

        for plot_counter, volume_factor in enumerate(bounding_lines):
            dv_list = []
            for chamber_pressure in self.pressure_plot_range:
                chamber_volume_plot = volume_factor*chamber_volume
                dv_list.append(self.dV_from_PV(chamber_pressure, chamber_volume_plot))
            plt.plot(self.pressure_plot_range/1e5, dv_list,
                     label=f'isochoric, V =  {chamber_volume_plot*1E3} L',
                     color=line_colors[plot_counter], linestyle='--')

    def show(self, title=''):
        plt.title(title)
        plt.xlabel('Chamber Pressure [bar]')
        plt.ylabel('Delivered Delta V [m/s]')
        plt.legend()
        plt.show()

if __name__ == '__main__':
    PB = WaterChamber(volume_reservoir=1.5/1000, pressure_reservoir=320E5)
    PB.plot_isentropic([1, 4, 8, 12, 16, 20])
    PB.plot_isochoric(5E-3)
    PB.show('Paintball Tank')

    # PBR.plot_isoenthalpic([.5, 1, 2])
    # PBR.plot_isochoric(5E-3)
    # PBR.show('Paintball Tank')

    CC = WaterChamber(volume_reservoir=.5/1000, pressure_reservoir=80E5)
    CC.plot_isentropic([1, 1, 1])
    CC.plot_isochoric(5E-3)
    CC.show('Combustion')

    # PB.plot_isoenthalpic([.5, 1, 2])
    # PB.plot_isochoric(5E-3)
    # PB.show('Combustion')







