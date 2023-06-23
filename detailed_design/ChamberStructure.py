"""
Chamber structural sizing class

General functional flow is the determination of a minimal chamber thickness under various loading conditions

Wrapped to optimise chamber shape
"""
from WaterChamber import WaterChamber
import matplotlib.pyplot as plt
import numpy as np

# <editor-fold desc="CONSTANTS">
MATERIAL_E_MODULUS = 70E9  # Pa
MATERIAL_DENSITY = 2710.  # kg/m3
MATERIAL_RATIO_POISSON = 1/3  # -
MATERIAL_STRESS_ALLOW = 870E6 / 2  # Pa

RESERVOIR_VOLUME = 2.032 / 1000  # m3
RESERVOIR_PRESSURE = 310.2641E5  # Pa

CHAMBER_VOLUME = 6.5 / 1000  # m3
CHAMBER_PRESSURE = 25E5  # Pa

SUBMARINE_MAX_THRUST = 300  # N

MIN_PRESSURE_REG = 20  # Bar

MIN_DV_LAUNCH = 21  # m/s
MIN_THRUST_LAUNCH = 850  # N

CHAMBER_RADIUS = 357 / 2 / 1000  # m
# </editor-fold>


class Chamber:
    def __init__(self, pressure_chamber=CHAMBER_PRESSURE, chamber_volume=CHAMBER_VOLUME, chamber_radius=CHAMBER_RADIUS,
                 jet_thrust=MIN_THRUST_LAUNCH, submarine_axial_thrust=SUBMARINE_MAX_THRUST,
                 youngs_mod=MATERIAL_E_MODULUS, r_poisson=MATERIAL_RATIO_POISSON):
        self.pressure_chamber = pressure_chamber
        self.chamber_volume = chamber_volume

        self.jet_thrust = jet_thrust
        self.chamber_exit_radius = np.sqrt((self.jet_thrust / 2 / self.pressure_chamber) / np.pi)
        self.chamber_radius = chamber_radius

        self.submarine_axial_thrust = submarine_axial_thrust

        self.youngs_mod = youngs_mod
        self.r_poisson = r_poisson

    @property
    def chamber_taper_length(self):
        return self.chamber_volume * 3 / np.pi \
               / (self.chamber_radius ** 2 + self.chamber_radius * self.chamber_exit_radius
                  + self.chamber_exit_radius ** 2)

    @staticmethod
    def plot_chamber_sizing():
        _ch = WaterChamber(pressure_reservoir=RESERVOIR_PRESSURE, volume_reservoir=RESERVOIR_VOLUME)

        # _ch.plot_isentropic([1, 2])
        dv, pres = _ch.get_isentropic(1)
        plt.fill_between(pres, dv, 200., color='lightcoral')

        dv2, pres2 = _ch.get_isentropic(2)
        plt.fill_between(pres, dv, dv2, color='w', hatch='/', label='No Redundancy', fc='gold')

        plt.fill_between(pres, dv2, 21., color='w', hatch='x', label='Design w/ Redundancy', fc='green')

        # p_min
        plt.fill_betweenx([MIN_DV_LAUNCH, 100], MIN_PRESSURE_REG, color='lightcoral', label='Invalid')

        # dV_min
        plt.fill_between([0, 300], 21., color='lightcoral')

        # plot lines
        plt.axvline(MIN_PRESSURE_REG, label=f'p_c_min = {MIN_PRESSURE_REG}', color='red', linestyle=':', linewidth=2.)
        plt.axhline(MIN_DV_LAUNCH, label=f'dV_min = {MIN_DV_LAUNCH}', color='red', linestyle=':', linewidth=2.)
        plt.plot(pres, dv, label=f'Max. dV, single expulsion', linestyle='-.', color='red', linewidth=2.)
        plt.plot(pres2, dv2, label=f'Max. dV, two expulsions', linestyle='--', color='green')

        # design point
        _ch.plot_isochoric(CHAMBER_VOLUME, bounding_lines=[1.])
        plt.axvline(25., label=f'p_c = {25.}', color='black')

        _ch.show()

    def thickness_buckling_cone(self, gamma=1/3):
        """
        Min. thickness mitigating buckling of an axially loaded truncated cone NASA method
        :param: gamma [-], dimensionless constant. Unity for ideal, 1/3 for conservative estimate
        :return t: thickness [m]
        """
        a = np.arctan((self.chamber_radius-self.chamber_exit_radius)/self.chamber_taper_length)
        P_cr = self.submarine_axial_thrust

        t = np.sqrt((P_cr * np.sqrt(3*(1 - self.r_poisson**2))) / (gamma * 2 * np.pi * self.youngs_mod * np.cos(a)**2))

        return t

    def thickness_hoop(self):
        ...


if __name__ == '__main__':
    instance = Chamber()
    instance.plot_chamber_sizing()
    print(f'Chamber radius {instance.chamber_radius*1000} [mm]')
    print(f'Exit radius {instance.chamber_exit_radius*1000} [mm]')
    print(f'Taper length {instance.chamber_taper_length*1000} [mm] \n ...')
    # print(instance.thickness_buckling_cone())

    print(f'Compression thickness {MIN_THRUST_LAUNCH/(MATERIAL_STRESS_ALLOW*2*np.pi*instance.chamber_exit_radius)*1000} [mm]')
    print(f'Buckling thickness {instance.thickness_buckling_cone()*1000} [mm]')
    print(f'Hoop stress thickness {CHAMBER_PRESSURE*CHAMBER_RADIUS/(2*MATERIAL_STRESS_ALLOW)*1000} [mm]')