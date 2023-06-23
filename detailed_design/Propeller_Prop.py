"""
1D energy solver for an object experiencing a variable drag force, optimised thrust,
 and a constant body force (e.g. gravity or buoyancy)
~ Jan
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimise


class VerticalPropulsion:
    def __init__(self, m: float, Cd: float, S: float, rho: float, force_vertical, v_intl,
                 dt=0.01, force_max_try=1000):
        """
        Object collecting parameters of 1 dimensional linear propulsion
        :param m: system mass for EOM [kg]
        :param Cd: drag coefficient [-]
        :param S: reference area [m2]
        :param rho: density [kg/m3]
        :param force_vertical: buoyancy or gravity force, positive inline with thrust [kgm/s2]
        :param v_intl: initial velocity [m/s]
        :param dt: sim timestep [s]
        :param force_max_try: max force in optimiser [N], high values fail to converge to optimum
        """
        self.specific_energy = 200*3600  # J/kg for energy
        self.specific_power = 20000  # W/kg for power

        self.m = m
        self.Cd = Cd
        self.S = S
        self.rho = rho
        self.force_buoyancy = force_vertical
        self.v0 = v_intl
        self.dt = dt
        self.max = force_max_try

        self.dV = None     # tracker for optimisation

        self.force_prop = None  # placeholder var for optimisation
        self.energy = None      # placeholder var for calc
        self.power = None       # placeholder var for calc
        self.mass_energy = None  # placeholder var for calc
        self.mass_power = None  # placeholder var for calc

    def drag(self, v: float) -> float:
        """
        Drag at velocity v
        :param v: aero-hydrodynamic velocity [m/s]
        return: drag [kgm/s2]
        """
        return .5*self.Cd*self.rho*self.S*v**2

    def energy_submerged_try(self, force_thrust: float) -> float:
        """
        1D submerged EOM resolution to energy requirement for buoyancy assisted underwater acceleration from v0 to v0+dV
        :param force_thrust: thrust force [N]
        return:
        """
        if force_thrust <= (self.drag(self.dV)-self.force_buoyancy):
            print(self.drag(self.dV))
            print(f'DEBUG: tried {force_thrust} [N] for {self.dV} [m/s], thrust too low, disqualifying')
            return 1e63  # equilibrium will never be reached

        v_0 = self.v0  # initial velocity
        v_i = v_0  # velocity tracker
        v_f = self.v0 + self.dV  # target velocity

        E_tot = 0.  # energy tracker

        while v_i < v_f:
            dv_dt = (1/self.m) * (force_thrust+self.force_buoyancy-self.drag(v_i))
            v_i += dv_dt*self.dt
            E_tot += force_thrust*v_i*self.dt  # energy per timestep, force * distance (velocity*time)

        return E_tot

    def energy_power_mass(self, dV):
        """
        Optimises energy use in 1D eom resolution for given dV
        Runs an optimiser to find the ideal thrust to accelerate to a final velocity with minimal energy
        :param dV: target velocity
        """
        self.dV = dV  # track

        opt = optimise.minimize_scalar(self.energy_submerged_try, bounds=(0., self.max), method='bounded')
        self.energy = opt.fun
        self.force_prop = opt.x
        self.power = self.dV*self.force_prop

        # cost
        self.mass_energy = self.energy/self.specific_energy
        self.mass_power = self.power/self.specific_power

        return self.energy, self.force_prop, self.power, self.mass_energy+self.mass_power


if __name__ == '__main__':
    """
    Main solution loop, initialising and plotting the 1D solver class for a submarine and airborne case
    """
    # Dict of initialisation parameters for the two cases
    Cd_shared = 0.05
    m_shared = 16.
    S_shared = 1.25
    rho_water = 998.
    rho_air = 1.
    env = {
        # 'Submerged': {
        #     'm': m_shared,
        #     'Cd': Cd_shared,
        #     'S': S_shared,
        #     'rho': rho_water,
        #     'force_vertical': 5.,  # Buoyancy
        #     'v_intl': 0.,
        #     'force_max_try': 100000.
        # },
        'Airborne': {
            'm': m_shared,
            'Cd': Cd_shared,
            'S': S_shared,
            'rho': rho_air,
            'force_vertical': -m_shared*9.81,  # Gravity
            'v_intl': 0.,
            'force_max_try': 1000.
        }
    }

    plot_params = {
        'Optimised Energy Consumption': [0, 'Energy [J]', 'dV [m/s]'],
        'Optimised Propulsive Force': [1, 'Force [N]', 'dV [m/s]'],
        'Optimised Propulsive Power': [2, 'Power [W]', 'dV [m/s]'],
        'Additional System Mass': [3, 'Mass [kg]', 'dV [m/s]']
    }
    for title, pars in plot_params.items():
        for name, params in env.items():
            study = VerticalPropulsion(**params)
            dV_range = np.arange(0.1, 30.1, .1)
            results_array = []
            for dV_plot in dV_range:
                idx = pars[0]
                results_array.append(study.energy_power_mass(dV_plot)[pars[0]])
            plt.plot(dV_range, results_array, label=name)

        plt.title(title)
        plt.xlabel(pars[2])
        plt.xlim(xmin=0.)
        plt.ylabel(pars[1])
        plt.legend()
        plt.show()

    study = VerticalPropulsion(**env['Airborne'])
    print(study.drag(1.0))
    print(study.drag(5.0))
    print(study.drag(10.0))