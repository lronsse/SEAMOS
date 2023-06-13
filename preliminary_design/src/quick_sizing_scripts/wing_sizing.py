"""
Inputs: v_stall, v_cruise, configuration, mtom, n_hops, energy density battery, clmax, rate of climb, h_hop, n_hops, roc_vtol

required inputs: mtom, nhops, thop, rocvtol

assumptions: rocvtol = 1m/s, roccruise = 1

Outputs: Wing surface, battery mass, power/thrust required, lift over drag

todo = nhops plot, provide l/d to Jan, tell billie size watts or kg, add timeline with elli, miko: weight + nhops

"""
import matplotlib.pyplot as plt
import numpy as np
#import ElectricDroneSizing as eds


class Configuration:
    def __init__(self, name, v_stall, v_cruise, cl_max, cd_0, cd_vtol, eta_p, e, ar, h_cruise, h_ceiling, mtom, roc,
                 VTOL, n_prop, h_hop, n_hops, roc_vtol, power_ttc):
        """
        class of aircraft type

        :param name: name of the architecture
        :param v_stall: stall speed of the architecture
        :param v_cruise: cruise speed of the architecture
        :param cl_max: maximum lift coefficient of the architecture
        :param cd_0: zero-lift-drag coefficient
        :param cd_vtol: drag coefficient in case of vtol
        :param eta_p: propulsive efficiency
        :param e: oswald efficiency factor
        :param ar: aspect ratio
        :param h_cruise: cruise altitude
        :param h_ceiling: service ceiling altitude
        :param mtom: maxium take-off mass
        :param roc: rate off climb required during flight
        :param VTOL: True/False, whether you need vtol or not
        :param n_prop: amount of proppelors used (use one for general comparison)
        :param t_hop: TO + Land + Hover time per hop
        :param n_hops: amount of hops needed per mission
        :param roc_vtol: rate of climb required for vtol
        """
        self.name = name
        self.v_stall = v_stall
        self.v_cruise = v_cruise
        self.v_max = self.v_cruise * 1.25
        self.cl_max = cl_max
        self.cd_0 = cd_0
        self.cd_vtol = cd_vtol
        self.eta_p = eta_p
        self.e = e
        self.ar = ar
        self.h_cruise = h_cruise
        self.h_ceiling = h_ceiling
        self.k = 1 / (np.pi * e * ar)
        self.lift_drag_max = 10
        self.mtom = mtom
        self.v_take_off = self.v_stall * 1.15
        self.cl_take_off = self.cl_max / 1.21
        self.cd_take_off = self.cd_0 + self.k * self.cl_take_off ** 2
        self.roc = roc
        self.n_prop = n_prop
        self.s_ratio = 1.2
        self.t_cruise = d_flight / self.v_cruise
        self.h_hop = h_hop
        self.n_hops = n_hops
        self.roc_vtol = roc_vtol
        self.plot = True
        optimisation = True
        if optimisation:
            self.plot = False
        self.power_consumption_flight_controller = 20
        self.power_consumption_ttc = power_ttc
        self.power_consumption_payload = 35

        # self.t_mission = 4 * 3600
        plot_nhops = False

        def calc_isa(h):
            """
            function calculating the isa properties at a function

            :type h: float
            :arg h: altitude at which you want to fly
            :return T: temperature at altitude, p: pressure at altitude, rho: density at altitude
            """
            a = -0.0065
            p0 = 101325.0
            T0 = 288.15
            h0 = 0
            R = 287
            T = T0 + a * (h - h0)
            p = p0 * (T / T0) ** (-(g / (a * R)))
            rho = p / (R * T)
            return T, p, rho

        def wing_loading_stall(v_stall, cl_max):
            """
            function determining the wing loading required to prevent stalling before required stall speed

            :param v_stall: Minimum stall speed determined
            :param cl_max: Maximum lift coefficient
            :return: wing_loading: wing loading at stall
            """
            wing_loading = 0.5 * v_stall ** 2 * rho_0 * cl_max
            return wing_loading

        def required_power_flight(v, ws):
            """
            Function describing the required power to fly at max speed

            :param v_max: maximum speed of the plane
            :param ws: wing loading
            :return: weight over power required to fly at max speed
            """
            weight_power = eta_p / (0.5 * rho_0 * v ** 3 * cd_0 / ws + 2 * self.k * ws / (self.rho_cruise * v))
            return weight_power

        def required_power_ceiling(ws):
            """
            Function determining the weight over power required to fly at service ceiling

            :param ws: wing loading
            :return: weight over power required to fly at max speed
            """
            rho_ceiling = calc_isa(h_ceiling)[2]
            sigma_ceiling = rho_0 / rho_ceiling
            weight_power = sigma_ceiling / (np.sqrt(2 * ws / (rho_ceiling * np.sqrt(3 * cd_0 / self.k))) * 1.155 / (
                        self.lift_drag_max * eta_p))
            return weight_power

        def required_power_rate_of_climb(ws):
            weight_power = 1 / (
                        self.roc / eta_p + np.sqrt(2 * ws / (self.rho_cruise * np.sqrt(3 * cd_0 / self.k))) * 1.155 / (
                            self.lift_drag_max * eta_p))
            return weight_power

        if VTOL:
            def required_power_vtol(ws):
                """
                function determining the power required to take off vertically
                :param ws: wing sizing
                :return: power required to take off vertically
                """
                tw = 1.2 * (1 + cd_vtol * 0.5 * rho_0 * self.roc_vtol ** 2 * self.s_ratio / ws)
                treq = tw * self.mtom * g
                fm = 0.4742 * treq ** 0.0793
                dl = 3.2261 * self.mtom + 74.991
                sp = self.mtom * g / (dl * self.n_prop)
                vh = np.sqrt(treq / (2 * rho_0 * sp))
                vi = (-self.roc_vtol / (2 * vh) + np.sqrt((self.roc_vtol / (2 * vh)) ** 2 + 1)) * vh
                pwreq = 1 / (treq * vi / (fm * self.mtom * g))
                return pwreq

        def lift_drag(ws):
            tw = 0.5 * self.rho_cruise * v_cruise ** 2 * cd_0 / ws + self.k * ws / (
                        0.5 * self.rho_cruise * v_cruise ** 2)
            ld = 1 / tw
            return ld

        def operation_cost(n_hops, mission_time):
            fixed_cost = 2500
            n_workers = 2
            amount_farms = 40
            days_work = np.ceil(amount_farms / n_hops)
            salary = 80
            supervision_cost = n_workers * (days_work * salary * (mission_time / 3600))
            operation_cost = fixed_cost + supervision_cost
            return operation_cost

        def battery_cost(energy_wh):
            return energy_wh

        '''
        def maintenance_cost(n_hops, airframe_cost):
            n_flights = np.ceil(80 / n_hops)
            return n_flights * airframe_cost / 2
        '''

        self.rho_cruise = calc_isa(h_cruise)[2]
        self.wing_loading_max = wing_loading_stall(v_stall, cl_max)

        # -------- Plotting ------------
        wingloading = np.ones(5)
        wingloading = wingloading * self.wing_loading_max
        begin = 10
        end = self.wing_loading_max * 1.1
        ws = np.arange(begin, end)
        a = required_power_ceiling(ws)
        b = required_power_flight(self.v_max, ws)
        cruise = required_power_flight(self.v_cruise, ws)
        c = required_power_rate_of_climb(ws)
        if VTOL:
            d = required_power_vtol(ws)

        '''
        Design requirements:

        prop: design towards top right

        above ceiling req
        below vmax req
        left of stall req
        '''
        if self.plot:
            plt.plot(ws, a, label='ceiling req')
            plt.plot(ws, b, label='vmax req')
            plt.plot(ws, cruise, label='cruise power needed, not design limiting')
            plt.plot(ws, c, label='roc req')
            if VTOL:
                plt.plot(ws, d, label='vtol req')
            plt.plot(wingloading, np.linspace(0, 3, 5), label='stall req')
            plt.xlabel('W/S')
            plt.ylabel('W/P')
            plt.legend()
            plt.grid()
            plt.title(f'W/S-W/P of {self.name}')
            plt.show()

        # ------- Design Selection --------

        if optimisation:
            self.wing_loading_design = self.wing_loading_max
            self.power_loading_design_cruise = required_power_flight(self.v_cruise, self.wing_loading_design)
            if VTOL:
                self.power_loading_design_vtol = required_power_vtol(self.wing_loading_design)
        else:
            self.wing_loading_design = float(input('Choose design point W/S: '))
            self.power_loading_design_cruise = float(input('Choose design point W/P cruise: '))
            if VTOL:
                self.power_loading_design_vtol = float(input('Choose design point W/P vtol: '))

        self.wing_surface = self.mtom * g / self.wing_loading_design
        self.power_required_cruise = self.mtom * g / self.power_loading_design_cruise
        self.battery_mass_vtol = 0
        self.t_mission = self.t_cruise

        self.energy_required_payload = self.power_consumption_payload * self.t_mission / 3600
        self.energy_required_flight_controller = self.power_consumption_flight_controller * self.t_mission / 3600
        self.energy_required_cruise = self.power_required_cruise * self.t_cruise / 3600 * 1.25

        self.energy_required = self.energy_required_cruise + self.energy_required_payload + self.energy_required_flight_controller

        if VTOL:
            self.power_required_vtol = self.mtom * g / self.power_loading_design_vtol
            self.t_hop = self.h_hop / self.roc_vtol
            self.t_hops = self.n_hops * self.t_hop * 2.5
            self.energy_required_vtol = self.power_required_vtol * self.t_hops / 3600
            self.battery_mass_vtol = self.energy_required_vtol / energy_density
            self.t_mission += self.t_hops
            self.energy_required += self.energy_required_vtol

        self.battery_mass_cruise = self.energy_required_cruise / energy_density
        self.lift_drag_cruise = lift_drag(self.wing_loading_design)

        self.battery_mass_payload = self.energy_required_payload / energy_density

        self.battery_mass_flight_controller = self.energy_required_flight_controller / energy_density

        self.battery_mass_total = self.battery_mass_vtol + self.battery_mass_cruise + self.battery_mass_flight_controller

        if VTOL:
            if plot_nhops:
                plt.plot(n_hops,
                         self.battery_mass_vtol + self.battery_mass_payload + self.battery_mass_flight_controller,
                         label=f'{self.name}')
                plt.xlabel('nhops')
                plt.ylabel('battery mass')
                plt.legend()
                plt.title(f'battery mass in function of number of hubs for {self.name}')
                plt.grid()
                plt.show()

        self.operation_cost = operation_cost(self.n_hops, self.t_mission)
        self.battery_cost = battery_cost(self.energy_required)
        # self.maintenance_cost = maintenance_cost(self.n_hops, self.airframe_cost)

        self.total_cost = self.operation_cost


rho_0 = 1.225
g = 9.81
d_flight = 80000
energy_density = 200


def weight_jan():
    weight = [[], []]
    hops = []
    cost = [[], []]
    W_montailsitter = 0
    W_monpuffin = 0
    W_monbibrid = 0
    mtailsitter = 60  # originial mass estimation
    m_puffin = 10
    m_bibrid = 10
    n_hops_tailsitter = 15
    n_hops_bibrid = 10

    for i in range(1, 41):
        n_hops = i
        hops.append(i)
        for i in range(100):  # converge mtom
            tailsitter = Configuration('Tailsitter', 15, 30, 1.4, 0.02, 0.5, 0.8, 0.75, 10, 2000, 5000, mtailsitter, 2,
                                       True, 1, 15, n_hops, 0.65, 6.5)
            W_montailsitter = (
                                          tailsitter.battery_mass_payload + tailsitter.battery_mass_vtol + tailsitter.battery_mass_flight_controller) * 9.81
            mtailsitter = \
                eds.sizing(0.8, energy_density * 3600, d_flight, tailsitter.lift_drag_cruise, 0.05, g, 8 * 9.81,
                           W_montailsitter, 0.5 * 9.81, eds.our_a, eds.our_b, 1, False)[4] / 9.81
            costtailsitter = tailsitter.total_cost

            puffin = Configuration('Puffin', 15, 30, 1.4, 0.02, 0, 0.8, 0.75, 10, 2000, 5000, m_puffin, 2, False, 1, 0,
                                   1, 2, 6.5)
            W_monpuffin = (
                                  puffin.battery_mass_payload + puffin.battery_mass_vtol + puffin.battery_mass_flight_controller) * 9.81
            m_puffin = \
                eds.sizing(0.8, energy_density * 3600, d_flight, puffin.lift_drag_cruise, 0.05, g, 8 * 9.81,
                           W_monpuffin,
                           0.5 * 9.81, eds.our_a, eds.our_b, 1, False)[4] / 9.81
            costpuffin = puffin.total_cost

            bibrid = Configuration('Hybrid', 15, 30, 1.4, 0.04, 0.75, 0.8, 0.75, 10, 2000, 5000, m_bibrid, 2, True, 1,
                                   15, n_hops * 2, 0.65, 6.5)
            W_monbibrid = (
                                  bibrid.battery_mass_payload + bibrid.battery_mass_vtol + bibrid.battery_mass_flight_controller) * 9.81
            m_bibrid = \
                eds.sizing(0.8, energy_density * 3600, d_flight, bibrid.lift_drag_cruise, 0.05, g, 2.5 * 9.81,
                           W_monbibrid,
                           0.5 * 9.81, eds.our_a, eds.our_b, 1, False)[4] / 9.81
            cost_bibrid = bibrid.total_cost

        cost[0].append(costtailsitter)
        cost[1].append(cost_bibrid)
        weight[0].append(mtailsitter)
        weight[1].append(m_bibrid)

    # print(W_montailsitter / 9.81, W_monbibrid / 9.81)

    plt.plot(hops, weight[0], label='tail')
    plt.plot(hops, weight[1], label='bibrid')
    plt.legend()
    plt.show()

    return weight


if __name__ == '__main__':
    v_stall = np.linspace(10, 25)
    wing_surface = []
    for v in v_stall:
        puffin = Configuration('puffin', v, 20, 1.356, 0.02, 0.5, 0.8, 0.8, 12, 4000, 6000, 14.352, 1, False, 1, 15, 20, 1, 6.5)
        wing_surface.append(puffin.wing_surface)

    plt.plot(v_stall, wing_surface)
    plt.show()

