"""
Inputs: v_stall, v_cruise, configuration, mtom, n_hops, energy density battery, clmax, rate of climb, h_hop, n_hops, roc_vtol

required inputs: mtom, nhops, thop, rocvtol

assumptions: rocvtol = 1m/s, roccruise = 1

Outputs: Wing surface, battery mass, power/thrust required, lift over drag



"""
import matplotlib.pyplot as plt
import numpy as np

class Configuration:
    def __init__(self, name, v_stall, v_cruise, cl_max, cd_0, cd_vtol, eta_p, e, ar, h_cruise, h_ceiling, mtom, roc, VTOL, n_prop, h_hop, n_hops, roc_vtol):
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
            weight_power = eta_p / (0.5 * rho_0 * v ** 3 * cd_0 / ws + 2 * self.k * ws/ (self.rho_cruise * v))
            return weight_power

        def required_power_ceiling(ws):
            """
            Function determining the weight over power required to fly at service ceiling

            :param ws: wing loading
            :return: weight over power required to fly at max speed
            """
            rho_ceiling = calc_isa(h_ceiling)[2]
            sigma_ceiling = rho_0 / rho_ceiling
            weight_power = sigma_ceiling / (np.sqrt(2 * ws / (rho_ceiling * np.sqrt(3 * cd_0 / self.k))) * 1.155 / (self.lift_drag_max * eta_p))
            return weight_power

        def required_power_rate_of_climb(ws):
            weight_power = 1 / ( self.roc / eta_p + np.sqrt(2 * ws / (self.rho_cruise * np.sqrt(3 * cd_0 / self.k))) * 1.155 / (self.lift_drag_max * eta_p))
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
            tw = 0.5 * self.rho_cruise * v_cruise ** 2 * cd_0 / ws + self.k * ws / (0.5 * self.rho_cruise * v_cruise ** 2)
            ld = 1 / tw
            return ld



        self.rho_cruise = calc_isa(h_cruise)[2]
        self.wing_loading_max = wing_loading_stall(v_stall, cl_max)

        # -------- Plotting ------------
        wingloading = np.ones(5)
        wingloading = wingloading * self.wing_loading_max
        begin = 10
        end = self.wing_loading_max * 1.1
        ws = np.arange(begin, end)
        a = required_power_ceiling(ws)
        b = required_power_flight(self.v_max ,ws)
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
            plt.plot(wingloading, np.linspace(0, 1, 5), label='stall req')
            plt.xlabel('W/S')
            plt.ylabel('W/P')
            plt.legend()
            plt.grid()
            plt.title(f'W/S-W/P of {self.name}')
            plt.show()

        # ------- Design Selection --------



        if optimisation:
            self.wing_loading_design = self.wing_loading_max
            self.power_loading_design_cruise = required_power_flight(self.v_max, self.wing_loading_design)
            if VTOL:
                self.power_loading_design_vtol = required_power_vtol(self.wing_loading_design)
        else:
            self.wing_loading_design = float(input('Choose design point W/S: '))
            self.power_loading_design_cruise = float(input('Choose design point W/P cruise: '))
            if VTOL:
                self.power_loading_design_vtol = float(input('Choose design point W/P vtol: '))

        self.wing_surface = self.mtom * g / self.wing_loading_design
        self.power_required_cruise = self.mtom * g / self.power_loading_design_cruise
        if VTOL:
            self.power_required_vtol = self.mtom * g / self.power_loading_design_vtol
            self.t_hop = self.h_hop / self.roc_vtol
            self.t_hops = self.n_hops * self.t_hop * 2.5
            self.energy_required_vtol = self.power_required_vtol * self.t_hops / 3600
            self.battery_mass_vtol = self.energy_required_vtol / energy_density

        self.energy_required_cruise = self.power_required_cruise * self.t_cruise / 3600

        self.battery_mass_cruise = self.energy_required_cruise / energy_density
        self.lift_drag_cruise = lift_drag(self.wing_loading_design)





rho_0 = 1.225
g = 9.81
d_flight = 80000
energy_density = 150

tailsitter = Configuration('Tailsitter', 15, 30, 1.4, 0.02, 0.5, 0.8, 0.75, 10, 2000, 5000, 60, 2, True, 1, 20, 10, 0.65)
puffin = Configuration('Puffin', 15, 30, 1.4, 0.02, 0, 0.8, 0.75, 10, 2000, 5000, 37.73, 2, False, 1, 20, 0, 0.65)
hybrid = Configuration('Hybrid', 15, 30, 1.4, 0.04, 0.5, 0.8, 0.75, 10, 2000, 6000, 27.09, 2, True, 1, 20, 80, 0.65)


print(f'power vtol: tailsitter = {tailsitter.power_required_vtol}, puffin = {puffin.power_required_cruise}, Hybrid = {hybrid.power_required_vtol}')
print(f'battery mass: tailsitter= {tailsitter.battery_mass_vtol} (tot = {tailsitter.battery_mass_cruise + tailsitter.battery_mass_vtol}), puffin = {puffin.battery_mass_cruise}, hybrid = {hybrid.battery_mass_vtol} (tot = {hybrid.battery_mass_cruise + hybrid.battery_mass_vtol})')
print(f'L/D cruise: tailsitter= {tailsitter.lift_drag_cruise}, puffin = {puffin.lift_drag_cruise}, hybrid = {hybrid.lift_drag_cruise}')
#print(tailsitter.wing_surface, tailsitter.battery_mass_cruise, tailsitter.battery_mass_vtol, tailsitter.wing_loading_max)
#print(puffin.wing_surface, puffin.battery_mass_cruise, puffin.lift_drag_cruise)
#print('test')

'''
v_stall = 10 # m/s
v_max = 20
cl_max = 1.4
rho_0 = 1.225
cd_0 = 0.04
rho_water = 1000
mu_water = 1.3059 * 10 ** -3
L_fus = 1.5
D_fus = 0.5
V_fuselage = L_fus * np.pi * (D_fus / 2) ** 2
S_water = V_fuselage ** (2/3)
rho_c = 0.8 # density at cruise altitude tbd
e = 0.7
b = 1
#S = 1
S_ratio = 1.3 # ratio between total and wing area
eta_p = 0.6 # propulsive efficiency
sigma = rho_c / rho_0 # ratio between cruise and SL density
rc = 2 # rate of climb [m/s]
#ar = b ** 2 / S # aspect ratio
ar = 10 # aspect ratio

n_prop = 1
g = 9.81
m = 60
cd_vtol = 0.1
v_cruise = 20 # m /s
d_flight = 80000 # m
t_flight = d_flight / v_cruise
n_hops = 10
Re_water = L_fus * rc * rho_water / mu_water
c_f = 0.0735 / (Re_water) ** 0.2
K_2 = 8
cte = 8
c_r = 0.00789 / (cte)

def wing_loading_stall(tw):
    return 0.5 * v_stall ** 2 * rho_0 * cl_max + tw * 0 #tw*0 for graphing
print(f'wing loading requirement due to stall = {wing_loading_stall(0)} N/m^2')

wing_loading_req = wing_loading_stall(0)
S = m * g / wing_loading_req

print(f'Wing area = {S} m²')




def wp(ws):
    return eta_p / (0.5 * rho_0 * v_max ** 3 * cd_0 / ws + 2 * k * ws / (rho_c * sigma * v_max))


def thrust_weight_cruise(ws):
    tw = 0.5 * rho_c * v_cruise ** 2 * cd_0 / ws + k * ws / (0.5 * rho_c * v_cruise ** 2)
    return tw

def thrust_weight_vtol(ws):
    tw = 1.2 * (1 + cd_vtol * 0.5 * rho_0 * rc ** 2 * S_ratio / ws)
    pw = tw * rc / eta_p
    treq = tw * m * g
    fm = 0.4742 * treq ** 0.0793
    dl = 3.2261 * m + 74.991
    sp = m * g / (dl * n_prop)
    vh = np.sqrt(treq / (2 * rho_0 * sp))
    vi = (-rc / (2*vh) + np.sqrt((rc / (2*vh))**2 +1)) * vh
    pwreq = treq * vi / fm
    return tw, treq, pwreq


print(f'{1 / thrust_weight_cruise(wing_loading_req)}')

D_water = 0.5 * rho_water * S_water * rc **2 * (c_f + c_r)
Drag = thrust_weight_vtol(wing_loading_req)[1] + D_water

print(f'Power needed = {thrust_weight_vtol(wing_loading_req)[2]} W')
print(f'thrust per prop = {Drag / n_prop} N')

e_density_battery = 150

print(f'battery weight hop {thrust_weight_vtol(21.4375)[2] * 80 / 3600 /150} kg')
print(thrust_weight_vtol(wing_loading_req)[2] * 40 / 3600 / 150)
print(f'battery weight flight {thrust_weight_cruise(wing_loading_req) * m * g * d_flight / 3600 / 150 / 0.6} kg')

m_battery_hop = thrust_weight_vtol(wing_loading_req)[2] * 80 / 3600 / 150 / 0.6
m_battery_flight = thrust_weight_cruise(wing_loading_req) * m * g * d_flight / 3600 / 150

m_battery = m_battery_flight + n_hops * m_battery_hop

print(f'total battery mass = {m_battery}')






plot = False

if plot == True:
    begin = 1
    end = 90

    plt.plot(np.linspace(begin, end), thrust_weight_vtol(np.linspace(begin, end))[0], label='vtol thrust req')
    plt.plot(wing_loading_stall(np.linspace(0, 6)), np.linspace(0, 8), label='stall req')
    plt.plot(np.linspace(begin, end), thrust_weight_cruise(np.linspace(begin, end)), label='cruise req')
    plt.legend()
    plt.grid()
    plt.show()

'''