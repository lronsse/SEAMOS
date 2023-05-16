import matplotlib.pyplot as plt
import numpy as np

v_stall = 10 # m/s
v_max = 20
cl_max = 1.4
rho_0 = 1.225
cd_0 = 0.02
rho_water = 1000
mu_water = 1.3059 * 10 ** -3
L = 1.5
D = 0.5
V_fuselage = L * np.pi * (D/2)** 2
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
k = 1 / (np.pi * e * ar)
n_prop = 1
g = 9.81
m = 10
cd_vtol = 0.1
v_cruise = 20 # m /s
d_flight = 80000 # m
t_flight = d_flight / v_cruise
n_hops = 10
Re_water = L * rc * rho_water / mu_water
c_f = 0.0735 / (Re_water) ** 0.2
K_2 = 8
cte = 8
c_r = 0.00789 / (cte)

def calc_dynamic_pressure(rho, v):
    return 0.5 * rho * v ** 2


def wing_loading_stall(tw):
    return 0.5 * v_stall ** 2 * rho_0 * cl_max + tw * 0 #tw*0 for graphing
print(f'wing loading requirement due to stall = {wing_loading_stall(0)} N/m^2')

wing_loading_req = wing_loading_stall(0)
S = m * g / wing_loading_req

print(f'Wing area = {S} m²')

def thrust_weight_cruise(ws):
    tw = 0.5 * rho_c * v_cruise ** 2 * cd_0 / ws + k * ws / (0.5 * rho_c * v_cruise ** 2)
    return tw

def wp(ws):
    return eta_p / (0.5 * rho_0 * v_max ** 3 * cd_0 / ws + 2 * k * ws / (rho_c * sigma * v_max))


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


print(thrust_weight_cruise(wing_loading_req) * m * g)

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


begin = 1
end = 90
plot = False

if plot == True:
    plt.plot(np.linspace(begin, end), thrust_weight_vtol(np.linspace(begin, end))[0], label='vtol thrust req')
    plt.plot(wing_loading_stall(np.linspace(0, 6)), np.linspace(0, 8), label='stall req')
    plt.plot(np.linspace(begin, end), thrust_weight_cruise(np.linspace(begin, end)), label='cruise req')
    plt.legend()
    plt.grid()
    plt.show()

