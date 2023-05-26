"""
Sizing of puffin mass and energy
"""
import numpy as np
import scipy.optimize as optimize

ops_time = []
m_mon = []
t_cruise = 30*2

""" Parameters """
# Route
v_current = 0.65  # m/s, current velocity
l_unit = 100.  # m, longitudinal unit length
w_space = l_unit  # m, lateral unit spacing
n_units = 10  # -, no. of units surveyed
# Underwater drag
S_sub = 2.  # m2, Submarine reference area
rho_water = 998.  # kg/m3, water density
Cd_sub = .82  # -, submarine Cd
mu_prop = .8  # -, propulsive efficiency
# Sensor
P_sensors = 30.  # W, sensor average power
# Reference motor params
p_ref = 1035.  # W
m_ref = 0.133  # kg
# Battery
e_dens_batt = 150*60*60  # J/kg


""" Calculate optimal transfer velocity for this parameter combination"""
# Define optimisation function
ret = []


def E_like_min(v_try):
    def time_implicit(t_transfer, v_prop=v_try):
        s_trans = np.sqrt(w_space ** 2 + (l_unit + v_current * t_transfer) ** 2)
        t2 = s_trans/v_prop
        return abs(t2 - t_transfer)

    t_transfer = optimize.minimize_scalar(time_implicit, bounds=(0., 10000.), method='bounded').x
    s_trans = np.sqrt(w_space**2 + (l_unit+v_current * t_transfer)**2)
    E = v_try**2 * s_trans  # energy-like
    ret.append([t_transfer, s_trans, v_try])
    return E


v_trans_opt = optimize.minimize_scalar(E_like_min, bounds=(v_current, 200.), method='bounded').x
t_trans, s_trans, _ = ret[-1]  # please dont hate me for this implementation <3
print(f'Optimal underwater transfer velocity: {v_trans_opt} [m/s]')
print(f'Underwater transfer duration: {t_trans} [s]')
print(f'Underwater transfer distance: {s_trans} [m]')

""" Define velocity relative to unit """
v_unit = v_current
t_mission = l_unit/v_unit * n_units + t_trans * (n_units-1)
print(f'{n_units} unit survey mission duration: {t_mission/60} [min]')


""" Evaluate dragging velocity"""
# unit surveying
dv_unit = v_unit-v_current
# unit transfer
dv_trans = v_trans_opt

""" Evaluate propulsive energy demand"""


def E_prop(dv, s):
    return .5 * Cd_sub * rho_water * S_sub * dv**2 * s * mu_prop


E_survey = n_units * E_prop(dv_unit, l_unit)
E_trans = (n_units - 1) * E_prop(dv_trans, s_trans)
print(f'Monitoring energy demand: {E_survey/10**6} [MJ]')
print(f'Transition energy demand: {E_trans/10**6} [MJ]')

""" Evaluate sensor energy demand"""
print(f'Sensor energy demand: {t_mission*P_sensors/10**6} [MJ]   NOT PROPAGATED TO BATTERY WEIGHT')


""" Evaluate propulsive power demand"""
P_prop_req = E_prop(dv_trans, s_trans)/s_trans * dv_trans / mu_prop
print(f'Required propulsive power: {P_prop_req} [W]')
pmR = p_ref/m_ref
m_motor = P_prop_req/pmR
print(f'Extrapolated motor mass: {m_motor} [kg]')


""" Final payload contributions"""
m_batt_survey = (E_survey+E_trans)/e_dens_batt
print(f'Monitoring battery mass: {m_batt_survey} [kg]')
print(f'Addtl. monitoring mass: {m_motor+m_batt_survey} [kg]')
m_mon.append(m_batt_survey)
ops_time.append(t_mission/60+t_cruise)

print(m_mon)
print(ops_time)
