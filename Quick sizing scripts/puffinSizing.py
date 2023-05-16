"""
Sizing of puffin mass and energy
"""
import numpy as np

""" Payload Definition """
m_sens = .8  # kg, mass of avionics, addtl. sensors
p_sens = 45.  # W, power requirement of avionics, addtl. sensors
e_dens_batt = 150  # Wh/kg, battery energy density

""" Flying Route Definition"""
t_fly = 30*60  # s, flight duration, one way

"""Submarine Route Definition """
v_sub = 1.5  # m/s, speed of navigation

n_units = 8  # number of units monitored in one go
mu_route = .7  # -, efficiency of the route
l_unit = 100.  # m
s_route = n_units*l_unit/mu_route  # m

t_route = s_route/v_sub  # s
# current ignored, max current 1 m/s

""" Hydrodynamic Definition """
a = 1.5  # m
b = 0.2  # m
rho_w = 1000.  # kg/m3, water density
S_sub = np.pi*a*b  # m2 sub cross area, ellipse model, iterate pls
Cd_sub = 0.82  # -, Cd of a cylinder

"""Submarine Battery Sizing"""
E_sub = .5 * rho_w * v_sub**2 * Cd_sub * S_sub * s_route
m_batt_sub = E_sub/(60*60)/e_dens_batt
m_batt_sens = (t_route+2*t_fly)*p_sens/(60*60)/e_dens_batt
m_batt_0 = m_batt_sub + m_batt_sens

print(f'Req. battery mass for submarine navigation [kg]: {m_batt_sub}')
print(f'Battery mass excl. airborne propulsion [kg]: {m_batt_0}')
print(f'Route duration [min] {t_route/60} for {n_units} units')


"""Underwater propulsion sizing"""
# https://en.wikipedia.org/wiki/Power-to-weight_ratio
# ElectriFly GPMG5220 brushless DC motor[37]
# assume same PmR
p_ref = 1035.  # W
m_ref = 0.133  # kg
pmR = p_ref/m_ref
mu_prop = .8

p_sub = .5 * rho_w * v_sub**3 * Cd_sub * S_sub
m_prop_sub = p_sub/mu_prop/pmR

m_pyl = m_sens + m_batt_0 + m_prop_sub
print(f'Req. submarine prop. power [W]: {p_sub}')
print(f'Approx. submarine engine mass [kg]: {m_prop_sub}')
print('-------------------')
print(f'Airborne payload mass [kg]: {m_pyl}')

"""... continue class I weight for A/C, finally apply 30% contingency to mass"""









