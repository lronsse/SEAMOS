import numpy as np


#Estimate power, battery mass and cost for static  (ballast tanks) and dynamic (thrusters) buoyancy control systems

#def thrust
#def thrust battery mass
#def tank battery mass
#def cost/kg idk

rho=1000 # [kg/m^3]
g=9.81 # [m/s^2]
thruster_thrust = 1.3 * 9.81 # [N]
thruster_power = 80 #[W]
thruster_m=0.210 # [kg] (for one unit)
mission_t=1.5 #[h]
Ed=200 # [Wh/kg]
sys_m = 16 #[kg] System mass
a=0.03 #[m/s^2] depth acceleration
V=16/1000 #submerged volume

x_cg = np.cbrt(V) / 2 #[m]



# ------------ Thruster Calculations ----------------

def thrust_req(sys_m,a,rho,V):
    T = (sys_m*a - sys_m * g + rho * g * V)
    return T

def thruster_battery_mass(T,thruster_thrust,mission_t,Ed,thruster_m):
    n_t = np.ceil(T / thruster_thrust) #number of thrusters
    thrust_battery_m = (thruster_power*mission_t)/Ed
    total_m=thrust_battery_m + n_t*thruster_m
    return thrust_battery_m,total_m,n_t


print(thruster_battery_mass(thrust_req(sys_m,a,rho,V),thruster_thrust,mission_t,Ed,thruster_m))



# ------------- Ballast Tanks Calculations -----------

def ballast_volume_req(m, a, rho, V):
    aup = - a
    adown = a
    empty_V = (m * g - m * aup) / (rho * g)
    full_V =  (m * g - m * adown) / (rho * g)
    tank_size = max(np.abs(V - full_V), np.abs(V-empty_V))
    return empty_V * 1000, full_V * 1000, tank_size * 1000





print(ballast_volume_req(16, 0.3, 1000, V))

