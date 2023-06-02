import numpy as np


#Estimate power, battery mass and cost for static  (ballast tanks) and dynamic (thrusters) buoyancy control systems

#def thrust
#def thrust battery mass
#def tank battery mass
#def cost/kg idk

rho = 1023.6  # [kg/m^3]
g = 9.81  # [m/s^2]
thruster_thrust = 1.3 * 9.81  # [N]
thruster_power = 80  # [W]
thruster_m = 0.210 # [kg] (for one unit)
mission_t = 1.5  # [h]
Ed = 200  # [Wh/kg]
sys_m = 16  # [kg] System mass
a = 0.3  # [m/s^2] depth acceleration
V = 20 / 1000  # submerged volume
rho_tank = 1.2  # kg/dm^3
thickness = 0.003  # [mm]
pump_volt = 16
pump_amp = 2
pump_mass = 0.340  # [kg]

x_cg = np.cbrt(V) / 2  # [m]

#Verification,


# ------------ Thruster Calculations ----------------

def thrust_req(sys_m,a,rho,V):
    T = (sys_m*a - sys_m * g + rho * g * V)
    return T

def thruster_battery_mass(T,thruster_thrust,mission_t,Ed,thruster_m):
    n_t = np.ceil(T / thruster_thrust) #number of thrusters
    thrust_battery_m = (thruster_power*mission_t)/Ed
    total_m=thrust_battery_m + n_t*thruster_m
    thruster_per_kg=(T)/total_m
    return thrust_battery_m,total_m,n_t,T/9.81


print(thruster_battery_mass(thrust_req(sys_m,a,rho,V),thruster_thrust,mission_t,Ed,thruster_m))


# ------------- Ballast Tanks Calculations -----------

#Todo
#Power calc
#Cost calc
#Buoyant force calc

# propeller water take off


def ballast_volume_req(m, a, rho, V):
    aup = - a
    adown = a
    empty_V = (m * g - m * aup) / (rho * g)
    full_V = (m * g - m * adown) / (rho * g)
    tank_size = max(np.abs(V - full_V), np.abs(V-empty_V))
    dim = np.cbrt(tank_size)
    tank_mass = 6 * dim ** 2 * thickness * rho_tank * 1000  # [kg]
    tank_per_kg = (tank_size * 1 * 9.81) / tank_mass
    # pump_watt=pump_amp*pump_volt
    pump_watt= 15 # 15-99 W
    pump_battery_m = (pump_watt * (90 / 60)) / Ed
    total_m = pump_battery_m + tank_mass + pump_mass
    return empty_V * 1000, full_V * 1000, tank_size * 1000,tank_mass, pump_battery_m,total_m




print(ballast_volume_req(1, 1, 1, 1))






#High strength steel (HY80)
#rho = 7.86 [kg/dm3]
#yield strength = 550 [MPa]
#Tensile strength = 207 [GPa]
#Specifc strength = 70

#Acrylic (HY80)
#rho = 1.2 [kg/dm3]
#yield strength = 103 [MPa]
#Tensile strength = 3.1 [GPa]
#Specifc strength = 86