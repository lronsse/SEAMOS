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
thickness = 0.005  # [m]
pump_volt = 16
pump_amp = 2
pump_mass = 0.340  # [kg]

x_cg = np.cbrt(V) / 2  # [m]

#Verification,


# ------------ Thruster Calculations ----------------

def thrust_req(sys_m,a,rho,V):
    T = (sys_m*a - sys_m * g + rho * g * V)
    return T


def thruster_battery_mass(T,thruster_thrust,mission_t,Ed,thruster_m,thruster_power):
    n_t = np.ceil(T / thruster_thrust) #number of thrusters
    thrust_battery_m = (thruster_power*mission_t)/Ed
    total_m=thrust_battery_m + n_t*thruster_m
    thruster_per_kg=(T)/total_m
    return thrust_battery_m,n_t,thruster_per_kg,total_m


print(thruster_battery_mass(thrust_req(sys_m,a,rho,V),thruster_thrust,mission_t,Ed,thruster_m,thruster_power))


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
    # R_in = (3*tank_size/(4*np.pi))**(1/3)
    R_in=0.15
    dim = np.cbrt(tank_size)
    tank_mass = ((dim+thickness)**2-(dim**2))*rho_tank
    # tank_mass = 6 * dim ** 2 * thickness * rho_tank * 1000  # [kg]
    tank_mass=(4*np.pi/3)*((R_in+thickness)**3-R_in**3)*rho_tank
    tank_per_kg = (tank_size * 1 * 9.81) / tank_mass
    # pump_watt=pump_amp*pump_volt
    pump_watt= 15 # 15-99 W
    pump_battery_m = (pump_watt * (5 / 60)) / Ed
    total_m = pump_battery_m + tank_mass + pump_mass
    return empty_V * 1000, full_V * 1000, tank_size * 1000,tank_mass, pump_battery_m,total_m



print("ahhhhh")
print(ballast_volume_req(sys_m, a, rho, V))



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



#### Verification


def test_thrust_req():
    expected_T = 108.86
    calculated_T = thrust_req(10, 5, 2, 8)
    assert calculated_T == expected_T

def test_thruster_battery_mass():
    test_T=95
    thrust_battery_m,total_m,n_t,thruster_per_kg=thruster_battery_mass(test_T,10,5,20,6,9)
    expected_n_t=10
    expected_thrust_battery_m = 2.25
    expected_total_m =62.25
    expected_thruster_per_kg = 1.52610442
    assert expected_n_t==n_t
    assert expected_thrust_battery_m == thrust_battery_m
    assert expected_total_m == total_m
    assert expected_thruster_per_kg - thruster_per_kg < 0.01