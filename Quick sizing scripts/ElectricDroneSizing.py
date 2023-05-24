"""
Adapted from Soham's implementation
"""
import scipy.optimize as optimise


def sizing(n_p, E_d, R, LD, non_cB, g, payload_weight, monitoring_weight, autopilot_weight, a, b, con, verbose=True):
    """Class I sizing function, taking unzipped parameter dict"""
    WbW = (1 + non_cB) * ((g / (n_p * E_d)) * (R / LD))  # -, battery weight fraction

    # solve implicit equation for M_TO
    def minimise_func(W_TO_try):  # Try a W_TO and see the match
        WeW = a * W_TO_try * 0.224809 + b  # Get empty weight fraction from stat relations
        return abs(payload_weight + autopilot_weight + monitoring_weight + WbW*W_TO_try + WeW*W_TO_try - W_TO_try)

    W_TO = optimise.minimize_scalar(minimise_func, bounds=(0., 40000.), method='bounded').x
    WplW = payload_weight / W_TO
    WmbW = monitoring_weight / W_TO
    WaW = autopilot_weight / W_TO
    WeW = a * W_TO * 0.224809 + b
    W_TO *= con
    if verbose:
        print(f"   Payload Weight / Take-off Weight: {WplW:1.4f} [-]")
        print(f"   Battery Weight / Take-off Weight: {WbW:1.4f} [-]")
        print(f" Autopilot Weight / Take-off Weight: {WaW:1.4f} [-]")
        print(f"Monitoring Weight / Take-off Weight: {WmbW:1.4f} [-]")
        print(f"     Empty Weight / Take-off Weight: {WeW:1.4f} [-]")
        print(f"                    Take-Off Weight: {W_TO:3.2f} [N]\n")
        print(f"                       Takeoff Mass: {W_TO / 9.81:3.2f} [kg]")
        print(f"                       Payload Mass: {(W_TO*WplW) / 9.81:3.2f} [kg]")
        print(f"                       Battery Mass: {(W_TO*WbW) / 9.81:3.2f} [kg]")
        print(f"                     Autopilot Mass: {(W_TO*WaW) / 9.81:3.2f} [kg]")
        print(f"                    Monitoring Mass: {(W_TO*WmbW) / 9.81:3.2f} [kg]")
        print(f"                         Empty Mass: {(W_TO*WeW) / 9.81:3.2f} [kg]\n")
        consumed_energy = WbW * W_TO * E_d / 9.81
        print(f"          Cruise Energy Consumption: {consumed_energy / (10 ** 6):3.2f} [MJ] or {consumed_energy / (10 ** 3) / 3600:3.2f} [kWh]")
        print(f"                       Battery Cost: {(consumed_energy / 3600) * 1:3.0f} euro")
        print(f"            Monitoring Battery Cost: {(W_TO * WmbW * E_d / 9.81 / 3600) * 1:3.0f} euro")

    return WplW, WbW, WaW, WeW, W_TO


# Calculated from SkyLane UAV (https://sky-drones.com/skylane)
our_a = -0.0013636
our_b = 0.545

sizing_data = {
    'Puffin': {
        # Payload parameters
        'payload_weight': .5*9.81,  # N, payload weight
        'monitoring_weight': 3.6437102010738838*9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0.5*9.81,  # N, autonomous weight
        'n_p': 0.7,  # -, propulsive efficiency
        'LD': 12,  # -, lift to drag ratio
        'non_cB': 0.05,  # Non cruise energy consumption

        'R': 80000,  # m, flight range
        # Constants
        'E_d': 200*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': our_a,  # -, regression constant for small RC Uav
        'b': our_b,  # -, regression constant for small RC Uav

        # Other
        'con': 1.3  # -, Contingency for actuating mechanisms, other features
    },
    'MultiSystem': {
        # Payload parameters
        'payload_weight': 7.*9.81,  # N, payload weight
        'monitoring_weight': 0. * 9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0.5*9.81,  # N, autonomous weight
        'n_p': 0.6,  # -, propulsive efficiency
        'LD': 8,  # -, lift to drag ratio
        'non_cB': 0.2,  # Non cruise energy consumption

        'R': 80000,  # m, flight range
        # Constants
        'E_d': 200*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': our_a,  # -, regression constant for small RC Uav
        'b': our_b,  # -, regression constant for small RC Uav

        # Other
        'con': 1.  # -, Contingency for actuating mechanisms, other features
    },
    'MAV': {
        # Payload parameters
        'payload_weight': .5*9.81,  # N, payload weight
        'monitoring_weight': 0. * 9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0.5*9.81,  # N, autonomous weight
        'n_p': 0.7,  # -, propulsive efficiency
        'LD': 10,  # -, lift to drag ratio
        'non_cB': 0.3,  # Non cruise energy consumption

        'R': 80000,  # m, flight range
        # Constants
        'E_d': 200*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': our_a,  # -, regression constant for small RC Uav
        'b': our_b,  # -, regression constant for small RC Uav

        # Other
        'con': 1.  # -, Contingency for actuating mechanisms, other features
    }
}


def test_sizing_verification():
    test_WplW, test_WbW, test_WaW, test_WeW, test_W_TO = sizing(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0, 0.8, -0.9, 1.0, 1.1, verbose=False)
    assert abs(test_WbW - 33.75) < 1
    assert abs(test_W_TO - 183.538) < 1
    assert abs(test_WeW - -32.759) < 1


def test_sizing_validation():
    # Validation with SkyLane-250: https://sky-drones.com/skylane
    test_WplW, test_WbW, test_WaW, test_WeW, test_W_TO = sizing(n_p=0.7, E_d=250*3600, R=300000, LD=12, non_cB=0.05, g=9.81,
                                                                payload_weight=1.2*9.81, autopilot_weight=0, monitoring_weight=0,
                                                                a=-0.0013636, b=0.545, con=1, verbose=True)
    assert abs(test_W_TO - 15*9.81) < 15


if __name__ == '__main__':
    data = {
        'Puffin': [sizing_data['Puffin'].copy(), [62.56410256410256, 71.69656773715104, 80.82903291019953, 89.961498083248, 99.09396325629649, 108.22642842934496, 117.35889360239345, 126.49135877544192, 135.62382394849038, 144.75628912153888, 153.8887542945874, 163.02121946763583, 172.15368464068433, 181.28614981373283, 190.41861498678128, 199.55108015982978, 208.68354533287825, 217.81601050592673, 226.94847567897523, 236.08094085202367, 245.21340602507217, 254.34587119812065, 263.4783363711691, 272.6108015442176, 281.7432667172661, 290.87573189031457, 300.008197063363, 309.1406622364116, 318.27312740946, 327.4055925825085, 336.53805775555696, 345.67052292860546, 354.8029881016539, 363.93545327470235, 373.0679184477509, 382.20038362079936, 391.33284879384786, 400.4653139668963, 409.5977791399448, 418.73024431299325], [0.0, 0.3951923715225416, 0.7903847430450832, 1.1855771145676248, 1.5807694860901664, 1.9759618576127078, 2.3711542291352496, 2.766346600657791, 3.161538972180333, 3.5567313437028747, 3.9519237152254156, 4.347116086747958, 4.742308458270499, 5.137500829793041, 5.532693201315582, 5.927885572838124, 6.323077944360666, 6.718270315883207, 7.113462687405749, 7.50865505892829, 7.903847430450831, 8.299039801973374, 8.694232173495916, 9.089424545018458, 9.484616916540999, 9.87980928806354, 10.275001659586081, 10.670194031108622, 11.065386402631164, 11.460578774153706, 11.855771145676249, 12.25096351719879, 12.646155888721331, 13.041348260243874, 13.436540631766414, 13.831733003288956, 14.226925374811499, 14.62211774633404, 15.01731011785658, 15.412502489379122]],
        'MultiSystem': [sizing_data['MultiSystem'].copy(), [], []],
        'MAV': [sizing_data['MultiSystem'].copy(), [], []]
    }

    for name, data in data.items():
        template = data[0]
        t_mission_arr = data[1]
        m_monitor_arr = data[2]

        n_arr = []
        C_
        for i in enumerate(t_mission_arr):
            n_units = i[0]+1
            t_mission = i[1]
            m_monitor = m_monitor_arr[i[0]]
            template['monitoring_weight'] = m_monitor*9.81
            WplW, WbW, WaW, WeW, W_TO = sizing(**template)

            C_sys = ...
            C_var = ...
            C_ops = ...




