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
        'payload_weight': 1.5*9.81,  # N, payload weight
        'monitoring_weight': 2.67*9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0.5*9.81,  # N, autonomous weight
        'n_p': 0.8,  # -, propulsive efficiency
        'LD': 15.37,  # -, lift to drag ratio
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
        'payload_weight': 8*9.81,  # N, payload weight
        'monitoring_weight': 5.25727 * 9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0.5*9.81,  # N, autonomous weight
        'n_p': 0.8,  # -, propulsive efficiency
        'LD': 15,  # -, lift to drag ratio
        'non_cB': 0.05,  # Non cruise energy consumption

        'R': 80000,  # m, flight range
        # Constants
        'E_d': 200*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': our_a,  # -, regression constant for small RC Uav
        'b': our_b,  # -, regression constant for small RC Uav

        # Other
        'con': 1.05  # -, Contingency for actuating mechanisms, other features
    },
    'MAV': {
        # Payload parameters
        'payload_weight': 2.5*9.81,  # N, payload weight
        'monitoring_weight': 4.599* 9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0.5*9.81,  # N, autonomous weight
        'n_p': 0.8,  # -, propulsive efficiency
        'LD': 8.9,  # -, lift to drag ratio
        'non_cB': 0.05,  # Non cruise energy consumption

        'R': 80000,  # m, flight range
        # Constants
        'E_d': 200*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': our_a,  # -, regression constant for small RC Uav
        'b': our_b,  # -, regression constant for small RC Uav

        # Other
        'con': 1.05  # -, Contingency for actuating mechanisms, other features
    }
}
#
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
    for name, params in sizing_data.items():
        print(f'==== START {name} sizing ====\n\t----  INPUTS ----')
        for param in params:
            print(f"{param: >25}: {params[param]}")
        print(f'\t---- OUTPUTS ----')
        sizing(**params)
        print(f'====   END {name} sizing  ====')


