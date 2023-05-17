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
        print(f"  Payload Weight / Take-off Weight: {WplW:1.4f} [-]")
        print(f"  Battery Weight / Take-off Weight: {WbW:1.4f} [-]")
        print(f"Autopilot Weight / Take-off Weight: {WaW:1.4f} [-]")
        print(f"Monitoring Weight / Take-off Weight: {WmbW:1.4f} [-]")
        print(f"    Empty Weight / Take-off Weight: {WeW:1.4f} [-]")
        print(f"                   Take-Off Weight: {W_TO:3.2f} [N]\n")
        print(f"                      Takeoff Mass: {W_TO / 9.81:3.2f} [kg]")
        print(f"                      Payload Mass: {(W_TO*WplW) / 9.81:3.2f} [kg]")
        print(f"                      Battery Mass: {(W_TO*WbW) / 9.81:3.2f} [kg]")
        print(f"                    Autopilot Mass: {(W_TO*WaW) / 9.81:3.2f} [kg]")
        print(f"                        Empty Mass: {(W_TO*WeW) / 9.81:3.2f} [kg]\n")
        print(f"         Cruise Energy Consumption: {WbW * W_TO * E_d / 9.81 / (10 ** 6):3.2f} [MJ]")
        print(f"                                  : {(WbW * W_TO * E_d / 9.81 / (10 ** 3))/3600:3.2f} [kWh]")
    return WplW, WbW, WaW, WeW, W_TO


sizing_data = {
    'Puffin': {
        # Payload parameters
        'payload_weight': .5*9.81,  # N, payload weight
        'monitoring_weight': 3.6437102010738838*9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0*9.81,  # N, autonomous weight
        'n_p': 0.7,  # -, propulsive efficiency
        'LD': 12,  # -, lift to drag ratio
        'non_cB': 0.05,  # Non cruise energy consumption

        'R': 80000,  # m, flight range
        # Constants
        'E_d': 150*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': -0.00296,  # -, regression constant for small RC Uav
        'b': 0.87,  # -, regression constant for small RC Uav

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
        'E_d': 150*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': -0.00296,  # -, regression constant for small RC Uav
        'b': 0.87,  # -, regression constant for small RC Uav

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
        'E_d': 150*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': -0.00296,  # -, regression constant for small RC Uav
        'b': 0.87,  # -, regression constant for small RC Uav

        # Other
        'con': 1.  # -, Contingency for actuating mechanisms, other features
    }
}



if __name__ == '__main__':
    for name, params in sizing_data.items():
        print(f'==== START {name} sizing ====\n\t----  INPUTS ----')
        for param in params:
            print(f"{param: >25}: {params[param]}")
        print(f'\t---- OUTPUTS ----')
        sizing(**params)
        print(f'====   END {name} sizing  ====')

