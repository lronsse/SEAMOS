"""
Adapted from Soham's implementation
"""
import scipy.optimize as optimise


def sizing(n_p, E_d, R, LD, non_cB, g, W_pl, W_auto, a, b, con):
    """Class I sizing function, taking unzipped parameter dict"""
    WbW = (1 + non_cB) * ((g / (n_p * E_d)) * (R / LD))  # -, battery weight fraction

    # solve implicit equation for M_TO
    def minimise_func(W_TO_try):  # Try a W_TO and see the match
        WeW = a * W_TO_try * 0.224809 + b  # Get empty weight fraction from stat relations
        return abs(W_pl + W_auto + WbW*W_TO_try + WeW*W_TO_try - W_TO_try)

    W_TO = optimise.minimize_scalar(minimise_func, bounds=(0., 40000.), method='bounded').x
    WplW = W_pl / W_TO
    WaW = W_auto / W_TO
    WeW = a * W_TO * 0.224809 + b
    W_TO *= con
    print('Wto', W_TO)
    print('WplW', WplW)
    print('WbW', WbW)
    print('WaW', WaW)
    print('WeW', WeW)
    print('Takeoff Mass [kg]:', W_TO/9.81)
    print('Cruise Energy Consumtpion [MJ]:', WbW*W_TO*E_d/9.81/(10**6))


sizing_data = {
    'Puffin': {
        # Payload parameters
        'W_pl': 3.2*9.81,  # N, payload weight

        # Flight parameters
        'W_auto': 0*9.81,  # N, autonomous weight
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
        'W_pl': 10.*9.81,  # N, payload weight

        # Flight parameters
        'W_auto': 0.5*9.81,  # N, autonomous weight
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
        'W_pl': .5*9.81,  # N, payload weight

        # Flight parameters
        'W_auto': 0.5*9.81,  # N, autonomous weight
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
        print(f'--------- START {name} sizing ---------')
        sizing(**params)
        print(f'--------- END {name} sizing ---------')

