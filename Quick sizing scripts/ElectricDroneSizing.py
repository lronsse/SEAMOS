"""
Adapted from Soham's implementation
"""
import scipy.optimize as optimise
import matplotlib.pyplot as plt


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
    return W_TO


sizing_data = {
    'Puffin': {
        # Payload parameters
        'payload_weight': .5*9.81,  # N, payload weight
        'monitoring_weight': 5.72*9.81,  # N, monitoring weight

        # Flight parameters
        'autopilot_weight': 0*9.81,  # N, autonomous weight
        'n_p': 0.7,  # -, propulsive efficiency
        'LD': 15.5,  # -, lift to drag ratio
        'non_cB': 0.05,  # Non cruise energy consumption

        'R': 80000,  # m, flight range
        # Constants
        'E_d': 250*3600,  # J/kg, energy density batteries
        'g': 9.81,  # m/s2, gravity constant
        'a': -0.00296,  # -, regression constant for small RC Uav
        'b': 0.87,  # -, regression constant for small RC Uav

        # Other
        'con': 1.3  # -, Contingency for actuating mechanisms, other features
    }
    # 'MultiSystem': {
    #     # Payload parameters
    #     'payload_weight': 7.*9.81,  # N, payload weight
    #     'monitoring_weight': 0. * 9.81,  # N, monitoring weight
    #
    #     # Flight parameters
    #     'autopilot_weight': 0.5*9.81,  # N, autonomous weight
    #     'n_p': 0.6,  # -, propulsive efficiency
    #     'LD': 8,  # -, lift to drag ratio
    #     'non_cB': 0.2,  # Non cruise energy consumption
    #
    #     'R': 80000,  # m, flight range
    #     # Constants
    #     'E_d': 150*3600,  # J/kg, energy density batteries
    #     'g': 9.81,  # m/s2, gravity constant
    #     'a': -0.00296,  # -, regression constant for small RC Uav
    #     'b': 0.87,  # -, regression constant for small RC Uav
    #
    #     # Other
    #     'con': 1.  # -, Contingency for actuating mechanisms, other features
    # },
    # 'MAV': {
    #     # Payload parameters
    #     'payload_weight': .5*9.81,  # N, payload weight
    #     'monitoring_weight': 0. * 9.81,  # N, monitoring weight
    #
    #     # Flight parameters
    #     'autopilot_weight': 0.5*9.81,  # N, autonomous weight
    #     'n_p': 0.7,  # -, propulsive efficiency
    #     'LD': 10,  # -, lift to drag ratio
    #     'non_cB': 0.3,  # Non cruise energy consumption
    #
    #     'R': 80000,  # m, flight range
    #     # Constants
    #     'E_d': 150*3600,  # J/kg, energy density batteries
    #     'g': 9.81,  # m/s2, gravity constant
    #     'a': -0.00296,  # -, regression constant for small RC Uav
    #     'b': 0.87,  # -, regression constant for small RC Uav
    #
    #     # Other
    #     'con': 1.  # -, Contingency for actuating mechanisms, other features
    # }
}


nr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
wr = [0.08697885737100923, 0.38337313601291545, 0.6797674146548216, 0.9761616932967278, 1.2725559719386341, 1.5689502505805402, 1.8653445292224464, 2.1617388078643525, 2.4581330865062587, 2.754527365148165, 3.050921643790071, 3.3473159224319775, 3.6437102010738833, 3.94010447971579, 4.236498758357696, 4.532893036999602, 4.829287315641508, 5.125681594283415, 5.422075872925321, 5.718470151567227, 6.0148644302091325, 6.311258708851039, 6.607652987492946, 6.904047266134852, 7.200441544776758, 7.496835823418664, 7.793230102060571, 8.089624380702476, 8.386018659344384, 8.68241293798629, 8.978807216628196, 9.275201495270101, 9.571595773912009, 9.867990052553916, 10.164384331195821, 10.460778609837726, 10.757172888479634, 11.05356716712154, 11.349961445763446, 11.646355724405351]
resr = []

j = 0
res = []
for i in nr:
    w = wr[j]
    puff_template = sizing_data['Puffin'].copy()
    puff_template['monitoring_weight'] = w*9.81

    res.append(sizing(**puff_template)/9.81)
    j += 1

plt.plot(nr, res)
plt.show()

# if __name__ == '__main__':
#     for name, params in sizing_data.items():
#         print(f'==== START {name} sizing ====\n\t----  INPUTS ----')
#         for param in params:
#             print(f"{param: >25}: {params[param]}")
#         print(f'\t---- OUTPUTS ----')
#         sizing(**params)
#         print(f'====   END {name} sizing  ====')

