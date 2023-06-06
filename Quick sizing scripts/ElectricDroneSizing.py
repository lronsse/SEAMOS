"""
Adapted from Soham's implementation
"""
import scipy.optimize as optimise
import numpy as np
import matplotlib.pyplot as plt
import wing_sizing as ws

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

    return WplW, WbW+WmbW, WaW, WeW, W_TO


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


def system_cost(name, m_batt, W_TO):
    if name == 'Puffin':
        C_fix = 18220
    if name == 'MultiSystem':
        C_fix = 26330
    if name == 'MAV':
        C_fix = 24730

    if name == 'Puffin':
        PW = 0.4
    else:
        PW = 0.1
    C_prop = PW*W_TO*0.0426+209.15

    C_batt = 200*m_batt

    return C_prop + C_batt + C_fix

def ops_cost(t_mission, n_units):
    t_mission /= 60  # convert to hours
    CpH = 80  # Eur per hour
    t_add_day = 2  # addtl. hours per day
    t_overhead_day = 5  # charged per day
    T_workday = 10  # length of a workday

    n_missions = (40/n_units)
    max_missions_day = (T_workday-t_add_day)/t_mission

    # print(max_missions_day, t_mission)
    n_days = n_missions/max_missions_day

    t_tot_month = n_days*(t_add_day+t_overhead_day)+t_mission*n_missions
    t_tot_year = t_tot_month*12

    #print(t_tot_month)
    return(t_tot_year*CpH)


if __name__ == '__main__':
    data = {
        'Puffin': [sizing_data['Puffin'].copy(), [62.56410256410256, 71.69656773715104, 80.82903291019953, 89.961498083248, 99.09396325629649, 108.22642842934496, 117.35889360239345, 126.49135877544192, 135.62382394849038, 144.75628912153888, 153.8887542945874, 163.02121946763583, 172.15368464068433, 181.28614981373283, 190.41861498678128, 199.55108015982978, 208.68354533287825, 217.81601050592673, 226.94847567897523, 236.08094085202367, 245.21340602507217, 254.34587119812065, 263.4783363711691, 272.6108015442176, 281.7432667172661, 290.87573189031457, 300.008197063363, 309.1406622364116, 318.27312740946, 327.4055925825085, 336.53805775555696, 345.67052292860546, 354.8029881016539, 363.93545327470235, 373.0679184477509, 382.20038362079936, 391.33284879384786, 400.4653139668963, 409.5977791399448, 418.73024431299325], [0.0, 0.3951923715225416, 0.7903847430450832, 1.1855771145676248, 1.5807694860901664, 1.9759618576127078, 2.3711542291352496, 2.766346600657791, 3.161538972180333, 3.5567313437028747, 3.9519237152254156, 4.347116086747958, 4.742308458270499, 5.137500829793041, 5.532693201315582, 5.927885572838124, 6.323077944360666, 6.718270315883207, 7.113462687405749, 7.50865505892829, 7.903847430450831, 8.299039801973374, 8.694232173495916, 9.089424545018458, 9.484616916540999, 9.87980928806354, 10.275001659586081, 10.670194031108622, 11.065386402631164, 11.460578774153706, 11.855771145676249, 12.25096351719879, 12.646155888721331, 13.041348260243874, 13.436540631766414, 13.831733003288956, 14.226925374811499, 14.62211774633404, 15.01731011785658, 15.412502489379122]],
        'MultiSystem': [sizing_data['MultiSystem'].copy(), [70.625, 74.58333333333333, 78.54166666666667, 82.5, 86.45833333333333, 90.41666666666667, 94.375, 98.33333333333333, 102.29166666666667, 106.25, 110.20833333333333, 114.16666666666667, 118.125, 122.08333333333333, 126.04166666666667, 130.0, 133.95833333333334, 137.91666666666666, 141.875, 145.83333333333334, 149.79166666666666, 153.75, 157.70833333333334, 161.66666666666666, 165.625, 169.58333333333334, 173.54166666666666, 177.5, 181.45833333333334, 185.41666666666666, 189.375, 193.33333333333334, 197.29166666666666, 201.25, 205.20833333333334, 209.16666666666666, 213.125, 217.08333333333334, 221.04166666666666, 225.0], [0.2300160399776165, 0.460032079955233, 0.6900481199328494, 0.920064159910466, 1.1500801998880825, 1.3800962398656988, 1.6101122798433156, 1.840128319820932, 2.0701443597985483, 2.300160399776165, 2.5301764397537814, 2.7601924797313977, 2.9902085197090145, 3.2202245596866312, 3.450240599664247, 3.680256639641864, 3.9102726796194798, 4.1402887195970965, 4.370304759574713, 4.60032079955233, 4.830336839529946, 5.060352879507563, 5.290368919485179, 5.520384959462795, 5.750400999440412, 5.980417039418029, 6.210433079395645, 6.4404491193732625, 6.6704651593508775, 6.900481199328494, 7.130497239306111, 7.360513279283728, 7.5905293192613446, 7.8205453592389595, 8.050561399216578, 8.280577439194193, 8.51059347917181, 8.740609519149427, 8.970625559127043, 9.20064159910466], ws.weight_jan()[0]],
        'MAV': [sizing_data['MAV'].copy(), [68.41666666666667, 70.16666666666667, 71.91666666666667, 73.66666666666667, 75.41666666666667, 77.16666666666667, 78.91666666666667, 80.66666666666667, 82.41666666666667, 84.16666666666667, 85.91666666666667, 87.66666666666667, 89.41666666666667, 91.16666666666667, 92.91666666666667, 94.66666666666667, 96.41666666666667, 98.16666666666667, 99.91666666666667, 101.66666666666667, 103.41666666666667, 105.16666666666667, 106.91666666666667, 108.66666666666667, 110.41666666666667, 112.16666666666667, 113.91666666666667, 115.66666666666667, 117.41666666666667, 119.16666666666667, 120.91666666666667, 122.66666666666667, 124.41666666666667, 126.16666666666667, 127.91666666666667, 129.66666666666666, 131.41666666666666, 133.16666666666666, 134.91666666666666, 136.66666666666666], [0.39751427149749313, 0.7950285429949863, 1.1925428144924795, 1.5900570859899725, 1.9875713574874654, 2.385085628984959, 2.782599900482452, 3.180114171979945, 3.577628443477438, 3.9751427149749308, 4.372656986472424, 4.770171257969918, 5.16768552946741, 5.565199800964904, 5.962714072462396, 6.36022834395989, 6.757742615457382, 7.155256886954876, 7.552771158452369, 7.9502854299498615, 8.347799701447356, 8.745313972944848, 9.142828244442342, 9.540342515939836, 9.937856787437328, 10.33537105893482, 10.732885330432314, 11.130399601929808, 11.527913873427302, 11.925428144924792, 12.322942416422286, 12.72045668791978, 13.117970959417272, 13.515485230914765, 13.912999502412259, 14.310513773909753, 14.708028045407245, 15.105542316904739, 15.50305658840223, 15.900570859899723], ws.weight_jan()[1]]
    }

    for name, data in data.items():
        template = data[0]
        t_mission_arr = data[1]
        m_monitor_arr = data[2]
        if name != 'Puffin':
            m_TO = data[3]
        n_arr = []
        C_sys_arr = []
        C_ops_arr = []
        C_tot_arr = []
        for i in enumerate(t_mission_arr):
            n_units = i[0]+1
            t_mission = i[1]
            m_monitor = m_monitor_arr[i[0]]
            template['monitoring_weight'] = m_monitor*9.81
            WplW, WbW, WaW, WeW, W_TO = sizing(**template, verbose=False)

            if name != 'Puffin':
                W_TO = m_TO[i[0]] * 9.81

            m_bat = WbW*W_TO/9.81
            C_sys = system_cost(name, m_bat, W_TO)
            C_ops = ops_cost(t_mission, n_units)
            C_tot = C_sys+C_ops

            C_sys_arr.append(C_sys)
            C_ops_arr.append(C_ops)
            C_tot_arr.append(C_tot)
            n_arr.append(n_units)

        print(min(C_sys_arr))
        plt.plot(n_arr, C_sys_arr, label='System Cost', color='green')
        plt.plot(n_arr, C_ops_arr, label='Operational Cost', color='blue')
        plt.plot(n_arr, C_tot_arr, label='Total Cost', color='red')
        C_min = min(C_tot_arr)
        n_opt = C_tot_arr.index(C_min) + 1
        plt.axhline(C_min, color='black', alpha=0.4, label='Min. Total Cost')
        plt.axvline(n_opt, color='black', alpha=0.4, label='Opt. n Units')
        plt.xlabel('Number of Units [-]')
        plt.ylabel('Cost [€]')
        if name == 'MAV':
            name = 'Hybrid'
        plt.title(name)
        plt.legend()
        plt.show()




