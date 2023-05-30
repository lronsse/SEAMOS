from weather import weather_availability

import numpy as np


def operations_per_mission(n_units, time_per_unit, cruise_speed, turnaround_time):
    cruise_time = cruise_dist / cruise_speed
    total_mission_time = 2 * cruise_time + n_units * time_per_unit
    print(f"            Total Mission time: {np.floor(total_mission_time/3600):1.0f} h {(total_mission_time/3600 - np.floor(total_mission_time/3600)) * 60:2.0f} m")
    amount_of_missions_per_day = int(np.floor((8*3600) / (total_mission_time + turnaround_time)))
    print(f"    Amount of missions per day: {amount_of_missions_per_day}")
    print(f"           Total units per day: {amount_of_missions_per_day * n_units}")
    days_for_all_units = np.ceil(40 / (amount_of_missions_per_day * n_units))
    print(f"      Total days for all units: {int(days_for_all_units)}")
    weather_results = weather_availability((days_for_all_units / 30) * 100)
    print(f"Must design for currents up to: {weather_results[0]:2.2f} m/s")
    print(f"   Must design for waves up to: {weather_results[1]:2.2f} m")
    print(f"    Must design for wind up to: {weather_results[2]:2.2f} m/s")


cruise_dist = 40000  # [m]
n_units_per_mission = 4


operations_per_mission(n_units=n_units_per_mission,
                       cruise_speed=20,  # [m/s]
                       time_per_unit=4.5576*60,  # [s]
                       turnaround_time=20*60
                       )
