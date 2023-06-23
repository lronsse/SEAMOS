import matplotlib.pyplot as plt
import numpy as np

offshore = np.linspace(20000, 60000, 500)
amount_of_units = 40
cruise_speed = 20
time_per_unit = 8 * 1.25
total_time = 8 * 60

fixed_times = 2 + 2 + 1 + 2 + 5 + 5 + 3

unit_time = time_per_unit * amount_of_units - 2
cruise_time = ((offshore * 2) / cruise_speed) / 60 * 1.25

print((cruise_time + unit_time + fixed_times) * 1.2 / 60)



n_farms = np.round_((total_time - fixed_times - cruise_time) / (time_per_unit))
n_farms = np.where(n_farms > 40, 40, n_farms)
print(n_farms)

plt.plot(offshore, n_farms)
plt.title('Distance vs Number of Farms')
plt.ylabel('Number of Units')
plt.xlabel('Distance from Shore')
plt.grid()
plt.show()
