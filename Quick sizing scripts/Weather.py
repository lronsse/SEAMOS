import numpy as np
from scipy import stats

percentage_duration = float(input('Enter percentage of month (write digits e.g. 80)',))

#Standard Deviations and averages
SD_current = 0.1547684955
SD_waves = 0.7158412638
SD_wind = 3.483927377

mu_current = 0.6676470588
mu_waves = 1.285751425
mu_wind = 7.891210879

# alpha to critical
alpha = 1-(percentage_duration/100)
n_sided = 1 # 2-sided test
z_crit = stats.norm.ppf(1-alpha/n_sided)

#Result

x_current = z_crit*SD_current+mu_current
x_waves = z_crit*SD_waves+mu_waves
x_wind = z_crit*SD_wind+mu_wind

print('Must design for currents up to', x_current, 'm/s')
print('Must design for waves up to', x_waves, 'm')
print('Must design for wind up to', x_wind, 'm/s')



