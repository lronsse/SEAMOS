import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def weather_availability(percentage_duration):
    # Standard Deviations and averages
    sd_current = 0.1547684955
    sd_waves = 0.7158412638
    sd_wind = 3.483927377

    mu_current = 0.6676470588
    mu_waves = 1.285751425
    mu_wind = 7.891210879

    # alpha to critical
    alpha = 1 - (percentage_duration / 100)
    n_sided = 1  # 2-sided test
    z_crit = stats.norm.ppf(1 - alpha / n_sided)

    # Result
    x_current = z_crit * sd_current + mu_current
    x_waves = z_crit * sd_waves + mu_waves
    x_wind = z_crit * sd_wind + mu_wind

    return x_current, x_waves, x_wind

rc = []
rwave = []
rwind = []
num=[]

if __name__ == '__main__':


    # result_current, result_waves, result_wind = weather_availability(float(input('Enter percentage of month (write digits e.g. 80): ',)))

    # print('Must design for currents up to', result_current, 'm/s')
    # print('Must design for waves up to', result_waves, 'm')
    # print('Must design for wind up to', result_wind, 'm/s')

    for i in range(100):
        result_current, result_waves, result_wind = weather_availability(i)
        num.append(i)
        rc.append(result_current)
        rwave.append(result_waves)
        rwind.append(result_wind)

    plt.plot(num, rc, label='Currents')
    plt.plot(num, rwave, label='Waves')
    plt.plot(num, rwind, label='Wind')
    plt.legend()
    plt.xlabel("Percentage")
    plt.ylabel("m/s")
    plt.show()



