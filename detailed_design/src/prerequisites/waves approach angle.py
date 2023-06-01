import math
import numpy as np
import matplotlib.pyplot as plt

def landing_wave_angle(approach_angle):
    time=np.linspace(0,40,200)

    approach_thing=-(approach_angle)*np.pi/180
    print(approach_thing)
    wave=np.zeros(len(time))
    wave_diff=np.zeros(len(time))
    landing_wave_angle=np.zeros(len(time))
    approach_line = np.zeros(len(time))
    for i in range(len(time)):
        wave[i] = 2*math.sin(time[i]/(2*np.pi))
        wave_diff[i] = math.cos(time[i]/(2*np.pi))/np.pi
        approach_line[i]=approach_thing
        landing_wave_angle[i]=-math.atan(np.abs((approach_thing-wave_diff[i])/(1+wave_diff[i]*approach_thing)))*180/np.pi
    print(wave_diff)
    print(landing_wave_angle)
    #plt.plot(time,approach_line,label='approach line')
    plt.plot(time, wave,label='wave line')
    #plt.plot(time, wave_diff)
    plt.plot(time,landing_wave_angle,label='angle line')
    plt.legend(loc='lower left')
    plt.show()






#landing_wave_angle(45)
print(np.sqrt(12.5))
print(np.sqrt(0.015/3.53/np.pi))