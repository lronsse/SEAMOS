import numpy as np
import matplotlib.pyplot as plt
#from new_underwater_drag import RF_flat
#from new_underwater_drag import Hull_form_drag
import math



def force_vs_time(time, V, alpha, rho_water, rho_B, C_d, m):
    F_lst = []
    t_lst = []
    R_cone = 0.1
    H_cone = 0.1
    g = 9.81
    pressure_lst = []
    dt=0.001
    for t in np.arange(0, time, dt):
        #F = (V * ((4 / 3) * rho_water * 3 * (t ** 2) * (((4 / pi) * V * np.tan(pi / 6)) ** 3))) * sin(alpha)
        #F = -m*g + (pi/2)*rho_water*C_d*(R_cone**2)* (V**2)*np.tanh(pi/9) + pi*rho_water*g*(R_cone**2)*(V*t-(2/3)*0.5)
        #F = (0.5*rho_water*C_d*pi*((V*t*tan(pi/9))**2)) - (g*((pi*V*t*tan(pi/9))**2)*V*t*(rho_water-rho_B))
        F = (4*np.sqrt(2))*rho_water*(V**(5/2))*((2*H_cone/(R_cone**2))**(-3/2))*(t**0.5)
        r=np.sqrt(V*time/(2*H_cone/R_cone**2))
        surface=V*dt*2*np.pi*r
        pressure =F/surface
        F_lst.append(F)
        t_lst.append(t)
        pressure_lst.append(pressure)

    return F_lst, t_lst,pressure_lst

#F_lst_time, t_lst = force_vs_time(0.3, 4, 0.912, 1030, 228.6, 0.3, 16)
def Paraboloid(x):
    return 0.1*x**2+10


def Paraboloid2(x):
    return 0.03*x**2+10


def Velocity_change(v_initial,mass,body_length,nose_length,r0):
    rho_water=1030
    rho_W=rho_water
    volume=0
    h0=nose_length
    mu_water = 0.00126
    eta_hull=5
    velocity=v_initial
    velocity_tab=[]
    time=0
    for x in np.linspace(0.001,nose_length,200):
        buoyancy=volume*rho_water*9.81
        r=r0/np.sqrt(h0)*np.sqrt(x)
        volume=np.pi*x*r**2/2
        #drag=surface*velocity*length*length/diameter
        S_submerged = (np.pi/6)*(r/x**2)*((r**2 + 4*x**2)**1.5 - r**3)
        ReW = rho_W * velocity * x / mu_water
        #print(ReW,'reynolds')
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2*r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact = (4 * np.sqrt(2)) * rho_water * (velocity ** (5 / 2)) * ((2 * h0 / (r0 ** 2)) ** (-3 / 2)) * (time ** 0.5)
        drag=RF_flat+RF_form
        dt=nose_length/200/velocity
        time=time+dt
        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        print(x,buoyancy,drag,impact)
        print(velocity)
        velocity_tab.append(velocity)
    for x in np.linspace(nose_length,body_length,200):
        buoyancy=volume*rho_water*9.81
        r=r0/np.sqrt(h0)*np.sqrt(x)
        volume=np.pi*x*r**2/2
        #drag=surface*velocity*length*length/diameter
        S_submerged_nose = (np.pi/6)*(r0/h0**2)*((r0**2 + 4*h0**2)**1.5 - r0**3)
        S_sumberged = S_submerged_nose+2*np.pi*r0*(x-nose_length)
        ReW = rho_W * velocity * x / mu_water
        #print(ReW,'reynolds')
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2*r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact = (4 * np.sqrt(2)) * rho_water * (velocity ** (5 / 2)) * ((2 * h0 / (r0 ** 2)) ** (-3 / 2)) * (time ** 0.5)
        drag=RF_flat+RF_form
        dt=nose_length/200/velocity
        time=time+dt
        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        print(x,buoyancy,drag,impact)
        print(velocity)
        velocity_tab.append(velocity)
    print(velocity_tab)
    plt.plot(np.linspace(0.001,body_length,400),velocity_tab)

    plt.xlabel('x location')
    plt.ylabel('Velocity')
    plt.show()

    return velocity_tab


#def Force_versus_X(velocity_tab):


F_lst_time, t_lst,pressure_lst_time = force_vs_time(0.3, 4, 0.912, 1030, 228.6, 0.3, 16)
fig, (Ftime, Fangle) = plt.subplots(2)
#fig.suptitle('Diving load impacts')
Ftime.plot(t_lst, F_lst_time, 'tab:orange')
Ftime.set_xlabel('Time [s]')
Ftime.set_ylabel('Force [N]')
Ftime.set_title('Impact force vs. time')
Fangle.plot(t_lst, pressure_lst_time, 'tab:red')
Fangle.set_xlabel('Time [s]')
Fangle.set_ylabel('Pressure [Pa]')
Fangle.set_title('Pressure vs time')


plt.show()

#def Beam_deflection(angle,Inertia):

#y=0.03x^{2}
#ef Velocity_change(v_initial,mass,body_length,nose_length,r0):
print(np.arctan(0.6)*180/np.pi)
Velocity_change(15,17,0.85,0.18,0.09)
#dimatere 10cm water prop
#body length 85 cm
#11cm
#2cm +tail length
#18cm fuselage
#20 cm
#0.0105 m2 wing

#0.5986wing location
#1.5m