import numpy as np
import matplotlib.pyplot as plt
#from new_underwater_drag import RF_flat
#from new_underwater_drag import Hull_form_drag
import math


def Velocity_change(v_initial,mass,body_length,nose_length,r0,wing_location,wing_length,root_chord,airfoil_height,K,beginning_folded,airfoil_area,tail_length):
    rho_water=1030
    rho_W=rho_water
    volume=0
    h0=nose_length
    mu_water = 0.00126
    eta_hull=5
    velocity=v_initial
    velocity_tab=[]
    time=0
    drag_tab=[]
    impact_tab=[]
    buoyancy_tab=[]
    surface_tab=[]
    cone_rad=0.025
    volume_tab=[]
    water_drag_tab=[]
    blade_width_tab=[]

    impact_blade_tab=[]
    alpha=0
    alpha=30*np.pi/180
    print('nose')
    blade_width=0.015
    x_nose_to_fuselage=np.linspace(0.001,nose_length,200)
    for x in x_nose_to_fuselage:
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
        drag=RF_flat#+RF_form
        dt=nose_length/200/velocity
        time=time+dt
        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        r_blade=0.025+x/np.tan(alpha)

        impact_blade=4/3*rho_water*(4/np.pi*cone_rad)**2*velocity**2*12/np.pi
        cone_rad=0.025+x/np.tan(alpha)
        impact_blade=blade_width/(2*np.pi*r_blade)*impact_blade

        if x<0.21*np.sin(alpha):
            impact_blade_tab.append(impact_blade)
            blade_width_tab.append(blade_width)
            blade_width = blade_width - 0.015 / 145
            print(cone_rad, blade_width)
        #print(x,buoyancy,drag,impact)
        #print(velocity)
        water_drag_tab.append(drag)
        drag_tab.append(buoyancy+drag+impact)
        velocity_tab.append(velocity)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)
    print(len(impact_blade_tab))
    x_blade=np.linspace(0,0.21,len(impact_blade_tab))
    plt.plot(x_blade,impact_blade_tab)
    plt.xlabel('Location on the blade (0 being the connection point) [m]')
    plt.ylabel('Maximal impact force[N]')
    plt.title('Impact vs the location on the blade ')
    plt.show()
    print(np.max(impact_blade_tab),'max impact')
    return impact_blade_tab
K=0.935
z=Velocity_change(15,17,1.05,0.18,0.09,0.5986,1.5,0.357,0.12*0.357,K,0.5986-0.1,0.105,0.258)


def Triangular_beam_failure(length, base, height, loads):
    elastic_modulus = 230*1000000000
    moment_of_inertia = (base * height ** 3) / 36

    centroid_position = height / 3
    max_stress = 0
    i=0
    stress_tab=[]
    for x in np.linspace(0,0.21,len(loads)):
        bending_moment = loads[i] * (length - x)*np.sin(60/180*np.pi)
        i = i + 1
        stress = bending_moment * centroid_position / moment_of_inertia
        stress_tab.append(stress/10**6)
        if stress > max_stress:
            max_stress = stress
    print(max_stress)
    # Determine if the beam will fail
    allowable_stress = elastic_modulus / 1.5  # Assuming a safety factor of 1.5
    print(allowable_stress)
    print(max_stress/allowable_stress,'allowable stress %')
    stress_tab=stress_tab
    plt.plot(np.linspace(0,0.21,len(loads)),stress_tab)
    plt.xlabel('x location on the beam [m]')
    plt.ylabel('Maximal stress [MPa]')
    plt.title('Stress versus location on the beam')
    plt.show()

    if max_stress > allowable_stress:
        print('The beam will fail')
    else:
        return ('The beam will not fail')


# Example usage
beam_length = 0.21  # Length of the beam in meters
base_length = 0.015  # Length of the base of the triangular cross-section in meters
height = 0.003  # Height of the triangular cross-section in meters
loads = z # Varying loads at different locations along the beam in Newtons

result = Triangular_beam_failure(beam_length, base_length, height, loads)
print(result)


