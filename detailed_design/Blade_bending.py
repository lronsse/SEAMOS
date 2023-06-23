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
    volume_tab=[]
    water_drag_tab=[]
    print('nose')
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

        #print(x,buoyancy,drag,impact)
        #print(velocity)
        water_drag_tab.append(drag)
        drag_tab.append(buoyancy+drag+impact)
        velocity_tab.append(velocity)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)
    print(volume,'Volume nose')
    print('max nose', np.max(impact_tab))
    print('nose to box',x)
    beginning_folded=0.493
    x_nose_to_box=np.linspace(nose_length,beginning_folded,200)
    for x in x_nose_to_box:
        buoyancy=volume*rho_water*9.81
        r=r0/np.sqrt(h0)*np.sqrt(x)
        volume=np.pi*h0*r0**2/2+np.pi*r0**2*(x-nose_length)
        #drag=surface*velocity*length*length/diameter
        S_submerged_nose = (np.pi/6)*(r0/h0**2)*((r0**2 + 4*h0**2)**1.5 - r0**3)
        S_submerged = S_submerged_nose+2*np.pi*r0*(x-nose_length)
        ReW = rho_W * velocity * x / mu_water
        #print(ReW,'reynolds \n')
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2*r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact=0
        #impact = (4 * np.sqrt(2)) * rho_water * (velocity ** (5 / 2)) * ((2 * h0 / (r0 ** 2)) ** (-3 / 2)) * (time ** 0.5)
        drag=RF_flat#+RF_form
        dt=(beginning_folded-nose_length)/200/velocity
        time=time+dt

        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        #print(x,buoyancy,drag,impact)
        #print(velocity)
        water_drag_tab.append(drag)
        velocity_tab.append(velocity)
        drag_tab.append(drag+buoyancy+impact)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        S_nose_full=S_submerged_nose
        volume_tab.append(volume)
        xmax=x



    print('Box')

    # print(np.max(impact_wing_tab),'max_impact_force wings')

    #box_dimensions
    box_width=0.1
    box_height=0.12*0.357
    box_length=0.107
    a=0.0001/2
    b=box_width/2
    S_box=0
    volume_box=0
    velocity_impact=velocity
    x_box=np.linspace(beginning_folded,beginning_folded+box_length,200)
    for x in x_box:
        #print(np.sqrt(1 - (a/b) ** 2), 'e')
        buoyancy=volume*rho_water*9.81
        r=r0/np.sqrt(h0)*np.sqrt(x)
        volume_box=volume_box+a*b*np.pi*velocity*dt
        volume = np.pi * h0 * r0 ** 2 / 2 + np.pi * r0 ** 2 * (x - nose_length) + volume_box
        #drag=surface*velocity*length*length/diameter
        #S_submerged_nose = (np.pi/6)*(r0/h0**2)*((r0**2 + 4*h0**2)**1.5 - r0**3)
        S_submerged_nose=S_nose_full
        #Perimeter=np.pi*(3*a*b-np.sqrt((3*a+b)*(3*b+a)))
        Perimeter=np.pi*2*np.sqrt((a**2+b**2)/2)
        print(Perimeter,'Perimiter')
        S_box=S_box+velocity*dt*Perimeter
        S_submerged = S_submerged_nose + 2 * np.pi * r0 * (x - nose_length)+S_box

        ReW = rho_W * velocity * x / mu_water
        print(ReW,'reynolds',x,velocity)
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2*r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact =velocity**2/2*K*4/3*np.pi*rho_water*a*(b)
        a=a+(box_height/2)/200
        #b=b+(box_width/2)/200

        drag=RF_flat#+RF_form
        dt=(box_length)/200/velocity
        time=time+dt
        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        print(x,buoyancy,drag,impact)
        print(velocity)
        water_drag_tab.append(drag)
        velocity_tab.append(velocity)
        drag_tab.append(drag+buoyancy+impact)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)


    print('aaaaaa buoyancy',buoyancy)
    print('box to wings')
    print('volume box',volume_box)

    print('volume',volume)
    volume_box=0
    wing_location=0.85-0.1
    x_box_to_wings=np.linspace(beginning_folded+box_length,wing_location,200)
    for x in x_box_to_wings:
        print('volume', volume)
        #print(np.sqrt(1 - (a/b) ** 2), 'e')
        buoyancy=volume*rho_water*9.81
        r=r0/np.sqrt(h0)*np.sqrt(x)
        volume_box=volume_box+(x-beginning_folded-box_length)*box_width*box_height
        volume = np.pi * h0 * r0 ** 2 / 2 + np.pi * r0 ** 2 * (x - nose_length) + volume_box
        #drag=surface*velocity*length*length/diameter
        #S_submerged_nose = (np.pi/6)*(r0/h0**2)*((r0**2 + 4*h0**2)**1.5 - r0**3)
        #Perimeter=np.pi*(3*a*b-np.sqrt((3*a+b)*(3*b+a)))

        S_box=S_box+velocity*dt*Perimeter
        S_submerged = S_submerged_nose + 2 * np.pi * r0 * (x - nose_length)+S_box

        ReW = rho_W * velocity * x / mu_water
        print(ReW,'reynolds',x,velocity)
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2*r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact =0
        #a=a+(box_height/2)/200
        #b=b+(box_width/2)/200

        drag=RF_flat#+RF_form
        dt=(wing_location-beginning_folded)/200/velocity
        time=time+dt
        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        print(x,buoyancy,drag,impact)
        print(velocity)
        water_drag_tab.append(drag)
        velocity_tab.append(velocity)
        drag_tab.append(drag+buoyancy+impact)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)
    print('wings smooth bit')



    wing_bit_width=0.25
    a=0.001/2
    b=wing_bit_width/2
    wing_bit_length=0.1

    S_bit=0
    x_wings_smooth_bit=np.linspace(wing_location,wing_location+wing_bit_length,200)
    for x in x_wings_smooth_bit:
        #print('volume', volume)
        #print(np.sqrt(1 - (a/b) ** 2), 'e')
        buoyancy=volume*rho_water*9.81
        r=r0/np.sqrt(h0)*np.sqrt(x)
        #volume_box=volume_box+(x-beginning_folded-box_length)*box_width*box_height
        #volume = np.pi * h0 * r0 ** 2 / 2 + np.pi * r0 ** 2 * (x - nose_length) + volume_box
        #drag=surface*velocity*length*length/diameter
        #S_submerged_nose = (np.pi/6)*(r0/h0**2)*((r0**2 + 4*h0**2)**1.5 - r0**3)
        Perimeter=np.pi*2*np.sqrt((a**2+b**2)/2)
        S_bit=S_bit+velocity*dt*Perimeter
        S_submerged = S_submerged_nose + 2 * np.pi * r0 * (x - nose_length)+S_box+S_bit
        print('wings smooth',S_submerged_nose, 2 * np.pi * r0 * (x - nose_length),S_box,S_bit)
        ReW = rho_W * velocity * x / mu_water
        print(ReW,'reynolds',x,velocity)
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2*r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact =2*velocity**2/2*K*4/3*np.pi*rho_water*a*(b)
        a=a+(box_height/2)/200
        #b=b+(wing_bit_width/2)/200

        drag=RF_flat#+RF_form
        dt=(wing_bit_length)/200/velocity
        time=time+dt
        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        print(x,buoyancy,drag,impact)
        print(velocity)
        water_drag_tab.append(drag)
        velocity_tab.append(velocity)
        drag_tab.append(drag+buoyancy+impact)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)

    print('wings')


    tail_position=body_length*.5+0.872
    x_wings=np.linspace(wing_location+wing_bit_length,tail_position,200)
    for x in x_wings:
        #print('volume', volume)
        #print(np.sqrt(1 - (a/b) ** 2), 'e')
        buoyancy=volume*rho_water*9.81
        r=r0/np.sqrt(h0)*np.sqrt(x)
        #volume_box=volume_box+(x-beginning_folded-box_length)*box_width*box_height
        #volume = np.pi * h0 * r0 ** 2 / 2 + np.pi * r0 ** 2 * (x - nose_length) + volume_box
        volume = np.pi * h0 * r0 ** 2 / 2 + volume_box#+ np.pi * r0 ** 2 * ( wing_location+wing_bit_length- nose_length) +  #+ 2*(airfoil_area)*(x-wing_location-wing_bit_length)*0.7
        #drag=surface*velocity*length*length/diameter
        S_submerged_nose = (np.pi/6)*(r0/h0**2)*((r0**2 + 4*h0**2)**1.5 - r0**3)
        Perimeter=np.pi*2*np.sqrt((a**2+b**2)/2)
        #S_box=S_box+velocity*dt*Perimeter
        S_submerged = S_submerged_nose + 2 * np.pi * r0 * (wing_location+wing_bit_length- nose_length)+S_box+ 2*Perimeter*(x-wing_location-wing_bit_length)+S_bit

        ReW = rho_W * velocity * x / mu_water
        print(ReW,'reynolds',x,velocity)
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2*r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact =0
        a=a+(box_height/2)/200
        #b=b+(wing_bit_width/2)/200

        drag=RF_flat#+RF_form
        dt=(wing_location-beginning_folded)/200/velocity
        time=time+dt
        velocity=velocity-(drag+buoyancy+impact-mass*9.81)/mass*dt
        print(x,buoyancy,drag,impact)
        print(velocity)
        water_drag_tab.append(drag)
        velocity_tab.append(velocity)
        drag_tab.append(drag+buoyancy+impact)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)
        print('S_bit', S_bit, 'the other one', 2 * Perimeter * (x - wing_location), 'S_box', S_box)
        xmax=x




    print('tail')
    print(x)
    tail_length=0.258
    b_tail=tail_length/2
    a_tail=b_tail*0.12
    K_axial=2.02

    end_tail=(.258-.129)+tail_position-0.05
    x_tail=np.linspace(tail_position,end_tail,200)
    # tail
    for x in x_tail:
        # print('volume', volume)
        # print(np.sqrt(1 - (a/b) ** 2), 'e')
        buoyancy = volume * rho_water * 9.81
        r = r0 / np.sqrt(h0) * np.sqrt(x)
        # volume_box=volume_box+(x-beginning_folded-box_length)*box_width*box_height
        # volume = np.pi * h0 * r0 ** 2 / 2 + np.pi * r0 ** 2 * (x - nose_length) + volume_box
        volume = np.pi * h0 * r0 ** 2 / 2 + volume_box+0.0011 #+ np.pi * r0 ** 2 * (wing_location + wing_bit_length - nose_length) + 2 * (airfoil_area) * (   x - wing_location - wing_bit_length)
        # drag=surface*velocity*length*length/diameter
        S_submerged_nose = (np.pi / 6) * (r0 / h0 ** 2) * ((r0 ** 2 + 4 * h0 ** 2) ** 1.5 - r0 ** 3)
        Perimeter = np.pi * 2 * np.sqrt((a ** 2 + b ** 2) / 2)
        # S_box=S_box+velocity*dt*Perimeter
        S_submerged = S_submerged_nose + 2 * np.pi * r0 * (
                wing_location + wing_bit_length - nose_length) + S_box + 2 * Perimeter * (
                              x - wing_location - wing_bit_length) * 0.7 + S_bit

        ReW = rho_W * np.abs(velocity) * x / mu_water
        print(ReW, 'reynolds', x, velocity)
        print(b_tail, 'btail', a_tail)
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * (velocity) ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2 * r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * (velocity) ** 2 * S_submerged * Cp
        if velocity > 0:
            impact = 3 * velocity ** 2 * K_axial * 2 / 3 * np.pi * rho_water * a_tail ** 2
        else:
            impact = 0
        b_tail = b_tail-(.258-0.179)/200
        a_tail = 0.12 * b_tail
        # b=b+(wing_bit_width/2)/200

        drag = RF_flat #+ RF_form
        dt = (end_tail - tail_position) / 200 / velocity
        time = time + dt
        if velocity > 0:
            velocity = velocity - (drag + buoyancy + impact - mass * 9.81) / mass * dt

        else:
            velocity = velocity - (-drag + buoyancy + impact - mass * 9.81) / mass * dt
            break
        print('S_bit', S_bit, 'S_total', S_submerged, 'S_box', S_box)
        print('x:', x, 'buoyancy', buoyancy, 'drag', drag, 'impact', impact, 'velocity:', velocity, 'Area', S_submerged)
        print(velocity)

        water_drag_tab.append(drag)
        velocity_tab.append(velocity)
        drag_tab.append(drag + buoyancy + impact)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)
        xmax = x
    print(b_tail,a_tail)
    print('Propeller')
    print(end_tail,end_tail+0.1)
    x_propeller=np.linspace(end_tail,end_tail+0.1,200)
    impact_water_prop_tab=[]
    r_water=0.01
    for x in x_propeller:
        # print('volume', volume)
        # print(np.sqrt(1 - (a/b) ** 2), 'e')
        buoyancy = volume * rho_water * 9.81
        r = r0 / np.sqrt(h0) * np.sqrt(x)
        # volume_box=volume_box+(x-beginning_folded-box_length)*box_width*box_height
        # volume = np.pi * h0 * r0 ** 2 / 2 + np.pi * r0 ** 2 * (x - nose_length) + volume_box
        volume = np.pi * h0 * r0 ** 2 / 2 + volume_box + 0.0011#+ np.pi * r0 ** 2 * (x - nose_length)+ 2 * (airfoil_area) * (
                    #x - wing_location) * 0.7 *0.7
        # drag=surface*velocity*length*length/diameter
        S_submerged_nose = (np.pi / 6) * (r0 / h0 ** 2) * ((r0 ** 2 + 4 * h0 ** 2) ** 1.5 - r0 ** 3)
        # Perimeter=np.pi*(3*a*b-np.sqrt((3*a+b)*(3*b+a)))
        # S_box=S_box+velocity*dt*Perimeter
        S_submerged = S_submerged_nose + 2 * np.pi * r0 * (body_length - nose_length) + S_box + 2 * Perimeter * (
                    x - wing_location) + S_bit + 2 * 0.056

        ReW = rho_W * velocity * x / mu_water
        print(ReW, 'reynolds', x, velocity)
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        RF_flat = 0.5 * rho_W * S_submerged * velocity ** 2 * CF_flat
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((nose_length / (2 * r)) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * velocity ** 2 * S_submerged * Cp
        impact_water_prop=8/3*rho_water*r_water**2*velocity**2
        impact = 3 * velocity ** 2 * K_axial * 2 / 3 * np.pi * rho_water * a_tail ** 2 +3*impact_water_prop
        b_tail = b_tail - (0.05) / 200
        a_tail = 0.12 * b_tail
        if x<end_tail+0.03:r_water=r_water+0.03/200

        impact_water_prop_tab.append(impact_water_prop)
        drag = RF_flat #+ RF_form
        dt = (wing_location + wing_bit_length + 1.25-end_tail) / 200 / velocity
        time = time + dt
        velocity = velocity - (drag + buoyancy + impact - mass * 9.81) / mass * dt
        print(x, 'buoyancy', buoyancy, 'drag', drag, impact)
        print(velocity)
        water_drag_tab.append(drag)
        velocity_tab.append(velocity)
        drag_tab.append(drag + buoyancy + impact)
        impact_tab.append(impact)
        buoyancy_tab.append(buoyancy)
        surface_tab.append(S_submerged)
        volume_tab.append(volume)

        xmax = x

    print(np.max(impact_water_prop_tab),x)
    print('end bit wings')


    #print(water_drag_tab)
    print('nose_length:', nose_length, 'beginning of the box', beginning_folded, 'box_to wings',
          beginning_folded + box_length, 'wing trailing edge smooth bit', wing_location, 'start wings',
          wing_location + wing_bit_length, 'start tail', body_length - tail_length, 'end tail', body_length,
          'end wings', wing_location + 1.25)
    plt.plot(np.linspace(0.001, xmax, len(drag_tab)), water_drag_tab)
    plt.xlabel('x location')
    plt.ylabel('Water drag')
    plt.show()
    plt.plot(np.linspace(0.001, xmax, len(velocity_tab)), velocity_tab)
    plt.xlabel('x location')
    plt.title('Velocity vs x')
    plt.ylabel('Velocity')
    plt.show()
    plt.plot(np.linspace(0.001, xmax, len(drag_tab)), drag_tab)
    plt.xlabel('x location')
    plt.ylabel('Forces up[N]')
    plt.show()
    plt.plot(np.linspace(0.001, xmax, len(velocity_tab)), impact_tab)
    plt.title('Impact force vs x')
    plt.xlabel('x location')
    plt.ylabel('Impact[N]')
    plt.show()
    plt.plot(np.linspace(0.001, xmax, len(drag_tab)), surface_tab)
    plt.xlabel('x location')
    plt.ylabel('Surface')
    plt.show()
    plt.plot(np.linspace(0.001, xmax, len(drag_tab)), volume_tab)
    plt.xlabel('x location')
    plt.ylabel('Volume')
    plt.show()
    plt.plot(np.linspace(0.001, xmax, len(drag_tab)), buoyancy_tab)
    plt.xlabel('x location')
    plt.ylabel('Buoyancy')
    plt.show()
    plt.plot(np.linspace(0.001, 0.1, len(impact_water_prop_tab)), impact_water_prop_tab)
    plt.xlabel('x location')
    plt.ylabel('Impact propeller')
    plt.show()
    print(np.max(impact_water_prop_tab))
    return velocity_tab


#def Force_versus_X(velocity_tab):


#F_lst_time, t_lst,pressure_lst_time = force_vs_time(0.3, 4, 0.912, 1030, 228.6, 0.3, 16)
#fig, (Ftime, Fangle) = plt.subplots(2)
#fig.suptitle('Diving load impacts')
#Ftime.plot(t_lst, F_lst_time, 'tab:orange')
#Ftime.set_xlabel('Time [s]')
#Ftime.set_ylabel('Force [N]')
#Ftime.set_title('Impact force vs. time')
#Fangle.plot(t_lst, pressure_lst_time, 'tab:red')
#Fangle.set_xlabel('Time [s]')
#Fangle.set_ylabel('Pressure [Pa]')
#Fangle.set_title('Pressure vs time')


#plt.show()

#def Beam_deflection(angle,Inertia):

#y=0.03x^{2}
#ef Velocity_change(v_initial,mass,body_length,nose_length,r0):
print(np.arctan(0.6)*180/np.pi)
print(np.sqrt(1-(0.12*2)**2),'e')
K=0.935
print(0.357*0.12,'whatever')
#Velocidef Velocity_change(v_initial: Any,
# mass: {__mul__},
                    #body_length: Any,
                    #nose_length: {__truediv__},
                    #r0: {__truediv__, __pow__},
                    #wing_location: {__sub__},
                    #wing_length: Any,
                    #root_chord: Any,
                    #airfoil_height: Any,
                    #K: Any,
                    #beginning_folded: Any)
Velocity_change(15,17,1.05,0.18,0.09,0.5986,1.5,0.357,0.12*0.357,K,0.5986-0.1,0.105,0.258)

#Nose end 0.18
#Body_length



#dimatere 10cm water prop
#body length 85 cm
#11cm
#2cm +tail length
#18cm fuselage
#20 cm
#0.0105 m2 wing

#0.5986wing location
#1.5m
#0.357 12% height