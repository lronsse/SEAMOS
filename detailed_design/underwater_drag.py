import numpy as np
import math

##Drag estimation method from Submarine Hydrodynamics textbook page 107

#Constants
V_water=1 #velocity [m/s]
rho_W = 1023 #denisty of water [kg/m^3]
L_hull = 1.5 #Length of hull [m]
D_hull = 0.3 #Diameter of hull [m]
eta_hull = 5 #3-6 depending on type of hull
A_front = (4*np.pi*(D_hull/2)**2)/2 #Front area of hull (assuming semi-sphere) [m^2]
Cd_water=0.02 #from literature
# A_plan =1.25+(0.146*3) #Wing area + tail area
AR=12
S=0.77
taper_ratio=0.4
b=np.sqrt(AR*S)
c_root=(2*S)/((1+taper_ratio)*b)
overlapped_area=0.5*(2*c_root)*(b/2)-(0.5*b/2*(2*c_root-D_hull))
new_area=S-overlapped_area
A_plan =new_area+(0.146*3) #Swept overlapped Wing area + tail area
# mu_water=1.3 * 10 ** ( -6 )
mu_water = 0.00126
ReW = rho_W*V_water*L_hull/mu_water  #Reynolds number
S=2*np.pi*(D_hull/2)*(L_hull-(D_hull/2))+2*np.pi*((D_hull/2)**2)+A_front


#Frictional resistance from hull
def RF_flat(rho_W,L_hull,D_hull,V_water,ReW):
    CF_flat=0.0667/((math.log10(ReW)-2)**2)
    S_hull=2.25*L_hull*D_hull
    RF_flat=0.5*rho_W*S_hull*V_water**2*CF_flat
    return RF_flat


#Hull form drag
def Hull_form_drag(ReW,eta_hull,L_hull,D_hull,rho_W,V_water,A_front):
    CF_form=0.075/((math.log10(ReW)-2)**2)
    Kp=eta_hull*((L_hull/D_hull)**(-1.7))
    Cp=Kp*CF_form
    RF_form=0.5*rho_W*V_water**2*A_front*Cp
    return RF_form

#Wing+tail skin friction

def control_surface_friction(Cd_water,A_plan):
    R_surface=0.5*rho_W*V_water**2*A_plan*Cd_water
    return R_surface

#Total drag
def total_drag():
    RF_flat_drag=RF_flat(rho_W, L_hull, D_hull, V_water, ReW)
    Hull_form_drag_drag= Hull_form_drag(ReW, eta_hull, L_hull, D_hull, rho_W, V_water,A_front)
    wingtail_friction_drag=control_surface_friction(Cd_water, A_plan)
    drag = RF_flat_drag + Hull_form_drag_drag+ wingtail_friction_drag
    Cd = drag / (0.5 * rho_W * V_water ** 2 * S)
    return drag,Cd,RF_flat_drag,Hull_form_drag_drag,wingtail_friction_drag


print(total_drag())

#
# print(0.71*16**0.48)
# print((0.71*16**0.48)/8.5)
# print(ReW/1e7)
# print(drag,Cd)