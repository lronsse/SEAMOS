import numpy as np
import math
# from structural_sizing import Wing
import matplotlib.pyplot as plt

# wing=ss.Wing
##Drag estimation method from Submarine Hydrodynamics textbook page 107

#Constants
V_water=2 #velocity [m/s]
V_trans=2#velocity [m/s]
rho_W = 1023 #denisty of water [kg/m^3]
L_hull = 1.764#Length of hull [m]
slender=9.8 #L/D
D_hull = L_hull/slender #Diameter of hull [m]
eta_hull = 5 #3-6 depending on type of hull
A_front = (4*np.pi*(D_hull/2)**2)/2 #Front area of hull (assuming semi-sphere) [m^2]
Cd_water=0.02 #from literature
# A_plan =1.25+(0.146*3) #Wing area + tail area #diving in wings unfolded
AR=12
S=0.75
taper_ratio=0.4
b=np.sqrt(AR*S)
c_root=(2*S)/((1+taper_ratio)*b)
hinge_diameter=0
overlapped_area=0.5*(2*c_root)*(b/2)-(0.5*b/2*(2*c_root-hinge_diameter))
new_area=(3*0.12*0.265)-overlapped_area+(3*(0.09*0.2*0.29))
# A_plan =new_area +(0.05624*3) #Swept overlapped Wing area + tail area
A_plan= 2*(0.3571428571428572*0.12*0.26530612244897966)+(3*(0.09*0.2*0.29))
# A_trans= 2*np.pi*D_hull/2*L_hull+4*np.pi*(D_hull/2)**2
A_trans=(0.18*1.764)+(2*(0.258*0.29))+(0.12*0.26530612244897966*0.26530612244897966)
# mu_water=1.3 * 10 ** ( -6 )
mu_water = 0.00126
ReW = rho_W*V_water*L_hull/mu_water  #Reynolds number
ReWt = rho_W*V_trans*L_hull/mu_water  #Reynolds number

S=2*np.pi*(D_hull/2)*(L_hull-(D_hull/2))+2*np.pi*((D_hull/2)**2)+A_front
hull_transverse_area = (0.18 * 1.764) #Hull side
control_surface_side = (2.2 * (0.258 * 0.29)) + (0.12 * 0.26530612244897966 * 0.26530612244897966) #ccontrol surface side

#Frictional resistance from hull
def Hull_friction_flat(rho_W,L_hull,D_hull,V_water,ReW):
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

def RF_tv(rho_W,V_trans,ReWt,A_trans):
    CF_flat=0.0667/((math.log10(ReWt)-2)**2)
    RF_transverse= 0.5*rho_W*A_trans*V_trans**2*CF_flat
    return RF_transverse

def Hull_form_drag_tv(ReWt,eta_hull,L_hull,D_hull,rho_W,V_trans,A_trans):
    CF_form=0.075/((math.log10(ReWt)-2)**2)
    Kp=eta_hull*((L_hull/D_hull)**(-1.7))
    Cp=Kp*CF_form
    RF_form=0.5*rho_W*V_trans**2*A_trans*Cp
    return RF_form

def control_surface_friction_tv(Cd_water,control_surface_side):
    R_surface=0.5*rho_W*V_trans**2*control_surface_side*Cd_water*0.02
    return R_surface

#Total drag
def total_drag(V_water):
    RF_flat_drag=Hull_friction_flat(rho_W,L_hull,D_hull,V_water,ReW)
    Hull_form_drag_drag= Hull_form_drag(ReW, eta_hull, L_hull, D_hull, rho_W, V_water,A_front)
    wingtail_friction_drag=control_surface_friction(Cd_water, A_plan)
    drag = RF_flat_drag + Hull_form_drag_drag+ wingtail_friction_drag
    Cd = drag / (0.5 * rho_W * V_water ** 2 * S)
    return drag,Cd,RF_flat_drag,Hull_form_drag_drag,wingtail_friction_drag

def total_drag_tv(V_trans):
    RF_flat_drag_tv=RF_tv(rho_W,V_trans,ReWt,A_trans)
    Hull_form_drag_drag_tv= Hull_form_drag_tv(ReWt,eta_hull,L_hull,D_hull,rho_W,V_trans,A_trans)
    wingtail_friction_drag_tv=control_surface_friction_tv(Cd_water,control_surface_side)
    drag_tv = RF_flat_drag_tv + Hull_form_drag_drag_tv+ wingtail_friction_drag_tv
    Cd_tv = drag_tv / (0.5 * rho_W * V_water ** 2 * S)
    return drag_tv,Cd_tv,RF_flat_drag_tv,Hull_form_drag_drag_tv,wingtail_friction_drag_tv

# def resultant_drag():
#     drag, Cd, RF_flat_drag, Hull_form_drag_drag, wingtail_friction_drag = total_drag()
#     drag_tv, Cd_tv, RF_flat_drag_tv, Hull_form_drag_drag_tv, wingtail_friction_drag_tv = total_drag_tv()
#     res_drag = np.sqrt(drag**2 + drag_tv**2)
#     return res_drag, f"axial: {drag}", f"transverse: {drag_tv}"

# # Define a range of velocities
velocities = np.linspace(0.1, 10, 100)
# #
# Lists to store the drag values
total_drag_values = []
total_drag_tv_values = []
resultant_drag_values = []
Hull_drag=[]
Pressure_drag=[]
Wingtail_drag=[]

# Calculate drag for each velocity
for V_waters in velocities:
    drag, Cd, RF_flat_drag, Hull_form_drag_drag, wingtail_friction_drag = total_drag(V_waters)
    drag_tv, Cd_tv, RF_flat_drag_tv, Hull_form_drag_drag_tv, wingtail_friction_drag_tv = total_drag_tv(V_waters)
    res_drag =np.sqrt(drag**2 + drag_tv**2)
    total_drag_values.append(drag)
    total_drag_tv_values.append(drag_tv)
    resultant_drag_values.append(res_drag)
    Hull_drag.append(RF_flat_drag)
    Pressure_drag.append(Hull_form_drag_drag)
    Wingtail_drag.append(wingtail_friction_drag)

# Plot the results
plt.plot(velocities, total_drag_values, label='Longitudinal Drag')
plt.plot(velocities, total_drag_tv_values, label='Transverse Drag')
plt.plot(velocities, resultant_drag_values, label='Resultant Drag')

plt.xlabel('Velocity [m/s]')
plt.ylabel('Underwater Drag [N]')
plt.title('Underwater Drag vs. Velocity')
plt.legend()
plt.grid(True)
plt.show()


plt.plot(velocities, total_drag_values, label='Longitudinal Drag')
plt.plot(velocities, Hull_drag, label='Hull Friction Drag')
plt.plot(velocities, Pressure_drag, label='Pressure Drag')
plt.plot(velocities, Wingtail_drag, label='Control Surface Friction Drag')


plt.xlabel('Velocity [m/s]')
plt.ylabel('Underwater Drag [N]')
plt.title('Longitudinal  Drag vs. Velocity')
plt.legend()
plt.grid(True)
# plt.show()




print("drag")
print(total_drag(V_water))
print(total_drag_tv(V_trans))
# print(resultant_drag())