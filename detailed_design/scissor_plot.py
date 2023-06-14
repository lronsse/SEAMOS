import matplotlib.pyplot as plt
import numpy as np
from structural_sizing import Wing


## Constants
Vh_V = 1
x_ac = 0.26 #
Mach = 0.1
beta = np.sqrt(1-Mach**2)
MAC = Wing.mean_aerodynamic_chord

eta_tail = 0.95
eta_wing = 0.95

x_cg_fw = 0.45#0.5585
x_cg_aft = 0.55#0.5694

A = 12
Ah = 4 #(2/3)*A # Not sure
half_sweep_wing = 2 #deg
half_sweep_tail = 4.76 #deg
b_f = 0.3 #fuselage diameter
b = 3.87 #wingspan
S = 1.25 #surface area
Sh = 11.73 #tail area
S_net = S - 7.36 # Kinda sketch (I did S_covered = root chord * b_f :-/)
taper = 0.4 # wing

print("OMG")
print((2*S)/((1+taper)*b))
l_h = 0.45 # Not accurate *******
r = 2*l_h/b
m_b_2 = 4.483
m = m_b_2 * 2/b
B_p = 3 # Number of propeller blades per propeller

Cm_ac = -0.3 # Not accurate ******
# CL_h = -0.35*Ah**(1/3) # Fixed tail
CL_h = -1 # Full moving tail

def Sh_S_stability(x_cg, CL_a_h, CL_a_Aminh, downwash, l_h=l_h, MAC=MAC, Vh_V=Vh_V, x_ac=x_ac, safety=0.05):
    """x_cg is the free variable"""
    return (1 / (CL_a_h / CL_a_Aminh * (1 - downwash) * l_h/MAC * Vh_V**2) * x_cg -
            (x_ac - safety) / (CL_a_h / CL_a_Aminh * (1 - downwash) * l_h/MAC * Vh_V**2))


def Sh_S_controllability(x_cg, CL_h, CL_a_Aminh, l_h=l_h, MAC=MAC, Vh_V=Vh_V, x_ac=x_ac, Cm_ac=Cm_ac):
    return (1 / (CL_h / CL_a_Aminh * l_h/MAC * Vh_V**2) * x_cg +
            (Cm_ac/CL_a_Aminh - x_ac) / (CL_h / CL_a_Aminh * l_h/MAC * Vh_V**2))



## Parameter functions
def CL_a_DATCOM(A, eta, half_sweep, beta):
    return 2*np.pi*A / (2 + np.sqrt(4 + (A*beta/eta)**2 * (1 + (np.tan(np.deg2rad(half_sweep))/beta)**2 )))

def CL_a_Aminh(CL_a_w, bf, b, S_net, S):
    return CL_a_w * (1 + 2.15*b_f/b) * S_net/S + 0.5*np.pi*bf**2/S

def downwash(CL_a_w, A=A, taper=taper, r=r, m=m, B_p=B_p):
    return (1.75 * CL_a_w / (np.pi * A * (taper*r)**0.25 * (1+np.abs(m)))) * (1 - 0.012 * B_p)

def tail_area_S(CLah,Cla_dash_h,deda,lh,c,VhV,xac,SM,S):
    ShS=(1/((CLah/Cla_dash_h)*(1-deda)*(lh/c)*VhV**2))-((xac-SM)/((CLah/Cla_dash_h)*(1-deda)*(lh/c)*VhV**2))
    Sh=ShS*S
    return Sh

# CLah=CL_a_DATCOM(A,eta_wing,half_sweep_wing,beta)
# Cla_dash_h=CL_a_Aminh_stab = CL_a_Aminh(CLah, b_f, b, S_net, S)
# deda = downwash(Cla_dash_h)
SM=0.05
lh=1.5
xac=0.55
print("ahhhh")
# print(tail_area_S(CLah,Cla_dash_h,deda,lh,MAC,Vh_V,xac,SM,S))



CL_a_w = CL_a_DATCOM(A, eta_wing, half_sweep_wing, beta)
CL_a_w_low = CL_a_DATCOM(A, eta_wing, half_sweep_wing, np.sqrt(1-0.171**2)) # Entered take-off speed here (take-off is limiting for controllability, was too lazy to make new Mach variable)
CL_a_h = CL_a_DATCOM(Ah, eta_tail, half_sweep_tail, beta)
dw = downwash(CL_a_w)
CL_a_Aminh_stab = CL_a_Aminh(CL_a_w, b_f, b, S_net, S)
CL_a_Aminh_cont = CL_a_Aminh(CL_a_w_low, b_f, b, S_net, S)


print(CL_a_Aminh_cont)
print(CL_a_Aminh_stab)
print(CL_a_h)


def tail_area_S(CLah,Cla_dash_h,deda,lh,c,VhV,xac,SM,S):
    ShS=(1/((CLah/Cla_dash_h)*(1-deda)*(lh/c)*VhV**2))-((xac-SM)/((CLah/Cla_dash_h)*(1-deda)*(lh/c)*VhV**2))
    Sh=ShS*S
    return Sh

CLadashh=5.65
CLah=4.722
print("AHHHUSETHIS")
area_tail=tail_area_S(CLah,CLadashh,dw,lh,MAC,Vh_V,xac,SM,S)
print(area_tail)
span_tail=np.sqrt(8*(area_tail))
print(span_tail)
# print(Sh_S_stability(0.5,4.722,5.65,dw(4.722),1.5,MAC,VhV,0.55,0.05))


stability_params = {'CL_a_h': CL_a_h,
                    "CL_a_Aminh": CL_a_Aminh_stab,
                    'downwash': dw}

controllability_params = {'CL_h': CL_h,
                          "CL_a_Aminh": CL_a_Aminh_cont}

x_cg_lst = np.arange(0, 1.01, 0.01)
stab_line = [Sh_S_stability(x, **stability_params) for x in x_cg_lst]
stabb_line = [Sh_S_stability(x, **stability_params, safety=0) for x in x_cg_lst]
cont_line = [Sh_S_controllability(x, **controllability_params) for x in x_cg_lst]

plt.plot(x_cg_lst, stab_line, label='stability')
plt.plot(x_cg_lst, stabb_line, label='neutral stability')
plt.fill_between(x_cg_lst, stabb_line, facecolor="lightgray", hatch="X", edgecolor="darkgray", linewidth=0.0)

plt.plot(x_cg_lst, cont_line, label='controllability')
plt.fill_between(x_cg_lst, cont_line, facecolor="lightgray", hatch="X", edgecolor="darkgray", linewidth=0.0)

plt.axvline(x_cg_fw, linestyle='dashed', color='black')
plt.axvline(x_cg_aft, linestyle='dashed', color='black')

plt.axhline(Sh/S, color='black')

plt.ylim(0, 0.5)
plt.xlim(0, 1)
plt.grid()
plt.xlabel(r'$X_{cg}$ [MAC]', fontsize=15)
plt.ylabel(r'$S_h / S$', fontsize=15)
plt.gca().set_aspect(2)
plt.legend(fontsize=15)
plt.show()
# print(downwash(CL_a_w))
