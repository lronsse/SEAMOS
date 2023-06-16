import matplotlib.pyplot as plt
import numpy as np
import structural_sizing as ss

wing = ss.wing

## Constants
Vh_V = 0.65
x_ac = 0.5
Mach = 0.1
beta = np.sqrt(1-Mach**2)
beta_low=np.sqrt(1-0.0437318**2)
MAC = wing.mean_aerodynamic_chord

eta_tail = 0.95
eta_wing = 0.95

x_cg_fw = 0.45 #0.5585
x_cg_aft = 0.55#0.5694
dxcg=x_cg_aft-x_cg_fw


tail_taper_ratio=0.4
A = wing.aspect_ratio
Ah = 4 #tail AR
AR_tail=4
half_sweep_wing = np.arctan(np.tan(wing.sweep_quarter_chord) - (4 / A) * (50 - 25) / 100 * (
                    1 - wing.taper_ratio) / (1 + wing.taper_ratio))
half_sweep_tail = np.arctan(np.tan(0) - (4 / AR_tail) * (50 - 100) / 100 * (
                    1 - tail_taper_ratio) / (1 + tail_taper_ratio))
b_f = 0.2 #Fuselage diameter [m]
b =wing.wing_span #wingspan [m]
S =  0.75 #wing surface area [m^2]
Sh = wing.tail_area #wing surface area [m^2]
S_net = S - (wing.root_chord*b_f) # Kinda sketch (I did S_covered = root chord * b_f :-/)
taper = wing.taper_ratio # wing

l_h = wing.tail_arm # Not accurate *******
r = 2*l_h/b
m_b_2 = 0.2
m = m_b_2 * 2/b
B_p = 2 # Number of propeller blades per propeller
l_f=0.84
delta=wing.sweep_quarter_chord
CL_zero=0.2
cm0_airfoil=-0.068
CL_a_dash_h_lowv=(2*np.pi*A / (2 + np.sqrt(4 + (A*beta_low/eta_wing)**2 * (1 + (np.tan((half_sweep_wing))/beta_low)**2 ))))* (1 + 2.15*b_f/b) * S_net/S + 0.5*np.pi*b_f**2/S

cmac_wing=cm0_airfoil*(wing.aspect_ratio*np.cos(delta)**2/(wing.aspect_ratio+2*np.cos(delta)))
cmac_fuselage=-1.8*(1-(2.5*b_f/l_f))*(np.pi*b_f*b_f*l_f)/(4*S*wing.mean_aerodynamic_chord)*(CL_zero/CL_a_dash_h_lowv)


Cm_ac = cmac_wing+cmac_fuselage
CL_h = -0.8 #-0.35*Ah**(1/3) # Fixed tail
print(CL_h)
print("AHHHHH")



def Sh_S_stability(x_cg, CL_a_h, CL_a_Aminh, downwash, l_h=l_h, MAC=MAC, Vh_V=Vh_V, x_ac=x_ac, safety=0.05):
    """x_cg is the free variable"""
    return (1 / (CL_a_h / CL_a_Aminh * (1 - downwash) * l_h/MAC * Vh_V**2) * x_cg -
            (x_ac - safety) / (CL_a_h / CL_a_Aminh * (1 - downwash) * l_h/MAC * Vh_V**2))


def Sh_S_controllability(x_cg, CL_h, CL_a_Aminh, l_h=l_h, MAC=MAC, Vh_V=Vh_V, x_ac=x_ac, Cm_ac=Cm_ac):
    return (1 / (CL_h / CL_a_Aminh * l_h/MAC * Vh_V**2) * x_cg +
            (Cm_ac/CL_a_Aminh - x_ac) / (CL_h / CL_a_Aminh * l_h/MAC * Vh_V**2))



## Parameter functions
def CL_a_DATCOM(A, eta, half_sweep, beta):
    return 2*np.pi*A / (2 + np.sqrt(4 + (A*beta/eta)**2 * (1 + (np.tan((half_sweep))/beta)**2 )))

def CL_a_Aminh(CL_a_w, bf, b, S_net, S):
    return CL_a_w * (1 + 2.15*b_f/b) * S_net/S + 0.5*np.pi*bf**2/S

def downwash(CL_a_w, A=A, taper=taper, r=r, m=m, B_p=B_p):
    return (1.75 * CL_a_w / (np.pi * A * (taper*r)**0.25 * (1+np.abs(m)))) * (1 - 0.012 * B_p)


CL_a_w = CL_a_DATCOM(A, eta_wing, half_sweep_wing, beta)
CL_a_w_low = CL_a_DATCOM(A, eta_wing, half_sweep_wing, np.sqrt(1-0.0437318**2)) # Entered take-off speed here (take-off is limiting for controllability, was too lazy to make new Mach variable)
CL_a_h = CL_a_DATCOM(Ah, eta_tail, half_sweep_tail, beta)
dw = downwash(CL_a_w)
CL_a_Aminh_stab = CL_a_Aminh(CL_a_w, b_f, b, S_net, S)
CL_a_Aminh_cont = CL_a_Aminh(CL_a_w_low, b_f, b, S_net, S)

#
# def horizontal_tailplan():
#     ShS_hori=(dxcg)/

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
print(downwash(CL_a_w))