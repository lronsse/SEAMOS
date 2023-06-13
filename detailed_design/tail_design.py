import numpy as np
import math

#constants
AR_wing=12 #Wing aspect ratio
Vv=0.03 #Vertical tail volume coefficient (from literature)
b_wing=2.5 #Wingspan [m]
S_wing=1.25 #Wing surface area [m^2]
tail_arm=0.75 #Tail arm to cg [m]
c_wing=0.52 #MAC of wing [m]
Vh=0.35 #Horizontal tail volume coefficient (from literature)
tail_taper_ratio=0.5 #what do you think this is

def tail_sizing(AR_wing,b_wing,S_wing,tail_arm,c_wing,tail_taper_ratio):
    Vv = 0.03  # Vertical tail volume coefficient (from literature)
    Vh = 0.35  # Horizontal tail volume coefficient (from literature)
    AR_tail=(2/3)*AR_wing
    Sh=(0.85*Vh*c_wing*S_wing)/tail_arm
    Sv=(0.85*Vv*b_wing*S_wing)/tail_arm
    S_projected_v=0.33*Sv
    S_projected_h=Sv-S_projected_v
    tail_anhedral=np.degrees(np.arctan(np.sqrt(S_projected_v/Sh)))
    tail_area=0.5*(Sh/(np.cos(np.radians(tail_anhedral)))**2) #
    tail_span=np.sqrt(AR_tail*tail_area) #
    tail_root_chord=(2*tail_area)/((1+tail_taper_ratio)*tail_span) #
    tail_tip_chord=tail_root_chord*tail_taper_ratio #
    tail_mac=(2/3)*(tail_root_chord)*(1+tail_taper_ratio+tail_taper_ratio**2)/(1+tail_taper_ratio) #
    tail_qc_sweep=np.degrees(np.arctan((((np.tan(0))-(4/AR_tail)*((-75/100)*((1-tail_taper_ratio)/(1+tail_taper_ratio))))))) #
    return tail_area,tail_span,tail_root_chord,tail_tip_chord,tail_mac,tail_qc_sweep,tail_anhedral,tail_qc_sweep

