import aerosandbox as asb
import aerosandbox.numpy as np
# import Airfoil_Design as ad
import structural_sizing as ss
import matplotlib.pyplot as plt
from fuselage import make_fuselage

print('start')
wing = ss.wing

# Define your airfoil data

tail_chord_root = wing.tail_root_chord
tail_tip_chord = wing.tail_tip_chord
tail_span = wing.tail_span
actual_length_of_tail_surfaces = tail_span / 2 / np.cos(30 * np.pi / 180)
moment_arm = wing.l_opt
x_ac=0.671-wing.root_chord/4

print(tail_chord_root)
# Create the airfoil object
airfoil = asb.Airfoil(
    name="NACA2412",
    generate_polars=True
)

airfoil_tail = asb.Airfoil(
    name="NACA0012",
    generate_polars=True
)

# Define the geometry of the wing
wing = asb.Wing(
    xyz_le=[0, 0.5, 0.1],  # Position of the leading edge
    xsecs=[  # Define the cross sections of the wing
        asb.WingXSec(  # Root section
            xyz_le=[x_ac-0.5, 0, 0],  # Position of the leading edge
            chord=wing.root_chord,
            twist=0,  # In degrees
            airfoil=airfoil,
            num_chordwise=12,
            num_spanwise=12,
        ),
        asb.WingXSec(  # Tip section
            xyz_le=[x_ac + wing.le_tip - 0.5, wing.wing_span / 2, 0],  # Position of the leading edge
            chord=wing.tip_chord,
            twist=0,  # In degrees
            airfoil=airfoil,
            num_chordwise=12,
            num_spanwise=12,
        ),
    ],
    symmetric=True
)

tail = asb.Wing(
    xyz_le=[0, 0, 0],  # Position of the leading edge
    xsecs=[  # Define the cross sections of the wing
        asb.WingXSec(  # Root section
            xyz_le=[moment_arm - tail_chord_root, 0.1*np.sin(60 * np.pi / 180), -0.1- 0.1*np.sin(30*np.pi / 180)],  # Position of the leading edge
            chord=tail_chord_root,
            twist=0,  # In degrees
            airfoil=airfoil_tail,
            num_chordwise=12,
            num_spanwise=12,
        ),
        asb.WingXSec(  # Tip section
            xyz_le=[moment_arm-tail_tip_chord, tail_span / 2 + 0.1* np.sin(60 * np.pi / 180), -0.1- 0.1*np.sin(30*np.pi / 180) -actual_length_of_tail_surfaces * np.sin(30 * np.pi / 180)],  # Position of the leading edge
            chord=tail_tip_chord,
            twist=0,  # In degrees
            airfoil=airfoil_tail,
            num_chordwise=12,
            num_spanwise=12,
        ),
    ],
    symmetric=True
)

v_tail = asb.Wing(
    xyz_le=[0, 0, 0],  # Position of the leading edge
    xsecs=[  # Define the cross sections of the wing
        asb.WingXSec(  # Root section
            xyz_le=[moment_arm - tail_chord_root, 0,0],  # Position of the leading edge
            chord=tail_chord_root,
            twist=0,  # In degrees
            airfoil=airfoil_tail,
            num_chordwise=12,
            num_spanwise=12,
        ),
        asb.WingXSec(  # Tip section
            xyz_le=[moment_arm-tail_tip_chord, 0, actual_length_of_tail_surfaces],  # Position of the leading edge
            chord=tail_tip_chord,
            twist=0,  # In degrees
            airfoil=airfoil_tail,
            num_chordwise=12,
            num_spanwise=12,
        ),
    ],
    symmetric=False
)

center_fuse = make_fuselage(
    boom_length=moment_arm,
    nose_length=0.5,
    fuse_diameter=0.2,
    boom_diameter=0.2,
    fuse_resolution=10,
)
# Define the aircraft
aircraft = asb.Airplane(
    wings=[wing,tail,v_tail],
    fuselages=[center_fuse],
)

print("AHHHH")
print(asb.Airplane.aerodynamic_center(aircraft))

# aircraft.draw()
alpha_array = np.arange(-5, 15, 0.5)
beta_array = np.arange(-5, 12, 0.5)
r_array = np.arange(-2, 10, 0.5)
p_array = np.arange(-2, 10, 0.5)

cl_array = []
cd_array = []
cm_array = []
moments = []

#Sideslip derivatives
CYb_array=[]
Clb_array=[]
Cnb_array=[]
#Roll derivatives
CYp_array=[]
Clp_array=[]
Cnp_array=[]
#Yaw derivatives
CYr_array=[]
Clr_array=[]
Cnr_array=[]
#
# # Define the operating point
for alpha in alpha_array:
    op_point = asb.OperatingPoint(
        velocity=20,  # In m/s
        alpha=alpha,  # In degrees
        beta=0,  # In degrees
        p=0, q=0, r=0,  # In rad/s
    )

    # Perform the VLM analysis
    vlm = asb.VortexLatticeMethod(
        airplane=aircraft,
        op_point=op_point,
    )
    aero = asb.AeroBuildup(
        airplane=aircraft,
        op_point=op_point,
    )


    # aero_solve = aero.run_with_stability_derivatives(alpha=True,beta=True,p=True,q=True,r=True)
    aero_solve = aero.run_with_stability_derivatives(alpha=True)


    # Solve the VLM


    cl_array.append(aero_solve['CL'])
    cd_array.append(aero_solve['CD'])
    cm_array.append(aero_solve['Cm'])
    moments.append(aero_solve['m_b'])

    # # Sideslip derivatives
    # CYb_array.append(aero_solve['CYb'])
    # Clb_array.append(aero_solve['Clb'])
    # Cnb_array.append(aero_solve['Cnb'])
    # # Roll derivatives
    # CYp_array.append(aero_solve['CYp'])
    # Clp_array.append(aero_solve['Clp'])
    # Cnp_array.append(aero_solve['Cnp'])
    # # Yaw derivatives
    # CYr_array.append(aero_solve['CYr'])
    # Clr_array.append(aero_solve['Clr'])
    # Cnr_array.append(aero_solve['Cnr'])

    print(aero_solve)
    cla_array = []
    # Print the results
    print("Lift coefficient:", aero_solve['CL'])
    print("Drag coefficient:", aero_solve['CD'])
    print("Moment coefficient:", aero_solve['Cm'])
    cla_array.append(aero_solve['CLa'])

max_moment = np.max(abs(np.array(moments)))





plt.subplot(3, 2, 1)
plt.plot(alpha_array, cl_array)
plt.grid()
plt.title('Lift Curve')
plt.subplot(3, 2, 2)
plt.plot(alpha_array, cd_array)
plt.grid()
plt.title('Drag Curve')
plt.subplot(3, 2, 3)
plt.plot(cd_array, cl_array)
plt.grid()
plt.title('Drag Polar')
plt.subplot(3, 2, 4)
plt.plot(alpha_array, cm_array)
plt.grid()
plt.title('Moment Curve')
plt.subplot(3, 2, 5)
plt.plot(alpha_array, np.array(cl_array) / np.array(cd_array))
plt.grid()
plt.title('Lift over Drag Curve')

plt.tight_layout()
plt.suptitle('Aerodynamic Properties Full Puffin')
plt.show()
#
# plt.plot(alpha_array, cm_array)
# plt.title('Moment Curve')
# plt.xlabel("Alpha")
# plt.ylabel("Cm")
# plt.show()

'''



# Define the operating point
for beta in beta_array:
    op_point = asb.OperatingPoint(
        velocity=20,  # In m/s
        alpha=0,  # In degrees
        beta=beta,  # In degrees
        p=0, q=0, r=0,  # In rad/s
    )

    # Perform the VLM analysis
    vlm = asb.VortexLatticeMethod(
        airplane=aircraft,
        op_point=op_point,
    )
    aero = asb.AeroBuildup(
        airplane=aircraft,
        op_point=op_point,
    )


    # aero_solve = aero.run_with_stability_derivatives(alpha=True,beta=True,p=True,q=True,r=True)
    aero_solve = aero.run_with_stability_derivatives(beta=True)


    # Solve the VLM

    #
    cl_array.append(aero_solve['CL'])
    cd_array.append(aero_solve['CD'])
    cm_array.append(aero_solve['Cm'])
    #moments.append(aero_solve['m_b'])

    # Sideslip derivatives
    #CYb_array.append(aero_solve['CYb'])
    #Clb_array.append(aero_solve['Clb'])
    #Cnb_array.append(aero_solve['Cnb'])
    # # Roll derivatives
    # CYp_array.append(aero_solve['CYp'])
    # Clp_array.append(aero_solve['Clp'])
    # Cnp_array.append(aero_solve['Cnp'])
    # # Yaw derivatives
    # CYr_array.append(aero_solve['CYr'])
    # Clr_array.append(aero_solve['Clr'])
    # Cnr_array.append(aero_solve['Cnr'])


# Define the operating point
for r in r_array:
    op_point = asb.OperatingPoint(
        velocity=20,  # In m/s
        alpha=0,  # In degrees
        beta=0,  # In degrees
        p=0, q=0, r=r,  # In rad/s
    )

    # Perform the VLM analysis
    vlm = asb.VortexLatticeMethod(
        airplane=aircraft,
        op_point=op_point,
    )
    aero = asb.AeroBuildup(
        airplane=aircraft,
        op_point=op_point,
    )


    # aero_solve = aero.run_with_stability_derivatives(alpha=True,beta=True,p=True,q=True,r=True)
    aero_solve = aero.run_with_stability_derivatives(r=True)


    # Solve the VLM


    # cl_array.append(aero_solve['CL'])
    # cd_array.append(aero_solve['CD'])
    # cm_array.append(aero_solve['Cm'])
    # moments.append(aero_solve['m_b'])
    #
    # # Sideslip derivatives
    # CYb_array.append(aero_solve['CYb'])
    # Clb_array.append(aero_solve['Clb'])
    # Cnb_array.append(aero_solve['Cnb'])
    # # Roll derivatives
    # CYp_array.append(aero_solve['CYp'])
    # Clp_array.append(aero_solve['Clp'])
    # Cnp_array.append(aero_solve['Cnp'])
    # Yaw derivatives
    CYr_array.append(aero_solve['CYr'])
    Clr_array.append(aero_solve['Clr'])
    Cnr_array.append(aero_solve['Cnr'])


# Define the operating point
for p in p_array:
    op_point = asb.OperatingPoint(
        velocity=20,  # In m/s
        alpha=0,  # In degrees
        beta=0,  # In degrees
        p=p, q=0, r=0,  # In rad/s
    )

    # Perform the VLM analysis
    vlm = asb.VortexLatticeMethod(
        airplane=aircraft,
        op_point=op_point,
    )
    aero = asb.AeroBuildup(
        airplane=aircraft,
        op_point=op_point,
    )


    # aero_solve = aero.run_with_stability_derivatives(alpha=True,beta=True,p=True,q=True,r=True)
    aero_solve = aero.run_with_stability_derivatives(p=True)


    # Solve the VLM


    # cl_array.append(aero_solve['CL'])
    # cd_array.append(aero_solve['CD'])
    # cm_array.append(aero_solve['Cm'])
    # moments.append(aero_solve['m_b'])
    #
    # # Sideslip derivatives
    # CYb_array.append(aero_solve['CYb'])
    # Clb_array.append(aero_solve['Clb'])
    # Cnb_array.append(aero_solve['Cnb'])
    # # Roll derivatives
    CYp_array.append(aero_solve['CYp'])
    Clp_array.append(aero_solve['Clp'])
    Cnp_array.append(aero_solve['Cnp'])
    # Yaw derivatives
    # CYr_array.append(aero_solve['CYr'])
    # Clr_array.append(aero_solve['Clr'])
    # Cnr_array.append(aero_solve['Cnr'])

plt.plot(beta_array, CYb_array)
plt.xlabel("Beta [deg]")
plt.ylabel("CYb")
plt.show()

Clb_array=list(map(lambda Clb: Clb - 0.055, Clb_array))


plt.plot(beta_array, Clb_array)
plt.xlabel("Beta [deg]")
plt.ylabel("Clb")
plt.show()

plt.plot(beta_array, Cnb_array)
plt.xlabel("Beta [deg]")
plt.ylabel("Cnb")
plt.show()



CYp_array=list(map(lambda CYp: CYp - 0.055, CYp_array))

plt.plot(p_array, CYp_array)
plt.xlabel("Roll rate [rad/s]")
plt.ylabel("CYp")
plt.show()


Clp_array=list(map(lambda Clp: Clp - 0.5, Clp_array))

plt.plot(p_array, Clp_array)
plt.xlabel("Roll rate [rad/s]")
plt.ylabel("Clp")
plt.show()

Cnp_array=list(map(lambda Cnp: Cnp - 0.3, Cnp_array))


plt.plot(p_array, Cnp_array)
plt.xlabel("Roll rate [rad/s]")
plt.ylabel("Cnp")
plt.show()

CYr_array=list(map(lambda CYr: CYr*0.003, CYr_array))


plt.plot(r_array, CYr_array)
plt.xlabel("Yaw rate [rad/s]")
plt.ylabel("CYr")
plt.show()

plt.plot(r_array, Clr_array)
plt.xlabel("Yaw rate [rad/s]")
plt.ylabel("Clr")
plt.show()

plt.plot(r_array, Cnr_array)
plt.xlabel("Yaw rate [rad/s]")
plt.ylabel("Cnr")
plt.show()



#
# def asymmetric(alpha):
#
#     muc, mub, CL, CD, CX0, CZ0 =
#     V0 = 20
#     b=wing.wing_area
#     Db = (b / V0)
#
#     X = np.matrix([[Db * (CYbdot - 2 * mub), 0, 0, 0],
#                     [0, -0.5 * Db, 0, 0],
#                     [0, 0, -4 * mub * KX2 * Db, 4 * mub * KXZ * Db],
#                     [Cnbdot * Db, 0, 4 * mub * KXZ * Db, -4 * mub * KZ2 * Db]])
#
#     Y= np.matrix([[CYb, CL, CYp, (CYr - 4 * mub)],
#                     [0, 0, 1, 0],
#                     [Clb, 0, Clp, Clr],
#                     [Cnb, 0, Cnp, Cnr]])
#
#     Z = np.matrix([[CYda, CYdr],
#                     [0, 0],
#                     [Clda, Cldr],
#                     [Cnda, Cndr]])
#
#     dim_x=4
#     dim_y=2
#
#     A = -np.linalg.pinv(X) @Y
#     B = -np.linalg.pinv(X) @Z
#     C = np.eye(dim_x)
#     D = np.zeros((dim_x,dim_y))
#
#
#     sys=ctrl.ss(A,B,C,D)
#
#     return sys
'''