import aerosandbox as asb
import aerosandbox.numpy as np
import Airfoil_Design as ad
import structural_sizing as ss
import matplotlib.pyplot as plt
from fuselage import make_fuselage

print('start')
wing = ss.wing

# Define your airfoil data

alpha_data = ad.data.index
cl_data = ad.data['CL']
cd_data = ad.data['CD']
cm_data = ad.data['Cm']

tail_chord_root = wing.tail_root_chord
tail_tip_chord = wing.tail_tip_chord
tail_span = wing.tail_span
actual_length_of_tail_surfaces = tail_span / 2 / np.cos(30 * np.pi / 180)
moment_arm = wing.l_opt
x_ac=0.6

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
    xyz_le=[0, 0, 0],  # Position of the leading edge
    xsecs=[  # Define the cross sections of the wing
        asb.WingXSec(  # Root section
            xyz_le=[0-x_ac, 0, 0],  # Position of the leading edge
            chord=wing.root_chord,
            twist=0,  # In degrees
            airfoil=airfoil,
            num_chordwise=12,
            num_spanwise=12,
        ),
        asb.WingXSec(  # Tip section
            xyz_le=[moment_arm-wing.le_tip, wing.wing_span / 2, 0],  # Position of the leading edge
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
    wings=[wing, v_tail, tail],
    fuselages=[center_fuse],
)

aircraft.draw()
alpha_array = np.arange(-5, 15, 0.5)
cl_array = []
cd_array = []
cm_array = []
moments = []
# Define the operating point
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


    aero_solve = aero.run_with_stability_derivatives()



    # Solve the VLM



    cl_array.append(aero_solve['CL'])
    cd_array.append(aero_solve['CD'])
    cm_array.append(aero_solve['Cm'])
    moments.append(aero_solve['m_b'])

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
