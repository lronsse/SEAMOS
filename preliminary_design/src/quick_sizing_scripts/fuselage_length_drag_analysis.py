import numpy as np
import math
import matplotlib.pyplot as plt

# Constants
V_water = 1  # velocity [m/s]
rho_W = 1023  # density of water [kg/m^3]
slender_ratios = np.arange(5, 12, 1)  # list of slenderness ratios to test
lengths = np.arange(0.5, 3, 0.5)  # list of lengths to test
AR = 12
taper_ratio = 0.4

# Create empty 2D array to store drag values
drag_values = np.zeros((len(lengths), len(slender_ratios)))

# Drag estimation method
def calculate_drag(length, slender):
    D_hull = length / slender  # Diameter of hull [m]
    eta_hull = 5  # 3-6 depending on the type of hull
    A_front = (4 * np.pi * (D_hull / 2) ** 2) / 2  # Front area of hull (assuming semi-sphere) [m^2]
    Cd_water = 0.02  # from literature

    # Reynolds number calculation
    mu_water = 0.00126
    ReW = rho_W * V_water * length / mu_water  # Reynolds number

    # Surface area calculation
    S = 2 * np.pi * (D_hull / 2) * (length - (D_hull / 2)) + 2 * np.pi * ((D_hull / 2) ** 2) + A_front

    # Frictional resistance from hull
    def RF_flat(ReW):
        CF_flat = 0.0667 / ((math.log10(ReW) - 2) ** 2)
        S_hull = 2.25 * length * D_hull
        RF_flat = 0.5 * rho_W * S_hull * V_water ** 2 * CF_flat
        return RF_flat

    # Hull form drag
    def Hull_form_drag():
        CF_form = 0.075 / ((math.log10(ReW) - 2) ** 2)
        Kp = eta_hull * ((length / D_hull) ** (-1.7))
        Cp = Kp * CF_form
        RF_form = 0.5 * rho_W * V_water ** 2 * A_front * Cp
        return RF_form

    # Wing+tail skin friction
    def control_surface_friction():
        b = np.sqrt(AR * S)
        c_root = (2 * S) / ((1 + taper_ratio) * b)
        overlapped_area = 0.5 * (2 * c_root) * (b / 2) - (0.5 * b / 2 * (2 * c_root - D_hull))
        new_area = S - overlapped_area
        A_plan = new_area + (0.146 * 3)  # Swept overlapped Wing area + tail area
        R_surface = 0.5 * rho_W * V_water ** 2 * A_plan * Cd_water
        return R_surface

    # Total drag
    def total_drag():
        RF_flat_drag = RF_flat(ReW)
        Hull_form_drag_drag = Hull_form_drag()
        wingtail_friction_drag = control_surface_friction()
        drag = RF_flat_drag + Hull_form_drag_drag + wingtail_friction_drag
        Cd = drag / (0.5 * rho_W * V_water ** 2 * S)
        return drag, Cd, RF_flat_drag, Hull_form_drag_drag, wingtail_friction_drag

    return total_drag()

# Loop over different lengths and slenderness ratios
for i, length in enumerate(lengths):
    for j, slender in enumerate(slender_ratios):
        drag, Cd, RF_flat_drag, Hull_form_drag_drag, wingtail_friction_drag = calculate_drag(length, slender)

        # Store drag value in the array
        drag_values[i, j] = drag

        print(f"Length: {length}, Slenderness Ratio: {slender}")
        print(f"Total Drag: {drag}")
        print(f"Drag Coefficient: {Cd}")
        print(f"Hull Flat Plate Resistance: {RF_flat_drag}")
        print(f"Hull Form Drag Resistance: {Hull_form_drag_drag}")
        print(f"Wing+Tail Skin Friction Resistance: {wingtail_friction_drag}")
        print("-" * 20)

# Plotting heatmap with increased resolution
fig, ax = plt.subplots()
im = ax.imshow(drag_values, cmap='coolwarm', extent=[slender_ratios[0], slender_ratios[-1], lengths[0], lengths[-1]], aspect='auto', interpolation='bilinear')
cbar = fig.colorbar(im)
cbar.set_label('Drag')
ax.set_xlabel('Slenderness Ratio')
ax.set_ylabel('Length')
ax.set_title('Drag Heatmap for Different Lengths and Slenderness Ratios')

# Add diameter annotations to the heatmap
for i in range(len(lengths)):
    for j in range(len(slender_ratios)):
        diameter = lengths[i] / slender_ratios[j]
        ax.text(slender_ratios[j], lengths[i], f'{diameter:.2f}', ha='center', va='center', color='white')

plt.show()


