import subprocess
import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D




def download_airfoil(airfoil_number):
    url = f'http://airfoiltools.com/airfoil/seligdatfile?airfoil=naca{airfoil_number}-il'
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception if the GET request was unsuccessful
    except requests.exceptions.HTTPError:
        # If the first URL didn't work, try with the "n" prefix
        airfoil_number = f'n{airfoil_number}'
        url = f'http://airfoiltools.com/airfoil/seligdatfile?airfoil={airfoil_number}-il'
        response = requests.get(url)
        response.raise_for_status()

    with open(f'{airfoil_number}.dat', 'w') as f:
        f.write(response.text)
    return airfoil_number

def run_xfoil(naca_number, mach, reynolds, alpha_start, alpha_end, alpha_step):
    # Define the XFOIL command
    xfoil_command = r"C:\Users\mathi\OneDrive\Documenten\School\TU Delft\BSC3\Q4\DSE\Programs\xfoil.exe"  # Replace with your actual path

    # Construct the airfoil name and filename from the NACA number
    airfoil_name = f'naca{naca_number}'
    airfoil_file = f'{airfoil_name}.dat'

    # Download the airfoil file if it doesn't exist
    if not os.path.isfile(airfoil_file):
        print(f'Airfoil file {airfoil_file} not found. Downloading...')
        airfoil_name = download_airfoil(naca_number)
        airfoil_file = f'{airfoil_name}.dat'

    # Define the XFOIL commands
    xfoil_commands = f'''
    LOAD {airfoil_file}
    PANEL
    OPER
    VISC {reynolds}
    MACH {mach}
    ITER 200
    PACC
    polarnaca{naca_number}_{alpha_step}.txt

    ASEQ {alpha_start} {alpha_end} {alpha_step}
    PACC
    QUIT
    '''

    # Run XFOIL with the defined commands
    try:
        process = subprocess.Popen(xfoil_command, stdin=subprocess.PIPE, text=True)
        process.communicate(xfoil_commands)
    except subprocess.CalledProcessError as e:
        print(f'XFOIL process failed with error: {e}')


def read_data(filename):
    # Define the column names
    col_names = ["alpha", "CL", "CD", "CDp", "Cm", "Top Xtr", "Bot Xtr", "Cpmin", "Chinge", "XCp"]

    # Read the data from the file
    data = pd.read_csv(filename, sep="\s+", skiprows=12, names=col_names)

    # Return the data
    return data

def get_coefficients(alpha, data):
    '''
    if alpha < data.index.min() or alpha > data.index.max():
        return None
    '''

    lift_coefficient = data.loc[alpha, 'CL']
    drag_coefficient = data.loc[alpha, 'CD']
    moment_coefficient = data.loc[alpha, 'Cm']
    return lift_coefficient, drag_coefficient, moment_coefficient

def planform(mach, wing_area, aspect_ratio):
    if mach < 0.7:
        sweep_quarter_chord = np.degrees(np.arccos(1))
    else:
        sweep_quarter_chord = np.degrees(np.arccos(0.75 * 0.935 / (mach + 0.03)))
    taper_ratio = 0.2 * (2 - sweep_quarter_chord * np.pi / 180)
    wing_span = np.sqrt(wing_area * aspect_ratio)
    root_chord = 2 * wing_area / ((1 + taper_ratio) * wing_span)
    tip_chord = taper_ratio * root_chord
    return sweep_quarter_chord, taper_ratio, wing_span, root_chord, tip_chord

def plot_planform(sweep_quarter_chord, taper_ratio, wing_span, root_chord, tip_chord):
    n_points = 100
    y = np.linspace(0, wing_span / 2, n_points)
    sweep_quarter_chord_rad = np.radians(sweep_quarter_chord)
    c_r = np.linspace(0, root_chord, n_points)
    x_quarter_root = root_chord / 4
    x_quarter_tip = x_quarter_root + wing_span / 2 * np.tan(sweep_quarter_chord_rad)
    x_quarter = x_quarter_root + y * np.tan(sweep_quarter_chord_rad)
    le_tip = x_quarter_tip - 0.25 * tip_chord
    te_tip = x_quarter_tip + 0.75 * tip_chord
    c_t = np.linspace(le_tip, te_tip, n_points)
    leading_edge = le_tip / (wing_span / 2) * y
    trailing_edge = root_chord - (root_chord - te_tip) / (wing_span / 2) * y

    mean_aerodynamic_chord = root_chord * 2 / 3 * ((1 + taper_ratio + taper_ratio ** 2) / (1 + taper_ratio))
    y_lemac = wing_span / 2 * (root_chord - mean_aerodynamic_chord) / (root_chord - tip_chord)
    leading_edge_mean_aerodynamic_chord = le_tip / (wing_span / 2) * y_lemac
    x_mean_aerodynamic_chord = np.linspace(leading_edge_mean_aerodynamic_chord, leading_edge_mean_aerodynamic_chord + mean_aerodynamic_chord, n_points)


    plt.plot(c_r, np.zeros(n_points), color='black')
    plt.plot(c_t, np.ones(n_points) * wing_span / 2, color='black')
    plt.plot(leading_edge, y, color='black')
    plt.plot(trailing_edge, y, color= 'black')
    plt.plot(x_quarter, y, 'r--', label='Quarter Chord')
    plt.plot(x_mean_aerodynamic_chord, np.ones(n_points) * y_lemac, 'r', label='Mean Aerodynamic Chord')
    plt.plot()
    plt.xlabel('Chordwise position (m)')
    plt.ylabel('Spanwise position (m)')
    plt.title('Wing Planform')
    plt.legend()
    plt.axis('equal')
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    naca_number = '2412'
    mach = 0.1
    reynolds = 500000

    alphamin = -5
    alphamax = 25
    delta_alpha = 0.2

    filename = f'polarnaca{naca_number}_{delta_alpha}.txt'



    alpha = np.arange(alphamin+1, alphamax-1, delta_alpha)

    if not os.path.isfile(filename):
        run_xfoil(naca_number, mach, reynolds, alphamin, alphamax, delta_alpha)


    data = read_data(filename)
    data = data.set_index('alpha')
    data = data.interpolate('linear')

    plt.plot(data.index, data['CL'])
    plt.title(f'CL alpha for naca{naca_number}')
    plt.show()

    wing = planform(0.1, 1.25, 5)
    plot_planform(wing[0], wing[1], wing[2], wing[3], wing[4])
    print(wing)




