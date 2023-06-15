import subprocess
import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import aerosandbox as asb


"""
File getting aerodynamic airfoil data

Inputs: airfoil geometry (currently only naca 4 digit)
        Reynolds number
        Mach number
        
Outputs: Airfoil Polars (cl, cd, cm vs alpha, cl vs cd)
         Airfoil stability coefficients
         
# todo: adapt code for different airfoils, try to implement wing as well (maybe use aeropy?)
"""


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



    # Define the XFOIL commands
    xfoil_commands = f'''
    NACA {naca_number}
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






naca_number = '2412'
mach = 0.1
reynolds = 500000

alphamin = -5
alphamax = 20
delta_alpha = 0.5



# List of airfoils to compare
airfoils = ['0412', '2412', '4412', '6412', '8412']

# Create a dictionary to store the data for each airfoil
airfoil_data = {}

# Run XFOIL for each airfoil and store the data
for airfoil in airfoils:
    filename = f'polarnaca{airfoil}_{delta_alpha}.txt'
    if not os.path.isfile(filename):
        run_xfoil(airfoil, mach, reynolds, alphamin, alphamax, delta_alpha)
    data = read_data(filename)
    data = data.set_index('alpha')
    data = data.interpolate('linear')
    airfoil_data[airfoil] = data
    # Delete the files
    os.remove(filename)
    airfoil_name = f'naca{airfoil}'
    #os.remove(f'{airfoil_name}.dat')


# Plot the aerodynamic coefficients for each airfoil
for airfoil, data in airfoil_data.items():
    #plt.figure(figsize=(10, 6))
    #plt.subplot(3, 2, 1)
    plt.plot(data.index, data['CL'], label=f'NACA {airfoil}')
plt.xlabel('Alpha (degrees)')
plt.ylabel('Lift Coefficient (CL)')
plt.title('Lift Coefficient vs Alpha')
plt.legend()
plt.grid(True)
plt.show()

for airfoil, data in airfoil_data.items():
    #plt.figure(figsize=(10, 6))
    #plt.subplot(3, 2, 2)
    plt.plot(data.index, data['CD'], label=f'NACA {airfoil}')
plt.xlabel('Alpha (degrees)')
plt.ylabel('Drag Coefficient (CD)')
plt.title('Drag Coefficient vs Alpha')
plt.legend()
plt.grid(True)
plt.show()


for airfoil, data in airfoil_data.items():
    #plt.figure(figsize=(10, 6))
    #plt.subplot(3, 2, 3)
    plt.plot(data.index, data['Cm'], label=f'NACA {airfoil}')
plt.xlabel('Alpha (degrees)')
plt.ylabel('Moment Coefficient (Cm)')
plt.title('Moment Coefficient vs Alpha')
plt.legend()
plt.grid(True)
plt.show()

for airfoil, data in airfoil_data.items():
    #plt.figure(figsize=(10, 6))
    #plt.subplot(3, 2, 4)
    plt.plot(data['CD'], data['CL'], label=f'NACA {airfoil}')
plt.xlabel('Drag Coefficient (CD)')
plt.ylabel('Lift Coefficient (CL)')
plt.title('Lift Coefficient vs Drag Coefficient')
plt.legend()
plt.grid(True)
plt.show()

for airfoil, data in airfoil_data.items():
    #plt.figure(figsize=(10, 6))
    #plt.subplot(3, 2, 5)
    plt.plot(data.index, data['CL'] / data['CD'], label=f'NACA {airfoil}')
plt.xlabel('Alpha (degrees)')
plt.ylabel('Lift over Drag (CL/CD)')
plt.title('Lift over Drag vs Alpha')
plt.legend()
plt.grid(True)
plt.show()




# Delete the files
#os.remove(filename)
#airfoil_name = f'naca{naca_number}'
#os.remove(f'{airfoil_name}.dat')
