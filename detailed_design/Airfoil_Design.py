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

# ... your existing code ...

# Define the range of alpha for the quasi-linear part of the lift curve
alpha_linear_range = np.arange(0, 10, delta_alpha)  # Adjust this range as needed

# Round the alpha values in the DataFrame index and alpha_linear_range to match the precision of your data
precision = 2  # Adjust this value as needed
data.index = pd.Series(data.index).round(precision).values
alpha_linear_range = alpha_linear_range.round(precision)

# Extract the lift coefficients in the linear range
cl_linear_range = data.loc[alpha_linear_range, 'CL']

# Fit a line to the lift coefficients in the linear range
slope, intercept = np.polyfit(alpha_linear_range, cl_linear_range, 1)

# Print the lift curve slope
print(f'Lift curve slope: {slope} per degree')

# Interpolate the drag coefficient at zero lift
cd0 = np.interp(0, data.index, data['CD'])

# Print the zero-lift drag coefficient
print(f'Zero-lift drag coefficient (Cd0): {cd0}')




