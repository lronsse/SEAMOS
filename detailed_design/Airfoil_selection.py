import subprocess
import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import aerosandbox as asb


def run_xfoil(naca_number, mach, reynolds, alpha_start, alpha_end, alpha_step):
    # Define the XFOIL command
    xfoil_command = r"C:\Users\mathi\OneDrive\Documenten\School\TU Delft\BSC3\Q4\DSE\Programs\xfoil.exe"  # Replace with your actual path

    # Define the XFOIL commands
    xfoil_commands = [f'NACA {naca_number}', 'OPER', f'VISC {reynolds}', f'MACH {mach}', 'ITER 200', f'ASEQ {alpha_start} {alpha_end} {alpha_step}', 'CPWR polar.txt',
                      'QUIT']

    # Start Xfoil
    process = subprocess.Popen(xfoil_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)

    # Send the commands to Xfoil one at a time
    for command in xfoil_commands:
        process.stdin.write(command + '\n')

    output, _ = process.communicate()

    # Print the output
    print(output)

    # Parse the output to get the polar data
    lines = output.split('\n')
    polar_data = []
    for line in lines:
        if 'a =' in line and 'CL =' in line:
            data = line.split()
            alpha = float(data[data.index('=')+1])
            CL = float(data[data.index('CL =')+2])
            Cm = float(data[data.index('Cm =')+1])
            CD = float(data[data.index('CD')+3])
            CDf = float(data[data.index('CDf')+2])
            CDp = float(data[data.index('CDp')+2])
            polar_data.append([alpha, CL, CD, CDp, Cm, CDf])

    polar = pd.DataFrame(polar_data, columns=['alpha', 'CL', 'CD', 'CDp', 'CM', 'CDf'])

    # Return the polar data
    return polar


# List of airfoils to compare
airfoils = ['0012', '2412', '4412']

# Run Xfoil for each airfoil
for airfoil in airfoils:
    polar = run_xfoil(airfoil, mach=0.1, reynolds=1e6, alpha_start=-10, alpha_end=10, alpha_step=1)
    if polar is not None:
        plt.plot(polar['alpha'], polar['CL'], label=f'NACA {airfoil}')

# Show the plot
plt.xlabel('Angle of Attack [deg]')
plt.ylabel('Lift Coefficient')
plt.legend()
plt.grid(True)
plt.show()
