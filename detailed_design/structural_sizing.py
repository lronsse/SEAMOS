"""
File that calculates all stresses on the structure

Stresses calculated: - Cruise: Lift, Drag, Weight, Thrust
                     - Underwater: Pressure, Drag, Thrust, Lift, Weight, Buoyancy
                     - AWT: Impact on nose and wings
                     - WAT: Impulse due to T-O system


Inputs: Geometry: - Wing (planform + airfoil)
                  - Fuselage
                  - Tail
        Loads: see stresses calculated
        Material


Outputs: - Stresses due to all loads
         - Design requirements needed for bearing with loads
         - Deflections (if wanted, let Mathis know)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

g = 9.81
n_points = 1000

class Material:
    def __init__(self, rho, sig_yld, sig_ult, E, G, nu):
        """
        Initiate variable of type material
        :param rho: Density of the material [kg/m^3]
        :type rho: float
        :param sig_yld: Yield strength of the material [N/m^2]
        :type sig_yld: float
        :param sig_ult: Ultimate strength of the material [N/m^2]
        :type sig_ult: float
        :param E: The elastic modules [N/m^2]
        :type E: float
        :param G: Shear modules [N/m^2]
        :type G: float
        :param nu: The Poisson's ratio of the material
        :type nu: float
        """
        self.rho = rho
        self.sig_yld = sig_yld

        self.sig_ult = sig_ult
        self.E = E
        self.G = G
        self.nu = nu

    def mass(self, volume):
        """
        This function will return the mass (in kilograms) of an object given its volume (in meters cubed)

        :param volume: The object's volume [m^3]
        :type volume: float
        :return: The mass of the object [kg]
        :rtype: float
        """
        return self.rho * volume


class Wing:
    def __init__(self, b, cr, ct, sweep_le, spar_rear, spar_front):
        """
        Initiate variable of type planform

        :param b: Wingspan of the planform [m]
        :type b: float
        :param cr: Cord at the root [m]
        :type cr: float
        :param ct: Cord at the tip [m]
        :type ct: float
        :param sweep_le: Leading edge sweep [rad]
        :param spar_rear: The percentage of the cord where the rear spar is
        :type spar_rear: float
        :param spar_front: The percentage of the cord where the front spar is
        :type spar_front: float
        """
        self.b = b
        self.cr = cr
        self.ct = ct
        self.sweep_le = sweep_le
        self.spar_rear = spar_rear
        self.spar_front = spar_front
        self.spar_dif = spar_rear - spar_front
        self.mass = 1

    def area(self):
        """
        This function returns the area (in m^2) of the wing

        :return: Area of the planform [m^2]
        :rtype: float
        """
        return self.b * (self.cr + self.ct) * 0.5

    def chord(self, y):
        """
        This function returns the cord of the planform at any given distance away form the root

        :param y: Distance away from the root [m]
        :type y: float
        :return: Cord length [m]
        :rtype: float
        """
        return self.cr - 2 * y * (self.cr - self.ct) / self.b

    def naca4(self, number, chord_length=1, n=n_points, plot=False):
        # NACA 4 digit airfoil generator
        m = int(number[0]) / 100.0
        p = int(number[1]) / 10.0
        t = int(number[2:]) / 100.0

        x = np.linspace(0, 1, n)
        yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1015 * x ** 4)

        yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))

        xu = x - yt * np.sin(np.arctan(np.gradient(yc, x)))
        xl = x + yt * np.sin(np.arctan(np.gradient(yc, x)))
        yu = yc + yt * np.cos(np.arctan(np.gradient(yc, x)))
        yl = yc - yt * np.cos(np.arctan(np.gradient(yc, x)))

        # Scale the airfoil
        xu *= chord_length
        yu *= chord_length
        xl *= chord_length
        yl *= chord_length

        max_thickness = max(yu - yl)
        chord_length = max(xu) - min(xu)

        if plot==True:
            plt.plot(xu, yu, 'b')
            plt.plot(xl, yl, 'b')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        return xu, yu, xl, yl, max_thickness, chord_length

    def moment_of_inertia(self, chord_length):
        points = self.naca4('2412', chord_length, n_points, False)
        I_x = np.sum(self.mass / n_points * points[1] ** 2) + np.sum(self.mass / n_points * points[3] ** 2)
        I_y = np.sum(self.mass / n_points * points[0] ** 2) + np.sum(self.mass / n_points * points[2] ** 2)
        return I_x, I_y



    def calculate_internal_loads(self, mass_plane_nowings, mass_wings, plot):
        # Input parameters
        # Length of the beam
        L = self.b / 2

        # List of point loads (position, magnitude)
        point_loads = [(0.0, mass_plane_nowings * g), (L / 2, mass_wings * g)]  # downwards positive

        # Non-uniform load as a function of x
        def w(x):
            return -5

        # Positions along the beam
        x = np.linspace(0, L, n_points)

        V = np.zeros_like(x)
        M = np.zeros_like(x)

        # Calculate shear force and bending moment due to point loads
        for pos, mag in point_loads:
            V -= np.where(x < pos, mag, 0)
            M -= np.where(x < pos, mag * (pos - x), 0)  # Change the sign here

        # Add shear force and bending moment due to distributed load
        V -= w(x) * (L - x)
        M += np.array([quad(lambda xi: w(xi) * (L - xi), xi, L)[0] for xi in x])
        if plot:
            plt.figure(figsize=(12, 6))

            plt.subplot(2, 1, 1)
            plt.plot(x, V, label='Shear force')
            plt.xlabel('x')
            plt.ylabel('V(x)')
            plt.title('Shear Force Diagram')
            plt.grid(True)

            plt.subplot(2, 1, 2)
            plt.plot(x, M, label='Bending moment')
            plt.xlabel('x')
            plt.ylabel('M(x)')
            plt.title('Bending Moment Diagram')
            plt.grid(True)

            plt.tight_layout()
            plt.show()
        return V, M



class Fuselage:
    def __int__(self, length, radius, n_fuselages):
        self.length = length
        self.radius = radius
        self.n_fuselages = n_fuselages
        self.volume = self.length * np.pi * self.radius ** 2  # Todo: find more accurate fuselage volume formula

    def longitudinal_stress(self, delta_p, thickness):
        """
        Function that calculates the longitudinal stress in a pressurised cylinder

        :param delta_p: pressure difference between inside and outside fuselage [Pa]
        :param thickness: thickness of the fuselage [m]
        :return: Circumferential stress [Pa]
        """
        return delta_p * self.radius / (2 * thickness)

    def circumferential_stress(self, delta_p, thickness):
        """
        Function that calculates the circumferential stress in a pressurised cylinder

        :param delta_p: pressure difference between inside and outside fuselage [Pa]
        :param thickness: thickness of the fuselage [m]
        :return: Circumferential stress [Pa]
        """
        return delta_p * self.radius / thickness

    def x_cg(self):
        """
        Function calculating the x_cg location of the fuselage

        :return: x-location of centre of gravity [m] of fuselage
        """
        return self.length / 2


wing = Wing(5, 1.3, 1, 0, 0.25, 0.75)
wing.moment_of_inertia(1)
wing.naca4('2412', plot=False)
wing.calculate_internal_loads(16, 5, False)
