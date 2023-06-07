"""
File that calculates all stresses on the structure

Stresses calculated: - Cruise: Lift, Drag, Weight, Thrust, Torsion on wings due to moment
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
from scipy.integrate import quad, simps, cumtrapz

g = 9.81
n_points = 1000
rho_water = 1023.56
airfoil = '2412'




class Material:
    def __init__(self, rho, sig_yld, sig_ult, E, G, nu):  # todo: create dictionary
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
    def __init__(self, b, cr, ct, sweep_le, spar_rear, spar_front, airfoil, thickness, wing_mass):
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
        self.mass = wing_mass
        self.airfoil = airfoil
        self.thickness_wing_sheet = thickness

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

    def naca4(self, chord_length, n=n_points, plot=False):
        # NACA 4 digit airfoil generator
        number = self.airfoil
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

        if plot:
            plt.plot(xu, yu, 'b')
            plt.plot(xl, yl, 'b')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        return xu, yu, xl, yl, max_thickness, chord_length

    def first_moment_area(self, chordlength):
        t = self.thickness_wing_sheet
        xu, yu, xl, yl, dummy1, dummy2 = self.naca4(chord_length=chordlength)

        xm = np.concatenate((xu, xl[::-1]))
        ym = np.concatenate((yu, yl[::-1])) / 2

        # Calculate the length of each differential element along the perimeter
        dl = np.sqrt(np.gradient(xm) ** 2 + np.gradient(ym) ** 2)

        # Calculate the distance from the neutral axis to each differential element
        y = ym - np.mean(ym)

        # Calculate the first moment of area
        Q = simps(y * t * dl, xm)
        return Q

    def moment_of_inertia(self, chord_length):
        points = self.naca4(chord_length, n_points, False)
        I_x = np.sum(self.mass / n_points * points[1] ** 2) + np.sum(self.mass / n_points * points[3] ** 2)
        I_y = np.sum(self.mass / n_points * points[0] ** 2) + np.sum(self.mass / n_points * points[2] ** 2)
        return I_x, I_y

    def second_moment_of_area(self, chordlength, t):
        # Concatenate upper and lower coordinates
        xu, yu, xl, yl, dummy1, dummy2 = self.naca4(chord_length=chordlength)

        x = np.concatenate((xu, np.flip(xl)))
        y = np.concatenate((yu, np.flip(yl)))

        # Calculate length along the wall
        s = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)

        x_centroid = np.mean(x)
        y_centroid = np.mean(y)

        # Calculate distance from the x-axis
        y = (y[:-1] + y[1:]) / 2 - y_centroid

        # Calculate differential area
        dA = t * s

        # Calculate second moment of area
        I = np.sum(y ** 2 * dA)

        return I

    def polar_moment_of_inertia(self, chordlength):
        # Concatenate upper and lower coordinates
        xu, yu, xl, yl, dummy1, dummy2 = self.naca4(chord_length=chordlength)
        t = self.thickness_wing_sheet

        x = np.concatenate((xu, np.flip(xl)))
        y = np.concatenate((yu, np.flip(yl)))

        # Calculate length along the wall
        s = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)

        # Calculate centroid of the airfoil
        x_centroid = np.mean(x)
        y_centroid = np.mean(y)

        # Calculate distance from the centroid
        r = np.sqrt((x[:-1] - x_centroid) ** 2 + (y[:-1] - y_centroid) ** 2)

        # Calculate differential area
        dA = t * s

        # Calculate polar moment of inertia
        J = np.sum(r ** 2 * dA)

        return J

    def calculate_internal_loads(self, mass_plane_nowings, mass_wings, plot):
        # Input parameters
        # Length of the beam
        L = self.b / 2

        # List of point loads (position, magnitude)
        point_loads = [(0.0, mass_plane_nowings * g / 2), (L / 2, mass_wings * g * 0 / 2)]  # downwards positive

        # Non-uniform load as a function of x
        def w(x):
            L = (mass_wings + mass_plane_nowings)*g
            return - (4 * L / (np.pi * self.b) * np.sqrt(1-(2*x/self.b)**2)) #+ mass_wings * g / self.b * x

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

        E = 70e9  # Modulus of elasticity (example value)
        I = self.second_moment_of_area(chord, thickness)  # Moment of inertia
        v = cumtrapz(cumtrapz(M / (E * I), x, initial=0), x, initial=0)

        if plot:
            plt.figure(figsize=(12, 6))

            plt.subplot(2, 2, 1)
            plt.plot(x, V, label='Shear force')
            plt.xlabel('x')
            plt.ylabel('V(x)')
            plt.title('Shear Force Diagram')
            plt.grid(True)

            plt.subplot(2, 2, 2)
            plt.plot(x, M, label='Bending moment')
            plt.xlabel('x')
            plt.ylabel('M(x)')
            plt.title('Bending Moment Diagram')
            plt.grid(True)

            plt.subplot(2, 2, 3)
            plt.plot(x, v, label='Deflection')
            plt.xlabel('x')
            plt.ylabel('v(x)')
            plt.title('Deflection Diagram')
            plt.grid(True)

            plt.subplot(2, 2, 4)
            plt.plot(x, w(x), label='Lift Distribution')
            plt.xlabel('x')
            plt.ylabel('L(x)')
            plt.title('Lift Distribution')
            plt.grid(True)

            plt.tight_layout()
            plt.show()
        return V, M, v

    def calculate_bending_stress(self, bending_moment, moment_of_inertia, airfoil_thickness):
        return bending_moment * (airfoil_thickness / 2) / moment_of_inertia

    def calculate_shear_stress(self, shear_force, second_moment_of_area, thickness, first_moment_area):
        return shear_force * first_moment_area / (second_moment_of_area * thickness)

    def shear_stress_torsion(self, moment_aero, chord_length):
        a = self.naca4(chord_length)[4]
        b = chord_length
        polar_moment_inertia = np.pi * a ** 3 * b ** 3 / (a ** 2 + b ** 2)
        polar_moment_inertia = self.polar_moment_of_inertia(chord_length)
        distance_centroid_farpoint = b / 2
        return moment_aero * distance_centroid_farpoint / polar_moment_inertia

    def wing_mass(self, material):
        # Calculate the perimeter of the airfoil cross-section
        xu, yu, xl, yl, dummy1, dummy2 = self.naca4(self.cr)
        perimeter = (np.sum(np.sqrt(np.diff(xu) ** 2 + np.diff(yu) ** 2)) + np.sum(
            np.sqrt(np.diff(xl) ** 2 + np.diff(yl) ** 2)))/2
        # Calculate the surface area of the wing
        surface_area = perimeter * self.b

        # Calculate the volume of the wing
        volume = surface_area * self.thickness_wing_sheet

        # Calculate the mass of the wing
        mass = material.rho * volume

        return mass


class Fuselage:
    def __init__(self, length, radius, n_fuselages, thickness_fuselage):
        self.length = length
        self.radius = radius
        self.n_fuselages = n_fuselages
        self.volume = self.length * np.pi * self.radius ** 2  # Todo: find more accurate fuselage volume formula
        self.thickness = thickness_fuselage

    def pressure_difference(self, depth, area, p_internal):
        p = rho_water * depth * g
        return abs(p_internal-p)
    
    def longitudinal_stress(self, delta_p):
        """
        Function that calculates the longitudinal stress in a pressurised cylinder
        
        :param delta_p: pressure difference between inside and outside fuselage [Pa]
        :param thickness: thickness of the fuselage [m]
        :return: Circumferential stress [Pa]
        """
        return delta_p * self.radius / (2 * self.thickness)

    def circumferential_stress(self, delta_p):
        """
        Function that calculates the circumferential stress in a pressurised cylinder

        :param delta_p: pressure difference between inside and outside fuselage [Pa]
        :param thickness: thickness of the fuselage [m]
        :return: Circumferential stress [Pa]
        """
        return delta_p * self.radius / self.thickness

    def x_cg(self):
        """
        Function calculating the x_cg location of the fuselage

        :return: x-location of centre of gravity [m] of fuselage
        """
        return self.length / 2

    def mass(self, material):
        volume = self.thickness * 2 * np.pi * self.radius * self.length * self.n_fuselages
        return volume * material.rho


thickness = 1 * 10 ** -3
AR = 5
S = 1.25
b = np.sqrt(AR * S)
chord = S/b
print(b)
print(chord)
wing = Wing(b, 1.3, 1, 0, 0.25, 0.75, airfoil, thickness, 4)
wing.naca4(chord_length=chord, plot=True)
wing.calculate_internal_loads(12, 4, True)
#print(wing.first_moment_area(chord))

print(max(abs(wing.calculate_shear_stress(wing.calculate_internal_loads(12, 4, False)[0], wing.second_moment_of_area(chord, thickness), thickness, wing.first_moment_area(chord)))))
#print(wing.shear_stress_torsion(150, chord))

print(max(abs(wing.calculate_bending_stress(wing.calculate_internal_loads(12, 4, False)[1], wing.second_moment_of_area(chord, thickness), wing.naca4(chord)[-2]))))

cfrp = Material(1600, 120, 180, 70000, 15000, 1.2)

fuselage = Fuselage(1, 0.1, 1, 1*10**-3)
print(f'mass fuselage {fuselage.mass(cfrp)}')
print(fuselage.circumferential_stress(20000))
wing_mass = wing.wing_mass(cfrp)
print('start')
for i in range(2):
    wing = Wing(b, 1.3, 1, 0, 0.25, 0.75, airfoil, thickness, wing_mass)
    wing.calculate_internal_loads(3*wing_mass, wing_mass, True)
    print('shear', max(abs(wing.calculate_shear_stress(wing.calculate_internal_loads(3*wing_mass, wing_mass, False)[0],
                                              wing.second_moment_of_area(chord, thickness), thickness,
                                              wing.first_moment_area(chord)))))

    print('bending', max(abs(wing.calculate_bending_stress(wing.calculate_internal_loads(3*wing_mass, wing_mass, False)[1],
                                                wing.second_moment_of_area(chord, thickness), wing.naca4(chord)[-2]))))
    wing_mass = wing.wing_mass(cfrp)
    print(wing_mass)
