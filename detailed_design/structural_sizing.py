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
         - Mass of wing
         - Mass of fuselage
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d

g = 9.81
n_points = 1000
rho_water = 1023.56
airfoil = '2412'


class Material:
    def __init__(self, rho, sig_yld, sig_ult, E, G, nu):  # todo: create dictionary
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
    def __init__(self, wing_area, aspect_ratio, mach, airfoil, thickness, material, plane_mass):
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
        self.tail_area = None
        self.le_tip = None
        self.first_moment_of_area = None
        self.plane_mass = plane_mass
        self.chord_array = None
        self.bending_moment = None
        self.shear_force = None
        self.tip_chord = None
        self.root_chord = None
        self.wing_span = None
        self.sweep_quarter_chord = None
        self.taper_ratio = None
        self.wing_area = wing_area
        self.aspect_ratio = aspect_ratio
        self.mass = None
        self.mach = mach
        self.airfoil = airfoil
        self.thickness_wing_sheet = thickness
        self.material = material

    def planform(self):
        if self.mach < 0.7:
            self.sweep_quarter_chord = np.degrees(np.arccos(1))
        else:
            self.sweep_quarter_chord = np.degrees(np.arccos(0.75 * 0.935 / (self.mach + 0.03)))
        self.taper_ratio = 0.2 * (2 - self.sweep_quarter_chord * np.pi / 180)
        self.wing_span = np.sqrt(self.wing_area * self.aspect_ratio)
        self.root_chord = 2 * self.wing_area / ((1 + self.taper_ratio) * self.wing_span)
        self.tip_chord = self.taper_ratio * self.root_chord
        return

    def plot_planform(self):
        n_points = 100
        y = np.linspace(0, self.wing_span / 2, n_points)
        sweep_quarter_chord_rad = np.radians(self.sweep_quarter_chord)
        c_r = np.linspace(0, self.root_chord, n_points)
        x_quarter_root = self.root_chord / 4
        x_quarter_tip = x_quarter_root + self.wing_span / 2 * np.tan(sweep_quarter_chord_rad)
        x_quarter = x_quarter_root + y * np.tan(sweep_quarter_chord_rad)
        self.le_tip = x_quarter_tip - 0.25 * self.tip_chord
        te_tip = x_quarter_tip + 0.75 * self.tip_chord
        c_t = np.linspace(self.le_tip, te_tip, n_points)
        leading_edge = self.le_tip / (self.wing_span / 2) * y
        trailing_edge = self.root_chord - (self.root_chord - te_tip) / (self.wing_span / 2) * y

        self.mean_aerodynamic_chord = self.root_chord * 2 / 3 * ((1 + self.taper_ratio + self.taper_ratio ** 2) / (1 + self.taper_ratio))
        y_lemac = self.wing_span / 2 * (self.root_chord - self.mean_aerodynamic_chord) / (self.root_chord - self.tip_chord)
        leading_edge_mean_aerodynamic_chord = self.le_tip / (self.wing_span / 2) * y_lemac
        x_mean_aerodynamic_chord = np.linspace(leading_edge_mean_aerodynamic_chord, leading_edge_mean_aerodynamic_chord + self.mean_aerodynamic_chord, n_points)

        plt.plot(c_r, np.zeros(n_points), color='black')
        plt.plot(c_t, np.ones(n_points) * self.wing_span / 2, color='black')
        plt.plot(leading_edge, y, color='black')
        plt.plot(trailing_edge, y, color='black')
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

    def chord(self, y):
        """
        This function returns the cord of the planform at any given distance away form the root

        :param y: Distance away from the root [m]
        :type y: float
        :return: Cord length [m]
        :rtype: float
        """
        self.planform()
        return self.root_chord - 2 * y * (self.root_chord - self.tip_chord) / self.wing_span

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
    '''
    def moment_of_inertia(self, chord_length):
        points = self.naca4(chord_length, n_points, False)
        I_x = np.sum(self.mass / n_points * points[1] ** 2) + np.sum(self.mass / n_points * points[3] ** 2)  # todo: improve mass estimate
        I_y = np.sum(self.mass / n_points * points[0] ** 2) + np.sum(self.mass / n_points * points[2] ** 2)
        return I_x, I_y
    '''

    def second_moment_area(self, chordlengths):
        I_array = []

        for chordlength in chordlengths:
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
            dA = self.thickness_wing_sheet * s

            # Calculate second moment of area
            I = np.sum(y ** 2 * dA)

            I_array.append(I)

        return np.array(I_array)


    def polar_moment_inertia(self):
        J_array = []
        t = self.thickness_wing_sheet

        for chordlength in self.chord_array:
            # Concatenate upper and lower coordinates
            xu, yu, xl, yl, dummy1, dummy2 = self.naca4(chord_length=chordlength)

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

            J_array.append(J)

        return np.array(J_array)

    def calculate_internal_loads(self, mass_plane_nowings, mass_wings, wing_mass_function, plot):  # todo: adapt to new planform
        # Input parameters
        # Length of the beam
        L = self.wing_span / 2

        # List of point loads (position, magnitude)
        point_loads = [(0.0, mass_plane_nowings * g / 2), (L / 2, mass_wings * g *0 / 2)]  # downwards positive

        # Non-uniform load as a function of x
        def w(x):
            lift = (mass_wings + mass_plane_nowings)*g
            return - (4 * lift / (np.pi * self.wing_span) * np.sqrt(1-(2*x/self.wing_span) ** 2)) + wing_mass_function(x) * g

        # Positions along the beam
        x = np.linspace(0, L, n_points)

        self.shear_force = np.zeros_like(x)
        self.bending_moment = np.zeros_like(x)

        # Calculate shear force and bending moment due to point loads
        for pos, mag in point_loads:
            self.shear_force -= np.where(x < pos, mag, 0)
            self.bending_moment -= np.where(x < pos, mag * (pos - x), 0)  # Change the sign here

        # Add shear force and bending moment due to distributed load
        self.shear_force -= w(x) * (L - x)
        self.bending_moment += np.array([quad(lambda xi: w(xi) * (L - xi), xi, L)[0] for xi in x])
        '''
        E = 70e9  # Modulus of elasticity (example value)
        I = self.second_moment_of_area(chord, thickness)  # Moment of inertia
        v = cumtrapz(cumtrapz(M / (E * I), x, initial=0), x, initial=0)
        '''
        if plot:
            plt.figure(figsize=(12, 6))

            plt.subplot(2, 2, 1)
            plt.plot(x, self.shear_force, label='Shear force')
            plt.xlabel('x')
            plt.ylabel('V(x)')
            plt.title('Shear Force Diagram')
            plt.grid(True)

            plt.subplot(2, 2, 2)
            plt.plot(x, self.bending_moment, label='Bending moment')
            plt.xlabel('x')
            plt.ylabel('M(x)')
            plt.title('Bending Moment Diagram')
            plt.grid(True)
            '''
            plt.subplot(2, 2, 3)
            plt.plot(x, v, label='Deflection')
            plt.xlabel('x')
            plt.ylabel('v(x)')
            plt.title('Deflection Diagram')
            plt.grid(True)
            '''
            plt.subplot(2, 2, 4)
            plt.plot(x, w(x), label='Lift Distribution')
            plt.xlabel('x')
            plt.ylabel('L(x)')
            plt.title('Lift Distribution')
            plt.grid(True)

            plt.tight_layout()
            plt.show()
        return self.shear_force, self.bending_moment

    def calculate_bending_stress(self, second_moment_of_area, airfoil_thickness):
        return self.bending_moment * (airfoil_thickness / 2) / second_moment_of_area

    def calculate_shear_stress(self, second_moment_of_area, first_moment_area):
        return self.shear_force * first_moment_area / (second_moment_of_area * self.thickness_wing_sheet)

    def shear_stress_torsion(self, moment_aero):
        shear_stress_torsion_array = []
        for chord_length in self.chord_array:
            shear_stress_torsion_array.append(moment_aero * chord_length / 2 / self.polar_moment_of_inertia)

        return shear_stress_torsion_array

    def calc_wing_mass(self, chordlengths, material):
        # Calculate the differences between consecutive x and y values
        mass_array = []
        volume = 0
        dL = self.wing_span / 2 / n_points

        for chordlength in chordlengths:
            xu, yu, xl, yl, dummy1, dummy2 = self.naca4(chord_length=chordlength)

            x_values = np.concatenate((xu, np.flip(xl)))
            y_values = np.concatenate((yu, np.flip(yl)))

            dx = np.diff(x_values)
            dy = np.diff(y_values)

            # Calculate the distances between each pair of consecutive points
            distances = np.sqrt(dx ** 2 + dy ** 2)

            # Sum the distances to get the perimeter
            perimeter = np.sum(distances)
            area = perimeter * thickness
            d_volume = area * dL
            d_mass = d_volume * material.rho
            mass_array.append(d_mass)
            volume = volume + d_volume

        def interpolate_masses(wingspan_values, mass_values):
            # Create the interpolation function
            f = interp1d(wingspan_values, mass_values, kind='cubic')

            # Return the interpolation function
            return f

        wing_mass_function = interpolate_masses(np.linspace(0, self.wing_span / 2, n_points), mass_array)

        self.mass = np.sum(mass_array) * 2

        return self.mass, wing_mass_function

    def wing_total_volume(self, chordlengths):
        volume_array = []

        for chordlength in chordlengths:
            xu, yu, xl, yl, dummy1, dummy2 = self.naca4(chord_length=chordlength, plot=False)

            d_volume = float(np.trapz(yu-yl, xu))
            volume_array.append(d_volume)


        self.total_volume = np.average(volume_array) * self.wing_span
        return self.total_volume

    def tail_sizing(self, tail_arm, tail_taper_ratio):
        Vv = 0.03
        Vh = 0.35
        AR_tail = (2 / 3) * self.aspect_ratio
        Sh = (0.85 * Vh * self.mean_aerodynamic_chord * self.wing_area) / tail_arm
        Sv = (0.85 * Vv * self.wing_span * self.wing_area) / tail_arm
        S_projected_v = 0.33 * Sv
        S_projected_h = Sv - S_projected_v
        tail_anhedral = np.degrees(np.arctan(np.sqrt(S_projected_v / Sh)))
        self.tail_area = 0.5 * (Sh / (np.cos(np.radians(tail_anhedral))) ** 2)
        self.tail_span = np.sqrt(AR_tail * self.tail_area)
        self.tail_root_chord = (2 * self.tail_area) / ((1 + tail_taper_ratio) * self.tail_span)  #
        self.tail_tip_chord = self.tail_root_chord * tail_taper_ratio  #
        self.tail_mac = (2 / 3) * (self.tail_root_chord) * (1 + tail_taper_ratio + tail_taper_ratio ** 2) / (
                    1 + tail_taper_ratio)  #
        self.tail_qc_sweep = np.degrees(np.arctan(
            (((np.tan(0)) - (4 / AR_tail) * ((-75 / 100) * ((1 - tail_taper_ratio) / (1 + tail_taper_ratio)))))))  #
        return self.tail_area, self.tail_span, self.tail_root_chord, self.tail_tip_chord, self.tail_mac, self.tail_qc_sweep, tail_anhedral

    def plot_tail(self):
        n_points = 100
        y = np.linspace(0, self.tail_span / 2, n_points)
        sweep_quarter_chord_rad = np.radians(self.tail_qc_sweep)
        c_r = np.linspace(0, self.tail_root_chord, n_points)
        x_quarter_root = self.tail_root_chord / 4
        x_quarter_tip = x_quarter_root + self.tail_span / 2 * np.tan(sweep_quarter_chord_rad)
        x_quarter = x_quarter_root + y * np.tan(sweep_quarter_chord_rad)
        le_tip = x_quarter_tip - 0.25 * self.tail_tip_chord
        te_tip = x_quarter_tip + 0.75 * self.tail_tip_chord
        c_t = np.linspace(le_tip, te_tip, n_points)
        leading_edge = le_tip / (self.tail_span / 2) * y
        trailing_edge = self.tail_root_chord - (self.tail_root_chord - te_tip) / (self.tail_span / 2) * y

        y_lemac = self.tail_span / 2 * (self.tail_root_chord - self.tail_mac) / (self.tail_root_chord - self.tail_tip_chord)
        leading_edge_mean_aerodynamic_chord = le_tip / (self.tail_span / 2) * y_lemac
        x_mean_aerodynamic_chord = np.linspace(leading_edge_mean_aerodynamic_chord, leading_edge_mean_aerodynamic_chord + self.tail_mac, n_points)

        plt.plot(c_r, np.zeros(n_points), color='black')
        plt.plot(c_t, np.ones(n_points) * self.tail_span / 2, color='black')
        plt.plot(leading_edge, y, color='black')
        plt.plot(trailing_edge, y, color='black')
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


    def wing_main(self, plot=False):
        self.planform()
        self.plot_planform()
        self.chord_array = np.linspace(self.root_chord, self.tip_chord, n_points)
        self.mass = self.calc_wing_mass(self.chord_array, self.material)[0]

        wing_mass_function = wing.calc_wing_mass(self.chord_array, alu)[1]

        self.calculate_internal_loads(self.plane_mass - self.mass, self.mass, wing_mass_function, plot)

        self.first_moment_of_area = self.first_moment_area(self.chord_array)
        self.second_moment_of_area = self.second_moment_area(self.chord_array)
        self.polar_moment_of_inertia = self.polar_moment_inertia()

        self.normal_stress_bending = self.calculate_bending_stress(self.second_moment_of_area, self.naca4(self.chord_array)[4])

        self.shear_stress_torsion = self.shear_stress_torsion(moment)
        self.shear_stress_shear = self.calculate_shear_stress(self.second_moment_of_area, self.first_moment_of_area)
        self.shear_stress = self.shear_stress_shear + self.shear_stress_torsion
        self.total_volume = self.wing_total_volume(self.chord_array)
        self.tail_sizing(0.75, 0.5)
        self.plot_tail()
        return



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

AR = 12
S = 0.7
mach = 0.1
moment = 150
alu = Material(1600, 180, 250, 70, 70, 1.2)



wing = Wing(S, AR, mach, airfoil, thickness, alu, 16)
wing.wing_main(True)
print(wing.tip_chord)
print(wing.root_chord)
print(wing.wing_span)
print(wing.mass)
print(wing.total_volume)
print(wing.taper_ratio)
print(wing.le_tip)
print(wing.mean_aerodynamic_chord)
print(wing.chord(wing.wing_span / 2 - 0.868))
print(f'Tail: Span ({wing.tail_span}), Area ({wing.tail_area}), ')

