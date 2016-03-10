__author__ = 'jacob'
import csv
import os
import numpy as np
import math
from fortranformat import FortranRecordWriter

'''
FORTRAN equation code: (Might not be correct equations)
               DENNY(j,k) = ((1/SQRT(j*ROF3N))(1 + .25EXP(-((j*ROF3N)* &
                (j*ROF3N))/(2*0.05*0.05)))/ &
 -               (SQRT(2*pi)(.7j)))EXP(-(k*k/(4.2((j*ROF3N)*(j*ROF3N))))
 +               (SQRT(2*pi)(.7j)))EXP(-((k*ZOF3N)*(k*ZOF3N)/ &
 +               (4.2((j*ROF3N)*(j*ROF3N))))
                ANGGY(j,k) = SQRT((j*ROF3N)*2.1*1+1/(j*ROF3N))

c---------------------------------------------
c Filling array with Rossby-Wave model
c G=gravitational constant, 1?
c W= width of the density bump, .05?
c pi constant? if not, add
c
c velocity(radius) : SQRT(r*2.1*G+G/r)
c
c density(radius, z) : ((1/SQRT(r))(1 + .25EXP(-(r*r)/(2W*W)))/(SQRT(2*pi)(.7r)))EXP(-(z*z/(4.2(r*r)))
c  J at dr go to ITYPE 7
c just read first numbers, not read ANGGY or DENNY, two equations above do that
c thene it should work
c L = 1 because it is rotating the grid made by hscf.f
c---------------------------------------------
'''

#TODO Add Length Scale
def gaussian(x, mu, sig):
    '''
    Return the gaussian distribution
    :param x: Generally, r, the radius
    :param mu: Generall r_nought, the center of the bump
    :param sig: delta_r, the width of the bump
    :return: The gaussian distribution, not scaled for the amplitude
    '''
    distribution = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    scaled = (1.0 / (math.sqrt(2 * math.pi) * sig)) * distribution
    return scaled


def p_nought(amplitude, radius, r_nought, delta_r, alpha):
    gaussian_bump = amplitude * gaussian(radius, r_nought, delta_r)
    print("Gaussian Bump: " + str(gaussian_bump))
    b_r = amplitude * gaussian_bump
    print("B(r): " + str(b_r))
    return math.pow(radius, -alpha) * b_r


def big_h(h, r):
    print("H: " + str(h * r))
    return h * r


def p_nought_coefficient(h, z, r):
    print("little h: " + str(h))
    print("Z: " + str(z))
    print("P Nought: " + str((1.0 - (z / big_h(h, r)) ** 2)))
    return (1.0 - (z / big_h(h, r)) ** 2)


def surface_density_profile(amplitude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index, jmin):
    if radius > jmin and z_height <= big_h(h, radius):
        density_at_point = p_nought(amplitude, radius, r_nought, delta_r, alpha=alpha) * (p_nought_coefficient(h, z_height,
                                                                                                               radius)) ** polytropic_index
    else:
        density_at_point = 0.000
    return density_at_point


# Angular Momentum functions
def pressure(density, constant, polytropic_index):
    power = 1 + 1/polytropic_index
    return constant * density**power


def phi(radius, z, mass_star, g):
    coefficient = 1 - (z**2/(2*radius**2))
    gravity_and_mass = (-g*mass_star)/radius
    return gravity_and_mass*coefficient


def gradient_phi_r(radius, z, mass_star, g):
    top = g*mass_star*(2*radius**2-3*z**2)
    bottom = 2*radius**4
    return top/bottom


def gradient_phi_z(radius, z, mass_star, g):
    top = g*mass_star*z
    bottom = radius**3
    return top / bottom


def gradient_pressure_z(ampltiude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index):
    e_term = math.exp(-(radius - r_nought)**2/(2*delta_r**2))
    alpha_term = radius**(-alpha)
    ampltiude_term = -(2*polytropic_index*z_height*ampltiude)/((radius*h)**2)
    polytropic_term = (1-(z_height**2)/((radius*h)**2))**(polytropic_index-1)
    return e_term*alpha_term*ampltiude_term*polytropic_term


def gradient_pressure_r(ampltiude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index):
    exponential_term = -(radius-r_nought)**2/((2*delta_r**2))
    inside_polytropic = (1-(z_height**2)/((radius*h)**2))
    over_h_term = 2*ampltiude*polytropic_index*z_height**2*radius**(-alpha-3)*math.exp(exponential_term)* \
                  inside_polytropic**(polytropic_index-1)
    final_h_term = over_h_term/(h**2)
    no_denom_term = alpha*ampltiude*radius**(-alpha-1)*math.exp(exponential_term)*inside_polytropic**polytropic_index
    over_delta_r_term = ampltiude*radius**(-alpha)*(radius-r_nought)*math.exp(exponential_term)*inside_polytropic \
                                                                                                **polytropic_index
    final_delta_r_term = over_delta_r_term/(delta_r**2)
    final_term = final_h_term - no_denom_term - final_delta_r_term
    return final_term


def velocity_field(ampltiude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index, density, mass_star, g, jmin):
    pressure_gradient_grid = []
    gravitational_gradient_grid = []
    pressure_gradient_grid.append((-(1/density)*gradient_pressure_z(ampltiude, radius, r_nought, delta_r, h, z_height,
                                                                    alpha, polytropic_index),
                                   -(1/density)*gradient_pressure_r(ampltiude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index)))
    gravitational_gradient_grid.append((-gradient_phi_r(radius, mass_star=mass_star, z=z_height, g=g),
                                        -gradient_phi_z(radius, mass_star=mass_star, z=z_height, g=g)))


def generate_fort_2(polytropic_index, model, jmax, kmax, jout, kout, log_central_density, iteration, mass_star, xcut,
                    xwidth, xnorm, type):
    """
    Generate fort.2 model file for CHYMERA Code
    :param polytropic_index: Polytropic index of star
    :param model: The star/disk model, 100 for star/disk, -2.0 for disk with point mass star
    :param jmax: radial grid size
    :param kmax: vertical grid size
    :param jout: radial size of star+disk. negative for just the disk
    :param kout: polar radius of star
    :param log_central_density: log of central density
    :param iteration: iteration parameter for model
    :param mass_star: point mass in center of star
    :param xcut: width of star in model
    :param xwidth: percentage of mass in the disk
    :param xnorm: TODO: Figure these x things out
    :param type: Currently either 'RWI' for Rossby Wave Instability, anything else for general disk model
    :return: A model file for input into the CHYMERA computational fluid dynamics code to run
    """

    if type == 'RWI':
        print('Using Rossby Wave equations for disk')
        denny = [[0 for x in range(jmax + 2)] for x in range(jmax + 2)]
        anggy = [[0 for x in range(jmax + 1)] for x in range(jmax + 1)]

        # Get the density array
        for row in range(jmax + 2):
            for column in range(jmax + 2):
                denny[row][column] = surface_density_profile(1.4, row + 1, 100, 20, 0.14, column + 1, 0.5, 1.5, 2)

        # Get angular momentum array
        for row in range(jmax + 1):
            for column in range(jmax + 1):
                anggy[row][column] = 1
        with open("temp", 'w') as model_file:
            fortran_writer = FortranRecordWriter('8(1PE10.3,2X)')
            # Write header line
            # Fortran saves out arrays column first, so first row in file would be the first entry in each row in array
            # Each line is composed of 8 floating point with 22 spaces with 15 after the decimal place, then two spaces
            # TODO Use Equations for Density and Angaulaer Momentum
            # Writing density array to fort.2
            # Save the first 8 values of the column to a temp array then write it, repeat until end of denny
            temp_denny = []
            for row in range(jmax + 2):
                for column in range(jmax + 2):
                    temp_denny.append(denny[column][row])
                    ''' if len(temp_denny) == 8:
                         writer.writerow(temp_denny)
                         temp_denny = []
                     if row == jmax + 1 and column == jmax + 1:
                         writer.writerow(temp_denny)
                         temp_denny = []'''
            output_text = fortran_writer.write(temp_denny)
            model_file.write(output_text)
            # Writing the angular momentum array to fort.2
            # Repeat what was done for the denny array
            temp_anggy = []
            for row in range(jmax + 1):
                for column in range(jmax + 1):
                    temp_anggy.append(anggy[column][row])
                    if len(temp_denny) == 8:
                        writer.writerow(temp_anggy)
                        temp_anggy = []
                        # TODO Input coordinates of points into RWI equations and output them with 8 density points per line
                        # TODO Then for specific angular momentum, same thing
    else:
        print('Using normal equations for disk')


generate_fort_2(1.5, 2.0, 256, 256, 10, 10, -10, 50, 1, 30, 30, 30, 'RWI')
