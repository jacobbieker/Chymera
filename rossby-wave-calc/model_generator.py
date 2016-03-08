__author__ = 'jacob'
import csv
import os
import numpy as np
import math

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
        with open("fort.2", 'w') as model_file:
            writer = csv.writer(model_file, delimiter=" ")
            # Write header line
            writer.writerow([str(polytropic_index)] + [str(jmax)] + [str(kmax)])
            # Fortran saves out arrays column first, so first row in file would be the first entry in each row in array
            # Each line is composed of 8 floating point with 22 spaces with 15 after the decimal place, then two spaces
            # TODO Use Equations for Density and Angaulaer Momentum
            # Writing density array to fort.2
            # Save the first 8 values of the column to a temp array then write it, repeat until end of denny
            temp_denny = []
            for row in range(jmax + 2):
                for column in range(jmax + 2):
                    temp_denny.append(denny[column][row])
                    if len(temp_denny) == 8:
                        writer.writerow(temp_denny)
                        temp_denny = []
                    if row == jmax + 2 and column == jmax + 2:
                        writer.writerow(temp_denny)
                        temp_denny = []
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
