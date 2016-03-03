__author__ = 'jacob'
import csv
import os
import numpy as np


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
        denny = [[0 for x in range(jmax+2)] for x in range(jmax+2)]
        anggy = [[0 for x in range(jmax+1)] for x in range(jmax+1)]
        with open("fort.2", 'w', newline='') as model_file:
            writer = csv.writer(model_file, delimiter="  ")
            # Write header line
            writer.writerow([str(polytropic_index)] + TODO: Stuff + [str(jmax)] + [str(kmax)])
            # Fortran saves out arrays column first, so first row in file would be the first entry in each row in array
            # Each line is composed of 8 floating point with 22 spaces with 15 after the decimal place, then two spaces

            # Writing density array to fort.2
            # Save the first 8 values of the column to a temp array then write it, repeat until end of denny
            temp_denny = []
            for row in range(len(denny)):
                for column in range(jmax+2):
                    temp_denny.append(denny[column][row])
                    if len(temp_denny) == 8:
                        writer.writerow(temp_denny)
                        temp_denny = []
            # Writing the angular momentum array to fort.2
            # Repeat what was done for the denny array
            temp_anggy = []
            for row in range(len(anggy)):
                for column in range(jmax+1):
                    temp_anggy.append(anggy[column][row])
                    if len(temp_denny) == 8:
                        writer.writerow(temp_anggy)
                        temp_anggy = []
                        #TODO Input coordinates of points into RWI equations and output them with 8 density points per line
                        # TODO Then for specific angular momentum, same thing
    else:
        print('Using normal equations for disk')
