__author__ = 'jacob'
import csv
import os


'''
FORTRAN equation code: (Might not be correct equations)
               DENNY(j,k) = ((1/SQRT(j*ROF3N))(1 + .25EXP(-((j*ROF3N)* &
                (j*ROF3N))/(2*0.05*0.05)))/ &
 -               (SQRT(2*pi)(.7j)))EXP(-(k*k/(4.2((j*ROF3N)*(j*ROF3N))))
 +               (SQRT(2*pi)(.7j)))EXP(-((k*ZOF3N)*(k*ZOF3N)/ &
 +               (4.2((j*ROF3N)*(j*ROF3N))))
                ANGGY(j,k) = SQRT((j*ROF3N)*2.1*1+1/(j*ROF3N))
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
        with open("fort.2", 'w', newline='') as model_file:
            writer = csv.writer(model_file, delimiter=" ")
            writer.writerow([str(polytropic_index)] + TODO: Stuff + [str(jmax)] + [str(kmax)])
            #TODO INput coordinates of points into RWI equations and output them with 8 density points per line
            # TODO Then for specific angular momentum, same thing
    else:
        print('Using normal equations for disk')
