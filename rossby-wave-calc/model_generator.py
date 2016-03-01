__author__ = 'jacob'

def generate_fort_2(polytropic_index, model, jmax, kmax, jout, kout, log_central_density, iteration, mass_star, xcut, xwidth, xnorm):
    '''

    :param polytropic_index: Polytropic index of star
    :param model: The star/disk model, 100 for star/disk, -2.0 for disk with point mass star
    :param jmax: radial grid size
    :param kmax: vertical grid size
    :param jout: radial size of star+disk
    :param kout: polar radius of star
    :param log_central_density: log of central density
    :param iteration: iteration parameter for model
    :param mass_star: point mass in center of star
    :param xcut: width of star in model
    :param xwidth: percentage of mass in the disk
    :param xnorm: TODO: Figure these x things out
    :return: A model file for input into the CHYMERA computational fluid dynamics code to run
    '''