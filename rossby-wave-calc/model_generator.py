__author__ = 'jacob'
import csv
import os
import numpy as np
import math
from fortranformat import FortranRecordWriter, FortranRecordReader

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

R∈[0.4,1.6]r0,Z∈[0,Zs] =[0,0.9], andα= 0.5 for the power-law part of the surfacedensity profile. The bump radius,
amplitude and widthare set to r0= 1,A= 1.4 and ∆r= 0.05r0,
'''

# RWI constants
r_nought = 100
delta_r = 0.05*r_nought
amplitude = 1.4
alpha = 0.5

# Constants taken from papers, etc.
h = 0.14
polytropic_index = 1.5
g = 1
mass_star = 1
r_size = 256
z_size = 256
jout = 0.4*r_nought
kout = 10
xnorm = 30
xcut = 30
xwidth = 30
iteration = 50
lcd = -10

'''
Potential new equations:

#Value of density_profile at r_nought
def density_proflie_nought(amplitude):
  return #Integrate from -inf to inf rho(r_nought,z) dz

#Density profile to be used in later equations
def density_profile(r, r_nought, amplitude, delta_r, density_profile_nought, alpha):
  return (density_profile_nought) * (r/r_nought)**(-alpha) * (1+(amplitude-1)*exp((r-r_nought)**2/(2*delta_r**2)))

#Value of big_h at r_nought  
def big_h_nought(h,r):
  return h * r

#General equation for big_h 
def big_h(big_H_nought, density_profile, amplitude, polytropic_index, r, r_nought):
  return big_h_nought * (density_profile/(density_profile_nought*amplitude))**(1/(2n+1)) * (r/r_nought)**(3n/(2n+1))

#Density at r_nought, don't remember what k is but I think it's an arbitrary constant that we used somewhere else
def rho_nought(g, mass_star, big_h, k, polyatropic_index)
  return (g*mass_star*big_h**2)/(2k*(1+n)*r^3)

#Density equation
def rho(rho_nought, z, big_h):
  return rho_nought*(1-(z**2/big_h**2))
  
#With these equations instead of the old equations I think the angular momentum should come out correctly because it will use the new
#value for big_h, but I wanted to just leave them commented out until we actually talked about it.
'''

# TODO Add Length Scale
def gaussian(x, mu, sig):
    '''
    Return the gaussian distribution
    :rtype : float
    :param x: Generally, r, the radius
    :param mu: Generall r_nought, the center of the bump
    :param sig: delta_r, the width of the bump
    :return: The gaussian distribution, not scaled for the amplitude
    '''
    distribution = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.0)))
    scaled = (1.0 / (math.sqrt(2 * math.pi) * sig)) * distribution
    return scaled


def p_nought(amplitude, radius, r_nought, delta_r, alpha):
    '''
    Calculates the P_o, which gives most of the bump for the model
    :rtype : float
    :param amplitude: Amplitude of the bump
    :param radius: Radius of the point, in radius * rof3n
    :param r_nought: Center point of bump, in r_nought * rof3n
    :param delta_r: Width of bump, in delta_r * rof3n
    :param alpha: coefficient
    :return: P_o
    '''
    gaussian_bump = gaussian(radius, r_nought, delta_r)
    # print("Gaussian Bump: " + str(gaussian_bump))
    b_r = amplitude * gaussian_bump
    # print("B(r): " + str(b_r))
    return math.pow(radius, -alpha) * b_r


def big_h(h, r):
    '''
    Calculates the H(r) value for the model
    :rtype : float
    :param h: little h value
    :param r: radius to calculate for
    :return: H(r) value for the model
    '''
    # print("H: " + str(h * r))
    return h * r


def p_nought_coefficient(h, z, r):
    '''
    Calculates the coefficent for P_o in
    :rtype : float
    :param h: little h value
    :param z: height of point
    :param r: radius from center of point
    :return: Coefficient of P_o
    '''
    # print("little h: " + str(h))
    # print("Z: " + str(z))
    # print("P Nought: " + str((1.0 - (z / big_h(h, r)) ** 2)))
    return 1.0 - (z / big_h(h, r)) ** 2


def surface_density_profile(amplitude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index, jmin, rof3n,
                            zof3n):
    '''
    Calculates the surface density for the grid point in terms of the desired Gaussian density bump
    :rtype : float
    :param amplitude: Amplitude of teh bump
    :param radius: radius in grid units
    :param r_nought: center of density bump
    :param delta_r: width of density bump
    :param h: disk thickness
    :param z_height: height of grid point
    :param alpha: coefficient given in papers
    :param polytropic_index: polytropic index of model
    :param jmin: min distance in which to calculate the angular momentum
    :param rof3n: length scale for r grid points
    :param zof3n: length scale for z grid points
    :return: Value for density at the given point
    '''
    z = z_height / big_h(h, radius * rof3n)
    if radius > jmin and z_height <= big_h(h, radius):
        density_at_point = p_nought(amplitude, radius * rof3n, r_nought * rof3n, delta_r * rof3n, alpha=alpha) \
                           * (p_nought_coefficient(h, z * zof3n, radius * rof3n)) ** polytropic_index
    else:
        density_at_point = 0.000
    return density_at_point


# Angular Momentum functions
def pressure(density, constant, polytropic_index):
    """
    Calculate pressure
    :rtype : float
    :param density: density at point
    :param constant: coefficient
    :param polytropic_index: polytropic index of model
    :return: Pressure at that density
    """
    power = 1 + 1 / polytropic_index
    return constant * density ** power


def phi(radius, z, mass_star, g):
    """
    Calculate phi
    :rtype : float
    :param radius: radius, in radius * length scale
    :param z: z height, in z * zof3n
    :param mass_star: mass of star
    :param g: graviational constant
    :return: Phi
    """
    coefficient = 1 - (z ** 2 / (2 * radius ** 2))
    gravity_and_mass = (-g * mass_star) / radius
    return gravity_and_mass * coefficient


def gradient_phi_r(radius, z, mass_star, g):
    """
    Calculate gradient of phi with respect to radius
    :rtype : float
    :param radius: Radius, in terms of radius * rof3n
    :param z: z, in terms of z * zof3n
    :param mass_star: mass of star
    :param g: gravitational constant
    :return: Gradient Phi
    """
    top = g * mass_star * (2 * radius ** 2 - 3 * z ** 2)
    bottom = 2 * radius ** 4
    return top / bottom


def gradient_phi_z(radius, z, mass_star, g):
    """
    Calculate gradient of phi with respect to height
    :rtype : float
    :param radius: Radius, in terms of radius * rof3n
    :param z: z, in terms of z * zof3n
    :param mass_star: mass of star
    :param g: gravitational constant
    :return: Gradient Phi
    """
    top = g * mass_star * z
    bottom = radius ** 3
    return top / bottom


def gradient_pressure_z(ampltiude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index):
    """
    Calculate the gradient of pressure in height
    :rtype : float
    :param ampltiude: Amplitude of bump
    :param radius: radius, in terms of radius * rof3n
    :param r_nought: center of bump, in terms of r_nought * rof3n
    :param delta_r: width of bump, in terms of delta_r * rof3n
    :param h: little h
    :param z_height: height of point, in terms of z * zof3n
    :param alpha: coefficient
    :param polytropic_index: polytropic index of model
    :return: Gradient in z direction of pressure
    """
    e_term = math.exp(-(radius - r_nought) ** 2 / (2 * delta_r ** 2))
    alpha_term = radius ** (-alpha)
    ampltiude_term = -(2 * polytropic_index * z_height * ampltiude) / ((radius * h) ** 2)
    polytropic_term = (1 - (z_height ** 2) / ((radius * h) ** 2)) ** (polytropic_index - 1)
    return e_term * alpha_term * ampltiude_term * polytropic_term


def gradient_pressure_r(ampltiude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index):
    """
    Calculate the gradient of pressure in radius
    :rtype : float
    :param ampltiude: Amplitude of bump
    :param radius: radius, in terms of radius * rof3n
    :param r_nought: center of bump, in terms of r_nought * rof3n
    :param delta_r: width of bump, in terms of delta_r * rof3n
    :param h: little h
    :param z_height: height of point, in terms of z * zof3n
    :param alpha: coefficient
    :param polytropic_index: polytropic index of model
    :return: Gradient in r direction of pressure
    """
    exponential_term = -(radius - r_nought) ** 2 / ((2 * delta_r ** 2))
    inside_polytropic = (1 - (z_height ** 2) / ((radius * h) ** 2))
    over_h_term = 2 * ampltiude * polytropic_index * z_height ** 2 * radius ** (-alpha - 3) * math.exp(
        exponential_term) * \
                  inside_polytropic ** (polytropic_index - 1)
    final_h_term = over_h_term / (h ** 2)
    no_denom_term = alpha * ampltiude * radius ** (-alpha - 1) * math.exp(
        exponential_term) * inside_polytropic ** polytropic_index
    over_delta_r_term = ampltiude * radius ** (-alpha) * (radius - r_nought) * math.exp(
        exponential_term) * inside_polytropic \
                            ** polytropic_index
    final_delta_r_term = over_delta_r_term / (delta_r ** 2)
    final_term = final_h_term - no_denom_term - final_delta_r_term
    return final_term


def velocity_field(ampltiude, radius, r_nought, delta_r, h, z_height, alpha, polytropic_index, density, mass_star, g):
    """
    Calculate velocity field
    :param ampltiude: Amplitude of density bump
    :param radius: radius in terms of radius * rof3n
    :param r_nought: center of bump, in r_nought * rof3n
    :param delta_r: width of bump, in delta_r * rof3n
    :param h: little h
    :param z_height: height, in terms of z * zof3n
    :param alpha: coefficient
    :param polytropic_index: polytopic index of model
    :param density: density at the point
    :param mass_star: mass of the star
    :param g: graviational constant
    :return: Value for the velocity field at that point
    """
    pressure_gradient_grid = []
    gravitational_gradient_grid = []
    pressure_gradient_grid.append(
        (-(1 / density) * gradient_pressure_z(ampltiude, radius, r_nought, delta_r, h, z_height,
                                              alpha, polytropic_index),
         -(1 / density) * gradient_pressure_r(ampltiude, radius, r_nought, delta_r, h, z_height, alpha,
                                              polytropic_index)))
    gravitational_gradient_grid.append((-gradient_phi_r(radius, mass_star=mass_star, z=z_height, g=g),
                                        -gradient_phi_z(radius, mass_star=mass_star, z=z_height, g=g)))


def angular_velocity(radius, g, mass_star):
    """
    Return angular velocity at a given radius, assumed that velocity is only along cylindrical coordinates
    :rtype : float
    :param radius: radius, in terms of radius * rof3n
    :param g: graviational constant
    :param mass_star: mass of star
    :return: Angular momentum of radius
    """
    l = math.sqrt((g * mass_star) / (radius) ** 3)
    return l


def unit_mass(rof3n, zof3n, density):
    """
    Calculate unit mass for the grid
    :rtype : float
    :param rof3n: length scale of radius
    :param zof3n: length scale of height
    :param density: density at a point
    :return: Unit mass
    """
    mass = rof3n * zof3n * density
    return mass


def velocity_squared(g, mass_star, h, radius):
    """
    Calculate v^2 for a given radius
    :rtype : float
    :param g: graviational constant
    :param mass_star: mass of the star
    :param h: little h
    :param radius: radius, in terms of radius * rof3n
    :return: v^2 for a point
    """
    ''' v^2 = (GM/r)*(1 - 3H^2/2r^2 + H/r * dH/dr)'''
    star_influence = (g * mass_star) / (radius)
    second_half = 1 - (3 * (h * radius) ** 2) / (2 * (radius) ** 2) + (h * radius) / (
        radius) * h
    return star_influence * second_half


def angular_velocity_1(radius, rof3n, z, density, g, mass_star):
    """
    Calculate angular velocity
    :rtype : float
    :param radius:
    :param rof3n:
    :param z: height, in terms of z * zof3n
    :param density: Density at point
    :param g: graviational constant
    :param mass_star: mass of star
    :return: Angular velocty for a given point
    """
    #TODO Figure out what this does
    velocity = radius * rof3n
    omega = math.sqrt((g * mass_star) / (radius * rof3n) ** 2)
    return velocity * omega


def angular_momentum(radius, rof3n, z, zof3n, density, g, mass_star, h):
    '''
    Calculates the angular momentum for the grid point (r, z) from the density from surface_density_profile function
    :param radius: Radius (in grid units)
    :param rof3n: Length scale factor for grid in r direction
    :param z: height (in grid units)
    :param zof3n: length scale factor for grid in z direction
    :param density: Density at the grid point (r, z)
    :param g: gravitational constant
    :param mass_star: mass of the star
    :param jmin: minimum radius in which to compute the angular momentum
    :return: The angular momentum for that grid point
    '''
    return (radius * rof3n) * unit_mass(rof3n=rof3n, zof3n=rof3n, density=density) * math.sqrt(
        velocity_squared(g=g, mass_star=mass_star, h=h, radius=radius * rof3n))
    # return angular_velocity_1(radius, rof3n, z, zof3n, density, g, mass_star) * (radius * rof3n) * unit_mass(radius,
    #                                                                                                        rof3n,
    #                                                                                                       z,
    #                                                                                                      zof3n,
    #                                                                                                     density)


def generate_fort_2(polytropic_index, jmax, kmax, jout, kout, log_central_density, iteration, mass_star, xcut,
                    xwidth, xnorm, type):
    """
    Generate fort.2 model file for CHYMERA Code
    :param polytropic_index: Polytropic index of star
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
    jmax2 = jmax + 2
    jmax1 = jmax + 1
    if type == 'RWI':
        print('Using Rossby Wave equations for disk')
        denny = [[0 for x in range(jmax2)] for x in range(jmax2)]
        anggy = [[0 for x in range(jmax1)] for x in range(jmax1)]

        constants_array = []
        header_line = ""
        with open("fort.2", "r") as example_file:
            header_line = example_file.readline()
            print(header_line)
        constants_array = header_line.split(" ")
        # Remove the first empty parts that come from the '3X' Fortran Formatting
        del constants_array[0]
        del constants_array[0]
        del constants_array[0]
        del constants_array[0]
        print(constants_array)
        # Convert to floats
        for element in range(len(constants_array)):
            constants_array[element] = float(constants_array[element])
        print(constants_array)
        '''
        Format of first line in fort.2 file
        PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,ROF3N,ZOF3N,
&        A1NEWZ,JREQ,KZPOL
        '''

        # Get the density array
        for column in range(jmax2):
            for row in range(jmax2):
                denny[row][column] = surface_density_profile(amplitude, row + 1, r_nought, delta_r, h, column + 1, alpha,
                                                             constants_array[0], jout, constants_array[6],
                                                             constants_array[7])

        # Get angular momentum array
        for column in range(jmax1):
            for row in range(jmax1):
                anggy[row][column] = angular_momentum(row + 1, constants_array[6], column + 1, constants_array[7],
                                                      denny[row][column], g, mass_star, h=h)
        print("Length of Anggy: " + str(len(anggy)))

        with open("temp", 'w') as model_file:
            model_file.write(header_line)
            fortran_writer = FortranRecordWriter('8(1PE10.3,2X)')
            # Write header line
            # Fortran saves out arrays column first, so first row in file would be the first entry in each row in array
            # Each line is composed of 8 floating point with 22 spaces with 15 after the decimal place, then two spaces
            # Writing density array to fort.2
            temp_denny = []
            for column in range(jmax2):
                for row in range(jmax2):
                    temp_denny.append(denny[column][row])
                    ''' if len(temp_denny) == 8:
                        writer.writerow(temp_denny)
                        temp_denny = []
                    if row == jmax1 and column == jmax1:
                        writer.writerow(temp_denny)
                        temp_denny = []'''
            output_text = fortran_writer.write(temp_denny)
            output_text += "\n"
            model_file.write(output_text)
            # Writing the angular momentum array to fort.2
            # Repeat what was done for the denny array
            temp_anggy = []
            for column in range(jmax1):
                for row in range(jmax1):
                    temp_anggy.append(anggy[column][row])
            print("Length of Temp Anggy: " + str(len(temp_anggy)))
            output_text = fortran_writer.write(temp_anggy)
            output_text += "\n"
            model_file.write(output_text)
            # TODO Input coordinates of points into RWI equations and output them with 8 density points per line
            # TODO Then for specific angular momentum, same thing
    else:
        print('Using normal equations for disk')


generate_fort_2(polytropic_index=polytropic_index, jmax=r_size, kmax=z_size, jout=jout, kout=kout,
                log_central_density=lcd, iteration=iteration, mass_star=mass_star, xcut=xcut, xnorm=xnorm,
                xwidth=xwidth, type='RWI')
