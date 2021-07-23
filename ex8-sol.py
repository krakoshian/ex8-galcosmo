

import scipy.integrate as integrate
import scipy.optimize as root
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
plt.rc('text', usetex=True)
plt.rc('font',family = 'serif', size=20)

#*************************************************************************
# TO DO:
# - Fill in variables for Omega_m0, Omega_X0, h (line 26)
# - Change redshift bounds to those given in question (line 37)
# - Add the four models (line 42)
# - Some integrals are over 1/E(z)dz, but others have extra terms. 
#   You will need to write a new function to integrate over for this case.
#   1/E(z)dz is already given for you.
# - Add in  equations for:
#    a) Age of the Universe
#    b) Luminosity distance
#    c) Co-moving volume
#    d) Alcock-Paczynski test (note this is inverted compared to Ex6!)
#   (starting line 68)
#*************************************************************************

#FILL IN VARIABLES
Omega_m0 = 0.3
Omega_X0 = 0.7
h = 0.7

#Hubble constant at present time in Gyr
H_0 = (h * 100. * u.km /u.s /u.Mpc).to(1/u.Gyr).value
#Speed of light in Gpc/Gyr
c = c.c.to(u.Gpc /u.Gyr).value

#Array of redshift values
lower_bound = 0.
upper_bound = 5.
redshift = np.linspace(lower_bound, upper_bound,100)

#FILL IN CASES
W = [
    lambda z: -1.,
    lambda z: -1./3.,
    lambda z: -0.5 + 0.1*z,
    lambda z: -0.5 - 0.05*z
    ]

def IntGeneric(function, lower, upper, parameters):
    """Integration function"""
    result, error =  integrate.quad(function, lower, upper, args=parameters)
    return result

def f(z,w):
    """Bracketed section of Eq.1 of the exercise sheet"""
    return (1. + w(z)) / (1. + z)

def E_inverse(z, Omega_m0, Omega_X0, w):
    """Returns E^-0.5 from Eq.1 of the exercise sheet, where H = H_0 * E """
    E_squared = Omega_m0 * (1. + z)**3 + Omega_X0 * np.exp(3. * IntGeneric(f, 0., z, w))
    
    return E_squared**-0.5
    
def E_inverse_2(z, Omega_m0, Omega_X0, w):
    """Returns E^-0.5/(1+z) from Eq.1 of the exercise sheet, where H = H_0 * E """
    E_squared = Omega_m0 * (1. + z)**3 + Omega_X0 * np.exp(3. * IntGeneric(f, 0., z, w))
    
    return E_squared**-0.5 / (1. + z)
    
#************
# FILL IN
#************

def Age(z, arguments):
    """Should calculate age of the Universe"""
    return 1. / H_0 * IntGeneric(E_inverse_2, z, np.inf, arguments)

def D_lum(z, arguments):
    """Should calculate luminosity distance"""
    return c / H_0 *(1.+z) * IntGeneric(E_inverse, 0., z, arguments)

def V_com(z, arguments):
    """Should calculate co-moving volume"""
    return 4. * np.pi / 3. * (c /H_0 * IntGeneric(E_inverse, 0., z, arguments))**3.
    
def AP_test(z, arguments):
    """Should calculate delta z/ delta theta"""
    Omega_m0, Omega_X0, w = arguments
    return 1/(E_inverse(z, Omega_m0, Omega_X0, w)) *IntGeneric(E_inverse, 0., z, arguments)

#*************************************
# The below should not need changing  
#*************************************

#Labels for plotting    
line_label = [
        r'$w=-1$',
        r'$w=-1/3$',
        r'$w=-0.5+0.1z$',
        r'$w=-0.5-0.05z$',
        ]
        
y_label = [
        'Time [Gyr]',
        'Luminosity distance [Gpc]',
        'Co-moving volume [Gpc$^3$]',
        r'$\Delta z / \Delta \theta$'
        ]
        
title = [
        'Age of the Universe',
        'Luminosity distance',
        'Co-moving volume',
        'Alcock-Paczynski test'
        ]
    
def plot(x, y, y_label, title):
    """Plotting function"""

    #Initialise plot
    fig, ax = plt.subplots(1, 1)
    
    #Colours for lines
    colors = ['#008141','#fdc513','#8b0a50','#4f7687']
    
    for y, c, l in zip(y, colors, line_label):
            ax.plot(x, y, c=c, label=l)
    
    #Add legend
    plt.legend()
    
    #Set formatting (bounds, title, etc)
    ax.set_xlabel(r'$z$')
    ax.set_xlim([lower_bound, upper_bound])
    ax.set_ylabel(y_label)
    ax.set_title(title)

    fig.set_size_inches(8, 6)
    fig.tight_layout()
    
    fig.savefig('{}.pdf'.format(title))
    print('saved')
    #plt.show()
    

#Collects data for each model over required redshift range

#Lists for each measurement
A = []
B = []
C = []
D = []

#Loop over models
for w in W:
    AA = []
    BB = []
    CC = []
    DD = []
    
    #Loop over redshifts
    for z in redshift:
        arguments = (Omega_m0, Omega_X0, w)

        AAA = Age(z, arguments)
        BBB = D_lum(z, arguments)
        CCC = V_com(z, arguments)
        DDD = AP_test(z, arguments)
        
        AA.append(AAA)
        BB.append(BBB)
        CC.append(CCC)
        DD.append(DDD)
        
    A.append(AA)
    B.append(BB)
    C.append(CC)
    D.append(DD)
    
#Plot each figure and save it    
plot(redshift, A, y_label[0], title[0])    
plot(redshift, B, y_label[1], title[1])
plot(redshift, C, y_label[2], title[2])
plot(redshift, D, y_label[3], title[3])

