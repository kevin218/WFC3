import numpy as np

#Constants
k     = 1.38065e-23                     #Boltzman's constant     (kg m^2 s^-2 K^-1)
g     = 6.67384e-11                     #gravitational constant  (m^3 kg^-1 s^-2)
#mu    = 2.2 * 1.0078e0 * 1.6605e-27     #mean molecular weight assuming H2/He   (kg)
r_jup = 7.1492e7                        #Radius of Jupiter  (m)
m_jup = 1.89813e27                      #Mass of Jupiter    (kg)
r_sun = 6.955e8                         #Radius of Sun      (m)

def calc_gravtiy(scale_height, temperature, mu=2.2):
    '''
    
    INPUTS
    ------
    scale_height    : Scale Height
    temperature     : Planet temperature in K
    
    RETURNS
    -------
    surface_gravity : Gravity
    '''
    mu *= 1.0078e0 * 1.6605e-27
    
    #scale_height    = k * temperature / surface_gravity / mu / 1000    #km
    surface_gravity   = k * temperature / scale_height / mu / 1000    #km
    
    return surface_gravity

def calc_scaleheight(mass, radius, temperature, mu=2.2, returnGravity=False):
    '''
    
    INPUTS
    ------
    mass            : Planet mass in M_jup
    radius          : Planet radius in R_jup
    temperature     : Planet temperature in K
    
    RETURNS
    -------
    scale_height    : Planet scale height in km
    '''
    mu *= 1.0078e0 * 1.6605e-27
    
    surface_gravity = g * (mass*m_jup) / (radius*r_jup)**2
    scale_height    = k * temperature / surface_gravity / mu / 1000    #km
    
    if returnGravity == False:
        return scale_height
    else:
        return scale_height, surface_gravity

def depth2scaleheight(depth, zeropt, Rs, mass, radius, temperature, mu=2.2):
    '''
    
    INPUTS
    ------
    depth           : Transit depths in %
    zeropt          : Transit depth zero reference point in %
    Rs              : Stellar radius in units of r_sun
    mass            : Planet mass in M_jup
    radius          : Planet radius in R_jup
    temperature     : Planet temperature in K
    
    RETURNS
    -------
    H               : Transits depth in units of scale heights
    '''
    
    #rprs  = np.sqrt(zeropt/100.)
    #H     = 0.5*(depth-zeropt)/100.*Rs*r_sun/rprs/1000   #km
    H     = 0.5*(depth-zeropt)/100.*(Rs*r_sun)**2/radius/r_jup/1000   #km
    H    /= calc_scaleheight(mass, radius, temperature, mu) #unitless
    
    return H
    
def scaleheight2depth(H, Rs, mass, radius, temperature, mu=2.2):
    '''
    
    INPUTS
    ------
    H               : Number of scale heights
    Rs              : Stellar radius in units of r_sun
    mass            : Planet mass in M_jup
    radius          : Planet radius in R_jup
    temperature     : Planet temperature in K
    
    RETURNS
    -------
    ddepth          : Variation in transit depth in %
    '''
    H       = np.copy(H)
    H      *= calc_scaleheight(mass, radius, temperature, mu)*1000  #m
    ddepth  = 2*H*radius*r_jup/(Rs*r_sun)**2
    
    return ddepth*100
