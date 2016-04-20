# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 00:28:47 2016

@author: robouser
"""

import mulens_class
from astropy.time import Time, TimeDelta
import numpy as np
from astropy import constants
import matplotlib.pyplot as plt

def generate_mag_err_relations():
    """Function to generate datafiles of the magnitude error relations"""
    
    event = mulens_class.MicrolensingEvent()
    event.u_min = 0.0
    event.u_offset = 0.001
    event.t_E = TimeDelta((1.0 * 24.0 * 3600.0),format='sec')
    event.phi = ( 0.0 * np.pi ) / 180.0
    event.rho = 0.001
    event.M_L = constants.M_sun * 0.3
    event.D_L = constants.pc * 3000.0
    event.D_S = constants.pc * 8000.0
    event.RA = '17:57:34.0'
    event.Dec = '-29:13:15.0'
    event.t_o = Time('2015-06-15T15:00:00', format='isot', scale='utc')
    event.get_earth_perihelion()
    
    exp_time = 200.0
    
    output = open('/home/robouser/mag_err_relation.data','w')
    output.write('# Column 1: magnitude\n')
    output.write('# Column 2: magnitude uncertainty (1m telescope on Earth)\n')
    output.write('# Column 3: magnitude uncertainty (Swift)\n')
    mags = []
    earth_merr = []
    swift_merr = []
    for mag in np.arange(12,18.0,0.01):
        event.mag_base = mag
        mags.append( mag )
        earth_merr.append( event.sim_mag_error( exp_time, mag, precision_model='1m') )
        swift_merr.append( event.sim_mag_error( exp_time, mag, precision_model='swift') )
        output.write( str(mag) + ' ' + str(earth_merr[-1]) + ' ' + str(swift_merr[-1]) + '\n' )
    output.close()
    
    mags = np.array(mags)
    earth_merr = np.array( earth_merr )
    swift_merr = np.array( swift_merr )
    
    fig = plt.figure(1,(12,12)) 
    plt.plot( mags, earth_merr, 'r.', label='Earth 1m')
    plt.plot( mags, swift_merr, 'b+', label='Swift')
    plt.xlabel( 'Mag' )
    plt.ylabel( 'Mag uncertainty' )
    plt.yscale('log')
    plt.legend(loc='best',frameon=False)
    plt.grid()
    plt.savefig('/home/robouser/mag_err_relation.png')
    
if __name__ == '__main__':
    generate_mag_err_relations()
    