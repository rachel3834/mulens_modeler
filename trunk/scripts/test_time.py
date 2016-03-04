# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 00:09:24 2016

@author: robouser
"""
import mulens_class
from astropy.time import TimeDelta
from astropy import constants
import numpy as np

def test_time():
    """Function to test the timeline calculations"""
    
    # Example lensing event
    event = mulens_class.MicrolensingEvent()
    event.u_o = 0.1
    event.t_E = TimeDelta((50.0 * 24.0 * 3600.0),format='sec')
    event.phi = ( 0.0 * np.pi ) / 180.0
    event.mag_base = 15.0
    event.rho = 0.005
    event.M_L = constants.M_sun * 0.3
    event.D_L = constants.pc * 4000.0
    event.D_S = constants.pc * 8000.0
    event.RA = '17:57:34.0'
    event.Dec = '-29:13:15.0'
    event.t_o = Time('2015-01-04T16:00:00', format='isot', scale='utc')
    event.t_p = Time('2015-01-04T06:37:00', format='isot', scale='utc')
    
    event.gen_event_timeline(debug=True)
    
if __name__ == '__main__':
    test_time()
    