# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 19:28:29 2016

@author: robouser
"""

import mulens_class
from astropy.time import TimeDelta
from astropy import constants
import numpy as np

def test_parallax():
    """Function to test the microlensing parallax calculations"""
    
    # Example lensing event
    event = mulens_class.MicrolensingEvent()
    event.u_o = 0.1
    event.t_E = TimeDelta((100.0 * 24.0 * 3600.0),format='sec')
    event.phi = ( 0.0 * np.pi ) / 180.0
    event.mag_base = 15.0
    event.rho = 0.02
    event.M_L = constants.M_sun * 0.3
    event.D_L = constants.pc * 6000.0
    event.D_S = constants.pc * 8000.0
    event.RA = '17:57:34.0'
    event.Dec = '-29:13:15.0'
    event.t_o = Time('2015-01-04T16:00:00', format='isot', scale='utc')
    event.t_p = Time('2015-01-04T06:37:00', format='isot', scale='utc')
    
    # Calculate PSPL lightcurve
    event.calc_D_lens_source()
    event.calc_einstein_radius()
    event.gen_event_timeline()
    event.calc_source_lens_rel_motion()
    event.calc_proj_observer_pos(parallax=False,satellite=False)
    event.calc_pspl_impact_param()
    event.calc_magnification(model='fspl')
    
    # Calculate PSPL+parallax lightcurve
    event.calc_proj_observer_pos(parallax=True,satellite=False)
    event.calc_pspl_parallax_impact_param()
    event.calc_magnification(model='fspl_parallax')
    
    event.plot_diff_models('parallax_test.png', \
                                model_list = ['fspl', 'fspl_parallax'])
    
    
if __name__ == '__main__':
    test_parallax()
    