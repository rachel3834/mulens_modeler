# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 05:46:14 2016

@author: robouser
"""
import mulens_class
import numpy as np

def calc_umin_grid():
    """Function to determine the umin grid necessary"""
    
    # Example lensing event
    event = mulens_class.MicrolensingEvent()
    event.t_E = TimeDelta((1.0 * 24.0 * 3600.0),format='sec')
    event.phi = ( 0.0 * np.pi ) / 180.0
    event.mag_base = 12.0
    event.rho = 0.001
    event.M_L = constants.M_sun * 0.3
    event.D_L = constants.pc * 4000.0
    event.D_S = constants.pc * 8000.0
    event.RA = '17:57:34.0'
    event.Dec = '-29:13:15.0'
    event.t_o = Time('2015-06-15T16:00:00', format='isot', scale='utc')
    event.t_p = Time('2015-06-15T06:37:00', format='isot', scale='utc')
    
    fileobj = open('umin_grid.dat', 'w')
    for umin in np.arange(0.00001, 0.2, 0.001):
        event.u_min = umin
        
        event.calc_D_lens_source()
        event.calc_einstein_radius()
        event.gen_event_timeline()
        event.calc_source_lens_rel_motion()
        event.calc_proj_observer_pos(parallax=True,satellite=False)
        event.calc_parallax_impact_param()
        
        fileobj.write('u_min = ' + str(event.u_min) + \
                            ' -> u_o = ' + str(event.u_o) + '\n' )
        fileobj.flush()
    fileobj.close()
    
if __name__ == '__main__':
    calc_umin_grid()