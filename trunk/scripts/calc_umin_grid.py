# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 05:46:14 2016

@author: robouser
"""
import mulens_class
import numpy as np
from astropy import constants, coordinates, units
from astropy.time import Time, TimeDelta
from swift_grid_simulator import parse_input_file

def calc_umin_grid(grid_dir):
    """Function to determine the umin grid necessary"""

    # Get the parameters of this grid-point:
    par_file = path.join( grid_dir, 'swift_grid.par' )
    params = parse_input_file( par_file )
    
    # Calculate the full range of parameter space covered by this 
    # grid point definition, EXCEPT umin:
    grid = construct_grid( params )
    
    # Loop over the full parameter space, recording those umin values 
    # which result in low u_o values:    
    fileobj = open(path.join(grid_dir,'umin_grid.dat'), 'w')
    for grid_point in grid:
        for offset in [ 0.0, 0.2 ]:
            
            print 'OFFSET = ',offset
            # Example lensing event
            event = mulens_class.MicrolensingEvent()
            event.u_min = 0.0
            event.u_offset= offset
            event.t_E = TimeDelta((grid_point[2] * 24.0 * 3600.0),format='sec')
            event.phi = ( grid_point[3] * np.pi ) / 180.0
            event.mag_base = 12.0
            event.rho = 0.001
            event.M_L = constants.M_sun * grid_point[0]
            event.D_L = constants.pc * grid_point[1]
            
            event.D_S = constants.pc * 8000.0
            event.RA = '17:57:34.0'
            event.Dec = '-29:13:15.0'
            event.t_o = Time('2015-06-15T16:00:00', format='isot', scale='utc')
            
            event.get_earth_perihelion()
        
            
            event.calc_D_lens_source()
            event.calc_einstein_radius()
            event.gen_event_timeline(lc_length=grid_point[2])
            event.calc_source_lens_rel_motion()
            event.calc_proj_observer_pos(parallax=True,satellite=False)
            event.calc_parallax_impact_param(set_uo=True)
            A = event.calc_pspl_A( event.u_o )
            
            #if event.u_o < 0.01:
            fileobj.write( event.summary(inc_uo=True) + '\n' )
            fileobj.flush()
            
            event.plot_lens_plane_motion()
            
            cont = False
            if cont == True:
                event.get_earth_perihelion()
                event.calc_D_lens_source()
                event.calc_einstein_radius()
                event.gen_event_timeline(lc_length=0.01)
                event.calc_source_lens_rel_motion()
                event.calc_proj_observer_pos(parallax=True,satellite=False)
                
            
    fileobj.close()
    
def construct_grid( params ):
    """Function to return a list of gridpoints.  Each list entry consists
    of a list of grid parameters:
    [ML, DL, tE, phi]
    """
        
    (temin, temax, teincr) = params['te_range']
    (phimin, phimax, phiincr) = params['phi_range']
    (rhomin, rhomax, rhoincr)= params['rho_range']
    (dlmin, dlmax, dlincr)= params['lens_distance_range']
    
    grid = []
    for ml in params['lens_mass_list']:
        for dl in np.arange( dlmin, dlmax, dlincr ):
            for te in np.arange( temin, temax, teincr ):
                for phi in np.arange( phimin, phimax, phiincr ):
                    grid.append( [ml,dl,te,phi] )
    return grid
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        grid_dir = raw_input( 'Please enter the path to the grid directory: ' )
    else:
        grid_dir = argv[1]
    calc_umin_grid( grid_dir )