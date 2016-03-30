# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 03:17:37 2016

@author: robouser
"""

import glob
import numpy as np
import swift_grid_simulator

##############################################################################
#               SWIFT GRID CHECKER
##############################################################################

def grid_checker( params ):

    grid = construct_grid( params )
    n_grid = str(len(grid))

    uvalues = np.array( params['uofflist_list'] )
    
    fileobj = open('grid_locale.dat','w')
    
    for g, grid_point in enumerate(grid):
        
        root_name = file_search_string( grid_point )
        
        file_list = glob.glob( root_name )
        file_list.sort()
        
        if len(file_list) == 4:
            for i,u in enumerate(params['uofflist_list']):
                fileobj.write(str(u) + ' ' + file_list[i]) + '\n'
        else:
            for i,u in enumerate(params['uofflist_list']):
                fileobj.write(str(u) + ' ERROR: Missing output\n'
            
    fileobj.close()
  
def file_search_string( grid_point ):
    """Method to return a string summarizing the basic FSPL parameters
    for use as a root file name"""
    
    te = grid_point[3] / ( 60.0 * 60.0 * 24.0 )   
    te = round( te, 1 )
    phi = round( grid_point[4], 3 )
    mag = round( grid_point[5], 1 )
    ml = round( grid_point[0]/constants.M_sun.value, 6 )
    dl = round( grid_point[1]/constants.pc.value, 4 )
    #except TypeError:
    #    uo = round( self.u_min, 4 )
    file_name = 'lc_*_' + str(te) + \
                    '_' + str(phi) + '_' + str(mag) +\
                    '_' + str(self.rho) + '_' + str(ml) +\
                    '_' + str(dl) + '_earth.dat'
    return file_name

def construct_grid( params ):
    """Function to return a list of gridpoints.  Each list entry consists
    of a list of grid parameters:
    [u0, tE, phi, Vbas, rho]
    """

    if 'uoffset_range' in params.keys():    
        (um_min, um_max,um_incr) = params['uoffset_range']
        uoffset_list = np.arange( um_min, um_max, um_incr )
    else:
        uoffset_list = params['uoffset_list']
        
    (temin, temax, teincr) = params['te_range']
    (phimin, phimax, phiincr) = params['phi_range']
    (vmin, vmax, vincr) = params['v_range']
    (rhomin, rhomax, rhoincr)= params['rho_range']    
    (dlmin, dlmax, dlincr)= params['lens_distance_range']
    
    grid = []
    for ml in params['lens_mass_list']:
        for dl in np.arange( dlmin, dlmax, dlincr ):
            for te in np.arange( temin, temax, teincr ):
                for phi in np.arange( phimin, phimax, phiincr ):
                    for Vbase in np.arange( vmin, vmax, vincr ):
                        for rho in np.arange( rhomin, rhomax, rhoincr ):
                            grid.append( [ml,dl,um,te,phi,Vbase,rho] )
    return grid

#############################################
# COMMANDLINE RUN SECTION
if __name__ == '__main__':
    
    if len(argv) == 1:
        print """Call sequence:
                python grid_checker.py [parameter_file]
              """
    else:
        file_path = argv[1]
        params = swift_grid_simulator.parse_input_file( file_path )
        grid_checker( params )