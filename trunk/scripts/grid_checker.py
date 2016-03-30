# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 03:17:37 2016

@author: robouser
"""

import glob
import numpy as np
import swift_grid_simulator
from sys import argv, exit
from astropy import constants
from os import path

##############################################################################
#               SWIFT GRID CHECKER
##############################################################################

def grid_checker( params ):

    uvalues = np.array( params['uoffset_list'] )
    nu = len(params['uoffset_list'])
    
    grid = construct_grid( params )
    n_grid = len(grid) * nu
    
    n_output = len( glob.glob( path.join( params['data_dir'], '*_earth.dat' ) ) )
    
    print 'Expecting ' + str(n_grid) + ' grid points, found output for ' + str(n_output)
    print 'Missing ' + str( n_grid - n_output ) + ' files'    
    
    fileobj = open(path.join(params['data_dir'],'grid_locale.dat'),'w')
    fileobj.write('# u_0  t_E  phi  mag  rho  M_L  D_L  file_name\n')
    for g, grid_point in enumerate(grid):
        
        root_name = file_search_string( grid_point )
        par_string = root_name.replace('lc_*_','').replace('_earth.dat','').replace('_',' ')
        
        file_list = glob.glob( path.join(params['data_dir'],root_name) )
        file_list.sort()
        
        if len(file_list) == nu:
            for i,u in enumerate(params['uoffset_list']):
                fileobj.write(str(u) + ' ' + par_string + ' ' + \
                        file_list[i] + '\n')
        else:
            for i,u in enumerate(params['uoffset_list']):
                fileobj.write(str(u) + ' ' + par_string + ' ' + \
                        ' ERROR: Missing output\n')
            
    fileobj.close()
  
def file_search_string( grid_point ):
    """Method to return a string summarizing the basic FSPL parameters
    for use as a root file name"""
    
    te = grid_point[2]
    te = round( te, 1 )
    phi = round( grid_point[3], 3 )
    mag = round( grid_point[4], 1 )
    ml = round( grid_point[0], 6 )
    dl = round( grid_point[1], 4 )
    rho = grid_point[5]
    
    file_name = 'lc_*_' + str(te) + \
                    '_' + str(phi) + '_' + str(mag) +\
                    '_' + str(rho) + '_' + str(ml) +\
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
                            grid.append( [ml,dl,te,phi,Vbase,rho] )
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
        params['data_dir'] = path.dirname(file_path)
        grid_checker( params )