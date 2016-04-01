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
    print 'Expecting ' + str(n_grid) + ' grid points, found output for ' + \
                            str(n_output)
    print 'Missing ' + str( n_grid - n_output ) + ' files'    
    
    # If available, read back the statistics file to include it in the output. 
    # Returns an empty dictionary if none available. 
    stats_data = read_grid_stats( params )
    
    fileobj = open(path.join(params['data_dir'],'grid_locale.dat'),'w')
    fileobj.write('# u_0  t_E  phi  mag  rho  M_L  D_L  dchisq  chisq_swift  chisq_swift_earth  lc_file_root  stats_file_prefix\n')
    for g, grid_point in enumerate(grid):
        
        # Find all the output files with these parameters, but include the
        # complete set of simulated u0 values together:
        root_name = file_search_string( grid_point )
        par_string = root_name.replace('lc_*_','').replace('_earth.dat','').replace('_',' ')
        
        # Ensure the file list is sorted into alphanumerical order:
        file_list = glob.glob( path.join(params['data_dir'],root_name) )
        file_list.sort()
        
        if len(file_list) == nu:
            for i,u in enumerate(params['uoffset_list']):
                prefix = path.basename(file_list[i]).replace('lc_','').replace('_earth.dat','')
                if len(stats_data) > 0:
                    (chisq_string, file_name) = stats_data[prefix]
                else:
                    chisq_string = '-9999.99999  -9999.99999  -9999.99999'
                    file_name = 'none_available'
                lc_file_root = path.basename(file_list[i]).replace('_earth.dat','')
                
                fileobj.write(str(u) + ' ' + par_string + ' ' + chisq_string + \
                        ' ' + lc_file_root + ' ' + file_name + '\n')
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
    phi = (float(grid_point[3])*np.pi) / 180.0
    phi = round( phi, 3 )
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
    [u0, tE, phi, Vbase, rho]
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

def read_grid_stats( params ):
    """Function to read back the grid statistics file"""
    
    stats_data = {}
    stats_file = path.join( params['data_dir'], 'grid_locale_stats.dat' )
    if path.isfile(stats_file) == True:
        file_lines = open(stats_file,'r').readlines()
        
        for line in file_lines:
            if line[0:1] != '#':
                entries = line.replace('\n','').split()
                par_string = '_'.join(entries[0:7])
                chisq_string = ' '.join(entries[7:10])
                file_name = entries[11]
                stats_data[par_string] = [ chisq_string, file_name ]
    
    return stats_data
    
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