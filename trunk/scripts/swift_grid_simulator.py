# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:12:44 2016

@author: rstreet
"""


################################################################################
#     	      	      	MICROLENSING EVENT SIMULATION
# 
# Program to test the implementation of Finite Source, Point Lens events
################################################################################

###########################
# IMPORT MODULES
from os import path
from sys import argv, exit
import mulens_class
import numpy as np
from astropy import constants
from astropy import coordinates
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
import copy
import glob
import log_utilities

def simulate_grid_models( params ):
    """Function to drive a simulation of a grid of microlensing models
    spanning user-defined ranges in u0, tE, phi, Vbase, rho).
    
    For each grid point, the simulation generates two lightcurves: 
    one FSPL including annual parallax
    one FSPL including annual + satellite parallax
    with datapoints which reflect the photometric precision likely from 
    a 1m telescope on Earth and the Swift satellite.
    """
    
    log = log_utilities.start_day_log( params, 'grid_sim' )
    log.info( 'Starting grid simulation' )
    
    grid = construct_grid( params )
    n_grid = str(len(grid))

    for g, grid_point in enumerate(grid):                    
        event = mulens_class.MicrolensingEvent()
        event.u_o = grid_point[0]
        event.t_E = TimeDelta((grid_point[1] * 24.0 * 3600.0),format='sec')
        event.phi = ( grid_point[2] * np.pi ) / 180.0
        event.mag_base = grid_point[3]
        event.rho = grid_point[4]
        event.M_L = constants.M_sun * params['lens_mass']
        event.D_L = constants.pc * params['lens_distance']
        event.D_S = constants.pc * params['source_distance']
        event.RA = '17:57:34.0'
        event.Dec = '-29:13:15.0'
        event.t_o = Time('2015-01-04T16:00:00', format='isot', scale='utc')
        event.t_p = Time('2015-01-04T06:37:00', format='isot', scale='utc')
        
        # Check for pre-existing output and skip if found:
        file_list = glob.glob( path.join( params['output_path'], \
                        event.root_file_name()+'*' ) )
                        
        if len(file_list) == 0:
        
            log.info( 'Computing for grid point parameters (' + str(g+1) + \
                    ' out of ' + n_grid + '):' )
            log.info( '-> ' + event.summary() )
            
            # Compute lens essential parameters
            event.calc_D_lens_source()
            log.info( '-> calculated the projected separation of lens and source' )
            event.calc_einstein_radius()
            log.info( '-> computed the Einstein radius' )
            event.gen_event_timeline(cadence=params['cadence'], \
                                        lc_length=params['lc_length'])
            log.info( '-> generated the event timeline' )
            event.calc_source_lens_rel_motion()
            log.info( '-> built lensing event object' )
            
            # For ease of handling later, a copy of the basic event
            # is taken and will be used to compute the same event
            # as seen from Swift:
            swift_event = copy.copy( event )
            log.info( '-> copied to Swift event object' )
            
            # Ground-based observer:        
            # Calculate the model lightcurve and datapoints for an FSPL 
            # event including annual parallax:
            event.calc_proj_observer_pos(parallax=True,satellite=False)
            log.info( '-> calculated the projected observer position' )
            event.calc_pspl_impact_param()
            log.info( '-> calculated the PSPL impact parameter' )
            event.calc_magnification(model='fspl')
            log.info( '-> calculated the magnification as a function of time' )
            event.simulate_data_points(model='fspl', phot_precision='1m')
            log.info( '-> Simulated ground-based model and data' )
            
            # Swift observer: 
            swift_event.swift_t = event.t[0]
            swift_event.calc_proj_observer_pos(parallax=True,satellite=True)
            log.info( '-> calculated the projected observer position' )
            swift_event.calc_pspl_impact_param()
            log.info( '-> calculated the PSPL impact parameter' )
            swift_event.calc_magnification(model='fspl')
            log.info( '-> calculated the magnification as a function of time' )
            swift_event.simulate_data_points(model='fspl', \
                                phot_precision='swift', window=0.75, interval=1.5)
            log.info( '-> Simulated Swift model and data' )
            
            # Output data lightcurves:
            file_path = path.join( params['output_path'], \
                            event.root_file_name()+'_earth.dat' )
            event.output_data( file_path )
            file_path = path.join( params['output_path'], \
                            event.root_file_name()+'_swift.dat' )
            swift_event.output_data( file_path )
            
            # Output model lightcurves:
            file_path = path.join(params['output_path'], \
                            event.root_file_name()+'_earth.model' )
            event.output_model( file_path, model='fspl' )
            file_path = path.join( params['output_path'], \
                            event.root_file_name()+'_swift.model' )
            swift_event.output_model( file_path, model='fspl' )
            log.info( '-> Completed output' )
        
        else:
            log.info( ' XX> Found existing output for grid point parameters:' )
            log.info( ' XX> ' + event.summary() )
            log.info( 'Skipping.' )
            
    log.info( 'Completed simulation' )
    
    log_utilities.end_day_log( log )
    
def construct_grid( params ):
    """Function to return a list of gridpoints.  Each list entry consists
    of a list of grid parameters:
    [u0, tE, phi, Vbas, rho]
    """

    (umin, umax, uincr) = params['u0_range']
    (temin, temax, teincr) = params['te_range']
    (phimin, phimax, phiincr) = params['phi_range']
    (vmin, vmax, vincr) = params['v_range']
    (rhomin, rhomax, rhoincr)= params['rho_range']    
    
    grid = []
    for u0 in np.arange( umin, umax, uincr ):
        for te in np.arange( temin, temax, teincr ):
            for phi in np.arange( phimin, phimax, phiincr ):
                for Vbase in np.arange( vmin, vmax, vincr ):
                    for rho in np.arange( rhomin, rhomax, rhoincr ):
                        grid.append( [u0,te,phi,Vbase,rho] )
    return grid

def parse_input_file( file_path ):
    """Function to parse the input file of simulation parameters into a 
    dictionary of the required format.
    Parameters:
        u0_range  min  max  incr    [units of RE]
        te_range  min  max  incr    [days]
        phi_range min  max  incr    [deg]
        v_range   min  max  incr    [mag]
        rho_range min  max  incr    [units of RE]
        cadence   float             [mins]
        lc_length float             [days]
        lens_mass float             [Msol]
        lens_distance float         [pc]
        source_distance float       [pc]
   
    """
    
    if path.isfile( file_path ) == False:
        print 'Error: Cannot find input parameter file ' + file_path
        exit()
    
    lines = open( file_path, 'r' ).readlines()
    params = {}
    for line in lines: 
        entries = line.split()
        key = str( entries[0] ).lower()
        if 'range' in key:
            rmin = float( entries[1] )
            rmax = float( entries[2] )
            incr = float( entries[3] )
            value = [ rmin, rmax, incr ]
        else:
            try:
                value = float( entries[1] )
            except ValueError:
                value = str( entries[1] )
        params[ key ] = value
    return params

#############################################
# COMMANDLINE RUN SECTION
if __name__ == '__main__':
    
    if len(argv) == 1:
        print """Call sequence:
                python swift_grid_simulation.py [parameter_file]
              """
    else:
        file_path = argv[1]
        params = parse_input_file( file_path )
        simulate_grid_models( params )