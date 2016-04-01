# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 18:33:13 2016

@author: rachel3834
"""

##############################################################################
#                           GRID STATISTICS
##############################################################################

from sys import argv, exit
from os import path
import numpy as np
import glob

class GridPoint:
    """Class describing the parameters of a point in a grid of microlensing 
    models covering ranges in a number of parameters"""
    
    def __init__( self ):
        self.u0 = None
        self.te = None
        self.phi = None
        self.rho = None
        self.mag_base = None
        self.chi2 = None
        self.earth_chi2 = None
        self.swift_earth_chi2 = None
        self.red_chi2 = None
        self.dchi2 = None
        self.n_pts = None
        self.ml = None
        self.dl = None
        
    def set_par_from_file( self, file_name ):
        """Method to set the parameters of the grid point from the filename
        of the model output"""
        
        pars = path.basename( str( file_name ) ).split('_')
        self.u0 = float( pars[1] )
        self.te = float( pars[2] )
        self.phi = float( pars[3] )
        self.mag_base = float( pars[4] )
        self.rho = float( pars[5] )        
        self.ml = float( pars[6] )
        self.dl = float( pars[7] )
        
    def summary( self ):
        """Method to output a summary of all grid-point parameters as a string"""
        
        line = str( self.u0 ) + ' ' + str( self.te ) + ' ' + \
                str( self.phi ) + ' ' + \
                str( self.mag_base ) + ' ' + str( self.rho ) + ' ' + \
                str( self.ml ) + ' ' + str( self.dl ) + ' ' + \
                str( self.dchi2 ) + ' ' + str( self.swift_chi2 ) + ' ' + \
                str( self.swift_earth_chi2 )
        return line
        
def grid_locale_stats( grid_dir ):
    """Function to loop over all models within a grid locale sub-directory.
    Each set of models represents the data for a single grid-point, and there
    may be several grid-point models within one locale directory.    
    """

    # Produce a list of the grid-points within this locale by selecting just
    # one of the datafiles for each grid-point.  Files are systematically 
    # named with the grid-point parameters.
    earth_file_list = glob.glob( path.join( grid_dir, 'lc_*earth.dat' ) )
    
    # Open output file:
    output_file = path.join( grid_dir, 'grid_locale_stats.dat' )
    output = open(output_file, 'w')    
    output.write('# u0     tE[d]    phi[rads]   mag_base[mag]   rho  ML[MSun]    DL[pc]   dchi_sq  swift_chi2 swift-earth_chi2\n')
    
    # Loop over each file, recording the computed difference statistic between
    # the Earth- and Swift-lightcurves for the same grid point model.
    grid_stats = []
    for earth_lc_file in earth_file_list:
        gp = GridPoint()
        gp.set_par_from_file( earth_lc_file )
        
        swift_lc_file = str(earth_lc_file).replace('earth','swift_new')
        
        earth_data = read_lc_file( earth_lc_file )
        swift_data = read_lc_file( swift_lc_file )
        
        swift_model_file = str(earth_lc_file).replace('earth.dat','swift.model')
        earth_model_file = str(earth_lc_file).replace('earth.dat','earth.model')
        
        earth_model = read_lc_file( earth_model_file, model=True )
        swift_model = read_lc_file( swift_model_file, model=True )
        
        #( gp.chi2, gp.red_chi2, gp.n_pts ) = diff_lcs( earth_data, swift_data )
        
        (gp.dchi2, gp.swift_chi2, gp.swift_earth_chi2) = \
            calc_dchi2( earth_data, earth_model, swift_data, swift_model )
        
        print gp.summary()
        file_name = path.basename( earth_lc_file ).replace('_earth.dat','')
        output.write( gp.summary() + ' ' + file_name + '\n' )
        output.flush()
        
    output.close()
    
def read_lc_file( lc_file, model=False ):
    """Function to read the lightcurve output datafiles"""

    if path.isfile( lc_file ) == None:
        print 'ERROR: Cannot find lightcurve ' + lc_file
        exit()
    
    file_lines = open( lc_file, 'r' ).readlines()
    
    lc_data = []
    for line in file_lines:
        if line[0:1] != '#':
            entries = str( line ).split()
            if model == False:
                lc_data.append( [ float( entries[0] ), \
                                float( entries[1] ), \
                                    float( entries[2] ) ] )
            else:
                
                lc_data.append( [ float( entries[0] ), \
                                float( entries[1] ) ] )
    
    lc_data = np.array( lc_data )
    return lc_data 

def calc_dchi2( earth_lc, earth_model, swift_lc, swift_model ):
    """Function to calculate the delta chi 2 between the swift-swift model and 
    Earth lightcurve"""
    
    # Extract the indices of the lightcurve timestamps common to both
    # datasets. Swift dataset has the fewest points, so for consistency, 
    # use this dataset to identify which points to include:
    idx = common_ts( swift_model, swift_lc )
    
    # Compute chi2 of the Swift lightcurve with its model:
    swift_chi2 = calc_chi2( swift_model, swift_lc, idx )

    # Identify the timestamps in common between the Swift and Earth 
    # lightcurves.  As the Swift lightcurve is the least sampled this
    # selects the same datapoints as for the chi2 test above. 
    idx = common_ts( earth_model, swift_lc )
    swift_earth_chi2 = calc_chi2( earth_model, swift_lc, idx )
    
    dchi2 = swift_chi2 - swift_earth_chi2
    return dchi2, swift_chi2, swift_earth_chi2

def calc_chi2( lc1, lc2, idx, use_err1=False, debug=False ):
    """Calculate chi2 between lightcurves, where idx is an array giving
    the common timestamps"""
        
    merr2 = lc2[:,2]
    if use_err1 == True:
        merr1 = lc1[idx,2]
        sigma = np.sqrt( merr1*merr1 + merr2*merr2 )
    else:
        sigma = merr2
        
    if debug == True:
        print 'Sigma = ', sigma
        
    diff = lc1[idx,1] - lc2[:,1]
    diff = diff / sigma
    if debug == True:
        print 'Weighted diff = ', diff
    
    chi2 = ( ( diff * diff ).sum() )
    if debug == True:
        print 'Chi^2: ',chi2
    
    return chi2

def common_ts( lc1, lc2, debug=False ):
    """Swift has different sampling from Earth, but the datapoints in both
    lightcurves are taken at the same timestamps - the Swift lightcurve
    is just missing some timestamps. So the first task is to select those
    parts of the Earth lightcurve which have correspondig Swift datapoints
    we can compare with. 

    The Earth lightcurve with the superset of datapoints should be given
    as lc1."""
    
    # Extract the indices of the lightcurve timestamps common to both
    # datasets. 
    idx = []
    ts1 = lc1[:,0]
    for i,ts in enumerate( lc2[:,0] ):
        if ts in ts1: 
            j = np.where( ts1 == ts )
            idx.append(j[0][0])
    
    if debug == True:
        print lc1[idx,0], lc2[:,0]
    
    return idx

def diff_lcs( lc1, lc2, debug=False ):
    """Function to compute basic statistics to estimate the degree of variation 
    between two lightcurves. 
    
    Swift has different sampling from Earth, but the datapoints in both
    lightcurves are taken at the same timestamps - the Swift lightcurve
    is just missing some timestamps. So the first task is to select those
    parts of the Earth lightcurve which have correspondig Swift datapoints
    we can compare with. 

    The Earth lightcurve with the superset of datapoints should be given
    as lc1.
    """
    
    # Extract the indices of the lightcurve timestamps common to both
    # datasets. 
    idx = common_ts( lc1, lc2 )
    
    # Compute the chi2:
    diff = lc1[idx,1] - lc2[:,1]
    if debug == True:
        print 'Diff = ', diff
    
    merr1 = lc1[idx,2]
    merr2 = lc2[:,2]
    sigma = np.sqrt( merr1*merr1 + merr2*merr2 )
    if debug == True:
        print 'Sigma = ', sigma
    
    diff = diff / sigma
    if debug == True:
        print 'Weighted diff = ', diff
    
    chi2 = ( ( diff * diff ).sum() )
    red_chi2 = chi2 / len(idx)
    if debug == True:
        print 'Chi2 = ', chi2
    
    return chi2, red_chi2, len(idx)
    
    
##############################################################################
if __name__ == '__main__':
    
    if len(argv) == 1:
        grid_dir = raw_input( 'Please enter the path to the grid directory: ' )
    else:
        grid_dir = argv[1]
        
    grid_locale_stats( grid_dir )
    