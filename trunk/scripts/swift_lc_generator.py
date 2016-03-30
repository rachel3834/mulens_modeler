# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 17:29:41 2016

@author: robouser

################################################################################
#     	      	      SWIFT LIGHTCURVE GENERATOR
# 
# Program to re-generate Swift lightcurves for a given model lightcurve
################################################################################
"""

from os import path
from sys import exit, argv
import numpy as np
import mulens_class
import glob

def process_grid_locale( data_dir ):
    """Function to loop over all Swift models in a given directory"""
    
    data_files = glob.glob( path.join( data_dir, 'lc_*_swift.dat' ) )
    exp_time = 200.0
    
    for dat_file in data_files:
        model_file = dat_file.replace('_swift.dat','_swift.model')
        new_dat_file = dat_file.replace('_swift.dat','_swift_new.dat')
        
        out_file = open(new_dat_file,'w')
        
        (data,header) = read_file( dat_file, 3 )
        (model,tmp) = read_file( model_file, 2 )
        
        out_file = open(new_dat_file,'w')
        for line in header:
            out_file.write(line)
        
        for i in range(len(data)):
            ts = data[i,0]
            idx = np.where( model[:,0] == ts )[0][0]
            m = model[idx,1]
            sig = data[i,2]
            mag = np.random.normal( m, sig )
            merr = sim_mag_error( exp_time, mag, precision_model='swift')
            
            out_file.write(str(ts) + '  ' + str(mag) + '  ' + str(merr) + '\n')
        
        out_file.close()
        print 'Processed ' + path.basename(dat_file) + ' -> ' + \
                        path.basename(new_dat_file)
    
def read_file( file_path, ncols ):
    """Function to read ASCII lightcurve data and model files"""
    
    d = []
    header = []
    if path.isfile(file_path) == False:
        print 'ERROR: Missing file ' + file_path
    else:
        file_lines = open(file_path, 'r').readlines()
        for line in file_lines:
            if line[0:1] != '#':
                entries = line.split()
                if ncols == 2:
                    d.append( [ float(entries[0]), float(entries[1]) ] )
                elif ncols == 3:
                    d.append( [ float(entries[0]), float(entries[1]), \
                                    float(entries[2]) ] )
            else:
                header.append(line)
    d = np.array(d)
    return d,header
    
def sim_mag_error( exp_time, mag, precision_model='1m', debug=False):
    """Method to approximate the photometric precision possible
    for a given telescope"""
    
    # Select simulation parameters depending on the class
    # of telescope selected:
    if precision_model == '1m':
        ZP = 25.0
        G = 2.0
        aperradius = 8.0
        RDN = 2.5
        skybkgd = 500.0
        airmass = 1.5
        telheight = 2000.0
        teldiam = 1.0
        exptime = 200.0
        scintillation_noise = True
    elif precision_model == 'swift':
        exptime = 200.0
        ZP = 23.0
        G = 1.0
        aperradius = 5.0
        RDN = 0.0
        skybkgd = 5.0
        teldiam = 0.3
        scintillation_noise = False
    height_o = 8000.0
    
    flux = ( 10**( ( mag - ZP ) / -2.5 ) ) * G
    logfactor = 2.5 * (1.0 / flux) * np.log10(np.exp(1.0))
    if debug == True and precision_model == 'swift': 
        print flux, mag
    npix_aper = np.pi*aperradius*aperradius
    sig_Read = np.sqrt(RDN*RDN*npix_aper)*logfactor
    var_Read = sig_Read*sig_Read
    invvar = 1.0/var_Read
    readnoise = 1.0/np.sqrt( invvar )
    if debug == True and precision_model == 'swift': 
        print 'read: ',sig_Read, var_Read, invvar, readnoise
    var_Sky = skybkgd * G * npix_aper
    sig_Sky = np.sqrt(var_Sky)*logfactor
    var_Sky = sig_Sky*sig_Sky
    invvar = 1.0/var_Sky
    skynoise = 1.0/np.sqrt( invvar )
    if debug == True and precision_model == 'swift': 
        print 'sky: ',sig_Sky, var_Sky, invvar, skynoise
        
    sig_Star = np.sqrt(flux)*logfactor
    var_Star = sig_Star * sig_Star
    invvar = 1.0/var_Star
    starnoise = 1.0/np.sqrt( invvar )
    
    if debug == True and precision_model == 'swift': 
        print 'star: ', sig_Star, var_Star, invvar, starnoise
        
    if scintillation_noise == True:
        sig_Scint = 0.09 * ( (teldiam*100.0)**-0.67) * \
                        (airmass**1.5) * np.exp(-telheight/height_o) * \
                        ((2.0*exptime)**-0.5)
        var_Scint = sig_Scint * sig_Scint
        invvar = 1.0/var_Scint
        scintnoise = 1.0/np.sqrt( invvar )
    
    err_sum = (readnoise*readnoise) + (skynoise*skynoise) + \
	      	(starnoise*starnoise)
    if scintillation_noise == True:
        err_sum = err_sum + (scintnoise*scintnoise)
    merr = np.sqrt( err_sum )

    if debug == True and precision_model == 'swift': 
        print 'Noise: ',readnoise, skynoise, starnoise, merr
        #if merr > 10.0:
         #   exit()
        
    return merr

if __name__ == '__main__':
    
    if len(argv) == 1:
        grid_dir = raw_input( 'Please enter the path to the grid directory: ' )
    else:
        grid_dir = argv[1]
        
    process_grid_locale( grid_dir )
    