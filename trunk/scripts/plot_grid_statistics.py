# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 17:29:47 2016

@author: robouser
"""

##############################################################################
#                       PLOT GRID STATISTICS
##############################################################################

from os import path
from sys import argv, exit
import numpy as np

def plot_heat_map_3d( stats_data ):
    """Function to plot a 3D heat map of the chi squared of each grid model
    Format of input data array:
        stats_data[ uo, te, phi, mag_base, rho ]
    """
    
    fig = plt.figure(1,(12,12), projection='3d') 
    ax = fig.add_axes([0.15, 0.95, 0.9, 0.9])   # [left, bottom, width, height]	   
    
    # Identify data points for a given pair of rho and phi:
    rx = np.where( stats_data[:,4] == 0.001 ) 
    px = np.where( stats_data[:,2] == 0.0 )
    idx =     
    
    u0 = stats_data[rx,0]
    te = stats_data[rx,1]
    mag = stats_data[rx,3]
    
    
    ax.scatter( stats_data[:,:,0,mag,rho], \
                stats_data[:,:,0,mag,rho], \
                stats_data[:,:,0,mag,rho] )

def read_stats_grid_data( stats_file ):
    """Function to read the file containing the statistical data from a grid of 
    models"""
    
    if path.isfile( stats_file ) == False:
        print 'Error: Cannot find input file ' + stats_file
        exit()
    
    with open( stats_file, 'r' ) as textFile:
        file_lines = textFile.readlines()
    
    data_list = []
    for line in file_lines:
        if line.lstrip()[0:1] != '#':
            entries = line.replace('\n','').split()
            line_data = []
            for value in entries:
                line_data.append( float(value) )
            data_list.append( line_data )
    
    stats_data = np.array(data_list)
    
    return stats_data

