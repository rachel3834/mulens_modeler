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
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def plot_heat_map_3d( stats_file ):
    """Function to plot a 3D heat map of the chi squared of each grid model
    Format of input data array:
        stats_data[ uo, te, phi, mag_base, rho ]
    """
    
    stats_data = read_stats_grid_data( stats_file )
    
    fig = plt.figure(1,(12,12)) 
    ax = fig.add_subplot(111, projection='3d')   # [left, bottom, width, height]	   
    
    # Identify data points for a given pair of rho and phi:
    rx = np.where( stats_data[:,4] == 0.001 ) 
    px = np.where( stats_data[:,2] == 0.0 )
    idx = np.intersect1d(rx,px)
    
    u0 = stats_data[idx,0]
    te = stats_data[idx,1]
    mag = stats_data[idx,3]
    chi2 = stats_data[idx,5]
    
    # Scaling the colourmap:
    nshade = 20.0
    cincr = ( chi2.max() - chi2.min() )
    cmap = ( chi2 - chi2.min() ) / cincr
    
    scatter_plot = ax.scatter( u0, te, mag, c=cmap, edgecolors='face', s=100 )
    ax.set_xlabel( '$u_{O}$', fontdict={'fontsize': 18 } )
    ax.set_ylabel( '$t_{E}$', fontdict={'fontsize': 18 } )
    ax.set_zlabel( 'Star mag (baseline)', fontdict={'fontsize': 18 } )
    ax.invert_zaxis()

    cbar = plt.colorbar(scatter_plot, shrink=0.5, aspect=20)
    cbar.set_ticks
    cticks = np.arange(0.0,1.1,0.1).tolist()
    cincr = ( chi2.max() - chi2.min() ) / 11.0
    clabels = []
    for c in np.arange(chi2.min(), chi2.max(), cincr):
        clabels.append( round(c,0) )
    cbar.set_ticks(cticks)
    cbar.set_ticklabels(clabels)
    
    data_dir = path.dirname( stats_file )
    plt.savefig( path.join( data_dir,'stats_heat_map.png' ) )

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

if __name__ == '__main__':
    
    if len(argv) == 1:
        stats_file = raw_input( 'Please enter path to the statistics file: ')
    else:
        stats_file = argv[1]
        
    plot_heat_map_3d( stats_file )