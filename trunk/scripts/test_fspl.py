# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 11:39:07 2016

@author: rstreet
"""

################################################################################
#     	      	      	SWIFT MICROLENSING EVENT SIMULATION
# 
# Program to simulate point-source, point-lens microlensing events as seen
# simultaneously from Earth and Swift.
################################################################################

###########################
# IMPORT MODULES
import mulens_class
import numpy as np
from astropy import constants
from astropy import coordinates
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt

def generate_fspl():
    """Function to generate a test case of an FSPL event lightcurve"""
    
    # Essential parameters of the event
    event = mulens_class.MicrolensingEvent()
    event.u_o = 0.00005
    event.rho = 0.001
    event.t_E = TimeDelta((100.0 * 24.0 * 3600.0),format='sec')	  # Must be in JD or seconds
    event.M_L = constants.M_sun * 0.3 	    # ~K star
    event.D_L = 6000.0 * constants.pc 	    # ~halfway to the Galactic centre
    event.D_S = 8000.0 * constants.pc 	    # Source in the Bulge
    event.phi = ( 00.0 * np.pi ) / 180.0 
    event.calc_D_lens_source()
    event.calc_einstein_radius()
    event.RA = '17:57:34.0'
    event.Dec = '-29:13:15.0'
    event.t_o = Time('2015-01-04T16:00:00', format='isot', scale='utc')
    event.t_p = Time('2015-01-04T06:37:00', format='isot', scale='utc')
    
    
    # Generate the event timeline = timestamp array for the lightcurve:
    event.gen_event_timeline()
    
    # Calculate the motion of the source relative to the lens (fixed at
    # the centre of the reference frame).  Straight line motion:
    event.calc_source_lens_rel_motion()
    
    
    ###################################
    # PSPL EVENT
    
    # Calculate the position of the observer projected to the lens plane in units of R_E:
    event.calc_proj_observer_pos(parallax=False,satellite=False)
    
    # Compute the PSPL lightcurve of this event as seen from a stationary Earth. 
    #event.calc_pspl_impact_param()
    event.calc_pspl_parallax_impact_param()
    event.calc_magnification(model='pspl')
    
    # Diagnostic plot of relative motion in the lens plane:
    event.plot_lens_plane_motion()
    
    # Begin plot:
    legend_text = '$u_{O}$='+str(event.u_o)						   
    fig = plt.figure(1,(12,12)) 
    
    # Single panel						   
    ax = fig.add_axes([0.15, 0.125, 0.775, 0.45])   #  [left, bottom, width, height]	   
    
    # Plot PSPL event seen from stationary Earth:
    plt.plot(event.t-2450000.0,event.A_t_pspl,'k-',label='PSPL')
    
    
    
    
    # Complete lightcurve plot (top panel):
    dt = 2.0
    xmin = event.t_o.jd - 2450000.0 - dt
    xmax = event.t_o.jd - 2450000.0 + dt
    plt.xlabel('JD-2450000.0', fontsize=18)						   
    plt.ylabel('A(t)', fontsize=18)		
    plt.axis([xmin,xmax,0.0,4000.0])					   
    plt.legend(loc='upper right',frameon=False) 
    ax.yaxis.grid() #vertical lines
    ax.xaxis.grid() #horizontal lines
    					   
    plt.savefig('pspl_curve.png')		
    

###########################
# COMMANDLINE RUN SECTION
if __name__ == '__main__':
    
    generate_fspl()
