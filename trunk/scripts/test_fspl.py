# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 11:39:07 2016

@author: rstreet
"""

################################################################################
#     	      	      	MICROLENSING EVENT SIMULATION
# 
# Program to test the implementation of Finite Source, Point Lens events
################################################################################

###########################
# IMPORT MODULES
import mulens_class
import numpy as np
from astropy import constants
from astropy import coordinates
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
from sys import exit, argv

def generate_fspl():
    """Function to generate a test case of an FSPL event lightcurve"""
    
    # Essential parameters of the event
    event = mulens_class.MicrolensingEvent()
    event.u_min = 0.0001
    event.rho = 0.001
    event.t_E = TimeDelta((1.0 * 24.0 * 3600.0),format='sec')	  # Must be in JD or seconds
    event.M_L = constants.M_sun * 0.3 	    # ~K star
    event.D_L = 6000.0 * constants.pc 	    # ~halfway to the Galactic centre
    event.D_S = 8000.0 * constants.pc 	    # Source in the Bulge
    event.phi = ( 90.0 * np.pi ) / 180.0 
    event.calc_D_lens_source()
    event.calc_einstein_radius()
    event.RA = '17:57:34.0'
    event.Dec = '-29:13:15.0'
    event.t_o = Time('2015-06-01T16:00:00', format='isot', scale='utc')
    event.t_p = Time('2015-06-01T06:37:00', format='isot', scale='utc')
    event.mag_base = 13.0
    
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
    #event.calc_pspl_impact_param()
    event.calc_parallax_impact_param()
    print 'u0 = ',event.u_o, ' cf u_min = ', event.u_min
    event.calc_magnification(model='pspl')
    print 'Amax PSPL= ',event.A_t_pspl.max()
    print '   cf Amax PSPL = ',event.calc_pspl_A( event.u_t.min() )
    event.calc_magnification(model='fspl')
    print 'Amax FSPL = ',event.A_t_fspl.max()
    i = np.where( event.A_t_fspl == event.A_t_fspl.max() )
    print '    at uo = ',event.u_t[i[0][0]]
    
    # Diagnostic plot of relative motion in the lens plane:
    event.plot_lens_plane_motion()
    
    # Begin plot:
    legend_text = '$u_{O}$='+str(event.u_o)						   
    fig = plt.figure(1,(12,12)) 
    
    # Single panel						   
    ax = fig.add_axes([0.15, 0.125, 0.775, 0.8])   #  [left, bottom, width, height]	   
    
    # Plot PSPL event seen from stationary Earth:
    plt.plot(event.t-2450000.0,event.A_t_pspl,'k-',label='PSPL')
    A_fspl = event.A_t_fspl
    plt.plot(event.t-2450000.0,A_fspl,'b--',label='FSPL')
    
    
    # Complete lightcurve plot (top panel):
    dt = 1.5 * ( event.t_E.value / ( 24.0 * 3600.0 ) )
    (xmin,xmax,ymin,ymax) = plt.axis()
    xmin = event.t_o.jd - 2450000.0 - dt
    xmax = event.t_o.jd - 2450000.0 + dt
    plt.xlabel('JD-2450000.0', fontsize=18)						   
    plt.ylabel('A(t)', fontsize=18)		
    plt.axis([xmin,xmax,ymin,ymax])					   
    plt.legend(loc='upper right',frameon=False) 
    ax.yaxis.grid() #vertical lines
    ax.xaxis.grid() #horizontal lines
    					   
    plt.savefig('pspl_curve.png')		
    
    # Output to datafile:
    fileobj = open('pspl_fspl_curve.dat','w')
    fileobj.write('# HJD   PSPL A    FSPL A  PSPL mag  FSPL mag\n')
    mag_pspl = event.calc_mag(model='pspl')
    mag_fspl = event.calc_mag(model='fspl')
    for i,t in enumerate(event.t):
        fileobj.write( str(t) + ' ' + str(event.A_t_pspl[i]) + \
                        ' ' + str(event.A_t_fspl[i]) + \
                        ' ' + str(mag_pspl[i]) + ' ' + str(mag_fspl[i]) + '\n' )
    fileobj.close()
    
###########################
# COMMANDLINE RUN SECTION
if __name__ == '__main__':
    
    generate_fspl()
