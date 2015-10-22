#!/usr/bin/python
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

###########################
# SIMULATOR 
def swift_pspl_simulator():
    '''Program to simulate point-source, point-lens microlensing events as seen
    simultaneously from Earth and Swift.'''
    
    # Instantiate a new instance of a microlensing event:
    event = mulens_class.MicrolensingEvent()
    event.u_o = 0.00005
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
    
    # t_c = time when the Earth is closest to the Sun-source line
    # For the purposes of this simulation, this is taken to be the same as t_o
    #event.t_c = Time('2015-06-25T12:00:00', format='isot', scale='utc') 
    
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
    
    # Top panel						   
    ax = fig.add_axes([0.15, 0.5, 0.775, 0.45])   #  [left, bottom, width, height]	   
    
    # Plot PSPL event seen from stationary Earth:
    plt.plot(event.t-2450000.0,event.A_t_pspl,'k-',label='PSPL')
    
    ###################################
    # PSPL + PARALLAX EVENT
    
    # Calculate the position of the observer projected to the lens plane in units of R_E:
    event.calc_proj_observer_pos(parallax=True,satellite=False)

    # Calculate the impact parameter taking Earth's orbital motion into account:
    event.calc_pspl_parallax_impact_param()
    
    # Calculate the magnification curve:
    event.calc_magnification(model='pspl_parallax')	      # event.A_t now contains PSPL+parallax lightcurve
    
    # Diagnostic plot of relative motion in the lens plane:
    event.plot_lens_plane_motion()
    
    # Plot PSPL+parallax event seen from moving Earth:
    plt.plot(event.t-2450000.0,event.A_t_pspl_parallax,'r-',label='PSPL+parallax')
    
    ###################################
    # PSPL + PARALLAX + SATELLITE EVENT
    
    satellite = True
    if satellite == True: 
        # Setting the reference time of the Swift orbit to the starting timestamp of the event:
        event.swift_t = event.t[0]
    
        # Calculate the position of the observer projected to the lens plane in units of R_E:
        event.calc_proj_observer_pos(parallax=True,satellite=True)
    
        # Calculate the impact parameter taking Earth's orbital motion into account:
        event.calc_pspl_parallax_impact_param()
    
        # Calculate the magnification curve:
        event.calc_magnification(model='pspl_parallax_satellite')	      # event.A_t now contains PSPL+parallax lightcurve
    
        # Diagnostic plot of relative motion in the lens plane:
        event.plot_lens_plane_motion()
    
        # Plot PSPL+parallax event seen from moving Earth:
        plt.plot(event.t-2450000.0,event.A_t_pspl_parallax_satellite,'b-',label='PSPL+parallax+satellite')
    
    
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
    
    # Add plot of residuals:			   
    ax = fig.add_axes([0.15, 0.125, 0.775, 0.3])   #  [left, bottom, width, height]	
    
    # Add plot of residuals PSPL+parallax - PSPL+parallax+satellite, if calculated:
    if satellite == True:
        label='A(PSPL+parallax+satellite) - A(PSPL+parallax)'
      	dmag1 = 2.5 * np.log10( event.A_t_pspl_parallax_satellite )
      	dmag2 = 2.5 * np.log10( event.A_t_pspl_parallax )
	dmag = dmag1 - dmag2
        plt.ylabel('$\Delta$(2.5 log(A)) [mag]', fontsize=18)	
        plt.plot(event.t-2450000.0,dmag,'b-',label=label)				   
        #plt.ylabel('A$_{PSPL+parallax+satellite}$ - A$_{PSPL+parallax}$', fontsize=18)	
    
    # Otherwise, plot the PSPL-PSPL+parallax curve instead:
    else:
        plt.plot(event.t-2450000.0,(event.A_t_pspl_parallax - event.A_t_pspl),'r-',label='A(PSPL+parallax)$ - A(PSPL)')					   
        plt.ylabel('A$_{PSPL+parallax}$ - A$_{PSPL}$', fontsize=18)	
    
    # Complete second plot and save to file:
    plt.axis([xmin,xmax,-0.001,0.001])	
    plt.xlabel('JD-2450000.0', fontsize=18)				   
    plt.legend(loc='upper right',frameon=False) 					   
    ax.yaxis.grid() #vertical lines
    ax.xaxis.grid() #horizontal lines						   
    plt.savefig('pspl_curve.png')								   

    # Print summary to screen:
    print 'Amax = ',(event.A_t_pspl_parallax).max(), ' u_t(min) = ',event.u_t.min(), ' u_o = ',event.u_o, ' t_o = ',event.t_o.jd


###########################
# COMMANDLINE RUN SECTION
if __name__ == '__main__':
    
    swift_pspl_simulator()
