## Program to simulate short-tE event observing sequence.

from astropy.time import Time, TimeDelta
import mulens_class
from astropy import constants
import matplotlib.pyplot as plt
from sys import argv, exit
from os import path
import numpy as np
import copy

def sim_short_te( params ):
    """Main driver function to simulate observing sequences for short-tE events
    """
    
    # Create event and configure its parameters:
    event = mulens_class.MicrolensingEvent()
    for key, value in params.items():
        setattr(event,key,value)
    
    # Compute the microlensing lightcurve:
    event.calc_D_lens_source()
    event.calc_einstein_radius()
    event.gen_event_timeline()
    event.calc_pspl_curve()
    mag = params['I_base'] - 2.5 * np.log10( event.A_t_pspl )
    
    # Simulate photometry from observing sequence.  
    # This is done by taking a copy of the event, then resetting the timestamps array, 
    # and re-calculating the lightcurve for just those timestamps:
    I_event = copy.copy(event)
    time_start = params['t_obs_start']
    time_end = time_start + TimeDelta( (I_event.t_E.value / 2.0), format='sec' )
    visit_duration = TimeDelta( (30.0*60.0), format='sec' )
    visit_cadence = TimeDelta( 355.0, format='sec' )
    I_event.t = gen_time_sequence( time_start, time_end, params['exposure_sequences'], params['visit_intervals'] )
    I_event.calc_pspl_curve()
    I_mag = params['I_base'] - 2.5 * np.log10( I_event.A_t_pspl )
    
    V_event = copy.copy(event)
    time_start = params['t_obs_start'] + visit_duration
    time_end = time_start + TimeDelta( (I_event.t_E.value / 2.0), format='sec' )
    visit_duration = TimeDelta( (30.0*60.0), format='sec' )
    visit_cadence = TimeDelta( 355.0, format='sec' )
    V_event.t = gen_time_sequence( time_start, time_end, params['exposure_sequences'], params['visit_intervals'] )
    V_event.calc_pspl_curve()
    V_mag = params['I_base'] - 2.5 * np.log10( V_event.A_t_pspl )
    
    # Plot event lightcurve:
    def select_plot_data( time_stamps, data, t_min, t_max ):
        """Function to select from the arrays given those datapoints within the time range.
	Timestamps should have been corrected for any plotting offsets before
	the function is called (e.g. ts-2450000.0).
	"""
	
        i = np.where( time_stamps >= t_min )
        j = np.where( time_stamps <= t_max )
        idx = np.intersect1d( i, j )
        
	return time_stamps[idx], data[idx]
   
    font = 22
    fig = plt.figure(1,(12,12))
    ax = fig.add_axes([0.15, 0.55, 0.775, 0.35])   #  [left, bottom, width, height]
    dt = 1.0
    t_min = event.t_o.jd - 2450000.0 - dt
    t_max = event.t_o.jd - 2450000.0 + dt
    (I_xplot, I_yplot) = select_plot_data( I_event.t-2450000.0, I_mag, t_min, t_max )
    (V_xplot, V_yplot) = select_plot_data( V_event.t-2450000.0, V_mag, t_min, t_max )
    (model_xplot, model_yplot) = select_plot_data( event.t-2450000.0, mag, t_min, t_max )
    plt.plot( model_xplot, model_yplot, 'k-', label='PSPL' )
    plt.plot( I_xplot, I_yplot, 'rd',label='I data' )
    plt.plot( V_xplot, V_yplot, 'bv',label='V data' )
    (xmin,xmax,ymin,ymax) = plt.axis()
    plt.xlabel('JD-2450000.0', fontsize=font)						   
    plt.ylabel('Magnitude', fontsize=font)
    plt.title('Simulated lightcurve of a $t_{E}$=' + \
             str( round( ( event.t_E.value/(24.0*60.0*60.0) ), 1) ) + 'd event, ', fontsize=font)				   
    plt.legend(loc='upper right',frameon=False, numpoints=1) 
    ax.tick_params(labelsize=font)
    plt.axis([xmin,xmax,ymax,ymin])						   
    
    # Plot time differential of event lightcurve
    # Note the gradient of the lightcurve is calculated this way and not in normal Python array difference
    # style because the interval between visits is variable.  Normal array subtraction ends up differencing
    # datapoints from very different points in the lightcurve, and stepping over the array would be complicated
    # due to the variable visit duration. 
    dt = 1.0
    def calc_lc_gradient( event, mag, interval ):
        grad = []
	ts = []
        for i in range( 0, len( event.t )-2, 2 ):
            if event.t[i] - event.t[i+1] < interval:
	        grad.append( ( mag[i] - mag[i+1] ) / ( event.t[i] - event.t[i+1] ) )
		ts.append( event.t[i] )
	return np.array( grad ), np.array( ts )
	
    ( grad_I, ts_I ) = calc_lc_gradient(I_event, I_mag, params['visit_intervals'][0].value)
    ( grad_V, ts_V ) = calc_lc_gradient(V_event, V_mag, params['visit_intervals'][0].value)
    ax = fig.add_axes([0.15, 0.1, 0.775, 0.35])   #  [left, bottom, width, height]
    plt.plot( ts_I-2450000.0, grad_I, 'rd')
    plt.plot( ts_V-2450000.0, grad_V, 'bv')
    plt.xlabel('JD-2450000.0', fontsize=font)		   
    plt.ylabel('Gradient [mag/d]', fontsize=font)
    plt.title('Lightcurve rate of change', fontsize=font)
    ax.tick_params(labelsize=font)
    (xmin,xmax,ymin,ymax) = plt.axis()
    xmin = I_event.t_o.jd - 2450000.0 - dt
    xmax = I_event.t_o.jd - 2450000.0 + dt
    plt.axis([xmin,xmax,ymax,ymin])	
    
    
    plt.savefig(params['plot_file'])
    
    
def read_event_file( file_path ):
    """Function to read the event parameters from a file"""
    
    if path.isfile( file_path ) == False:
        print 'Error: Cannot find input file ' + file_path
    
    file_lines = open( file_path, 'r' ).readlines()
    
    # Parse parameters to convert to standard units:
    print 'Input parameters:'
    params = {}
    for line in file_lines:
        ( key, value ) = line.replace('\n','').split( '=' )
	key = key.lstrip().rstrip()
	value = str( value ).lstrip().rstrip()
	if key in [ 'u_o', 't_E', 'M_L', 'D_L', 'D_S', 'phi', 'I_base' ]: 
	    value = float(value)
	if key == 't_E':
	    value = TimeDelta((value * 24.0 * 3600.0),format='sec')
	elif key == 'M_L':
	    value = constants.M_sun * value
	elif key == 'D_L' or key == 'D_S':
	    value = value * constants.pc
	elif key == 'phi':
	    value = ( value * np.pi ) / 180.0 
	elif key in [ 't_o', 't_obs_start' ]:
	    value = Time( value , format='isot', scale='utc' )
	elif key == 'visit_intervals':
	    tlist = value.split(',')
	    value = []
	    for entry in tlist: value.append( TimeDelta( (float( entry )*3600.0), format='sec' ) )
	elif key == 'exposure_sequences':
	    tlist1 = value.split(',')
	    value = []
	    for entry in tlist1: 
	        tlist2 = entry.split(':')
		sequence = []
		for exp in tlist2:
		    sequence.append( TimeDelta( float( exp ), format='sec' )  )
		value.append( sequence )
	params[key] = value
        print key, value
	
    return params

def gen_time_sequence(time_start, time_end, exposure_sequences, visit_intervals ):
    """Function to generate timestamps of simulated data"""
    
    # Generate an array of incremental timestamps throughout the event in JD:
    ts = []
    t = time_start
    v = -1
    while t <= time_end:
        
	# 1 visit consists of a set of sequential exposures:
	v = v + 1
	if v >= len(exposure_sequences):
	    sequence = exposure_sequences[-1]
	else:
	    sequence = exposure_sequences[v]
	
	for exptime in sequence:
	    t = t + exptime
	    ts.append(t.jd)
	    
	# Then there is a gap of length visit_interval before the next visit:
	if v >= len(visit_intervals):
	    interval = visit_intervals[-1]
	else:
	    interval = visit_intervals[v]
	t = t + interval
    
    ts = np.array(ts)
    
    return ts
    
#################################################
if __name__ == '__main__':

    help_text = """Simulator for short-tE event observations.
    Useage:
    > python observe_shortte.py path_parameter_file
    """
    
    if len(argv) > 1:
        file_path = argv[1]
	params = read_event_file( file_path )
	sim_short_te(params)
    else:
        print help_text
	
