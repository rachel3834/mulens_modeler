################################################################################
#     	      	      	    MICROLENSING EVENT CLASS

##########################
# IMPORTED MODULES
from sys import exit
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants, coordinates, units
from astropy.time import Time, TimeDelta
from pyslalib import slalib

##########################
# PHYSICAL AND CALCULATED CONSTANTS

# The eccentricity of Earth's orbit
earth_orbit_e = 0.017

G = constants.G.value
M_earth = constants.M_earth.value
R_earth = constants.R_earth.value

class MicrolensingEvent():
    """Class to describe a microlensing event."""
    
    # Initialise object instance:
    def __init__(self):
        
        # Projected minimum separation of lens and source 
        # [units R_E, lens plane]
        self.u_o = None
        	
        # Time of minimum projected separation of lens and source:
        self.t_o = None
        	
        # Einstein crossing time, [d]:
        self.t_E = None
        	
        # Projected separation of lens and source [units R_E, lens plane] 
        #as a function of time
        self.u_t = None
        	
        # Angular size of the source star in units of R_E
        self.rho = None
        
        # Source-lens relative velocity (assumed to be straight line motion) 
        # [ms^-1]
        self.v_SL = None
        	
        # Arrays of the lens-source relative motion in the lens plane as 
        # functions of time:
        self.x_lens = None
        self.y_lens = None
        	
        # Angle between the source's trajectory and ecliptic north. 
        # For the time being, assuming all motion to be in the plane
        # of the ecliptic, so phi = 90deg:
        self.phi = 0.0
        	
        # Distances from the observer->lens and source and the distance 
        # between lens and source:
        self.D_L = None
        self.D_S = None
        self.D_LS = None
        	
        # Lens physical parameters:
        self.M_L  = None
        	
        # Sky coordinates and ecliptic latitude of source:
        self.RA = None
        self.Dec = None
        self.beta = None
        
        # Source baseline magnitude:
        self.mag_base = None
        self.mag = None
        self.merr = None
        
        # Lens Einstein radius:
        self.R_E = None
        	
        # Array of the magnification of the event as a function of time:
        self.A_t_pspl = None
        self.A_t_fspl = None
        self.A_t_pspl_parallax = None
        self.A_t_pspl_parallax_satellite = None
            
        # Parameterization of Earth's orbit around the Sun:
        # t_p = time of perihelion
        # t_c = time when the Earth is closest to the Sun-source line
        # omega = time of maximum light as an angular frequency
        # big_omega = parameter of orbital motion of Earth around the Sun 
        # projected into the lens plane
        self.t_c = None
        self.t_p = None
        self.omega = None
        self.big_omega_o = ( 2.0 * np.pi ) / 365.25       # [rads day^-1] 
        self.omega_dt = None
        self.omega_sq_dt_sq = None
        self.big_omega_dtc = None
        self.x_earth = None
        self.y_earth = None
        self.earth_helio_position = []
        self.earth_helio_radius = None
        	
        # Projected position of the observer in the lens plane, units of R_E:
        self.x_obs = 0.0
        self.y_obs = 0.0
        	
        # Parameterization of Swift's orbit around Earth:
        # Radius of Swift's orbit from Earth centre (converting nautical 
        # mile altitude->m):
        self.swift_orbit_radius = ( 324.0 * 1852.0 ) + R_earth
        self.swift_orbit_period = np.sqrt( ( 4.0 * np.pi**2 * \
                               self.swift_orbit_radius**3 ) / ( G * M_earth ) )
        self.omega_swift = 1.0 / self.swift_orbit_period
        self.swift_t = None
	
    ###############################
    # CALCULATE MAGNIFICATION
    def calc_magnification(self,model='pspl'):
        """Method to calculate the magnification during a microlensing event.
        PSPL expression taken from Mao & Pacinzski
        FSPL computed according to the method of Lee, C-H et al., 2009, ApJ, 
        695, 200.
        """
        dbg = False
        
        def calc_f( theta, u, rho ):
            """Function to evaluate the value of f(theta) where
            f(theta) = u2(theta)sqrt(u2(theta)**2 + 4) -
                   u1(theta)*sqrt(u1(theta)**2 + 4)
            """
            # Calculate integration limits:
            if u <= rho:
                u1 = 0.0
                
                u2 = ( u * np.cos( theta ) ) + \
                       ( np.sqrt( rho*rho - u*u*(np.sin(theta))**2 ) )
                  
            else:
                if theta <= np.arcsin( rho / u ):
                    u1 = ( u * np.cos( theta ) ) - \
                        ( np.sqrt( rho*rho - u*u*(np.sin(theta))**2 ) )
                  
                    u2 = ( u * np.cos( theta ) ) + \
                        ( np.sqrt( rho*rho - u*u*(np.sin(theta))**2 ) )
                else:
                    u1 = 0.0
                    u2 = 0.0
                 
            f = ( u2 * np.sqrt( u2*u2 + 4.0 ) ) - \
                       ( u1 * np.sqrt( u1*u1 + 4.0 ) )
            
            return f
              
        if model == 'pspl':
            A_t = ( self.u_t * self.u_t + 2 ) \
                / ( self.u_t * np.sqrt( self.u_t * self.u_t + 4 ) )
          
        elif model == 'fspl':
            n = 50.0
            A_t = np.zeros( [len(self.u_t)] )
            for i,u in enumerate(self.u_t):
                up = u + self.rho
                udiff = u - self.rho
                if u <= self.rho:
                    factors = ( 1.0 / (np.pi * self.rho * self.rho) ) * \
                              ( np.pi / ( 2.0 * n ) )
                    term1 = ( ( up * np.sqrt( up*up + 4.0 ) ) \
                             - ( udiff * np.sqrt( udiff*udiff + 4.0 ) ) ) \
                             / 3.0
                      
                    term2 = 0.0
                    for k in np.arange(1.0,n-1.0,1.0):
                        theta = ( 2.0 * k * np.pi ) / ( 2.0 * n )
                        f = calc_f( theta, u, self.rho )
                        term2 = term2 + f
                    term2 = term2 * (2.0/3.0)
                    
                    term3 = 0.0
                    for k in np.arange(1.0,n,1.0):
                        theta = ( ( 2 * k - 1.0 ) * np.pi ) / ( 2.0 * n )
                        f = calc_f( theta, u, self.rho )
                        term3 = term3 + f
                    term3 = term3 * ( 4.0 / 3.0 )
                      
                    A_t[i] = factors * ( term1 + term2 + term3 ) 
                    if dbg == True:                    
                        print 'u ('+str(u)+') <= rho ('+str(self.rho)+\
                        '), A(t) = '+str(A_t[i])
                        print factors, term1, term2, term3
                else:
                    factors = ( 1.0 / (np.pi*self.rho*self.rho) ) \
                                * ( np.arcsin( self.rho / u ) / n )
                    term1 = ( ( up * np.sqrt( up*up + 4.0 ) ) \
                             - ( udiff * np.sqrt( udiff*udiff + 4.0 ) ) ) \
                              / 3.0
                      
                    term2 = 0.0
                    for k in np.arange(1.0, ((n/2.0)-1.0), 1.0):
                        theta = ( 2.0 * k * np.arcsin( self.rho / u ) ) / n
                        f = calc_f( theta, u, self.rho )
                        term2 = term2 + f
                    term2 = term2 * (2.0/3.0)
                      
                    term3 = 0.0
                    for k in np.arange( 1.0, (n/2.0), 1.0 ):
                        theta = ( ( (2.0 * k) - 1.0 ) * \
                                   np.arcsin( self.rho/u ) ) / n
                        f = calc_f( theta, u, self.rho )
                        term3 = term3 + f
                    term3 = term3 * (4.0/3.0)
                      
                    A_t[i] = factors * ( term1 + term2 + term3 )
                    if dbg == True:                    
                        print 'u ('+str(u)+') > rho ('+str(self.rho)+\
                        ') A(t) = '+str(A_t[i])
        
        if model == 'pspl': self.A_t_pspl = A_t
        elif model == 'fspl': self.A_t_fspl = A_t
        elif model == 'pspl_parallax': self.A_t_pspl_parallax = A_t
        elif model == 'pspl_parallax_satellite': 
            self.A_t_pspl_parallax_satellite = A_t
            	
    ###############################
    # CALCULATE EINSTEIN RADIUS
    def calc_einstein_radius(self):
        """Method to calculate the Einstein radius of the lens in the lens 
        plane in physical units.
        Returns R_E in [m]"""
	
        nom = ( 4.0 * constants.G * self.M_L * self.D_L * self.D_LS )
        denom = ( constants.c * constants.c * self.D_S )
        self.R_E = np.sqrt( nom / denom )
    

    ################################
    # CALCULATE LENS-SOURCE SEPARATION
    def calc_D_lens_source(self):
        """Method to calculate the physical separation distance of lens and 
        source."""
	
        if self.D_S != None and self.D_L != None:
            self.D_LS = self.D_S - self.D_L

    #################################
    # CALCULATE SOURCE ECLIPTIC LATITUDE
    def calc_beta(self):
        """Method to calculate the ecliptic latitude of the source from 
        its coordinates."""
	
        sky_coords = coordinates.ICRS(ra=self.RA,dec=self.Dec, \
                         unit=(units.hourangle, units.degree))
        self.beta = sky_coords.galactic.b
        #print 'Beta = ',self.beta, ' cos(beta) = ',(np.cos(self.beta.value))
	
    ################################
    # CALCULATE EARTHS ORBITAL MOTION
    def calc_earth_orbit(self):
        """Method to calculate the parameterization of Earths orbital 
        motion"""
	
    	# Code from original implementation using parameterized Earth orbit. 
    	# Kept here for reference.  
    	# Calculate omega = Einstein crossing time as an angular frequency in units of days:
    	#self.omega  = 1.0/ ( self.t_E.value / ( 60.0 * 60.0 * 24.0 ) )
    	#print 'Time of maximum light as an angular frequency (omega) = ',
         #self.omega,' t_E = ',self.t_E
    	
    	# t_c = time when observer is closest to the Sun-source line
    	# t_p = time of Earth perihelion
            # Orbital motion of Earth around the Sun projected into the lens plane
    	# big_omega terms pre-calculated in arrays of the orbital motion parameter for 
    	# each timestep in the lightcurve:
    	#dtc = self.t - self.t_c.jd
    	#dtp = self.t - self.t_p.jd
    	#dt = self.t - self.t_o.jd
    	#self.omega_dt = self.omega * dt
    	#self.omega_sq_dt_sq =  self.omega * self.omega * dt * dt
    	#self.big_omega_dtc = self.big_omega_o * dtc + 2.0 * earth_orbit_e * np.sin( self.big_omega_o * dtp )
    	#print 'big_omega_dtc: ',self.big_omega_o, dtp, np.sin( self.big_omega_o * dtp )
    	
    	# Calculate projected position of Earth in lens plane for all timestamps
    	# during the event:
    	#self.x_earth = self.alpha * np.sin(self.beta.value) * np.cos(self.big_omega_dtc) 
    	#self.y_earth = self.alpha * np.cos(self.beta.value) * np.sin(self.big_omega_dtc) 
    	#print ' beta: ',np.sin(self.beta.value), np.cos(self.beta.value)
    	
    	#print 'x,y_earth: ',self.x_earth, self.y_earth
    	
    	# Calculate the position of Earth from its known orbital elements for all timestamps
    	# in the event timeline. SLALIB does not handle np.arrays, only taking the first entry, 
    	# thus a loop is necessary here.  The result is then converted from units of AU -> m.
        R_earth = []
        for ts in self.t:
            ( ph, vh, pb, vb ) = slalib.sla_epv(ts-2400000.5)
            self.earth_helio_position.append(ph)
            R_earth.append( np.sqrt( ( ph * ph ).sum() ) )
        self.earth_helio_radius = R_earth * constants.au
    	#print self.earth_helio_radius
	
    ################################
    # CALCULATE PROJECTED OBSERVER POSITION
    def calc_proj_observer_pos(self,parallax=False,satellite=False):
        """Method to calculate the position of the observer projected 
        to the lens plane, factoring in both annual parallax and 
        satellite parallax."""
	
    	# By default, the observer is assumed to be at the centre of the 
        # Earth, 1AU away from the Sun-source line:
        if parallax == False:
            r_obs = np.array( [ constants.au.value ] * len(self.t) )
	
    	# Including Earth's true orbital motion:
        if parallax == True:
    	    # Calculate the Earth's radius from the Sun during the event:
    	    # Technically this currently assumes that Earth-bound observers 
    	    # are at the centre of the planet
            self.calc_earth_orbit()
    	
    	    # Radius of the observer from the centre of the Sun, factoring 
            # in the orbit of the satellite if requested:
            r_obs = self.earth_helio_radius
	
    	# Including the additional altitude of a satellite orbiting the Earth:
        if satellite == True: 
            r_obs = r_obs.value + self.swift_orbit_radius + \
                     ( self.omega_swift * ( self.t - self.swift_t ) )
	
    	# Identify the time of closest approach of the observer to the 
    	# Sun-source line:
        i = np.where( r_obs == r_obs.min() )
        self.t_c = self.t[i][0]
    	#print 'Time of closest approach of observer to Sun-source line: ', \
         #self.t_c
    	
    	# Projecting this radius to the lens plane:
        self.alpha = ( r_obs * ( 1.0 - ( self.D_L / self.D_S ) ) ) \
                     / self.R_E
    	#print 'alpha = ',self.alpha.min()
          	#self.alpha = ( self.alpha * np.pi ) / 180.0
            
    	# Calculate omega = Einstein crossing time as an angular frequency 
        # in units of days:
        dt = self.t - self.t_o.jd
        self.omega  = 1.0/ ( self.t_E.value / ( 60.0 * 60.0 * 24.0 ) )
        self.omega_dt = self.omega * dt
        self.omega_sq_dt_sq =  self.omega * self.omega * dt * dt
    	
    	# Calculate Omega(t-tc) for all timestamps during the event:
        dtc = self.t - self.t_c
        dtp = self.t - self.t_p.jd
        self.big_omega_dtc = self.big_omega_o * dtc + 2.0 * \
                 earth_orbit_e * np.sin( self.big_omega_o * dtp )
    	
    	# Factoring in the ecliptic latitude, beta, to arrive at the 
        # projected observer location in the lens plane as a function 
        # of time during the event:
        self.calc_beta()
        self.x_obs = self.alpha.value * np.sin(self.beta) * \
                     np.cos(self.big_omega_dtc)
        self.y_obs = self.alpha.value * np.cos(self.beta) * \
                     np.sin(self.big_omega_dtc)
        self.r_obs = np.sqrt( ( self.x_obs*self.x_obs ) + \
                    ( self.y_obs*self.y_obs ) )
	
    ###############################
    # POINT-SOURCE, POINT-LENS
    
    # Calculate source-lens relative motion for a generic PSPL event:
    def calc_source_lens_rel_motion(self):
        """Method to calculate the source-lens relative motion for a 
        generic PSPL event"""
	
    	# Calculate the relative source-lens velocity [ ms^-1] 
        # based on the given t_E:
        self.v_SL = ( 2.0 * self.R_E.value ) / self.t_E.value
    	#print 'Lens-source relative velocity = ',self.v_SL,'m/s'
    	#print 'Einstein radius = ',self.R_E, ' crossing time = ',
        #(self.t_E/ ( 60.0 * 60.0 * 24.0 )),'d'
    	
    	# Assume the motion is such that the components of it in the x, 
        # y directions in the lens plane are equal.  Convert the 
        # velocity units into R_E/d:
        nu = 2.0 / ( self.t_E.value / ( 60.0 * 60.0 * 24.0 ) ) 
        self.v_SL_dx = nu * np.cos(self.phi)
        self.v_SL_dy = nu * np.sin(self.phi)
    	
    	# Calculate the time-dependent position of the source relative 
    	# to the lens in the lens plane in units of the Einstein radius:
        x_incr = ( 1.0 + 1.0 )/len(self.t)
        self.x_lens = np.arange(-1.0,1.0,x_incr)
        self.y_lens = ( self.x_lens * np.tan(self.phi) ) + \
                 ( self.u_o * np.cos(self.phi) )
        self.r_lens = np.sqrt(self.x_lens*self.x_lens + \
                 self.y_lens*self.y_lens)
          	
    ################################
    # EVENT TIMELINE
    def gen_event_timeline(self,cadence=None):
        """Method to generate the event timeline, an array of timestamps 
        for datapoints on the event lightcurve, 
        spanning the range t_o +/- (t_E/2).
        Optional argument:
            cadence float [mins]
        """
	
    	# Compute the start and end times of the event, note these are in 
        # UTC format:
    	# Plots extend to +/- 2 * t_E
        half_t_E = TimeDelta( ( self.t_E.value * 2.0 ), format='sec' )
        event_start = self.t_o-half_t_E
        event_end = self.t_o+half_t_E
        if cadence == None:
            t_incr = TimeDelta( ( self.t_E.value / 5000.0 ) , format='sec')
        else:
            t_incr = TimeDelta( ( cadence * 60.0 ), format='sec')
    	
    	# Generate an array of incremental timestamps throughout the event 
        # in JD:
        ts = []
        t = event_start
        while t <= event_end:
            t = t + t_incr
            ts.append(t.jd)
            self.t = np.array(ts)
    	
    def calc_pspl_impact_param(self):
        """Method to plot (and calculate if necessary) the point-source, 
        	point-lens lightcurve to file."""
	
    	# Calculate the impact parameter as a function of time
    	# throughout the event, assuming uniform-speed, straight-line 
    	# relative motion of lens and source:
        self.calc_source_lens_rel_motion()
        dt = ( self.t - self.t_o.jd ) / self.t_E
        self.u_t = np.sqrt( ( dt.value * dt.value ) + \
                     ( self.u_o * self.u_o ) )
    	#print self.u_t

    def calc_pspl_curve(self):
        """Method to compute the magnification of the event for a given set of timestamps.
        Requires the event instance to be fully configured.
        """
	
    	# Re-calculate the event lightcurve with these timestamps:
        self.calc_source_lens_rel_motion()
        self.calc_pspl_impact_param()
        self.calc_magnification()
	
    ###############################
    # PSPL WITH PARALLAX DUE TO EARTH'S ORBITAL MOTION
    
    # Calculate PSPL lightcurve with parallax:
    def calc_pspl_parallax_impact_param(self):
        """Method to calculate a point-source, point-lens lightcurve given 
        the lens parameters, including the parallax due to Earths orbital 
        motion.
        Formulation derived from Buchalter, A., Kamionkowski, M., 1997, 
        ApJ, 483, 782.
        """
        dbg = False
   
    	# Pre-calculate some repeated factors for convenience and speed:
        alpha_sq = self.alpha.value * self.alpha.value
	
    	# Calculate impact parameter as a function of time:
        term1 =  self.u_o*self.u_o
        term2 = self.omega_sq_dt_sq
        term3 = alpha_sq * ( np.sin( self.big_omega_dtc ) )**2 
        term4 = alpha_sq * ( ( np.sin( self.beta.value ) )**2 ) * \
                     ( ( np.cos(self.big_omega_dtc) )**2 ) 
        term5 = 2.0 * self.alpha.value * np.sin(self.big_omega_dtc) * \
                         ( self.omega_dt * np.sin(self.phi) + \
                             self.u_o * np.cos(self.phi) ) 
        term6 = 2.0 * self.alpha.value * np.sin(self.beta.value) * \
                     np.cos(self.big_omega_dtc) * \
                     ( self.u_o * np.sin(self.phi) - self.omega_dt * \
                     np.cos(self.phi) )
    
        usq_t = term1 + term2 + term3 + term4 + term5 + term6
    	
    	# Print these terms at t_o?  Or find when they are all minimised? 
        # Term 3 largest at t=t_o
        if dbg == True:
            i = np.where(self.t == self.t_o.jd)
            print 'u(t)^2 terms: ',term1, \
    	      	      	term2[i], \
    			term3[i], \
    			term4[i], \
    			term5[i],\
    			term6[i]
    	
    	#dx = self.x_lens - self.x_obs
    	#dy = self.y_lens - self.y_obs
    	
    	#usq_t = np.sqrt( dx*dx + dy*dy )
    	
    	#print usq_t
    	
        self.u_t = np.sqrt(usq_t) 
    	
        #print self.u_t

    ##################################
    # OBSERVATION SIMULATION
    def simulate_data_points(self,model='pspl', phot_precision='1m'):
        """Method to extract a series of datapoints from the current model
        with a given photometric cadence and precision."""
        
        if model == 'pspl': 
            A_t = self.A_t_pspl
        if model == 'fspl': 
            A_t = self.A_t_fspl
        if model == 'pspl_parallax': 
            A_t = self.A_t_pspl_parallax
        if model == 'pspl_parallax_satelite': 
            A_t = self.A_t_pspl_parallax_satellite
        
        for A in A_t:
            mag = self.mag_base + A
            merr = self.sim_mag_error( mag, phot_precision )
        self.mag = np.array( mag )
        self.merr = np.array( merr )

    def sim_mag_error( self, mag, precision_model='1m'):
        """Method to approximate the photometric precision possible
        for a given telescope"""
        
        # Select simulation parameters depending on the class
        # of telescope selected:
        if precision_model == '1m':
            ZP = 23.0
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
            ZP = 19.0
            G = 1.0
            aperradius = 8.0
            RDN = 0.0
            skybkgd = 50.0
            teldiam = 0.3
            scintillation_noise = False
        height_o = 8000.0
        
        flux = ( 10**( ( mag - ZP ) / -2.5 ) ) * G
        logfactor = 2.5 * (1.0 / flux) * np.log10(np.exp(1.0))
        npix_aper = np.pi*aperradius*aperradius
        sig_Read = np.sqrt(RDN*RDN*npix_aper)*logfactor
        var_Read = sig_Read*sig_Read
        invvar = 1.0/var_Read
        readnoise = 1.0/np.sqrt( invvar )
        var_Sky = skybkgd * G * npix_aper
        sig_Sky = np.sqrt(var_Sky)*logfactor
        var_Sky = sig_Sky*sig_Sky
        invvar = 1.0/var_Sky
        skynoise = 1.0/np.sqrt( invvar )
        sig_Star = np.sqrt(flux)*logfactor
        var_Star = sig_Star * sig_Star
        invvar = 1.0/var_Star
        starnoise = 1.0/np.sqrt( invvar )
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
        
        return merr
        
    ##################################
    # PLOT LENS PLANE MOTION
    def plot_lens_plane_motion(self):
        """Method to plot the motion of the lens as a function of time in 
        the lens plane geometry"""
	
    	# Function to describe the circle representing the Einstein ring:
        def circle(r):
            phi=np.arange(0,6.28,0.01)	    # 0 -> 2.pi
            xcir = r*np.cos(phi)
            ycir = r*np.sin(phi)
            return xcir, ycir
    	
        # Set up figure:				   
        fig = plt.figure(2,(10,12)) 
        ax = fig.add_axes([0.1, 0.3, 0.85, 0.65])   #  [left, bottom, width, height]
	
    	# Plot the trajectory of the lens relative to the Sun-source 
        # line during the event:
        plt.plot(self.x_lens,self.y_lens,'r-')
	
    	# Plot the trajectory of the observer in the lens plane during 
        # the event:
        plt.plot(self.x_obs,self.y_obs,'b-')
	
        plt.xlabel('x [R$_{E}$]', fontsize=18)
        plt.ylabel('y [R$_{E}$]', fontsize=18)
        (x_E, y_E) = circle(1.0)      # Plot in units of the Einstein radius
        plt.plot(x_E,y_E,'k-.')
        plt.axis([-2.0,2.0,-2.0,2.0])
        plt.box('off')
        plt.grid(True)
        plt.axhline(y=0, color='black')
        plt.axvline(x=0, color='black')
    
        # Plot lens and observer separation as functions of time:
        ax = fig.add_axes([0.1, 0.05, 0.85, 0.2])   #  [left, bottom, width, height]
        plt.plot(self.t-2450000.0,(self.y_lens),'r-')
        plt.plot(self.t-2450000.0,(self.y_obs),'b-')
        plt.plot(np.array( [self.t_o.jd-2450000.0]*2 ), \
                         np.array([-2.0,2.0]),'r-.')
        plt.xlabel('JD-2450000.0', fontsize=18) 				     
        plt.ylabel('y$_{lens}$ [R$_{E}$]', fontsize=18) 
        ax.yaxis.grid() #vertical lines
        ax.xaxis.grid() #horizontal lines		    
    
        plt.savefig('lens_plane_motion.png')
        plt.close(2)
