################################################################################
#     	      	      	    MICROLENSING EVENT CLASS

##########################
# IMPORTED MODULES
from sys import exit
from os import path
from shutil import move
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants, coordinates, units
from astropy.time import Time, TimeDelta
from pyslalib import slalib
import timesubs
import utilities

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
        # u_min is the minimum projected separation of the lens from the 
        # Sun-Source line.  
        # For a PSPL event without parallax (i.e. Earth is static)
        # u_o = u_min
        # but for events with parallax, u_o (observed) may be different
        self.u_o = None
        self.u_min = None
        self.u_offset = None
        
        # Time of minimum projected separation of lens and source:
        self.t_o = None
        self.t_o_jd = None
        
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
        # functions of time, and at t0:
        self.x_lens = None
        self.y_lens = None
        self.x_lens_o = None
        self.y_lens_o = None
        	
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
        
        # Source baseline magnitude and simulated data points. 
        # Timestamps are by default in HJD, JD recording for testing
        self.mag_base = None
        self.ts = None
        self.ts_jd = None
        self.mag = None
        self.merr = None
        
        # Lens Einstein radius:
        self.R_E = None
        	
        # Array of the magnification of the event as a function of time:
        self.A_t_pspl = None
        self.A_t_fspl = None
        self.A_t_pspl_parallax = None
        self.A_t_fspl_parallax = None
        self.A_t_pspl_parallax_satellite = None
        self.A_t_fspl_parallax_satellite = None
            
        # Parameterization of Earth's orbit around the Sun:
        # t_p = time of perihelion
        # t_c = time when the Earth is closest to the Sun-source line
        # omega = time of maximum light as an angular frequency
        # big_omega = parameter of orbital motion of Earth around the Sun 
        # projected into the lens plane
        self.t_c = None
        self.t_p = None
        self.t_p_jd = None
        self.omega = None
        self.big_omega_o = ( 2.0 * np.pi ) / 365.25       # [rads day^-1] 
        self.omega_dt = None
        self.omega_sq_dt_sq = None
        self.big_omega_dtc = None
        self.x_earth = None
        self.y_earth = None
        self.earth_helio_position = []
        self.earth_helio_radius = None
        self.earth_helio_radius_o = None
        	
        # Projected position of the observer in the lens plane, units of R_E, 
        # as functions of time and at t0:
        self.x_obs = 0.0
        self.y_obs = 0.0
        self.x_obs_o = 0.0
        self.y_obs_o = 0.0
        	
        # Parameterization of Swift's orbit around Earth:
        # Radius of Swift's orbit from Earth centre (converting nautical 
        # mile altitude->m):
        self.swift_orbit_radius = ( 324.0 * 1852.0 ) + R_earth
        self.swift_orbit_period = np.sqrt( ( 4.0 * np.pi**2 * \
                               self.swift_orbit_radius**3 ) / ( G * M_earth ) )
        self.omega_swift = 1.0 / self.swift_orbit_period
        self.swift_t = None
    
    def summary(self, inc_uo=False):    
        """Method to print a summary of parameters"""
        
        ml = round( self.M_L.value/constants.M_sun.value, 6 )
        dl = round( self.D_L.value/constants.pc.value, 4 )
        if inc_uo == False:
            output = 't_E=' + str(self.t_E/(60.0*60.0*24.0)) + 'd, t_0=' + str(self.t_o) + \
                ', u_min=' + str(self.u_min) + ', rho=' + str(self.rho) + \
                ', mag_base=' + str(self.mag_base) + ', phi=' + \
                str(self.phi) + ', M_L=' + str(ml) + 'Msol, D_L=' + \
                str(dl) + 'pc'
        else:
            uo = round( self.u_o, 6 )
            uoffset = round( self.u_offset, 6 )
            output = 't_E=' + str(self.t_E/(60.0*60.0*24.0)) + 'd, t_0=' + str(self.t_o) + \
                ', u_min=' + str(self.u_min) + ', u_o=' + str(uo) + \
                ', u_offset=' + str(uoffset) + ', rho=' + str(self.rho) + \
                ', mag_base=' + str(self.mag_base) + ', phi=' + \
                str(self.phi) + ', M_L=' + str(ml) + 'Msol, D_L=' + \
                str(dl) + 'pc, R_E=' + str(self.R_E)
        return output
    
    def root_file_name( self, params=None ):
        """Method to return a string summarizing the basic FSPL parameters
        for use as a root file name"""
        
        if params == None:
            te = self.t_E.value / ( 60.0 * 60.0 * 24.0 )   
            te = round( te, 1 )
            phi = round( self.phi, 3 )
            mag = round( self.mag_base, 1 )
            uo = round( self.u_offset, 6 )
            ml = round( self.M_L.value/constants.M_sun.value, 6 )
            dl = round( self.D_L.value/constants.pc.value, 4 )
            rho = round( self.rho, 3 )
        else:
            te = round( float(params['te']),1 )
            phi = round( float(params['phi']), 3 )
            mag = round( float(params['mag_base']),1 )
            uo = round( float(params['uo']), 6 )
            ml = round( float(params['ml']), 6 )
            dl = round( float(params['dl']), 4 )
            rho = round( float(params['rho']), 3 )
            
        file_name = 'lc_' + str(uo) + '_' + str(te) + \
                        '_' + str(phi) + '_' + str(mag) +\
                        '_' + str(rho) + '_' + str(ml) +\
                        '_' + str(dl)
        return file_name
        
    ###############################
    # CALCULATE MAGNIFICATION
    def calc_pspl_A( self, u_t ):
        A_t = ( u_t * u_t + 2 ) \
            / ( u_t * np.sqrt( u_t * u_t + 4 ) )
        return A_t
        
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
        
        
        # The ratio of u_o to rho is used to judge whether finite
        # source effects will ever be significant in this lightcurve.  
        z_o = self.u_o/self.rho
        if 'pspl' in model or z_o > 10.0:
            A_t = self.calc_pspl_A( self.u_t ) 
            
        elif 'fspl' in model:
            z = self.u_t/self.rho
            idx = np.where( z > 10.0 )
            jdx = np.where( z <= 10.0 )
            n = 500.0
            A_t = np.zeros( [len(self.u_t)] )
            
            # For points in the lightcurve where z > 10.0 and finite source
            # effects should be neglible, use a PSPL model:
            A_t[idx] = self.calc_pspl_A( self.u_t[idx] )
            
            # For points in the lightcurve where z <= 10.0, use Lee's
            # approximation method:
            for j in jdx[0].tolist():
                u = self.u_t[j]
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
                      
                    A_t[j] = factors * ( term1 + term2 + term3 ) 
                    if dbg == True:                    
                        print 'u ('+str(u)+') <= rho ('+str(self.rho)+\
                        '), A(t) = '+str(A_t[j])
                        print 'Factors: ',factors, term1, term2, term3
                else:
                    factors = ( 1.0 / (np.pi*self.rho*self.rho) ) \
                                * ( np.arcsin( self.rho / u ) / n )
                    factors = factors.value
                    
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
                    
                    A_t[j] = factors * ( term1 + term2 + term3 )
                    if dbg == True:
                        print 'u ('+str(u)+') > rho ('+str(self.rho)+\
                        ') A(t) = '+str(A_t[j])
                        print 'Factors: ',factors, term1, term2, term3
        
        if model == 'pspl': 
            self.A_t_pspl = A_t
        elif model == 'fspl': self.A_t_fspl = A_t
        elif model == 'pspl_parallax': self.A_t_pspl_parallax = A_t
        elif model == 'fspl_parallax': self.A_t_fspl_parallax = A_t
        elif model == 'pspl_parallax_satellite': 
            self.A_t_pspl_parallax_satellite = A_t
        elif model == 'fspl_parallax_satellite': 
            self.A_t_fspl_parallax_satellite = A_t
        
        
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
        def calc_earth_helio_position(t):
            ( ph, vh, pb, vb ) = slalib.sla_epv(ts-2400000.5)
            r_orb_earth = np.sqrt( ( ph * ph ).sum() )
            return ph,r_orb_earth
            
        R_earth = []
        for ts in self.t:
            (ph,r) = calc_earth_helio_position(ts)
            self.earth_helio_position.append(ph)
            R_earth.append( r )
        self.earth_helio_radius = R_earth * constants.au
        
        (ph,r) = calc_earth_helio_position(self.t_o)
        self.earth_helio_radius_o = r * constants.au
        
    	#print self.earth_helio_radius
    
    def get_earth_perihelion( self ):
        """Method to return the appropriate date of Earth's perihelion 
        for the year of the to
        As this is derived from tabular data provided by
        aa.usno.navy.mil/data/docs/EarthSeasons.php
        it only covers years 2000+        
        """
        
        perihelion_table = { \
            2000: Time('2000-01-03T05:18:00', format='isot', scale='utc'),
            2001: Time('2001-01-04T08:52:00', format='isot', scale='utc'),
            2002: Time('2002-01-02T14:09:00', format='isot', scale='utc'),
            2003: Time('2003-01-04T05:02:00', format='isot', scale='utc'),
            2004: Time('2004-01-04T17:42:00', format='isot', scale='utc'),
            2005: Time('2005-01-02T00:35:00', format='isot', scale='utc'),
            2006: Time('2006-01-04T15:30:00', format='isot', scale='utc'),
            2007: Time('2007-01-03T19:43:00', format='isot', scale='utc'),
            2008: Time('2008-01-02T23:51:00', format='isot', scale='utc'),
            2009: Time('2009-01-04T15:30:00', format='isot', scale='utc'),
            2010: Time('2010-01-03T00:09:00', format='isot', scale='utc'),
            2011: Time('2011-01-03T18:32:00', format='isot', scale='utc'),
            2012: Time('2012-01-05T00:32:00', format='isot', scale='utc'),
            2013: Time('2013-01-02T04:38:00', format='isot', scale='utc'),
            2014: Time('2014-01-04T11:59:00', format='isot', scale='utc'),
            2015: Time('2015-01-04T06:36:00', format='isot', scale='utc'),
            2016: Time('2016-01-02T22:49:00', format='isot', scale='utc')
            }
    
        year = int(self.t_o.jyear)
        self.t_p = perihelion_table[year]
        
    ################################
    # CALCULATE PROJECTED OBSERVER POSITION
    def calc_proj_observer_pos(self,parallax=False,satellite=False,debug=False):
        """Method to calculate the position of the observer projected 
        to the lens plane, factoring in both annual parallax and 
        satellite parallax."""
	
        # By default, the observer is assumed to be at the centre of the 
        # Earth, 1AU away from the Sun-source line:
        if parallax == False:
            r_obs = np.array( [ constants.au.value ] * len(self.t) )
            r_obs_o = constants.au.value
            
        # Including Earth's true orbital motion:
        if parallax == True:
            # Calculate the Earth's radius from the Sun during the event:
            # Technically this currently assumes that Earth-bound observers 
            # are at the centre of the planet
            self.calc_earth_orbit()
    	
            # Radius of the observer from the centre of the Sun, factoring 
            # in the orbit of the satellite if requested:
            r_obs = self.earth_helio_radius
            r_obs_o = self.earth_helio_radius_o
            if debug==True:
                print 'Earths orbital motion at t_o: ',r_obs_o
                
        # Including the additional altitude of a satellite orbiting the Earth:
        def radius_satellite(t, sat_orbit_radius, sat_omega, sat_t, robs):
            r = robs + sat_orbit_radius + \
                     ( sat_omega * ( t - sat_t ) )
            return r
            
        if satellite == True: 
            #r_obs = r_obs.value + self.swift_orbit_radius + \
            #         ( self.omega_swift * ( self.t - self.swift_t ) )
            r_obs = radius_satellite(self.t, self.swift_orbit_radius, \
                            self.omega_swift, self.swift_t, r_obs.value)
            r_obs_o = radius_satellite(self.t_o, self.swift_orbit_radius, \
                            self.omega_swift, self.swift_t, r_obs_o.value)
            
        # Identify the time of closest approach of the observer to the 
        # Sun-source line:
        i = np.where( r_obs == r_obs.min() )
        self.t_c = self.t[i]
    	
        # Projecting this radius to the lens plane:
        # self.alpha = ( r_obs * ( 1.0 - ( self.D_L / self.D_S ) ) ) \
        #             / self.R_E
        def calc_alpha(robs, DL, DS, RE):
            alpha = ( robs/ RE ) * ( 1.0 - ( DL / DS ) ) 
            return alpha
            
        self.alpha = calc_alpha(r_obs, self.D_L, self.D_S, self.R_E)
        alpha_o = calc_alpha(r_obs_o, self.D_L, self.D_S, self.R_E)
        if debug == True:
            print 'alpha = ',self.alpha.min(), alpha_o, r_obs_o, \
            self.D_L, self.D_S, self.R_E
            
        # Calculate omega = Einstein crossing time as an angular frequency 
        # in units of days:
        dt = self.t - self.t_o
        self.omega  = 1.0/ ( self.t_E.value / ( 60.0 * 60.0 * 24.0 ) )
        self.omega_dt = self.omega * dt
        self.omega_sq_dt_sq =  self.omega * self.omega * dt * dt
    	
        # Calculate Omega(t-tc) for all timestamps during the event:
        def calc_big_omega(big_omega_o,dtc,earth_orbit_e, dtp):
            big_omega = big_omega_o * dtc + 2.0 * \
                 earth_orbit_e * np.sin( big_omega_o * dtp )
            return big_omega
        dtc = self.t - self.t_c
        dtp = self.t - self.t_p
        #self.big_omega_dtc = self.big_omega_o * dtc + 2.0 * \
        #         earth_orbit_e * np.sin( self.big_omega_o * dtp )
        self.big_omega_dtc = calc_big_omega(self.big_omega_o,dtc, \
                                        earth_orbit_e, dtp)
        dtc_o = self.t_o - self.t_c
        dtp_o = self.t_o - self.t_p
        big_omega_dtc_o = calc_big_omega(self.big_omega_o,dtc_o, \
                                        earth_orbit_e, dtp_o)
        if debug == True:
            print 'OMEGA: ',(self.big_omega_o* dtc), self.t_c,self.t, \
                earth_orbit_e, np.sin( self.big_omega_o * dtp ), \
                self.big_omega_dtc, big_omega_dtc_o
                
        # Factoring in the ecliptic latitude, beta, to arrive at the 
        # projected observer location in the lens plane as a function 
        # of time during the event:
        # self.x_obs = self.alpha.value * np.sin(self.beta) * \
        #             np.cos(self.big_omega_dtc)
        # self.y_obs = self.alpha.value * np.cos(self.beta) * \
        #             np.sin(self.big_omega_dtc)
        # self.r_obs = np.sqrt( ( self.x_obs*self.x_obs ) + \
        #            ( self.y_obs*self.y_obs ) )
        self.calc_beta()
        def calc_obs_position(alpha, beta, big_omega_dtc):
            x_obs = alpha * np.sin(beta) * \
                     np.cos(big_omega_dtc)
            y_obs = alpha * np.cos(beta) * \
                     np.sin(big_omega_dtc)
            r_obs = np.sqrt( ( x_obs*x_obs ) + ( y_obs*y_obs ) )
            return x_obs,y_obs,r_obs
            
        (self.x_obs, self.y_obs, self.r_obs) = \
            calc_obs_position(self.alpha.value, self.beta, self.big_omega_dtc)
            
        (self.x_obs_o, self.y_obs_o, r_obs_o) = \
            calc_obs_position(alpha_o.value, self.beta, big_omega_dtc_o)
           
        idx = np.where( abs(self.t-self.t_o) == abs(self.t-self.t_o).min() )
        self.x_obs_o = self.x_obs[idx]
        self.y_obs_o = self.y_obs[idx]
        self.it_o = idx
        
        if debug == True:
            print 'Obs positions: ',self.x_obs,self.y_obs
            print 'Obs position at t0: ',self.x_obs_o, self.y_obs_o, r_obs_o
            imid = int(float(len(self.t))/2.0)
            print 'Obs position in middle of lightcurve: ', \
                    self.x_obs[imid],self.y_obs[imid]
        
        # Normalizing the observer's position so that the lens position
        # closest approach is to the observer, offset to the current
        # uoffset value:
        if debug == True:        
            print 'X: ',self.alpha.value, (np.sin(self.beta)), \
                (np.cos(self.big_omega_dtc)), self.x_obs
            print 'Y: ',self.alpha.value, (np.cos(self.beta)), \
                (np.sin(self.big_omega_dtc)), self.y_obs
        
        #self.x_obs = self.x_obs - self.x_obs_o[0] + ( self.u_offset / np.sqrt(2.0) )
        #self.x_obs = self.x_obs - self.x_obs_o[0]
        #self.y_obs = self.y_obs - self.y_obs_o[0]
        #self.x_obs_o = self.x_obs_o[0] - self.x_obs_o[0]
        #self.x_obs_o = self.y_obs_o[0] - self.y_obs_o[0]
        
        if debug == True:
            print 'Obs position at t0: ',self.x_obs_o, self.y_obs_o, r_obs_o
            print 'Obs position at middle of x, y position arrays: ',\
                    self.x_obs[imid],self.y_obs[imid]
            idx = np.where( abs(self.t-self.t_o) == abs(self.t-self.t_o).min() )
            print 'Minimum in t-to at : ',idx, self.t[idx],self.t_o, \
                    self.x_obs[idx], self.y_obs[idx], (self.t-self.t_o)[idx]
            
    ###############################
    # POINT-SOURCE, POINT-LENS
    
    # Calculate source-lens relative motion for a generic PSPL event:
    def calc_source_lens_rel_motion(self, debug=False):
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
    	
    	# Calculate the time-dependent position of the lens relative 
    	# to the source in the lens plane in units of the Einstein radius:
        #x_incr = ( 1.0 + 1.0 )/len(self.t)
        #self.x_lens = np.arange(-1.0,1.0,x_incr)
        #self.y_lens = ( self.x_lens * np.tan(self.phi) ) + \
        #         ( self.u_min * np.cos(self.phi) )
        #self.r_lens = np.sqrt(self.x_lens*self.x_lens + \
        #         self.y_lens*self.y_lens)
        def calc_lens_position(u_min, phi,v_SL_dx,v_SL_dy,t,t_o):
            x_lens = ( -1.0 * u_min * np.sin(phi) ) + \
                        ( v_SL_dx * ( t - t_o ) )
            y_lens = ( u_min * np.cos(phi) ) + \
                        ( v_SL_dy * (t - t_o ) )
            r_lens = np.sqrt(x_lens*x_lens + y_lens*y_lens)
            return x_lens,y_lens,r_lens
            
        (self.x_lens,self.y_lens,self.r_lens) = calc_lens_position(self.u_min, \
                        self.phi,self.v_SL_dx,self.v_SL_dy,self.t,self.t_o)
                
        (self.x_lens_o,self.y_lens_o,r_lens_o) = calc_lens_position(self.u_min, \
                        self.phi,self.v_SL_dx,self.v_SL_dy,self.t_o,self.t_o)
        
        if debug == True:
            print 'X_Lens: ',self.x_lens
            print 'Y_Lens: ',self.y_lens
            print 'X,Y Lens at t0: ',self.x_lens_o, self.y_lens_o
            
        
    ################################
    # EVENT TIMELINE
    def gen_event_timeline(self,cadence=None,lc_length=None, \
                            force_t0_obs=False, debug=False):
        """Method to generate the event timeline, an array of timestamps 
        for datapoints on the event lightcurve, 
        spanning the range t_o +/- (t_E/2).
        Optional argument:
            cadence float [mins]
        """
        
        # Generate an array of incremental timestamps throughout the event 
        # in JD:
        def gen_timestamps(t_start, t_end, t_incr, target_position, debug):
            ts = []
            jd = []
            t = t_start
            while t <= t_end:
                t = t + t_incr
                jd.append( t.jd )
                hjd = self.jd_to_hjd(t, target_position, debug=debug)
                ts.append(hjd)
            return ts, jd
            
        # Compute the start and end times of the event, note these are in 
        # UTC format:
        # Plots extend to +/- 2 * t_E by default, or this can be fixed
        if lc_length == None:        
            half_lc = TimeDelta( ( self.t_E.value * 2.0 ), format='sec' )
        else:
            half_lc = TimeDelta( ( (lc_length*60.0*60.0*24.0) / 2.0 ), 
                                        format='sec' )
        
        # Convert coordinates for future use:
        (ra_rads, dec_rads) = utilities.sex2rads(self.RA, self.Dec)
        target_position = slalib.sla_dcs2c( ra_rads, dec_rads )
        
        if cadence == None:
            t_incr = TimeDelta( ( self.t_E.value / 5000.0 ) , format='sec')
        else:
            t_incr = TimeDelta( ( cadence * 60.0 ), format='sec')
            
        # Store the Time version of t_o in JD for future reference:
        if self.t_o_jd == None:
            self.t_o_jd = self.t_o
        
        if force_t0_obs == False:
            event_start = self.t_o_jd - half_lc
            event_end = self.t_o_jd + half_lc
                
            # Generate an array of incremental timestamps throughout the event 
            # in JD:
            (ts, jd) = gen_timestamps(event_start, event_end, t_incr, \
                            target_position, debug)

        else:
            
            # Build the lightcurve in two halves, centered at t0:
            event_start = self.t_o_jd
            event_end = self.t_o_jd + half_lc
            
            (ts1, jd1) = gen_timestamps(event_start, event_end, t_incr, \
                            target_position, debug)
            event_start = self.t_o_jd - half_lc
            event_end = self.t_o_jd
            
            (ts2, jd2) = gen_timestamps(event_start, event_end, t_incr, \
                            target_position, debug)
            ts = ts1 + ts2
            jd = jd1 + jd2
            
        self.t = np.array(ts)
        self.ts_jd = np.array(jd)
        #print 'TO (JD) = ',self.t_o.jd
        self.t_o = self.jd_to_hjd(self.t_o_jd, target_position, debug=debug)
        #print 'TO (HJD) = ',self.t_o
    
        #print 'T_P (JD) = ',self.t_p.jd
        if self.t_p_jd == None:
            self.t_p_jd = self.t_p
            
        self.t_p = self.jd_to_hjd(self.t_p_jd, target_position, debug=debug)
        #print 'T_P (HJD) = ',self.t_p
            
    
    def jd_to_hjd(self, t, target_position, debug=False):
        """Calculate the HJD timestamp corresponding to the JD given for the
        current event"""
        
        if debug == True:
            print 'TIME JD: ',t, t.jd
            
        # Calculate the MJD (UTC) timestamp:
        mjd_utc = t.jd - 2400000.5
        if debug == True:
            print 'TIME MJD_UTC: ',mjd_utc
        
        # Correct the MJD to TT:
        mjd_tt = timesubs.mjd_utc2mjd_tt(mjd_utc)
        if debug == True:
            print 'TIME MJD_TT: ',mjd_tt, t.tt.jd
        
        # Calculate Earth's position and velocity, both heliocentric
        # and barycentric for this date
        (earth_helio_position, vh, pb, vb) = slalib.sla_epv( mjd_tt )
        if debug == True:
            print 'Earth Cartesian position: ',earth_helio_position
            print 'Target cartesian position: ', target_position
        
        # Calculate the light travel time delay from the target to 
        # the Sun:
        dv = slalib.sla_dvdv( earth_helio_position, target_position )            
        tcorr = dv * ( constants.au.value / constants.c.value )
        
        if debug == True:
            print 'TIME tcorr: ', tcorr, 's', (tcorr/60.0),'mins'
        
        # Calculating the HJD:
        hjd = mjd_tt + tcorr/86400.0 + 2400000.5
        if debug == True:
            print 'TIME HJD: ',hjd,'\n'

        return hjd
    
    
    def calc_static_observer_impact_param(self):
        """Method to plot (and calculate if necessary) the point-source, 
        	point-lens lightcurve to file."""
	
    	# Calculate the impact parameter as a function of time
    	# throughout the event, assuming uniform-speed, straight-line 
    	# relative motion of lens and source:
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
    def calc_parallax_impact_param(self, set_uo=False):
        """Method to calculate a point-source, point-lens lightcurve given 
        the lens parameters, including the parallax due to Earths orbital 
        motion.
        Formulation derived from Buchalter, A., Kamionkowski, M., 1997, 
        ApJ, 483, 782.
        """
        dbg = False
        method = 'offset'
        
    	# Pre-calculate some repeated factors for convenience and speed:
        alpha_sq = self.alpha.value * self.alpha.value
	
    	# Calculate impact parameter as a function of time:
        # The projected separation of lens-source from the observer's 
        # perspective. 
        if method=='first_principles':
            term1 =  self.u_min*self.u_min
            term2 = self.omega_sq_dt_sq
            term3 = alpha_sq * ( np.sin( self.big_omega_dtc ) )**2 
            term4 = alpha_sq * ( ( np.sin( self.beta.value ) )**2 ) * \
                        ( ( np.cos(self.big_omega_dtc) )**2 ) 
            term5 = -2.0 * self.alpha.value * np.sin(self.big_omega_dtc) * \
                             ( self.omega_dt * np.sin(self.phi) + \
                                 self.u_min * np.cos(self.phi) ) 
            term6 = 2.0 * self.alpha.value * np.sin(self.beta.value) * \
                         np.cos(self.big_omega_dtc) * \
                         ( self.u_min * np.sin(self.phi) - self.omega_dt * \
                         np.cos(self.phi) )
        
            usq_t = term1 + term2 + term3 + term4 + term5 + term6
        	
        	# Print these terms at t_o?  Or find when they are all minimised? 
            # Term 3 largest at t=t_o
            if dbg == True:
                print 'u(t)^2 terms: ',term1, \
        	      	      	term2[0], \
        			term3[0], \
        			term4[0], \
        			term5[0],\
        			term6[0], \
                       alpha_sq, ( np.sin( self.beta.value ) )**2, \
                       ( np.cos(self.big_omega_dtc) )**2
        
        else:
            
            # Offset the trajectory of the lens so that it passes
            # close to the position of the source at t0. 
            (x_lens, y_lens) = self.offset_lens()
            
            dx = x_lens - self.x_obs
            dy = y_lens - self.y_obs
            usq_t = dx*dx + dy*dy
            
            if dbg == True:
                idx = self.it_o
                print 'LENS at t[idx]: ',self.x_lens[self.it_o], self.y_lens[self.it_o]
                print 'OBS at t[idx]: ',self.x_obs[self.it_o], self.y_obs[self.it_o]
                print 'Delta pos: ',dx[idx],dy[idx]
                u = np.sqrt( dx[idx]*dx[idx] + dy[idx]*dy[idx] )
                print 'u [idx]: ',u
            
    	#dx = self.x_lens - self.x_obs
    	#dy = self.y_lens - self.y_obs
    	
    	#usq_t = np.sqrt( dx*dx + dy*dy )
    	
    	#print usq_t
    	
        self.u_t = np.sqrt(usq_t)
        if set_uo == True:
            self.u_o = self.u_t.min()
        
        #print 'u(t) min: ',self.u_t.min()

    def offset_lens(self):
        uoffset = self.u_offset / np.sqrt(2.0)
        x_lens = self.x_lens + self.x_obs_o + uoffset
        y_lens = self.y_lens + self.y_obs_o + uoffset
        return x_lens, y_lens
        
    ##################################
    # OBSERVATION SIMULATION
    def calc_mag( self, model='pspl' ):
        """Function to calculate the magnitude of an event as a function of
        time given its baseline magnitude and magnification"""
        
        if model == 'pspl': 
            A_t = self.A_t_pspl
        if model == 'fspl': 
            A_t = self.A_t_fspl
        if model == 'pspl_parallax': 
            A_t = self.A_t_pspl_parallax
        if model == 'pspl_parallax_satelite': 
            A_t = self.A_t_pspl_parallax_satellite
        if model == 'fspl_parallax': 
            A_t = self.A_t_fspl_parallax
        if model == 'fspl_parallax_satelite': 
            A_t = self.A_t_fspl_parallax_satellite
        
        mag = self.mag_base - 2.5 * np.log10( A_t )
        return mag

    def data_within_interval(self,obs_start,obs_end,A_t,exp_time,phot_precision):
        
        idx1 = np.where( self.t >= obs_start )
        idx2 = np.where( self.t <= obs_end )
        idx = np.intersect1d(idx1,idx2)
        
        ts = self.t[idx].tolist()
        mag = []
        merr = []
        for i in idx:
            mag_mean = self.mag_base - 2.5 * np.log10( A_t[i] ) 
            merr_est = self.sim_mag_error( exp_time, mag_mean, phot_precision ) 
            mag.append( np.random.normal( mag_mean, merr_est ) )
            merr.append( self.sim_mag_error( exp_time, mag[-1], phot_precision ) )
        
        return ts,mag,merr

    def select_data_at_intervals(self,A_t,obs_start,obs_end,\
                            interval, window, exp_time, phot_precision, \
                            force_t0_obs=False):
        """Method to select datapoints from the existing model at intervals 
        according to the observational parameters"""
        ts = []
        mag = []
        merr = []
        first_group = True
        
        def test_continue(interval, obs_start, obs_end, interval_start):
            if interval < 0:
                if interval_start > obs_start:
                    return True
                else:
                    return False
            
            else:
                if interval_start < obs_end:
                    return True
                else:
                    return False
                
        
        # If forcing observations around t0, one interval is split over t0, 
        # so the first group is a half-group:
        if force_t0_obs == True:
            if interval > 0:
                interval_start = obs_start
                interval_end = obs_start + (window/2.0)
            else:
                interval_start = obs_end - (window/2.0)
                interval_end = obs_end
            (t,m,me) = self.data_within_interval(interval_start,interval_end,\
                                        A_t,exp_time,phot_precision)
                                            
            ts = ts + t
            mag = mag + m
            merr = merr + me
            
            interval_start = interval_start + interval
            interval_end = interval_start + window
        else:
            interval_start = obs_start
            interval_end = interval_start + window
            
        cont = test_continue(interval, obs_start, obs_end, interval_start)
        while cont == True:
            (t,m,me) = self.data_within_interval(interval_start,interval_end,A_t,\
                                exp_time, phot_precision)
            ts = ts + t
            mag = mag + m
            merr = merr + me
            interval_start = interval_start + interval
            interval_end = interval_start + window
            cont = test_continue(interval, obs_start, obs_end, interval_start)
            
        return ts, mag, merr

    def simulate_data_points(self,model='pspl', phot_precision='1m', \
            window=None, interval=None,force_t0_obs=False,log=None):
        """Method to extract a series of datapoints from the current model
        with a given photometric cadence and precision."""

        # Exposure time set to 200s for consistency, observing cadence
        # expected to be 220s:
        exp_time = 200.0
        if log != None:
            log.info('Simulating data points with window='+str(window)+\
                    ' and interval='+str(interval))
        if model == 'pspl': 
            A_t = self.A_t_pspl
        if model == 'fspl': 
            A_t = self.A_t_fspl
        if model == 'pspl_parallax': 
            A_t = self.A_t_pspl_parallax
        if model == 'pspl_parallax_satelite': 
            A_t = self.A_t_pspl_parallax_satellite
        
        if window != None:
            window = window / 24.0
            interval = interval / 24.0
        
        # Continuous sampling:
        if window == None:
            obs_start = self.t[0]
            obs_end = self.t[-1]
            (ts,mag,merr) = self.data_within_interval(obs_start,obs_end,A_t,\
                                exp_time, phot_precision)
                                
        # Observations at periodic intervals:
        else:
            if force_t0_obs == False:
                obs_start = self.t[0]
                obs_end = self.t[-1]
                (ts, mag, merr) = self.select_data_at_intervals(A_t,obs_start,\
                                obs_end,interval, window, exp_time, phot_precision,\
                                force_t0_obs=force_t0_obs)
            else:
                obs_start = self.t[0]
                obs_end = self.t_o
                (t1, m1, me1) = self.select_data_at_intervals(A_t,obs_start,\
                                obs_end,-1.0*interval, window, exp_time, phot_precision,\
                                force_t0_obs=force_t0_obs)
                
                obs_start = self.t_o
                obs_end = self.t[-1]
                (t2, m2, me2) = self.select_data_at_intervals(A_t,obs_start,\
                                obs_end,interval, window, exp_time, phot_precision,\
                                force_t0_obs=force_t0_obs)
                ts = t1 + t2
                mag = m1 + m2
                merr = me1 + me2
                print len(t1), len(t2), len(ts), interval
            
        self.ts = np.array( ts )
        self.mag = np.array( mag )
        self.merr = np.array( merr )
        
    def sim_mag_error( self, exp_time, mag, precision_model='1m', debug=False):
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
    
    ##################################
    # FILE OUTPUT 
    def output_data( self, file_name ):
        """Function to output a lightcurve file containing the data
        from the current model"""

        if path.isfile( file_name ) == True:
            bkup = file_name.replace('.dat', '_old.dat')
            move(file_name, bkup)
        
        te = self.t_E / ( 24.0 * 60.0 * 60.0 )
        fileobj = open( file_name, 'w' )
        fileobj.write( '# u_o = ' + str(self.u_o) + '\n' )
        fileobj.write( '# u_min = ' + str(self.u_min) + '\n' )
        fileobj.write( '# t_E = ' + str(te) + 'd\n' )
        fileobj.write( '# phi = ' + str(self.phi) + 'rads\n' )
        fileobj.write( '# rho = ' + str(self.rho) + '\n' )
        fileobj.write( '# Baseline mag = ' + str(self.mag_base) + '\n' )
        fileobj.write( '# HJD       mag         merr\n')
        n_data = len( self.mag ) - 1
        for i in range(0,n_data,1):
            fileobj.write( str(self.ts[i]) + '  ' + str(self.mag[i]) + \
                     '  ' + str(self.merr[i]) + '\n' )
        fileobj.close()
    
    def output_model( self, file_name, model='pspl' ):
        """Function to output a lightcurve file containing the data
        from the current model"""
        
        if model == 'pspl': 
            A_t = self.A_t_pspl
        if model == 'fspl': 
            A_t = self.A_t_fspl
        if model == 'pspl_parallax': 
            A_t = self.A_t_pspl_parallax
        if model == 'pspl_parallax_satelite': 
            A_t = self.A_t_pspl_parallax_satellite
        
        te = self.t_E / ( 24.0 * 60.0 * 60.0 )
        fileobj = open( file_name, 'w' )
        fileobj.write( '# u_o = ' + str(self.u_o) + '\n' )
        fileobj.write( '# u_min = ' + str(self.u_min) + '\n' )
        fileobj.write( '# t_E = ' + str(te) + 'd\n' )
        fileobj.write( '# phi = ' + str(self.phi) + 'rads\n' )
        fileobj.write( '# rho = ' + str(self.rho) + '\n' )
        fileobj.write( '# Baseline mag = ' + str(self.mag_base) + '\n' )
        fileobj.write( '# HJD       mag\n')
        n_data = len( A_t ) - 1
        for i in range(0,n_data,1):
            mag = self.mag_base - 2.5 * np.log10( A_t[i] )
            fileobj.write( str(self.t[i]) + '  ' + str(mag) + '\n' )
        fileobj.close()
        
    ##################################
    # PLOT LENS PLANE MOTION
    def plot_lens_plane_motion(self,params=None):
        """Method to plot the motion of the lens as a function of time in 
        the lens plane geometry"""
	
    	# Function to describe the circle representing the Einstein ring:
        def circle(r):
            phi=np.arange(0,6.28,0.01)	    # 0 -> 2.pi
            xcir = r*np.cos(phi)
            ycir = r*np.sin(phi)
            return xcir, ycir
    	
        # Set up figure:				   
        fig = plt.figure(2,(16,10)) 
        ax = fig.add_axes([0.1, 0.45, 0.4, 0.4])   #  [left, bottom, width, height]
	
    	  # Plot the trajectory of the lens relative to the Sun-source 
        # line during the event:
        (x_lens,y_lens) = self.offset_lens()
        plt.plot(x_lens,y_lens,'r-')
        plt.plot(self.x_lens_o,self.y_lens_o,'r.')
	
         # Plot the trajectory of the observer in the lens plane during 
        # the event:
        plt.plot(self.x_obs,self.y_obs,'b-')
        plt.plot(self.x_obs_o,self.y_obs_o,'b.')
        
        plt.xlabel('x [R$_{E}$]', fontsize=18)
        plt.ylabel('y [R$_{E}$]', fontsize=18)
        (x_E, y_E) = circle(1.0)      # Plot in units of the Einstein radius
        plt.plot(x_E+self.x_obs_o,y_E+self.y_obs_o,'k-.')
        #plt.axis([-2.0,2.0,-2.0,2.0])
        plt.box('off')
        plt.grid(True)
        plt.axhline(y=0, color='black')
        plt.axvline(x=0, color='black')
    
        # Plot a zoom-in of the trajectory of the lens relative to the Sun-source 
        # line during the event:
        ax = fig.add_axes([0.55, 0.45, 0.4, 0.4])   #  [left, bottom, width, height]
        plt.plot(x_lens,y_lens,'r-')
        plt.plot(self.x_lens_o,self.y_lens_o,'r.')
	
         # Plot the trajectory of the observer in the lens plane during 
        # the event:
        plt.plot(self.x_obs,self.y_obs,'b.')
        plt.plot(self.x_obs_o,self.y_obs_o,'b.')
        
        plt.xlabel('x [R$_{E}$]', fontsize=18)
        plt.ylabel('y [R$_{E}$]', fontsize=18)
        (x_E, y_E) = circle(1.0)      # Plot in units of the Einstein radius
        plt.plot(x_E+self.x_obs_o,y_E+self.y_obs_o,'k-.')
        xmin = self.x_obs_o - 2.0
        xmax = self.x_obs_o + 2.0
        ymin = self.y_obs_o - 2.0
        ymax = self.y_obs_o + 2.0
        plt.axis([xmin,xmax,ymin,ymax])
        plt.box('off')
        plt.grid(True)
        plt.axhline(y=0, color='black')
        plt.axvline(x=0, color='black')
    
    
        # Plot lens and observer separation as functions of time:
        ax = fig.add_axes([0.1, 0.1, 0.3, 0.2])   #  [left, bottom, width, height]
        plt.plot(self.t-2450000.0,(x_lens),'r-')
        plt.plot(self.t-2450000.0,(self.x_obs),'b-')
        (xmin,xmax,ymin,ymax) = plt.axis()
        plt.plot(np.array( [self.t_o-2450000.0]*2 ), \
                         np.array([ymin,ymax]),'r-.')
        plt.xlabel('HJD-2450000.0', fontsize=18) 				     
        plt.ylabel('x$_{lens}$ [R$_{E}$]', fontsize=18) 
        ax.yaxis.grid() #vertical lines
        ax.xaxis.grid() #horizontal lines		    
        
        ax = fig.add_axes([0.55, 0.1, 0.3, 0.2])   #  [left, bottom, width, height]
        plt.plot(self.t-2450000.0,(y_lens),'r-')
        plt.plot(self.t-2450000.0,(self.y_obs),'b-')
        (xmin,xmax,ymin,ymax) = plt.axis()
        plt.plot(np.array( [self.t_o-2450000.0]*2 ), \
                         np.array([ymin,ymax]),'r-.')
        plt.xlabel('HJD-2450000.0', fontsize=18) 				     
        plt.ylabel('y$_{lens}$ [R$_{E}$]', fontsize=18) 
        ax.yaxis.grid() #vertical lines
        ax.xaxis.grid() #horizontal lines		    
        
        if params != None:
            plot_file = path.join( params['output_path'], 'lens_plane_motion.png')
        else:
            plot_file = 'lens_plane_motion.png'
        plt.savefig(plot_file)
        plt.close(2)

    def plot_lightcurve(self, plot_file, model_list=['pspl']):
        """Method to plot a lightcurve of the event"""
        
            
        # Set up figure:				   
        fig = plt.figure(2,(10,12)) 
        ax = fig.add_axes([0.1, 0.3, 0.85, 0.65])   #  [left, bottom, width, height]
        
        plt_style = { 'pspl': 'r-', 'fspl': 'b-', \
                        'pspl_parallax': 'm-', 'pspl_parallax_satelite': 'k-' }        
        
        # Plot the list of models
        for model in model_list:
            mag = self.get_mag_lc( self, model )
            plt.plot( self.t, mag, plt_style[model] )
        
        
        plt.xlabel('HJD-2450000.0', fontsize=18) 				     
        plt.ylabel('Mag', fontsize=18) 
        ax.yaxis.grid() #vertical lines
        ax.xaxis.grid() #horizontal lines		    
    
        plt.savefig(plot_file)
        plt.close(2)
    
    def get_mag_lc( self, model ):

        ZP = 25.0
        available_models = ['pspl', 'pspl_parallax', 'pspl_parallax_satelite', \
                            'fspl', 'fspl_parallax', 'fspl_parallax_satelite']
        if model not in available_models:
            print 'Error: Unrecognised model requested ' + model
            print 'Need one of ', available_models
            exit()
        
        A_t = getattr(self, 'A_t_'+model)
        
        if A_t == None:
            print 'Error: magnification for ' + model + ' not yet calculated'
            exit()
            
        mag = ZP - 2.5 * np.log10( A_t )
        
        return mag
        
    def plot_diff_models(self, plot_file, model_list=['pspl']):
        """Method to plot lightcurve of the event including different effects, 
        and to plot a panel differencing the models"""
        
        legend_text = '$u_{O}$='+str(self.u_o)	
        plt_style = { 'pspl': 'r-', 'fspl': 'b--', \
                        'pspl_parallax': 'm-', 'pspl_parallax_satelite': 'k-', \
                        'fspl_parallax': 'm--', 'fspl_parallax_satelite': 'k--'}
                        
        fig = plt.figure(1,(12,12)) 
    
        # Top panel						   
        ax = fig.add_axes([0.15, 0.5, 0.775, 0.45])   #  [left, bottom, width, height]	   
    
        # Plot PSPL event seen from stationary Earth:
        lc_mag = {}
        for model in model_list:
            lc_mag[model] = self.get_mag_lc( model )
            plt.plot(self.t-2450000.0,lc_mag[model],plt_style[model],label=model)
            
        dt = 2.0
        (xmin,xmax,ymin,ymax) = plt.axis()
        #xmin = self.t_o.jd - 2450000.0
        #xmax = self.t_o.jd - 2450000.0
        plt.xlabel('HJD-2450000.0', fontsize=18)						   
        plt.ylabel('Mag', fontsize=18)		
        plt.axis([xmin,xmax,ymax,ymin])
        plt.legend(loc='upper right',frameon=False) 
        ax.yaxis.grid() #vertical lines
        ax.xaxis.grid() #horizontal lines
        
        # Add plot of residuals:			   
        ax = fig.add_axes([0.15, 0.125, 0.775, 0.3])   #  [left, bottom, width, height]
        diff = lc_mag[model_list[0]] - lc_mag[model_list[1]]
        label = model_list[0]+' - '+model_list[1]
        plt.plot(self.t-2450000.0,diff,'r-',label=label)					   
        plt.ylabel(label, fontsize=18)	
        
        # Complete second plot and save to file:
        (xmin,xmax,ymin,ymax) = plt.axis()
        plt.axis([xmin,xmax,ymax,ymin])	
        plt.xlabel('HJD-2450000.0', fontsize=18)				   
        plt.legend(loc='upper right',frameon=False) 					   
        ax.yaxis.grid() #vertical lines
        ax.xaxis.grid() #horizontal lines						   
        plt.savefig(plot_file)