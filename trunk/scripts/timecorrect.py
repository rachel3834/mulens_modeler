#!/usr/bin/env python -Wignore::DeprecationWarning                                                                      
import pyfits
import slalib as S
import warnings
from math import radians,degrees, cos, sin
from correct_subs import *
from numpy import array
from sys import argv, exit

debug = False
#  Astronomical unit in km and light-seconds(TDB), sun's radius in AU,
#  Shapiro constant (2*G*M_sun/c^3 in TDB seconds)
aukm = 149597870.691
ausec = 499.0047838061
sunrau = 0.00465247263738,
shapiro = 9.8509818966176734e-06

# Declare version and corrections performed
catalog_version = 'LCOGT-FITS-SRCCAT-1.4.3'
version = 'timecorrect.py $Revision: 7319 $'
corrections = 'ROMER,CLOCK,EINSTEIN'
#corrections = 'ROMER,CLOCK,EINSTEIN,SHAPIRO' # Not yet...

if len(argv) == 2:
    fitsimfile = argv[1]
else:
#   fitsimfile = 'bpl1m001-kb72-20110930-1026-e00_bp_ff.fits'
#   fitscatfile = 'catalogue_f1026.fits.1.0'
    fitsimfile = 'bpl1m001-kb72-20111028-0080-e90.fits'

# Name of the source catalog file to read in based on the image filename we were
# passed on the cmdline.
if fitsimfile.rfind('e90') > 1:
    fitscatfile = fitsimfile.replace('e90', 'e90_cat')
elif fitsimfile.rfind('s90') > 1:
    fitscatfile = fitsimfile.replace('s90', 's90_cat')
elif fitsimfile.rfind('e10') > 1:
    fitscatfile = fitsimfile.replace('e10', 'e10_cat')
elif fitsimfile.rfind('s10') > 1:
    fitscatfile = fitsimfile.replace('s10', 's10_cat')
else:
    print "Unrecognized filetype"
    exit(-1)



# Open source detection catalog
hdulist = pyfits.open(fitscatfile)

if debug: hdulist.info()
prihdr = hdulist[0].header

tbdata = hdulist[1].data
cols = hdulist[1].columns
extnhdr = hdulist[1].header
  
# Get needed header items
hdrs_needed = obtain_hdrpars(prihdr, debug)

# Compute TT midpoint of the exposure (in MJD)
hdrs_needed = compute_tt_mid(hdrs_needed, debug)

# Compute UT1
hdrs_needed = compute_ut1(hdrs_needed)
if debug: print hdrs_needed

# Need *geocentric* position of the observatory for the TDB clock 
# correction so we have to convert from the geodetic position given
# by a GPS.
(r,z) = sla_geoc_iers2003(hdrs_needed['OBS_LAT'], hdrs_needed['OBS_HGT'])

# Calculate relativistic clock correction (TDB-TT) and put into dictionary
tdb_tt = S.sla_rcc(hdrs_needed['MJD_TT_mid'], hdrs_needed['UT1'], 
    -hdrs_needed['OBS_ELONG'], r, z)
tdb_utc = hdrs_needed['TT-UTC'] + tdb_tt
if debug: print 'TDB-TT(RCC)=', tdb_tt
if debug: print 'TDB-UTC    =', tdb_utc
new_dict = { 'MJD_TDB_mid' : hdrs_needed['MJD_TT_mid'] + (tdb_tt/86400.0) }
hdrs_needed.update(new_dict)

#  Compute star-independent parameters in preparation for
#  conversions between mean place and geocentric apparent place.
amprms = S.sla_mappa(2000.0, hdrs_needed['MJD_TT_mid'])

# Reform into precession/nutation matrix
irmat = amprms[12:21].reshape(3,3)

# Earth's barycentric position (in AU) from sla_evp (max error 6900 km ~=23ms)
#pos_earth_bary = array([amprms[1],amprms[2],amprms[3]])

# Earth's barycentric position (in AU) from sla_epv (max error 13.4 km ~=45us)
(epos_hel, evel_hel, pos_earth_bary, evel_bary) = S.sla_epv(hdrs_needed['MJD_TDB_mid'])

# Position of observer
hdrs_needed = compute_local_st(hdrs_needed, debug)
if debug: print 'LAST=', hdrs_needed['LST_App_mid']
pv = array([ r*cos(hdrs_needed['LST_App_mid']), 
    	     r*sin(hdrs_needed['LST_App_mid']),
	     z])
pvnew = S.sla_dmxv(irmat, pv)
if debug: print 'Geocenter->observer (mean) (AU)=',pvnew/aukm
# Compute barycenter->topocentric observer position vector
bary_topo = pos_earth_bary+(pvnew/aukm)

# Loop over all stars in table computing apparent RA, Dec 
airmass_list = []
bjd_list = []
for i in tbdata:
    if debug == False:
	ra_degs = i['ALPHA_J2000']
	dec_degs =  i['DELTA_J2000']
    else:
# Hardwired RA, Dec for testing 
	ra_degs = (20.0+(24.0/60.0)+(29.740/3600.0))*15.0
	dec_degs = 16.0 + (45.0/60.0) + (43.79/3600.0) 
    ra_rads = radians(ra_degs)
    dec_rads = radians(dec_degs)

# Form Cartesian vector from mean J2000.0 spherical RA, Dec
#
    v = S.sla_dcs2c(ra_rads,dec_rads)

# Light travel time (Romer delay)
#
    ltt = S.sla_dvdv(bary_topo, v) * ausec
    if debug: print 'LTT(secs)=', ltt
    delta_t = ltt + tdb_utc # +shapiro_corr
    if debug: print 'BJD_TDB-JD_UTC= ', delta_t
    bjd_list.append(delta_t)
    

bjds = array(bjd_list)
