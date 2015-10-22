import pyfits
import slalib as S
from math import degrees, radians, sin, cos, sqrt, modf

class FITSHdrException(Exception):
    '''Raised when a required FITS header keyword is missing'''
    
    def __init__(self, keyword):
    	self.keyword = keyword
    
    def __str__(self):
    	return "Required keyword '" + self.keyword + "' missing"

def obtain_hdrpars(hdr,debug=False):
    if debug == -1:
      	print "DANGER! Using assumed header values"
    	headers = {
	    'MJD_UTC' : 55834.1117707361,
	    'DATEOBS_UTC' : '2011-09-30T02:40:56.943',
	    'EXPSTART_UTC' : '02:40:56.943',
	    'EXPTIME' : 30.0,
	    'OBS_ELONG' : radians(-119.863103),
	    'OBS_LAT' : radians(34.433161),
	    'OBS_HGT' : 7.0,
	    'UT1-UTC' : -0.3204774, #0.0,
	    'XP' : radians(0.179833/3600.0), # 0.0,
	    'YP' : radians(0.378683/3600.0), #0.0,
	    'HUMIDITY' : 79.8/100.0,
	    'TEMP_K' : 17.1+273.15,
	    'PRES_MB' : 1013.1169501,
	    'WAVELENGTH' : 550.0/1000.0,
	    'LAPSERATE' : 0.0065
	}
    else:
# Read from passed hdr    
    	headers = {}
# Dictionary of concepts and header keywords to look for, with default values
# to use if not found. -1 means default can't be assumed and we should throw
# an exception if we don't find the keyword in the header
    	hdr_keys = { 	    
	    'MJD_UTC' : ('MJD-OBS', -1),
	    'DATEOBS_UTC' : ('DATE-OBS', -1),
	    'EXPTIME' : ('EXPTIME', -1),
	    'EXPSTART_UTC' : ('UTSTART', -1),
            'DAYOBS_UTC' : ('DAY-OBS', -1),
	    'OBS_ELONG' : ('LONGITUD', -1),
	    'OBS_LAT' : ('LATITUDE', -1),
	    'OBS_HGT' : ('HEIGHT', -1),
	    'UT1-UTC' : ('UT1-UTC', 0.0),
	    'XP' :  	('POLARMOX', 0.0),
	    'YP' :  	('POLARMOY', 0.0),
	    'EOP_SRC' : ('EOPSRC', 'timecorrect.py'),
# Default values for humidity and temp. need to be in % and degrees C as in 
# header, as they will be written and then converted to correct units for SLALIB 
	    'HUMIDITY' : ('WMSHUMID', 55.0), 
	    'TEMP_K' :  ('WMSTEMP', 15.0),
	    'PRES_MB' : ('WMSPRES', 1013.25),
	    'WAVELENGTH' : ('WAVELEN', 0.55),
	    'LAPSERATE' : ('LAPSRATE', 0.0065)
	}

# Search through header for all the values we need.
	for i in hdr_keys.keys():
	    fits_keyword = hdr_keys[i][0]
	    default = hdr_keys[i][1]
    	    newitem = {}
	    if fits_keyword in hdr:
# Found, extract value
    	    	value = hdr[fits_keyword]
		if value == 'UNKNOWN':
		    if debug: print 'UNKNOWN value found for ' + fits_keyword
    	    	    if default == -1:
# We found an UNKNOWN for the value and we *really* need a proper value, throw exception
		      	raise FITSHdrException(fits_keyword)
		    else:
# Store original value
                        origitem = { 'original_'+fits_keyword : value }
                        headers.update(origitem)
		      	value = default
      	      	      	hdr.update(fits_keyword, default, 'Assumed value rewritten by timecorrect.py')
# After we have written it to the header, convert the value to the right units
                        value = convert_values(fits_keyword, value)
      	      	else:
		    value = convert_values(fits_keyword, value)
    	    	newitem = { i : value }
	    else:
# Not found, check if we need it
    	    	if default == -1:
# We, like, *really* need it man, throw exception
		    raise FITSHdrException(fits_keyword)
		else:
# Not essential, use default and write keyword and assumed value to FITS header
		    newitem = { i : default }
		    hdr.update(fits_keyword, default, 'Assumed value written by timecorrect.py')
# Append to dictionary
	    headers.update(newitem)

# Verify headers for MJD-OBS != DATE-OBS and DATE/MJD-OBS in TDB problems
    headers = verify_mjd_utc(headers,debug)

    return headers

def convert_values(keyword, value, debug=False):
    '''Handles unit conversion for FITS header keywords to the values expected
    by timecorrect/SLALIB.
    Unrecognized keywords are passed through with value unchanged'''
    
    newvalue = value
    if keyword == 'LATITUDE' or keyword == 'LONGITUD':
# Convert from degrees to radians    
    	newvalue = radians(value)
    elif  keyword == 'WMSTEMP':
# Convert from deg C to deg K
    	newvalue = value + 273.15
    elif keyword == 'POLARMOX' or keyword == 'POLARMOY':
# Convert from arcseconds to radians
    	newvalue = radians(value/3600.0)
    elif keyword == 'WMSHUMID':
# Convert from percentage to 0.0..1.0 fraction
    	newvalue = value / 100.0

    if debug: print "Keyword, Old, new values=",keyword, ' ', value, ' ', newvalue
    return newvalue

def copy_new_headers(imhdr, cathdr, debug=False):
    '''Copies over header items from <imhdr> to <cathdr>, modifying <cathdr> 
    in place.
    Keywords in the skip_headers[] list are not copied from <imhdr>; those
    in the del_headers[] list are deleted from <cathdr>'''

# Keywords in the Image header to skip
    skip_headers = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'EXTEND',
    	    	    'BSCALE', 'BZERO', 'HDUCLAS1', 'HDUCLAS2', 'HDSTYPE', 'LBOUND1', 'LBOUND2',
		    'CTYPE1', 'CTYPE2', 'CD1_1', 'CD1_2','CD2_1', 'CD2_2',  
		    'PC1_1', 'PC1_2','PC2_1', 'PC2_2',
		    'QV1_0', 'QV1_1', 'QV1_7', 'QV1_9', 'QV1_17', 'QV1_21', 'QV1_31', 'QV1_37',
		    'QV2_0', 'QV2_1', 'QV2_7', 'QV2_9', 'QV2_17', 'QV2_21', 'QV2_31', 'QV2_37',
		     ]
# Keywords in the Catalog header to delete
    del_headers = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2',
    	    	   'CROTA1', 'CROTA2', 'EPOCH',]
 
    cathdr.add_comment('The following are a copy of the image headers', after='extend')
    startptr = cathdr.ascardlist().keys().index('ORIGIN')
    if debug: print startptr
    for item in imhdr.ascardlist():
    	if item.key in skip_headers:
	    action =  "Skipping"
	elif item.key in del_headers:
	    action = "Deleting"
	else:
	    action =  "To be copied"
    	    if item.key == 'COMMENT':
	    	cathdr.add_comment(item.value, after=startptr)
	    elif item.key == 'HISTORY':
# Don't copy HISTORY keywords right now
		startptr -= 1
#	    	cathdr.add_history(item.value, after=startptr)
	    elif item.key == '':
# Don't copy blank keywords right now
		startptr -= 1
#	    	cathdr.add_blank(after=startptr)
    	    elif item.key in cathdr:
# Keyword already exists, update in place
	    	cathdr.update(item.key, item.value, item.comment)
		startptr -= 1
	    else:
# Keyword does not exist, add it
	    	cathdr.update(item.key, item.value, item.comment, after=startptr)
 	    startptr += 1
   	if debug: print "%8s %s:-> %s" % ( item.key, item.value, action )

# Remove un-needed header keywords
    for delitem in del_headers:    
	if delitem in cathdr:
	    if debug: print "Deleting " + delitem
	    del cathdr[delitem]

    cathdr.add_comment('Catalog extraction headers:', before='fitsfile')

    return cathdr

def show_header(header):
    '''Basic pretty printer of a passed pyfits Header class. Prints
    keyword, value and comment, one per line.'''
     
    for item in header.ascardlist():
    	print "%-8s= %20s %s" % ( item.key, item.value, item.comment)
    
    return

def safe_add_keyword(header, keyword, value, comment, **position):
    '''Safely adds a new keyword to the header'''
    return
    
def add_missing_WCS(header,debug=False):
    '''Adds missing WCS keywords back to the header'''
    
    CD_string = 'WCS CD transformation matrix'
    QV_string = 'Projection parameter (nominal WCS)'
    new_WCS_keys = [
        ['CD1_2',  0.0, CD_string], 
        ['CD2_1',  0.0, CD_string], 
        ['QV1_0',  0.0, QV_string], 
        ['QV1_1',  1.0, QV_string], 
        ['QV1_7',  0.0, QV_string], 
        ['QV1_9',  0.0, QV_string], 
        ['QV1_17', 0.0, QV_string], 
        ['QV1_21', 0.0, QV_string], 
        ['QV1_31', 0.0, QV_string], 
        ['QV1_37', 0.0, QV_string], 
        ['QV2_0',  0.0, QV_string], 
        ['QV2_1',  1.0, QV_string], 
        ['QV2_7',  0.0, QV_string], 
        ['QV2_9',  0.0, QV_string], 
        ['QV2_17', 0.0, QV_string], 
        ['QV2_21', 0.0, QV_string], 
        ['QV2_31', 0.0, QV_string], 
        ['QV2_37', 0.0, QV_string], 
    ]
    prior_keyword = 'CD1_1'
    for new_item in new_WCS_keys:
    	keyword = new_item[0]
    	if keyword not in header:
    	    if debug: print "Inserting", keyword, new_item[1], new_item[2]
    	    header.update(keyword, new_item[1], new_item[2], after=prior_keyword)

        prior_keyword = keyword
# First QV keyword needs to come after last CD entry (CD2_2) not prior item in 
# list (CD2_1)
        if keyword == 'CD2_1':
          prior_keyword = 'CD2_2'
    return header

def update_sexaperd_keyword(header, debug=False):
    '''Updates the passed <header>'s SEXAPERD keyword (if present),
    replacing it with SEXAPED1-4 with the correct aperture sizes'''

# Fixed Aperture sizes in arcsec
    apers_arcsec = [1.0, 3.0, 5.0, 7.0]

# Check for required keywords
    if 'sexaperd' in header and 'pixscale' in header:
# Obtain pixelscale (in arcsec/pixel)
        platescale = header['pixscale']
        if debug: print "Platescale=", platescale
        if platescale != 'UNKNOWN':
# For each aperture, convert arcsec->pixels, round to 2 d.p and add to list
            num_ap = 1
            for aperture in apers_arcsec:
                aper_value = "%.2f" % (aperture / platescale )
                aper_value = float(aper_value)
                aper_keyword = 'SEXAPED' + str(num_ap)
                aper_comment = 'FIXED APERTURE DIAMETER #' + str(num_ap) + ' (PIXELS)'
                if debug: print aper_keyword,'=',aper_value,'/',aper_comment
# Need to treat first aperture differently as we need to update old SEXAPERD
# keyword to new values, then rename it
                if num_ap == 1:
                    header.update('sexaperd', aper_value, comment=aper_comment)
                    header.rename_key('sexaperd', aper_keyword)
                else:
                    header.update(aper_keyword, aper_value, comment=aper_comment, after=prior_keyword)
		num_ap = num_ap+1
                prior_keyword = aper_keyword

    return header

def verify_mjd_utc(hdrs,debug=False):
    '''Checks whether the UTSTART value from the header (as
    hdrs['EXPSTART_UTC']) matches the fractional part of MJD-OBS (to within 1ms).
    If it does not, then hdrs['MJD_UTC'] is replaced with the new value from 
    UTSTART, preserving the old value as hdrs['original_MJD_UTC']. 
    hdrs['DATE-OBS'] is also modified, pasting in the hdrs['EXPSTART_UTC']
    string after the '<yyyy>-<mm>-<dd>T'. The original hdrs['DATE-OBS'] is
    preserved as hdrs['original_DATE-OBS'].
    This works around 2 problems:
    1) There is a loss of precision somewhere in the Instrument Agent/pubsub 
    that means MJD-OBS differs from DATE-OBS/UTSTART by a few seconds.
    2) An AST bug that causes DATE-OBS and MJD-OBS in the output FITS to be 
    in ~TDB and not the requested UTC (Issue #4389)'''
 
    if 'DATEOBS_UTC' in hdrs and 'EXPSTART_UTC' in hdrs and 'MJD_UTC' in hdrs and 'DAYOBS_UTC' in hdrs:
    	if debug: print "Needed headers found"
# Replace colons with spaces
	iptr = 1
	utc_start_str = hdrs['EXPSTART_UTC'].replace(':',' ')
# Parse string and convert to radians
    	(iptr, utc_start , status)=S.sla_dafin(utc_start_str, iptr)
	if status == 0:
# If successful, convert to days
	    utc_start = degrees(utc_start)/24.0
	    mjd = modf(hdrs['MJD_UTC'])
# Compare with fractional part of MJD-OBS to a millisecond
	    if debug: print "UTC_START, MJD_OBS", utc_start, mjd[0]
	    if abs(utc_start-mjd[0]) >= (1e-3/86400.0):
	    	if debug: print "UTC_START != MJD_OBS", (utc_start-mjd[0])*86400.0
# Update MJD_UTC (MJD-OBS in header), preserve old value
    	    	new_mjd = mjd[1] + utc_start
                if (utc_start-mjd[0])*86400.0 > 85000.0:
                    print "TDB vs UTC wrap"
                    new_mjd = new_mjd - 1.0
		print "WARN: Modifying MJD_UTC to:", new_mjd
		hdrs['original_MJD_UTC'] = hdrs['MJD_UTC']
		hdrs['MJD_UTC'] = new_mjd
# Update DATEOBS_UTC (DATE-OBS in header), preserve old value
		loc_t = hdrs['DATEOBS_UTC'].find('T')
		if loc_t != -1:
    	    	    hdrs['original_DATEOBS_UTC'] = hdrs['DATEOBS_UTC']
# Check if the munged DATE-OBS(in TDB) has wrapped over midnight while the original UTC hadn't.
# If so we must use the original DAY-OBS keyword value to join back with the UTSTART
# XXX This code assumes there will not be more than 1400 leapseconds...
                    if (utc_start-mjd[0])*86400.0 > 85000.0:
                        print "TDB vs UTC wrap"
                        new_datestr =  "%4s-%2s-%2sT" % ( hdrs['DAYOBS_UTC'][0:4], 
                                                          hdrs['DAYOBS_UTC'][4:6], 
                                                          hdrs['DAYOBS_UTC'][6:8])
                        new_dateobs_utc = new_datestr+hdrs['EXPSTART_UTC']
                    else:
		        new_dateobs_utc = hdrs['DATEOBS_UTC'][0:loc_t+1]+hdrs['EXPSTART_UTC']
		    hdrs['DATEOBS_UTC'] = new_dateobs_utc
		    print "WARN: Modifying DATEOBS_UTC to:",new_dateobs_utc
	else:
	    if debug: print "Problems in conversion"
    return hdrs
    
def compute_tt_mid(hdrs,debug=False):
    '''Correct a UTC MJD start time to a TT MJD mid time. Also returns the
    TT-UTC offset for later use.'''

    UTC_mid_MJD = hdrs['MJD_UTC'] + (hdrs['EXPTIME'] / 86400.0 / 2.0)
    dtt = S.sla_dtt(UTC_mid_MJD)
    TT_mid_MJD = UTC_mid_MJD + dtt/86400.0
    if debug: print S.sla_dtt(UTC_mid_MJD)
    new_dict = { 'MJD_UTC_mid' : UTC_mid_MJD, 
    	    	 'MJD_TT_mid' : TT_mid_MJD,
		 'TT-UTC' : dtt }
    hdrs.update(new_dict)
    
    return hdrs

def compute_ut1(hdrs,debug=False):
    '''Computes UT1'''
    
    if hdrs.has_key('MJD_UTC_mid') and hdrs.has_key('UT1-UTC'):
    	ut1 = (hdrs['MJD_UTC_mid'] + (hdrs['UT1-UTC']/86400.0)) % 1.0
    	new_dict = { 'UT1' : ut1 }
	hdrs.update(new_dict)

    return hdrs

def sla_geoc_iers2003(p, h):
    '''Convert geodetic position to geocentric.
    *  Given:
    *     p     dp     latitude (geodetic, radians)
    *     h     dp     height above reference spheroid (geodetic, metres)
    *
    *  Returned:
    *     r     dp     distance from Earth axis (km)
    *     z     dp     distance from plane of Earth equator (km)
    *  Notes:
    *
    *  1  Geocentric latitude can be obtained by evaluating ATAN2(Z,R).
    *
    *  2  This version is an update of the original sla_geoc (which used 
    * 	  IAU 1976 constants) to use IERS2003 constants.
    *
    *  3  This version returns km for direct use with sla_RCC unlike the AU of the
    *     original'''

# Earth equatorial radius (metres)
    a0=6378137.0

# Reference spheroid flattening factor and useful function
    f=1.0/298.257223563
    b=(1.0-f)**2

# Geodetic to geocentric conversion
#
    sp = sin(p)
    cp = cos(p)
    c = 1.0 / sqrt(cp * cp + b * sp * sp )
    s = b * c
    r = ((a0*c + h) * cp)/1000.0
    z = ((a0*s + h) * sp)/1000.0
    
    return (r, z)
 
def compute_local_st(hdrs,debug=False):
    '''Compute the Local Apparent Sidereal Time. A TDB (ideally) or TT
    midtime must be provided as part of the <hdrs> dictionary otherwise
    an exception is raised.'''

# Compute the Greenwich mean sidereal time. This requires the full UT1 date
# as a MJD, not just the fraction.
    gmst = S.sla_gmst(hdrs['MJD_UTC_mid'] + (hdrs['UT1-UTC']/86400.0))

# Now local apparent sidereal time by adding longitude and the equation of the
# equinoxes (which needs the TDB date)
    if hdrs.has_key('MJD_TDB_mid'):
    	date = hdrs['MJD_TDB_mid']
	if debug: print "Using TDB time:", date
    elif hdrs.has_key('MJD_TT_mid'):
    	date = hdrs['MJD_TT_mid']
	if debug: print "Using TT time:", date
    else:
    	raise TimeError('Need a TDB or TT time available')
    
    last = gmst + hdrs['OBS_ELONG'] + S.sla_eqeqx(date)
    new_dict = { 'LST_App_mid' : last }
    hdrs.update(new_dict)

    return hdrs
    
def radec2strings(ra_radians, dec_radians, seperator=' '):
    '''Format an (RA, Dec) pair (in radians) into a tuple of strings, with 
    configurable seperator (defaults to <space>)'''

    ra_format = "%02.2d%c%02.2d%c%02.2d.%02.2d"
    dec_format = "%s%02.2d%c%02.2d%c%02.2d.%d"

    (rsign, ra ) = S.sla_dr2tf(2, ra_radians) 
    (dsign, dec) = S.sla_dr2af(1, dec_radians)

    ra_str = ra_format % ( ra[0], seperator, ra[1], seperator, ra[2],  ra[3] )
    dec_str = dec_format % ( dsign, dec[0], seperator, dec[1], seperator, dec[2], dec[3] )
    
    return (ra_str, dec_str)

