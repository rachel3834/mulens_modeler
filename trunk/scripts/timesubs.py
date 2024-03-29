############################################################
# TIMESUBS
#
# Routines taken wholesale from Tim Lister's new_NEO_Code package
############################################################

from datetime import datetime,timedelta
from pyslalib import slalib as S
import math

def datetime2mjd_utc(d):
    
# Compute MJD for UTC
    (mjd, status) = S.sla_cldj(d.year, d.month, d.day)
    if status != 0:
    	return None
    (fday, status ) = S.sla_dtf2d(d.hour, d.minute, d.second+(d.microsecond/1e6))
    if status != 0:
    	return None
    mjd_utc = mjd + fday

    return mjd_utc
    
def mjd_utc2mjd_tt(mjd_utc, dbg=False):
    '''Converts a MJD in UTC (MJD_UTC) to a MJD in TT (Terrestial Time) which is
    needed for any position/ephemeris-based calculations.
    UTC->TT consists of: UTC->TAI = 10s offset + 24 leapseconds (last one 2009 Jan 1.)
    	    	    	 TAI->TT  = 32.184s fixed offset'''
# UTC->TT offset
    tt_utc = S.sla_dtt(mjd_utc)
    if dbg: print 'TT-UTC(s)=', tt_utc

# Correct MJD to MJD(TT)
    mjd_tt = mjd_utc + (tt_utc/86400.0)
    if dbg: print 'MJD(TT)  =  ', mjd_tt

    return mjd_tt

def datetime2mjd_tdb(date, obsvr_long, obsvr_lat, obsvr_hgt, dbg=False):
    
    auinkm = 149597870.691
# Compute MJD_UTC from passed datetime
    mjd_utc = datetime2mjd_utc(date)
    if mjd_utc == None: return None

# Compute MJD_TT
    mjd_tt = mjd_utc2mjd_tt(mjd_utc, dbg)

# Compute TT->TDB

# Convert geodetic position to geocentric distance from spin axis (r) and from 
# equatorial plane (z)
    (r, z) = S.sla_geoc(obsvr_lat, obsvr_hgt)

    ut1 = compute_ut1(mjd_utc, dbg)
    if dbg: print "UT1=", ut1

# Compute relativistic clock correction TDB->TT    
    tdb_tt = S.sla_rcc(mjd_tt, ut1, -obsvr_long, r*auinkm, z*auinkm)
    if dbg: print "(TDB-TT)=", tdb_tt
    if dbg: print "(CT-UT)=", S.sla_dtt(mjd_utc)+tdb_tt

    mjd_tdb = mjd_tt + (tdb_tt/86400.0)
    
    return mjd_tdb
    
def ut1_minus_utc(mjd_utc, dbg=False):
    '''Compute UT1-UTC (in seconds), needed for tasks that require the Earth's orientation.
    UT1-UTC can be had from IERS Bulletin A (http://maia.usno.navy.mil/ser7/ser7.dat)
    but only for a short timespan and in arrears requiring continual downloading.
    Really need to get and read ftp://maia.usno.navy.mil/ser7/finals.all
    to get accurate UT1 value. Exercise for the reader...
    Currently we fake it by asuming 0.0. This will be wrong by at most +/- 0.9s 
    until they do away with leapseconds.'''

    dut = 0.0
    return dut

def compute_ut1(mjd_utc, dbg=False):
    '''Compute UT1 (as fraction of a day), needed for tasks that require the Earth's orientation.
    Currently we fake it by taking the fractional part of the day. This is good
    to +/- 0.9s until they do away with leapseconds.'''

    dut = ut1_minus_utc(mjd_utc)
    if dbg: print "DUT=", dut
    ut1 = (mjd_utc - int(mjd_utc)) + (dut/86400.0)

    return ut1

def parse_neocp_date(neocp_datestr, dbg=False):
    '''Parse dates from the NEOCP (e.g. '(Nov. 16.81 UT)' ) into a datetime object and 
    return this. No sanity checking of the input is done'''
    month_map = { 'Jan' : 1,
    	    	  'Feb' : 2,
		  'Mar' : 3,
		  'Apr' : 4,
		  'May' : 5,
		  'Jun' : 6,
		  'Jul' : 7,
		  'Aug' : 8,
		  'Sep' : 9,
		  'Oct' : 10,
		  'Nov' : 11,
		  'Dec' : 12 } 
		  
    chunks = neocp_datestr.split(' ')
    if dbg: print chunks
    if len(chunks) != 3: return None
    month_str = chunks[0].replace('(', '').replace('.', '')
    day_chunks = chunks[1].split('.') 
    if dbg: print day_chunks
    neocp_datetime = datetime(year=datetime.utcnow().year, month=month_map[month_str[0:3]],
    	day=int(day_chunks[0]))
	
    decimal_day = float('0.' + day_chunks[1].split()[0])
    neocp_datetime = neocp_datetime + timedelta(days=decimal_day)
    
    return neocp_datetime
    
def round_datetime(date_to_round, round_mins=10, round_up=False):
    '''Rounds the passed datetime object, <date_to_round>, to the 
    'floor' (default) or the 'ceiling' (if [roundup=True]) of
    the nearest passed amount (which defaults to 10min)'''
    
    correct_mins = 0
    if round_up:
        correct_mins = round_mins
    date_to_round = date_to_round - timedelta(minutes=(date_to_round.minute % round_mins)-correct_mins,
                	seconds=date_to_round.second,
                	microseconds=date_to_round.microsecond)

    return date_to_round
