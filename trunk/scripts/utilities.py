# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 23:52:35 2016

@author: robouser
"""

from numpy import pi

#####################
# SEX2DECDEG
def sex2decdeg(ra_str,dec_str):
    '''Function to convert an RA and Dec in sexigesimal format to decimal degrees'''
    
    ra_hrs = sexig2dec(ra_str)
    ra_deg = ra_hrs * 15.0
    dec_deg = sexig2dec(dec_str)
    
    return (ra_deg, dec_deg)

#####################
# SEX2RADS
def sex2rads(ra_str,dec_str):
    '''Function to convert an RA and Dec in sexigesimal format to decimal radians'''
    
    (ra_deg, dec_deg) = sex2decdeg(ra_str,dec_str)
    ra_rads = deg2rads(ra_deg)
    dec_rads = deg2rads(dec_deg)
    
    return (ra_rads, dec_rads)

#####################
# DEG2RAD
def deg2rads(angle_deg):
    """Function to convert decimal degrees to decimal radians"""
    
    angle_rad = ( pi * angle_deg ) / 180.0
    return angle_rad

#####################
# Function: SEXIG2DEC
def sexig2dec(CoordStr):
    '''Function to convert a sexigesimal coordinate string into a decimal float, returning a value in
        the same units as the string passed to it.'''
    
    # Ensure that the string is separated by ':':
    CoordStr = CoordStr.lstrip().rstrip().replace(' ',':')
    
    # Strip the sign, if any, from the first character, making a note if its negative:
    if CoordStr[0:1] == '-':
        Sign = -1.0
        CoordStr = CoordStr[1:]
    else:
        Sign = 1.0
    
    # Split the CoordStr on the ':':
    CoordList = CoordStr.split(':')
    
    # Assuming this presents us with a 3-item list, the last two items of which will
    # be handled as minutes and seconds:
    try:
        Decimal = float(CoordList[0]) + (float(CoordList[1])/60.0) + (float(CoordList[2])/3600.0)
        Decimal = Sign*Decimal
    except:
        Decimal = 0
    
    # Return with the decimal float:
    return Decimal
