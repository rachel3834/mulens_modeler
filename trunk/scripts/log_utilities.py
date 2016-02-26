# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 00:00:05 2016

@author: robouser
"""

#############################################################################
#                       LOGGING UTILITIES
#############################################################################

import logging
from os import path
from astropy.time import Time
import glob

def start_day_log( config, log_name ):
    """Function to initialize a new log file.  The naming convention for the
    file is [log_name]_[UTC_date].log.  A new file is automatically created 
    if none for the current UTC day already exist, otherwise output is appended
    to an existing file.
    This function also configures the log file to provide timestamps for 
    all entries.  
    """

    ts = Time.now()    
    ts = ts.iso.split()[0]
    
    log_file = path.join( config['log_directory'], \
                    config['log_root_name'] + '_' + ts + '.log' )

    # Look for previous logs and rollover the date if the latest log
    # isn't from the curent date:


    # To capture the logging stream from the whole script, create
    # a log instance together with a console handler.  
    # Set formatting as appropriate.
    log = logging.getLogger( log_name )
    
    if len(log.handlers) == 0:
        log.setLevel( logging.INFO )
        file_handler = logging.FileHandler( log_file )
        file_handler.setLevel( logging.INFO )
    
        console_handler = logging.StreamHandler()
        console_handler.setLevel( logging.INFO )
    
        formatter = logging.Formatter( fmt='%(asctime)s %(message)s', \
                                    datefmt='%Y-%m-%dT%H:%M:%S' )
        file_handler.setFormatter( formatter )
        console_handler.setFormatter( formatter )
    
        log.addHandler( file_handler )
        log.addHandler( console_handler )
    
    log.info( '\n------------------------------------------------------\n')
    return log
    
def end_day_log( log ):
    """Function to cleanly shutdown logging functions with last timestamped
    entry"""
    
    log.info( 'Processing complete\n' )
    logging.shutdown()

