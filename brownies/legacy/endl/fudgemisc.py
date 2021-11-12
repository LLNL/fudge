# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains useful top-level fudge routines that do not fit into any other module.
"""

import os
import sys
from . import fudgeParameters


def readOnly( ) :
    """Sets Fudge's internal readonly flag which causes fudge to not create working directories 
    for the isotopes read in. Note, you cannot unset this flag once readOnly has been called."""

    fudgeParameters.ReadOnly = 1

def verbose( mode = 1 ) :
    """Sets Fudge's internal verbose flag to mode. Mode = 0 suppress all informational messages."""

    fudgeParameters.VerboseMode = int(mode)

def checkMessagesToString( message, indentation = '  ', subIndentation = '  ' ) :

    s = []
    if( type( message ) == type( '' ) ) :
        s += [ indentation + message ]
    else :
        for m in message :
            if( type( m ) == type( '' ) ) :
                s += [ indentation + m  ]
            elif( ( type( m ) == type( [] ) ) or ( type( m ) == type( [] ) ) ) :
                if( len( m ) != 2 ) : raise Exception( 'Error in checkMessagesToString: len( m ) = %d != 2 )' % len( m ) )
                s += [ indentation + m[0] ]
                s += checkMessagesToString( m[1], indentation = indentation + subIndentation  )
            else :
                raise Exception( 'Error in checkMessagesToString: invalid message type = %s' % type( m ) )
    return( s )

def printWarning( s ) :

    sys.stderr.write( s )
    if( s[-1:] != '\n' ) : sys.stderr.write( '\n' )

def findPythonFile( file ) :
    """For internal use only."""

    if( os.path.exists( file ) ) : return file
    if( 'PYTHONPATH' in os.environ ) :
        ps = os.environ['PYTHONPATH'].split( ":" )
        for p in ps :
            f = os.path.join( p, file )
            if( os.path.exists( f ) ) : return f
    f = os.path.join( os.path.dirname( __file__ ), file )
    if( os.path.exists( f ) ) : return f
    f = os.path.join( "/usr/apps/fudge/current/bin", file )
    if( os.path.exists( f ) ) : return f
    return None

def stringWithPrefixSuffix( list, Prefix = "", Suffix = "" ) :
    """For internal use only."""

    PS = "\n"
    if( len( list ) > 0 ) :
        list[0] = Prefix + list[0]
        list[-1] = list[-1] + Suffix + "\n"
        PS = Suffix + "\n" + Prefix
    list = PS.join( list )
    return list

def getFormat( self ) :
    """Checks if self has a getFormat method. If it does not, then None is returned. If it does then: 
    1) if getFormat returns a string or None then it is returned else getFormat must return an 
    integer, n1. If n1 is 12 (default for legacy ENDL) then None is returned else returns string 
    '%n1.n2e' where n2 = n1 - 7."""

    try :
        format = self.getFormat( )
        if( ( format is None ) or ( type( format ) == type( '' ) ) ) : return( format )
        if( format == 12 ) :
            format = None
        else :
            format = '%%%d.%de' % ( format, format - 7 )
    except :
        format = None
    return( format )
