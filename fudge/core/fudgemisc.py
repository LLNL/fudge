# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
This module contains useful top-level fudge routines that do not fit into any other module.
"""

import os
import sys
from fudge import fudgeParameters

def readOnly( ) :
    """Sets Fudge's internal readonly flag which causes fudge to not create working directories 
    for the isotopes read in. Note, you cannot unset this flag once readOnly has been called."""

    fudgeParameters.ReadOnly = 1

def verbose( mode = 1 ) :
    """Sets Fudge's internal verbose flag to mode. Mode = 0 suppress all informational messages."""

    fudgeParameters.VerboseMode = int( mode )

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
    "For internal use only."

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
    "For internal use only."

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
