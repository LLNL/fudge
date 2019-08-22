# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""
This module contains useful fudge math routines that do not fit into any other module.
"""

import math

try :
    import numpy
    numpyFloat64 = numpy.float64( 1. )
except :
    numpyFloat64 = 1.

from pqu import PQU as PQUModule
import xData.standards as standardsModule

from fudge.core.utilities import brb

__metaclass__ = type

def thickenXYList( list, tester, biSectionMax = 6, interpolation = standardsModule.interpolation.linlinToken ) :
    """
    This functions takes a list of (x,y) points and a function, tester.evaluateAtX, and bi-sectionally adds points to 
    obtain linear-linear tolerance of the returned list and tester.evaluateAtX to tester.relativeTolerance. At most 
    biSectionMax bi-sections are performed between each consecutive pair of inputted points. It is assumed that the 
    initial list of points and the function tester.evaluateAtX agree to tolerance tester.relativeTolerance. The instance
    tester must contain the members relativeTolerance and absoluteTolerance as well as the method evaluateAtX. The
    method evaluateAtX takes an x-value and returns its y-value.
    """

    def thickenXYList2( interpolation, xl, yl, xu, yu, newList, tester, level ) :

        if( level == biSectionMax ) : return
        level += 1
        if( ( interpolation == standardsModule.interpolation.linlinToken ) or ( interpolation == standardsModule.interpolation.loglinToken ) ) :
            xMid = 0.5  * ( xl + xu )
        else :
            xMid = math.sqrt( xl * xu );

        if( ( interpolation == standardsModule.interpolation.linlinToken ) or ( interpolation == standardsModule.interpolation.linlogToken ) ) :
            yMid = 0.5  * ( yl + yu )
        else :
            yMid = math.sqrt( yl * yu )

        y = tester.evaluateAtX( xMid )

        dy = abs( y - yMid )
        if( ( dy > abs( y * tester.relativeTolerance ) ) and ( dy > tester.absoluteTolerance ) ) :
            newList.append( [ xMid, y ] )
            thickenXYList2( interpolation, xl, yl, xMid, y, newList, tester, level )
            thickenXYList2( interpolation, xMid, y, xu, yu, newList, tester, level )

    if( len( list ) < 2 ) : raise Exception( "len( list ) = %d < 2" % len( list ) )
    newList = []
    for i1, xy in enumerate( list ) :
        x2, y2 = xy
        if( i1 > 0 ) : thickenXYList2( interpolation, x1, y1, x2, y2, newList, tester, 0 )
        newList.append( [ x2, y2 ] )
        x1, y1 = x2, y2
    newList.sort( )
    return( newList )

def checkForNaN( v, str, printErrors = True, indentation = "", messages = None ) :

    if( v != v ) :
        s = "%s: value is nan" % str
        if( type( messages ) != type( None ) ) : messages.append( s )
        if( printErrors ) : printWarning( '\n'.join( checkMessagesToString( s, indentation = indentation ) ) ) 

def isNumber( n ) :

    if( type( 1. ) == type( n ) ) : return( True )
    if( type( 1 ) == type( n ) ) : return( True )
    if( type( numpyFloat64 ) == type( n ) ) : return( True )
    return( False )

def toInt( value ) :

    try :
        return( int( value ) )
    except :
        raise TypeError( 'value must be a convertible to a int; it is of %s: str( value )[:128] = "%s"' % ( brb.getType( value ), str( value )[:128] ) )

def toFloat( value ) :

    try :
        return( float( value ) )
    except :
        raise TypeError( 'value must be a convertible to a float; it is of %s: str( value )[:128] = "%s"' % ( brb.getType( value ), str( value )[:128] ) )

def checkNumber( v, str = "", printErrors = True, indentation = "", messages = None, maxAbsFloatValue = None ) :

    checkForNaN( v, str, printErrors = printErrors, indentation = indentation, messages = messages )
    if( not( maxAbsFloatValue is None ) ) :
        if( abs( v ) > maxAbsFloatValue ) :
            s = "%s: abs of number %s exceeds maxAbsFloatValue = %s" % ( str, v, maxAbsFloatValue )
            if( not( messages is None ) ) : messages.append( s )
            if( printErrors ) : printWarning( '\n'.join( checkMessagesToString( s, indentation = indentation ) ) )

def getValue( n ) :

    if( isNumber( n ) ) : return( n )
    if( isinstance( n, PQUModule.PQU ) ) : return( n.getValue( ) )
    raise Exception( 'Invalue number object = %s' % brb.getType( n ) )
