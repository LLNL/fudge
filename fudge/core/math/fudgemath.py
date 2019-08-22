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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
This module contains useful fudge math routines that do not fit into any other module.
"""

from pqu import physicalQuantityWithUncertainty
from fudge.core.utilities import brb
try :
    import numpy
    numpyFloat64 = numpy.float64( 1. )
except :
    numpyFloat64 = 1.

__metaclass__ = type

def thickenXYList( list, tester, biSectionMax = 6 ) :
    """This functions takes a list of (x,y) points and a function, tester.evaluateAtX, and bi-sectionally adds points to 
    obtain linear-linear tolerance of the returned list and tester.evaluateAtX to tester.relativeTolerance. At most 
    biSectionMax bi-sections are performed between each consecutive pair of inputted points. It is assumed that the 
    initial list of points and the function tester.evaluateAtX agree to tolerance tester.relativeTolerance. The instance
    tester must contain the members relativeTolerance and absoluteTolerance as well as the method evaluateAtX. The
    method evaluateAtX takes an x-value and returns its y-value."""

    def thickenXYList2( xl, yl, xu, yu, newList, tester, level ) :

        if( level == biSectionMax ) : return
        level += 1
        xMid = 0.5  * ( xl + xu )
        yMid = 0.5  * ( yl + yu )
        y = tester.evaluateAtX( xMid )
        dy = abs( y - yMid )
        if( ( dy > abs( y * tester.relativeTolerance ) ) and ( dy > tester.absoluteTolerance ) ) :
            newList.append( [ xMid, y ] )
            thickenXYList2( xl, yl, xMid, y, newList, tester, level )
            thickenXYList2( xMid, y, xu, yu, newList, tester, level )

    if( len( list ) < 2 ) : raise Exception( "len( list ) = %2 < 2" % len( list ) )
    newList = []
    for i, xy in enumerate( list ) :
        x2, y2 = xy
        if( i > 0 ) : thickenXYList2( x1, y1, x2, y2, newList, tester, 0 )
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
    if( isinstance( n, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) : return( n.getValue( ) )
    raise Exception( 'Invalue number object = %s' % brb.getType( n ) )
