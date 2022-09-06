# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from numericalFunctions import pointwiseXY_C

accuracy = 1e-2
biSectionMax = 0.

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

CPATH = '../../../../Test/Flat/binaryMath'

os.system( 'cd %s; make -s clean; ./flatMath -v > v' % CPATH )

def skipBlankLines( ls ) :

    i = 0
    for i, l in enumerate( ls ) :
        if( l.strip( ) != '' ) : break
    ls = ls[i:]
    if( ( len( ls ) == 1 ) and ( ls[0].strip( ) == '' ) ) : ls = []
    return( ls )

def getIntegerValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = int( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def getDoubleValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = float( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.12e' % v1, '%.12e' % v2
    sv1, sv2 = '%.7e' % float( sv1 ), '%.7e' % float( sv2 )
    if( sv1 != sv2 ) : print( '<%s> <%s>' % ( sv1, sv2 ) )
    if( sv1 != sv2 ) : raise Exception( '%s: values %e and %e diff by %e at %d for label = %s' % ( __file__, v1, v2, v2 - v1, i, label ) )

def getXYData( ls, biSectionMax, accuracy ) :

    ls, length = getIntegerValue( 'length', ls )
    data = [ list( map( float, ls[i].split( ) ) ) for i in range( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10, biSectionMax = biSectionMax, accuracy = accuracy, safeDivide = True, interpolation = "flat" )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    return( ls, data )

def getCommand( ls ) :

    s = ls[0].split( )
    if( len( s ) != 2 ) : raise Exception( 'Invalid command = "%s"' % ls[0][:-1] )
    if( s[0] != "#" ) : raise Exception( 'Invalid command = "%s"' % ls[0][:-1] )
    return( ls[1:], s[1] )

def compareXYs( XYs1, XYs2, label ) :

    if( len( XYs1 ) != len( XYs2 ) ) : raise Exception( 'for %s: len( XYs1 ) = %s != len( XYs2 ) = %s' % ( label, len( XYs1 ), len( XYs2 ) ) )
    for i, xy in enumerate( XYs1 ) :
        compareValues( "x division " + label, count, xy[0], XYs2[i][0] )
        compareValues( "y division " + label, count, xy[1], XYs2[i][1] )

def mathParse( count, ls ) :

    ls, command = getCommand( ls ) 
    if( command == 'double' ) :
        ls = doubleCheck( count, ls )
    elif( command == 'all_double' ) :
        ls = allDoubleCheck( count, ls )
    elif( command == 'binary_add_sub' ) :
        ls = binaryAddSubCheck( count, ls )
    elif( command == 'binary_mul_div' ) :
        ls = binaryMulDivCheck( count, ls )
    else :
        raise  Exception( 'Invalid command = "%s"' % command )
    return( ls )

def doubleCheck( count, ls ) :

    ls, d = getDoubleValue( 'double', ls )
    ls, o = getCommand( ls )
    if( o not in '+-=*/\\' ) : raise Exception( 'Unknown operator "%s"' % o )
    ls, XYs = getXYData( ls, biSectionMax, accuracy )
    ls, resultsC = getXYData( ls, biSectionMax, accuracy )
    if(   o == '+'  ) : results = XYs + d
    elif( o == '-'  ) : results = XYs - d
    elif( o == '='  ) : results = d - XYs
    elif( o == '*'  ) : results = XYs * d
    elif( o == '/'  ) : results = XYs / d
    elif( o == '\\' ) : results = d / XYs
    compareXYs( resultsC, results, "doubleCheck %s" % o )
    return( ls )

def allDoubleCheck( count, ls ) :

    ls, d = getDoubleValue( 'double', ls )
    ls, XYs = getXYData( ls, biSectionMax, accuracy )
    ls, resultsC = getXYData( ls, biSectionMax, accuracy )
    results = ( ( d * ( XYs + d ) ) - d ) / d
    results = ( ( ( d * results ) + d ) / d ) - d
    compareXYs( resultsC, results, "allDoubleCheck" )
    return( ls )

def binaryAddSubCheck( count, ls ) :

    ls, XYs1 = getXYData( ls, biSectionMax, accuracy )
    ls, XYs2 = getXYData( ls, biSectionMax, accuracy )
    ls, resultsC = getXYData( ls, biSectionMax, accuracy )
    results = XYs1 + XYs2
    compareXYs( resultsC, results, "binaryAddSubCheck" )
    ls, dummy = getXYData( ls, biSectionMax, accuracy )
    ls, resultsC = getXYData( ls, biSectionMax, accuracy )
    results = results - XYs2
    compareXYs( resultsC, results, "binaryAddSubCheck" )
    return( ls )

def binaryMulDivCheck( count, ls ) :

    ls, XYs1 = getXYData( ls, biSectionMax, accuracy )
    ls, XYs2 = getXYData( ls, biSectionMax, accuracy )
    ls, resultsC = getXYData( ls, biSectionMax, accuracy )
    results = XYs1 * XYs2
    compareXYs( resultsC, results, "binaryMulDivCheck" )
    ls, dummy = getXYData( ls, biSectionMax, accuracy )
    ls, resultsC = getXYData( ls, biSectionMax, accuracy )
    results = results / XYs2
    compareXYs( resultsC, results, "binaryMulDivCheck" )
    return( ls )

f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

ls, accuracy = getDoubleValue( 'accuracy', ls )

count = 0
while( len( ls ) ) :
    count += 1
    ls = mathParse( count, ls )
