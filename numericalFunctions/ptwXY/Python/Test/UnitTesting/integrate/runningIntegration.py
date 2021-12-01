# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from numericalFunctions import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

CPATH = '../../../../Test/UnitTesting/integrate'

os.system( 'cd %s; make -s clean; ./runningIntegration -v > v' % CPATH )

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

def compareValues( i, v1, v2 ) :

    d = v1 - v2
    r = d
    if( v2 != 0 ) : r /= v2
    if( abs( r ) > 1e-10 ) : raise Exception( '%s: values %e and %e diff by %e or relatively %e at %d' % ( __file__, v1, v2, d, r, i ) )

def getXYData( ls ) :

    ls = skipBlankLines( ls )
    ls, length = getIntegerValue( 'length', ls )
    data = [ list( map( float, ls[i].split( ) ) ) for i in range( length ) ]
    XYs, sums = [], []
    for x, y, s in data :
        XYs.append( ( x, y ) )
        sums.append( s )
    XYs = pointwiseXY_C.pointwiseXY_C( XYs, initialSize = len( data ), overflowSize = 10 )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    return( ls, XYs, sums )

def checkIntegration( ls ) :

    ls, XYs, sums = getXYData( ls )
    runningSum = XYs.runningIntegral( )
    for i, s in enumerate( runningSum ) :
        compareValues( i, s, sums[i] )
    return( ls )

f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

while( len( ls ) ) : ls = checkIntegration( ls )
