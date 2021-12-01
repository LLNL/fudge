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

CPATH = '../../../../Test/UnitTesting/Others'

os.system( 'cd %s; make -s clean; ./trim -v > v' % CPATH )

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

def getXYData( ls ) :

    ls, length = getIntegerValue( 'length', ls )
    data = [ list( map( float, ls[i].split( ) ) ) for i in range( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10 )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    return( ls, data )

def checkClipping( count, ls ) :

    ls, data = getXYData( ls )
    ls, trimmedC = getXYData( ls )
    trimmed = data.trim( )
    if( len( trimmed ) != len( trimmedC ) ) : raise Exception( '%s: at %d len( trimmed ) = %d != len( trimmedC ) = %d' % \
        ( __file__, count, len( trimmed ), len( trimmedC ) ) )
    for i, xy in enumerate( trimmed ) :
        compareValues( "x trimmed", count, xy[0], trimmedC[i][0] )
        compareValues( "y trimmed", count, xy[1], trimmedC[i][1] )
    return( ls )

f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

count = 0
while( len( ls ) ) :
    count += 1
    if( count == 6 ) :
        ls, dummy = getXYData( ls )
        ls, dummy = getXYData( ls )
    ls = checkClipping( count, ls )
