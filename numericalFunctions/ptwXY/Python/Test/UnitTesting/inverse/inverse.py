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

CPATH = '../../../../Test/UnitTesting/inverse'

os.system( 'cd %s; make -s clean; ./inverse -v > v' % CPATH )

def skipBlankLines( ls ) :

    i = 0
    for i, l in enumerate( ls ) :
        if( l.strip( ) != '' ) : break
    ls = ls[i:]
    if( ( len( ls ) == 1 ) and ( ls[0].strip( ) == '' ) ) : ls = []
    return( ls )

def getStringValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = ls[0].split( '=' )[1].strip( )
    return( ls[1:], value )

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

def compare( count, xy1, xy2 ) :

    if( len( xy1 ) != len( xy2 ) ) : raise Exception( "count = %d: len( xy1 ) = %d != len( xy2 ) = %d" % ( count, len( xy1 ), len( xy2 ) ) )

def getXYData( ls ) :

    ls, length = getIntegerValue( 'length', ls )
    ls, interpolation = getIntegerValue( 'interpolation', ls )
    ls, interpolationStr = getStringValue( 'interpolation string', ls )
    data = [ list( map( float, ls[i].split( ) ) ) for i in range( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), interpolation = interpolationStr )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    return( ls, data )

def inverse( count, n1, i1 ) :

    if( len( n1 ) != len( i1 ) ) : raise Exception( "count = %d: len( n1 ) = %d != len( i1 ) = %d" % ( count, len( n1 ), len( i1 ) ) )
    i2 = n1.inverse( )
    compare( count, i1, i2 )

def checkInverse( count, ls ) :

    ls, norm1 = getXYData( ls )
    ls, inverse1 = getXYData( ls )
    ls, norm2 = getXYData( ls )

    inverse( count, norm1, inverse1 )
    inverse( count, inverse1, norm2 )
    compare( count, norm1, norm2 )
    return( ls )

f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

count = 0
while( len( ls ) ) :
    count += 1
    ls = checkInverse( count, ls )
