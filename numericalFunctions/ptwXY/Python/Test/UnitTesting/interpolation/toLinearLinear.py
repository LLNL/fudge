# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from numericalFunctions import pointwiseXY_C

options = ''
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )

CPATH = '../../../../Test/UnitTesting/interpolation'

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

def getXYData( ls, accuracy, interpolation = 'lin-lin' ) :

    ls, length = getIntegerValue( 'length', ls )
    data = [ list( map( float, ls[i].split( ) ) ) for i in range( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10, accuracy = accuracy, interpolation = interpolation )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    return( ls, data )

def checkIntepolation( xlog, ylog ) :

    xlogs = { 0 : 'lin', 1 : 'log' }
    ylogs = { 0 : 'lin', 1 : 'log' }
    f = open( os.path.join( CPATH, 'curve_sparse.dat' ) )
    ls = f.readlines( )
    f.close( )
    ls, xlog_ = getIntegerValue( 'xlog', ls )
    if( xlog != xlog_ ) : raise exception( 'xlog %s != xlog_ = %s' % ( xlog, xlog_ ) )
    ls, ylog_ = getIntegerValue( 'ylog', ls )
    if( ylog != ylog_ ) : raise exception( 'ylog %s != ylog_ = %s' % ( ylog, ylog_ ) )
    ls, accuracy = getDoubleValue( 'accuracy', ls )
    ls, sparseData = getXYData( ls, accuracy, interpolation = '%s-%s' % ( ylogs[ylog], xlogs[xlog] ) )

    f = open( os.path.join( CPATH, 'curve_dense.dat' ) )
    ls = f.readlines( )
    f.close( )
    ls, denseData = getXYData( ls, accuracy, interpolation = '%s-%s' % ( ylogs[ylog], xlogs[xlog] ) )

    f = open( os.path.join( CPATH, 'curve_linear.dat' ) )
    ls = f.readlines( )
    f.close( )
    ls, linearData = getXYData( ls, accuracy )

    n = sparseData.changeInterpolation( 'lin-lin', accuracy = accuracy )

    for x, y in denseData : 
        yn = n.evaluate( x )
        if( abs( y - yn ) > accuracy * ( abs( y ) + abs( yn ) ) ) : 
            raise Exception( '%s: at %s values %e and %e diff by %e' % ( __file__, x, y, yn , y - yn ) )

xlogs = { 0 : '', 1 : '-xlog' }
ylogs = { 0 : '', 1 : '-ylog' }
for xlog in xlogs :
    for ylog in ylogs :
        if( '-e' in options ) : print( __file__, xlogs[xlog], ylogs[ylog] )
        os.system( 'cd %s; make -s clean; ./toLinearLinear %s %s' % ( CPATH, xlogs[xlog], ylogs[ylog] ) )
        checkIntepolation( xlog, ylog )
