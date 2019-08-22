#! /bin/env python

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

import sys
sys.path.insert( 0, '../../../../../lib' )

import os
import pointwiseXY_C

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

def getXYData( ls, accuracy, interpolation = 'linear,linear' ) :

    ls, length = getIntegerValue( 'length', ls )
    data = [ map( float, ls[i].split( ) ) for i in xrange( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10, accuracy = accuracy, interpolation = interpolation )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    return( ls, data )

def checkIntepolation( xlog, ylog ) :

    xlogs = { 0 : 'linear', 1 : 'log' }
    ylogs = { 0 : 'linear', 1 : 'log' }
    f = open( os.path.join( CPATH, 'curve_sparse.dat' ) )
    ls = f.readlines( )
    f.close( )
    ls, xlog_ = getIntegerValue( 'xlog', ls )
    if( xlog != xlog_ ) : raise exception( 'xlog %s != xlog_ = %s' % ( xlog, xlog_ ) )
    ls, ylog_ = getIntegerValue( 'ylog', ls )
    if( ylog != ylog_ ) : raise exception( 'ylog %s != ylog_ = %s' % ( ylog, ylog_ ) )
    ls, accuracy = getDoubleValue( 'accuracy', ls )
    ls, sparseData = getXYData( ls, accuracy, interpolation = '%s,%s' % ( xlogs[xlog], ylogs[ylog] ) )

    f = open( os.path.join( CPATH, 'curve_dense.dat' ) )
    ls = f.readlines( )
    f.close( )
    ls, denseData = getXYData( ls, accuracy, interpolation = '%s,%s' % ( xlogs[xlog], ylogs[ylog] ) )

    f = open( os.path.join( CPATH, 'curve_linear.dat' ) )
    ls = f.readlines( )
    f.close( )
    ls, linearData = getXYData( ls, accuracy )

    n = sparseData.changeInterpolation( 'linear,linear', accuracy = accuracy )

    for x, y in denseData : 
        yn = n.getValue( x )
        if( abs( y - yn ) > accuracy * ( abs( y ) + abs( yn ) ) ) : 
            raise Exception( '%s: at %s values %e and %e diff by %e' % ( __file__, x, y, yn , y - yn ) )

xlogs = { 0 : '', 1 : '-xlog' }
ylogs = { 0 : '', 1 : '-ylog' }
for xlog in xlogs :
    for ylog in ylogs :
        if( '-e' in options ) : print __file__, xlogs[xlog], ylogs[ylog]
        os.system( 'cd %s; make -s clean; ./toLinearLinear %s %s' % ( CPATH, xlogs[xlog], ylogs[ylog] ) )
        checkIntepolation( xlog, ylog )
