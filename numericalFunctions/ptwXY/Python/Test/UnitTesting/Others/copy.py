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

def compareValues( label, i, v1, v2 ) :

    if( v1 != v2 ) : raise Exception( '%s: values %e and %e diff by %e at %d for label = %s' % ( __file__, v1, v2, v2 - v1, i, label ) )

def compareXYs( xys1, xys2 ) :

    if( len( xys1 ) != len( xys2 ) ) : raise Exception( "len( xys1 ) = %s != len( xys2 ) = %d" % ( len( xys1 ), len( xys2 ) ) )
    for i, xy1 in enumerate( xys1 ) :
        xy2 = xys2[i]
        compareValues( 'x', i, xy1[0], xy2[0] )
        compareValues( 'x', i, xy1[1], xy2[1] )

xs = [ x + .1 for x in range( 2, 33 ) ]
ys = [ x * x for x in xs ]
xys = list( zip( xs, ys ) )
xylist = []
for i, x in enumerate( xs ) :
    xylist.append( x )
    xylist.append( ys[i] )

xys1 = pointwiseXY_C.pointwiseXY_C( data = xys )
xys2 = pointwiseXY_C.pointwiseXY_C( data = xys, dataForm = 'XYs' )
compareXYs( xys1, xys2 )

xys2 = pointwiseXY_C.pointwiseXY_C( data = [ xs, ys ], dataForm = 'xsandys' )
compareXYs( xys1, xys2 )

xys2 = pointwiseXY_C.pointwiseXY_C( data = [ xs, ys ], dataForm = 'XsAndYs' )
compareXYs( xys1, xys2 )

xys2 = pointwiseXY_C.pointwiseXY_C( data = xylist, dataForm = 'List' )
compareXYs( xys1, xys2 )

xys2 = xys1.copy( )
compareXYs( xys1, xys2 )

xys2 = xys1.copyDataToXYs( )
compareXYs( xys1, xys2 )

xs2, ys2 = xys1.copyDataToXsAndYs( )
if( xs != xs2 ) : raise Exception( "xs != xs2" )
if( ys != ys2 ) : raise Exception( "ys != ys2" )

xScale = 3.14
yScale = 2.3e-6

for i in range( len( xs ) ) :
    xs[i] *= xScale
    ys[i] *= yScale
xys3 = pointwiseXY_C.pointwiseXY_C( data = [ xs, ys ], dataForm = 'XsAndYs' )

xys2 = xys1.copyDataToXYs( xScale = xScale, yScale = yScale )
compareXYs( xys3, xys2 )

xs2, ys2 = xys1.copyDataToXsAndYs( xScale = xScale, yScale = yScale )
if( xs != xs2 ) : raise Exception( "xs != xs2" )
if( ys != ys2 ) : raise Exception( "ys != ys2" )
