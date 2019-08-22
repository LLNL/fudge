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

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.7g' % v1, '%.7g' % v2
    if( sv1 != sv2 ) : raise Exception( 'Values %s %s diff at %d for label = %s' % ( v1, v2, i, label ) )

def compareData( data1, data2 ) :

    if( len( data1 ) != len( data2 ) ) : raise Exception( 'len( data1 ) = %s != len( data2 ) = %s' % ( len( data1 ), len( data2 ) ) )
    for i, xy1 in enumerate( data1 ) :
        xy2 = data1[i]
        compareValues( 'x-value', i, xy1[0], xy2[0] )
        compareValues( 'y-value', i, xy1[1], xy2[1] )

def compareAxis( label, data1, data2 ) :

    if( len( data1 ) != len( data2 ) ) : raise Exception( 'len( data1 ) = %s != len( data2 ) = %s' % ( len( data1 ), len( data2 ) ) )
    for i, xy1 in enumerate( data1 ) :
        compareValues( label, i, xy1, data2[i] )


import random
random.seed( 314159 )
xMin, xMax = -100, 50
yMin, yMax = -1., 10.
n = 100

xs, ys, xys = [], [], []
for i in xrange( n ) :
    x, y = ( xMax - xMin ) * random.random() + xMin, ( yMax - yMin ) * random.random() + yMin
    xys.append( [ x, y ] )
xys.sort( )                 # Hopefully, no two x-values are the same.
for x, y in xys :
    xs.append( x )
    ys.append( y )

xysPy = pointwiseXY_C.pointwiseXY_C( data = xys )

xysPy2 = pointwiseXY_C.pointwiseXY_C( data = [ xs, ys ], dataForm = 'xsAndYs' )
compareData( xysPy, xysPy2 )

xysPy2 = pointwiseXY_C.pointwiseXY_C( )
xysPy2.setData( xys )
compareData( xysPy, xysPy2 )

xysPy2 = pointwiseXY_C.pointwiseXY_C( )
xysPy2.setDataFromXsAndYs( xs, ys )
compareData( xysPy, xysPy2 )

copiedXYs = xysPy.copyDataToXYs( )
compareData( xys, copiedXYs )

copiedXs, copiedYs = xysPy.copyDataToXsAndYs( )
compareAxis( 'x-value', xs, copiedXs )
compareAxis( 'y-value', ys, copiedYs )
