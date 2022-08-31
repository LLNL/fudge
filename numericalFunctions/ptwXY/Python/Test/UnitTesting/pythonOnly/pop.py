# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import math
import random
import sys
sys.path.insert( 0, '../../Utilities' )

from numericalFunctions import pointwiseXY_C

import utilities

options = utilities.getOptions( __file__ )

xData = [ 0.5 * i - 0.145 for i in range( 51 ) ]
yData = [ math.sin( x ) for x in xData ]

xys = pointwiseXY_C.pointwiseXY_C( [ xData, yData ], initialSize = len( xData ), overflowSize = 10, dataForm = 'XsAndYs' )

n1 = len( xData )
for i1 in range( n1 ) :
    n2 = len( xData )
    i2 = random.randint( 0, n2 - 1 )
    x1 = xData.pop( i2 )
    y1 = yData.pop( i2 )
    x2, y2 = xys.pop( i2 )
    if( x1 != x2 ) : raise Exception( 'x1 = %.17e != x2 = %.17e, length left = %d' % ( x1, x2, len( xData ) ) )
    if( y1 != y2 ) : raise Exception( 'y1 = %.17e != y2 = %.17e, length left = %d' % ( y1, y2, len( xData ) ) )
    if( ( len( xData ) != len( yData ) ) or ( len( xData ) != len( xys ) ) ) : 
        raise Exception( 'len( xData ) = %d != len( yData ) = %d != len( xys ) = %d' % ( len( xData ), len( yData ), len( xys ) ) )

if( len( xData ) != 0 ) : raise Exception( 'len( xData ) = %d != 0' % len( xData ) )
if( len( yData ) != 0 ) : raise Exception( 'len( yData ) = %d != 0' % len( yData ) )
if( len( xys ) != 0 ) : raise Exception( 'len( xys ) = %d != 0' % len( xys ) )
