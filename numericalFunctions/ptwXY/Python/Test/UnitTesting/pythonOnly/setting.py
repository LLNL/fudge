# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import math
import sys
sys.path.insert( 0, '../../Utilities' )

from numericalFunctions import pointwiseXY_C

import utilities

options = utilities.getOptions( __file__ )

xData = [ 0.5 * i - 0.145 for i in range( 51 ) ]
yData = [ math.sin( x ) for x in xData ]
XYData = [ [ x, yData[i] ] for i, x in enumerate( xData ) ]
listData = []
for i, x in enumerate( xData ) :
    listData.append( x )
    listData.append( yData[i] )

dataXYs = pointwiseXY_C.pointwiseXY_C( XYData, initialSize = len( XYData ), overflowSize = 10 )
dataXYs2 = pointwiseXY_C.pointwiseXY_C( XYData, initialSize = len( XYData ), overflowSize = 10, dataForm = 'xys' )
dataXsAndYs = pointwiseXY_C.pointwiseXY_C( [ xData, yData ], initialSize = len( XYData ), overflowSize = 10, dataForm = 'XsAndYs' )
dataList = pointwiseXY_C.pointwiseXY_C( listData, initialSize = len( XYData ), overflowSize = 10, dataForm = 'List' )

def cmp( p1, p2 ) :

    status = 0
    d = p1 - p2
    if( d.rangeMin( ) != 0 ) : status = 1
    if( d.rangeMax( ) != 0 ) : status = 1
    return( status )

status = cmp( dataXYs, dataXYs2 )
status += cmp( dataXYs, dataXsAndYs )
status += cmp( dataXYs, dataList )
if( status ) : raise Exception( '%s: %d sets not the same' % ( __file__, status ) )
