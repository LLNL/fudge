#! /bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

import sys, math
sys.path.insert( 0, '../../Utilities' )
sys.path.insert( 0, '../../../../../lib' )

import pointwiseXY_C
import utilities
options = utilities.getOptions( __file__ )

xData = [ 0.5 * i - 0.145 for i in xrange( 51 ) ]
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
    if( d.yMin( ) != 0 ) : status = 1
    if( d.yMax( ) != 0 ) : status = 1
    return( status )

status = cmp( dataXYs, dataXYs2 )
status += cmp( dataXYs, dataXsAndYs )
status += cmp( dataXYs, dataList )
if( status ) : raise Exception( '%s: %d sets not the same' % ( __file__, status ) )
