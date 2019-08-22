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
