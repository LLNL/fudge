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
rangeMin = min( [ y for x, y in xys ] )
rangeMax = max( [ y for x, y in xys ] )
range1 = [ rangeMin, rangeMax ]
range2 = xys.range( )
range3 = [ xys.rangeMin( ), xys.rangeMax( ) ]

if( range1 != range2 ) : raise Exception( 'xys.range( ) failed "%s" vs. "%s"' % ( range1, range2 ) )
if( range1 != range3 ) : raise Exception( 'xys.rangeMin( ) and/or xys.rangeMax( ) failed "%s" vs. "%s"' % ( range1, range3 ) )
