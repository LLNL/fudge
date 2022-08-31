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
domain = [ xData[0], xData[-1] ]
domain2 = xys.domain( )
domain3 = [ xys.domainMin( ), xys.domainMax( ) ]

if( domain != domain2 ) : raise Exception( 'xys.domain( ) failed "%s" vs. "%s"' % ( domain, domain2 ) )
if( domain != domain3 ) : raise Exception( 'xys.domainMin( ) and/or xys.domainMax( ) failed "%s" vs. "%s"' % ( domain, domain3 ) )
