# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys

from numericalFunctions import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

status = 0

def check( xys, x, correctIndex ) :

    global status

    index = xys.lowerIndexBoundingX( x )
    if( correctIndex != index ) :
        if( status == 0 ) : print( 'Error from %s: correctIndex = %d for x = %e and got index = %d' % ( __file__, correctIndex, x, index ) )
        status += 1

data = [ [ -1.1, 0 ], [ 0.7, 1 ], [ 10, 0 ], [ 1001.1, 1 ] ]
xys = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10 )

i1 = -1
x1 = data[0][0] - 1
for i2, ( x2, y ) in enumerate( data ) :
    x = 0.5 * ( x1 + x2 )
    check( xys, x, i1 )
    i1 = i2
    x1 = x2

i1 = check( xys, x1 + 1, -1 )
sys.exit( status )
