# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import copy

from numericalFunctions import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

def checkSlicing( XYs_base, i1, i2, doPrint = False ) :

    XYs = copy.copy( XYs_base )
    pXYs = pointwiseXY_C.pointwiseXY_C( data = XYs, initialSize = 20 )
    sXYs = XYs[i1:i2]
    if( doPrint ) :
        print( )
        print( i1, i2 )
        for xy in XYs : print( xy )
        print( )
        print( sXYs )
    spXYs = pXYs[i1:i2]
    if( doPrint ) : spXYs
    if( len( sXYs ) != len( spXYs ) ) : raise Exception( "%s: len( sXYs ) = %d != len( spXYs ) = %d: index1 = %d, i2 = %d" % \
        ( __file__, len( sXYs ), len( spXYs ), i1, i2 ) )

    for i, xy in enumerate( sXYs ) :
        if( xy[0] != spXYs[i][0] ) : raise Exception( "%s: difference at index = %d: %e %e" % ( __file__, i, xy[0], spXYs[i][0] ) )

    slice = pXYs.getslice( i1, i2 )
    for i, xy in enumerate( sXYs ) :
        if( xy[0] != slice[i][0] ) : raise Exception( "%s: getslice failed for difference at index = %d: %e %e" % ( __file__, i, xy[0], slice[i][0] ) )

XYs = [ [ float( x ), float( x )**2 ] for x in range( 12 ) ]

checkSlicing( XYs, 4, 8 )
checkSlicing( XYs, -4, 8 )
checkSlicing( XYs, -4, -2 )
checkSlicing( XYs, 4, -2 )
checkSlicing( XYs, 4, -8 )
checkSlicing( XYs, -4 - 7 * len( XYs ), 8 )
checkSlicing( XYs, 3, -4 - 7 * len( XYs ) )
