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

value = 2.4
orig = pointwiseXY_C.pointwiseXY_C( [ [ 0, 1 ], [ 1, 3 ], [ 3, -1 ] ] )
errs = []

def check( operator, other ) :

    global orig
    m1 = eval( 'orig %s other' % operator )
    m2 = orig.copy( )
    m3 = m2
    exec( 'm2 %s= other' % operator )
    d1 = m2 - m1
    if( d1.rangeMin( ) != 0. ) : errs.append( 'failure for operator "%s": d1.rangeMin( ) = %e != 0.' % ( operator, d1.rangeMin( ) ) )
    if( d1.rangeMax( ) != 0. ) : errs.append( 'failure for operator "%s": d1.rangeMax( ) = %e != 0.' % ( operator, d1.rangeMax( ) ) )
    if( m2 is not m3 ) : errs.append( 'failure for operator "%s": m2 is not m3' % operator )

check( '+', value )
check( '-', value )
check( '*', value )
check( '/', value )

other = pointwiseXY_C.pointwiseXY_C( [ [ 0, 1 ], [ 1, 4 ], [ 3, 3 ] ] )
check( '+', other )
check( '-', other )
check( '*', other )
check( '/', other )

if( len( errs ) > 0 ) :
    for err in errs : print( '   ', err )
    raise Exception( '%s: in place binary operator failed' % __file__ )
