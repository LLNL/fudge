#! /bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys
sys.path.insert( 0, '../../Utilities' )

from numericalFunctions import pointwiseXY_C

import utilities

options = utilities.getOptions( __file__ )

CPATH = '../../../../Test/UnitTesting/interpolation'
fIn = open( os.path.join( CPATH, 'strings.out' ) )
strings = []
for line in fIn.readlines( ) :
    if( line == '\n' ) : continue
    strings.append( line.split( '=' )[1].strip( ) )
fIn.close( )

def compareStrings( stringC, stringPy ) :

    _stringC = stringC[1:-1]
    if( _stringC != stringPy ) : raise Exception( '_stringC = "%s" not equal to stringPy = "%s"' % ( _stringC, stringPy ) )

def check( ptwXY2, interpolation, index ) :

    ptwXY = pointwiseXY_C.pointwiseXY_C( [ [ 1, 1 ], [ 10, 10 ] ], interpolation = interpolation )
    string = ptwXY.getInterpolation( )
    compareStrings( strings[index], string )
    if( '-v' in options ) : print( string )

    ptwXY2 = ptwXY.copy( )
    string = ptwXY.getInterpolation( )
    compareStrings( strings[index+1], string )
    if( '-v' in options ) : print( string )
    if( '-v' in options ) : print( string )

    return( ptwXY )

ptwXY2 = pointwiseXY_C.pointwiseXY_C( [ [ 1, 1 ], [ 10, 10 ] ], interpolation = 'charged-particle' )
ptwXY2 = check( ptwXY2, 'lin-lin', 0 )
ptwXY2 = check( ptwXY2, 'log-lin', 3 )
ptwXY2 = check( ptwXY2, 'lin-log', 6 )
ptwXY2 = check( ptwXY2, 'log-log', 9 )
ptwXY2 = check( ptwXY2, 'flat', 12 )
ptwXY2 = check( ptwXY2, 'charged-particle', 15 )
