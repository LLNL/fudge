# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

#
# This routine check that a raise in the applyFunction's callback function gets propagated up the stack.
#

import os
import sys

from numericalFunctions import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

CPATH = '../../../../Test/UnitTesting/toFunction'

accuracy = 1e-3

def func( x, parameters ) :

    self = parameters[0]
    index = self.lowerIndexBoundingX( x )
    self[len( self )]               # This should cause a python raise.
    return( x * x )

linear = pointwiseXY_C.pointwiseXY_C( [ [ -2, -2 ], [ 1, 1 ], [ 2, 2 ], [ 5, 5 ] ] )
parameters = [ linear ]
try :
    squared = linear.applyFunction( func, parameters, accuracy = accuracy, biSectionMax = 10 )
except :
    sys.exit( 0 )
raise Exception( 'The callback raise did not propagate properly.' )
