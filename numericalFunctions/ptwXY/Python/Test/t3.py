#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This file compares the results of evaluate from an endl2dmathClasses.endl2dmath and pointwiseXY_C
instances for the same data at random x values. If a relative difference greater than eps is found
an error is printed. This should be repeatable.
"""

import time
import random

from numericalFunctions import pointwiseXY_C

import endl2dmathClasses, endlmisc

eps = sys.float_info.epsilon
data = endl2dmathClasses.endl2dmath( endlmisc.read2dDataFile( 'Data/t3.py.in' ) )

f = pointwiseXY_C.pointwiseXY_C( initialSize = 40, overflowSize = 10 )
t0 = time.clock( )
f.setData( data.data )
print( '# time =', time.clock( ) - t0 )
print( 'len( data ) =', len( data ) )
print( 'len( f ) =', len( f ) )

random.seed( 314159 )
n = 10000
np = n / 10
xMin, xMax = data.domainMin( ), data.domainMax( )
r = random.random()
x = xMin * r + xMax * ( 1. - r )
print( '# time =', time.clock( ) - t0 )
for i in range( n ) :
    r = random.random( )
    x = xMin * r + xMax * ( 1. - r )
    y1 = data.evaluate( x )
    y2 = f.evaluate( x )
    if( abs( y2 - y1 ) > eps * ( y1 + y2 ) ) : print( "#ERROR:", i, x, y1, y2 )
    if( ( i + 1 ) % np == 0 ) : print( "%s of %s" % ( i + 1, n ) )
print( '# time =', time.clock( ) - t0 )
