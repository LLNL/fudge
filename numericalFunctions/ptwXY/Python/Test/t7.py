#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This is a test of the method copy. Creates a pointwiseXY_C object and makes a copy. The outputs should be the same.
"""

import time
import random

from numericalFunctions import pointwiseXY_C

d = pointwiseXY_C.pointwiseXY_C( initialSize = 100000, overflowSize = 3100 )

xMin, xMax = 1, 100
yMin, yMax = -1., 0.
n = 100000
t0 = time.process_time( )
random.seed( 314159 )
for i in range( n ) :
    r = random.random()
    x = xMin * r + xMax * ( 1. - r )
    r = random.random()
    y = yMin * r + yMax * ( 1. - r )
    d.setValue( x, y )

print( '# time', time.process_time( ) - t0 )
d2 = d.copy( )
print( '# time', time.process_time( ) - t0 )
print( len( d ), len( d2 ) )

print( '# time', time.process_time( ) - t0 )
f = open( 'Temp/t7.orig.out', 'w' )
for x, y in d : f.write( '%.12f %.12f\n' % ( x, y ) )
f.close( )

f = open( 'Temp/t7.copy.out', 'w' )
for x, y in d2 : f.write( '%.12f %.12f\n' % ( x, y ) )
f.close( )
