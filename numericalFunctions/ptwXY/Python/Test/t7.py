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

"""
This is a test of the method copy. Creates a pointwiseXY_C object and makes a copy. The outputs should be the same.
"""
import sys
sys.path.append( '../../../lib' )

import time
import random
import pointwiseXY_C

d = pointwiseXY_C.pointwiseXY_C( initialSize = 100000, overflowSize = 3100 )

xMin, xMax = 1, 100
yMin, yMax = -1., 0.
n = 100000
t0 = time.clock( )
random.seed( 314159 )
for i in xrange( n ) :
    r = random.random()
    x = xMin * r + xMax * ( 1. - r )
    r = random.random()
    y = yMin * r + yMax * ( 1. - r )
    d.setValue( x, y )

print '# time', time.clock( ) - t0
d2 = d.copy( )
print '# time', time.clock( ) - t0
print len( d ), len( d2 )

print '# time', time.clock( ) - t0
f = open( 'Temp/t7.orig.out', 'w' )
for x, y in d : f.write( '%.12f %.12f\n' % ( x, y ) )
f.close( )

f = open( 'Temp/t7.copy.out', 'w' )
for x, y in d2 : f.write( '%.12f %.12f\n' % ( x, y ) )
f.close( )
