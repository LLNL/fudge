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
This file compares the results of getValue from an endl2dmathClasses.endl2dmath and pointwiseXY_C
instances for the same data at random x values. If a relative difference greater than eps is found
an error is printed. This should be repeatable.
"""
import sys
sys.path.append( '../../../lib' )

import time
import random
import endl2dmathClasses, endlmisc
import pointwiseXY_C

eps = sys.float_info.epsilon
data = endl2dmathClasses.endl2dmath( endlmisc.read2dDataFile( 'Data/t3.py.in' ) )

f = pointwiseXY_C.pointwiseXY_C( initialSize = 40, overflowSize = 10 )
t0 = time.clock( )
f.setData( data.data )
print '# time =', time.clock( ) - t0
print 'len( data ) =', len( data )
print 'len( f ) =', len( f )

random.seed( 314159 )
n = 10000
np = n / 10
xMin, xMax = data.xMin( ), data.xMax( )
r = random.random()
x = xMin * r + xMax * ( 1. - r )
print '# time =', time.clock( ) - t0
for i in xrange( n ) :
    r = random.random( )
    x = xMin * r + xMax * ( 1. - r )
    y1 = data.getValue( x )
    y2 = f.getValue( x )
    if( abs( y2 - y1 ) > eps * ( y1 + y2 ) ) : print "#ERROR:", i, x, y1, y2
    if( ( i + 1 ) % np == 0 ) : print "%s of %s" % ( i + 1, n )
print '# time =', time.clock( ) - t0
