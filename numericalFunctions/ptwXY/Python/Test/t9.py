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
This stuff should be repeatable.
"""
import sys
sys.path.append( '../../../lib' )

import pointwiseXY_C
import random

class t( pointwiseXY_C.pointwiseXY_C ) :

    def __init__( self, n, overflowSize = 10, xMin = -1, xMax = 1, yMin = 0, yMax = 100 ) :

        pointwiseXY_C.pointwiseXY_C.__init__( self, initialSize = n, overflowSize = 10 )
        for i in xrange( n ) :
            r = random.random( )
            x = xMin * r + xMax * ( 1. - r )
            r = random.random( )
            y = yMin * r + yMax * ( 1. - r )
            self.setValue( x, y )

random.seed( 314159 )
a = t( 5 )
b = t( 6 )

print a
print b
if( a[0][0] < b[0][0] ) :
    b.setValue( a[0][0], 0 )
else :
    a.setValue( b[0][0], 0 )

if( a[-1][0] > b[-1][0] ) :
    b.setValue( a[-1][0], 0 )
else :
    a.setValue( b[-1][0], 0 )

print a
print b
print a + b
