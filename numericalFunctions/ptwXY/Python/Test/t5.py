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
Some simple stuff that is not repeatable as pointers are printed.
"""
import sys
sys.path.append( '../../../lib' )

import random
import pointwiseXY_C

f = pointwiseXY_C.pointwiseXY_C( initialSize = 7, overflowSize = 3 )

xMin, xMax = -100, 10
yMin, yMax = -1., 0.
n = 100
i = 0
random.seed( 314159 )
while( i < n ) :
    r = random.random()
    x = xMin * r + xMax * ( 1. - r )
    r = random.random()
    y = yMin * r + yMax * ( 1. - r )
    print 'f.setValue( %.18e, %.18e )' % ( x, y )
    f.setValue( x, y )
    i += 1
    xm = xMin - 1
    for x, y in f :
        if( x < xm ) : raise Exception( 'x = %16e < xm = %16e: i = %d' % ( x, xm, i ) )
        xm = x

length = len( f )
f.showInteralStructure( printPointersAsNull = True )
