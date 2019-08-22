#! /usr/bin/env python

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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
This file compares the results of evaluate from an endl2dmathClasses.endl2dmath and pointwiseXY_C
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
xMin, xMax = data.domainMin( ), data.domainMax( )
r = random.random()
x = xMin * r + xMax * ( 1. - r )
print '# time =', time.clock( ) - t0
for i in xrange( n ) :
    r = random.random( )
    x = xMin * r + xMax * ( 1. - r )
    y1 = data.evaluate( x )
    y2 = f.evaluate( x )
    if( abs( y2 - y1 ) > eps * ( y1 + y2 ) ) : print "#ERROR:", i, x, y1, y2
    if( ( i + 1 ) % np == 0 ) : print "%s of %s" % ( i + 1, n )
print '# time =', time.clock( ) - t0
