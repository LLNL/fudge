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
This stuff should be repeatable. Simple test of operators '+', '-', '*', 'neg', 'abs', 'pow' and union method.
"""
import sys
sys.path.append( '../../../lib' )

import random
import pointwiseXY_C

d = pointwiseXY_C.pointwiseXY_C( initialSize = 100, overflowSize = 31 )

d.setValue( 1, 11. )
d.setValue( 3.14, 3.12 )
d.setValue( 7.66, -1.5 )
d.setValue( 3.66, -1.6 )
d.setValue( 2.66,  1.7 )

print '---- d ----'
print d

e = d
del e[2]

print '---- del d[2] ----'
print d

d = pointwiseXY_C.pointwiseXY_C( initialSize = 100, overflowSize = 31 )

d.setValue( 1, 11. )
d.setValue( 3.14, 3.12 )
d.setValue( 7.66, -1.5 )
d.setValue( 3.66, -1.6 )
d.setValue( 2.66,  1.7 )
e[2] = d[2]

print '---- e ----'
print e

e = pointwiseXY_C.pointwiseXY_C( initialSize = 5, overflowSize = 3 )
xMin, xMax = 0.5, 8
yMin, yMax = -1.2, 10.
random.seed( 314159 )
for i in xrange( 7 ) :
    r = random.random( )
    x = xMin * r + xMax * ( 1. - r )
    r = random.random( )
    y = yMin * r + yMax * ( 1. - r )
    e.setValue( x, y )

print '---- d ----'
print d

print '---- neg of d ----'
print -d

print '---- abs of d ----'
print abs( d )

print '---- d.pow( 3 ) ----'
print d.__pow__( 3 )

print '---- e ----'
print e

print '---- d.union( e ) ----'
u = d.union( e )
print u

print '---- union and map e onto d ----'
u = d.union( e, True )
print u

print '---- union and map d onto e ----'
u = e.union( d, True )
print u

print '---- union and map and trim of e onto d ----'
u = d.union( e, True, True )
print u

e[ 0] = [ e[ 0][0], 0. ]        # needed to make e's domain mutual with d
e[-1] = [ e[-1][0], 0. ]

print '---- e ----'
print e

print '---- d + e ----'
print d + e

print '---- e + d ----'
print e + d

print '---- d + 1.41 ----'
print d + 1.41

print '----  1.41 + d ----'
print 1.41 + d

print '---- 1.41 - d ----'
print 1.41 - d

print '---- d - 1.41 ----'
print d - 1.41

print '---- d * 1.41 ----'
print d * 1.41

print '---- 1.41 * d ----'
print 1.41 * d

print '---- d + d ----'
print d + d

print '---- d - d ----'
print d - d

print '---- d * d ----'
print d * d

print '---- d * d * d * d ----'
print d * d * d * d

print '---- d**4 ----'
print d**4

print '---- d**-5 ----'
print d**-5
