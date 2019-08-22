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

import sys, math, random
sys.path.insert( 0, '../../Utilities' )
sys.path.insert( 0, '../../../../../lib' )

import pointwiseXY_C
import utilities
options = utilities.getOptions( __file__ )

xData = [ 0.5 * i - 0.145 for i in xrange( 51 ) ]
yData = [ math.sin( x ) for x in xData ]

xys = pointwiseXY_C.pointwiseXY_C( [ xData, yData ], initialSize = len( xData ), overflowSize = 10, dataForm = 'XsAndYs' )

n1 = len( xData )
for i1 in xrange( n1 ) :
    n2 = len( xData )
    i2 = random.randint( 0, n2 - 1 )
    x1 = xData.pop( i2 )
    y1 = yData.pop( i2 )
    x2, y2 = xys.pop( i2 )
    if( x1 != x2 ) : raise Exception( 'x1 = %.17e != x2 = %.17e, length left = %d' % ( x1, x2, len( xData ) ) )
    if( y1 != y2 ) : raise Exception( 'y1 = %.17e != y2 = %.17e, length left = %d' % ( y1, y2, len( xData ) ) )
    if( ( len( xData ) != len( yData ) ) or ( len( xData ) != len( xys ) ) ) : 
        raise Exception( 'len( xData ) = %d != len( yData ) = %d != len( xys ) = %d' % ( len( xData ), len( yData ), len( xys ) ) )

if( len( xData ) != 0 ) : raise Exception( 'len( xData ) = %d != 0' % len( xData ) )
if( len( yData ) != 0 ) : raise Exception( 'len( yData ) = %d != 0' % len( yData ) )
if( len( xys ) != 0 ) : raise Exception( 'len( xys ) = %d != 0' % len( xys ) )
