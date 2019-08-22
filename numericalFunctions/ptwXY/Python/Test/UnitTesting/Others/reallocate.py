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

import sys
sys.path.insert( 0, '../../../../../lib' )

import os
import pointwiseXY_C

verbose = False

options = []
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )
for argv in sys.argv[1:] : options += [ argv ]

if( '-e' in options ) : print __file__
if( '-v' in options ) : verbose = True

def printIfVerbose( ) :

    if( not( verbose ) ) : return
    print
    print '# length =' , len( data ), '  allocated size =', data.allocatedSize( )
    print '# overflow length =', data.overflowLength( ), '  overflow allocated size =', data.overflowAllocatedSize( )

if( verbose ) :
    doc = pointwiseXY_C.pointwiseXY_C.reallocatePoints.__doc__.split( '\n' )
    for d in doc : print "#", d
    print "#"
    doc = pointwiseXY_C.pointwiseXY_C.reallocateOverflowPoints.__doc__.split( '\n' )
    for d in doc : print "#", d

data = [ [ i, i * 10 + 3 ] for i in xrange( 51 ) ]
data = pointwiseXY_C.pointwiseXY_C( data, initialSize = 10, overflowSize = 10 )

printIfVerbose( )

data.reallocatePoints( 1000 )
data.reallocateOverflowPoints( 25 )
if( data.allocatedSize( ) != 1000 ) : raise Exception( '%s: data.allocatedSize( ) = %s != 1000' % ( __file__, data.allocatedSize( ) ) )
if( data.overflowAllocatedSize( ) != 25 ) : raise Exception( '%s: data.overflowAllocatedSize( ) = %s != 25' % i\
    ( __file__, data.overflowAllocatedSize( ) ) )

printIfVerbose( )

for i in xrange( 103 ) : data.setValue( i + 0.5, i * 10 + 3.5 )
printIfVerbose( )
