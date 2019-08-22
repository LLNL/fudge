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

import sys, os
sys.path.insert( 0, '../../Utilities' )
sys.path.insert( 0, '../../../../../lib' )

import pointwiseXY_C
import utilities
options = utilities.getOptions( __file__ )

CPATH = '../../../../Test/UnitTesting/interpolation'
fIn = open( os.path.join( CPATH, 'strings.out' ) )
strings = []
for line in fIn.readlines( ) :
    if( line == '\n' ) : continue
    strings.append( line.split( '=' )[1].strip( ) )
fIn.close( )

def compareStrings( stringC, stringPy ) :

    _stringC = stringC[1:-1]
    if( _stringC != stringPy ) : raise Exception( '_stringC = "%s" not equal to stringPy = "%s"' % ( _stringC, stringPy ) )

def check( ptwXY2, interpolation, index ) :

    ptwXY = pointwiseXY_C.pointwiseXY_C( [ [ 1, 1 ], [ 10, 10 ] ], interpolation = interpolation )
    string = ptwXY.getInterpolation( )
    compareStrings( strings[index], string )
    if( '-v' in options ) : print string

    ptwXY2 = ptwXY.copy( )
    string = ptwXY.getInterpolation( )
    compareStrings( strings[index+1], string )
    if( '-v' in options ) : print string
    if( '-v' in options ) : print string

    return( ptwXY )

ptwXY2 = pointwiseXY_C.pointwiseXY_C( [ [ 1, 1 ], [ 10, 10 ] ], interpolation = 'charged-particle' )
ptwXY2 = check( ptwXY2, 'lin,lin', 0 )
ptwXY2 = check( ptwXY2, 'lin,log', 3 )
ptwXY2 = check( ptwXY2, 'log,lin', 6 )
ptwXY2 = check( ptwXY2, 'log,log', 9 )
ptwXY2 = check( ptwXY2, 'flat', 12 )
ptwXY2 = check( ptwXY2, 'charged-particle', 15 )
