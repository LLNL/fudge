#! /bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
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
ptwXY2 = check( ptwXY2, 'lin-lin', 0 )
ptwXY2 = check( ptwXY2, 'log-lin', 3 )
ptwXY2 = check( ptwXY2, 'lin-log', 6 )
ptwXY2 = check( ptwXY2, 'log-log', 9 )
ptwXY2 = check( ptwXY2, 'flat', 12 )
ptwXY2 = check( ptwXY2, 'charged-particle', 15 )
