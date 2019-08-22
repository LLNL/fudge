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
sys.path.insert( 0, '../../Utilities' )

import os, math
import pointwiseXY_C
import utilities

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

CPATH = '../../../../Test/UnitTesting/toFunction'

os.system( 'cd %s; createFromFunction -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def getData( ls ) :

    ls = utilities.skipBlankLines( ls )
    if( len( ls ) == 0 ) : return( None, None, None, None, None, None )
    label, accuracy = None, None
    label, ls = ls[0].strip( ), ls[1:]
    if( label != '# Xs' ) : raise Exception( '%s: missing line "# Xs": found "%s"' % ( __file__, label ) )
    ls, Xs = utilities.getXData( ls, '=' )
    ls, accuracy = utilities.getDoubleValue( 'accuracy', ls )
    ls, biSectionMax = utilities.getDoubleValue( 'biSectionMax', ls )
    ls, data = utilities.getXYData( ls )
    return( ls, label, accuracy, biSectionMax, Xs, pointwiseXY_C.pointwiseXY_C( data, accuracy = accuracy ) )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.7g' % v1, '%.7g' % v2
    if( sv1 != sv2 ) : raise Exception( '%s: values %s %s diff at %d for label = %s' % ( __file__, v1, v2, i, label ) )

def f_x_Sin_xx( x, args ) :

       return( x * math.sin( x * x ) )

def createFromFunction( label, accuracy, biSectionMax, Xs, original ) :

    data = pointwiseXY_C.createFromFunction( Xs, f_x_Sin_xx, None, accuracy, biSectionMax, checkForRoots = 1 )
    if( len( data ) != len( original ) ) : raise Exception( '%s: len( data ) = %d != len( original ) = %d for label = "%s"' % \
        ( __file__, len( data ), len( original ), label ) )

while( 1 ) :
    ls, label, accuracy, biSectionMax, Xs, original = getData( ls )
    if( ls is None ) : break
    createFromFunction( label, accuracy, biSectionMax, Xs, original )
