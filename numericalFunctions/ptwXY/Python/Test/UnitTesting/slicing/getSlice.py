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

import os, copy
import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

def checkSlicing( XYs_base, i1, i2, doPrint = False ) :

    XYs = copy.copy( XYs_base )
    pXYs = pointwiseXY_C.pointwiseXY_C( data = XYs, initialSize = 20 )
    sXYs = XYs[i1:i2]
    if( doPrint ) :
        print
        print i1, i2
        for xy in XYs : print xy
        print
        print sXYs
    spXYs = pXYs[i1:i2]
    if( doPrint ) : spXYs
    if( len( sXYs ) != len( spXYs ) ) : raise Exception( "%s: len( sXYs ) = %d != len( spXYs ) = %d: index1 = %d, i2 = %d" % \
        ( __file__, len( sXYs ), len( spXYs ), i1, i2 ) )
    for i, xy in enumerate( sXYs ) :
        if( xy[0] != spXYs[i][0] ) : raise Exception( "%s: difference at index = %d: %e %e" % ( __file__, i, xy[0], spXYs[i][0] ) )
    
XYs = [ [ float( x ), float( x )**2 ] for x in xrange( 12 ) ]

checkSlicing( XYs, 4, 8 )
checkSlicing( XYs, -4, 8 )
checkSlicing( XYs, -4, -2 )
checkSlicing( XYs, 4, -2 )
checkSlicing( XYs, 4, -8 )
checkSlicing( XYs, -4 - 7 * len( XYs ), 8 )
checkSlicing( XYs, 3, -4 - 7 * len( XYs ) )
