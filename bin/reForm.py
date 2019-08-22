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

import sys, os
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )
import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule

f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

f = open( sys.argv[2], 'w' )

END66Characters = endfFormatsModule.endfSENDLine( 0, 0 )[:66]
fillZero = '%11d' % 0

index = 1
ls = ls[1:]
f.write( '                                                                     1 0  0    0\n' )
newLines = []
for i, l in enumerate( ls ) :
    if( l[75:80] == '99999' ) :
        newLines.append( l[:-1] )
        index = 0
    elif( l[70:75] == ' 0  0' ) :
        newLines.append( l[:-1] )
        index = 0
    else :
        newLines.append( "%s%5d" % ( l[:75], index ) )
    index += 1

for l in newLines :
    if( int( l[66:70] ) == -1 ) :
        l = END66Characters + l[66:]
    if( int( l[72:75] ) == 0 ) :
        l = END66Characters + l[66:]
    elif( l[71:75] != '1451' ) :
        if( l.find( '.' ) < 0 ) :
            s = l[:66].rstrip( )
            if( len( s ) < 56 ) :
                n = ( len( s ) + 10 ) / 11
                s += ( n * 11 - len( s ) ) * ' '
                for i in xrange( n, 6 ) : s += fillZero
                l = s + l[66:]
    f.write( l + '\n' )
f.close( )
