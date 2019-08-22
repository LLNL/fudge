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

import sys, os
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )
from fudge.legacy.converting import endfFormats

f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

f = open( sys.argv[2], 'w' )

END66Characters = endfFormats.endfSENDLine( 0, 0 )[:66]
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
