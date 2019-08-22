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
sys.path.insert(0, os.path.dirname( binDir ) )
from fudge.legacy.converting import endfFileToGNDMisc, endfFormats

f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

f = open( sys.argv[2], 'w' )

r6 = xrange( 0, 66, 11 )
s = ' 0.9999'
n = len( s )
for i, l in enumerate( ls ) :
    MF = int( l[70:72] )
    if( MF == 1 ) :
        f.write( l )
    else :
        for i6 in r6 :
            if( l[i6:i6+n] == s ) :
                nl = ''
                for i6 in r6 :
                    datum = l[:11]
                    if( datum[:n] == s ) : datum = endfFormats.floatToFunky( endfFileToGNDMisc.funkyFloatStringToFloat( 0, datum ) )
                    nl += datum
                    l = l[11:]
                l = nl + l
        f.write( l )
f.close( )
