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

import sys
f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

f = open( sys.argv[2], 'w' )

MF = int( sys.argv[3] )
if( sys.argv[4] == '-' ) :
    MT = None
else :
    MT = int( sys.argv[4] )

skipLines, extraSkipLines = 2, 0
if( ( MF == 2 ) and ( MT == 151 ) ) : skipLines = 3
for l in ls :
    try :
        fileMF, fileMT = int( l[70:72] ), int( l[72:75] )
    except :
        f.write( l )
        continue
    if( MF == fileMF ) :
        if( MT is None ) :
            extraSkipLines = skipLines
        elif( MT == fileMT ) :
            extraSkipLines = skipLines
    if( extraSkipLines > 0 ) :
        extraSkipLines -= 1
    else :
        f.write( l )
f.close( )
