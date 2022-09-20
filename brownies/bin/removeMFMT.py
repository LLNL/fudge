#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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
