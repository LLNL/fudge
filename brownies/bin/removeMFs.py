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

MFs = [ int( MF ) for MF in sys.argv[3:] ]

extraSkipLines = 0
for l in ls :
    try :
        fileMF = int( l[70:72] )
    except :
        f.write( l )
        continue
    if( fileMF in MFs ) : extraSkipLines = 2
    if( extraSkipLines > 0 ) :
        extraSkipLines -= 1
    else :
        f.write( l )
f.close( )
