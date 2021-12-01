#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import brownies.legacy.toENDF6.endfFormats as endfFormatsModule
from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc

f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

f = open( sys.argv[2], 'w' )

r6 = range( 0, 66, 11 )
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
                    if( datum[:n] == s ) : datum = endfFormatsModule.floatToFunky(endfFileToGNDSMisc.funkyFloatStringToFloat(0, datum))
                    nl += datum
                    l = l[11:]
                l = nl + l
        f.write( l )
f.close( )
