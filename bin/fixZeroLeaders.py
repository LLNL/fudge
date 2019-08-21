# <<BEGIN-copyright>>
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
