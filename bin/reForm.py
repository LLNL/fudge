# <<BEGIN-copyright>>
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
