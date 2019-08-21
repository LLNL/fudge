# <<BEGIN-copyright>>
# <<END-copyright>>

import sys
f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

f = open( sys.argv[2], 'w' )

MTs = [ int( MT ) for MT in sys.argv[3:] ]

skipLines, extraSkipLines = 2, 0
priorMF = None
for l in ls :
    try :
        fileMF, fileMT = int( l[70:72] ), int( l[72:75] )
    except :
        f.write( l )
        continue
    if( priorMF == fileMF ) :
        if( fileMT in MTs ) : extraSkipLines = skipLines
    if( extraSkipLines > 0 ) :
        extraSkipLines -= 1
    else :
        f.write( l )
    priorMF = fileMF
f.close( )
