# <<BEGIN-copyright>>
# <<END-copyright>>

import sys
f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

f = open( sys.argv[2], 'w' )
for l in ls : f.write( l[:75] + '\n' )
f.close( )
