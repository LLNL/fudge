# <<BEGIN-copyright>>
# <<END-copyright>>

import sys, os

file = sys.argv[1]

status = os.system( 'diff v %s.dat > /dev/null' % file )
if( status != 0 ) : print "ERROR for %s" % file
