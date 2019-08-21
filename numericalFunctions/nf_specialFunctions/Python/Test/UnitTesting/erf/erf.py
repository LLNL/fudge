#
# <<BEGIN-copyright>>
# <<END-copyright>>
#
import sys, os
sys.path.insert( 0, '../../../../../lib' )
import nf_specialFunctions_C

options = []
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )
for argv in sys.argv[1:] : options += [ argv ]

if( '-e' in options ) : print __file__

x = 1e-4
while( x < 5 ) :
    erf  = nf_specialFunctions_C.erf( x )
    erfc = nf_specialFunctions_C.erf( x, True )
    print x, erf, erfc, erf+ erfc
    x *= 1.2
