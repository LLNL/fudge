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

f = open( '../../../../Test/UnitTesting/gammaFunctions/gammaTest.dat' )
ls = f.readlines( )
f.close( )

errors = 0
for l in ls :
    if( ';' in l ) : break
    if( ',' not in l ) : continue
    x, f = l.split( '{' )[1].split( '}' )[0].split( ',' )
    x = float( x )
    f = float( f )
    gamma = nf_specialFunctions_C.gamma( x )
    r = 1
    if( f != 0 ) : r = gamma / f - 1
    if( abs( r ) > 3e-15 ) :
        print x, f, gamma, r
        errors += 1

if( errors > 0 ) : raise Exception( "FAILED: %s" % __file__ )
