import sys
sys.path.insert( 0, '../../' )

from pqu import PQU

def add( pqu1, pqu2 ) :

    print
    print pqu1
    print pqu1.info( significantDigits = 15 )
    print pqu2
    print pqu2.info( significantDigits = 15 )
    sum = pqu1 + pqu2
    print sum
    print sum.info( significantDigits = 15 )

def loop( v1, v2 ) :

    pqu1 = PQU.PQU( v1 + 'm' )
    for prefix in PQU._prefixes :
        add( pqu1, PQU.PQU( v2 + prefix[0] + 'm' ) )
        add( PQU.PQU( v2 + prefix[0] + 'm' ), pqu1 )

loop( '100.1(4)', '100.1(4)' )
loop( '550.1(4)', '10.1(4)' )
loop( '9999.1(4)', '10.1(4)' )
loop( '9999.1(4)', '99.1(1)' )
