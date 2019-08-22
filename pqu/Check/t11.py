import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import floatToShortestString

def f( s ) :

    print s, floatToShortestString( float( s ), includeSign = True, significantDigits = 12, favorEFormBy = -12 )

f( '-inf' )
f( 'inf' )
f( '+inf' )
f( 'nan' )
f( -1e10 )
f( 1e10 )
