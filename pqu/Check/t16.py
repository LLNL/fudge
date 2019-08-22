import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

c = PQU( "0.0000 +/- 3e-4" )
print
print c
print c.info( significantDigits = 15 )

c = PQU( "0.0000(3)" )
print
print c
print c.info( significantDigits = 15 )
