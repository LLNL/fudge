import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

t = PQU( "314159.", 's')
print t.info( significantDigits = 15 )
s = t.inUnitsOf( 'd', 'h', 'min', 's' )
for i in s : print '\n', i, '\n', i.info( significantDigits = 15 )

