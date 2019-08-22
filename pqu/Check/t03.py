import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

print PQU( "1.23", "m", "12e-2" )
print PQU( "1.23 m", uncertainty = "12e-2" )
print
print PQU( "1.23", "m", "12%" )
print PQU( "1.23 m", uncertainty = "12%" )
print PQU( "1.23 +/- 12%", "m" )
print PQU( "1.23+/-12% m" )
