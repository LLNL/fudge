import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import pqu_float, pqu_uncertainty, PQU

a = pqu_float( 1.234, 4 )
b = pqu_float( 1.234, 4, True )
print a
print a.info( significantDigits = 15 )
print
print b
print b.info( significantDigits = 15 )
print
u = pqu_uncertainty( pqu_uncertainty.pqu_uncertaintyStylePlusMinus, .12, isPercent = True )
print u
print u.info( significantDigits = 15 )

c = PQU( "1e-12 +/- 1e-16" )
print
print c
print c.info( significantDigits = 15 )

c = PQU( "1.0000e-12(1)" )
print
print c
print c.info( significantDigits = 15 )
