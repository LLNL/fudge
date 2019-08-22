import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

print
a2 = PQU( "2.300000000000(1) MeV" )
print a2
print a2.info( significantDigits = 15 )
a2 = PQU( "2.300000000001(1) MeV" )
print a2
print a2.info( significantDigits = 15 )

print
a2 = PQU( "2.300000000003(1)", "MeV" )
print a2
print a2.info( significantDigits = 15 )


print
l = PQU( 10., 'm' )
big_l = PQU( 10., 'km' )
sum_l = big_l + l
print sum_l
print l.info( significantDigits = 15 )
print big_l.info( significantDigits = 15 )
print sum_l.info( significantDigits = 15 )

print
E = PQU( 1000, 'MeV/c**2' )                 # This is similar to the prior one.
print E
print E.info( significantDigits = 15 )
kg = E.inUnitsOf( 'kg' )
print kg
print kg.info( significantDigits = 15 )
