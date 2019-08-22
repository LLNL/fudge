from pqu.PQU import PQU

v = PQU( '1.234(56) eV' )
print v

u = v.getValueAs( 'MeV' )
print type( u ), u

u = v.getValueAs( 'MeV', asPQU = True )
print type( u ), u

u = v.getUncertaintyValueAs( )
print type( u ), u

u = v.getUncertaintyValueAs( 'MeV' )
print type( u ), u

u = v.getUncertaintyValueAs( asPQU = True )
print type( u ), u

u = v.getUncertaintyValueAs( 'MeV', asPQU = True )
print type( u ), u
