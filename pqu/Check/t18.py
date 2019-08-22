import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import pqu_uncertainty, PQU

def changeUnc( s ) :

    print '============'
    b = PQU( s )
    print b
    print b.info( significantDigits = 15 )
    print
    b.changeUncertaintyPercent( not( b.uncertainty.isPercent( ) ) )
    print b
    print b.info( significantDigits = 15 )
    print
    b.changeUncertaintyPercent( not( b.uncertainty.isPercent( ) ) )
    print b
    print b.info( significantDigits = 15 )

changeUnc( "2.300000 +/- 4e-6 MeV" )
changeUnc( "2.30e-12 +/- 0.41% MeV" )
changeUnc( "8.30e+12 +/- 0.041% MeV" )
