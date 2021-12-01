# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import pqu_uncertainty, PQU

def changeUnc( s ) :

    print('============')
    b = PQU( s )
    print(b)
    print(b.info( significantDigits = 15 ))
    print()
    b.changeUncertaintyStyle( pqu_uncertainty.pqu_uncertaintyStyleParenthesis )
    print(b)
    print(b.info( significantDigits = 15 ))
    print()
    b.changeUncertaintyStyle( pqu_uncertainty.pqu_uncertaintyStylePlusMinus )
    print(b)
    print(b.info( significantDigits = 15 ))

changeUnc( "2.300098(4) MeV" )
changeUnc( "2.300000 +/- 4e-6 MeV" )
changeUnc( "2.300e12(4) MeV" )
changeUnc( "2.30e-12 +/- 4e-18 MeV" )
changeUnc( "2.300e-12(4)%" )
changeUnc( "(2.3e-12 +/- 4e-18)%" )


styles = [ pqu_uncertainty.pqu_uncertaintyStyleNone, pqu_uncertainty.pqu_uncertaintyStylePlusMinus, pqu_uncertainty.pqu_uncertaintyStyleParenthesis ]
pqus = { pqu_uncertainty.pqu_uncertaintyStyleNone : PQU( "2.300098" ),
        pqu_uncertainty.pqu_uncertaintyStylePlusMinus : PQU( "2.300098 +/- 4e-6 MeV" ),
        pqu_uncertainty.pqu_uncertaintyStyleParenthesis : PQU( "2.300098(4) MeV" ) }

for style in styles :
    pqu = pqus[style]
    for toStyle in styles :
        try :
            pqu.changeUncertaintyStyle( toStyle )
            if( style != toStyle ) :
                if( ( style == pqu_uncertainty.pqu_uncertaintyStyleNone ) or ( toStyle == pqu_uncertainty.pqu_uncertaintyStyleNone ) ) :
                    print('========= FAILED ========= 1) changeUncertaintyStyle: %s to %s' % ( style, toStyle ))
        except :
            if( style == toStyle ) : print('========= FAILED ========= 2) changeUncertaintyStyle: %s to %s' % ( style, toStyle ))
