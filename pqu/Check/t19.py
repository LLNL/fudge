# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU_float, PQU

def f( s, unc = None ) :

    print()
    print(s)
    print(unc)
    a = PQU( s, uncertainty = unc )
    print(a)
    print(a.info( significantDigits = 15 ))
    a.truncate( True, True )
    print(a)
    print(a.info( significantDigits = 15 ))

f( '14850. +/- 1000. eV/mm' )
f( '14850. +/- 100. eV/mm' )
f( '1485. +/- 10. eV/mm' )
f( '1485. +/- 1. eV/mm' )
f( '1485. +/- 1.1 eV/mm' )
f( '1485. +/- 0.1 eV/mm' )
f( '1485. +/- 0.11 eV/mm' )
f( '1485. +/- 0.01 eV/mm' )
f( '1485. +/- 0.011 eV/mm' )
f( '14850. +/- 1780. eV/mm' )
f( '14850. +/- 178. eV/mm' )
f( '14850. +/- 17.8 eV/mm' )
f( '14850. +/- 1.78 eV/mm' )
f( '14850. +/- 0.178 eV/mm' )

print()
f( '1485. +/- 1e1 eV/mm' )
f( '1485. eV/mm', unc = PQU_float( 10, 1 ) )
