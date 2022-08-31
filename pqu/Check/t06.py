# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU, Parsers

def f( s ) :

    print()
    print('-------------')
    a = Parsers.parsePQUString( s )
    print(s)
    v = PQU( a[0], a[1], a[2] )
    print(v)
    print(v.info( significantDigits = 15 ))
    print()
    p = 3.2 * v
    print(p)
    print(p.info( significantDigits = 15 ))

    print()
    p = v * "3.2 MeV"
    print(p)
    print(p.info( significantDigits = 15 ))

    print()
    p = v * "3.2%"
    print(p)
    print(p.info( significantDigits = 15 ))

f( '12' )
f( '12%' )
f( '12 J' )
f( '12(3)' )
f( '12(3) J' )
f( '12 +/- 3' )
f( '12 +/- 3 J' )
f( '12(3)%' )
f( '12 +/- 3%' )
f( '(12 +/- 3)%' )
f( '12 +/- 3% J' )
