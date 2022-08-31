# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

def f( s, shouldRaise = False ) :

    print()
    try :
        v = PQU( s )
    except :
        if( not( shouldRaise ) ) :
            print(s)
            print('========= FAILED =========')
        return
    if( shouldRaise ) :
        print(s)
        print('========= FAILED =========')
    print(v)
    print(float( v ), v.getUncertaintyValueAs( ))

f( '333 +/- 1' )
f( '333 +/- 1%' )
f( '333% +/- 1', True )   # Should not be allowed
f( '333% +/- 1%', True )  # Should not be allowed
f( '333%(1)', True )      # Should not be allowed
f( '333(1)%' )

print()
a = PQU( '10.0(2) m' )
b = PQU( '11.0(3) m' )
r = b / a
print(a)
print(a.info( significantDigits = 15 ))
print(b)
print(b.info( significantDigits = 15 ))
print(r)
print(r.info( significantDigits = 15 ))
