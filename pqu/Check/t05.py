# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import valueOrPQ, _findUnit, PQU, floatToShortestString

doRaise = False

def a( v, unitFrom, unitTo ) :

    def b( v, unitFrom, unitTo, asPQU ) :

        vv = valueOrPQ( v, unitFrom, unitTo, asPQU )

        if( isinstance( vv, PQU ) ):
            print(vv)
        else:
            print(floatToShortestString(vv))

    print()
    print('---------')
    try :
        b( v, unitFrom, unitTo, False )
    except :
        if( doRaise ) : raise
        print("FAILED:", v, "'%s'" % unitFrom, "'%s'" % unitTo)
        return
    b( v, unitFrom, unitTo, True )

dl = _findUnit( "" )
mm = _findUnit( "mm" )
pqu_m = PQU( "3.14 m" )

a( 3.14, None, None )
a( 3.14, "",   None )
a( 3.14, None, "" )
a( 3.14, "",   "" )

a( 3.14, dl,   None )
a( 3.14, None, dl )
a( 3.14, dl,   dl )

a( 3.14, "m", None )
a( 3.14, "m", "cm" )
a( 3.14, "cm", "m" )

doRaise = True
a( pqu_m, None, None )
a( pqu_m, None, "cm" )
a( pqu_m, None, "km" )
try :
    a( pqu_m, "", None )
    print('========= FAILED =========')
except :
    pass
