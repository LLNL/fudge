# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import floatToShortestString

def f( s ) :

    print(s, floatToShortestString( float( s ), includeSign = True, significantDigits = 12, favorEFormBy = -12 ))

f( '-inf' )
f( 'inf' )
f( '+inf' )
f( 'nan' )
f( -1e10 )
f( 1e10 )
