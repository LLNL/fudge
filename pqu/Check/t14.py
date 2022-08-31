# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import floatToShortestString

value = 100
for trimZeros in [ 0, 1 ] :
    for keepPeriod in [ 0, 1 ] :
        for favorEFormBy in [ -5, 5 ] :
            for significantDigits in range( 3 ) :
                a = floatToShortestString( value, significantDigits, trimZeros, keepPeriod, favorEFormBy )
                print(trimZeros, keepPeriod, "%2s" % favorEFormBy, significantDigits, "<%s>" % a)
print(floatToShortestString( 1.234e-9, 12 ))
