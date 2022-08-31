# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU_float, PQU_uncertainty, PQU

a = PQU_float( 1.234, 4 )
b = PQU_float( 1.234, 4, True )
print(a)
print(a.info( significantDigits = 15 ))
print()
print(b)
print(b.info( significantDigits = 15 ))
print()
u = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStylePlusMinus, .12, isPercent = True )
print(u)
print(u.info( significantDigits = 15 ))

c = PQU( "1e-12 +/- 1e-16" )
print()
print(c)
print(c.info( significantDigits = 15 ))

c = PQU( "1.0000e-12(1)" )
print()
print(c)
print(c.info( significantDigits = 15 ))
