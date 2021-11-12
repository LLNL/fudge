# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

a = PQU( "2.31(1) MeV / cm" )
print(a)
print()
print(a**6)           # This prints "152.(3) MeV**6/cm**6" which is correct.
print((a**6).info( significantDigits = 15 ))
print()
print(a**7)           # This prints "351.(1) MeV**7/cm**7" which should be "351.(9) MeV**7/cm**7" as the uncertainty is 8.88.
print((a**7).info( significantDigits = 15 ))
print()
print(a**12)
print((a**12).info( significantDigits = 15 ))

a2 = PQU( 2.3, "MeV" )
a = PQU( "2.3000000000002(1) MeV" )
print(a2.info( significantDigits = 15 ))

print(a == a2)
print(a != a2)
print(a <= a2)
print(a < a2)
print(a >= a2)
print(a > a2)

print()
a = PQU( "2.31(1) MeV / cm" )
print(a)
a3 = a**7
print(a3)
print(a**6)

print()

a = PQU( 2.3, "MeV" )
a1 = PQU( "2.3 +/- 1e-14 MeV" )
a2 = PQU( "2.30000000000003(1)", "MeV" )
b = PQU( 3e2, "m" )
c = PQU( 1, "J" )
print('a = ', a)
print('a1 = ', a1)
print('a2 = ', a2)
print('b = ', b)
print('c = ', c)
print(c.inBaseUnits( ))
print()
print(a.inBaseUnits( ))
print(a.getValueAs( 'm**2*kg/s**2' ))
print(a * b)
print(a)
a.convertToUnit( 'J' )
print(a)
print(a * b)
print(a.inBaseUnits( ))
print(b.inBaseUnits( ))

print()
c = PQU( 342, "", ".12%" )
print(c)
print(c.info( significantDigits = 15 ))
print(100 + c)
