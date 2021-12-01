# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

a = PQU( '10.0(2) m' )
b = PQU( '22.0(3) m' )
r = b / a
print()
print(a)
print(a.info( significantDigits = 15 ))
print()
print(b)
print(b.info( significantDigits = 15 ))
print()
print(r)
print(r.info( significantDigits = 15 ))
