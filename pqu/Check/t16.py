# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

c = PQU( "0.0000 +/- 3e-4" )
print()
print(c)
print(c.info( significantDigits = 15 ))

c = PQU( "0.0000(3)" )
print()
print(c)
print(c.info( significantDigits = 15 ))
