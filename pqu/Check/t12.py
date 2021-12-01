# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

t = PQU( "314159.", 's')
print(t.info( significantDigits = 15 ))
s = t.inUnitsOf( 'd', 'h', 'min', 's' )
for i in s : print('\n%s \n%s' % (i, i.info( significantDigits = 15 )))

