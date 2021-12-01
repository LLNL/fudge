# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

print("2.32(12)")
p1 = PQU( "2.32(12)" )
print(p1)

print()
print("2.32(123)")
try :
    p1 = PQU( "2.32(123)" )
except :
    print(p1)
    print(p1.info( significantDigits = 15 ))
    print('========= FAILED =========')

print()
print("2.3237(1234)")
try :
    p1 = PQU( "2.3237(1234)" )  # print "2.3237(12)" and not "2.32(12)".
except :
    print(p1)
    print(p1.info( significantDigits = 15 ))
    print('========= FAILED =========')

print()
print("2.32374321(1234)")
try :
    p1 = PQU( "2.32374321(1234)" )  # print "2.3237(12)" and not "2.32(12)".
except :
    print(p1)
    print(p1.info( significantDigits = 15 ))
    print('========= FAILED =========' )
