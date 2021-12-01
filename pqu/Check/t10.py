# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import PQU

def func( a1, a2 ) :

    def func2( op, a1, a2 ) :

        print('\na1 %s a2' % op)
        try :
            b = eval( 'a1 %s a2' % op )
            print(b)
            print(b.info( significantDigits = 15 ))
        except ZeroDivisionError :
            print('Divide by zero')
        except :
            raise

    print('\n')
    print(a1)
    print(a1.info( significantDigits = 15 ))
    if( a1 is not a2 ) :
        print(a2)
        print(a2.info( significantDigits = 15 ))
    func2( '+', a1, a2 )
    func2( '-', a1, a2 )
    func2( '*', a1, a2 )
    func2( '/', a1, a2 )

a1 = PQU( 0, 'eV', .13 )
a2 = PQU( 1.2, 'eV', .2 )
print(1 / a2)
func( a1, a1 )
func( a2, a2 )
func( a1, a2 )
func( a2, a1 )
