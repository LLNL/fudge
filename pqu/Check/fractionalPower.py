# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU

def f(factor, unit):

    print()
    print(factor, unit)
    quantity = PQU.PQU(factor, unit)
    print(quantity)
    print(PQU.PQU(quantity))
    print(quantity**(1/2))
    print(quantity * PQU.PQU(3, 'MeV'))
    print(quantity * PQU.PQU(3, 'MeV**(3/2)'))
    print(quantity * PQU.PQU(3, 'MeV**2'))

f(10, 'b * MeV**(-1/2)')
f(10, '(b * MeV**(-1/2))**(1/3)')
f(10, 'b * MeV**(1/2)')
f(10, 'b * MeV**2')
f(10, '(b * MeV**2)**(1/3)')
