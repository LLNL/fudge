# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU
PQU.Debugging = True

def Print(quantity):

    print()
    print('   ', quantity)
    print('   ', quantity.unit.symbols)
    print('   ', quantity.unit.powers)

def f(factor, unit, power):

    print()
    print(factor, unit, power)
    quantity = PQU.PQU(factor, unit)
    Print(quantity)
    Print(quantity**power)
    Print(quantity * PQU.PQU(1, 'MeV**(1/2)'))
    Print(quantity / PQU.PQU(1, 'MeV**(1/2)'))

f(10, 'b * MeV**(-1/2)', 2)
f(10, '(b * MeV**(-1/2))**(1/3)', 6)
f(10, 'b * MeV**(1/2)', 2)
f(10, 'b * MeV', 1)

a = PQU.PQU(1, 'MeV**(1/2)')
b = PQU.PQU(1, 'MeV')

print(a, b, a.isPhysicalUnitOf(a, b))
print(a**2, b, a.isPhysicalUnitOf(a**2, b))
print(a, b**(1/2), a.isPhysicalUnitOf(a, b**(1/2)))
