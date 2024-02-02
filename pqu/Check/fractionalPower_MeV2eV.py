# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU
PQU.Debugging = True

def to_eV(quantity):

    s = str(quantity.unit).replace('MeV', 'eV')
    qq = PQU.PQU(quantity).convertToUnit(s)
    print('to_eV:')
    print('    quantity:  ',  quantity.toString(significantDigits=12))
    print('    unit in eV:', s)
    print('    in eV:     ', qq.toString(significantDigits=12))

def f(factor, unit):

    print()
    print(factor, unit)
    quantity = PQU.PQU(factor, unit)
    print(quantity.toString(significantDigits=12))
    to_eV(quantity)
    to_eV(quantity**(1/2))

f(10, 'b * MeV**(-1/2)')
f(10, '(b * MeV**(-1/2))**(1/3)')
f(10, 'b * MeV**(1/2)')
f(10, 'b * MeV**2')
