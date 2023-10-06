#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

summaryDocStringPoPs = '''Parses the bit components of an intid value and prints results.'''

description = '''Parses the bit components of an intid value and prints results.'''

import argparse

from PoPs import intId as intIdModule

parser = argparse.ArgumentParser(description=description)
parser.add_argument('intids', type=int, nargs='*',          help='''The list particle intids to parse.''')

if __name__ == '__main__':
    args = parser.parse_args()

    intidDB = intIdModule.IntidDB()

    for intid in args.intids:
        print()
        print('intid:', intid)

        print('Bits:')
        print('     3         2         1')
        print('    10987654321098765432109876543210')
        print('    {0:32b}'.format(intid))

        family, isAnti, *others = intIdModule.parseIntid(intid)
        print('Parts:')
        print('   family:', family)
        print('   isAnti:', isAnti)
        print('   others:', *others)
        print('id:', intIdModule.idFromIntid(intid, intidDB))
