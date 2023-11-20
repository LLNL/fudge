#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

summaryDocStringPoPs = '''Builds TNSL target ids and intids to an intid file.'''

description = '''
    Builds TNSL target ids and intids to an intid file.

    If the option "--inputPath" is present, that file is read in as an intid file and its contents are used to initialize the built intid file.
'''

import pathlib
import argparse

from PoPs import intId as intIdModule

parser = argparse.ArgumentParser(description=description)
parser.add_argument('outputPath', type=pathlib.Path,                    help='The path of the output file.')
parser.add_argument('--inputPath', type=pathlib.Path,                   help='The path of the input file.')
parser.add_argument('ids', nargs='*',                                   help='List of TNSL target ids.')

if __name__ == '__main__':
    args = parser.parse_args()

    intidDB = intIdModule.IntidDB(args.inputPath)

    intidsNext = 0
    for intid in intidDB.intids:
       intidsNext = max(intidsNext, intIdModule.getBits(intid, 0, 20) + 1)
       family, isAnti, *others = intIdModule.parseIntid(intid)
       if family == intIdModule.Family.TNSL:
            intidsNext = max(intidsNext, others[0] + 1)

    for id in args.ids:
        intidDB.add(id, intIdModule.setTNSL_intid(intidsNext))
        intidsNext += 1

    intidDB.saveToFile(args.outputPath)
