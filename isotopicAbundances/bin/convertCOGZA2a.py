#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import argparse
import pathlib

from PoPs.chemicalElements import misc as miscModule
from isotopicAbundances import isotopicAbundances as isotopicAbundancesModule

description ='''
Converts an LLNL COG isotopic abundance file into an IsotopicAbundancesByChemicalElement instance and write to a file.
Also, prints some stats about the data.
'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('COG_file', type=pathlib.Path,              help='COG isotopic abundance file.')
parser.add_argument('output', type=pathlib.Path,                help='Name of output file.')
parser.add_argument('evaluation', type=str,                     help='Evaluation name for the outputted file.')

args = parser.parse_args()

isotopicAbundancesByChemicalElement = isotopicAbundancesModule.IsotopicAbundancesByChemicalElement(args.evaluation)

missing = []
atomFractionSums = []
with args.COG_file.open() as fIn:
    for line in fIn.readlines():
        values = line.split()

        Z = int(values[0]) // 1000
        symbol = miscModule.symbolFromZ[Z]

        count = int(values[1])
        if count == 0:
            if Z <= 92:
                missing.append((symbol, Z))
            continue

        chemicalElement = isotopicAbundancesModule.ChemicalElement(symbol)
        for index in range(2, 2 * (count + 1),2):
            ZA = int(values[index])
            atomFraction = float(values[index+1])
            id = miscModule.idFromZA(ZA)
            chemicalElement.isotopes.add(isotopicAbundancesModule.Isotope(id, atomFraction, -1.0))
        isotopicAbundancesByChemicalElement.chemicalElements.add(chemicalElement)
        atomFractionSums.append((symbol, Z, chemicalElement.atomFractionSum()))

isotopicAbundancesByChemicalElement.saveToFile(args.output)

if len(missing) > 0:
    print('Not data for elements:', file=sys.stderr)
    for symbol, Z in missing:
        print('    %-4s %s' % (Z, symbol), file=sys.stderr)

print('\nSums:', file=sys.stderr)
for Z, symbol, atomFractionSum in atomFractionSums:
    if abs(atomFractionSum - 1.0) > 1e-14:
        print('    %-4s %-3s %.10s' % (Z, symbol, atomFractionSum), file=sys.stderr)
