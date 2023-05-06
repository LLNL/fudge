#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib

from PoPs.chemicalElements import misc as miscModule
from isotopicAbundances import isotopicAbundances as isotopicAbundancesModule

summaryDocStringIAM = '''This script prints the difference of two isotopic abundance files.'''

description = '''This script prints the difference of two isotopic abundance files 
(i.e., isotopicAbundancesByChemicalElement node files).'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('files', type=pathlib.Path, nargs=2,                       help='The two isotopic abundance files to diff.')

args = parser.parse_args()

isotopicAbundance1 = isotopicAbundancesModule.read(args.files[0])
chemicalElements1 = isotopicAbundance1.chemicalElements
isotopicAbundance2 = isotopicAbundancesModule.read(args.files[1])
chemicalElements2 = isotopicAbundance2.chemicalElements

symbols = set([chemicalElement.symbol for chemicalElement in chemicalElements1] +
              [chemicalElement.symbol for chemicalElement in chemicalElements2])
Zs_symbols = sorted([[miscModule.ZFromSymbol[symbol], symbol] for symbol in symbols])

for Z, symbol in Zs_symbols:
    print('    %s:' % symbol, end='')
    if symbol not in chemicalElements1:
        print(' chemical elements not in first isotopic abundance file.')
    elif symbol not in chemicalElements2:
        print(' chemical elements not in second isotopic abundance file.')
    else:
        print()
        isotopes1 = chemicalElements1[symbol].isotopes
        isotopes2 = chemicalElements2[symbol].isotopes
        isotopes = set([isotope.id for isotope in isotopes1] + [isotope.id for isotope in isotopes2])
        ZAs_id = sorted([[miscModule.ZAInfo_fromString(id)[2], id] for id in isotopes])
        for ZA, id in ZAs_id:
            print('        %s:' % id, end='')
            if id not in isotopes1:
                print(' isotope not in first isotopic abundance file.')
            elif id not in isotopes2:
                print(' isotope not in second isotopic abundance file.')
            else:
                atomFraction1 = isotopes1[id].atomFraction
                atomFraction2 = isotopes2[id].atomFraction
                diff = atomFraction1 - atomFraction2
                if diff != 0:
                    print(' difference %10.2e: %s versus %s' % (diff, atomFraction1, atomFraction2))
                else:
                    print()
