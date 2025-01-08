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

from PoPs import database as databaseModule
from PoPs.chemicalElements import misc as miscModule
from fudge.core.utilities import abundance as abundanceModule

from isotopicAbundances import isotopicAbundances as isotopicAbundancesModule

description ='''
Converts a legacy LLNL FUDGE isotopic abundance file into an IsotopicAbundancesByChemicalElement instance and write to a file.
Also, prints some stats about the data.
'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('pops', type=pathlib.Path,                  help='A GNDS pops file.')
parser.add_argument('output', type=pathlib.Path,                help='Name of output file.')
parser.add_argument('evaluation', type=str,                     help='Evaluation name for the outputted file.')

args = parser.parse_args()

pops = databaseModule.read(args.pops)

isotopicAbundancesByChemicalElement = isotopicAbundancesModule.IsotopicAbundancesByChemicalElement(args.evaluation)

missing = []
atomFractionSums = []
for chemicalElement in pops.chemicalElements:
    if chemicalElement.Z > 92:
        continue
    try:
        isotopes = abundanceModule.getElementsNaturalIsotopes(chemicalElement.Z)
    except:
        missing.append((chemicalElement.symbol, chemicalElement.Z))
        continue

    chemicalElement2 = isotopicAbundancesModule.ChemicalElement(chemicalElement.symbol)
    for isotope in isotopes:
        atomFraction, uncertainty = isotopes[isotope]
        id = miscModule.idFromZAndA(chemicalElement.Z, isotope, True)
        chemicalElement2.isotopes.add(isotopicAbundancesModule.Isotope(id, atomFraction / 100., uncertainty / 100.))
    if len(chemicalElement2.isotopes) > 0:
        isotopicAbundancesByChemicalElement.chemicalElements.add(chemicalElement2)
        atomFractionSums.append((chemicalElement.symbol, chemicalElement.Z, chemicalElement2.atomFractionSum()))

isotopicAbundancesByChemicalElement.saveToFile(args.output)

if len(missing) > 0:
    print('Not data for elements:', file=sys.stderr)
    for symbol, Z in missing:
        print('    %-4s %s' % (Z, symbol), file=sys.stderr)

print('\nSums:', file=sys.stderr)
for Z, symbol, atomFractionSum in atomFractionSums:
    if abs(atomFractionSum - 1.0) > 1e-14:
        print('    %-4s %-3s %.10s' % (Z, symbol, atomFractionSum), file=sys.stderr)
