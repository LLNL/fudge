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

summaryDocStringIAM = '''Reads in an isotopic abundance file (i.e., an isotopicAbundancesByChemicalElement node) and prints the isotopic abundance sum for each of its chemical element.'''

description = summaryDocStringIAM

parser = argparse.ArgumentParser(description=description)
parser.add_argument('file', type=pathlib.Path,                  help='The isotopic abundance file.')

args = parser.parse_args()

isotopicAbundance = isotopicAbundancesModule.read(args.file)

for chemicalElement in isotopicAbundance.chemicalElements:
    sum = chemicalElement.atomFractionSum()
    print('    %-2s %-20s %s' % (chemicalElement.symbol, sum, 1.0 - sum))
