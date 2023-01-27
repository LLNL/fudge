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
Converts isotopic abundance data from a "special" "Commission on Isotopic Abundances and Atomic Weights" (CIAAW) file into an 
IsotopicAbundancesByChemicalElement instance and write to a file.
'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('CIAAW_file', type=pathlib.Path,            help='CIAAW isotopic abundance file.')
parser.add_argument('output', type=pathlib.Path,                help='Name of output file.')
parser.add_argument('evaluation', type=str,                     help='Evaluation name for the outputted file.')

args = parser.parse_args()

with args.CIAAW_file.open() as fIn:
    lines = fIn.readlines()

isotopicAbundancesByChemicalElement = isotopicAbundancesModule.IsotopicAbundancesByChemicalElement(args.evaluation)
isotopicAbundancesByChemicalElement.documentation.body.body = '''
Data from Commission on Isotopic Abundances and Atomic Weights (CIAAW) 2021. 
These data are from the web-site https://www.ciaaw.org/isotopic-abundances.htm.'''

for line in lines:
    Z = int(line.split()[0])
    if Z <= 92:                                 # Check if a new chemical element by looking for a chemcial element name matching Z also on the line.
        name = miscModule.nameFromZ[Z].lower()
        if name in line:
            symbol = line.split(name)[0].split()[1]
            print(Z, name, symbol)
            chemicalElement = isotopicAbundancesModule.ChemicalElement(symbol)
            isotopicAbundancesByChemicalElement.chemicalElements.add(chemicalElement)
            line = line.split(name)[-1]
    A = int(line.split()[0])
    line = ''.join(line.split()[1:])
    uncertainty = 0
    if '[' in line:
        rangeMin, rangeMax = line.split(',')
        rangeMin = float(rangeMin.split('[')[1])
        rangeMax = float(rangeMax.split(']')[0])
        atomFraction = 0.5 * (rangeMin + rangeMax)
        uncertainty = (rangeMax - rangeMin) / 2
    else:
        if '(' in line:
            atomFraction, uncertainty = line.split('(')
            uncertainty = uncertainty.split(')')[0]
            power10 = len(uncertainty)
            uncertainty = atomFraction + uncertainty
            atomFraction = float(atomFraction)
            uncertainty = 10**power10 * (float(uncertainty) - atomFraction)
        else:
            atomFraction = float(line)
    id = '%s%s' % (symbol, A)
    chemicalElement.isotopes.add(isotopicAbundancesModule.Isotope(id, atomFraction, uncertainty))

isotopicAbundancesByChemicalElement.saveToFile(args.output)
