#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import warning as warningModule
from fudge import GNDS_file as GNDS_fileModule

summaryDocStringFUDGE = '''Reads GNDS files and runs all of FUDGE physics tests on each.'''

description1 = '''Read one or more GNDS files into Fudge and run all physics tests.
    Sample use: python checkGNDS.py n-001_H_001.xml n-001_H_002.xml ...
    If file n-001_H_001-cov.xml (or -covar.xml) exists, covariances will automatically be read and checked.
'''

__doc__ = description1

parser = argparse.ArgumentParser(description1)
parser.add_argument('gnds', nargs='+',                                  help='GNDS and/or PoPs file(s) to check.')
parser.add_argument('-e', '--ebalance', action='store_true',            help='Include energy balance warnings.')
parser.add_argument('--normTolerance', type=float, default=1e-5,        help='Include energy balance warnings.')
parser.add_argument('-f', '--failOnException', action='store_true',     help='Be strict!')

if __name__ == '__main__':
    args = parser.parse_args()

    for fileName in args.gnds:

        gnds = GNDS_fileModule.read(fileName)

        warnings = gnds.check(checkEnergyBalance=args.ebalance, failOnException=args.failOnException, normTolerance=args.normTolerance)
        if args.ebalance:
            print(warnings)
        else:
            print(warnings.filter(exclude=[warningModule.EnergyImbalance]))

        covariances = []
        if hasattr(gnds, 'loadCovariances'):
            covariances = gnds.loadCovariances()

        for covarianceSuite in covariances:
            print('\nChecking covariance file %s' % covarianceSuite.sourcePath)
            covWarnings = covarianceSuite.check()
            print(covWarnings)
