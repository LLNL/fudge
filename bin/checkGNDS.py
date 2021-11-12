#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import argparse

from LUPY import GNDSType as GNDSTypeModule
from fudge import warning
from PoPs import database as databaseModule

description1 = """Read one or more GNDS files into Fudge and run all physics tests.
    Sample use: python checkGNDS.py n-001_H_001.xml n-001_H_002.xml ...
    If file n-001_H_001-cov.xml (or -covar.xml) exists, covariances will automatically be read and checked.
"""

__doc__ = description1

parser = argparse.ArgumentParser( description1 )
parser.add_argument( 'gnds', nargs = "+", help = "GNDS and/or PoPs file(s) to check" )
parser.add_argument( '-e', '--ebalance', action = 'store_true', help = "include energy balance warnings" )
parser.add_argument( '-f', '--failOnException', action = 'store_true', help = "Be strict!" )

if __name__ == '__main__' :
    args = parser.parse_args()

    for fileName in args.gnds :

        covariances = []
        name, dummy = GNDSTypeModule.type( fileName )
        if( name == databaseModule.database.moniker ) :
            gnds = GNDSTypeModule.read( fileName )
        else :

            gnds = GNDSTypeModule.read( fileName )
            if hasattr(gnds, 'loadCovariances'):
                covariances = gnds.loadCovariances()

        warnings = gnds.check( checkEnergyBalance = args.ebalance, failOnException = args.failOnException )
        if args.ebalance:
            print( warnings )
        else:
            print( warnings.filter( exclude=[warning.energyImbalance] ) )

        for covarianceSuite in covariances:
            print( "\nChecking covariance file %s" % covarianceSuite.sourcePath)
            covWarnings = covarianceSuite.check()
            print( covWarnings )

