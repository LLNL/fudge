#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
from brownies.legacy.converting.endfFileToGNDS import endfFileToGNDS

# Process command line options
parser = argparse.ArgumentParser(description='Check and ENDF file')
parser.add_argument('inFile', type=str, help='The ENDF file you want to translate.' )
args = parser.parse_args()

# Now translate
rce = endfFileToGNDS( args.inFile, toStdOut=True, skipBadData=True )
myEval, myCov = rce['reactionSuite'], rce['covarianceSuite']
print( '\n\n' )

# Check the evaluation
print( "Checking evaluation for "+args.inFile )
print( "------------------------------------------------" )
print( myEval.check() )

print( '\n' )

# Check the covariance
print( "Checking covariances for "+args.inFile )
print( "------------------------------------------------" )
print( myCov.check() )
