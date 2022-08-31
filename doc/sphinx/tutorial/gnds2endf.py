#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse, os
from fudge import reactionSuite
from fudge.covariances import covarianceSuite

# Process command line options
parser = argparse.ArgumentParser(description='Translate GNDS into ENDF')
parser.add_argument('inFilePrefix', type=str, help='The prefix of the GNDS files you want to translate.' )
parser.add_argument('-o', dest='outFile', default=None, help='Specify the output file' )
parser.add_argument('--style', type=str, help='Data style to translate back to ENDF-6', default='eval' )
args = parser.parse_args()

# Compute input file names
inEvalFile = args.inFilePrefix + '.gnds.xml'
inCovFile = args.inFilePrefix + '.gndsCov.xml'

# Compute the output file name
if args.outFile == None: outFile = args.inFilePrefix + '.endf'
else:                    outFile = args.outFile
    
# Read in XML files
myEval = reactionSuite.ReactionSuite.readXML_file( inEvalFile )
if os.path.exists( inCovFile ): myCov = covarianceSuite.CovarianceSuite.readXML_file( inCovFile, reactionSuite=myEval )
else:                           myCov = None

# Now translate
open( outFile, mode='w' ).write( myEval.toENDF6( args.style, {'verbosity':0}, covarianceSuite=myCov ) )

