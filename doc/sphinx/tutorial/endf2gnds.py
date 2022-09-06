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
parser = argparse.ArgumentParser(description='Translate ENDF into GNDS')
parser.add_argument('inFile', type=str, help='The ENDF file you want to translate.' )
parser.add_argument('-o', dest='outFilePrefix', default=None, help='''Specify the output file's prefix to be ``outFilePrefix``.  The outputted files have extensions ".gnds.xml" and ".gndsCov.xml" for the GNDS main evaluations and covariance files.''' )
args = parser.parse_args()

# Compute output file names
if args.outFilePrefix != None:
    outEvalFile = args.outFilePrefix + '.gnds.xml'
    outCovFile = args.outFilePrefix + '.gndsCov.xml'
else:
    outEvalFile = args.inFile.replace( '.endf', '.gnds.xml' )
    outCovFile = args.inFile.replace( '.endf', '.gndsCov.xml' )
    
# Now translate
rce = endfFileToGNDS( args.inFile, toStdOut=True, skipBadData=True )
myEval, myCov = rce['reactionSuite'], rce['covarianceSuite']
open( outEvalFile, mode='w' ).writelines( line+'\n' for line in myEval.toXML_strList( ) )
if myCov is not None:
     open( outCovFile, mode='w' ).writelines( line+'\n' for line in myCov.toXML_strList( ) )
