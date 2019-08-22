#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

import os
import argparse
from fudge.legacy.converting.endfFileToGND import endfFileToGND

# Process command line options
parser = argparse.ArgumentParser(description='Translate ENDF into GND')
parser.add_argument('inFile', type=str, help='ENDF-6 file to be translated' )
parser.add_argument('-o', dest='outFilePrefix', default=None,
    help="""Specify the output file's prefix to be ``outFilePrefix``.
Resulting files have extensions ".gnd.xml" (main evaluation) and ".gndCov.xml" (covariance data)""" )
parser.add_argument('-r', dest='reconstruct', default=False, action='store_true',
    help='Reconstruct pointwise cross sections from resonance parameters after translating?' )
args = parser.parse_args()

# Compute output file names
inFile = os.path.split( args.inFile )[-1]
if args.outFilePrefix != None:
    outEvalFile = args.outFilePrefix + '.gnd.xml'
    outCovFile = args.outFilePrefix + '.gndCov.xml'
else:
    outEvalFile = inFile.replace( '.endf', '.gnd.xml' )
    outCovFile = inFile.replace( '.endf', '.gndCov.xml' )

# Now translate
translated = endfFileToGND( args.inFile, toStdOut=True, skipBadData=False )

if translated['errors']:
    print( "Warning! The following errors occured during translation:" )
    for warn in translated['errors']:
        print( warn )

if args.reconstruct:
    translated['reactionSuite'].reconstructResonances(0.001)

with open( outEvalFile, mode='w' ) as fout:
    fout.writelines( line+'\n' for line in translated['reactionSuite'].toXMLList( {'verbosity':0} ) )
if translated['covarianceSuite'] is not None:
    with open( outCovFile, mode='w' ) as fout:
        fout.writelines( line+'\n' for line in translated['covarianceSuite'].toXMLList( {'verbosity':0} ) )
