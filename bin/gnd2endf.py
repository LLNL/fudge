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

import argparse, os
from fudge.gnd import reactionSuite, covariances

# Process command line options
parser = argparse.ArgumentParser(description='Translate GND into ENDF')
parser.add_argument('inFilePrefix', type=str, help='The prefix of the GND files you want to translate.' )
parser.add_argument('-o', dest='outFile', default=None, help='Specify the output file' )
args = parser.parse_args()

# Compute input file names
inEvalFile = args.inFilePrefix + '.gnd.xml'
inCovFile = args.inFilePrefix + '.gndCov.xml'

# Compute the output file name
if args.outFile == None: outFile = args.inFilePrefix + '.endf'
else:                    outFile = args.outFile
    
# Read in XML files
myEval = reactionSuite.readXML( inEvalFile )
if os.path.exists( inCovFile ): myCov = covariances.readXML( inCovFile, reactionSuite=myEval )
else:                           myCov = None

# Now translate
with open( outFile, mode='w' ) as fout:
    fout.write( myEval.toENDF6( {'verbosity':0}, covarianceSuite=myCov ) )

