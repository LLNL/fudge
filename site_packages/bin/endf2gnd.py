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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

import sys, os
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert(0, os.path.dirname( binDir ) )

from fudge.legacy.converting import endfFileToGND, endfFormats
from fudge.legacy.converting.ENDFToGND import endfFileToGNDMisc
from fudge.processing import processingInfo

def process_args():
    # see http://docs.python.org/lib/optparse-tutorial.html
    import argparse
    parser = argparse.ArgumentParser(description="translates an ENDF file to the new GND format")
    parser.add_argument("inFile", help="input file")
    parser.add_argument("-v", action="store_true", dest="verbose", help="enable verbose output")
    parser.add_argument("-q", action="store_false", dest="verbose", help="disable verbose output")
    parser.add_argument("--skipBadData", action="store_true", default=False, help="skip bad data, rather than throw an exception, when reading an ENDF file")
    return parser.parse_args()


args = process_args()

inFile = args.inFile
outFile = inFile+'.gnd.xml'
outCovFile = inFile+'.gndCov.xml'

try:
    results = endfFileToGND.endfFileToGND( inFile, toStdOut = args.verbose, skipBadData = args.skipBadData )
    x = results['reactionSuite']
    c = results['covarianceSuite']
except Exception as err:
    sys.stderr.write( 'WARNING: ENDF READ HALTED BECAUSE '+str(err) )
    exit()
        
f = open( outFile, 'w' )
try:
    f.write( '\n'.join( x.toXMLList( ) + [ '' ] ) )
except Exception as err:
    sys.stderr.write( 'WARNING: MAIN ENDF WRITE HALTED BECAUSE '+str(err) )
    exit()
f.close( )
if c:
    f = open( outCovFile, 'w' )
    try:
        f.write( '\n'.join( c.toXMLList( ) + [ '' ] ) )
    except Exception as err:
        sys.stderr.write( 'WARNING: COVARIANCE ENDF WRITE HALTED BECAUSE '+str(err) )
        exit()
    f.close()

