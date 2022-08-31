#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys, os, glob
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert(0, os.path.dirname( binDir ) )

from brownies.legacy.converting import endfFileToGNDS

def process_args():
    # see http://docs.python.org/lib/optparse-tutorial.html
    import argparse
    parser = argparse.ArgumentParser(description="translates an ENDF files to the new GNDS format")
    parser.add_argument("-v", action="store_true", dest="verbose", help="enable verbose output")
    parser.add_argument("-q", action="store_false", dest="verbose", help="disable verbose output")
    parser.add_argument("-k", action="store_true", dest="keepGoing", default=False, help="keep going, namely skip to next file on error")
    parser.add_argument("--noSkipBadData", dest='skipBadData', action="store_false", default=True, help="do not skip bad data when reading an ENDF file, throw an exception instead")
    return parser.parse_args()


args = process_args()


for inFile in glob.glob('*.endf'):
    outFile = inFile.replace('.endf','.gnds.xml')
    outCovFile = inFile.replace('.endf','.gndsCov.xml')

    x=None
    c=None

    if args.verbose:
        print()
        print()
        print(50*'-')
        print("Translating ENDF file into GNDS: ", inFile)
        print(50 * '-')

    try:
        results = endfFileToGNDS.endfFileToGNDS( inFile, toStdOut = args.verbose, skipBadData = args.skipBadData, verbose = 0 )
        x = results['reactionSuite']
        c = results['covarianceSuite']
    except Exception as err:
        sys.stderr.write( 'WARNING: ENDF READ HALTED BECAUSE "'+str(err)+'" for file %s\n'%inFile )
        if not args.keepGoing: exit()

    f = open( outFile, 'w' )
    try:
        f.write( '\n'.join( x.toXML_strList(  ) ) )
    except Exception as err:
        sys.stderr.write( 'WARNING: MAIN ENDF WRITE HALTED BECAUSE "'+str(err)+'" for file %s\n'%inFile )
        if not args.keepGoing: exit()
    f.close( )
    if c:
        f = open( outCovFile, 'w' )
        try:
            f.write( '\n'.join( c.toXML_strList(  ) ) )
        except Exception as err:
            sys.stderr.write( 'WARNING: COVARIANCE ENDF WRITE HALTED BECAUSE  "'+str(err)+'" for file %s\n'%inFile )
            if not args.keepGoing: exit()
        f.close()
