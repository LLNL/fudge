#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import sys
import traceback


def readEvaluation(filename, verbose=True, skipBadData=True, continuumSpectraFix=False):
    """
    Read in an evaluation in either Fudge/GNDS or ENDF and return result as Fudge classes
    """

    firstline = open(filename).readline()

    # Is the file a GNDS file?
    if firstline.startswith("<reactionSuite ") or firstline.startswith("<?xml"):
        import fudge
        RS = fudge.reactionSuite.ReactionSuite.readXML_file(filename)
        try:
            CS = fudge.covariances.covarianceSuite.CovarianceSuite.readXML_file(filename.replace('.gnds.', '.gndsCov.'))
        except:
            CS = None
        return {'reactionSuite': RS, 'covarianceSuite': CS, 'errors': []}

    # Maybe its an ENDF file?
    elif firstline.endswith(' 0  0    0\n') or \
            firstline.endswith(' 0  0    0\r') or \
            firstline.endswith(' 0  0    0\r\n') or \
            filename.endswith('.endf'):
        from brownies.legacy.converting import endfFileToGNDS
        return endfFileToGNDS.endfFileToGNDS(filename, toStdOut=verbose, skipBadData=skipBadData,
                                             continuumSpectraFix=continuumSpectraFix)

    # Failed!
    else:
        print("WARNING: Unknown file type, not reading %s" % filename)


def process_args():
    """Process command line options"""
    parser = argparse.ArgumentParser(description='Check an ENDF or GNDS file')
    parser.add_argument('inFile', type=str, help='The ENDF or GNDS file you want to check.')
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Enable verbose output')
    parser.add_argument('--skipCov', default=False, action='store_true', help='skip covariance checks')
    parser.add_argument('--skipEnergyBalance', default=False, action='store_true', help='skip energy balance checks')
    parser.add_argument('--skipRes', default=False, action='store_true',
                        help='skip resonance reconstruction and associated checks')
    parser.add_argument('--traceback', default=False, action='store_true',
                        help='on exception, print the python traceback')
    parser.add_argument("--continuumSpectraFix", default=False, action="store_true",
                        help="Skip unnormalizeable continuum gamma distributions")
    return parser.parse_args()


if __name__ == "__main__":
    args = process_args()

    if args.skipRes:
        raise NotImplementedError("write skipRes option")

    # Read the evaluation
    print("\nErrors encountered on read of " + args.inFile)
    print("------------------------------------------------")
    result = readEvaluation(args.inFile, verbose=args.verbose, skipBadData=True,
                            continuumSpectraFix=args.continuumSpectraFix)
    myEval = result['reactionSuite']
    myCov = result['covarianceSuite']
    print('\n'.join(result['errors']))
    print('\n\n')

    # Check the evaluation
    try:
        print("Checking evaluation for " + args.inFile)
        print("------------------------------------------------")
        print(myEval.check(checkEnergyBalance=not args.skipEnergyBalance))
        print('\n')
    except Exception as err:
        print("Checking evaluation failed, got", str(err))
        print('\n')
        if args.traceback:
            traceback.print_exc(file=sys.stdout)

    # Check the covariance
    if not args.skipCov and myCov is not None:
        try:
            print("Checking covariances for " + args.inFile)
            print("------------------------------------------------")
            print(myCov.check())
            print('\n')
        except Exception as err:
            print("Checking covariance failed, got", str(err))
            print('\n')
            if args.traceback:
                traceback.print_exc(file=sys.stdout)
