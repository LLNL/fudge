#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import argparse, traceback, sys

def readEvaluation( filename, verbose=True, skipBadData=True, continuumSpectraFix=False ):
    '''
    Read in an evaluation in either Fudge/GNDS or ENDF and return result as Fudge classes
    '''

    firstline = open(filename).readline()

    # Is the file a GNDS file?
    if firstline.startswith( "<?xml" ):
        import fudge.gnds
        RS = fudge.gnds.reactionSuite.readXML( filename)
        try:
            CS = fudge.gnds.covariances.covarianceSuite.readXML( filename.replace('.gnds.', '.gndsCov.'))
        except:
            CS = None
        return {'reactionSuite':RS, 'covarianceSuite':CS, 'errors':[]}

    # Maybe its an ENDF file?
    elif firstline.endswith(' 0  0    0\n') or \
            firstline.endswith(' 0  0    0\r') or \
            firstline.endswith(' 0  0    0\r\n') or \
            filename.endswith('.endf'):
        from fudge.legacy.converting import endfFileToGNDS
        return endfFileToGNDS.endfFileToGNDS( filename, toStdOut=verbose, skipBadData=skipBadData, continuumSpectraFix=continuumSpectraFix )

    # Failed!
    else: print "WARNING: Unknown file type, not reading %s"% filename

def process_args():
    '''Process command line options'''
    parser = argparse.ArgumentParser(description='Check an ENDF or GNDS file')
    parser.add_argument('inFile', type=str, help='The ENDF or GNDS file you want to check.' )
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Enable verbose output' )
    parser.add_argument('--skipCov',default=False,action='store_true',help='skip covariance checks' )
    parser.add_argument('--skipEnergyBalance',default=False,action='store_true',help='skip energy balance checks' )
    parser.add_argument('--skipRes',default=False,action='store_true',help='skip resonance reconstruction and associated checks' )
    parser.add_argument('--traceback',default=False, action='store_true',help='on exception, print the python traceback')
    parser.add_argument( "--continuumSpectraFix", default = False, action = "store_true", help = "Skip unnormalizeable continuum gamma distributions" )
    return parser.parse_args()

if __name__=="__main__":
    args=process_args()

    if args.skipRes: raise NotImplementedError( "write skipRes option" )

    # Read the evaluation
    print "\nErrors encountered on read of "+args.inFile
    print "------------------------------------------------"
    result = readEvaluation( args.inFile, verbose=args.verbose, skipBadData=True, continuumSpectraFix=args.continuumSpectraFix )
    myEval=result['reactionSuite']
    myCov=result['covarianceSuite']
    print '\n'.join( result['errors'] )
    print '\n\n'


    # Check the evaluation
    try:
        print "Checking evaluation for "+args.inFile
        print "------------------------------------------------"
        print myEval.check(checkEnergyBalance=not args.skipEnergyBalance)
        print '\n'
    except Exception, err:
        print "Checking evaluation failed, got", str(err)
        print '\n'
        if args.traceback:
            traceback.print_exc(file=sys.stdout)


    # Check the covariance
    if not args.skipCov and myCov is not None:
        try:
            print "Checking covariances for "+args.inFile
            print "------------------------------------------------"
            print myCov.check()
            print '\n'
        except Exception, err:
            print "Checking covariance failed, got", str(err)
            print '\n'
            if args.traceback:
                traceback.print_exc(file=sys.stdout)
