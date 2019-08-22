#! /usr/bin/env python
import argparse
from fudge.legacy.converting.endfFileToGND import endfFileToGND

# Process command line options
parser = argparse.ArgumentParser(description='Check an ENDF file')
parser.add_argument('inFile', type=str, help='The ENDF file you want to check.' )
parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Enable verbose output' )
parser.add_argument('-o', dest='outFile', default=None, type=str, help='Output file prefix for gnd files (.gnd.xml and .gndCov.xml)' )

args = parser.parse_args()

# Now translate
result = endfFileToGND( args.inFile, toStdOut=args.verbose, skipBadData=True )
myEval=result['reactionSuite']
myCov=result['covarianceSuite']
print '\n\n'

# Check the evaluation
print "Errors encountered on read of "+args.inFile
print "------------------------------------------------"
print '\n'.join( result['errors'] )
print '\n'

# Check the evaluation
try:
    print "Checking evaluation for "+args.inFile
    print "------------------------------------------------"
    print myEval.check()
    print '\n'
except Exception, err: 
    print "Checking evaluation failed, got", str(err)
    print '\n'


# Check the covariance
try:
    print "Checking covariances for "+args.inFile
    print "------------------------------------------------"
    print myCov.check()
    print '\n'
except Exception, err: 
    print "Checking covariance failed, got", str(err)
    print '\n'

if args.outFile != None:
    myEval.saveToFile( args.outFile+'.gnd.xml' )
    myCov.saveToFile( args.outFile+'.gndCov.xml' )