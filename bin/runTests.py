#! /usr/bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

import subprocess, argparse

testList = '''
fudge/core/math/test/testFudgeMath.py
fudge/core/math/test/testXYs.py
fudge/core/utilities/test/testUtilities.py
fudge/gnd/reactionData/test/test_crossSection.py
fudge/gnd/test/testCovariances.py
fudge/legacy/endl/test/test_endlProject.py
pqu/test/testPhysicalQuantityWithUncertainty.py
xmlrunner/tests/testsuite.py
xmlrunner/tests/testsuite_cases.py
'''.split('\n')[1:-1]

parser = argparse.ArgumentParser( "Unit test driver script" )
parser.add_argument( '-l', dest='listTests', default=False, action='store_true', help="List the unit test files" )
parser.add_argument( '-r', dest='runTests', default=False, action='store_true', help="Run all the unit tests" )
parser.add_argument( '-e', dest='editTests', default=False, action='store_true', help="Edit all the unit tests using the editor defined by the alias 'edit'" )

if __name__=="__main__":
    args = parser.parse_args()
    
    if args.listTests:
        print "Files to test:" 
        for t in testList:
            if t.startswith( '#' ): continue
            print '   ',t
    elif args.runTests:
        for t in testList:
            if t.startswith( '#' ): continue
            print '\n\n\n\n\nRunning', t, '...'
            subprocess.call( ['python', t ] )
    elif args.editTests:
        print ' '.join(  [ 'edit' ] + testList )
        #subprocess.call( [ 'edit' ] + testList )