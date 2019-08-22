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