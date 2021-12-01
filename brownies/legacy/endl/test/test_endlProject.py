# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import unittest
import os

from brownies.legacy.endl import fudgeParameters
from brownies.legacy.endl import __path__
from brownies.legacy.endl.endlProject import endlProject

fudgeParameters.VerboseMode = False

class testEndlProject( unittest.TestCase ): 
    
    def testInit( self ):
        e = endlProject( os.sep.join( __path__ + [ 'test', 'testdb' ] ) )
        self.assertEqual( e.ZAList(), [ 'za001001' ] )

    def testReadZA( self ): 
        e = endlProject( os.sep.join( __path__ + [ 'test', 'testdb' ] ) )
        self.assertEqual( e.ZAList(), [ 'za001001' ] )
        za = e.readZA( 1001 )
        za.read()
        self.assertEqual( [ x.C for x in za.findDatas( I=0 ) ], [ 1, 10, 46 ] )

if __name__== "__main__":
    unittest.main()
