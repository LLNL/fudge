# <<BEGIN-copyright>>
# <<END-copyright>>

import unittest
import os
from fudge.legacy.endl import endlProject, __path__

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
