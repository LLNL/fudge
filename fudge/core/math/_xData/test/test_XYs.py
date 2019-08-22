if __name__ == "__main__": 

    import unittest
    from fudge.core.math._xData.XYs import XYs
    import fudge.core.math._xData.axes as axes
    
    class Test_XYs( unittest.TestCase ): 

        def setUp( self ): 
            vl1 = axes.axes( )
            vl1[0] = axes.axis( 'energy_in', 0, 'eV', interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
            vl1[1] = axes.axis( 'crossSection', 1, 'b' )
            self.a = XYs( vl1, data=[], accuracy=1e-3 )
    
        def test_init( self ):
            self.assertNotEqual( self.a, None )

        def test_getOtherInterpolationFunctionAndInfo( self ):
            '''This function is only defined in the compiled extension of XYs (pointwiseXYs?).  If it is not an attribute, we are not using the compiled version'''
            self.assertTrue( hasattr( self.a, "getOtherInterpolationFunctionAndInfo" ), "Missing attribute getOtherInterpolationFunctionAndInfo" )
            
    unittest.main()
