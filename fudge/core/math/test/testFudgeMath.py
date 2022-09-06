#! /usr/bin/env python3
# encoding: utf-8

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
test fudge/core/math
cmattoon, 3/24/2011
"""

import math
import unittest
from xData.series1d import LegendreSeries, Legendre


class TestLegendre(unittest.TestCase):

    def test0(self): 
        """Simple test, just check a couple of values for n = 0,1,2"""
        self.assertEqual( Legendre( 0, 0.0 ), 1.0 )
        self.assertEqual( Legendre( 1, 1.0 ), 1.0 )
        self.assertEqual( Legendre( 2, 1.0 ), 1.0 )
        self.assertEqual( Legendre( 1, -1.0 ), -1.0 )
        
    def test1(self):
        """Check that all odd order polynomials evaluate to zero at the origin"""
        for n in range( 50 ):
            self.assertEqual( Legendre( 2*n+1, 0.0 ), 0.0 )
        
    def test2(self): 
        """Check that all odd order polynomials are in fact odd and even ones are in fact even"""
        for n in range( 50 ):
            self.assertEqual( Legendre( 2*n+1, 0.5 ), -Legendre( 2*n+1, -0.5 ) )
            self.assertEqual( Legendre( 2*n, 0.5 ), Legendre( 2*n, -0.5 ) )

    def test3(self): 
        """Check some more values, these are taken from Abramowitz and Stegun, 8.15 Example 1"""
        self.assertEqual( Legendre( 0, 0.3141592654 ), 1.0 )
        self.assertEqual( Legendre( 1, 0.3141592654 ), 0.3141592654 )
        self.assertAlmostEqual( Legendre( 2, 0.3141592654 ), -0.3519559340 )
        self.assertAlmostEqual( Legendre( 3, 0.3141592654 ), -0.3937232064 )
        self.assertAlmostEqual( Legendre( 4, 0.3141592654 ), 0.0475063122 )
        self.assertAlmostEqual( Legendre( 5, 0.3141592654 ), 0.3418427517 )
        self.assertAlmostEqual( Legendre( 6, 0.3141592654 ), 0.1572986975 )
        self.assertAlmostEqual( Legendre( 7, 0.3141592654 ), -0.2012339354 )
        self.assertAlmostEqual( Legendre( 8, 0.3141592654 ), -0.2561729328 )
        self.assertEqual( Legendre( 0, 2.6, False ), 1.0 )
        self.assertEqual( Legendre( 1, 2.6, False ), 2.6 )
        self.assertAlmostEqual( Legendre( 2, 2.6, False ), 9.64 )
        self.assertAlmostEqual( Legendre( 3, 2.6, False ), 40.04 )
        self.assertAlmostEqual( Legendre( 4, 2.6, False ), 174.952 )
        self.assertAlmostEqual( Legendre( 5, 2.6, False ), 786.74336 )
        self.assertAlmostEqual( Legendre( 6, 2.6, False ), 3604.350016 )
        self.assertAlmostEqual( Legendre( 7, 2.6, False ), 16729.51005, 5 )
        self.assertAlmostEqual( Legendre( 8, 2.6, False ), 78402.55522, 4 )

    def test4(self): 
        """Check range testing of parameters"""
        self.assertRaises( ValueError, Legendre, n=-1, mu=0.0 )
        self.assertRaises( ValueError, Legendre, n=1, mu=10.0 )
        
    def test5(self): 
        """Check recursion relation"""
        def lhs( n, x ): return (n+1)*Legendre(n+1,x)
        def rhs( n, x ): return (2*n+1)*x*Legendre(n,x) - n*Legendre(n-1,x)
        for n in range( 1, 50 ):
            self.assertAlmostEqual( lhs( n, 0.5 ), rhs( n, 0.5 ) )
        

if __name__ == '__main__':
    unittest.main()

