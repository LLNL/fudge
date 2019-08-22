#!/usr/bin/env python
# encoding: utf-8

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

"""
test fudge/core/math
cmattoon, 3/24/2011
"""

import math
import unittest
from fudge.core.math.xData import LegendreSeries

__metaclass__ = type

class testLegendre(unittest.TestCase):

    def test0(self): 
        '''Simple test, just check a couple of values for n = 0,1,2'''
        self.assertEqual( LegendreSeries.Legendre( 0, 0.0 ), 1.0 )
        self.assertEqual( LegendreSeries.Legendre( 1, 1.0 ), 1.0 )
        self.assertEqual( LegendreSeries.Legendre( 2, 1.0 ), 1.0 )
        self.assertEqual( LegendreSeries.Legendre( 1, -1.0 ), -1.0 )
        
    def test1(self):
        '''Check that all odd order polynomials evaluate to zero at the origin'''
        for n in range( 50 ):
            self.assertEqual( LegendreSeries.Legendre( 2*n+1, 0.0 ), 0.0 )           
        
    def test2(self): 
        '''Check that all odd order polynomials are in fact odd and even ones are in fact even'''
        for n in range( 50 ):
            self.assertEqual( LegendreSeries.Legendre( 2*n+1, 0.5 ), -LegendreSeries.Legendre( 2*n+1, -0.5 ) )           
            self.assertEqual( LegendreSeries.Legendre( 2*n, 0.5 ), LegendreSeries.Legendre( 2*n, -0.5 ) )           

    def test3(self): 
        '''Check some more values, these are taken from Abramowitz and Stegun, 8.15 Example 1'''
        self.assertEqual( LegendreSeries.Legendre( 0, 0.3141592654 ), 1.0 )
        self.assertEqual( LegendreSeries.Legendre( 1, 0.3141592654 ), 0.3141592654 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 2, 0.3141592654 ), -0.3519559340 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 3, 0.3141592654 ), -0.3937232064 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 4, 0.3141592654 ), 0.0475063122 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 5, 0.3141592654 ), 0.3418427517 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 6, 0.3141592654 ), 0.1572986975 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 7, 0.3141592654 ), -0.2012339354 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 8, 0.3141592654 ), -0.2561729328 )
        self.assertEqual( LegendreSeries.Legendre( 0, 2.6, False ), 1.0 )
        self.assertEqual( LegendreSeries.Legendre( 1, 2.6, False ), 2.6 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 2, 2.6, False ), 9.64 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 3, 2.6, False ), 40.04 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 4, 2.6, False ), 174.952 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 5, 2.6, False ), 786.74336 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 6, 2.6, False ), 3604.350016 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 7, 2.6, False ), 16729.51005, 5 )
        self.assertAlmostEqual( LegendreSeries.Legendre( 8, 2.6, False ), 78402.55522, 4 )

    def test4(self): 
        '''Check range testing of parameters'''
        self.assertRaises( ValueError, LegendreSeries.Legendre, n=-1, x=0.0 )
        self.assertRaises( ValueError, LegendreSeries.Legendre, n=1, x=10.0 )
        
    def test5(self): 
        '''Check recursion relation'''
        def lhs( n, x ): return (n+1)*LegendreSeries.Legendre(n+1,x)
        def rhs( n, x ): return (2*n+1)*x*LegendreSeries.Legendre(n,x) - n*LegendreSeries.Legendre(n-1,x)
        for n in range( 1, 50 ):
            self.assertAlmostEqual( lhs( n, 0.5 ), rhs( n, 0.5 ) )
        

if __name__ == '__main__':
    unittest.main()

