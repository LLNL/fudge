#!/usr/bin/env python
# encoding: utf-8

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

"""
test fudge/core/math
cmattoon, 3/24/2011
"""

import math
import unittest
from xData.series1d import LegendreSeries, Legendre

__metaclass__ = type

class testLegendre(unittest.TestCase):

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

