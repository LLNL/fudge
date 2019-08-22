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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
test fudge/gnd/covariances/
dbrown, 12/5/2012
"""

import unittest, copy
from pqu import PQU
from fudge.gnd.covariances import base
defaultAccuracy = 0.001

from fudge.core.utilities.xmlNode import xmlNode
from xml.etree import cElementTree

def covFromString( s ):
    return base.covarianceMatrix.parseXMLNode( xmlNode( cElementTree.fromstring( s ), xmlNode.etree ) )
         
class Test_covariance_baseClass( unittest.TestCase ): 

    def setUp(self):
        self.cov_const = covFromString( '''    <covarianceMatrix type="relative">
      <axes>
        <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="2"> 1e-5 2e7</axis>
        <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>
        <axis index="2" label="matrix_elements" unit=""/></axes>
      <matrix rows="1" columns="1" form="symmetric" precision="6">
         0.0144</matrix></covarianceMatrix>''' )

    def test_toXMLList(self):
        self.assertEqual( self.cov_const.toXMLList(), [
            '<covarianceMatrix type="relative">',
            '  <axes>',
            '    <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="2"> 1e-5 2e7</axis>',
            '    <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>',
            '    <axis index="2" label="matrix_elements" unit=""/></axes>',
            '  <matrix rows="1" columns="1" form="symmetric" precision="6">',
            '     1.440000e-02</matrix></covarianceMatrix>'] )
            
    def test_convertAxesToUnits(self): pass
    
    def test_toCovarianceMatrix(self): pass
    
    def test_toCorrelationMatrix(self): pass
    
    def test_getRowBounds(self): 
        self.assertEquals( self.cov_const.getRowBounds(), (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))

    def test_getColumnBounds(self):
        self.assertEquals( self.cov_const.getColumnBounds(), (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
    
    def test_toAbsolute_and_toRelative(self): 
        '''
        The constant covariance in self.cov_const is 0.01 everywhere and relative.  So, in absolute covariance
        it is :: 
        
            cov(E,E') = sigma(E)*0.0144*sigma(E')
        
        With a constant cross section of 1.5 b, the covariance should be (1.5 b)^2*0.0144 = 3.24e-2 b^2
        '''
        import fudge.gnd.reactionData.crossSection 
        pointwise = fudge.gnd.reactionData.crossSection.pointwise
        ptwise = pointwise(pointwise.defaultAxes(), [ [1e-5,1.5], [20.0e6,1.5] ], 1e-4 )
        original = copy.deepcopy(self.cov_const) 

        # start as relative, so first call should do nothing
        self.assertEqual( self.cov_const.toRelative().toXMLList(), original.toXMLList() )

        # call to absolute better change things
        self.assertEqual(
            self.cov_const.toAbsolute(ptwise).toXMLList(),
            [ x.replace('relative','absolute').replace('1.44','3.24') for x in original.toXMLList() ] )
    
        # call to relative better bring us back
        self.assertEqual( self.cov_const.toRelative(ptwise).toXMLList(), original.toXMLList() )
    
    def test_check(self): 
        self.assertEqual( self.cov_const.check({
            'checkUncLimits':False, 
            'negativeEigenTolerance':1e-8, 
            'eigenvalueRatioTolerance':1e8}), [] )
    
    def test_fix(self): pass
    
    def test_group(self): pass
    
    def test_removeExtraZeros(self): pass
    
    def test_getUncertaintyVector(self): 
        # Test it as-is
        self.assertEqual(
            [xy for xy in self.cov_const.getUncertaintyVector()],
            [[1e-05, 0.12], [20000000.0, 0.12]] )

        # Try making it absolute
        import fudge.gnd.reactionData.crossSection 
        pointwise = fudge.gnd.reactionData.crossSection.pointwise
        ptwise = pointwise( pointwise.defaultAxes(), [ [1e-5,1.5], [20.0e6,1.5] ], 1e-4 )
        self.assertEqual(
            [xy for xy in self.cov_const.getUncertaintyVector(ptwise, relative=False)],
            [[1e-05, 0.18], [20000000.0, 0.18]] )
    
    def test_toENDF6(self): pass

if __name__=="__main__":
    unittest.main()
