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
test fudge/gnd/covariances/
dbrown, 12/5/2012
"""

import unittest, cStringIO, os
from fudge.gnd.covariances import *
from fudge.gnd.covariances.covarianceSuite import readXML as CovReadXML
from fudge.gnd.reactionSuite import readXML as RxnReadXML

TEST_DATA_PATH, this_filename = os.path.split(__file__)
FeEvaluation =  RxnReadXML( open(TEST_DATA_PATH+os.sep+'n-026_Fe_056-endfbvii.1.endf.gnd.xml') )
FeCovariance =  CovReadXML( open(TEST_DATA_PATH+os.sep+'n-026_Fe_056-endfbvii.1.endf.gndCov.xml'), reactionSuite=FeEvaluation )


class TestCaseBase( unittest.TestCase ):
    def assertXMLListsEqual(self,x1,x2):
        x1List = []
        for line in x1:
            if '\n' in line: x1List += [x.strip() for x in line.split('\n')]
            else: x1List.append(line.strip())
        x2List = []
        for line in x2:
            if '\n' in line: x2List += [x.strip() for x in line.split('\n')]
            else: x2List.append(line.strip())
        return self.assertItemsEqual(x1List,x2List)


class Test_covarianceSuite( TestCaseBase ):
    
    def setUp(self): 
        self.covSuite = FeCovariance
    
    def test_readXML(self):
        '''Also tests __init__ & parseXMLNode'''
        self.assertIsInstance( self.covSuite, covarianceSuite.covarianceSuite )
        
    def test__getitem__(self):
        answer = """<section label="nonelastic">
    <rowData ENDF_MFMT="33,3" xlink:href="/reactionSuite/sums/crossSections/crossSectionSum[@label='nonelastic']/crossSection/resonancesWithBackground[@label='eval']/regions1d"/>
    <mixed label="eval">
      <sum label="0" domainMin="1e-05" domainMax="862270.0" domainUnit="eV">
        <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->
        <summand ENDF_MFMT="33,102" coefficient="1.0" xlink:href="/covarianceSuite/section[@label='Fe57 + photon']"/></sum>
      <covarianceMatrix label="1" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="4">1e-5 862270 4e6 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="4,4" compression="diagonal">
            <values length="4">0 4e-4 9e-4 0</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="2" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="9">1e-5 862270 1e6 2e6 4e6 6e6 8e6 1.4e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="9,9" compression="diagonal">
            <values length="9">0 1.584e-3 0.025344 1.584e-3 1.584e-3 1.584e-3 1.584e-3 1.584e-3 0</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="3" type="relative" ENDFconversionFlag="LB=8">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="9">1e-5 862270 1e6 2e6 4e6 6e6 8e6 1.4e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="9,9" compression="diagonal">
            <values length="9">0 4.84e-10 1.357e-5 9.8057e-6 2.6915e-5 3.3173e-5 3.2769e-5 3.0954e-5 0</values></array></gridded2d></covarianceMatrix></mixed></section>""".split('\n')
        self.maxDiff = None
        self.assertXMLListsEqual( self.covSuite[2].toXMLList(), answer )

    def test__len__(self):
        self.assertEqual( len(self.covSuite),40)

    def test_addSection(self):
        pass

    def test_addModelParameterCovariance(self):
        pass
    
    def test_addReactionSum(self):
        pass

    def test_addExternalReaction(self):
        pass

    def test_saveToOpenedFile( self):
        pass

    def test_saveToFile( self):
        pass
    
    def test_check( self ):
        pass

    def test_fix( self ):
        pass
    
    def test_removeExtraZeros(self):
        pass
     
    def test_toXMLList(self):
        '''Tested already in test__getitem__'''
        pass

    def test_toENDF6(self):
        pass

if __name__=="__main__":
    unittest.main()
