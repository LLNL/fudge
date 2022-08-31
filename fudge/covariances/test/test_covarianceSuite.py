#! /usr/bin/env python3
# encoding: utf-8

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
test fudge/covariances/
dbrown, 12/5/2012
"""

import unittest, os

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge.covariances import *
from fudge import reactionSuite
from fudge.covariances import covarianceSuite

TEST_DATA_PATH, this_filename = os.path.split(__file__)
FeEvaluation =  reactionSuite.ReactionSuite.readXML_file( TEST_DATA_PATH+os.sep+'n-026_Fe_056-endfbvii.1.endf.gnds.xml')
FeCovariance =  covarianceSuite.CovarianceSuite.readXML_file( TEST_DATA_PATH+os.sep+'n-026_Fe_056-endfbvii.1.endf.gndsCov.xml', reactionSuite=FeEvaluation )


class TestCaseBase( unittest.TestCase ):

    def assertXMLListsEqual(self,x1,x2):

        listAssertEqual = getattr( self, 'assertCountEqual', None )
        if( listAssertEqual is None ) : listAssertEqual = getattr( self, 'assertItemsEqual' )

        x1List = []
        for line in x1:
            if '\n' in line: x1List += [x.strip() for x in line.split('\n')]
            else: x1List.append(line.strip())
        x2List = []
        for line in x2:
            if '\n' in line: x2List += [x.strip() for x in line.split('\n')]
            else: x2List.append(line.strip())
        return listAssertEqual( x1List, x2List )

class Test_covarianceSuite( TestCaseBase ):

    def setUp(self):
        self.covSuite = FeCovariance

    def test_readXML(self):
        '''Also tests __init__ & parseNode'''
        self.assertIsInstance( self.covSuite, covarianceSuite.CovarianceSuite )

    def test__getitem__(self):
        answer = """<covarianceSection label="nonelastic">
    <rowData ENDF_MFMT="33,3" href="$reactions#/reactionSuite/sums/crossSectionSums/crossSectionSum[@label='nonelastic']/crossSection/resonancesWithBackground[@label='eval']"/>
    <mixed label="eval">
      <sum label="0" domainMin="1e-05" domainMax="862270.0" domainUnit="eV">
        <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->
        <summand ENDF_MFMT="33,102" coefficient="1.0" href="/covarianceSuite/covarianceSections/covarianceSection[@label='Fe57 + photon']/mixed[@label='eval']"/></sum>
      <covarianceMatrix label="1" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 862270 4e6 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="boundaries">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="3,3" compression="diagonal">
            <values>0 4e-4 9e-4</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="2" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 862270 1e6 2e6 4e6 6e6 8e6 1.4e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="boundaries">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="8,8" compression="diagonal">
            <values>0 1.584e-3 0.025344 1.584e-3 1.584e-3 1.584e-3 1.584e-3 1.584e-3</values></array></gridded2d></covarianceMatrix>
      <shortRangeSelfScalingVariance label="3" type="absolute" dependenceOnProcessedGroupWidth="inverse">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 862270 1e6 2e6 4e6 6e6 8e6 1.4e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="boundaries">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="8,8" compression="diagonal">
            <values>0 4.84e-10 1.357e-5 9.8057e-6 2.6915e-5 3.3173e-5 3.2769e-5 3.0954e-5</values></array></gridded2d></shortRangeSelfScalingVariance></mixed></covarianceSection>""".split('\n')
        self.maxDiff = None
        self.assertXMLListsEqual( self.covSuite.covarianceSections[2].toXML_strList(formatVersion=GNDS_formatVersionModule.gnds2_0), answer )

    def test__len__(self):
        self.assertEqual( len(self.covSuite.covarianceSections),60)

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

    def test_toXML_strList(self):
        '''Tested already in test__getitem__'''
        pass

    def test_toENDF6(self):
        pass

if __name__=="__main__":
    unittest.main()
