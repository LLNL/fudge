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

Note, the tests below operate on the one and only copy of the covariance data loaded below.  Make sure
each test only fools with one section otherwise you'll be chasing ghost bugs:

    * run 1: test A changes file, test B depends on results from test A

    * run 2 (different order): test B fails because test A didn't change the files in an "expected" way

"""

from fudge import reactionSuite
from fudge.covariances import covarianceSuite
from fudge.covariances import enums as covarianceEnumsModule

import unittest, os, copy

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge import covariances

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import link as linkModule
from xData import gridded as griddedModule
from xData import values as valuesModule
from xData import xDataArray as arrayModule

TEST_DATA_PATH, this_filename = os.path.split(__file__)
FeEvaluation = reactionSuite.ReactionSuite.readXML_file(
    TEST_DATA_PATH + os.sep + 'n-026_Fe_056-endfbvii.1.endf.gnds.xml')
FeCovariance = covarianceSuite.CovarianceSuite.readXML_file(
    TEST_DATA_PATH + os.sep + 'n-026_Fe_056-endfbvii.1.endf.gndsCov.xml', reactionSuite=FeEvaluation)


class TestCaseBase(unittest.TestCase):

    def assertXMLListsEqual(self, x1, x2):

        listAssertEqual = getattr(self, 'assertCountEqual', None)
        if listAssertEqual is None: listAssertEqual = getattr(self, 'assertItemsEqual')

        x1List = []
        for line in x1:
            if '\n' in line:
                x1List += [x.strip() for x in line.split('\n')]
            else:
                x1List.append(line.strip())
        x2List = []
        for line in x2:
            if '\n' in line:
                x2List += [x.strip() for x in line.split('\n')]
            else:
                x2List.append(line.strip())
        return listAssertEqual(x1List, x2List)


class Test_mixed(TestCaseBase):

    def setUp(self):
        # The COMMARA-2.0 33 group structure
        self.__groupBoundaries = [
            1.9640E+07, 1.0000E+07, 6.0653E+06, 3.6788E+06, 2.2313E+06, 1.3534E+06,
            8.2085E+05, 4.9787E+05, 3.0197E+05, 1.8316E+05, 1.1109E+05, 6.7380E+04,
            4.0868E+04, 2.4788E+04, 1.5034E+04, 9.1188E+03, 5.5308E+03, 3.3546E+03,
            2.0347E+03, 1.2341E+03, 7.4852E+02, 4.5400E+02, 3.0433E+02, 1.4863E+02,
            9.1661E+01, 6.7904E+01, 4.0169E+01, 2.2603E+01, 1.3710E+01, 8.3153E+00,
            4.0000E+00, 5.4000E-01, 1.0000E-01, 1e-5]
        self.__groupBoundaries.reverse()
        self.__groupUnit = 'eV'

        # ...................... example matrix 'a' ......................
        axes = axesModule.Axes(3, labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'),
                                               2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit,
                                  style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values([1.0000E-07, 1.1109E-01, 1.3534E+00, 1.9640E+01]))
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit, style=xDataEnumsModule.GridStyle.none,
                                  values=linkModule.Link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.Full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.Symmetry.lower)
        self.a = covariances.covarianceMatrix.CovarianceMatrix('eval',
                                                               matrix=griddedModule.Gridded2d(axes, myMatrix),
                                                               type=covarianceEnumsModule.Type.relative)

        # ...................... example matrix 'b' ......................
        axes = axesModule.Axes(3, labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'),
                                               2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit,
                                  style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values([1.0e-5, 0.100, 1.0, 20.0]))
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit, style=xDataEnumsModule.GridStyle.none,
                                  values=linkModule.Link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.Full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.Symmetry.lower)
        self.b = covariances.covarianceMatrix.CovarianceMatrix('eval', matrix=griddedModule.Gridded2d(axes, myMatrix),
                                                               type=covarianceEnumsModule.Type.relative)

        # ...................... example matrix 'c' ......................
        axes = axesModule.Axes(3, labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'),
                                               2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit,
                                  style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values([1.0000E-07, 6.7380E-02, 1.1109E-01, 1.3534E+00]))
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit, style=xDataEnumsModule.GridStyle.none,
                                  values=linkModule.Link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.Full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.Symmetry.lower)
        self.c = covariances.covarianceMatrix.CovarianceMatrix('eval', matrix=griddedModule.Gridded2d(axes, myMatrix),
                                                               type=covarianceEnumsModule.Type.relative)

        self.abc = covariances.mixed.MixedForm(components=[self.a, self.b, self.c])

    def test__getitem__(self):
        self.maxDiff = None
        self.assertXMLListsEqual(
            FeCovariance.covarianceSections[36]['eval'].toXML_strList(
                formatVersion=GNDS_formatVersionModule.version_1_10), '''<mixed label="eval">
      <covarianceMatrix label="0" type="absolute">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="2,2" compression="diagonal">
            <values>0 9e-8</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="1" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="2,2" compression="diagonal">
            <values>0 8e-2</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="2" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="5,5" compression="diagonal">
            <values>0 0.072 0.072 0.072 0.072</values></array></gridded2d></covarianceMatrix>
      <shortRangeSelfScalingVariance label="3" type="absolute" dependenceOnProcessedGroupWidth="inverse">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="5,5" compression="diagonal">
            <values>0 7.2e-14 1.5335e-14 2.9892e-12 1.9177e-9</values></array></gridded2d></shortRangeSelfScalingVariance></mixed>'''.split(
                '\n'))

    def test__len__(self):
        self.assertEqual(len(FeCovariance.covarianceSections), 60)

    def test_toXML_strList_11(self):
        """This covariance in the Fe56 file is of mixed form and is made of 4 diagonal subspaces."""
        self.maxDiff = None
        self.assertXMLListsEqual(
            FeCovariance.covarianceSections[11].toXML_strList(formatVersion=GNDS_formatVersionModule.version_1_10), '''  <section label="n + (Fe56_e5 -> Fe56 + photon)">
    <rowData ENDF_MFMT="33,55" href="$reactions#/reactionSuite/reactions/reaction[@label='n + (Fe56_e5 -> Fe56 + photon)']/crossSection/XYs1d[@label='eval']"/>
    <mixed label="eval">
      <covarianceMatrix label="0" type="absolute">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 3013400 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="2,2" compression="diagonal">
            <values>0 2e-6</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="1" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 3013400 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="2,2" compression="diagonal">
            <values>0 5e-3</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="2" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 3013400 4e6 5e6 6e6 8e6 1e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="7,7" compression="diagonal">
            <values>0 4.5e-3 4.5e-3 4.5e-3 4.5e-3 4.5e-3 4.5e-3</values></array></gridded2d></covarianceMatrix>
      <shortRangeSelfScalingVariance label="3" type="absolute" dependenceOnProcessedGroupWidth="inverse">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 3013400 4e6 5e6 6e6 8e6 1e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="7,7" compression="diagonal">
            <values>0 3.5972e-6 7.369e-6 3.8185e-6 1.376e-6 1.353e-7 1.1602e-8</values></array></gridded2d></shortRangeSelfScalingVariance></mixed></section>'''.split(
                '\n'))

    def test_toXML_strList_36(self):
        self.maxDiff = None
        self.assertXMLListsEqual(
            FeCovariance.covarianceSections[36].toXML_strList(formatVersion=GNDS_formatVersionModule.version_1_10),
            '''<section label="H3 + Mn54 [inclusive]">
    <rowData ENDF_MFMT="33,105" href="$reactions#/reactionSuite/reactions/reaction[@label='H3 + Mn54 [inclusive]']/crossSection/regions1d[@label='eval']"/>
    <mixed label="eval">
      <covarianceMatrix label="0" type="absolute">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="2,2" compression="diagonal">
            <values>0 9e-8</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="1" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="2,2" compression="diagonal">
            <values>0 8e-2</values></array></gridded2d></covarianceMatrix>
      <covarianceMatrix label="2" type="relative">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="5,5" compression="diagonal">
            <values>0 0.072 0.072 0.072 0.072</values></array></gridded2d></covarianceMatrix>
      <shortRangeSelfScalingVariance label="3" type="absolute" dependenceOnProcessedGroupWidth="inverse">
        <gridded2d>
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values>1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link href="../../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="5,5" compression="diagonal">
            <values>0 7.2e-14 1.5335e-14 2.9892e-12 1.9177e-9</values></array></gridded2d></shortRangeSelfScalingVariance></mixed></section>'''.split(
                '\n'))

    def test_check(self):

        listAssertEqual = getattr(self, 'assertCountEqual', None)
        if (listAssertEqual is None): listAssertEqual = getattr(self, 'assertItemsEqual')

        listAssertEqual(FeCovariance.covarianceSections[36].check(
            {'checkUncLimits': False, 'negativeEigenTolerance': 1e-8, 'eigenvalueRatioTolerance': 1e8}), [])

    def test_addComponent(self):
        """Test adding components the two standard ways"""
        ABC = covariances.mixed.MixedForm()
        ABC.addComponent(self.a)
        ABC.addComponent(self.b)
        ABC.addComponent(self.c)
        self.assertXMLListsEqual(self.abc.toXML_strList(formatVersion=GNDS_formatVersionModule.version_1_10),
                                 ABC.toXML_strList(formatVersion=GNDS_formatVersionModule.version_1_10))

    def test_constructArray(self):

        listAssertEqual = getattr(self, 'assertCountEqual', None)
        if (listAssertEqual is None): listAssertEqual = getattr(self, 'assertItemsEqual')

        listAssertEqual(self.a.matrix.array.constructArray().ravel(),
                        self.abc.components[0].matrix.array.constructArray().ravel())

    @unittest.skip("FIXME")
    def test_group(self):
        """test grouping fun on a MixedForm"""

        # build the mixed matrix & make sure is kosher
        self.assertEqual(list(map(str, self.abc.check(
            {'checkUncLimits': False, 'negativeEigenTolerance': 0.0001, 'eigenvalueRatioTolerance': 0.0001}))), [])
        self.assertXMLListsEqual(self.abc.toXML_strList(),
                                 """<mixed label="None">
  <covarianceMatrix label="eval" type="relative">
    <gridded2d>
      <axes>
        <grid index="2" label="row_energy_bounds" unit="MeV" style="boundaries">
          <values>1e-7 0.11109 1.3534 19.64</values></grid>
        <grid index="1" label="column_energy_bounds" unit="MeV" style="link">
          <link href="../../grid[@index='2']/values"/></grid>
        <axis index="0" label="matrix_elements" unit="b**2"/></axes>
      <array shape="3,3" symmetry="lower">
        <values>4 1 9 0 0 25</values></array></gridded2d></covarianceMatrix>
  <covarianceMatrix label="eval" type="relative">
    <gridded2d>
      <axes>
        <grid index="2" label="row_energy_bounds" unit="MeV" style="boundaries">
          <values>1e-5 0.1 1 20</values></grid>
        <grid index="1" label="column_energy_bounds" unit="MeV" style="link">
          <link href="../../grid[@index='2']/values"/></grid>
        <axis index="0" label="matrix_elements" unit="b**2"/></axes>
      <array shape="3,3" symmetry="lower">
        <values>4 1 9 0 0 25</values></array></gridded2d></covarianceMatrix>
  <covarianceMatrix label="eval" type="relative">
    <gridded2d>
      <axes>
        <grid index="2" label="row_energy_bounds" unit="MeV" style="boundaries">
          <values>1e-7 0.06738 0.11109 1.3534</values></grid>
        <grid index="1" label="column_energy_bounds" unit="MeV" style="link">
          <link href="../../grid[@index='2']/values"/></grid>
        <axis index="0" label="matrix_elements" unit="b**2"/></axes>
      <array shape="3,3" symmetry="lower">
        <values>4 1 9 0 0 25</values></array></gridded2d></covarianceMatrix></mixed>""".split('\n'))

        abc_c = self.abc.toCovarianceMatrix()
        self.assertEqual(list(map(str, abc_c.check(
            {'checkUncLimits': False, 'negativeEigenTolerance': 0.0001, 'eigenvalueRatioTolerance': 0.0001}))), [])
        self.assertXMLListsEqual(abc_c.toXML_strList(),
                                 """<covarianceMatrix label="composed" type="relative">
  <gridded2d>
    <axes>
      <grid index="2" label="row_energy_bounds" unit="MeV" style="boundaries">
        <values>1e-7 1e-5 0.06738 0.1 0.11109 1 1.3534 19.64 20</values></grid>
      <grid index="1" label="column_energy_bounds" unit="MeV" style="link">
        <link href="../../grid[@index='2']/values"/></grid>
      <axis index="0" label="matrix_elements" unit="b**2"/></axes>
    <array shape="8,8" symmetry="lower">
      <values>8 8 12 5 9 17 5 6 14 22 1 2 2 10 43 1 1 1 1 34 59 0 0 0 0 0 25 50 0 0 0 0 0 25 25 25</values></array></gridded2d></covarianceMatrix>""".split(
                                     '\n'))

        self.maxDiff = None
        g = abc_c.group(groupBoundaries=(self.__groupBoundaries, self.__groupBoundaries),
                        groupUnit=(self.__groupUnit, self.__groupUnit))
        g.convertAxesToUnits(('b**2', 'MeV', 'MeV'))
        self.assertTrue(g.matrix.axes[0].unit, 'MeV')
        self.assertTrue(g.matrix.axes[1].unit, 'MeV')
        self.assertTrue(g.matrix.axes[2].unit, 'b**2')
        self.assertEqual(list(map(str, g.check(
            {'checkUncLimits': False, 'negativeEigenTolerance': 0.0001, 'eigenvalueRatioTolerance': 0.0001}))), [])
        if True:
            self.assertXMLListsEqual(g.toXML_strList(),
                                     """<covarianceMatrix label="composed" type="relative">
  <gridded2d>
    <axes>
      <grid index="2" label="row_energy_bounds" unit="MeV" style="boundaries">
        <values>1e-7 0.11109 1.3534 19.64</values></grid>
      <grid index="1" label="column_energy_bounds" unit="MeV" style="link">
        <link href="../../grid[@index='2']/values"/></grid>
      <axis index="0" label="matrix_elements" unit="b**2"/></axes>
    <array shape="33,33" symmetry="lower">
      <values>2.09046853425614e-31 1.50619069577192e-15 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868
        10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868
        10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 1.50619069577192e-15 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 10.8521624451868 3.17404481019479e-16 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 40.6309282309796
        3.17404481019479e-16 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 40.6309282309796 40.6309282309796 3.17404481019479e-16 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 40.6309282309796 40.6309282309796 40.6309282309796 3.17404481019479e-16 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 40.6309282309796 40.6309282309796 40.6309282309796 40.6309282309796 3.17404481019479e-16 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969 2.28691160987969
        2.28691160987969 40.6309282309796 40.6309282309796 40.6309282309796 40.6309282309796 40.6309282309796 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 50 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 50 50 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 50 50 50 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 50 50 50 50 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 7.11175149519846 50 50 50 50 50</values></array></gridded2d></covarianceMatrix>""".split(
                                         '\n'))
        if False:
            # self.abc.plot(xlog=True, ylog=True, title='abc')
            # self.abc.plot(xlog=False, ylog=False, title='abc')
            self.a.plot(title="a")
            self.b.plot(xlog=True, ylog=True, title='b')
            self.c.plot(title='c')
            g.plot(xlog=True, ylog=True, title='abc grouped')

    def test_getMatchingComponent_0(self):
        self.assertXMLListsEqual(FeCovariance.covarianceSections[33]['eval'].getMatchingComponent(
            rowBounds=(1.e-5, 2.e7),
            columnBounds=(1.e-5, 2.e7)).toXML_strList(formatVersion=GNDS_formatVersionModule.version_1_10),
                                 ['<covarianceMatrix label="1" type="relative">',
                                  '<gridded2d>',
                                  '<axes>',
                                  '<grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">',
                                  '<values>1e-5 8.5e5 2e6 3e6 4e6 5e6 6e6 7e6 1e7 2e7</values></grid>',
                                  '<grid index="1" label="column_energy_bounds" unit="eV" style="link">',
                                  '<link href="../../grid[@index=\'2\']/values"/></grid>',
                                  '<axis index="0" label="matrix_elements" unit=""/></axes>',
                                  '<array shape="9,9" compression="diagonal">',
                                  '<values>0 0.0396 0.0891 0.0891 0.0891 0.0891 0.0891 0.0891 0.3782</values></array></gridded2d></covarianceMatrix>'])

    def test_getMatchingComponent_0_stripped(self):
        x = FeCovariance.covarianceSections[33]['eval'].getMatchingComponent(
            rowBounds=(1.e-5, 2.e7),
            columnBounds=(1.e-5, 2.e7))
        x.removeExtraZeros()
        self.assertXMLListsEqual(x.toXML_strList(formatVersion=GNDS_formatVersionModule.version_2_0),
                                 ['<covarianceMatrix label="1" type="relative">',
                                  '<gridded2d>',
                                  '<axes>',
                                  '<grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">',
                                  '<values>8.5e5 2e6 3e6 4e6 5e6 6e6 7e6 1e7 2e7</values></grid>',
                                  '<grid index="1" label="column_energy_bounds" unit="eV" style="boundaries">',
                                  '<link href="../../grid[@index=\'2\']/values"/></grid>',
                                  '<axis index="0" label="matrix_elements" unit=""/></axes>',
                                  '<array shape="8,8" symmetry="lower">',
                                  '<values>0.0396 0 0.0891 0 0 0.0891 0 0 0 0.0891 0 0 0 0 0.0891 0 0 0 0 0 0.0891 0 0 0 0 0 0 0.0891 0 0 0 0 0 0 0 0.3782</values></array></gridded2d></covarianceMatrix>'])

    def test_getMatchingComponent_1(self):
        """Fails because there is no component of section 1 with such bounds"""
        self.assertRaises(ValueError,
                          FeCovariance.covarianceSections[1]['eval'].getMatchingComponent,
                          rowBounds=(1.e-5, 2.e7),
                          columnBounds=(1.e-5, 2.e7))

    def test_getMatchingComponent_2(self):
        self.maxDiff = None
        self.assertXMLListsEqual(FeCovariance.covarianceSections[1]['eval'].getMatchingComponent(
            rowBounds=(1.e-5, 850636),
            columnBounds=(1.e-5, 850636)).toXML_strList(formatVersion=GNDS_formatVersionModule.version_1_10),
                                 ['<covarianceMatrix label="1" type="relative">',
                                  '<gridded2d>',
                                  '<axes>',
                                  '<grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">',
                                  '<values>1e-5 20 3e2 33830.5 76772 102950 152676 204928 238299 446595 483428 518904 553281 850636</values></grid>',
                                  '<grid index="1" label="column_energy_bounds" unit="eV" style="link">',
                                  '<link href="../../grid[@index=\'2\']/values"/></grid>',
                                  '<axis index="0" label="matrix_elements" unit=""/></axes>',
                                  '<array shape="13,13" symmetry="lower">',
                                  '<values>1.6e-3 2.4e-3 3.6e-3 0 0 3.176954e-3 0 0 3.468456e-3 0.01514682 0 0 2.606664e-3 5.691679e-3 8.554986e-3 0 0 2.975089e-3 6.496139e-3 4.882072e-3 0.0111442 0 0 2.995769e-3 6.541293e-3 4.916007e-3 5.610834e-3 0.01129967 0 0 3.582893e-3 7.823286e-3 5.879469e-3 6.710471e-3 6.757115e-3 0.01616281 0 0 3.334617e-3 7.281172e-3 5.472051e-3 6.245469e-3 6.288881e-3 7.521405e-3 0.01400042 0 0 3.325916e-3 7.262173e-3 5.457774e-3 6.229173e-3 6.272472e-3 7.50178e-3 6.981945e-3 0.01392745 0 0 3.059465e-3 6.680374e-3 5.020531e-3 5.730131e-3 5.769962e-3 6.900785e-3 6.422595e-3 6.405837e-3 0.01178528 0 0 3.185348e-3 6.955242e-3 5.227104e-3 5.965901e-3 6.00737e-3 7.184721e-3 6.686857e-3 6.669409e-3 6.135098e-3 0.01277506 0 0 2.919549e-3 6.374865e-3 4.790931e-3 5.468079e-3 5.506088e-3 6.585196e-3 6.128875e-3 6.112884e-3 5.623158e-3 5.854526e-3 0.010732</values></array></gridded2d></covarianceMatrix>'])

    def test_makeSafeBounds(self):
        cc = copy.copy(FeCovariance.covarianceSections[0]['eval'])
        self.assertEqual(
            [c.rowBounds() for c in cc.components],
            [(1.e-5, 862270.), (1e-5, 2.e7)])
        cc.makeSafeBounds()
        self.assertEqual(
            [c.rowBounds() for c in cc.components],
            [(1.e-5, 862270.), (862270., 2.e7)])

    def test_getRowBounds(self):
        self.assertEqual(FeCovariance.covarianceSections[0]['eval'].rowBounds('eV'), (1.e-5, 2.e7))
        self.assertEqual(FeCovariance.covarianceSections[1]['eval'].rowBounds('eV'), (1.e-5, 2.e7))
        self.assertEqual(FeCovariance.covarianceSections[2]['eval'].rowBounds('eV'), (1.e-5, 2.e7))

    def test_getColumnBounds(self):
        self.assertEqual(FeCovariance.covarianceSections[0]['eval'].columnBounds('eV'), (1.e-5, 2.e7))
        self.assertEqual(FeCovariance.covarianceSections[1]['eval'].columnBounds('eV'), (1.e-5, 2.e7))
        self.assertEqual(FeCovariance.covarianceSections[2]['eval'].columnBounds('eV'), (1.e-5, 2.e7))

    @unittest.skip("FIXME")
    def test_getUncertaintyVector(self):
        self.assertEqual(
            repr(FeCovariance.covarianceSections[0]['eval'].getUncertaintyVector(
                theData=FeEvaluation.getReaction('elastic').crossSection)),
            '''   1.00000000e-05   6.72446355e-02
   4.99999995e-01   6.72446355e-02
   5.00000005e-01   1.17323442e-01
   1.99999998e+01   1.17323442e-01
   2.00000002e+01   1.25557915e-01
   2.99999997e+02   1.25557915e-01
   3.00000003e+02   1.23861794e-01
   1.14003999e+03   1.23861794e-01
   1.14004001e+03   1.14460950e-01
   1.75099998e+03   1.14460950e-01
   1.75100002e+03   1.35933822e-01
   3.38304997e+04   1.35933822e-01
   3.38305003e+04   1.74493180e-01
   3.54794996e+04   1.74493180e-01
   3.54795004e+04   1.58269027e-01
   5.36474995e+04   1.58269027e-01
   5.36475005e+04   1.44941267e-01
   7.67719992e+04   1.44941267e-01
   7.67720008e+04   1.20067219e-01
   7.79314992e+04   1.20067219e-01
   7.79315008e+04   1.15763669e-01
   8.93689991e+04   1.15763669e-01
   8.93690009e+04   1.28532700e-01
   9.45609991e+04   1.28532700e-01
   9.45610009e+04   1.42054307e-01
   1.02949999e+05   1.42054307e-01
   1.02950001e+05   1.50892810e-01
   1.04507999e+05   1.50892810e-01
   1.04508001e+05   1.58779848e-01
   1.52675998e+05   1.58779848e-01
   1.52676002e+05   1.66531078e-01
   2.03064998e+05   1.66531078e-01
   2.03065002e+05   1.36276807e-01
   2.04927998e+05   1.36276807e-01
   2.04928002e+05   1.53083337e-01
   2.38298998e+05   1.53083337e-01
   2.38299002e+05   1.45849642e-01
   2.61279997e+05   1.45849642e-01
   2.61280003e+05   3.03438379e-01
   4.46594996e+05   3.03438379e-01
   4.46595004e+05   3.08426588e-01
   4.83427995e+05   3.08426588e-01
   4.83428005e+05   3.04934075e-01
   5.18903995e+05   3.04934075e-01
   5.18904005e+05   3.06552720e-01
   5.53280994e+05   3.06552720e-01
   5.53281006e+05   2.81911263e-01
   7.94898992e+05   2.81911263e-01
   7.94899008e+05   0.00000000e+00
   8.62269991e+05   0.00000000e+00
   8.62270009e+05   9.94987437e-03
   9.99999990e+05   9.94987437e-03
   1.00000001e+06   7.45922248e-02
   1.99999998e+06   7.45922248e-02
   2.00000002e+06   3.97994975e-02
   4.99999995e+06   3.97994975e-02
   5.00000005e+06   2.98496231e-02
   9.99999990e+06   2.98496231e-02
   1.00000001e+07   1.98997487e-02
   1.49999998e+07   1.98997487e-02
   1.50000001e+07   1.98997487e-02
   2.00000000e+07   1.98997487e-02
''')

    @unittest.skip("FIXME")
    def test_toCovarianceMatrix(self):
        self.maxDiff = None
        # The values below are correct, I checked them by hand (what a pain) DAB 13 Jun 2016
        self.assertXMLListsEqual(self.abc.toCovarianceMatrix().toXML_strList(), """<covarianceMatrix label="composed" type="relative">
  <gridded2d>
    <axes>
      <grid index="2" label="row_energy_bounds" unit="MeV" style="boundaries">
        <values>1e-7 1e-5 0.06738 0.1 0.11109 1 1.3534 19.64 20</values></grid>
      <grid index="1" label="column_energy_bounds" unit="MeV" style="link">
        <link href="../../grid[@index='2']/values"/></grid>
      <axis index="0" label="matrix_elements" unit="b**2"/></axes>
    <array shape="8,8" symmetry="lower">
      <values>8 8 12 5 9 17 5 6 14 22 1 2 2 10 43 1 1 1 1 34 59 0 0 0 0 0 25 50 0 0 0 0 0 25 25 25</values></array></gridded2d></covarianceMatrix>""".split(
            '\n'))
        self.abc.plot()

    @unittest.skip("FIXME")
    def test_toAbsolute(self):
        l = FeCovariance.covarianceSections[0]['eval'].toAbsolute()
        self.assertXMLListsEqual(l.toXML_strList(), """<mixed label="eval">
  <covarianceMatrix label="composed" type="absolute">
    <gridded2d>
      <axes>
        <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
          <values>1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899 850636</values></grid>
        <grid index="1" label="row_energy_bounds" unit="eV" style="link">
          <link href="../../grid[@index='2']/values"/></grid>
        <axis index="0" label="matrix_elements" unit="b**2"/></axes>
      <array shape="25,25" symmetry="lower">
        <values>4.521841e-3 7.561844e-3 0.01376479 8.361844e-3 0.01456479 0.01576479 5.961844e-3 0.01216479 0.01216479 0.015341744 0 0 0 3.176954e-3 0.013101309 0 0 0 3.176954e-3 9.338387e-3 0.018478004 0 0 0 3.468456e-3 9.629889e-3 0.018769506 0.03044787 0 0 0 3.468456e-3 8.425108e-3 9.623028e-3 0.021301392 0.025049085 0 0 0 3.468456e-3 7.281858e-3 8.203479e-3 0.019881843 0.018955976 0.021007971 0 0 0 2.606664e-3 6.420066e-3 7.341687e-3 0.010426702 9.500835e-3 0.01155283 0.014416137 0 0 0 2.606664e-3 6.074221e-3 6.912258e-3 9.997273e-3 9.155375e-3 8.356475e-3 0.011219782 0.013401227 0 0 0 2.606664e-3 7.052283e-3 8.126697e-3 0.011211712 0.010132348 9.10811e-3 0.011971417 0.011661574 0.016520655 0 0 0 2.606664e-3 7.977071e-3 9.274988e-3 0.012360003 0.011056106 9.818803e-3 0.01268211 0.012307813 0.013366339 0.020179426 0 0 0 2.975089e-3 8.345496e-3 9.643413e-3 0.013164463 0.011860566 0.010623263
          9.009196e-3 8.634899e-3 9.693425e-3 0.016506512 0.02276864 0 0 0 2.975089e-3 8.882803e-3 0.010310575 0.013831625 0.012397274 0.01103618 9.422113e-3 9.010368e-3 0.010174798 0.011275803 0.017537931 0.02521104 0 0 0 2.995769e-3 9.381029e-3 0.010924214 0.014469738 0.012919443 0.011448326 9.82304e-3 9.378011e-3 0.010636568 0.011826572 0.012521399 0.013212797 0.0277326 0 0 0 2.995769e-3 7.243324e-3 8.269871e-3 0.011815395 0.010784119 9.805513e-3 8.180227e-3 7.884188e-3 8.721396e-3 9.513002e-3 0.010207829 0.010667755 0.016765365 0.018571368 0 0 0 3.582893e-3 7.830448e-3 8.856995e-3 0.013097388 0.012066112 0.011087506 9.143689e-3 8.84765e-3 9.684858e-3 0.010476464 0.011307466 0.011767392 0.01222281 0.014028813 0.023434508 0 0 0 3.334617e-3 7.582172e-3 8.608719e-3 0.012555274 0.011523998 0.010545392 8.736271e-3 8.440232e-3 9.27744e-3 0.010069046 0.010842464 0.01130239 0.011754576 0.013560579 0.014793103 0.021272118 0 0 0 3.334617e-3 0.017252577 0.020616257 0.024562812 0.021183632 0.017977032 0.016167911
          0.015197887 0.017941161 0.020535021 0.021308439 0.022815479 0.024198321 0.018202461 0.019434985 0.025914 0.09207485 0 0 0 3.325916e-3 0.017519686 0.020950026 0.024886283 0.021440143 0.018170003 0.016365604 0.015376347 0.018173984 0.020819244 0.021590643 0.023127553 0.024536822 0.018422142 0.01965145 0.019131615 0.046792765 0.09512696 0 0 0 3.059465e-3 0.017253235 0.020683575 0.024304484 0.020858344 0.017588204 0.015928361 0.014939104 0.017736741 0.020382001 0.021091601 0.022628511 0.024034312 0.017919632 0.019050455 0.018572265 0.046233415 0.087605347 0.09298479 0 0 0 3.185348e-3 0.017379118 0.020809458 0.024579352 0.021133212 0.017863072 0.016134934 0.015145677 0.017943314 0.020588574 0.021327371 0.022864281 0.02427172 0.01815704 0.019334391 0.018836527 0.046497677 0.087868919 0.087334608 0.09397457 0 0 0 2.919549e-3 0.016961699 0.020355389 0.023810705 0.020401375 0.017166165 0.015582231 0.01460355 0.017371301 0.019988301 0.020665449 0.022185939 0.023575328 0.017525978 0.018605086 0.018148765 0.045514415 0.046278934 0.045789208 0.046020576 0.09020596
          0 0 0 2.919549e-3 0.017258749 0.020724229 0.024179545 0.020698095 0.017394445 0.015810511 0.014811121 0.017637431 0.020309791 0.020986939 0.022539589 0.023957568 0.017780238 0.018859346 0.018403025 0.046347575 0.047128604 0.046638878 0.046870246 0.05130957 0.0936039</values></array></gridded2d></covarianceMatrix>
  <covarianceMatrix label="1" type="absolute">
    <gridded2d>
      <axes>
        <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
          <values>862270 1e6 2e6 5e6 1e7 1.5e7 2e7</values></grid>
        <grid index="1" label="column_energy_bounds" unit="eV" style="link">
          <link href="../../grid[@index='2']/values"/></grid>
        <axis index="0" label="matrix_elements" unit="b**2"/></axes>
      <array shape="6,6" symmetry="lower">
        <values>9.9e-5 0 5.564e-3 0 0 1.584e-3 0 0 0 8.91e-4 0 0 0 0 3.96e-4 0 0 0 0 0 3.96e-4</values></array></gridded2d></covarianceMatrix></mixed>""".split(
            '\n'))

    @unittest.skip("FIXME")
    def test_toRelative(self):
        l = FeCovariance.covarianceSections[0]['eval'].toRelative()
        self.assertXMLListsEqual(l.toXML_strList(), """<mixed label="eval">
  <covarianceMatrix label="composed" type="relative">
    <gridded2d>
      <axes>
        <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
          <values>1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899 850636</values></grid>
        <grid index="1" label="row_energy_bounds" unit="eV" style="link">
          <link href="../../grid[@index='2']/values"/></grid>
        <axis index="0" label="matrix_elements" unit=""/></axes>
      <array shape="25,25" symmetry="lower">
        <values>4.521841e-3 7.561844e-3 0.01376479 8.361844e-3 0.01456479 0.01576479 5.961844e-3 0.01216479 0.01216479 0.015341744 0 0 0 3.176954e-3 0.013101309 0 0 0 3.176954e-3 9.338387e-3 0.018478004 0 0 0 3.468456e-3 9.629889e-3 0.018769506 0.03044787 0 0 0 3.468456e-3 8.425108e-3 9.623028e-3 0.021301392 0.025049085 0 0 0 3.468456e-3 7.281858e-3 8.203479e-3 0.019881843 0.018955976 0.021007971 0 0 0 2.606664e-3 6.420066e-3 7.341687e-3 0.010426702 9.500835e-3 0.01155283 0.014416137 0 0 0 2.606664e-3 6.074221e-3 6.912258e-3 9.997273e-3 9.155375e-3 8.356475e-3 0.011219782 0.013401227 0 0 0 2.606664e-3 7.052283e-3 8.126697e-3 0.011211712 0.010132348 9.10811e-3 0.011971417 0.011661574 0.016520655 0 0 0 2.606664e-3 7.977071e-3 9.274988e-3 0.012360003 0.011056106 9.818803e-3 0.01268211 0.012307813 0.013366339 0.020179426 0 0 0 2.975089e-3 8.345496e-3 9.643413e-3 0.013164463 0.011860566 0.010623263
          9.009196e-3 8.634899e-3 9.693425e-3 0.016506512 0.02276864 0 0 0 2.975089e-3 8.882803e-3 0.010310575 0.013831625 0.012397274 0.01103618 9.422113e-3 9.010368e-3 0.010174798 0.011275803 0.017537931 0.02521104 0 0 0 2.995769e-3 9.381029e-3 0.010924214 0.014469738 0.012919443 0.011448326 9.82304e-3 9.378011e-3 0.010636568 0.011826572 0.012521399 0.013212797 0.0277326 0 0 0 2.995769e-3 7.243324e-3 8.269871e-3 0.011815395 0.010784119 9.805513e-3 8.180227e-3 7.884188e-3 8.721396e-3 9.513002e-3 0.010207829 0.010667755 0.016765365 0.018571368 0 0 0 3.582893e-3 7.830448e-3 8.856995e-3 0.013097388 0.012066112 0.011087506 9.143689e-3 8.84765e-3 9.684858e-3 0.010476464 0.011307466 0.011767392 0.01222281 0.014028813 0.023434508 0 0 0 3.334617e-3 7.582172e-3 8.608719e-3 0.012555274 0.011523998 0.010545392 8.736271e-3 8.440232e-3 9.27744e-3 0.010069046 0.010842464 0.01130239 0.011754576 0.013560579 0.014793103 0.021272118 0 0 0 3.334617e-3 0.017252577 0.020616257 0.024562812 0.021183632 0.017977032 0.016167911
          0.015197887 0.017941161 0.020535021 0.021308439 0.022815479 0.024198321 0.018202461 0.019434985 0.025914 0.09207485 0 0 0 3.325916e-3 0.017519686 0.020950026 0.024886283 0.021440143 0.018170003 0.016365604 0.015376347 0.018173984 0.020819244 0.021590643 0.023127553 0.024536822 0.018422142 0.01965145 0.019131615 0.046792765 0.09512696 0 0 0 3.059465e-3 0.017253235 0.020683575 0.024304484 0.020858344 0.017588204 0.015928361 0.014939104 0.017736741 0.020382001 0.021091601 0.022628511 0.024034312 0.017919632 0.019050455 0.018572265 0.046233415 0.087605347 0.09298479 0 0 0 3.185348e-3 0.017379118 0.020809458 0.024579352 0.021133212 0.017863072 0.016134934 0.015145677 0.017943314 0.020588574 0.021327371 0.022864281 0.02427172 0.01815704 0.019334391 0.018836527 0.046497677 0.087868919 0.087334608 0.09397457 0 0 0 2.919549e-3 0.016961699 0.020355389 0.023810705 0.020401375 0.017166165 0.015582231 0.01460355 0.017371301 0.019988301 0.020665449 0.022185939 0.023575328 0.017525978 0.018605086 0.018148765 0.045514415 0.046278934 0.045789208 0.046020576 0.09020596
          0 0 0 2.919549e-3 0.017258749 0.020724229 0.024179545 0.020698095 0.017394445 0.015810511 0.014811121 0.017637431 0.020309791 0.020986939 0.022539589 0.023957568 0.017780238 0.018859346 0.018403025 0.046347575 0.047128604 0.046638878 0.046870246 0.05130957 0.0936039</values></array></gridded2d></covarianceMatrix>
  <covarianceMatrix label="1" type="relative">
    <gridded2d>
      <axes>
        <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
          <values>862270 1e6 2e6 5e6 1e7 1.5e7 2e7</values></grid>
        <grid index="1" label="column_energy_bounds" unit="eV" style="link">
          <link href="../../grid[@index='2']/values"/></grid>
        <axis index="0" label="matrix_elements" unit=""/></axes>
      <array shape="6,6" symmetry="lower">
        <values>9.9e-5 0 5.564e-3 0 0 1.584e-3 0 0 0 8.91e-4 0 0 0 0 3.96e-4 0 0 0 0 0 3.96e-4</values></array></gridded2d></covarianceMatrix></mixed>""".split(
            '\n'))


if __name__ == "__main__":
    unittest.main()
