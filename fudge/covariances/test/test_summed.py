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

import unittest
from pqu import PQU
from fudge import covariances
from fudge.covariances import enums as covarianceEnumsModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import link as linkModule
from xData import gridded as griddedModule
from xData import values as valuesModule
from xData import xDataArray as arrayModule

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


class Test_summed( TestCaseBase ):

    def setUp(self):
        # ...................... example matrix 'a' ......................
        axes = axesModule.Axes(3, labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'), 2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit, style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values([1.0000E-07, 1.1109E-01, 1.3534E+00, 1.9640E+01]))
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit, style=xDataEnumsModule.GridStyle.none,
                                  values=linkModule.Link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.Full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.Symmetry.lower)
        self.a = covariances.covarianceMatrix.CovarianceMatrix('eval', matrix=griddedModule.Gridded2d(axes, myMatrix),
                                         type=covarianceEnumsModule.Type.relative)

        # ...................... example matrix 'b' ......................
        axes = axesModule.Axes(3, labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'), 2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit, style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values([1.0e-5, 0.100, 1.0, 20.0]))
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit, style=xDataEnumsModule.GridStyle.none,
                                  values=linkModule.Link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.Full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.Symmetry.lower)
        self.b = covariances.covarianceMatrix.CovarianceMatrix('eval', matrix=griddedModule.Gridded2d(axes, myMatrix),
                                         type=covarianceEnumsModule.Type.relative)

        # ...................... example matrix 'c' ......................
        axes = axesModule.Axes(3, labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'), 2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit, style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values([1.0000E-07, 6.7380E-02, 1.1109E-01, 1.3534E+00]))
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit, style=xDataEnumsModule.GridStyle.none,
                                  values=linkModule.Link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.Full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.Symmetry.lower)
        self.c = covariances.covarianceMatrix.CovarianceMatrix('eval', matrix=griddedModule.Gridded2d(axes, myMatrix),
                                         type=covarianceEnumsModule.Type.relative)

        # ...................... combine them for example matrix 'abc' ......................
        abc=covariances.mixed.MixedForm(components=[self.a, self.b, self.c])
        # FIXME: learn how to add abc to a section & to a covarianceSuite!, sumabc is built wrong!

        # ...................... 'sumabc' is just a way to exercise the summed class ......................
        bds=abc.rowBounds()
        self.sumabc=covariances.summed.SummedCovariance(label='test', domainMin=float(bds[0]), domainMax=float(bds[1]),
            domainUnit=abc[0].matrix.axes[-1].unit, summands=[covariances.summed.Summand(link=abc, path=abc.toXLink(), coefficient=1.0)])

    @unittest.skip("FIXME")
    def test__getitem__(self):
        self.assertXMLListsEqual( self.sumabc.summands[0].toXML_strList(),['<link coefficient="1.0" xlink:href="/mixed"/>'] )

    def test__len__(self):
        self.assertEqual( len(self.sumabc.summands), 1 )

    @unittest.skip("FIXME")
    def test_toXML_strList(self):
        print('\n'.join(self.sumabc.toXML_strList()))
        self.assertXMLListsEqual( self.sumabc.toXML_strList(),"""<sum lowerBound="1e-7 MeV" upperBound="20 MeV">
  <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->
  <link coefficient="1.0" xlink:href="/mixed"/></sum>""".split('\n') )

    @unittest.skip("FIXME")
    def test_getReferredCovariance(self):
        pass

    def test_getRowBounds(self):
        bounds = tuple([PQU.PQU(bound, self.sumabc.domainUnit) for bound in self.sumabc.rowBounds()])
        self.assertEqual( bounds, (PQU.PQU( "1.e-7 MeV" ), PQU.PQU( "20. MeV" )))

    def test_getColumnBounds(self):
        bounds = tuple([PQU.PQU(bound, self.sumabc.domainUnit) for bound in self.sumabc.columnBounds()])
        self.assertEqual( bounds, (PQU.PQU( "1.e-7 MeV" ), PQU.PQU( "20. MeV" )))

    @unittest.skip("FIXME")
    def test_toCovarianceMatrix(self):
        self.assertXMLListsEqual(self.sumabc.toCovarianceMatrix().toXML_strList(),"""<covarianceMatrix label="composed" type="relative">
  <gridded2d>
    <axes>
      <grid index="2" label="row_energy_bounds" unit="MeV" style="boundaries">
        <values length="9">1e-7 1e-5 0.06738 0.1 0.11109 1 1.3534 19.64 20</values></grid>
      <grid index="1" label="column_energy_bounds" unit="MeV" style="link">
        <link xlink:href="../../grid[@index='2']/values"/></grid>
      <axis index="0" label="matrix_elements" unit="b**2"/></axes>
    <array shape="8,8" symmetry="lower">
      <values length="36">8 8 12 5 9 17 5 6 14 22 1 2 2 10 43 1 1 1 1 34 59 0 0 0 0 0 25 50 0 0 0 0 0 25 25 25</values></array></gridded2d></covarianceMatrix>""".split('\n'))

    @unittest.skip("needs real data, so the getUncertainty test in mixed does double duty and tests summed too")
    def test_getUncertaintyVector(self):
        pass


if __name__=="__main__":
    unittest.main()
