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

Note, the tests below operate on the one and only copy of the covariance data loaded below.  Make sure
each test only fools with one section otherwise you'll be chasing ghost bugs:

    * run 1: test A changes file, test B depends on results from test A

    * run 2 (different order): test B fails because test A didn't change the files in an "expected" way

"""

import unittest
from pqu import PQU
from fudge.gnd import covariances
import xData.axes as axesModule
import xData.link as linkModule
import xData.gridded as griddedModule
import xData.values as valuesModule
import xData.array as arrayModule

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


class Test_summed( TestCaseBase ):

    def setUp(self):
        # ...................... example matrix 'a' ......................
        axes = axesModule.axes(
            labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'),
                         2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.grid(axes[2].label, axes[2].index, axes[2].unit,
                                  style=axesModule.boundariesGridToken,
                                  values=valuesModule.values([1.0000E-07, 1.1109E-01, 1.3534E+00, 1.9640E+01]))
        axes[1] = axesModule.grid(axes[1].label, axes[1].index, axes[1].unit,
                                  style=axesModule.linkGridToken,
                                  values=linkModule.link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.symmetryLowerToken)
        self.a = covariances.covarianceMatrix('eval', matrix=griddedModule.gridded2d(axes, myMatrix),
                                         type=covariances.tokens.relativeToken)

        # ...................... example matrix 'b' ......................
        axes = axesModule.axes(
            labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'),
                         2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.grid(axes[2].label, axes[2].index, axes[2].unit,
                                  style=axesModule.boundariesGridToken,
                                  values=valuesModule.values([1.0e-5, 0.100, 1.0, 20.0]))
        axes[1] = axesModule.grid(axes[1].label, axes[1].index, axes[1].unit,
                                  style=axesModule.linkGridToken,
                                  values=linkModule.link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.symmetryLowerToken)
        self.b = covariances.covarianceMatrix('eval', matrix=griddedModule.gridded2d(axes, myMatrix),
                                         type=covariances.tokens.relativeToken)

        # ...................... example matrix 'c' ......................
        axes = axesModule.axes(
            labelsUnits={0: ('matrix_elements', 'b**2'), 1: ('column_energy_bounds', 'MeV'),
                         2: ('row_energy_bounds', 'MeV')})
        axes[2] = axesModule.grid(axes[2].label, axes[2].index, axes[2].unit,
                                  style=axesModule.boundariesGridToken,
                                  values=valuesModule.values([1.0000E-07, 6.7380E-02, 1.1109E-01, 1.3534E+00]))
        axes[1] = axesModule.grid(axes[1].label, axes[1].index, axes[1].unit,
                                  style=axesModule.linkGridToken,
                                  values=linkModule.link(link=axes[2].values, relative=True))
        myMatrix = arrayModule.full((3, 3), [4.0, 1.0, 9.0, 0.0, 0.0, 25.0], symmetry=arrayModule.symmetryLowerToken)
        self.c = covariances.covarianceMatrix('eval', matrix=griddedModule.gridded2d(axes, myMatrix),
                                         type=covariances.tokens.relativeToken)

        # ...................... combine them for example matrix 'abc' ......................
        abc=covariances.mixed.mixedForm(components=[self.a, self.b, self.c])
        # FIXME: learn how to add abc to a section & to a covarianceSuite!, sumabc is built wrong!

        # ...................... 'sumabc' is just a way to exercise the summed class ......................
        bds=abc.getRowBounds()
        self.sumabc=covariances.summed.summedCovariance(label='test', domainMin=float(bds[0]), domainMax=float(bds[1]),
            domainUnit=abc[0].matrix.axes[-1].unit, pointerList=[linkModule.link(link=abc, path=abc.toXLink(), coefficient=1.0)])

    @unittest.skip("FIXME")
    def test__getitem__(self):
        self.assertXMLListsEqual( self.sumabc.pointerList[0].toXMLList(),['<link coefficient="1.0" xlink:href="/mixed"/>'] )

    def test__len__(self): 
        self.assertEqual( len(self.sumabc.pointerList), 1 )
        
    @unittest.skip("FIXME")
    def test_toXMLList(self):
        print '\n'.join(self.sumabc.toXMLList())
        self.assertXMLListsEqual( self.sumabc.toXMLList(),"""<sum lowerBound="1e-7 MeV" upperBound="20 MeV">
  <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->
  <link coefficient="1.0" xlink:href="/mixed"/></sum>""".split('\n') )

    @unittest.skip("FIXME")
    def test_getReferredCovariance(self):
        pass

    def test_getRowBounds(self):
        bounds = tuple([PQU.PQU(bound, self.sumabc.domainUnit) for bound in self.sumabc.getRowBounds()])
        self.assertEqual( bounds, (PQU.PQU( "1.e-7 MeV" ), PQU.PQU( "20. MeV" )))

    def test_getColumnBounds(self):
        bounds = tuple([PQU.PQU(bound, self.sumabc.domainUnit) for bound in self.sumabc.getColumnBounds()])
        self.assertEqual( bounds, (PQU.PQU( "1.e-7 MeV" ), PQU.PQU( "20. MeV" )))

    @unittest.skip("FIXME")
    def test_toCovarianceMatrix(self):
        self.assertXMLListsEqual(self.sumabc.toCovarianceMatrix().toXMLList(),"""<covarianceMatrix label="composed" type="relative">
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
