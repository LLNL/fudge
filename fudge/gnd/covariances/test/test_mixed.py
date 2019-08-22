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
from pqu import PQU
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


class Test_mixed( TestCaseBase ):
    
    def test__getitem__(self):
        self.assertXMLListsEqual(
            FeCovariance[36]['eval'].toXMLList(), '''<mixed label="eval">
      <covarianceMatrix index="0" type="absolute">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="3">1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="3,3" compression="diagonal">
            <values length="3">0 9e-8 0</values></array></gridded></covarianceMatrix>
      <covarianceMatrix index="1" type="relative">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="3">1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="3,3" compression="diagonal">
            <values length="3">0 8e-2 0</values></array></gridded></covarianceMatrix>
      <covarianceMatrix index="2" type="relative">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="6">1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="6,6" compression="diagonal">
            <values length="6">0 0.072 0.072 0.072 0.072 0</values></array></gridded></covarianceMatrix>
      <covarianceMatrix index="3" type="relative" ENDFconversionFlag="LB=8">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="6">1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="6,6" compression="diagonal">
            <values length="6">0 7.2e-14 1.5335e-14 2.9892e-12 1.9177e-9 0</values></array></gridded></covarianceMatrix></mixed>'''.split('\n') )

    def test__len__(self): 
        self.assertEqual( len(FeCovariance), 40 )
        
    def test_toXMLList(self):
        self.assertXMLListsEqual( FeCovariance[36].toXMLList(), '''<section label="36" id="H3 + Mn54_s">
    <rowData ENDF_MFMT="33,105" xlink:href="/reactionSuite/reactions/reaction[@label='33']/crossSection/regions[@label='eval']"/>
    <mixed label="eval">
      <covarianceMatrix index="0" type="absolute">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="3">1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit="b**2"/></axes>
          <array shape="3,3" compression="diagonal">
            <values length="3">0 9e-8 0</values></array></gridded></covarianceMatrix>
      <covarianceMatrix index="1" type="relative">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="3">1e-5 1.2143e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="3,3" compression="diagonal">
            <values length="3">0 8e-2 0</values></array></gridded></covarianceMatrix>
      <covarianceMatrix index="2" type="relative">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="6">1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="6,6" compression="diagonal">
            <values length="6">0 0.072 0.072 0.072 0.072 0</values></array></gridded></covarianceMatrix>
      <covarianceMatrix index="3" type="relative" ENDFconversionFlag="LB=8">
        <gridded dimension="2">
          <axes>
            <grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">
              <values length="6">1e-5 1.2143e7 1.3e7 1.45e7 1.75e7 2e7</values></grid>
            <grid index="1" label="column_energy_bounds" unit="eV" style="link">
              <link xlink:href="../grid[@index='2']/values"/></grid>
            <axis index="0" label="matrix_elements" unit=""/></axes>
          <array shape="6,6" compression="diagonal">
            <values length="6">0 7.2e-14 1.5335e-14 2.9892e-12 1.9177e-9 0</values></array></gridded></covarianceMatrix></mixed></section>'''.split('\n') )

    def test_check(self):
        self.assertItemsEqual( FeCovariance[36].check({
            'checkUncLimits':False,
            'negativeEigenTolerance':1e-8,
            'eigenvalueRatioTolerance':1e8}), [] )

    def test_fix(self): pass
    
    def test_plot(self):pass
    
    def test_addComponent(self): pass
    
    def test_getMatchingComponent_0(self):
        self.assertXMLListsEqual( FeCovariance[33]['eval'].getMatchingComponent(
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" ))).toXMLList(), 
                         [  '<covarianceMatrix index="1" type="relative">',
                            '<gridded dimension="2">',
                            '<axes>',
                            '<grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">',
                            '<values length="10">1e-5 8.5e5 2e6 3e6 4e6 5e6 6e6 7e6 1e7 2e7</values></grid>',
                            '<grid index="1" label="column_energy_bounds" unit="eV" style="link">',
                            '<link xlink:href="../grid[@index=\'2\']/values"/></grid>',
                            '<axis index="0" label="matrix_elements" unit=""/></axes>',
                            '<array shape="10,10" compression="diagonal">',
                            '<values length="10">0 0.0396 0.0891 0.0891 0.0891 0.0891 0.0891 0.0891 0.3782 0</values></array></gridded></covarianceMatrix>'] )

    def test_getMatchingComponent_0_stripped(self):
        x = FeCovariance[33]['eval'].getMatchingComponent(
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        x.removeExtraZeros()
        self.assertXMLListsEqual( x.toXMLList(),
                        [  '<covarianceMatrix index="1" type="relative">',
                            '<gridded dimension="2">',
                            '<axes>',
                            '<grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">',
                            '<values length="10">1e-5 8.5e5 2e6 3e6 4e6 5e6 6e6 7e6 1e7 2e7</values></grid>',
                            '<grid index="1" label="column_energy_bounds" unit="eV" style="link">',
                            '<link xlink:href="../grid[@index=\'2\']/values"/></grid>',
                            '<axis index="0" label="matrix_elements" unit=""/></axes>',
                            '<array shape="8,8" symmetry="lower">',
                            '<values length="36">0.0396 0 0.0891 0 0 0.0891 0 0 0 0.0891 0 0 0 0 0.0891 0 0 0 0 0 0.0891 0 0 0 0 0 0 0.0891 0 0 0 0 0 0 0 0.3782</values></array></gridded></covarianceMatrix>'] )

    def test_getMatchingComponent_1(self):
        '''Fails because there is no component of section 1 with such bounds'''
        self.assertRaises(  ValueError,
                            FeCovariance[1]['eval'].getMatchingComponent,
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )) )
    
    def test_getMatchingComponent_2(self):
        self.maxDiff=None
        self.assertXMLListsEqual( FeCovariance[1]['eval'].getMatchingComponent(
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "850636 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "850636 eV" ))).toXMLList(), 
                         [  '<covarianceMatrix index="1" type="relative">',
                            '<gridded dimension="2">',
                            '<axes>',
                            '<grid index="2" label="row_energy_bounds" unit="eV" style="boundaries">',
                            '<values length="14">1e-5 20 3e2 33830.5 76772 102950 152676 204928 238299 446595 483428 518904 553281 850636</values></grid>',
                            '<grid index="1" label="column_energy_bounds" unit="eV" style="link">',
                            '<link xlink:href="../grid[@index=\'2\']/values"/></grid>',
                            '<axis index="0" label="matrix_elements" unit=""/></axes>',
                            '<array shape="13,13" symmetry="lower">',
                            '<values length="91">1.6e-3 2.4e-3 3.6e-3 0 0 3.176954e-3 0 0 3.468456e-3 0.01514682 0 0 2.606664e-3 5.691679e-3 8.554986e-3 0 0 2.975089e-3 6.496139e-3 4.882072e-3 0.0111442 0 0 2.995769e-3 6.541293e-3 4.916007e-3 5.610834e-3 0.01129967 0 0 3.582893e-3 7.823286e-3 5.879469e-3 6.710471e-3 6.757115e-3 0.01616281 0 0 3.334617e-3 7.281172e-3 5.472051e-3 6.245469e-3 6.288881e-3 7.521405e-3 0.01400042 0 0 3.325916e-3 7.262173e-3 5.457774e-3 6.229173e-3 6.272472e-3 7.50178e-3 6.981945e-3 0.01392745 0 0 3.059465e-3 6.680374e-3 5.020531e-3 5.730131e-3 5.769962e-3 6.900785e-3 6.422595e-3 6.405837e-3 0.01178528 0 0 3.185348e-3 6.955242e-3 5.227104e-3 5.965901e-3 6.00737e-3 7.184721e-3 6.686857e-3 6.669409e-3 6.135098e-3 0.01277506 0 0 2.919549e-3 6.374865e-3 4.790931e-3 5.468079e-3 5.506088e-3 6.585196e-3 6.128875e-3 6.112884e-3 5.623158e-3 5.854526e-3 0.010732</values></array></gridded></covarianceMatrix>'] )
    
    def test_shrinkToBounds(self):
        self.maxDiff=None
        self.assertXMLListsEqual(
            FeCovariance[0]['eval'].shrinkToBounds((PQU.PQU( "862270 eV" ), PQU.PQU( "1.5e7 eV" )),).toXMLList(),
            ['<mixed>',
             '  <covarianceMatrix index="1" type="relative">',
             '    <axes>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="6"> 862270 1e6 2e6 5e6 1e7 1.5e7</axis>',
             '      <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>',
             '      <axis index="2" label="matrix_elements" unit=""/></axes>',
             '    <matrix rows="6" columns="6" form="diagonal" precision="6"> 9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix></mixed>'])
    
    def test_makeSafeBounds(self):
        self.assertEqual(
            [c.getRowBounds() for c in FeCovariance[0]['eval'].components],
            [(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "862270. eV" )), (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" ))])
        FeCovariance[0]['eval'].makeSafeBounds()
        self.assertEqual(
            [c.getRowBounds() for c in FeCovariance[0]['eval'].components],
            [(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "862270. eV" )), (PQU.PQU( "862270. eV" ), PQU.PQU( "1.5e7 eV" ))])

    def test_getRowBounds(self):
        self.assertEqual( FeCovariance[0]['eval'].getRowBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( FeCovariance[1]['eval'].getRowBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( FeCovariance[2]['eval'].getRowBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
    
    def test_getColumnBounds(self):
        self.assertEqual( FeCovariance[0]['eval'].getColumnBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( FeCovariance[1]['eval'].getColumnBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( FeCovariance[2]['eval'].getColumnBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
    
    def test_getUncertaintyVector(self): 
        self.assertEqual( 
            repr(FeCovariance[0]['eval'].getUncertaintyVector(
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
        
    def test_toCovarianceMatrix(self):
        self.assertXMLListsEqual( FeCovariance[36]['eval'].toCovarianceMatrix().toXMLList(),
            ['<covarianceMatrix type="relative">',
             '  <axes>',
             '    <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="32"> 1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899 862270 1e6 2e6 5e6 1e7 1.5e7 2e7</axis>',
             '    <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="32"> 1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899 862270 1e6 2e6 5e6 1e7 1.5e7 2e7</axis>',
             '    <axis index="2" label="matrix_elements" unit=""/></axes>',
             '  <matrix rows="31" columns="31" form="symmetric">',
             '    0.004521841',
             '    0.007561844 0.01376479',
             '    0.008361844 0.01456479 0.01576479',
             '    0.005961844 0.01216479 0.01216479 0.015341744',
             '    0.0 0.0 0.0 0.003176954 0.013101309',
             '    0.0 0.0 0.0 0.003176954 0.009338387 0.018478004',
             '    0.0 0.0 0.0 0.003468456 0.009629889 0.018769506 0.03044787',
             '    0.0 0.0 0.0 0.003468456 0.008425108 0.009623028 0.021301392 0.025049085',
             '    0.0 0.0 0.0 0.003468456 0.007281858 0.008203479 0.019881843 0.018955976 0.021007971',
             '    0.0 0.0 0.0 0.002606664 0.006420066 0.007341687 0.010426702 0.009500835 0.01155283 0.014416137',
             '    0.0 0.0 0.0 0.002606664 0.006074221 0.006912258 0.009997273 0.009155375 0.008356475 0.011219782 0.013401227',
             '    0.0 0.0 0.0 0.002606664 0.007052283 0.008126697 0.011211712 0.010132348 0.00910811 0.011971417 0.011661574 0.016520655',
             '    0.0 0.0 0.0 0.002606664 0.007977071 0.009274988 0.012360003 0.011056106 0.009818803 0.01268211 0.012307813 0.013366339 0.020179426',
             '    0.0 0.0 0.0 0.002975089 0.008345496 0.009643413 0.013164463 0.011860566 0.010623263 0.009009196 0.008634899 0.009693425 0.016506512 0.02276864',
             '    0.0 0.0 0.0 0.002975089 0.008882803 0.010310575 0.013831625 0.012397274 0.01103618 0.009422113 0.009010368 0.010174798 0.011275803 0.017537931 0.02521104',
             '    0.0 0.0 0.0 0.002995769 0.009381029 0.010924214 0.014469738 0.012919443 0.011448326 0.00982304 0.009378011 0.010636568 0.011826572 0.012521399 0.013212797 0.0277326',
             '    0.0 0.0 0.0 0.002995769 0.007243324 0.008269871 0.011815395 0.010784119 0.009805513 0.008180227 0.007884188 0.008721396 0.009513002 0.010207829 0.010667755 0.016765365 0.018571368',
             '    0.0 0.0 0.0 0.003582893 0.007830448 0.008856995 0.013097388 0.012066112 0.011087506 0.009143689 0.00884765 0.009684858 0.010476464 0.011307466 0.011767392 0.01222281 0.014028813 0.023434508',
             '    0.0 0.0 0.0 0.003334617 0.007582172 0.008608719 0.012555274 0.011523998 0.010545392 0.008736271 0.008440232 0.00927744 0.010069046 0.010842464 0.01130239 0.011754576 0.013560579 0.014793103 0.021272118',
             '    0.0 0.0 0.0 0.003334617 0.017252577 0.020616257 0.024562812 0.021183632 0.017977032 0.016167911 0.015197887 0.017941161 0.020535021 0.021308439 0.022815479 0.024198321 0.018202461 0.019434985 0.025914 0.09207485',
             '    0.0 0.0 0.0 0.003325916 0.017519686 0.020950026 0.024886283 0.021440143 0.018170003 0.016365604 0.015376347 0.018173984 0.020819244 0.021590643 0.023127553 0.024536822 0.018422142 0.01965145 0.019131615 0.046792765 0.09512696',
             '    0.0 0.0 0.0 0.003059465 0.017253235 0.020683575 0.024304484 0.020858344 0.017588204 0.015928361 0.014939104 0.017736741 0.020382001 0.021091601 0.022628511 0.024034312 0.017919632 0.019050455 0.018572265 0.046233415 0.087605347 0.09298479',
             '    0.0 0.0 0.0 0.003185348 0.017379118 0.020809458 0.024579352 0.021133212 0.017863072 0.016134934 0.015145677 0.017943314 0.020588574 0.021327371 0.022864281 0.02427172 0.01815704 0.019334391 0.018836527 0.046497677 0.087868919 0.087334608 0.09397457',
             '    0.0 0.0 0.0 0.0 0.01404215 0.01743584 0.01743584 0.01402651 0.0107913 0.0107913 0.009812619 0.01258037 0.01519737 0.01519737 0.01671786 0.01806924 0.01201989 0.01201989 0.01201989 0.03938554 0.04016605 0.04016605 0.04016605 0.07947396',
             '    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0',
             '    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 9.9e-05',
             '    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.005564',
             '    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.001584',
             '    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.000891',
             '    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.000396',
             '    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.000396</matrix></covarianceMatrix>'] )
    
    def test_toAbsolute(self): 
        self.assertXMLListsEqual(
            FeCovariance[0]['eval'].toAbsolute(
                rowData=FeEvaluation.getReaction('elastic').crossSection.toPointwise_withLinearXYs(1e-8,1e-8)).toXMLList(),
            ['<mixed>',
             '  <covarianceMatrix type="absolute">',
             '    <axes>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="25"> 1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899</axis>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="25"> 1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899</axis>',
             '      <axis index="2" label="matrix_elements" unit=""/></axes>',
             '    <matrix rows="24" columns="24" form="symmetric">',
             '      0.004521841',
             '      0.007561844 0.01376479',
             '      0.008361844 0.01456479 0.01576479',
             '      0.005961844 0.01216479 0.01216479 0.015341744',
             '      0.0 0.0 0.0 0.003176954 0.013101309',
             '      0.0 0.0 0.0 0.003176954 0.009338387 0.018478004',
             '      0.0 0.0 0.0 0.003468456 0.009629889 0.018769506 0.03044787',
             '      0.0 0.0 0.0 0.003468456 0.008425108 0.009623028 0.021301392 0.025049085',
             '      0.0 0.0 0.0 0.003468456 0.007281858 0.008203479 0.019881843 0.018955976 0.021007971',
             '      0.0 0.0 0.0 0.002606664 0.006420066 0.007341687 0.010426702 0.009500835 0.01155283 0.014416137',
             '      0.0 0.0 0.0 0.002606664 0.006074221 0.006912258 0.009997273 0.009155375 0.008356475 0.011219782 0.013401227',
             '      0.0 0.0 0.0 0.002606664 0.007052283 0.008126697 0.011211712 0.010132348 0.00910811 0.011971417 0.011661574 0.016520655',
             '      0.0 0.0 0.0 0.002606664 0.007977071 0.009274988 0.012360003 0.011056106 0.009818803 0.01268211 0.012307813 0.013366339 0.020179426',
             '      0.0 0.0 0.0 0.002975089 0.008345496 0.009643413 0.013164463 0.011860566 0.010623263 0.009009196 0.008634899 0.009693425 0.016506512 0.02276864',
             '      0.0 0.0 0.0 0.002975089 0.008882803 0.010310575 0.013831625 0.012397274 0.01103618 0.009422113 0.009010368 0.010174798 0.011275803 0.017537931 0.02521104',
             '      0.0 0.0 0.0 0.002995769 0.009381029 0.010924214 0.014469738 0.012919443 0.011448326 0.00982304 0.009378011 0.010636568 0.011826572 0.012521399 0.013212797 0.0277326',
             '      0.0 0.0 0.0 0.002995769 0.007243324 0.008269871 0.011815395 0.010784119 0.009805513 0.008180227 0.007884188 0.008721396 0.009513002 0.010207829 0.010667755 0.016765365 0.018571368',
             '      0.0 0.0 0.0 0.003582893 0.007830448 0.008856995 0.013097388 0.012066112 0.011087506 0.009143689 0.00884765 0.009684858 0.010476464 0.011307466 0.011767392 0.01222281 0.014028813 0.023434508',
             '      0.0 0.0 0.0 0.003334617 0.007582172 0.008608719 0.012555274 0.011523998 0.010545392 0.008736271 0.008440232 0.00927744 0.010069046 0.010842464 0.01130239 0.011754576 0.013560579 0.014793103 0.021272118',
             '      0.0 0.0 0.0 0.003334617 0.017252577 0.020616257 0.024562812 0.021183632 0.017977032 0.016167911 0.015197887 0.017941161 0.020535021 0.021308439 0.022815479 0.024198321 0.018202461 0.019434985 0.025914 0.09207485',
             '      0.0 0.0 0.0 0.003325916 0.017519686 0.020950026 0.024886283 0.021440143 0.018170003 0.016365604 0.015376347 0.018173984 0.020819244 0.021590643 0.023127553 0.024536822 0.018422142 0.01965145 0.019131615 0.046792765 0.09512696',
             '      0.0 0.0 0.0 0.003059465 0.017253235 0.020683575 0.024304484 0.020858344 0.017588204 0.015928361 0.014939104 0.017736741 0.020382001 0.021091601 0.022628511 0.024034312 0.017919632 0.019050455 0.018572265 0.046233415 0.087605347 0.09298479',
             '      0.0 0.0 0.0 0.003185348 0.017379118 0.020809458 0.024579352 0.021133212 0.017863072 0.016134934 0.015145677 0.017943314 0.020588574 0.021327371 0.022864281 0.02427172 0.01815704 0.019334391 0.018836527 0.046497677 0.087868919 0.087334608 0.09397457',
             '      0.0 0.0 0.0 0.0 0.01404215 0.01743584 0.01743584 0.01402651 0.0107913 0.0107913 0.009812619 0.01258037 0.01519737 0.01519737 0.01671786 0.01806924 0.01201989 0.01201989 0.01201989 0.03938554 0.04016605 0.04016605 0.04016605 0.07947396</matrix></covarianceMatrix>',
             '  <covarianceMatrix index="1" type="absolute">',
             '    <axes>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="8"> 1e-5 862270 1e6 2e6 5e6 1e7 1.5e7 2e7</axis>',
             '      <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>',
             '      <axis index="2" label="matrix_elements" unit=""/></axes>',
             '    <matrix rows="7" columns="7" form="diagonal" precision="6"> 0.000000e+00  9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix></mixed>'] )

    def test_toRelative(self): 
        self.assertXMLListsEqual( FeCovariance[0]['eval'].toRelative().toXMLList(),
            ['<mixed>',
             '  <covarianceMatrix type="relative">',
             '    <axes>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="25"> 1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899</axis>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="25"> 1e-5 0.5 20 3e2 1140.04 1751 33830.5 35479.5 53647.5 76772 77931.5 89369 94561 102950 104508 152676 203065 204928 238299 261280 446595 483428 518904 553281 794899</axis>',
             '      <axis index="2" label="matrix_elements" unit=""/></axes>',
             '    <matrix rows="24" columns="24" form="symmetric">',
             '      0.004521841',
             '      0.007561844 0.01376479',
             '      0.008361844 0.01456479 0.01576479',
             '      0.005961844 0.01216479 0.01216479 0.015341744',
             '      0.0 0.0 0.0 0.003176954 0.013101309',
             '      0.0 0.0 0.0 0.003176954 0.009338387 0.018478004',
             '      0.0 0.0 0.0 0.003468456 0.009629889 0.018769506 0.03044787',
             '      0.0 0.0 0.0 0.003468456 0.008425108 0.009623028 0.021301392 0.025049085',
             '      0.0 0.0 0.0 0.003468456 0.007281858 0.008203479 0.019881843 0.018955976 0.021007971',
             '      0.0 0.0 0.0 0.002606664 0.006420066 0.007341687 0.010426702 0.009500835 0.01155283 0.014416137',
             '      0.0 0.0 0.0 0.002606664 0.006074221 0.006912258 0.009997273 0.009155375 0.008356475 0.011219782 0.013401227',
             '      0.0 0.0 0.0 0.002606664 0.007052283 0.008126697 0.011211712 0.010132348 0.00910811 0.011971417 0.011661574 0.016520655',
             '      0.0 0.0 0.0 0.002606664 0.007977071 0.009274988 0.012360003 0.011056106 0.009818803 0.01268211 0.012307813 0.013366339 0.020179426',
             '      0.0 0.0 0.0 0.002975089 0.008345496 0.009643413 0.013164463 0.011860566 0.010623263 0.009009196 0.008634899 0.009693425 0.016506512 0.02276864',
             '      0.0 0.0 0.0 0.002975089 0.008882803 0.010310575 0.013831625 0.012397274 0.01103618 0.009422113 0.009010368 0.010174798 0.011275803 0.017537931 0.02521104',
             '      0.0 0.0 0.0 0.002995769 0.009381029 0.010924214 0.014469738 0.012919443 0.011448326 0.00982304 0.009378011 0.010636568 0.011826572 0.012521399 0.013212797 0.0277326',
             '      0.0 0.0 0.0 0.002995769 0.007243324 0.008269871 0.011815395 0.010784119 0.009805513 0.008180227 0.007884188 0.008721396 0.009513002 0.010207829 0.010667755 0.016765365 0.018571368',
             '      0.0 0.0 0.0 0.003582893 0.007830448 0.008856995 0.013097388 0.012066112 0.011087506 0.009143689 0.00884765 0.009684858 0.010476464 0.011307466 0.011767392 0.01222281 0.014028813 0.023434508',
             '      0.0 0.0 0.0 0.003334617 0.007582172 0.008608719 0.012555274 0.011523998 0.010545392 0.008736271 0.008440232 0.00927744 0.010069046 0.010842464 0.01130239 0.011754576 0.013560579 0.014793103 0.021272118',
             '      0.0 0.0 0.0 0.003334617 0.017252577 0.020616257 0.024562812 0.021183632 0.017977032 0.016167911 0.015197887 0.017941161 0.020535021 0.021308439 0.022815479 0.024198321 0.018202461 0.019434985 0.025914 0.09207485',
             '      0.0 0.0 0.0 0.003325916 0.017519686 0.020950026 0.024886283 0.021440143 0.018170003 0.016365604 0.015376347 0.018173984 0.020819244 0.021590643 0.023127553 0.024536822 0.018422142 0.01965145 0.019131615 0.046792765 0.09512696',
             '      0.0 0.0 0.0 0.003059465 0.017253235 0.020683575 0.024304484 0.020858344 0.017588204 0.015928361 0.014939104 0.017736741 0.020382001 0.021091601 0.022628511 0.024034312 0.017919632 0.019050455 0.018572265 0.046233415 0.087605347 0.09298479',
             '      0.0 0.0 0.0 0.003185348 0.017379118 0.020809458 0.024579352 0.021133212 0.017863072 0.016134934 0.015145677 0.017943314 0.020588574 0.021327371 0.022864281 0.02427172 0.01815704 0.019334391 0.018836527 0.046497677 0.087868919 0.087334608 0.09397457',
             '      0.0 0.0 0.0 0.0 0.01404215 0.01743584 0.01743584 0.01402651 0.0107913 0.0107913 0.009812619 0.01258037 0.01519737 0.01519737 0.01671786 0.01806924 0.01201989 0.01201989 0.01201989 0.03938554 0.04016605 0.04016605 0.04016605 0.07947396</matrix></covarianceMatrix>',
             '  <covarianceMatrix index="1" type="relative">',
             '    <axes>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="8"> 1e-5 862270 1e6 2e6 5e6 1e7 1.5e7 2e7</axis>',
             '      <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>',
             '      <axis index="2" label="matrix_elements" unit=""/></axes>',
             '    <matrix rows="7" columns="7" form="diagonal" precision="6"> 0.000000e+00  9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix></mixed>'] )
    

if __name__=="__main__":
    unittest.main()
