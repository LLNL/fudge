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

Fe56CovTestData = '''<?xml version="1.0" encoding="UTF-8"?>
<covarianceSuite projectile="n" target="Fe56" format="gnd version 1.0" xmlns:xlink="http://www.w3.org/1999/xlink">
  <styles>
    <style name="evaluated" version="7.1.4" library="ENDF/B"></style></styles>
  <section label="0" id="total" nativeData="mixed">
    <rowData xlink:type="simple" xlink:href="/reactionSuite/summedReaction[@label='37']/crossSection" ENDF_MFMT="33,1"/>
    <mixed>
      <sum index="0" lowerBound="1e-5 eV" upperBound="862270 eV">
        <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->
        <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label='1']" coefficient="1.0" ENDF_MFMT="33,2"/>
        <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label='33']" coefficient="1.0" ENDF_MFMT="33,102"/></sum>
      <covarianceMatrix index="1" type="relative">
        <axes>
          <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="8"> 1e-5 862270 1e6 2e6 5e6 1e7 1.5e7 2e7</axis>
          <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>
          <axis index="2" label="matrix_elements" unit=""/></axes>
        <matrix rows="7" columns="7" form="diagonal" precision="6"> 0.000000e+00  9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix></mixed></section>
  <section label="1" id="n + Fe56" nativeData="mixed">
    <rowData xlink:type="simple" xlink:href="/reactionSuite/reaction[@label='0']/crossSection" ENDF_MFMT="33,2"/>
    <mixed>
      <sum index="0" lowerBound="862270 eV" upperBound="2e7 eV">
        <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->
        <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label='0']" coefficient="1.0" ENDF_MFMT="33,1"/>
        <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label='33']" coefficient="-1.0" ENDF_MFMT="33,102"/></sum>
      <covarianceMatrix index="1" type="relative">
        <axes>
          <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="14"> 1e-5 20 3e2 33830.5 76772 102950 152676 204928 238299 446595 483428 518904 553281 850636</axis>
          <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>
          <axis index="2" label="matrix_elements" unit=""/></axes>
        <matrix rows="13" columns="13" form="symmetric" precision="6">
           1.600000e-03
           2.400000e-03  3.600000e-03
           0.000000e+00  0.000000e+00  3.176954e-03
           0.000000e+00  0.000000e+00  3.468456e-03  1.514682e-02
           0.000000e+00  0.000000e+00  2.606664e-03  5.691679e-03  8.554986e-03
           0.000000e+00  0.000000e+00  2.975089e-03  6.496139e-03  4.882072e-03  1.114420e-02
           0.000000e+00  0.000000e+00  2.995769e-03  6.541293e-03  4.916007e-03  5.610834e-03  1.129967e-02
           0.000000e+00  0.000000e+00  3.582893e-03  7.823286e-03  5.879469e-03  6.710471e-03  6.757115e-03  1.616281e-02
           0.000000e+00  0.000000e+00  3.334617e-03  7.281172e-03  5.472051e-03  6.245469e-03  6.288881e-03  7.521405e-03  1.400042e-02
           0.000000e+00  0.000000e+00  3.325916e-03  7.262173e-03  5.457774e-03  6.229173e-03  6.272472e-03  7.501780e-03  6.981945e-03  1.392745e-02
           0.000000e+00  0.000000e+00  3.059465e-03  6.680374e-03  5.020531e-03  5.730131e-03  5.769962e-03  6.900785e-03  6.422595e-03  6.405837e-03  1.178528e-02
           0.000000e+00  0.000000e+00  3.185348e-03  6.955242e-03  5.227104e-03  5.965901e-03  6.007370e-03  7.184721e-03  6.686857e-03  6.669409e-03  6.135098e-03  1.277506e-02
           0.000000e+00  0.000000e+00  2.919549e-03  6.374865e-03  4.790931e-03  5.468079e-03  5.506088e-03  6.585196e-03  6.128875e-03  6.112884e-03  5.623158e-03  5.854526e-03  1.073200e-02</matrix></covarianceMatrix></mixed></section>
  <section label="33" id="Fe57 + gamma" nativeData="mixed">
    <rowData xlink:type="simple" xlink:href="/reactionSuite/reaction[@label='29']/crossSection" ENDF_MFMT="33,102"/>
    <mixed>
      <covarianceMatrix index="0" type="relative">
        <axes>
          <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="17"> 1e-5 0.5 1140.04 1751 35479.5 53647.5 77931.5 89369 94561 104508 152676 203065 261280 446595 553281 794899 850636</axis>
          <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>
          <axis index="2" label="matrix_elements" unit=""/></axes>
        <matrix rows="16" columns="16" form="symmetric" precision="6">
           2.921841e-03
           5.961844e-03  1.216479e-02
           0.000000e+00  0.000000e+00  9.924355e-03
           0.000000e+00  0.000000e+00  6.161433e-03  1.530105e-02
           0.000000e+00  0.000000e+00  4.956652e-03  6.154572e-03  9.902265e-03
           0.000000e+00  0.000000e+00  3.813402e-03  4.735023e-03  3.809156e-03  5.861151e-03
           0.000000e+00  0.000000e+00  3.467557e-03  4.305594e-03  3.463696e-03  2.664796e-03  4.846241e-03
           0.000000e+00  0.000000e+00  4.445619e-03  5.520033e-03  4.440669e-03  3.416431e-03  3.106588e-03  7.965669e-03
           0.000000e+00  0.000000e+00  5.370407e-03  6.668324e-03  5.364427e-03  4.127124e-03  3.752827e-03  4.811353e-03  1.162444e-02
           0.000000e+00  0.000000e+00  5.907714e-03  7.335486e-03  5.901135e-03  4.540041e-03  4.128296e-03  5.292726e-03  6.393731e-03  1.406684e-02
           0.000000e+00  0.000000e+00  6.385260e-03  7.928445e-03  6.378150e-03  4.907033e-03  4.462004e-03  5.720561e-03  6.910565e-03  7.601963e-03  1.643293e-02
           0.000000e+00  0.000000e+00  4.247555e-03  5.274102e-03  4.242826e-03  3.264220e-03  2.968181e-03  3.805389e-03  4.596995e-03  5.056921e-03  5.465695e-03  7.271698e-03
           0.000000e+00  0.000000e+00  1.391796e-02  1.728164e-02  1.390246e-02  1.069586e-02  9.725836e-03  1.246911e-02  1.506297e-02  1.657001e-02  1.790944e-02  1.191358e-02  7.807443e-02
           0.000000e+00  0.000000e+00  1.419377e-02  1.762411e-02  1.417797e-02  1.090783e-02  9.918573e-03  1.271621e-02  1.536147e-02  1.689838e-02  1.826435e-02  1.214967e-02  3.981082e-02  8.119951e-02
           0.000000e+00  0.000000e+00  1.404215e-02  1.743584e-02  1.402651e-02  1.079130e-02  9.812619e-03  1.258037e-02  1.519737e-02  1.671786e-02  1.806924e-02  1.201989e-02  3.938554e-02  4.016605e-02  7.947396e-02
           0.000000e+00  0.000000e+00  1.433920e-02  1.780468e-02  1.432323e-02  1.101958e-02  1.002019e-02  1.284650e-02  1.551886e-02  1.707151e-02  1.845148e-02  1.227415e-02  4.021870e-02  4.101572e-02  4.057757e-02  8.287190e-02</matrix></covarianceMatrix>
      <covarianceMatrix index="1" type="relative">
        <axes>
          <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="10"> 1e-5 8.5e5 2e6 3e6 4e6 5e6 6e6 7e6 1e7 2e7</axis>
          <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>
          <axis index="2" label="matrix_elements" unit=""/></axes>
        <matrix rows="9" columns="9" form="diagonal" precision="6"> 0.000000e+00  3.960000e-02  8.910000e-02  8.910000e-02  8.910000e-02  8.910000e-02  8.910000e-02  8.910000e-02  3.782000e-01 </matrix></covarianceMatrix></mixed></section>
</covarianceSuite>
'''
Fe56RxnTestData = '''<?xml version="1.0" encoding="UTF-8"?>
<reactionSuite projectile="n" target="Fe56" format="gnd version 1.2" temperature="0 K" xmlns:xlink="http://www.w3.org/1999/xlink">
  <styles>
    <style name="evaluated" version="7.1.4" library="ENDF/B"></style></styles>
  <documentations>
    <documentation name="endfDoc"><![CDATA[
 26-Fe- 56 LANL,ORNL  EVAL-SEP96 M.B.Chadwick,P.G.Young,C.Y.Fu
]]></documentation></documentations>
  <particles>
    <particle name="gamma" genre="photon" mass="0 amu" transportable="true"/>
    <particle name="n" genre="nucleus" mass="1.00866491574 amu" transportable="true"/>
    <particle name="Fe56" genre="nucleus" mass="55.934504237446 amu">
      <level name="Fe56_e0" label="0" energy="0 eV" spin="0"/></particle>
    <particle name="Fe57" genre="nucleus" mass="56.935393969 amu"/></particles>
  <reaction label="0" outputChannel="n + Fe56" date="1996-09-01" ENDF_MT="2">
    <crossSection nativeData="linear">
      <linear xData="XYs" length="2" accuracy="0.001">
        <axes>
          <axis index="0" label="energy_in" unit="eV" interpolation="lin,lin" frame="lab"/>
          <axis index="1" label="crossSection" unit="b" frame="lab"/></axes>
        <data> 1e-5 1.0 2e7 1.0 </data></linear></crossSection>
    <outputChannel genre="twoBody" Q="0.0 eV">
      <product name="n" label="n" multiplicity="1">
        <distributions nativeData="none"></distributions></product>
      <product name="Fe56" label="Fe56" multiplicity="1">
        <distributions nativeData="none"></distributions></product></outputChannel></reaction>
  <reaction label="29" outputChannel="Fe57 + gamma" date="1996-09-01" ENDF_MT="102">
    <crossSection nativeData="linear">
      <linear xData="XYs" length="2" accuracy="0.001">
        <axes>
          <axis index="0" label="energy_in" unit="eV" interpolation="lin,lin" frame="lab"/>
          <axis index="1" label="crossSection" unit="b" frame="lab"/></axes>
        <data> 1e-5 1.0 2e7 1.0 </data></linear></crossSection>
    <outputChannel genre="twoBody" Q="7646090 eV">
      <product name="gamma" label="gamma" multiplicity="1">
        <distributions nativeData="none"></distributions></product>
      <product name="Fe57" label="Fe57" multiplicity="1">
        <distributions nativeData="none"></distributions></product></outputChannel></reaction>
  <summedReaction label="37" name="total" Q="0 eV" date="1996-09-01" ENDF_MT="1">
    <crossSection nativeData="linear">
      <linear xData="XYs" length="2" accuracy="0.001">
        <axes>
          <axis index="0" label="energy_in" unit="eV" interpolation="lin,lin" frame="lab"/>
          <axis index="1" label="crossSection" unit="b" frame="lab"/></axes>
        <data> 1e-5 2.0 2e7 2.0 </data></linear></crossSection></summedReaction>
</reactionSuite>'''

import unittest, cStringIO #, sys
from pqu import PQU
from fudge.gnd.covariances import readXML as CovReadXML
from fudge.gnd.reactionSuite import readXML as RxnReadXML

#sys.setrecursionlimit(75)

class Test_mixed( unittest.TestCase ): 
    
    def setUp(self): 
        self.rxnSuite = RxnReadXML(cStringIO.StringIO(Fe56RxnTestData))
        self.covSuite = CovReadXML(cStringIO.StringIO(Fe56CovTestData),reactionSuite=self.rxnSuite)

    def test__getitem__(self):
        self.assertEqual( 
            self.covSuite[1].getNativeData()[0].toXMLList(), 
            ['<sum index="0" lowerBound="862270 eV" upperBound="2e7 eV">',
             '  <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->',
             '  <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label=\'0\']" coefficient="1.0" ENDF_MFMT="33,1"/>',
             '  <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label=\'33\']" coefficient="-1.0" ENDF_MFMT="33,102"/></sum>'])

    def test__len__(self): 
        self.assertEqual( len(self.covSuite[1].getNativeData()), 2 )
        
    def test_toXMLList(self):
        self.assertEqual( 
            self.covSuite[0].getNativeData().toXMLList(), 
            ['<mixed>',
             '  <sum index="0" lowerBound="1e-5 eV" upperBound="862270 eV">',
             '    <!-- The matrix for this reaction equals the weighted sum of the following matrices: -->',
             '    <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label=\'1\']" coefficient="1.0" ENDF_MFMT="33,2"/>',
             '    <summand xlink:type="simple" xlink:href="/covarianceSuite/section[@label=\'33\']" coefficient="1.0" ENDF_MFMT="33,102"/></sum>',
             '  <covarianceMatrix index="1" type="relative">',
             '    <axes>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="8"> 1e-5 862270 1e6 2e6 5e6 1e7 1.5e7 2e7</axis>',
             '      <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>',
             '      <axis index="2" label="matrix_elements" unit=""/></axes>',
             '    <matrix rows="7" columns="7" form="diagonal" precision="6"> 0.000000e+00  9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix></mixed>'])

    def test_check(self): pass
    
    def test_fix(self): pass
    
    def test_plot(self):pass
    
    def test_addComponent(self): pass
    
    def test_getMatchingComponent_0(self):
        self.assertEqual( self.covSuite[0].getNativeData().getMatchingComponent(
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" ))).toXMLList(), 
                         ['<covarianceMatrix index="1" type="relative">', 
                          '  <axes>', 
                          '    <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="8"> 1e-5 862270 1e6 2e6 5e6 1e7 1.5e7 2e7</axis>', 
                          '    <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>', 
                          '    <axis index="2" label="matrix_elements" unit=""/></axes>', 
                          '  <matrix rows="7" columns="7" form="diagonal" precision="6"> 0.000000e+00  9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix>'] )

    def test_getMatchingComponent_0_stripped(self):
        x = self.covSuite[0].getNativeData().getMatchingComponent(
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        x.removeExtraZeros()
        self.assertEqual( x.toXMLList(), 
                         ['<covarianceMatrix index="1" type="relative">', 
                          '  <axes>', 
                          '    <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="6"> 862270 1e6 2e6 5e6 1e7 1.5e7</axis>', 
                          '    <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>', 
                          '    <axis index="2" label="matrix_elements" unit=""/></axes>', 
                          '  <matrix rows="6" columns="6" form="diagonal" precision="6"> 9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix>'] )

    @unittest.expectedFailure
    def test_getMatchingComponent_1(self):
        '''Fails because there is no component of section 1 with such bounds'''
        self.assertFails( self.covSuite[1].getNativeData().getMatchingComponent(
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" ))).toXMLList(), [] )
    
    def test_getMatchingComponent_2(self):
        self.assertEqual( self.covSuite[1].getNativeData().getMatchingComponent(
                            rowBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "850636 eV" )),
                            columnBounds = (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "850636 eV" ))).toXMLList(), 
                         ['<covarianceMatrix index="1" type="relative">', 
                          '  <axes>', 
                          '    <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="14"> 1e-5 20 3e2 33830.5 76772 102950 152676 204928 238299 446595 483428 518904 553281 850636</axis>', 
                          '    <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>', 
                          '    <axis index="2" label="matrix_elements" unit=""/></axes>', 
                          '  <matrix rows="13" columns="13" form="symmetric" precision="6">', 
                          '     1.600000e-03', 
                          '     2.400000e-03  3.600000e-03', 
                          '     0.000000e+00  0.000000e+00  3.176954e-03', 
                          '     0.000000e+00  0.000000e+00  3.468456e-03  1.514682e-02', 
                          '     0.000000e+00  0.000000e+00  2.606664e-03  5.691679e-03  8.554986e-03', 
                          '     0.000000e+00  0.000000e+00  2.975089e-03  6.496139e-03  4.882072e-03  1.114420e-02', 
                          '     0.000000e+00  0.000000e+00  2.995769e-03  6.541293e-03  4.916007e-03  5.610834e-03  1.129967e-02', 
                          '     0.000000e+00  0.000000e+00  3.582893e-03  7.823286e-03  5.879469e-03  6.710471e-03  6.757115e-03  1.616281e-02', 
                          '     0.000000e+00  0.000000e+00  3.334617e-03  7.281172e-03  5.472051e-03  6.245469e-03  6.288881e-03  7.521405e-03  1.400042e-02', 
                          '     0.000000e+00  0.000000e+00  3.325916e-03  7.262173e-03  5.457774e-03  6.229173e-03  6.272472e-03  7.501780e-03  6.981945e-03  1.392745e-02', 
                          '     0.000000e+00  0.000000e+00  3.059465e-03  6.680374e-03  5.020531e-03  5.730131e-03  5.769962e-03  6.900785e-03  6.422595e-03  6.405837e-03  1.178528e-02', 
                          '     0.000000e+00  0.000000e+00  3.185348e-03  6.955242e-03  5.227104e-03  5.965901e-03  6.007370e-03  7.184721e-03  6.686857e-03  6.669409e-03  6.135098e-03  1.277506e-02', 
                          '     0.000000e+00  0.000000e+00  2.919549e-03  6.374865e-03  4.790931e-03  5.468079e-03  5.506088e-03  6.585196e-03  6.128875e-03  6.112884e-03  5.623158e-03  5.854526e-03  1.073200e-02</matrix></covarianceMatrix>'] )
    
    def test_shrinkToBounds(self):
        self.assertEqual( 
            self.covSuite[0].getNativeData().shrinkToBounds((PQU.PQU( "862270 eV" ), PQU.PQU( "1.5e7 eV" )),).toXMLList(),
            ['<mixed>',
             '  <covarianceMatrix index="1" type="relative">',
             '    <axes>',
             '      <axis index="0" label="row_energy_bounds" unit="eV" interpolation="lin,flat" length="6"> 862270 1e6 2e6 5e6 1e7 1.5e7</axis>',
             '      <axis index="1" label="column_energy_bounds" unit="eV" interpolation="lin,flat" mirror_row_energy_bounds="true"/>',
             '      <axis index="2" label="matrix_elements" unit=""/></axes>',
             '    <matrix rows="6" columns="6" form="diagonal" precision="6"> 9.900000e-05  5.564000e-03  1.584000e-03  8.910000e-04  3.960000e-04  3.960000e-04 </matrix></covarianceMatrix></mixed>'])
    
    def test_makeSafeBounds(self):
        self.assertEqual(
            [c.getRowBounds() for c in self.covSuite[0].getNativeData().components],
            [(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "862270. eV" )), (PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" ))])
        self.covSuite[0].getNativeData().makeSafeBounds()
        self.assertEqual(
            [c.getRowBounds() for c in self.covSuite[0].getNativeData().components],
            [(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "862270. eV" )), (PQU.PQU( "862270. eV" ), PQU.PQU( "1.5e7 eV" ))])

    def test_getRowBounds(self):
        self.assertEqual( self.covSuite[0].getNativeData().getRowBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( self.covSuite[1].getNativeData().getRowBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( self.covSuite[2].getNativeData().getRowBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
    
    def test_getColumnBounds(self):
        self.assertEqual( self.covSuite[0].getNativeData().getColumnBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( self.covSuite[1].getNativeData().getColumnBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
        self.assertEqual( self.covSuite[2].getNativeData().getColumnBounds(),(PQU.PQU( "1.e-5 eV" ), PQU.PQU( "2.e7 eV" )))
    
    def test_getUncertaintyVector(self): 
        self.assertEqual( 
            repr(self.covSuite[0].getNativeData().getUncertaintyVector(
                theData=self.rxnSuite.getReaction('elastic').crossSection.getNativeData())),
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
        self.assertEqual( self.covSuite[0].getNativeData().toCovarianceMatrix().toXMLList(), 
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
        self.assertEqual( 
            self.covSuite[0].getNativeData().toAbsolute(
                rowData=self.rxnSuite.getReaction('elastic').crossSection.getNativeData().toPointwise_withLinearXYs(1e-8,1e-8)).toXMLList(), 
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
        self.assertEqual( self.covSuite[0].getNativeData().toRelative().toXMLList(), 
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
