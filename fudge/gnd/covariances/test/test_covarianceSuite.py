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

import unittest, cStringIO
from fudge.gnd.covariances import *
from fudge.gnd.reactionSuite import reactionSuite

Fe56TestData = '''<?xml version="1.0" encoding="UTF-8"?>
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

class Test_covarianceSuite( unittest.TestCase ):
    
    def setUp(self): 
        self.covSuite = readXML(cStringIO.StringIO(Fe56TestData)) #,reactionSuite=reactionSuite('n','Fe56'))
    
    def test_readXML(self):
        '''Also tests __init__ & parseXMLNode'''
        self.assertIsInstance( self.covSuite, covarianceSuite )
        
    def test__getitem__(self):
        self.assertEqual( '\n'.join(self.covSuite[2].toXMLList()), '''<section label="33" id="Fe57 + gamma" nativeData="mixed">
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
      <matrix rows="9" columns="9" form="diagonal" precision="6"> 0.000000e+00  3.960000e-02  8.910000e-02  8.910000e-02  8.910000e-02  8.910000e-02  8.910000e-02  8.910000e-02  3.782000e-01 </matrix></covarianceMatrix></mixed></section>'''    )

    def test__len__(self):
        self.assertEqual( len(self.covSuite),3)

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
