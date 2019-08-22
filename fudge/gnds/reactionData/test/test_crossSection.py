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
test fudge/gnds/reactionData/
dbrown, 12/5/2012
"""

import unittest
from fudge.gnds.reactionData.crossSection import *
defaultAccuracy = 0.001
from xData import XYs

class TestPointwise( unittest.TestCase ):

    def test_basic_stuff(self):
        ptwise_const=XYs1d(axes=XYs1d.defaultAxes(labelsUnits={
            XYs.yAxisIndex : ( 'crossSection', 'b' ),
            XYs.xAxisIndex : ( 'energy_in', 'eV' ) }), data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.assertEqual(ptwise_const.domainUnit, 'eV' )
        self.assertEqual(ptwise_const.toXMLList(), [
            '<XYs1d>',
            '  <axes>',
            '    <axis index="1" label="energy_in" unit="eV"/>',
            '    <axis index="0" label="crossSection" unit="b"/></axes>',
            '  <values>',
            '    1.00000000e-05 1.00000000e+00 2.00000000e+07 1.00000000e+00</values></XYs1d>'] )
            
    def test_group(self):
        ptwise=XYs1d(axes=XYs1d.defaultAxes(), data=[
            [1e-5,1.0], 
            [1.0,1.0], 
            [10.0,2.0], 
            [100.0,3.0], 
            [1000.0,1.0], 
            [20.0e6,1.0] ] )
#        ptwise.plot()
        ptwise.group( [1e-5,1.0,1000.0,20.e6] )#.plot()

    def test_evaluate(self):
        ptwise_const=XYs1d(axes=XYs1d.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( 'crossSection', 'b' ), \
                XYs.xAxisIndex : ( 'energy_in', 'eV' ) }), data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.assertEqual(ptwise_const.evaluate( 10.0e6 ), 1.0 )


class TestCrossSection( unittest.TestCase ): 

    def setUp(self):
    
        # A constant (kind of like (n,el))
        self.xs_const = component()
        self.ptwise_const=XYs1d(label='eval',axes=XYs1d.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( 'crossSection', 'b' ), \
                XYs.xAxisIndex : ( 'energy_in', 'eV' ) }), \
            data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.xs_const.add(self.ptwise_const)

        # A triangle (kind of like (n,2n), if you squint real hard)
        self.xs_triangle = component()
        self.ptwise_triangle=XYs1d(label='eval',axes=XYs1d.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( 'crossSection', 'b' ), \
                XYs.xAxisIndex : ( 'energy_in', 'eV' ) }), \
            data=[ [1e6,0.0], [10.0e6,1.0], [20.0e6,0.0] ] )
        self.xs_triangle.add(self.ptwise_triangle)

    def test_basic_stuff(self):
        self.assertTrue( isinstance( self.xs_const, component ) )
        self.assertEqual( self.xs_const.attribute, None )
        self.assertEqual( self.xs_const.moniker, 'crossSection' )
        self.assertEqual( self.xs_const.domainUnit, 'eV' )
        
    def test_xlink(self):
        self.assertEqual( self.xs_const.toRelativeXLink(), '' )
        
    def test_parent(self):
        '''
        Tests of traversing up GNDS hierarchy
        '''
        self.assertEqual( self.xs_const.ancestor, None )
        self.assertEqual( self.xs_const.getRootAncestor(), self.xs_const )
        self.assertEqual( self.xs_const.ancestor, None )
        self.xs_const.setAncestor( 'fred' )
        self.assertEqual( self.xs_const.ancestor, 'fred' )

    def test_ancestry(self):
        '''
        FIXME: findAttributeInAncestry
        FIXME: findClassInAncestry
        '''
        pass

    def test_forms(self):
        self.assertItemsEqual(
            [x.moniker for x in self.xs_const.allowedClasses],
            ['Ys1d','XYs1d','regions1d','gridded1d','resonancesWithBackground','reference'] )
        self.xs_const.remove('gridded')

    def test_evaluated(self):
        self.assertEqual( self.xs_const.evaluated.label, 'eval' )
        self.assertEqual( self.xs_const.evaluated.moniker, 'XYs1d' )
        self.assertEqual( self.xs_const.evaluated.evaluate(1.4e5), 1.0 )
        
    def test_domains(self):
        for form in self.xs_const:
            self.assertEqual(form.domainMax, 20e6 )
            self.assertEqual(form.domainMax, 20.0e6 )
            self.assertEqual(form.domainMin, 1.0e-5 )
            self.assertEqual(form.domainUnitConversionFactor('MeV'), 1e-6 )
            self.assertEqual(form.domainMin, 1e-05 )
            self.assertEqual(form.domainMax, 20000000.0 )
            self.assertEqual(form.domainUnit, 'eV' )
        self.assertEqual( self.xs_const.domainMax, 20.0e6 )
        self.assertEqual( self.xs_const.domainMin, 1.0e-5 )
        self.assertEqual( self.xs_const.domainUnitConversionFactor('MeV'), 1e-6 )
        self.assertEqual( self.xs_const.domainMin, 1e-05 )
        self.assertEqual( self.xs_const.domainMax, 20000000.0 )
        self.assertEqual( self.xs_const.domainUnit, 'eV' )

    def test_check(self):
        from pqu import PQU
        info = {'Q':1.0,'kinematicFactor':1.0,'dThreshold':PQU.PQU('1e-4 b'),
                'crossSectionEnergyMin':PQU.PQU('1.0e-5 eV'), 'crossSectionEnergyMax':PQU.PQU('20.0 MeV'),
                'CoulombChannel':False}
        self.assertEqual( self.xs_const.check(info), [] )

    def test_xml_output(self):
        self.assertEqual( self.xs_const.toXLink(), '/crossSection' )
        self.assertEqual( self.xs_const.toXMLList(), [
            '<crossSection>',
            '  <XYs1d label="eval">',
            '    <axes>',
            '      <axis index="1" label="energy_in" unit="eV"/>',
            '      <axis index="0" label="crossSection" unit="b"/></axes>',
            '    <values>',
            '      1.00000000e-05 1.00000000e+00 2.00000000e+07 1.00000000e+00</values></XYs1d></crossSection>'] )
    


if __name__=="__main__":
    unittest.main()
