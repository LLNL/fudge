#! /usr/bin/env python3
# encoding: utf-8

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
test fudge/reactionData/
dbrown, 12/5/2012
"""

import unittest
from fudge.reactionData.crossSection import *
defaultAccuracy = 0.001
from xData import XYs1d as XYs1dModule

class TestPointwise( unittest.TestCase ):

    def test_basic_stuff(self):
        ptwise_const=XYs1d(axes=XYs1d.defaultAxes(labelsUnits={
            XYs1dModule.yAxisIndex : ( 'crossSection', 'b' ),
            XYs1dModule.xAxisIndex : ( 'energy_in', 'eV' ) }), data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.assertEqual(ptwise_const.domainUnit, 'eV' )
        self.assertEqual(ptwise_const.toXML_strList(), [
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
        ptwise_const=XYs1d(axes=XYs1d.defaultAxes(labelsUnits={
                XYs1dModule.yAxisIndex : ( 'crossSection', 'b' ),
                XYs1dModule.xAxisIndex : ( 'energy_in', 'eV' ) }), data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.assertEqual(ptwise_const.evaluate( 10.0e6 ), 1.0 )

    def test_applyFunction(self):
        """
        Simple test of applyFunction() function.  The function we'll apply is exp(-x) so this is equivalent to the test_exp tests above

        Test function: y(x)=x/1e6+1.0 (the 1e6 converts eV to MeV)

        So, z(x)=exp(-y(x))=exp(-x/1e6-1.0)

        Expected results:

              x     y(x)      z(y(x))
            -------------------------------
             1.0e6   2.0      0.1353352832366127
            15.0e6  16.0      1.1253517471925912e-07
        """

        ptwise_linear = XYs1d(axes=XYs1d.defaultAxes(labelsUnits={
            XYs1dModule.yAxisIndex: ('crossSection', 'b'),
            XYs1dModule.xAxisIndex: ('energy_in', 'eV')}), data=[[1e-5, 1.0], [20.0e6, 21.0]])

        self.assertAlmostEqual(ptwise_linear.evaluate(15.0e6), 16.0)
#        self.assertAlmostEqual(ptwise_linear.applyFunction(lambda x, y: math.exp(-x), None).evaluate(15.0e6), math.exp(-16.0))  # This should work, but fails
        self.assertAlmostEqual(ptwise_linear.evaluate(1.0e6), 2.0)
#        self.assertAlmostEqual(ptwise_linear.applyFunction(lambda x, y: math.exp(-x), None).evaluate(1.0e6), math.exp(-2.0))    # This should work, but fails
        self.assertAlmostEqual(ptwise_linear.applyFunction(lambda x, y: math.exp(-ptwise_linear.evaluate(x)), None).evaluate(1.0e6), math.exp(-2.0), 3)  # This should absolutely fail and is the wrong way to do it


class TestCrossSection( unittest.TestCase ):

    def setUp(self):
    
        # A constant (kind of like (n,el))
        self.xs_const = Component()
        self.ptwise_const=XYs1d(label='eval',axes=XYs1d.defaultAxes(labelsUnits={
            XYs1dModule.yAxisIndex : ( 'crossSection', 'b' ),
            XYs1dModule.xAxisIndex : ( 'energy_in', 'eV' ) }),
            data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.xs_const.add(self.ptwise_const)

        # A triangle (kind of like (n,2n), if you squint real hard)
        self.xs_triangle = Component()
        self.ptwise_triangle=XYs1d(label='eval',axes=XYs1d.defaultAxes(labelsUnits={
            XYs1dModule.yAxisIndex : ( 'crossSection', 'b' ),
            XYs1dModule.xAxisIndex : ( 'energy_in', 'eV' ) }),
            data=[ [1e6,0.0], [10.0e6,1.0], [20.0e6,0.0] ] )
        self.xs_triangle.add(self.ptwise_triangle)

    def test_basic_stuff(self):
        self.assertTrue( isinstance( self.xs_const, Component ) )
        self.assertEqual( self.xs_const.moniker, 'crossSection' )
        self.assertEqual( self.xs_const.domainUnit, 'eV' )
        
    def test_xlink(self):
        self.assertEqual( self.xs_const.toRelativeXLink(), '' )
        
    def test_parent(self):
        """
        Tests of traversing up GNDS hierarchy
        """
        self.assertEqual( self.xs_const.ancestor, None )
        self.assertEqual( self.xs_const.rootAncestor, self.xs_const )
        self.assertEqual( self.xs_const.ancestor, None )
        self.xs_const.setAncestor( 'fred' )
        self.assertEqual( self.xs_const.ancestor, 'fred' )

    def test_ancestry(self):
        """
        FIXME: findAttributeInAncestry
        FIXME: findClassInAncestry
        """
        pass

    def test_forms(self):

        listAssertEqual = getattr( self, 'assertCountEqual', None )
        if listAssertEqual is None: listAssertEqual = getattr( self, 'assertItemsEqual')

        listAssertEqual( [ x.moniker for x in self.xs_const.allowedClasses ],
                [ 'Ys1d', 'XYs1d', 'regions1d', 'gridded1d', 'resonancesWithBackground', 'reference',
                  'URR_probabilityTables1d', 'CoulombPlusNuclearElastic', 'thermalNeutronScatteringLaw1d' ] )
        self.xs_const.remove('gridded') # ??? why does this test have a side effect?

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
        info = {
            'Q': 1.0,
            'kinematicFactor': 1.0,
            'dThreshold': PQU.PQU('1e-4 b'),
            'crossSectionEnergyMin': PQU.PQU('1.0e-5 eV'),
            'crossSectionEnergyMax': PQU.PQU('20.0 MeV'),
            'CoulombChannel': False,
            'ContinuumOutputChannel': False
            }
        self.assertEqual( self.xs_const.check(info), [] )

    def test_xml_output(self):
        self.assertEqual( self.xs_const.toXLink(), '/crossSection' )
        self.assertEqual( self.xs_const.toXML_strList(), [
            '<crossSection>',
            '  <XYs1d label="eval">',
            '    <axes>',
            '      <axis index="1" label="energy_in" unit="eV"/>',
            '      <axis index="0" label="crossSection" unit="b"/></axes>',
            '    <values>',
            '      1.00000000e-05 1.00000000e+00 2.00000000e+07 1.00000000e+00</values></XYs1d></crossSection>'] )
    


if __name__=="__main__":
    unittest.main()
