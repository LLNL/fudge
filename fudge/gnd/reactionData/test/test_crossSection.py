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
test fudge/gnd/reactionData/
dbrown, 12/5/2012
"""

import unittest
from fudge.gnd.reactionData.crossSection import *
defaultAccuracy = 0.001
from xData import XYs

class TestPointwise( unittest.TestCase ):

    def test_basic_stuff(self):
        ptwise_const=pointwise(axes=pointwise.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( 'crossSection', 'b' ), \
                XYs.xAxisIndex : ( 'energy_in', 'eV' ) }), data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.assertEqual(ptwise_const.domainUnit( ), 'eV' )
        self.assertEqual(ptwise_const.toXMLList(), [\
            '<XYs>',\
            '  <axes>',\
            '    <axis index="1" label="energy_in" unit="eV"/>',\
            '    <axis index="0" label="crossSection" unit="b"/></axes>',\
            '  <values length="4">1e-5 1 2e7 1</values></XYs>'] )
            
    def test_group(self):
        ptwise=pointwise(axes=pointwise.defaultAxes(), data=[
            [1e-5,1.0], 
            [1.0,1.0], 
            [10.0,2.0], 
            [100.0,3.0], 
            [1000.0,1.0], 
            [20.0e6,1.0] ] )
#        ptwise.plot()
        ptwise.group( [1e-5,1.0,1000.0,20.e6] )#.plot()

class TestCrossSection( unittest.TestCase ): 

    def setUp(self):
    
        # A constant (kind of like (n,el))
        self.xs_const = component()
        self.ptwise_const=pointwise(axes=pointwise.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( 'crossSection', 'b' ), \
                XYs.xAxisIndex : ( 'energy_in', 'eV' ) }), \
            data=[ [1e-5,1.0], [20.0e6,1.0] ] )
        self.xs_const.add(self.ptwise_const)

        # A triangle (kind of like (n,2n), if you squint real hard)
        self.xs_triangle = component()
        self.ptwise_triangle=pointwise(axes=pointwise.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( 'crossSection', 'b' ), \
                XYs.xAxisIndex : ( 'energy_in', 'eV' ) }), \
            data=[ [1e6,0.0], [10.0e6,1.0], [20.0e6,0.0] ] )
        self.xs_triangle.add(self.ptwise_triangle)

    def test_basic_stuff(self):
        self.assertTrue( isinstance( self.xs_const, component ) )
        self.assertEqual( self.xs_const.attribute, None )
        self.assertEqual( self.xs_const.moniker, 'crossSection' )
        self.assertEqual( self.xs_const.domainUnit(), 'eV' )
        
    def test_xlink(self):
        self.assertEqual( self.xs_const.toRelativeXLink(), '' )
        
    def test_parent(self):
        '''
        Tests of traversing up GND hierarchy
        '''
        self.assertEqual( self.xs_const.ancestor, None )
        self.assertEqual( self.xs_const.getRootAncestor(), self.xs_const )
        self.assertEqual( self.xs_const.getAncestor(), None )
        self.xs_const.setAncestor( 'fred' )
        self.assertEqual( self.xs_const.getAncestor(), 'fred' )

    def test_ancestry(self):
        '''
        FIXME: findAttributeInAncestry
        FIXME: findClassInAncestry
        '''
        pass

    def test_forms(self):
        self.assertItemsEqual( [x.moniker for x in self.xs_const.allowedClasses], ['XYs','regions','multiGroup','resonanceRegion','resonancesWithBackground', 'weightedPointwise', 'reference'] )
        self.assertEqual( self.xs_const.getFormTokens(), ['pointwise'] )
        self.assertEqual( self.xs_const.getFormByToken('pointwise'), self.ptwise_const )
        self.xs_const.removeStyle('pointwise',force=True)
        self.assertEqual( self.xs_const.getFormTokens(), [] )
                
    def test_genre(self):
        self.assertEqual( self.xs_const.genre, 'crossSection' )
        self.assertEqual( self.xs_const.getGenre(), 'crossSection' )
        
    def test_nativeData(self):
        print dir(self.xs_const)
        self.assertEqual( self.xs_const.getNativeDataToken(), 'pointwise' )
        self.assertEqual( self.xs_const.nativeData, 'pointwise' )
        self.assertEqual( self.xs_const.getNativeData(), self.ptwise_const )
        self.assertEqual( self.xs_const.getNativeData().getValue(1.4e5), 1.0 )
        
    def test_domains(self):
        self.assertEqual( self.xs_const.domainMax(), 20.0e6 )
        self.assertEqual( self.xs_const.domainMin(), 1.0e-5 )
        self.assertEqual( self.xs_const.domainUnitConversionFactor('MeV'), 1e-6 )
        self.assertEqual( self.xs_const.domain(), (1e-05, 20000000.0) )
        self.assertEqual( self.xs_const.domainUnit(), 'eV' )

    def test_check(self):
        info = {'Q':1.0,'kinematicFactor':1.0,'dThreshold':PQU.PQU('1e-4 b'),'crossSectionEnergyMax':PQU.PQU('20.0 MeV')}
        self.assertEqual( self.xs_const.check(info), [] )

    def test_xml_output(self):
        self.assertEqual( self.xs_const.toXLink(), '/crossSection' )
        self.assertEqual( self.xs_const.toXMLList(), [
            '<crossSection>',
            '  <XYs>',
            '    <axes>',
            '      <axis index="1" label="energy_in" unit="eV"/>',
            '      <axis index="0" label="crossSection" unit="b"/></axes>',
            '    <values length="4">1e-5 1 2e7 1</values></XYs></crossSection>'] )
    


if __name__=="__main__":
    unittest.main()
