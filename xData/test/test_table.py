# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import unittest

import xData.table as table
from xml.etree import cElementTree as parser


class TestTable(unittest.TestCase):

    def setUp(self):
        self.tt = table.table( columns=[
            table.columnHeader( 0, 'energy', 'eV' ),
            table.columnHeader( 1, 'spin', '' ),
            table.columnHeader( 2, 'neutronWidth', 'eV' ),
            table.columnHeader( 3, 'captureWidth', 'eV' ), ] )
        for dat in ([-1.7, table.blank(), 3.2, 0.089], [3.4, 0.5, 4.2, 0.072], [5.6, 1.5, 2.76, 0.064]):
            self.tt.addRow( dat )
        self.xmlstring = '\n'.join( self.tt.toXMLList() )
        element = parser.fromstring( self.xmlstring )
        self.tt2 = table.table.parseXMLNode( element, xPath = [], linkData = {'conversionTable' : {'index':int} } )
        self.xmlstring2 = '\n'.join( self.tt2.toXMLList() )


    def test_table(self):
        #print(self.xmlstring)
        self.assertEqual(self.xmlstring, self.xmlstring2)

    def test_getColumn(self):
        # extract a column, converting eV -> MeV
        self.assertListEqual(self.tt2.getColumn('captureWidth', 'MeV'), [8.899999999999999e-08, 7.2e-08, 6.4e-08])
        # THIS^ WILL NOT WORK WITH blank() OBJECT TYPES, ONLY TESTS IF COLUMN IS EXTRACTED AND UNITS ARE CONVERTED

    def test_addrow(self):
        """ Test if the ValueError raises when adding row of unequal length"""
        with self.assertRaises(ValueError):
            self.tt.addRow([1,1,1])


# to run unit test
if __name__=="__main__":
    unittest.main()
