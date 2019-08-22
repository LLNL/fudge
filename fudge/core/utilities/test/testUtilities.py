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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
test fudge/core/utilities
cmattoon, 3/24/2011
"""

import unittest
from fudge.core import utilities

class fudgeZA_test(unittest.TestCase):
    def setUp(self):
        self.ZA = utilities.fudgeZA
    
    def test1(self):
        # check some random Z values:
        for z,sym in ((20,"Ca"),(94,"Pu"),(3,"Li")):
            self.assertEqual( self.ZA.ZToSymbol(z), sym )
            self.assertEqual( self.ZA.SymbolToZ(sym), z )
        
        for z,el in ((111,"Roentgenium"),(77,"Iridium"),(43,"Technetium")):
            self.assertEqual( self.ZA.ZToLabel(z), el )
            self.assertEqual( self.ZA.LabelToZ(el), z )
        
    def test2(self):
        # z out of bounds
        max_z_value = self.ZA.nZs
        self.assertEqual( self.ZA.ZToSymbol(max_z_value+1), None )
        self.assertEqual( self.ZA.ZToSymbol(-1), None )
        self.assertEqual( self.ZA.SymbolToZ('noSuchSym'), None )
        self.assertEqual( self.ZA.LabelToZ('noSuchSym'), None )

class brb_test(unittest.TestCase):
    def setUp(self):
        self.brb = utilities.brb
    
    def test1(self):
        # ???
        pass

class testReactionStrings(unittest.TestCase):
    def setUp(self):
        self.RS = utilities.reactionStrings
    
    def test1(self):
        for reac in ("n + Th232 -> Th233 + gamma",
            "n + Th232 -> n + (Th232_e2 -> He4 + (Ra228 -> He4 + Rn224))",
            "H2 + H3 -> n + He4",
            "n + C11 -> (C12 -> He4 + (Be8 -> He4[multiplicity:'2']))",
            "n + Fe56 -> n + Fe56 [compound]",
            "gamma + Pu239 -> n[multiplicity:'energyDependent', decayRate:'1.45210000e-02'] + Pu238",
            "n + Cl35 -> n + Cl35_s",
            "n + Pu239 -> n[multiplicity:'energyDependent', emissionMode:'prompt'] + n[multiplicity:'energyDependent', emissionMode:'delayed', decayRate:'1.32710000e-02'] + n[multiplicity:'energyDependent', emissionMode:'delayed', decayRate:'3.08810000e-02'] + n[multiplicity:'energyDependent', emissionMode:'delayed', decayRate:'1.13370000e-01'] + n[multiplicity:'energyDependent', emissionMode:'delayed', decayRate:'2.92500000e-01'] + n[multiplicity:'energyDependent', emissionMode:'delayed', decayRate:'8.57490000e-01'] + n[multiplicity:'energyDependent', emissionMode:'delayed', decayRate:'2.72970000e+00'] + gamma[multiplicity:'energyDependent'] [total fission]"
            ):
            self.assertEqual( reac, str(self.RS.parseReaction(reac)) )
    
    def test2(self):
        # should handle extra spaces, etc:
        string = "H2+Rb87  ->n[multiplicity: 'energyDependent'][fake string]"
        parsed = self.RS.parseReaction(string)
        self.assertEqual( str(parsed.target), "Rb87" )
        self.assertEqual( parsed.info, "fake string" )
        self.assertEqual( str(parsed), "H2 + Rb87 -> n[multiplicity:'energyDependent'] [fake string]" )
        pass
if __name__ == '__main__':
    unittest.main()

