#! /usr/bin/env python

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

import random

from PoPs import database as databaseModule

from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import lepton as leptonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nuclide as nuclideModule

from PoPs.groups import isotope as isotopeModule
from PoPs.groups import chemicalElement as chemicalElementModule

database = databaseModule.database( 'test', '1.2.3' )

#
# Test adding an isotope to database.
#
photon = gaugeBosonModule.particle( 'photon' )
database.add( photon )

#
# Test adding leptons to database.
#
lepton = leptonModule.particle( 'e-', generation = 'electronic' )
database.add( lepton )

lepton = leptonModule.particle( 'mu-_anti', generation = 'muonic' )
database.add( lepton )

#
# Test adding baryons to database.
#
baryon = baryonModule.particle( 'p' )
database.add( baryon )

baryon = baryonModule.particle( 'n_anti' )
database.add( baryon )

#
# Test adding a nuclear level to database.
#
level = nuclideModule.particle( 'O16_e12' )
database.add( level )

#
# Test adding an isotope to database.
#
isotope = isotopeModule.isotope( 'N15', 15 )
database.add( isotope )
level = nuclideModule.particle( 'N15_e5' )
database.add( level )

#
# Test adding a nuclear level to an isotope.
#
level = nuclideModule.particle( 'N15_e7' )
database.add( level )

#
# Test adding a chemical element to database.
#
chemicalElement = chemicalElementModule.chemicalElement( 'Pu', 94, 'Plutonium' )
database.add( chemicalElement )
isotope = isotopeModule.isotope( 'Pu238', 238 )
database.add( isotope )
level = nuclideModule.particle( 'Pu238' )
database.add( level )
level = nuclideModule.particle( 'Pu238_e2' )
database.add( level )

#
# Test adding an isotope to a chemical element.
#
isotope = isotopeModule.isotope( 'Pu237', 237 )
database.add( isotope )

#
# Test adding a nuclear level to a chemical element.
#
level = nuclideModule.particle( 'Pu239_e4' )
database.add( level )

xmld1 = database.toXML( )
print xmld1
database2 = database.parseXMLStringAsClass( xmld1 )
if( xmld1 != database2.toXML( ) ) : raise Exception( 'Fix me.' )
