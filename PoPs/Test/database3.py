#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import random

from PoPs import database as databaseModule

from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import lepton as leptonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nuclide as nuclideModule

from PoPs.chemicalElements import isotope as isotopeModule
from PoPs.chemicalElements import chemicalElement as chemicalElementModule

database = databaseModule.Database( 'test', '1.2.3' )

#
# Test adding an isotope to database.
#
photon = gaugeBosonModule.Particle( 'photon' )
database.add( photon )

#
# Test adding leptons to database.
#
lepton = leptonModule.Particle( 'e-', generation = 'electronic' )
database.add( lepton )

lepton = leptonModule.Particle( 'mu-_anti', generation = 'muonic' )
database.add( lepton )

#
# Test adding baryons to database.
#
baryon = baryonModule.Particle( 'p' )
database.add( baryon )

baryon = baryonModule.Particle( 'n_anti' )
database.add( baryon )

#
# Test adding a nuclear level to database.
#
level = nuclideModule.Particle( 'O16_e12' )
database.add( level )

#
# Test adding an isotope to database.
#
isotope = isotopeModule.Isotope( 'N15', 15 )
database.add( isotope )
level = nuclideModule.Particle( 'N15_e5' )
database.add( level )

#
# Test adding a nuclear level to an isotope.
#
level = nuclideModule.Particle( 'N15_e7' )
database.add( level )

#
# Test adding a chemical element to database.
#
chemicalElement = chemicalElementModule.ChemicalElement( 'Pu', 94, 'Plutonium' )
database.add( chemicalElement )
isotope = isotopeModule.Isotope( 'Pu238', 238 )
database.add( isotope )
level = nuclideModule.Particle( 'Pu238' )
database.add( level )
level = nuclideModule.Particle( 'Pu238_e2' )
database.add( level )

#
# Test adding an isotope to a chemical element.
#
isotope = isotopeModule.Isotope( 'Pu237', 237 )
database.add( isotope )

#
# Test adding a nuclear level to a chemical element.
#
level = nuclideModule.Particle( 'Pu239_e4' )
database.add( level )

xmld1 = database.toXML( )
print( xmld1 )
database2 = database.parseXMLString(xmld1)
if( xmld1 != database2.toXML( ) ) : raise Exception( 'Fix me.' )
