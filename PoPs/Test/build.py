#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import random
import fractions

from PoPs import database as databaseModule
from PoPs import alias as aliasModule
from PoPs import misc as miscModule
from PoPs import IDs as IDsModule

from PoPs.quantities import spin as spinModule

from PoPs.families import lepton as leptonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nucleus as nucleusModule
from PoPs.families import nuclide as nuclideModule

oneHalf = fractions.Fraction( '1/2' )
one = fractions.Fraction( '1' )
spinUnit = spinModule.baseUnit

database = databaseModule.Database( 'test', '1.2.3' )

particle = miscModule.buildParticleFromRawData( leptonModule.Particle, IDsModule.electron,
        mass = ( 0.0005485801, 'amu' ),         spin = ( oneHalf, spinUnit ),                     parity = ( 1, '' ),
        charge = ( -1, 'e' ),                   halflife = ( 'stable', 's' ),               generation = 'electronic' )
database.add( particle )

particle = miscModule.buildParticleFromRawData( leptonModule.Particle, IDsModule.electron + '_anti',
        mass = ( 0.0005485801, 'amu' ),         spin = ( oneHalf, spinUnit ),                     parity = ( -1, '' ),
        charge = (  1, 'e' ),                   halflife = ( 'stable', 's' ),               generation = 'electronic' )
database.add( particle )

particle = miscModule.buildParticleFromRawData( leptonModule.Particle, 'mu',
        mass = ( 0.1134289267, 'amu' ),         spin = ( oneHalf, spinUnit ),                     parity = ( 1, '' ),
        charge = ( -1, 'e' ),                   halflife = ( 2.1969811e-6, 's' ),           generation = 'muonic' )
database.add( particle )

database.add(aliasModule.Alias('electron', 'e-'))
database.add(aliasModule.Alias('e+', 'e-_anti'))
database.add(aliasModule.Alias('positron', 'e-_anti'))

particle = miscModule.buildParticleFromRawData( baryonModule.Particle, 'n',
        mass = ( 1.00866491588, 'amu' ),         spin = ( oneHalf, spinUnit ),                    parity = ( 1, '' ),
        charge = (  0, 'e' ),                   halflife = ( 881., 's' ) )
database.add( particle )

particle = miscModule.buildParticleFromRawData( baryonModule.Particle, 'p',
        mass = ( 1.007276466812, 'amu' ),         spin = ( oneHalf, spinUnit ),                   parity = ( 1, '' ),
        charge = (  1, 'e' ),                   halflife = ( 'stable', 's' ) )
database.add( particle )

nucleus = miscModule.buildParticleFromRawData( nucleusModule.Particle, 'O16', index = 0, energy = ( 0.0, 'eV' ) )
particle = miscModule.buildParticleFromRawData( nuclideModule.Particle, 'O16',
        mass = ( 15.994913988, 'amu' ), energy = ( 0.0, 'eV' ) )
database.add( particle )

nucleus = miscModule.buildParticleFromRawData( nucleusModule.Particle, 'o16_e3', index = 3, energy = ( 6917100.0, 'eV' ) )
particle = miscModule.buildParticleFromRawData( nuclideModule.Particle, 'O16_e3', nucleus = nucleus )
database.add( particle )

xmld1 = database.toXML( )
print( xmld1 )
database2 = database.parseXMLString(xmld1)

if( xmld1 != database2.toXML( ) ) : raise Exception( 'Fix me.' )
