#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import halflife as halflifeModule

from PoPs.families import lepton as leptonModule

electron = leptonModule.particle( 'e-', generation = 'electronic' )

mass = massModule.double( 'base', 5.48579909070e-4, quantityModule.stringToPhysicalUnit( 'amu' ) )
electron.mass.add( mass )

charge = chargeModule.integer( 'base', -1, quantityModule.stringToPhysicalUnit( 'e' ) )
electron.charge.add( charge )

halflife = halflifeModule.string( 'base', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
electron.halflife.add( halflife )

spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( '1/2' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
electron.spin.add( spin )

parity = parityModule.integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
electron.parity.add( parity )

print( electron.toXML( ) )
print()
electron2 = leptonModule.particle.parseXMLStringAsClass( electron.toXML( ) )
if( electron.toXML( ) != electron2.toXML( ) ) : raise Exception( 'Fix me' )

suite = leptonModule.suite( )
suite.add( electron )
print( suite.toXML( ) )
suite2 = leptonModule.suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )
