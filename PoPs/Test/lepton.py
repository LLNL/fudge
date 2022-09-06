#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

electron = leptonModule.Particle( 'e-', generation = 'electronic' )

mass = massModule.Double( 'base', 5.48579909070e-4, quantityModule.stringToPhysicalUnit( 'amu' ) )
electron.mass.add( mass )

charge = chargeModule.Integer( 'base', -1, quantityModule.stringToPhysicalUnit( 'e' ) )
electron.charge.add( charge )

halflife = halflifeModule.String( 'base', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
electron.halflife.add( halflife )

spin = spinModule.Fraction( 'base', spinModule.Fraction.toValueType( '1/2' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
electron.spin.add( spin )

parity = parityModule.Integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
electron.parity.add( parity )

print( electron.toXML( ) )
print()
electron2 = leptonModule.Particle.parseXMLString(electron.toXML())
if( electron.toXML( ) != electron2.toXML( ) ) : raise Exception( 'Fix me' )

suite = leptonModule.Suite( )
suite.add( electron )
print( suite.toXML( ) )
suite2 = leptonModule.Suite.parseXMLString( suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )
