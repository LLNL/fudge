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

from PoPs.families import gaugeBoson as gaugeBosonModule

photon = gaugeBosonModule.particle( 'photon' )

mass = massModule.double( 'base', 0, quantityModule.stringToPhysicalUnit( 'amu' ) )
photon.mass.add( mass )

charge = chargeModule.integer( 'base', 0, quantityModule.stringToPhysicalUnit( 'e' ) )
photon.charge.add( charge )

halflife = halflifeModule.double( 'base', 1e100, quantityModule.stringToPhysicalUnit( 's' ) )
photon.halflife.add( halflife )

spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( '1' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
photon.spin.add( spin )

parity = parityModule.integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
photon.parity.add( parity )

print( photon.toXML( ) )
print()
photon2 = gaugeBosonModule.particle.parseXMLStringAsClass( photon.toXML( ) )
if( photon.toXML( ) != photon2.toXML( ) ) : raise Exception( 'Fix me' )

suite = gaugeBosonModule.suite( )
suite.add( photon )
print( suite.toXML( ) )
suite2 = gaugeBosonModule.suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )
