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

from PoPs.families import gaugeBoson as gaugeBosonModule

photon = gaugeBosonModule.Particle( 'photon' )

mass = massModule.Double( 'base', 0, quantityModule.stringToPhysicalUnit( 'amu' ) )
photon.mass.add( mass )

charge = chargeModule.Integer( 'base', 0, quantityModule.stringToPhysicalUnit( 'e' ) )
photon.charge.add( charge )

halflife = halflifeModule.Double( 'base', 1e100, quantityModule.stringToPhysicalUnit( 's' ) )
photon.halflife.add( halflife )

spin = spinModule.Fraction( 'base', spinModule.Fraction.toValueType( '1' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
photon.spin.add( spin )

parity = parityModule.Integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
photon.parity.add( parity )

print( photon.toXML( ) )
print()
photon2 = gaugeBosonModule.Particle.parseXMLString(photon.toXML())
if( photon.toXML( ) != photon2.toXML( ) ) : raise Exception( 'Fix me' )

suite = gaugeBosonModule.Suite( )
suite.add( photon )
print( suite.toXML( ) )
suite2 = gaugeBosonModule.Suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )
