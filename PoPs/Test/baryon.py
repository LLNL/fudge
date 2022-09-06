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

from PoPs.families import baryon as baryonModule

proton = baryonModule.Particle( 'p' )

mass = massModule.Double( 'base', 1.007276466812, quantityModule.stringToPhysicalUnit( 'amu' ) )
proton.mass.add( mass )

charge = chargeModule.Integer( 'base', 1, quantityModule.stringToPhysicalUnit( 'e' ) )
proton.charge.add( charge )

halflife = halflifeModule.String( 'base', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
proton.halflife.add( halflife )

spin = spinModule.Fraction( 'base', spinModule.Fraction.toValueType( '1/2' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
proton.spin.add( spin )

parity = parityModule.Integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
proton.parity.add( parity )

print( proton.toXML( ) )
print()
proton2 = baryonModule.Particle.parseXMLString(proton.toXML())
if( proton.toXML( ) != proton2.toXML( ) ) : raise Exception( 'Fix me' )

suite = baryonModule.Suite( )
suite.add( proton )
print( suite.toXML( ) )
suite2 = baryonModule.Suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )
