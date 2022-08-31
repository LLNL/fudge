#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from PoPs.chemicalElements import misc as chemicalElementMiscModule

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from PoPs.families import nuclide as nuclideModule 

index = 3
name = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( 'O16', index )

level = nuclideModule.Particle( name )
energy = nuclearEnergyLevelModule.Double( 'base', 6917100, quantityModule.stringToPhysicalUnit( 'eV' ) )
level.nucleus.energy.add( energy )

mass = massModule.Double( 'base', 15.99491461956, quantityModule.stringToPhysicalUnit( 'amu' ) )
level.mass.add( mass )

charge = chargeModule.Integer( 'base', 0, quantityModule.stringToPhysicalUnit( 'e' ) )
level.charge.add( charge )

halflife = halflifeModule.Double( 'base', 1e100, quantityModule.stringToPhysicalUnit( 's' ) )
level.halflife.add( halflife )

spin = spinModule.Fraction( 'base', spinModule.Fraction.toValueType( '5/2' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
level.spin.add( spin )

parity = parityModule.Integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
level.parity.add( parity )

xmll1 = level.toXML( )
print( xmll1 )
level2 = nuclideModule.Particle.parseXMLString(level.toXML())
if( xmll1 != level2.toXML( ) ) : raise Exception( 'Fix me.' )

print( level.mass[0].pqu( ) )
print( level.mass.pqu( ) )
print( level.mass[0].pqu( 'MeV/c**2' ) )
print( level.mass.pqu( 'MeV/c**2' ) )
print( PQUModule.floatToShortestString( level.mass[0].float( 'MeV/c**2' ), 12 ) )
print( PQUModule.floatToShortestString( level.mass.float( 'MeV/c**2' ), 12 ) )
