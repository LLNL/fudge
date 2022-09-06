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
from PoPs.chemicalElements import isotope as isotopeModule

A = 16
isotopeID = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( 'O', A )
isotope = isotopeModule.Isotope( isotopeID, A )

data = [ [ 0, 15.99491461956, 0,       None, None, None, None ],
         [ 1,           None, 6049400, None, None, None, None ],
         [ 2,           None, 6129893, None, None, None, None ],
         [ 3,           None, 6917100, None, None, None, None ] ]

for index, mass, energy, charge, halflife, spin, parity in data :
    name = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeID, index )
    level = nuclideModule.Particle( name )
    energy = nuclearEnergyLevelModule.Double( 'base', energy, quantityModule.stringToPhysicalUnit( 'eV' ) )
    level.nucleus.energy.add( energy )

    if( mass is not None ) :
        mass = massModule.Double( 'base', mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        level.mass.add( mass )

    if( charge is not None ) :
        charge = chargeModule.Integer( 'base', charge, quantityModule.stringToPhysicalUnit( 'e' ) )
        level.charge.add( charge )

    if( halflife is not None ) :
        halflife = halflifeModule.Double( 'base', halflife, quantityModule.stringToPhysicalUnit( 's' ) )
        level.halflife.add( halflife )

    if( spin is not None ) :
        spin = spinModule.Fraction( 'base', spinModule.Fraction.toValueType( spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
        level.spin.add( spin )

    if( parity is not None ) :
        parity = parityModule.Integer( 'base', parity, quantityModule.stringToPhysicalUnit( '' ) )
        level.parity.add( parity )

    isotope.add( level )

xmli1 = isotope.toXML( )
print( xmli1 )
isotope2 = isotopeModule.Isotope.parseXMLString(xmli1)
if( xmli1 != isotope2.toXML( ) ) : raise Exception( 'Fix me.' )

nuclide = isotope[isotopeID]
print( nuclide.mass.pqu( ) )
print( nuclide.mass.pqu( 'MeV/c**2' ) )
print( nuclide.mass.float( 'MeV/c**2' ) )

print()
print( nuclide.energy.pqu( ) )
print( nuclide.energy.pqu( 'MeV' ) )
print( PQUModule.floatToShortestString( nuclide.energy.float( 'MeV' ), 12 ) )
