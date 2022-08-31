#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import random

from PoPs import database as databaseModule

from PoPs.chemicalElements import misc as chemicalElementMiscModule

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import nuclide as nuclideModule

database = databaseModule.Database( 'test', '1.2.3' )

O16Data = { 0 : [ 15.99491461956, 0,       None, None, None, None ],
            1 : [           None, 6049400, None, None, None, None ],
            2 : [           None, 6129893, None, None, None, None ],
            3 : [           None, 6917100, None, None, None, None ] }


symbol = chemicalElementMiscModule.symbolFromZ[8]
for A, data in [ [ 17, O16Data ], [ 16, O16Data ] ] :
    isotopeID = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( symbol, A )
    keys = random.sample( [ key for key in data ], len( data ) )
    for index in keys :
        mass, energy, charge, halflife, spin, parity = data[index]

        name = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeID, index )
        level = nuclideModule.Particle( name )
        energy = nuclearEnergyLevelModule.Double( 'base', energy, quantityModule.stringToPhysicalUnit( 'keV' ) )
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

        database.add( level )

photon = gaugeBosonModule.Particle( 'photon' )

mass = massModule.Double( 'base', 0, quantityModule.stringToPhysicalUnit( 'amu' ) )
photon.mass.add( mass )

charge = chargeModule.Integer( 'base', 0, quantityModule.stringToPhysicalUnit( 'e' ) )
photon.charge.add( charge )

halflife = halflifeModule.String( 'base', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
photon.halflife.add( halflife )

spin = spinModule.Fraction( 'base', spinModule.Fraction.toValueType( '1' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
photon.spin.add( spin )

parity = parityModule.Integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
photon.parity.add( parity )

database.add( photon )

print( database['photon'] )
print( database.chemicalElements.getSymbol( 'O' ) )
print( database['O16'] )
print( database['O16_e2'] )
try :
    print( database['O15'] )
    testFailed = True
except :
    testFailed = False
if( testFailed ) : raise Exception( 'particle "O15" found when it should not have been' )
