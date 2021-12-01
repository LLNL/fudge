#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import random

from PoPs import database as databaseModule
from PoPs import alias as aliasModule

from PoPs.groups import misc as chemicalElementMiscModule

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import lepton as leptonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nuclide as nuclideModule

from PoPs.groups import isotope as isotopeModule
from PoPs.groups import chemicalElement as chemicalElementModule

database = databaseModule.database( 'test', '1.2.3' )

def nuclides( Z, A, data ) :

    symbol = chemicalElementMiscModule.symbolFromZ[Z]
    isotopeID = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( symbol, A )
    keys = random.sample( [ key for key in data ], len( data ) )
    for index in keys :
        name = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeID, index )
        nameLower = name[:1].lower( ) + name[1:]
        level = nuclideModule.particle( name )
        nucleus = level.nucleus

        atomicMass, nuclearMass, energy, charge, halflife, spin, parity = data[index]

        energy = nuclearEnergyLevelModule.double( 'base', energy, quantityModule.stringToPhysicalUnit( 'keV' ) )
        level.energy.add( energy )

        if( atomicMass is not None ) :
            mass = massModule.double( 'base', atomicMass, quantityModule.stringToPhysicalUnit( 'amu' ) )
            level.mass.add( mass )

        if( nuclearMass is not None ) :
            mass = massModule.double( 'base', nuclearMass, quantityModule.stringToPhysicalUnit( 'amu' ) )
            nucleus.mass.add( mass )

        if( charge is not None ) :
            charge = chargeModule.integer( 'base', charge, quantityModule.stringToPhysicalUnit( 'e' ) )
            nucleus.charge.add( charge )

        if( halflife is not None ) :
            if( halflife == 'stable' ) :
                halflife = halflifeModule.string( 'base', halflife, quantityModule.stringToPhysicalUnit( 's' ) )
            else :
                time, unit = halflife.split( )
                halflife = halflifeModule.double( 'base', float( time ), quantityModule.stringToPhysicalUnit( unit ) )
            nucleus.halflife.add( halflife )

        if( spin is not None ) :
            spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
            nucleus.spin.add( spin )

        if( parity is not None ) :
            parity = parityModule.integer( 'base', parity, quantityModule.stringToPhysicalUnit( '' ) )
            nucleus.parity.add( parity )

        database.add( level )

#                         atomMass    nuclearMass,     energy,   charge, halflife,   spin, parity
H1Data =   { 0 : [ 1.00782503207, 1.007276466879,          0,        1,  'stable',  '1/2',      1 ] }
nuclides( 1, 1, H1Data )

nuclearMass = 2.01410177785 - 0.0005485801
#                         atomMass    nuclearMass,     energy,   charge,  halflife,   spin, parity
H2Data =   { 0 : [ 2.01410177785,   nuclearMass,          0,        1,  'stable',   None,      1 ] }
nuclides( 1, 2, H2Data )

nuclearMass = 3.01604927767 - 0.0005485801
#                         atomMass    nuclearMass,     energy,   charge,  halflife,   spin, parity
H3Data =   { 0 : [ 3.01604927767,   nuclearMass,          0,        1, '12.32 yr',  '1/2',      1 ] }
nuclides( 1, 3, H3Data )

#                         atomMass    nuclearMass,     energy,   charge, halflife,   spin, parity
O16Data = { 0 : [ 15.99491461956,       None,             0,     None,     None,   None,   None ],
            1 : [           None,       None,      6049.400,     None,     None,   None,   None ],
            2 : [           None,       None,      6129.893,     None,     None,   None,   None ],
            3 : [           None,       None,      6917.100,     None,     None,   None,   None ] }
for A, data in [ [ 17, O16Data ], [ 16, O16Data ] ] : nuclides( 8, A, data )

#                    atomMass      nuclearMass,   energy, charge, halflife, spin, parity
Am242Data = { 0 : [ 242.059549159,      None,        0,   None, '16.02 h',    1,     -1 ],
              1 : [          None,      None,   44.092,   None,      None,    0,     -1 ],
              2 : [          None,      None,    48.60,   None,  '141 yr',    5,     -1 ],
              3 : [          None,      None,    52.70,   None,      None,    3,     -1 ] }
nuclides( 95, 242, Am242Data )
database.add( aliasModule.metaStable( 'Am242_m1', 'Am242_e2', 1 ) )

photon = gaugeBosonModule.particle( 'photon' )

mass = massModule.double( 'base', 0, quantityModule.stringToPhysicalUnit( 'amu' ) )
photon.mass.add( mass )

charge = chargeModule.integer( 'base', 0, quantityModule.stringToPhysicalUnit( 'e' ) )
photon.charge.add( charge )

halflife = halflifeModule.string( 'base', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
photon.halflife.add( halflife )

spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( '1' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
photon.spin.add( spin )

parity = parityModule.integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
photon.parity.add( parity )

database.add( photon )

leptons = [ [ 'e-',      5.48579909070e-4, -1,     'stable', '1/2',  1, 'electronic' ],
            [ 'e-_anti', 5.48579909070e-4,  1,     'stable', '1/2', -1, 'electronic' ],
            [ 'mu',      0.1134289267,     -1, 2.1969811e-6, '1/2',  1, 'muonic' ] ]

for id, mass, charge, halflife, spin, parity, generation in leptons :
    lepton = leptonModule.particle( id, generation = generation )

    mass = massModule.double( 'base', mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
    lepton.mass.add( mass )

    charge = chargeModule.integer( 'base', charge, quantityModule.stringToPhysicalUnit( 'e' ) )
    lepton.charge.add( charge )

    if( halflife == 'stable' ) :
        halflife = halflifeModule.string( 'base', halflife, quantityModule.stringToPhysicalUnit( 's' ) )
    else :
        halflife = halflifeModule.double( 'base', halflife, quantityModule.stringToPhysicalUnit( 's' ) )
    lepton.halflife.add( halflife )

    spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
    lepton.spin.add( spin )

    parity = parityModule.integer( 'base', parity, quantityModule.stringToPhysicalUnit( '' ) )
    lepton.parity.add( parity )

    database.add( lepton )

database.add( aliasModule.particle( 'electron', 'e-' ) )
database.add( aliasModule.particle( 'e+', 'e-_anti' ) )
database.add( aliasModule.particle( 'positron', 'e-_anti' ) )
database.add( aliasModule.particle( 'd', 'h2' ) )
database.add( aliasModule.particle( 't', 'h3' ) )

baryons = [ [ 'n', 1.00866491588,     0,    881.5, '1/2', 1 ],
            [ 'p', 1.007276466812,    1, 'stable', '1/2', 1 ] ]

for _id, _mass, _charge, _halflife, _spin, _parity in baryons :
    for anti in [ '', '_anti' ] :
        baryon = baryonModule.particle( _id + anti )

        mass = massModule.double( 'base', _mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        baryon.mass.add( mass )

        charge = chargeModule.integer( 'base', _charge, quantityModule.stringToPhysicalUnit( 'e' ) )
        baryon.charge.add( charge )

        if( _halflife == 'stable' ) :
            halflife = halflifeModule.string( 'base', _halflife, quantityModule.stringToPhysicalUnit( 's' ) )
        else :
            halflife = halflifeModule.double( 'base', _halflife, quantityModule.stringToPhysicalUnit( 's' ) )
        baryon.halflife.add( halflife )

        spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( _spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
        baryon.spin.add( spin )

        parity = parityModule.integer( 'base', _parity, quantityModule.stringToPhysicalUnit( '' ) )
        baryon.parity.add( parity )

        database.add( baryon )

xmld1 = database.toXML( )
print( xmld1 )
database2 = database.parseXMLStringAsClass( xmld1 )
if( xmld1 != database2.toXML( ) ) : raise Exception( 'Fix me.' )
