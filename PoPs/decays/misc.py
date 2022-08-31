# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Helper method for photon decay data
"""

from LUPY import enums as enumsModule

from .. import IDs as IDsModule
from ..families import nucleus as nucleusModule
from ..families import nuclide as nuclideModule

class Mode(enumsModule.Enum):

    electroMagnetic = enumsModule.auto()

def photonBranchingData( pops, id ) :
    """
    This function summarizes all photon branches for the specified nuclide. The returned object is a dictionary with 
    keys 'energy' and 'photons'. The value of the 'energy' key is the energy level of the nucleus. The value of the 
    'photons' key is a list of information about the the photons (i.e., gammas) emittied from the decay of the nucleus. 
    For each photo, the information is the tuple branching ratio, pqu of the photon energy and final nuclear state. 
    For example, for Fe56_e9 in ENDF/B-VIII.0, the data in the 'photons' key is

    [ ( 0.013, PQU( "0.265580000000000 MeV" ), 'Fe56_e7' ), 
      ( 0.987, PQU( "1.303445 MeV" ), 'Fe56_e2' ) ]

    These data state that the Fe56_e9 nuclear level decays to Fe56_e7 0.13% of time emitting a photon of energy 0.26558 MeV
    and it decays to Fe56_e2 98.7% of time emitting a photon of energy 1.303445 MeV.

    :return: dictionary with level energy and photon branching data for each nuclide in the isotope.
    """

    nuclide = pops[id]
    if( isinstance( nuclide, nucleusModule.Particle ) ) : nuclide = nuclide.nuclide
    if( not( isinstance( nuclide, nuclideModule.Particle ) ) ) : raise TypeError( 'id "%s" not a nuclide or nucleus' % id )

    branchingData = {}
    for nuclide in nuclide.isotope :
        photons = []
        energy = nuclide.energy[0].pqu( )
        for decayMode in nuclide.decayData.decayModes :
            if decayMode.mode == Mode.electroMagnetic:
                probability = decayMode.probability[0].value
                decayPath = decayMode.decayPath[0]
                products = [ decayProduct.pid for decayProduct in decayPath.products ]
                products.remove( IDsModule.photon )
                residual = products[0]
                photons.append( ( probability, energy - branchingData[residual]['energy'], residual ) )
        branchingData[nuclide.id] = { 'energy' : energy, 'photons' : photons }
    return( branchingData )
