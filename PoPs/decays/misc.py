# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module contains general enums and functions of use for decay data.
'''

from LUPY import enums as enumsModule

from .. import IDs as IDsModule
from ..families import nucleus as nucleusModule
from ..families import nuclide as nuclideModule

class Mode(enumsModule.Enum):

    electroMagnetic = enumsModule.auto()

def completePhotonBranchingData(initialState, branchingData, probability, completePhotons):
    '''
    For *initialState*, determines the multiplicity (i.e., total probability) that a photon from
    level *initialState* or any level it decays to is emitted. Data are added to *completePhotons*
    as needed. *completePhotons* is a dictionary with keys "inital:final" where "initial" and final

    represent the decay of initial to final leading to the emission of a photon. For example,
    for the decay from level "Fe56_e6" to level "Fe56_e3" the key is "Fe56_e6:Fe56_e3". The data
    for each key is the energy difference of the initial and final levels and the multiplicity (i.e., 
    totoal probability) that the photon was emitted. For example, for Fe56 in ENDF/B-VIII.0., the
    data for the key "Fe56_e6:Fe56_e3" are "[PQU( "462521.000000000 eV" ), 0.009888135544683001]".

    This function is mainly for use by the function **photonBranchingData** to fill in its "completePhotons" data.

    :param initialState:            The current state (nuclide) the decaying isotope is in.
    :param branchingData:           Branching data calculated by the function **photonBranchingData**.
    :param probability:             The probability of reaching level *initialState*.
    :param completePhotons:         The list of all emitted photon data.
    '''

    gammas = branchingData[initialState]['photons']
    for branchingRatio, gammaEnergy, finalState, photonEmissionProbability in gammas:
        key = '%s:%s' % (initialState, finalState)
        if key not in completePhotons:
            completePhotons[key] = [gammaEnergy, 0.0]
        completePhotons[key][1] += branchingRatio * probability * photonEmissionProbability
        completePhotonBranchingData(finalState, branchingData, branchingRatio * probability, completePhotons)

def photonBranchingData(pops, id):
    '''
    This function summarizes all photon branches for the specified nuclide. The returned object is a dictionary with 
    keys "energy", "photons" and "completePhotons". The value of the "energy" key is the energy level of the nucleus. 
    The value of the "photons" key is a list of information about the the photons (i.e., gammas) emittied from the decay
    of the nucleus.  For each photo, the information is the tuple branching ratio, **PQU** of the photon energy, final 
    nuclear state and photon emission probability. For example, for Fe56_e9 in ENDF/B-VIII.0, the data in the "photons" key is:

    [ ( 0.013, PQU( "0.265580000000000 MeV" ), "Fe56_e7", 0.9875568 ), 
      ( 0.987, PQU( "1.303445 MeV" ), "Fe56_e2", 0.999862 ) ]

    These data state that the Fe56_e9 nuclear level decays to Fe56_e7 0.13% of time emitting a photon of 
    energy 0.26558 MeV 98.75568% of the time, and it decays to Fe56_e2 98.7% of time emitting a photon of 
    energy 1.303445 MeV 99.9862% of the time.

    The value of the "completePhotons" key is a dictionary. Each key in the dictionary is comprised of the initial
    and final state ids for each emtted photon in the format "intial:final" (e.g., "Fe56_e6:Fe56_e3") which means 
    that it is the photon emitted from the transition from the initial level to the final level (e.g., from level 
    "Fe56_e6" to level "Fe56_e3"). The data for each key is the energy of the emitted photon and its multiplicity.
    The multiplicity is the total probability that this photon is emitted.
    These data are obtained by calling the function **completePhotonBranchingData**.

    :param pops:        A PoPs database.
    :param id:          The id of the requested nuclide.

    :return:            Dictionary with level energy and photon branching data for each nuclide in the isotope.

    :raises TypeError:  If id is in *pops* but not of family nuclide or nucleus.
    '''

    nuclide = pops[id]
    if isinstance(nuclide, nucleusModule.Particle):
        nuclide = nuclide.nuclide
    if not isinstance(nuclide, nuclideModule.Particle):
        raise TypeError('Id "%s" not a nuclide or nucleus' % id)

    branchingData = {}
    for nuclide in nuclide.isotope:
        photons = []
        energy = nuclide.energy[0].pqu()
        for decayMode in nuclide.decayData.decayModes:
            if decayMode.mode == Mode.electroMagnetic:
                probability = decayMode.probability[0].value
                photonEmissionProbability = 1
                if len(decayMode.photonEmissionProbabilities) > 0:
                    photonEmissionProbability = decayMode.photonEmissionProbabilities[0].value
                decayPath = decayMode.decayPath[0]
                products = [ decayProduct.pid for decayProduct in decayPath.products ]
                products.remove(IDsModule.photon)
                residual = products[0]
                photons.append((probability, energy - branchingData[residual]['energy'], residual, photonEmissionProbability))
        branchingData[nuclide.id] = {'energy': energy, 'photons': photons}

    for nuclideId, nuclideData in branchingData.items():
        completePhotons = {}
        completePhotonBranchingData(nuclideId, branchingData, 1.0, completePhotons)
        nuclideData['completePhotons'] = completePhotons

    return branchingData
