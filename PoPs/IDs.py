# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Define some common particle ids
"""

from . import misc as miscModule

# gauge bosons
photon = 'photon'

# leptons
neutrinoPrefix = 'nu_'

electron = 'e-'
electronAnti = miscModule.returnAntiParticleIDFromId( electron )
positron = electronAnti

electronNeutrino = neutrinoPrefix + electron
electronAntiNeutrino = miscModule.returnAntiParticleIDFromId( electronNeutrino )

muon = 'mu-'
muonAnti = miscModule.returnAntiParticleIDFromId( muon )

muonNeutrino = neutrinoPrefix + muon
muonNeutrinoAnti = miscModule.returnAntiParticleIDFromId( muonNeutrino )

tauon = 'tau-'
tauonAnti = miscModule.returnAntiParticleIDFromId( tauon )

tauonNeutrino = neutrinoPrefix + tauon
tauonNeutrinoAnti = miscModule.returnAntiParticleIDFromId( tauonNeutrino )

# baryons
neutron = 'n'
neutronAnti = miscModule.returnAntiParticleIDFromId( neutron )

proton = 'p'
protonAnti = miscModule.returnAntiParticleIDFromId( proton )

# Familiar nuclear particle names.
familiarPhoton = 'g'
familiarProton = proton
familiarDeuteron = 'd'
familiarTriton = 't'
familiarHelion = 'h'
familiarAlpha = 'a'
