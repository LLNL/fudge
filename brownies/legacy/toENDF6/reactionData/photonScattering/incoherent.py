# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherent as incoherentModule

#
# ScatteringFactor.Form
#
def toENDF6( self, MT, endfMFList, targetInfo ) :

    self.data.toENDF6( MT, endfMFList, targetInfo )

incoherentModule.ScatteringFactor.toENDF6 = toENDF6

#
# incoherent.Form
#
def toENDF6( self, MT, endfMFList, targetInfo ) :

    self.scatteringFactor.toENDF6( MT, endfMFList, targetInfo )

incoherentModule.Form.toENDF6 = toENDF6
