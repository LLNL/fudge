# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherent as incoherentModule

#
# incoherent.form
#
def toENDF6( self, MT, endfMFList, targetInfo ) :

    self.scatteringFunction.toENDF6( MT, endfMFList, targetInfo )

incoherentModule.form.toENDF6 = toENDF6
