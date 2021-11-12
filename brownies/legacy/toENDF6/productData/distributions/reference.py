# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import reference as referenceModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CPNElasticModule


def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    if isinstance(self.link, CPNElasticModule.form):
        self.link.toENDF6( MT, endfMFList, flags, targetInfo )
    else:
        pass

referenceModule.form.toENDF6 = toENDF6
