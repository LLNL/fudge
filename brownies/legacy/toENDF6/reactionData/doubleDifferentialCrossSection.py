# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import fudge.reactionData.doubleDifferentialCrossSection.doubleDifferentialCrossSection as doubleDifferentialCrossSectionModule
import fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic as CPElasticModule

def toENDF6( self, MT, endfMFList, targetInfo ) :
    """
    Convert self into ENDF format

    :param int MT: The ENDF reaction designator, MT
    :param endfMFList:
    :param targetInfo:
    """

    if( targetInfo['style'] in self ) :

        form = self[targetInfo['style']]
        if isinstance(form, CPElasticModule.CoulombPlusNuclearElastic.Form):
            return      # this gets translated when handling outgoing products
        form.toENDF6( MT, endfMFList, targetInfo )

doubleDifferentialCrossSectionModule.Component.toENDF6 = toENDF6
