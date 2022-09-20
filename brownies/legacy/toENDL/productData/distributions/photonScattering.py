# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.reactions import base as reactionBaseModule
from fudge.productData.distributions import photonScattering as photonScatteringModule

#
# CoherentPhotonScattering
#
def toENDL( self ) :

    reaction = self.findClassInAncestry( reactionBaseModule.Base_reaction )
    doubleDifferentialCrossSection = reaction.doubleDifferentialCrossSection[0]
    formFactor = doubleDifferentialCrossSection.formFactor.data
    data = [ formFactor[0][0] ]
    for xy in formFactor[1] : data.append( xy )
    return( { 941 : data } )

photonScatteringModule.CoherentPhotonScattering.Form.toENDL = toENDL

#
# IncoherentPhotonScattering
#
def toENDL( self ) :

    reaction = self.findClassInAncestry( reactionBaseModule.Base_reaction )
    incoherentPhotonScattering = reaction.doubleDifferentialCrossSection[0]
    formFactor = incoherentPhotonScattering.scatteringFactor
    data = [ formFactor[0][0] ]
    for xy in formFactor[1] : data.append( xy )
    return( { 942 : data } )

photonScatteringModule.IncoherentPhotonScattering.Form.toENDL = toENDL
