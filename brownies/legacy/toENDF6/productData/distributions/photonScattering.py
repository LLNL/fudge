# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import photonScattering as photonScatteringModule

#
# coherent
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    pass

photonScatteringModule.coherentPhotonScattering.form.toENDF6 = toENDF6

#
# incoherent
#
photonScatteringModule.incoherentPhotonScattering.form.toENDF6 = toENDF6
