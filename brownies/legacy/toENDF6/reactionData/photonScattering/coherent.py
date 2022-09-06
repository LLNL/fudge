# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import coherent as coherentModule


#
# coherent.Form
#
def toENDF6( self, MT, endfMFList, targetInfo ) :

    if( self.formFactor is not None ) : self.formFactor.data.toENDF6( MT, endfMFList, targetInfo )
    if( self.anomalousScatteringFactor_realPart is not None ) : self.anomalousScatteringFactor_realPart.data.toENDF6( MT, endfMFList, targetInfo )
    if( self.anomalousScatteringFactor_imaginaryPart is not None ) : self.anomalousScatteringFactor_imaginaryPart.data.toENDF6( MT, endfMFList, targetInfo )

coherentModule.Form.toENDF6 = toENDF6
