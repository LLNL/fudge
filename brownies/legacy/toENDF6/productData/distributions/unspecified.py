# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import unspecified as unspecifiedModule

from ... import gndsToENDF6 as gndsToENDF6Module

def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, 0, self.productFrame, [] )

unspecifiedModule.form.toENDF6 = toENDF6
