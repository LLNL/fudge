# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import reference as referenceModule


#
# form
#
def toENDL( self ) :

    return( self.link.toENDL( ) )

referenceModule.Form.toENDL = toENDL
