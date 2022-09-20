# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import unspecified as unspecifiedModule

#
# form
#
def toENDL( self ) :

    return( None )

unspecifiedModule.Form.toENDL = toENDL
