# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import branching3d as branching3dModule

#
# form
#
def toENDL( self ) :

    return( None )

branching3dModule.Form.toENDL = toENDL
