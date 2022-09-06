# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import energy as energyModule

#
# XYs1d
#
def toENDL( self ) :

    return( [ self.outerDomainValue, [ [ x, y ] for x, y in self ] ] )

energyModule.XYs1d.toENDL = toENDL

#
# XYs2d
#
def toENDL( self ) :

    return( [ function1d.toENDL( ) for function1d in self ] )

energyModule.XYs2d.toENDL = toENDL
