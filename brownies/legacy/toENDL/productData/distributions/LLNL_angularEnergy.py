# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import LLNL_angularEnergy as LLNL_angularEnergyModule

#
# XYs1d
#
def toENDL( self ) :

    return( [ self.outerDomainValue, [ [ x, y ] for x, y in self ] ] )

LLNL_angularEnergyModule.XYs1d.toENDL = toENDL

#
# XYs2d
#
def toENDL( self ) :

    return( [ self.outerDomainValue, [ function1d.toENDL( ) for function1d in self ] ] )

LLNL_angularEnergyModule.XYs2d.toENDL = toENDL

#
# XYs3d
#
def toENDL( self ) :

    return( [ function2d.toENDL( ) for function2d in self ] )

LLNL_angularEnergyModule.XYs3d.toENDL = toENDL

#
# LLNLAngularOfAngularEnergySubform
#
def toENDL( self ) :

    return( self.data.toENDL( ) )

LLNL_angularEnergyModule.LLNLAngularOfAngularEnergySubform.toENDL = toENDL

#
# LLNLAngularEnergyOfAngularEnergySubform
#
def toENDL( self ) :

    return( self.data.toENDL( ) )

LLNL_angularEnergyModule.LLNLAngularEnergyOfAngularEnergySubform.toENDL = toENDL

#
# LLNLAngularEnergyForm
#
def toENDL( self ) :

    I1And3 = { 1 : self.angularSubform.toENDL( ) }
    I1And3[3] = self.angularEnergySubform.toENDL( )
    return( I1And3 )

LLNL_angularEnergyModule.LLNLAngularEnergyForm.toENDL = toENDL
