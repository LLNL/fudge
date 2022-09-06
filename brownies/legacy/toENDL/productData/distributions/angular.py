# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge import outputChannel as outputChannelModule
from fudge.productData.distributions import angular as angularModule

#
# XYs1d
#
def toENDL( self ) :

    return( [ self.outerDomainValue, [ [ x, y ] for x, y in self ] ] )

angularModule.XYs1d.toENDL = toENDL

#
# XYs2d
#
def toENDL( self ) :

    return( [ function1d.toENDL( ) for function1d in self ] )

angularModule.XYs2d.toENDL = toENDL

#
# Isotropic2d
#
def toENDL( self ) :

    return( None )

angularModule.Isotropic2d.toENDL = toENDL

#
# Recoil
#
def toENDL( self ) :

    outputChannel = self.findClassInAncestry( outputChannelModule.OutputChannel )
    product = outputChannel.products[0]
    distribution = product.distribution[0].invert( ).toENDL( )
    return( distribution )

angularModule.Recoil.toENDL = toENDL

#
# TwoBody
#
def toENDL( self ) :

    return( { 1 : self.angularSubform.toENDL( ) } )

angularModule.TwoBody.toENDL = toENDL

#
# form
#
def toENDL( self ) :

    print( self.moniker )

angularModule.Form.toENDL = toENDL
