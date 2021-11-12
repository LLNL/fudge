# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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
# isotropic2d
#
def toENDL( self ) :

    return( None )

angularModule.isotropic2d.toENDL = toENDL

#
# recoil
#
def toENDL( self ) :

    outputChannel = self.findClassInAncestry( outputChannelModule.outputChannel )
    product = outputChannel.products[0]
    distribution = product.distribution[0].invert( ).toENDL( )
    return( distribution )

angularModule.recoil.toENDL = toENDL

#
# twoBodyForm
#
def toENDL( self ) :

    return( { 1 : self.angularSubform.toENDL( ) } )

angularModule.twoBodyForm.toENDL = toENDL

#
# form
#
def toENDL( self ) :

    print( self.moniker )

angularModule.form.toENDL = toENDL
