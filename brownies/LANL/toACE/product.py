# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the classes in the fudge.product module.
"""

from PoPs import IDs as IDsPoPsModule

from fudge import product as productModule

from fudge.productData.distributions import unspecified as unspecifiedModule, energyAngularMC as energyAngularMCModule
from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.productData.distributions import angularEnergyMC as angularEnergyMCModule
from fudge.productData.distributions import KalbachMann as KalbachMannModule


def toACE( self, cdf_style, MTData, MT, verbose ) :

    if( verbose > 2 ) : print('        %s: label = %s' % (self.id, self.label))
    if( self.id in [ IDsPoPsModule.neutron, IDsPoPsModule.photon ] ) :
        if( self.multiplicity.isConstant( ) ) :
            multiplicity = self.multiplicity.evaluate( 0 )
        else :
            multiplicity = self.multiplicity.evaluated

    if( self.id == IDsPoPsModule.neutron ) :
        angularData = None
        energyData = None
        evaluated = self.distribution.evaluated

        if( ( MT == 18 ) and isinstance( evaluated, unspecifiedModule.form ) ) :
            distribution = evaluated
        else :
            try :
                distribution = self.distribution[cdf_style.label]
            except :
                distribution = evaluated

            if( isinstance( distribution, angularModule.twoBodyForm ) ) :
                angularData = distribution.angularSubform
            elif( isinstance( distribution, uncorrelatedModule.form ) ) :
                angularData = distribution.angularSubform.data
                energyData = distribution.energySubform.data
            elif( isinstance( distribution, energyAngularMCModule.form ) ) :
                energyData = distribution
            elif( isinstance( distribution, angularEnergyMCModule.form ) ) :
                angularData = distribution.angular.data
                energyData = distribution.angularEnergy.data
            elif( isinstance( distribution, KalbachMannModule.form ) ) :
                energyData = distribution
            else :
                raise Exception( 'Unsupported distribution form = %s' % distribution.moniker )

        frame = distribution.productFrame
        MTData[IDsPoPsModule.neutron].append( { 'product' : self, 'frame' : frame, 'multiplicity' : multiplicity, 
            'angularData' : angularData, 'energyData' : energyData } )

productModule.product.toACE = toACE
