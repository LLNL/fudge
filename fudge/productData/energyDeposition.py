# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.regions as regionsModule
from xData import gridded as griddedModule

from fudge import abstractClasses as abstractClassesModule

__metaclass__ = type

averageProductEnergyToken = 'averageProductEnergy'

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( averageProductEnergyToken, 0, energyUnit )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

class baseEnergyDepositionForm( abstractClassesModule.form ) :

    pass

class XYs1d( baseEnergyDepositionForm, XYsModule.XYs1d ) :

    def __init__( self, **kwargs ) :

        baseEnergyDepositionForm.__init__( self )
        XYsModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print('%sMulti-grouping XYs1d average energy' % indent)

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self ) )

class regions1d( baseEnergyDepositionForm, regionsModule.regions1d ) :

    def __init__( self, **kwargs ) :

        baseEnergyDepositionForm.__init__( self )
        regionsModule.regions1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        linear = self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
        return( linear.processMultiGroup( style, tempInfo, indent ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class gridded1d( baseEnergyDepositionForm, griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded1d.__init__( self, **kwargs )
#
# energyDeposition component
#
class component( abstractClassesModule.component ) :

    moniker = averageProductEnergyToken

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( XYs1d, regions1d, gridded1d ) )
