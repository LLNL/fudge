# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule
from xData import gridded as griddedModule

from fudge import abstractClasses as abstractClassesModule

def defaultAxes( energyUnit ) :

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( Component.moniker, 0, energyUnit )
    axes[1] = axesModule.Axis( 'energy_in', 1, energyUnit )
    return( axes )

class BaseEnergyDepositionForm( abstractClassesModule.Form ) :

    keyName = 'label'

class XYs1d( BaseEnergyDepositionForm, XYs1dModule.XYs1d ) :

    def __init__( self, **kwargs ) :

        BaseEnergyDepositionForm.__init__( self )
        XYs1dModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print('%sMulti-grouping XYs1d average energy' % indent)

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( Gridded1d, style, tempInfo, self ) )

class Regions1d( BaseEnergyDepositionForm, regionsModule.Regions1d ) :

    def __init__( self, **kwargs ) :

        BaseEnergyDepositionForm.__init__( self )
        regionsModule.Regions1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        linear = self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
        return( linear.processMultiGroup( style, tempInfo, indent ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class Gridded1d(BaseEnergyDepositionForm, griddedModule.Gridded1d):

    def __init__(self, axes, array, **kwargs):

        BaseEnergyDepositionForm.__init__(self)
        griddedModule.Gridded1d.__init__(self, axes, array, **kwargs)
#
# averageProductEnergy component
#
class Component( abstractClassesModule.Component ) :

    moniker = 'averageProductEnergy'

    def __init__( self ) :

        abstractClassesModule.Component.__init__( self, ( XYs1d, Regions1d, Gridded1d ) )
