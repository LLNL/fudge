# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import gridded as griddedModule

from fudge import abstractClasses as abstractClassesModule


def defaultAxes( energyUnit ) :

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( 'energy_available', 0, energyUnit )
    axes[1] = axesModule.Axis( 'energy_in', 1, energyUnit )
    return( axes )

class BaseAvailableEnergyForm( abstractClassesModule.Form ) :

    pass
#
# availableEnergy forms
#
class XYs1d( BaseAvailableEnergyForm, XYs1dModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        XYs1dModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print('%sMulti-grouping XYs1d available energy' % indent)

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( Gridded1d, style, tempInfo, self ) )

class Gridded1d(BaseAvailableEnergyForm, griddedModule.Gridded1d):

    def __init__(self, axes, array, **kwargs ):

        BaseAvailableEnergyForm.__init__(self)
        griddedModule.Gridded1d.__init__(self, axes, array, **kwargs)
#
# availableEnergy component
#
class Component( abstractClassesModule.Component ) :

    moniker = 'availableEnergy'

    def __init__( self ) :

        abstractClassesModule.Component.__init__( self, ( XYs1d, Gridded1d, ) )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Reads an xml <availableEnergy> component element into fudge, including all forms in the component."""

        xPath.append( element.tag )

        aec = cls()
        for form in element:
            formClass = {
                    XYs1d.moniker : XYs1d,
                    Gridded1d.moniker : Gridded1d,
                }.get( form.tag )
            if( formClass is None ) : raise Exception( "encountered unknown availableEnergy form: %s" % form.tag )
            try :
                newForm = formClass.parseNodeUsingClass(form, xPath, linkData, **kwargs)
            except Exception as e:
                raise Exception( "availableEnergy/%s: %s" % (form.tag, e) )
            aec.add( newForm )

        xPath.pop()

        return aec
