# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import gridded as griddedModule

from fudge import abstractClasses as abstractClassesModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( 'energy_available', 0, energyUnit )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

class baseAvailableEnergyForm( abstractClassesModule.form ) :

    pass
#
# availableEnergy forms
#
class XYs1d( baseAvailableEnergyForm, XYsModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        XYsModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print('%sMulti-grouping XYs1d available energy' % indent)

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self ) )

class gridded1d( baseAvailableEnergyForm, griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded1d.__init__( self, **kwargs )
#
# availableEnergy component
#
class component( abstractClassesModule.component ) :

    moniker = 'availableEnergy'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( XYs1d, gridded1d, ) )

def parseXMLNode( element, xPath, linkData ) :
    """Reads an xml <availableEnergy> component element into fudge, including all forms in the component."""

    xPath.append( element.tag )

    aec = component()
    for form in element:
        formClass = {
                XYs1d.moniker : XYs1d,
                gridded1d.moniker : gridded1d,
            }.get( form.tag )
        if( formClass is None ) : raise Exception( "encountered unknown availableEnergy form: %s" % form.tag )
        try :
            newForm = formClass.parseXMLNode( form, xPath, linkData )
        except Exception as e:
            raise Exception( "availableEnergy/%s: %s" % (form.tag, e) )
        aec.add( newForm )
    xPath.pop()
    return aec
