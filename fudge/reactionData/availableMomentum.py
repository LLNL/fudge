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

def defaultAxes( energyUnit, momentumUnit ) : 

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( 'momentum_available', 0, momentumUnit )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

class baseAvailableMomentumForm( abstractClassesModule.form ) :

    pass
#
# availableMomentum forms
#
class XYs1d( baseAvailableMomentumForm, XYsModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        XYsModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print('%sMulti-grouping XYs1d available momentum' % indent)

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self ) )

class gridded1d( baseAvailableMomentumForm, griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded1d.__init__( self, **kwargs )
#
# availableMomentum forms
#
class component( abstractClassesModule.component ) :

    moniker = 'availableMomentum'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( XYs1d, gridded1d, ) )

def parseXMLNode( element, xPath, linkData ):
    """Reads an xml <availableMomentum> component element into fudge, including all forms in the component."""

    xPath.append( element.tag )
    amc = component()
    for form in element :
        formClass = {
            XYs1d.moniker : XYs1d,
            gridded1d.moniker : gridded1d
        }.get( form.tag )
        if( formClass is None ) : raise Exception( "encountered unknown availableMomentum form: %s" % form.tag )
        try :
            newForm = formClass.parseXMLNode( form, xPath, linkData )
        except Exception as e :
            raise Exception( "availableMomentum/%s: %s" % ( form.tag, e ) )
        amc.add( newForm )
    xPath.pop()
    return amc

def calculateMomentumPoints( style, massInE, EMin, EMax, energyUnit, accuracy = 1e-4 ) :
    """
    This is a temporary function (hopefully) to calculate momentum vs E.
    What is really needed is a function like rriint in mcfgen.
    """

    import math

    def insertPointIfNeeded( momentum, Ep1, Ep2, level ) :

        if( level == 16 ) : return

        E1, p1 = Ep1
        E2, p2 = Ep2
        s = ( p2 - p1 ) / ( E2 - E1 )
        Em = massInE / ( 2. * s * s )           # Point of largest difference to linear assumption between E1 and E2.
        pc = math.sqrt( 2 * massInE * Em )
        pi = s * ( Em - E1 ) + p1
        if( abs( pi / pc  - 1 ) > accuracy ) :
            m = [ Em, pc ]
            momentum.append( m )
            insertPointIfNeeded( momentum, Ep1, m, level + 1 )
            insertPointIfNeeded( momentum, m, Ep2, level + 1 )

    if( massInE == 0 ) :            # gammas
        momentum = [ [ EMin, EMin ], [ EMax, EMax ] ]
    else :
        momentum = [ [ EMin, math.sqrt( 2 * massInE * EMin ) ], [ EMax, math.sqrt( 2 * massInE * EMax ) ] ]
        insertPointIfNeeded( momentum, momentum[0], momentum[1], 0 )
    momentum.sort( )
    axes = defaultAxes( energyUnit, energyUnit + '/c' )
    return( XYs1d( data = momentum, axes = axes, label = style.label ) )
