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


def defaultAxes(energyUnit, momentumUnit=None):

    if momentumUnit is None:
        momentumUnit = energyUnit + '/c'

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( 'momentum_available', 0, momentumUnit )
    axes[1] = axesModule.Axis( 'energy_in', 1, energyUnit )
    return( axes )

class BaseAvailableMomentumForm( abstractClassesModule.Form ) :

    pass
#
# availableMomentum forms
#
class XYs1d( BaseAvailableMomentumForm, XYs1dModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        XYs1dModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print('%sMulti-grouping XYs1d available momentum' % indent)

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( Gridded1d, style, tempInfo, self ) )

class Gridded1d(BaseAvailableMomentumForm, griddedModule.Gridded1d):

    def __init__(self, axes, array, **kwargs):

        BaseAvailableMomentumForm.__init__(self)
        griddedModule.Gridded1d.__init__(self, axes, array, **kwargs)
#
# availableMomentum component.
#
class Component( abstractClassesModule.Component ) :

    moniker = 'availableMomentum'

    def __init__( self ) :

        abstractClassesModule.Component.__init__( self, ( XYs1d, Gridded1d, ) )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Reads an xml <availableMomentum> component element into fudge, including all forms in the component."""

        xPath.append( element.tag )

        amc = cls()
        for form in element :
            formClass = {
                XYs1d.moniker : XYs1d,
                Gridded1d.moniker : Gridded1d
            }.get( form.tag )
            if( formClass is None ) : raise Exception( "encountered unknown availableMomentum form: %s" % form.tag )
            try :
                newForm = formClass.parseNodeUsingClass(form, xPath, linkData, **kwargs)
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
