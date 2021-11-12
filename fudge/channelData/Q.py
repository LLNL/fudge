# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import axes as axesModule
from xData import constant as constantModule
from xData import XYs as XYsModule
from xData import gridded as griddedModule

from fudge import abstractClasses as abstractClassesModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( )
# BRB6 make 'energy_in' a token
    axes[1] = axesModule.axis( 'energy_in', 0, energyUnit )
    axes[0] = axesModule.axis( component.moniker, 1, energyUnit )
    return( axes )

#
# Q forms
#
class baseQForm( abstractClassesModule.form ) :

    def diff( self, other, diffResults ) :

        if( isinstance( other, constant1d ) ) :
            diffResults.append( '%s non-constant Q - 1' % other.moniker, '', self.toXLink( ), other.toXLink( ) )

class constant1d( baseQForm, constantModule.constant1d ) :

    def __init__( self, Q, domainMin, domainMax, axes, label = None ) :

        baseQForm.__init__( self )
        constantModule.constant1d.__init__( self, Q, domainMin, domainMax, axes = axes, label = label )

    def diff( self, other, diffResults ) :

        if( not( isinstance( other, constant1d ) ) ) :
            diffResults.append( '%s non-constant Q - 2' % other.moniker, '', self.toXLink( ), other.toXLink( ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """This method returns the Q-value as linear-linear XYs1d data which spans self's domain."""

        return( XYs1d( data = [ [ self.domainMin, self.value ], [ self.domainMax, self.value ] ], axes = self.axes ) )

class XYs1d( baseQForm, XYsModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        baseQForm.__init__( self )
        XYsModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self ) )

class gridded1d( baseQForm, griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded1d.__init__( self, **kwargs )
#
# Q component
#
class component( abstractClassesModule.component ) :

    moniker = 'Q'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( constant1d, XYs1d, gridded1d ) )

    def thresholdQAs( self, unit ) :

        for form in self :
            if( isinstance( form, constant1d ) ) :
                return( PQUModule.PQU( form.value, form.axes[0].unit ).getValueAs( unit ) )
        raise ValueError( 'Q-value is not a constant' )

    def getConstantAs( self, unit ) :

        return( self.thresholdQAs( unit ) )

    def getConstant( self ):

        for form in self :
            if( isinstance( form, constant1d ) ) :
                return( form.value )
            elif( isinstance( form, XYs1d ) ) :
                return( form.evaluate( form.domainMin ) )

        raise ValueError( 'Q-value is not a constant' )

    def evaluate( self, E ) :

        return( self.getEvaluated( ).evaluate( E ) )

def parseXMLNode( QElement, xPath, linkData ):
    """
    Reads an xml <Q> element into fudge, including all child forms
    """

    xPath.append( QElement.tag )
    kwargs = {}     # In case we need special interpolation rules for Q, see gnds/reactionData/crossSection.py

    Q = component( )
    for form in QElement :
        formClass = {
                constant1d.moniker     : constant1d,
                XYs1d.moniker          : XYs1d,
                gridded1d.moniker      : gridded1d      }.get( form.tag )
        if( formClass is None ) : raise Exception( "unknown Q form: %s" % form.tag )
        newForm = formClass.parseXMLNode( form, xPath, linkData, **kwargs )
        Q.add( newForm )
    xPath.pop( )
    return( Q )
