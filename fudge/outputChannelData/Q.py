# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import axes as axesModule
from xData import constant as constantModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule
from xData import series1d as series1dModule
from xData import gridded as griddedModule

from fudge import abstractClasses as abstractClassesModule

def defaultAxes( energyUnit ) :

    axes = axesModule.Axes(2)
# BRB6 make 'energy_in' a token
    axes[1] = axesModule.Axis( 'energy_in', 0, energyUnit )
    axes[0] = axesModule.Axis( Component.moniker, 1, energyUnit )
    return( axes )

#
# Q forms
#
class BaseQForm( abstractClassesModule.Form ) :

    keyName = 'label'

    def diff( self, other, diffResults ) :

        if( isinstance( other, Constant1d ) ) :
            diffResults.append( '%s non-constant Q - 1' % other.moniker, '', self.toXLink( ), other.toXLink( ) )

class Constant1d( BaseQForm, constantModule.Constant1d ) :

    def __init__(self, Q, domainMin, domainMax, axes=None, label=None, index=None, outerDomainValue=None):

        BaseQForm.__init__( self )
        constantModule.Constant1d.__init__(self, Q, domainMin, domainMax, axes=axes, label=label)

    def diff( self, other, diffResults ) :

        if( not( isinstance( other, Constant1d ) ) ) :
            diffResults.append( '%s non-constant Q - 2' % other.moniker, '', self.toXLink( ), other.toXLink( ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """This method returns the Q-value as linear-linear XYs1d data which spans self's domain."""

        return( XYs1d( data = [ [ self.domainMin, self.value ], [ self.domainMax, self.value ] ], axes = self.axes ) )

class XYs1d( BaseQForm, XYs1dModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        BaseQForm.__init__( self )
        XYs1dModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( Gridded1d, style, tempInfo, self ) )

class Regions1d(BaseQForm, regionsModule.Regions1d):

    def __init__(self, **kwargs):

        BaseQForm.__init__(self)
        regionsModule.Regions1d.__init__(self, **kwargs)

    def processMultiGroup(self, style, tempInfo, indent):

        raise NotImplementedError('Multi-group processing has not been implemented, please contact FUDGE developers.')

class Polynomial1d(BaseQForm, series1dModule.Polynomial1d):

    def __init__(self, **kwargs):

        BaseQForm.__init__(self)
        series1dModule.Polynomial1d.__init__(self, **kwargs)

    def processMultiGroup(self, style, tempInfo, indent):

        raise NotImplementedError('Multi-group processing has not been implemented, please contact FUDGE developers.')

class Gridded1d(BaseQForm, griddedModule.Gridded1d):

    def __init__(self, axes, array, **kwargs):

        BaseQForm.__init__(self)
        griddedModule.Gridded1d.__init__(self, axes, array, **kwargs)
#
# Q Component
#
class Component( abstractClassesModule.Component ) :

    moniker = 'Q'

    def __init__(self):

        abstractClassesModule.Component.__init__(self, (Constant1d, XYs1d, Regions1d, Polynomial1d, Gridded1d))

    def thresholdQAs( self, unit ) :

        for form in self :
            if( isinstance( form, Constant1d ) ) :
                return( PQUModule.PQU( form.value, form.axes[0].unit ).getValueAs( unit ) )
        raise ValueError( 'Q-value is not a constant' )

    def getConstantAs( self, unit ) :

        return( self.thresholdQAs( unit ) )

    def getConstant( self ):

        for form in self :
            if( isinstance( form, Constant1d ) ) :
                return( form.value )
            elif( isinstance( form, XYs1d ) ) :
                return( form.evaluate( form.domainMin ) )

        raise ValueError( 'Q-value is not a constant' )

    def evaluate( self, E ) :

        return( self.getEvaluated( ).evaluate( E ) )
