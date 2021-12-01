# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains a special type of Q-value, unique to fission reactions.
Fission releases energy several different ways (neutrons, gammas, etc.), and it's useful to subdivide the Q-value
into these different terms.
"""

# FIXME this is really a form of Q-value. Should it be defined in Q.py instead?

import xData.ancestry as ancestryModule
import xData.standards as standardsModule
import xData.axes as axesModule
import xData.series1d as series1dModule
import xData.XYs as XYsModule
import xData.gridded as griddedModule
import xData.uncertainties as uncertaintiesModule

from fudge import abstractClasses as abstractClassesModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( 'energy_out', 0, energyUnit )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

class polynomial1d( series1dModule.polynomial1d ):

    def __init__( self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None,
            index = None, valueType = standardsModule.types.float64Token, outerDomainValue = None, label = None, sep = ' ',
            coefficientUncertainties = None ):

        series1dModule.polynomial1d.__init__( self, coefficients = coefficients, domainMin = domainMin, domainMax = domainMax,
                    lowerIndex = lowerIndex, axes = axes, index = index, valueType = valueType, outerDomainValue = outerDomainValue, label = label, sep = sep )
        if coefficientUncertainties is not None:
            self.uncertainty = uncertaintiesModule.uncertainty( functional = series1dModule.polynomial1d( coefficientUncertainties, domainMin, domainMax, axes=axes ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class XYs1d( XYsModule.XYs1d ) :

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        multiGroup = miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self )
        multiGroup.label = None
        return( multiGroup )

class gridded1d( griddedModule.gridded1d ) :

    pass

class fissionEnergyReleaseTerm( ancestryModule.ancestry ):
    """ Base class for all types of fission energy release. """

    ancestryMembers = ( 'data', )

    def __init__( self, data ) :

        ancestryModule.ancestry.__init__( self )
        self.data = data

    @property
    def data(self):
        return self.__data

    @data.setter
    def data( self, value ) :

        if not isinstance( value, ( XYs1d, polynomial1d, gridded1d ) ) :
            raise TypeError( 'Invalid class "%s" for fissionEnergyReleaseTerm' % type(value) )
        value.setAncestor(self)
        self.__data = value

    def check(self, info):
        from fudge import warning

        warnings = []
        linearized = self.data.toPointwise_withLinearXYs(accuracy=1e-6)
        if linearized.rangeMin < 0:
            warnings.append( warning.badFissionEnergyRelease( worstCase=linearized.rangeMin, obj=self) )
        if not 0.1 < linearized.rangeMax / linearized.rangeMin < 10:
            warnings.append( warning.badFissionEnergyRelease( worstCase=linearized.rangeMax, obj=self) )
        return warnings

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.__class__( self.data.processMultiGroup( style, tempInfo, indent ) ) )

    def toXMLList( self, indent="", **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xmllist = ['%s<%s>' % (indent, self.moniker)]
        xmllist += self.data.toXMLList( indent2, **kwargs )
        xmllist[-1] += '</%s>' % self.moniker
        return xmllist

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        dataClass = {
            XYs1d.moniker: XYs1d,
            polynomial1d.moniker: polynomial1d,
            gridded1d.moniker: gridded1d
        }.get( element[0].tag, None )
        data = dataClass.parseXMLNode( element[0], xPath, linkData )
        FERT = cls( data )
        xPath.pop()
        return FERT

class promptProductKE( fissionEnergyReleaseTerm ):

    moniker = 'promptProductKE'

class promptNeutronKE( fissionEnergyReleaseTerm ):

    moniker = 'promptNeutronKE'

class delayedNeutronKE( fissionEnergyReleaseTerm ):

    moniker = 'delayedNeutronKE'

class promptGammaEnergy( fissionEnergyReleaseTerm ):

    moniker = 'promptGammaEnergy'

class delayedGammaEnergy( fissionEnergyReleaseTerm ):

    moniker = 'delayedGammaEnergy'

class delayedBetaEnergy( fissionEnergyReleaseTerm ):

    moniker = 'delayedBetaEnergy'

class neutrinoEnergy( fissionEnergyReleaseTerm ):

    moniker = 'neutrinoEnergy'

class nonNeutrinoEnergy( fissionEnergyReleaseTerm ):

    moniker = 'nonNeutrinoEnergy'

class totalEnergy( fissionEnergyReleaseTerm ):

    moniker = 'totalEnergy'

class field:
    """
    Descriptor to ensure ancestry is set when adding energy release terms
    to fissionEnergyRelease class.
    """

    def __init__(self, Class):
        self.Class = Class
        self.fname = '__' + Class.moniker

    def __get__(self, instance, owner):
        return getattr( instance, self.fname, None )

    def __set__(self, instance, value):
        if not isinstance( value, self.Class ):
            raise TypeError( "Incorrect type: expected %s, got %s" % (self.Class.moniker, type(value)) )
        value.setAncestor(instance)
        setattr( instance, self.fname, value )

class fissionEnergyRelease( ancestryModule.ancestry ) :
    """
    Store average energy released to different types of fission products.
    (prompt and delayed neutrons, prompt / delayed gammas, betas, neutrinos, etc.)
    Each term is currently (when translating from ENDF) stored as a polynomial expansion,
    although we expect to also see XYs1d representations in future evaluations
    """

    moniker = 'fissionEnergyRelease'

    ancestryMembers = ( 'promptProductKE',      'promptNeutronKE',      'delayedNeutronKE', 'promptGammaEnergy',
                        'delayedGammaEnergy',   'delayedBetaEnergy',    'neutrinoEnergy',   'nonNeutrinoEnergy',
                        'totalEnergy' )

    promptProductKE = field( promptProductKE )
    promptNeutronKE = field( promptNeutronKE )
    delayedNeutronKE = field( delayedNeutronKE )
    promptGammaEnergy = field( promptGammaEnergy )
    delayedGammaEnergy = field( delayedGammaEnergy )
    delayedBetaEnergy = field( delayedBetaEnergy )
    neutrinoEnergy = field( neutrinoEnergy )
    nonNeutrinoEnergy = field( nonNeutrinoEnergy )
    totalEnergy = field( totalEnergy )

    def __init__( self, label, **kwargs ) :

        ancestryModule.ancestry.__init__( self )
        self.label = label

        for key in list( kwargs.keys( ) ) :
            if( key in self.ancestryMembers ) : setattr( self, key, kwargs.pop( key, None ) )
        if( kwargs ) : raise TypeError( "fissionEnergyRelease received unexpected argument(s): '%s'" % list( kwargs.keys( ) ) )

    def __iter__( self ) :

        for term in self.ancestryMembers : yield getattr( self, term )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        for term in self.ancestryMembers :
            fwarnings = getattr( self, term ).convertUnits( unitMap )

    def check( self, info ) :

        from fudge import warning
        warnings = []
        for term in self.ancestryMembers:
            fwarnings = getattr(self, term).check( info )
            if fwarnings:
                warnings.append( warning.context('%s:' % term, fwarnings) )
        return warnings

    def processMultiGroup( self, style, tempInfo, indent ) :

        kwargs = {}
        for term in self.ancestryMembers :
            kwargs[term] = getattr( self, term ).processMultiGroup( style, tempInfo, indent )
        return( fissionEnergyRelease( style.label, **kwargs ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.promptProductKE.data.toPointwise_withLinearXYs( **kwargs ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent="", **kwargs ):

        indent2 = indent+"  "
        xmlList = ['%s<%s label="%s">' % ( indent, self.moniker, self.label )]
        for term in self.ancestryMembers:
            xmlList += getattr(self, term).toXMLList( indent2, **kwargs )
        xmlList[-1] += "</%s>" % self.moniker

        return xmlList

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        """Parse <fissionEnergyRelease> node from xml."""

        xPath.append( element.tag )

        children = {}
        childClasses = dict( [ ( cls.moniker, cls ) for cls in ( promptProductKE, promptNeutronKE, delayedNeutronKE, promptGammaEnergy,
                delayedGammaEnergy, delayedBetaEnergy, neutrinoEnergy, nonNeutrinoEnergy, totalEnergy ) ] )
        for child in element :
            children[ child.tag ] = childClasses[ child.tag ].parseXMLNode( child, xPath, linkData )
        fer = fissionEnergyRelease( label = element.get( "label" ), **children )

        xPath.pop( )
        return fer

class fissionEnergyReleases( abstractClassesModule.component ) :

    moniker = 'fissionEnergyReleases'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( fissionEnergyRelease, ) )
