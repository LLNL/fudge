# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import series1d as series1dModule
from xData import XYs1d as XYs1dModule
from xData import gridded as griddedModule
from xData import uncertainties as uncertaintiesModule
from xData import vector as vectorModule

from fudge import abstractClasses as abstractClassesModule
from fudge.processing import transporting as transportingModule


def defaultAxes( energyUnit ) :

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( 'energy_out', 0, energyUnit )
    axes[1] = axesModule.Axis( 'energy_in', 1, energyUnit )
    return( axes )

class Polynomial1d( series1dModule.Polynomial1d ):

    def __init__(self, coefficients, domainMin, domainMax, lowerIndex=0, axes=None,
            index=None, valueType=xDataEnumsModule.ValueType.float64, outerDomainValue=None, label=None, coefficientUncertainties=None):

        series1dModule.Polynomial1d.__init__( self, coefficients = coefficients, domainMin = domainMin, domainMax = domainMax,
                    lowerIndex = lowerIndex, axes = axes, index = index, valueType = valueType, outerDomainValue = outerDomainValue, label = label )
        if coefficientUncertainties is not None:
            self.uncertainty = uncertaintiesModule.Uncertainty( functional = series1dModule.Polynomial1d( coefficientUncertainties, domainMin, domainMax, axes=axes ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class XYs1d( XYs1dModule.XYs1d ) :

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        multiGroup = miscellaneousModule.groupFunctionCrossSectionAndFlux( Gridded1d, style, tempInfo, self )
        multiGroup.label = None
        return( multiGroup )

class Gridded1d( griddedModule.Gridded1d ) :

    pass

class FissionEnergyReleaseTerm( ancestryModule.AncestryIO ):
    """ Base class for all types of fission energy release. """

    ancestryMembers = ( 'data', )

    def __init__( self, data ) :

        ancestryModule.AncestryIO.__init__( self )
        self.data = data

    @property
    def data(self):
        return self.__data

    @data.setter
    def data( self, value ) :

        if not isinstance( value, ( XYs1d, Polynomial1d, Gridded1d ) ) :
            raise TypeError( 'Invalid class "%s" for fissionEnergyReleaseTerm' % type(value) )
        value.setAncestor(self)
        self.__data = value

    def check(self, info):
        from fudge import warning

        warnings = []
        linearized = self.data.toPointwise_withLinearXYs(accuracy=1e-6)
        if linearized.rangeMin < 0:
            warnings.append( warning.BadFissionEnergyRelease( worstCase=linearized.rangeMin, obj=self) )
        if not 0.1 < linearized.rangeMax / linearized.rangeMin < 10:
            warnings.append( warning.BadFissionEnergyRelease( worstCase=linearized.rangeMax, obj=self) )
        return warnings

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.__class__( self.data.processMultiGroup( style, tempInfo, indent ) ) )

    def toXML_strList( self, indent="", **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xmllist = ['%s<%s>' % (indent, self.moniker)]
        xmllist += self.data.toXML_strList( indent2, **kwargs )
        xmllist[-1] += '</%s>' % self.moniker
        return xmllist

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        dataClass = {
            XYs1d.moniker: XYs1d,
            Polynomial1d.moniker: Polynomial1d,
            Gridded1d.moniker: Gridded1d
        }.get( element[0].tag, None )
        data = dataClass.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        FERT = cls( data )
        xPath.pop()
        return FERT

class PromptProductKE( FissionEnergyReleaseTerm ):

    moniker = 'promptProductKE'

class PromptNeutronKE( FissionEnergyReleaseTerm ):

    moniker = 'promptNeutronKE'

class DelayedNeutronKE( FissionEnergyReleaseTerm ):

    moniker = 'delayedNeutronKE'

class PromptGammaEnergy( FissionEnergyReleaseTerm ):

    moniker = 'promptGammaEnergy'

class DelayedGammaEnergy( FissionEnergyReleaseTerm ):

    moniker = 'delayedGammaEnergy'

class DelayedBetaEnergy( FissionEnergyReleaseTerm ):

    moniker = 'delayedBetaEnergy'

class NeutrinoEnergy( FissionEnergyReleaseTerm ):

    moniker = 'neutrinoEnergy'

class NonNeutrinoEnergy( FissionEnergyReleaseTerm ):

    moniker = 'nonNeutrinoEnergy'

class TotalEnergy( FissionEnergyReleaseTerm ):

    moniker = 'totalEnergy'

class Field:
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

class FissionEnergyRelease( ancestryModule.AncestryIO ) :
    """
    Store average energy released to different types of fission products.
    (prompt and delayed neutrons, prompt / delayed gammas, betas, neutrinos, etc.)
    Each term is currently (when translating from ENDF) stored as a polynomial expansion,
    although we expect to also see XYs1d representations in future evaluations
    """

    moniker = 'fissionEnergyRelease'
    keyName = 'label'

    ancestryMembers = ( 'promptProductKE',      'promptNeutronKE',      'delayedNeutronKE', 'promptGammaEnergy',
                        'delayedGammaEnergy',   'delayedBetaEnergy',    'neutrinoEnergy',   'nonNeutrinoEnergy',
                        'totalEnergy' )

    promptProductKE = Field( PromptProductKE )
    promptNeutronKE = Field( PromptNeutronKE )
    delayedNeutronKE = Field( DelayedNeutronKE )
    promptGammaEnergy = Field( PromptGammaEnergy )
    delayedGammaEnergy = Field( DelayedGammaEnergy )
    delayedBetaEnergy = Field( DelayedBetaEnergy )
    neutrinoEnergy = Field( NeutrinoEnergy )
    nonNeutrinoEnergy = Field( NonNeutrinoEnergy )
    totalEnergy = Field( TotalEnergy )

    def __init__( self, label, **kwargs ) :

        ancestryModule.AncestryIO.__init__( self )
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
                warnings.append( warning.Context('%s:' % term, fwarnings) )
        return warnings

    def processMultiGroup( self, style, tempInfo, indent ) :

        kwargs = {}
        for term in self.ancestryMembers :
            kwargs[term] = getattr( self, term ).processMultiGroup( style, tempInfo, indent )
        return( FissionEnergyRelease( style.label, **kwargs ) )

    def multiGroupQ(self, multiGroupSettings):
        """
        Returns the multi-group Q-value.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        """
        
        Q = vectorModule.Vector()
        if multiGroupSettings.delayedNeutrons == transportingModule.DelayedNeutrons.on:
            Q += self.delayedNeutronKE.data.constructVector()
            Q += self.delayedGammaEnergy.data.constructVector()
            Q += self.delayedBetaEnergy.data.constructVector()

        return Q

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.promptProductKE.data.toPointwise_withLinearXYs( **kwargs ) )

    def toXML_strList( self, indent="", **kwargs ):

        indent2 = indent+"  "
        xmlList = ['%s<%s label="%s">' % ( indent, self.moniker, self.label )]
        for term in self.ancestryMembers:
            xmlList += getattr(self, term).toXML_strList( indent2, **kwargs )
        xmlList[-1] += "</%s>" % self.moniker

        return xmlList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Parse <fissionEnergyRelease> node from xml."""

        xPath.append( element.tag )

        children = {}
        childClasses = dict( [ ( cls.moniker, cls ) for cls in ( PromptProductKE, PromptNeutronKE, DelayedNeutronKE, PromptGammaEnergy,
                DelayedGammaEnergy, DelayedBetaEnergy, NeutrinoEnergy, NonNeutrinoEnergy, TotalEnergy ) ] )
        for child in element :
            children[ child.tag ] = childClasses[ child.tag ].parseNodeUsingClass(child, xPath, linkData, **kwargs)

        fer = cls( label = element.get( "label" ), **children )

        xPath.pop( )
        return fer

class FissionEnergyReleases( abstractClassesModule.Component ) :

    moniker = 'fissionEnergyReleases'

    def __init__( self ) :

        abstractClassesModule.Component.__init__( self, ( FissionEnergyRelease, ) )
