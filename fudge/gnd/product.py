# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""
This module contains the product class.

GND classes store reaction data in a hierarchical format. At the base is the product class. Next is the
channel that represents input and output channels. An input channel is a list of the inputcoming particles. 
As example, for a neutron hitting oxygen-16, the input channel is composed of the neutron and oxygen-16 
(i.e., 'n_1 + O_16'). An output channel is a list of the outgoing particles (some of which may decay into 
other particles). The third class in the hierarchy is the reaction class which is used to specify a 
reaction. A reaction is an input and output channel together, (e.g. n_1 + O_16  --> H_1 + N_16). Finally, the 
reactionSuite class is used to store a list of reaction for a given input channel. The following outline 
summarizes the classes.

:product:
    A particle distribution information.

:channel:
    A list of particles of class product.

:reaction:
    An input channel of class channel (e.g. n + Pu_239) and
    a single outgoing channel of sub-class of channel (e.g., 2n + Pu_239).
    Sub-classes of channel are twoBodyOutputChannel, NBodyOutputChannel and productionChannel.
    
:reactionSuite:
    An input channel of class channel (e.g. n + Pu_239)
    and a list of outgoing channels of class reaction with the same input channel instance::
    
        (e.g., n + Pu_239 --> n + Pu_239, n + Pu_239 --> n + Pu_239_e1, n + Pu_239 --> n + Pu_239_e2, n + Pu_239 --> 2n + Pu_239 ... )
"""
import string

from fudge.core.utilities import brb

import xData.ancestry as ancestryModule
from xData import standards as standardsModule

from fudge.gnd import suites as suitesModule

from .productData import distributions as distributionsModule
from .productData.distributions import distribution as distributionModule
from .productData.distributions import angular as angularModule
from .productData.distributions import energy as energyModule
from .productData.distributions import unspecified as unspecifiedModule
from .productData.distributions import reference as referenceModule
from .productData.distributions import KalbachMann as KalbachMannModule

from .productData import multiplicity as multiplicityModule
from .productData import energyDeposition as energyDepositionModule
from .productData import momentumDeposition as momentumDepositionModule
from . import channels as channelsModule
from . import productData

__metaclass__ = type

class product( ancestryModule.ancestry ) :
    """
    This is the class for a gnd particle. A gnd particle can decay (i.e., breakup), the decay
    formula is defined via the outputChannel member.
    """

    moniker = 'product'

    def __init__( self, particle, label = None, attributes = {}, outputChannel = None ) :
        """Creates a new product object."""

        ancestryModule.ancestry.__init__( self )
        self.particle = particle
        self.__label = label
        self.attributes = {}
        for q in attributes : self.addAttribute( q, attributes[q] )
        self.outputChannel = None
        if( outputChannel is not None ) : self.addOutputChannel( outputChannel )

        self.multiplicity = multiplicityModule.component( )
        self.multiplicity.setAncestor( self )

        self.energyDeposition = energyDepositionModule.component( )
        self.energyDeposition.setAncestor( self )

        self.momentumDeposition = momentumDepositionModule.component( )
        self.momentumDeposition.setAncestor( self )

        self.distribution = distributionModule.component( )
        self.distribution.setAncestor( self )

    def __cmp__( self, other ) :
        """Compares self to other."""

        if( self.name < other.name ) : return( -1 )
        if( self.name > other.name ) : return(  1 )
        if( self.outputChannel < other.outputChannel ) : return( -1 )
        if( self.outputChannel > other.outputChannel ) : return(  1 )
        return( 0 )

    def  __str__( self ) :
        """Converts product object to a string representation."""

        return( self.toString( simpleString = False ) )

    @property
    def name( self ) :
        """Returns self's name"""

        return( self.particle.name )

    @property
    def label( self ) :
        """Returns self's label."""

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = value


    def addOutputChannel( self, outputChannel ) :
        """Adds outputChannel to particle."""

        if( isinstance( outputChannel, channelsModule.channel ) ) :
            outputChannel.setAncestor( self )
            self.outputChannel = outputChannel
        else :
            raise Exception( 'Invalid decay channel = %s' % brb.getType( outputChannel ) )

    def addAttribute( self, name, value ) :
        """Add name and value to attribute list."""

        self.attributes[name] = value

    def checkProductFrame( self ) :
        """
        Calls checkProductFrame for self's distributions and if present for its outputChannel.
        """

        self.distribution.checkProductFrame( )
        if( self.outputChannel is not None ) : self.outputChannel.checkProductFrame( )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.multiplicity.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.multiplicity.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getLevelAsFloat( self, unit, default = 0. ) :

        if( hasattr( self.particle, 'getLevelAsFloat' ) ) : return( self.particle.getLevelAsFloat( unit, default = default ) )
        return( default )

    def getMass( self, unit ) :
        """Returns the mass of the particle if possible, otherwise None is returned."""

        return( self.particle.getMass( unit ) )

    def getAttribute( self, name ) :
        """Returns value for attribute name if it exists; otherwise, returns None."""

        if( name in self.attributes ) : return( self.attributes[name] )
        return( None )

    def calculateAverageProductData( self, style, indent, **kwargs ) :

        verbosity = kwargs['verbosity']
        indent2 = indent  + kwargs['incrementalIndent']
        indent3 = indent2 + kwargs['incrementalIndent']
        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']

        if( verbosity > 1 ) :
            print '%s%s: label = %s: calculating average product data' % ( indent, self.name, self.label )
        if( len( self.distribution ) > 0 ) :
            multiplicity = style.findFormMatchingDerivedStyles( self.multiplicity )
            kwargs['multiplicity'] = multiplicity.toPointwise_withLinearXYs( 1e-8, 1e-8 )
            kwargs['product'] = self
            energyData, momentumData = self.distribution.calculateAverageProductData( style, indent = indent2, **kwargs )
            if( energyData is not None ) :
                axes = energyDepositionModule.XYs1d.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
                if( len( energyData ) == 1 ) :
                    averageEnergy = energyDepositionModule.XYs1d( data = energyData[0], axes = axes, label = style.label, 
                            accuracy = energyAccuracy ) 
                else :
                    averageEnergy = energyDepositionModule.regions1d( axes = axes, label = style.label )
                    for energyDataRegion in energyData :
                        averageEnergyRegion = energyDepositionModule.XYs1d( data = energyDataRegion, axes = axes, 
                                accuracy = energyAccuracy )
                        averageEnergy.append( averageEnergyRegion )
                self.energyDeposition.add( averageEnergy )
            if( momentumData is not None ) :
                axes = momentumDepositionModule.XYs1d.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
                if( len( momentumData ) == 1 ) :
                    averageMomentum = momentumDepositionModule.XYs1d( data = momentumData[0], axes = axes, label = style.label, 
                            accuracy = energyAccuracy ) 
                else :
                    averageMomentum = momentumDepositionModule.regions1d( axes = axes, label = style.label ) 
                    for momentumDataRegion in momentumData :
                        averageMomentumRegion = momentumDepositionModule.XYs1d( data = momentumDataRegion, axes = axes, 
                                accuracy = momentumAccuracy ) 
                        averageMomentum.append( averageMomentumRegion )
                self.momentumDeposition.add( averageMomentum )
        if( self.outputChannel is not None ) : self.outputChannel.calculateAverageProductData( style, indent = indent3, **kwargs )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']
        status = 0

        tempInfo['workFile'].append( self.label )

        if( self.particle.name in style.transportables ) :
            if( verbosity > 1 ) :
                print '%s%s: label = %s: Sn processing' % ( indent, self.name, self.label )

            productMass = tempInfo['masses']['Product']             # Save to restore later
            tempInfo['masses']['Product'] = self.getMass( tempInfo['massUnit'] )
            tempInfo['product'] = self
            tempInfo['multiplicity'] = self.multiplicity

            self.multiplicity.processSnMultiGroup( style, tempInfo, indent )
            self.energyDeposition.processSnMultiGroup( style, tempInfo, indent )
            self.momentumDeposition.processSnMultiGroup( style, tempInfo, indent )
            try :
                self.distribution.processSnMultiGroup( style, tempInfo, indent )
            except :
                if( tempInfo['logFile'] is None ) :
                    raise
                else :
                    import traceback
                    tempInfo['logFile'].write( '\n' + self.toXLink() + ':\n' + traceback.format_exc( ) + '\n' )
                    status = 1
            tempInfo['masses']['Product'] = productMass

        if( self.outputChannel is not None ) : status += self.outputChannel.processSnMultiGroup( style, tempInfo, indent2 )
        del tempInfo['workFile'][-1]
        return( status )

    def check( self, info ):
        """ check product and distributions """
        from fudge.gnd import warning
        base = distributionsModule.base
        warnings = []

        multWarnings = self.multiplicity.check( info )
        if multWarnings:
            warnings.append( warning.context("Multiplicity:", multWarnings) )

        if( ( self.label in info['transportables'] ) and ( not self.distribution.hasData( ) ) ) :
            warnings.append( warning.missingDistribution( self.label, self ) )

        for form in self.distribution:

            if info['isTwoBody']:
                if( form.productFrame != standardsModule.frames.centerOfMassToken ) :
                    warnings.append( warning.wrong2BodyFrame( form ) )

                if form.moniker not in (angularModule.twoBodyForm.moniker,
                                        referenceModule.form.moniker,
                                        angularModule.CoulombExpansionForm.moniker,
                                        unspecifiedModule.form.moniker):
                    warnings.append( warning.wrongDistributionComponent( form.moniker, '2-body' ) )
            else:
                if form.moniker in (angularModule.twoBodyForm.moniker,
                                    angularModule.form.moniker,
                                    energyModule.form.moniker):
                    if not ( ( self.name =='gamma' ) and ( 'discrete' in self.attributes ) ) :
                        warnings.append( warning.wrongDistributionComponent( form.moniker, 'N-body' ) )

            def checkForm( form, uncorrelatedSubform='' ):
                distributionErrors = []
                try:
                    if hasattr(form, 'domain') and form.domain() != info['crossSectionDomain']:
                        for idx in range(2):    # check lower and upper edges: only warn if they disagree by > eps
                            ratio = form.domain()[idx] / info['crossSectionDomain'][idx]
                            if (ratio < 1-standardsModule.floats.epsilon or ratio > 1+standardsModule.floats.epsilon):
                                distributionErrors.append( warning.domain_mismatch(
                                        *(form.domain() + info['crossSectionDomain']), obj=form ) )
                                break
                except IndexError as err:
                    distributionErrors.append(warning.ExceptionRaised(err))
                    if info['failOnException']: raise
                if not hasattr(form,'check'):
                    warnings.append( warning.NotImplemented(form.moniker, form ) )
                    if info['failOnException']:
                        raise NotImplementedError("Checking distribution form '%s'" % form.moniker)
                else:
                    distributionErrors += form.check( info )
                if distributionErrors:
                    warnings.append( warning.context("Distribution %s %s - %s:"
                        % (form.moniker, uncorrelatedSubform, form.moniker), distributionErrors ) )

            if isinstance(form, distributionsModule.uncorrelated.form):
                for subformName in ('angularSubform','energySubform'):
                    subform = getattr(form, subformName ).data
                    checkForm( subform, subformName.replace('Subform','') )

            elif isinstance(form, KalbachMannModule.form):
                checkForm( form )

            else:
                for subform in form.subforms:
                    checkForm( subform )

        if self.outputChannel is not None:
            parentIsTwoBody = info['isTwoBody']
            info['isTwoBody'] = isinstance( self.outputChannel, channelsModule.twoBodyOutputChannel )
            for decayProduct in self.outputChannel:
                decayWarnings = decayProduct.check( info )
                if decayWarnings:
                    warnings.append( warning.context("Decay product: %s" % decayProduct.label, decayWarnings) )
            info['isTwoBody'] = parentIsTwoBody # reset to parent channel

        return warnings

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeString = ''
        for q in sorted(self.attributes) :
            if( q in [ 'discrete', 'decayRate', 'primary' ] ) :
                attributeString += ' %s="%s"' % ( q, self.attributes[q].toString( keepPeriod = False ) )
            else :
                attributeString += ' %s="%s"' % ( q, self.attributes[q] )
        xmlString = [ '%s<%s name="%s" label="%s"%s>' % ( indent, self.moniker, self.name, self.label, attributeString ) ]

        xmlString += self.multiplicity.toXMLList( indent2, **kwargs )
        xmlString += self.distribution.toXMLList( indent2, **kwargs )
        xmlString += self.energyDeposition.toXMLList( indent2, **kwargs )
        xmlString += self.momentumDeposition.toXMLList( indent2, **kwargs )

        if( not ( self.outputChannel is None ) ) :
            xmlString += self.outputChannel.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toString( self, simpleString = False, exposeGammaMultiplicity = False ) :
        """
        Returns a string representation of self. If simpleString is True, the string contains only the final 
        particles, and not any intermediate particles.
        """

        if( simpleString == True ) :
            s = self.name
            multiplicity = self.multiplicity.getXMLAttribute( )
            if( ( s != 'gamma' ) or exposeGammaMultiplicity ) :
                if( type( multiplicity ) == type( 1 ) ) :
                    if( multiplicity != 1 ) : s = "%s[multiplicity:'%s']" % ( s, multiplicity )
                else :
                    s = "%s[multiplicity:'%s']" % ( s, multiplicity )
        else :
            s = self.name
            qs = ''
            c = '['
            if( ( s != 'gamma' ) or exposeGammaMultiplicity ) :
                multiplicity = self.multiplicity.getXMLAttribute( )
                if( multiplicity != 1 ) :
                    qs += "%smultiplicity:'%s'" % ( c, multiplicity )
                    c = ', '
            for q in self.attributes :
                if( q == 'ENDFconversionFlag' ) : continue
                v = self.attributes[q]
                if( ( q == 'multiplicity' ) and ( v == 1 ) ) : continue
                qs += "%s%s:'%s'" % ( c, q, v )
                c = ', '
            if( len( qs ) > 0 ) : qs += ']'
            s = '%s%s' % ( s, qs )
            if( self.outputChannel is not None ) : s = '(%s -> %s)' % ( s, self.outputChannel )
        return( s )

    @staticmethod
    def parseXMLNode( productElement, xPath, linkData ):
        """Translate a <product> element from xml."""

        xPath.append( '%s[@label="%s"]' % ( productElement.tag, productElement.get( 'label' ) ) )
        attrs = dict( productElement.items( ) )
        particle = linkData['particles'].getParticle( attrs.pop( 'name' ) )
        prod = product( particle, label = attrs.pop( 'label' ) )
        prod.multiplicity.parseXMLNode( productElement.find( multiplicityModule.component.moniker ), xPath, linkData )
        prod.distribution.parseXMLNode( productElement.find( distributionModule.component.moniker ), xPath, linkData )
        outputChannel = productElement.find( channelsModule.outputChannelToken )
        if( outputChannel ) :
            prod.addOutputChannel( channelsModule.parseXMLNode( outputChannel, xPath, linkData ) )
        for attr in ( 'decayRate', 'primary', 'discrete' ) :
            if( attr in attrs ) :
                from pqu import PQU
                attrs[attr] = PQU.PQU( attrs[attr] )
        prod.attributes = attrs

        depositionEnergyToken = energyDepositionModule.component.moniker
        if( productElement.find( depositionEnergyToken ) is not None ) :
            prod.energyDeposition.parseXMLNode( productElement.find( depositionEnergyToken ), xPath, linkData )
        depositionMomentumToken = productData.momentumDeposition.component.moniker
        if( productElement.find( depositionMomentumToken ) is not None ) :
            prod.momentumDeposition.parseXMLNode( productElement.find( depositionMomentumToken ), xPath, linkData )
        xPath.pop( )
        return( prod )

class products( suitesModule.suite ) :

    moniker = 'products'

    def __init__( self ) :

        suitesModule.suite.__init__( self, ( product, ) )

    def uniqueLabel( self, product ) :
        """
        If product's label is the same as another product's label in self, construct a new unique label
        based on product's name appended with '__' and one or more lower case letters (i.e., 'a' to 'z').
        """

        def incrementSuffix( suffix ) :

            if( len( suffix ) == 0 ) : return( 'a' )
            index = string.ascii_lowercase.find( suffix[-1] ) + 1
            if( index != 26 ) : return( suffix[:-1] + string.ascii_lowercase[index] )
            return( incrementSuffix( suffix[:-1] ) + 'a' )

        if( product.label is None ) : product.__label = product.name
        if( product.label in self ) :
            name__ = product.name + '__'
            n1 = len( name__ )
            l1 = 0
            suffixes = []
            for _product in self :          # Find list of longest labels that start with name__.
                if( _product.label[:n1] == name__ ) :
                    suffix = _product.label[n1:]
                    if( not( suffix.islower( ) ) ) : continue       # Ignore non-standard labels.
                    l2 = len( suffix )
                    if( l2 < l1 ) : continue
                    if( l2 > l1 ) :
                        l1 = l2
                        suffixes = []
                    suffixes.append( suffix )
            if( len( suffixes ) == 0 ) :
                suffix = 'a'
            else :
                suffix = incrementSuffix( sorted( suffixes )[-1] )
            product.label = name__ + suffix
        return( product )
