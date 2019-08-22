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

GND classes store reaction data in a hierarchical format. At the top is the product class. Next is the
channel, representing the reaction output channels. An output channel contains the Q-value and a list
of outgoing products (some of which may decay into other particles).
The channel is stored inside the reaction class, which consists of an input and output channel together,
(e.g. n_1 + O_16  --> H_1 + N_16). Finally, the reactionSuite class is used to store a list of reaction
for a given input channel (that is, for the incoming particles that come together in the reaction).
The following outline summarizes the classes.

:product:
    Outgoing particle, including multiplicity and distribution information.

:channel:
    A Q-value plus a list of products.

:reaction:
    Contains a cross section (defined in the gnd.reactionData.crossSection module), along with
    an outgoing channel which must be one of the classes defined in gnd.channels (e.g., 2n + Pu239).
    
:reactionSuite:
    Contains all reactions sharing the same input channel (e.g. n + Pu239)
    
        (e.g., n + Pu239 --> n + Pu239, n + Pu239 --> n + Pu239_e1,
        n + Pu239 --> n + Pu239_e2, n + Pu239 --> 2n + Pu239 ... )
"""

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from fudge.core.utilities import brb

from xData import ancestry as ancestryModule
from xData import physicalQuantity as physicalQuantityModule

from fudge.gnd import suites as suitesModule

from .productData.distributions import distribution as distributionModule
from .productData.distributions import unspecified as unspecifiedModule
from .productData import multiplicity as multiplicityModule
from .productData import energyDeposition as energyDepositionModule
from .productData import momentumDeposition as momentumDepositionModule
from . import channels as channelsModule

__metaclass__ = type

class product( ancestryModule.ancestry ) :
    """
    This is the class for a gnd product. If the product can decay (e.g. for breakup reactions),
    the resulting decay information is defined in the product outputChannel
    """

    moniker = 'product'
    ancestryMembers = ( 'multiplicity', 'distribution', 'outputChannel',
                        'energyDeposition', 'momentumDeposition' )

    def __init__( self, id, label = None, attributes = None, outputChannel = None ) :
        """Creates a new product object."""

        ancestryModule.ancestry.__init__( self )

        self.__id = id
        self.__label = label
        self.__particle = None  # lazy evaluation
        self.attributes = {}
        if attributes is not None:
            for q in attributes : self.addAttribute( q, attributes[q] )

        self.outputChannel = None
        if( outputChannel is not None ) : self.addOutputChannel( outputChannel )

        self.__multiplicity = multiplicityModule.component( )
        self.__multiplicity.setAncestor( self )

        self.__energyDeposition = energyDepositionModule.component( )
        self.__energyDeposition.setAncestor( self )

        self.__momentumDeposition = momentumDepositionModule.component( )
        self.__momentumDeposition.setAncestor( self )

        self.__distribution = distributionModule.component( )
        self.__distribution.setAncestor( self )

    def __cmp__( self, other ) :
        """Compares self to other."""

        if( self.id < other.id ) : return( -1 )
        if( self.id > other.id ) : return(  1 )
        if( self.outputChannel < other.outputChannel ) : return( -1 )
        if( self.outputChannel > other.outputChannel ) : return(  1 )
        return( 0 )

    def  __str__( self ) :
        """Converts product object to a string representation."""

        return( self.toString( simpleString = False ) )

    @property
    def id( self ) :

        return( self.__id )

    @property
    def name( self ) :

        return( self.id )

    @property
    def particle(self):
        if self.__particle is None:
            self.__particle = self.findAttributeInAncestry('PoPs')[self.id]
        return self.__particle

    @property
    def label( self ) :
        """Returns self's label."""

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = value

    @property
    def multiplicity( self ) :

        return( self.__multiplicity )

    @property
    def energyDeposition( self ) :

        return( self.__energyDeposition )

    @property
    def momentumDeposition( self ) :

        return( self.__momentumDeposition )

    @property
    def distribution( self ) :

        return( self.__distribution )

    def addOutputChannel( self, outputChannel ) :
        """Adds outputChannel to particle."""

        if( isinstance( outputChannel, channelsModule.channel ) ) :
            outputChannel.setAncestor( self )
            self.outputChannel = outputChannel
        else :
            raise TypeError( 'Invalid decay channel = %s' % brb.getType( outputChannel ) )

    def addAttribute( self, name, value ) :
        """Add name and value to attribute list."""

        self.attributes[name] = value

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        if( 'originationLevel' in self.attributes ) :
            originationLevel = PQUModule.PQU( self.attributes['originationLevel'] )
            originationLevel = physicalQuantityModule.physicalQuantity( float( originationLevel ), originationLevel.unit )
            originationLevel.convertUnits( unitMap )
            self.attributes['originationLevel'] = originationLevel.toString( significantDigits = 12 )
        self.multiplicity.convertUnits( unitMap )
        self.energyDeposition.convertUnits( unitMap )
        self.momentumDeposition.convertUnits( unitMap )
        self.distribution.convertUnits( unitMap )
        if( self.outputChannel is not None ) : self.outputChannel.convertUnits( unitMap )

    def checkProductFrame( self ) :
        """
        Calls checkProductFrame for self's distributions and if present for its outputChannel.
        """

        self.distribution.checkProductFrame( )
        if( self.outputChannel is not None ) : self.outputChannel.checkProductFrame( )

    @property
    def domainMin( self ) :

        return( self.multiplicity.domainMin )

    @property
    def domainMax( self ) :

        return( self.multiplicity.domainMax )

    @property
    def domainUnit( self ) :

        return( self.multiplicity.domainUnit )

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

        reactionSuite = kwargs['reactionSuite']
        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = kwargs['momentumDepositionUnit']
        massUnit = kwargs['massUnit']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']

        kwargs['product'] = self
        try :
            kwargs['productMass'] = reactionSuite.PoPs[self.id].getMass( massUnit )
        except :                        # Can happend when mass is not needed for evaluation and hence not stored.
            kwargs['productMass'] = None

        if( verbosity > 1 ) :
            print '%s%s: label = %s: calculating average product data' % ( indent, self.id, self.label )
        if( len( self.distribution ) > 0 ) :
            multiplicity = style.findFormMatchingDerivedStyle( self.multiplicity )
            kwargs['multiplicity'] = multiplicity.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
            energyData, momentumData = self.distribution.calculateAverageProductData( style, indent = indent2, **kwargs )
            if( energyData is not None ) :
                axes = energyDepositionModule.defaultAxes( energyUnit = energyUnit )
                if( len( energyData ) == 1 ) :
                    averageEnergy = energyDepositionModule.XYs1d( data = energyData[0], axes = axes, label = style.label ) 
                else :
                    averageEnergy = energyDepositionModule.regions1d( axes = axes, label = style.label )
                    for energyDataRegion in energyData :
                        averageEnergyRegion = energyDepositionModule.XYs1d( data = energyDataRegion, axes = axes )
                        averageEnergy.append( averageEnergyRegion )
                self.energyDeposition.add( averageEnergy )
            if( momentumData is not None ) :
                axes = momentumDepositionModule.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
                if( len( momentumData ) == 1 ) :
                    averageMomentum = momentumDepositionModule.XYs1d( data = momentumData[0], axes = axes, label = style.label ) 
                else :
                    averageMomentum = momentumDepositionModule.regions1d( axes = axes, label = style.label ) 
                    for momentumDataRegion in momentumData :
                        averageMomentumRegion = momentumDepositionModule.XYs1d( data = momentumDataRegion, axes = axes ) 
                        averageMomentum.append( averageMomentumRegion )
                self.momentumDeposition.add( averageMomentum )
        if( self.outputChannel is not None ) : self.outputChannel.calculateAverageProductData( style, indent = indent3, **kwargs )

    def processMC( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']
        if( verbosity > 1 ) : print '%s%s: label = %s: MC processing' % ( indent, self.id, self.label )

        self.distribution.processMC( style, tempInfo, indent )
        if( self.outputChannel is not None ) : self.outputChannel.processMC( style, tempInfo, indent2 )

    def processMultiGroup( self, style, tempInfo, indent ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        tempInfo['workFile'].append( self.label )

        doIt = not( isinstance( self.distribution[0], unspecifiedModule.form ) )
        if( doIt and ( self.id in style.transportables ) ) :
            if( verbosity > 1 ) :
                print '%s%s: label = %s: multiGroup processing' % ( indent, self.id, self.label )

            productMass = tempInfo['masses']['Product']             # Save to restore later
            tempInfo['masses']['Product'] = self.getMass( tempInfo['massUnit'] )
            tempInfo['product'] = self
            tempInfo['multiplicity'] = self.multiplicity

            self.multiplicity.processMultiGroup( style, tempInfo, indent )
            self.energyDeposition.processMultiGroup( style, tempInfo, indent )
            self.momentumDeposition.processMultiGroup( style, tempInfo, indent )
            try :
                self.distribution.processMultiGroup( style, tempInfo, indent )
            except :
                if( tempInfo['logFile'] is None ) :
                    raise
                else :
                    import traceback
                    tempInfo['logFile'].write( '\n' + self.toXLink() + ':\n' + traceback.format_exc( ) + '\n' )
                    tempInfo['failures'] += 1
            tempInfo['masses']['Product'] = productMass

        if( self.outputChannel is not None ) : self.outputChannel.processMultiGroup( style, tempInfo, indent2 )
        del tempInfo['workFile'][-1]

    def check( self, info ):
        """ check product multiplicity, distribution and breakup products (if applicable) """
        from fudge.gnd import warning
        warnings = []

        multWarnings = self.multiplicity.check( info )
        if multWarnings:
            warnings.append( warning.context("Multiplicity:", multWarnings) )

        if( ( self.label in info['transportables'] ) and ( not self.distribution.hasData( ) ) ) :
            warnings.append( warning.missingDistribution( self.label, self ) )

        distributionWarnings = self.distribution.check( info )
        if distributionWarnings:
            warnings.append( warning.context("Distribution:", distributionWarnings) )

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
        xmlString = [ '%s<%s name="%s" label="%s"%s>' % ( indent, self.moniker, self.id, self.label, attributeString ) ]

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

        def multiplicityString( ) :

            multiplicity = ''
            if( ( self.id != IDsPoPsModule.photon ) or exposeGammaMultiplicity ) :
                _multiplicity = self.multiplicity.evaluated
                if( isinstance( _multiplicity, multiplicityModule.constant1d ) ) :
                    iValue = int( _multiplicity.constant )
                    if( iValue == _multiplicity.constant ) :
                        if( iValue != 1 ) : multiplicity = '%s' % iValue
                    else :
                        multiplicity = '%s ' % multiplicity.constant
                else :
                    multiplicity = 'm(E)*'
            return( multiplicity )

        multiplicity = multiplicityString( )
        if( simpleString == True ) :
            productName = '%s%s' % ( multiplicity, self.id )
        else :
            productName = self.id
            qs = ''
            c = '['
            for q in self.attributes :
                if( q == 'ENDFconversionFlag' ) : continue
                v = self.attributes[q]
                qs += "%s%s:'%s'" % ( c, q, v )
                c = ', '
            if( len( qs ) > 0 ) : qs += ']'
            productName = '%s%s%s' % ( multiplicity, productName, qs )
            if( self.outputChannel is not None ) : productName = '(%s -> %s)' % ( productName, self.outputChannel )
        return( productName )

    @staticmethod
    def parseXMLNode( productElement, xPath, linkData ):
        """Translate a <product> element from xml."""

        xPath.append( '%s[@label="%s"]' % ( productElement.tag, productElement.get( 'label' ) ) )
        attrs = dict( productElement.items( ) )
        prod = product( id = attrs.pop( 'name' ), label = attrs.pop( 'label' ) )
        prod.multiplicity.parseXMLNode( productElement.find( multiplicityModule.component.moniker ), xPath, linkData )
        prod.distribution.parseXMLNode( productElement.find( distributionModule.component.moniker ), xPath, linkData )
        outputChannel = productElement.find( channelsModule.outputChannelToken )
        if( outputChannel ) :
            prod.addOutputChannel( channelsModule.parseXMLNode( outputChannel, xPath, linkData ) )
        for attr in ( 'decayRate', 'primary', 'discrete' ) :
            if( attr in attrs ) :
                attrs[attr] = PQUModule.PQU( attrs[attr] )
        prod.attributes = attrs

        depositionEnergyToken = energyDepositionModule.component.moniker
        if( productElement.find( depositionEnergyToken ) is not None ) :
            prod.energyDeposition.parseXMLNode( productElement.find( depositionEnergyToken ), xPath, linkData )
        depositionMomentumToken = momentumDepositionModule.component.moniker
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

        if( product.label is None ) : product.label = product.id
        return( suitesModule.suite.uniqueLabel( self, product ) )
