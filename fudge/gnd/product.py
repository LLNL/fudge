# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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
from fudge.core.ancestry import ancestry
from fudge.core.utilities import brb
from .productData import distributions
from . import channels, productData, tokens

__metaclass__ = type

class product( ancestry ) :
    """This is the class for a gnd particle. A gnd particle can decay (i.e., breakup), the decay
    formula is defined via the decayChannel member."""

    def __init__( self, particle, ancestryName, label = None, attributes = {}, multiplicity = 1, 
            decayChannel = None ) :
        """Creates a new product object."""

        ancestry.__init__( self, ancestryName, None, attribute = 'label' )
        self.particle = particle
        self.label = label
        self.attributes = {}
        for q in attributes : self.addAttribute( q, attributes[q] )
        self.data = {}
        self.decayChannel = None
        if( decayChannel != None ) : self.addDecayChannel( decayChannel )
        if( isinstance( multiplicity, productData.multiplicity.component ) ) :
            self.multiplicity = multiplicity
        else:
            if( ( type( multiplicity ) == type( 1 ) ) or ( type( multiplicity ) == type( 1. ) ) ) :
                multiplicity = productData.multiplicity.constant( multiplicity )
            self.multiplicity = productData.multiplicity.component( multiplicity.getForm( ) )
            self.multiplicity.addForm( multiplicity )
        self.distributions = distributions.base.distribution( distributions.base.noneComponentToken )
        self.multiplicity.setParent( self )
        self.distributions.setParent( self )

    def __cmp__( self, other ) :
        """Compares self to other."""

        if( self.getName( ) < other.getName( ) ) : return( -1 )
        if( self.getName( ) > other.getName( ) ) : return(  1 )
        if( self.decayChannel < other.decayChannel ) : return( -1 )
        if( self.decayChannel > other.decayChannel ) : return(  1 )
        return( 0 )

    def __len__( self ) :
        """Returns the number of data sets for self."""

        return( len( self.data ) )

    def  __str__( self ) :
        """Converts product object to a string representation."""

        return( self.toString( simpleString = False ) )

    def addData( self, data ) :

        from productData import energyDeposition, momentumDeposition
        genre = data.genre
        if( genre not in [ energyDeposition.component.genre, momentumDeposition.component.genre ] ) : raise Exception( 'Unsupport genre = %s for data' % genre )
        if( genre in self.data ) : raise Exception( 'genre = %s already exists' % genre )
        self.data[genre] = data

    def addDecayChannel( self, decayChannel ) :
        """Adds decayChannel to particle."""

        if( isinstance( decayChannel, channels.channel ) ) :
            decayChannel.setParent( self )
            self.decayChannel = decayChannel
        else :
            raise Exception( 'Invalid decay channel = %s' % brb.getType( decayChannel ) )

    def addDistributionComponent( self, component ) : 

        self.distributions.addComponent( component )

    def addAttribute( self, name, value ) :
        """Add name and value to attribute list."""

        self.attributes[name] = value

    def checkProductFrame( self ) :
        """
        Calls checkProductFrame for self's distributions and if present for its decayChannel.
        """

        self.distributions.checkProductFrame( )
        if( self.decayChannel is not None ) : self.decayChannel.checkProductFrame( )

    def getDataTokens( self ) :
        """Returns the list of data for this particle."""

        return( self.data )

    def getDataByToken( self, dataToken ) :
        """Returns the data for dataToken."""

        if( dataToken not in self.data ) : return( None )
        return( self.data[dataToken] )

    def getDistributionNativeData( self ) :

        return( self.distributions.getNativeDataToken( ) )

    def getDistributionGenreNames( self ) : 

        return( self.distributions.keys( ) )

    def getDistributionComponentByToken( self, componentToken ) : 

        return( self.distributions[componentToken] )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.multiplicity.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.multiplicity.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getLevelAsFloat( self, unit, default = 0. ) :

        return( self.particle.getLevelAsFloat( unit, default = default ) )

    def getMass( self, unit ) :
        """Returns the mass of the particle if possible, otherwise None is returned."""

        return( self.particle.getMass( unit ) )

    def setMultiplicity( self, data ) :
        """Sets the multiplicity nativeData to be data."""

        self.multiplicity.addForm( data, asNativeData = True )

    def getName( self ) :
        """Returns self's name"""

        return( self.particle.getName( ) )

    def getToken( self ) :
        """Returns self's token, same as getName."""

        return( self.getName( ) )

    def getLabel( self ) :
        """Returns self's label."""

        return( self.label )

    def getAttribute( self, name ) :
        """Returns value for attribute name if it exists; otherwise, returns None."""

        if( name in self.attributes ) : return( self.attributes[name] )
        return( None )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        if( processInfo['verbosity'] >= 20 ) : print '%s%s: label = %s: calculating deposition data' % ( verbosityIndent, self.getName( ), self.label )
        if( self.distributions.nativeData not in [ distributions.base.noneComponentToken, distributions.base.unknownGenre ] ) :
            tempInfo['multiplicity'] = self.multiplicity
            tempInfo['product'] = self
            for newData in self.distributions.calculateDepositionData( processInfo, tempInfo ) :
                if( newData.getGenre( ) not in self.data ) :
                    self.data[newData.getGenre( )] = newData.getComponentsClass( )( newData.getForm( ) )
                self.data[newData.getGenre( )].addForm( newData )
        if( not ( self.decayChannel is None ) ) :
            for product in self.decayChannel : product.calculateDepositionData( processInfo, tempInfo, verbosityIndent + '    ' )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        doProcess = True
        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( not( processInfo.isProcessParticle( self.getName( ) ) ) ) : doProcess = False

        processInfo['workFile'].append( self.label )
        if( doProcess ) :
            if( processInfo['verbosity'] >= 20 ) : print '%s%s: label = %s: processing: (%s)' % \
                ( verbosityIndent, self.getName( ), self.label, self.distributions.nativeData )

            productMass = tempInfo['masses']['Product']             # Save to restore later
            tempInfo['masses']['Product'] = self.getMass( tempInfo['massUnit'] )

            if( self.distributions.nativeData not in [ distributions.base.noneComponentToken, distributions.base.unknownGenre ] ) :
                tempInfo['product'] = self
                tempInfo['multiplicity'] = self.multiplicity
                self.multiplicity.process( processInfo, tempInfo, verbosityIndent )
                try :
                    self.distributions.process( processInfo, tempInfo, verbosityIndent + '    ' )
                except :
                    if( processInfo['logFile'] is None ) :
                        raise
                    else :
                        import traceback
                        processInfo['logFile'].write( '\n' + self.toXLink() + ':\n' + traceback.format_exc( ) + '\n' )
                for genre in self.data : self.data[genre].process( processInfo, tempInfo, verbosityIndent + '    ' )
            tempInfo['masses']['Product'] = productMass

        if( self.decayChannel is not None ) :
            for product in self.decayChannel : product.process( processInfo, tempInfo, verbosityIndent + '    ' )
        del processInfo['workFile'][-1]

    def check( self, info ):
        """ check product and distributions """
        from fudge.gnd import warning
        base = distributions.base
        warnings = []

        multWarnings = self.multiplicity.check( info )
        if multWarnings:
            warnings.append( warning.context("Multiplicity:", multWarnings) )

        if self.getLabel() in info['transportables'] and self.distributions.nativeData == 'none':
            warnings.append( warning.missingDistribution( self.getLabel(), self ) )

        if info['isTwoBody']:
            if self.distributions.nativeData not in (base.angularComponentToken, base.referenceComponentToken,
                    base.CoulombElasticComponentToken, base.noneComponentToken):
                warnings.append( warning.wrongDistributionComponent( self.distributions.nativeData, '2-body' ) )
        else:
            if self.distributions.nativeData in (base.angularComponentToken,base.energyComponentToken):
                if not (self.getName()=='gamma' and 'discrete' in self.attributes):
                    warnings.append( warning.wrongDistributionComponent( self.distributions.nativeData, 'N-body' ) )

        for component in self.distributions.components.values():

            if info['isTwoBody']:
                nativeForm = component.forms[ component.nativeData ]
                from fudge.core.math.xData import axes
                if nativeForm.getProductFrame() != axes.centerOfMassToken :
                    warnings.append( warning.wrong2BodyFrame( nativeForm ) )

            def checkForm( form, uncorrelatedComponent='' ):
                distributionErrors = []
                if hasattr(form, 'getDomain') and form.getDomain() != info['crossSectionDomain']:
                    distributionErrors.append( warning.domain_mismatch( *(form.getDomain() + info['crossSectionDomain']),
                        obj=form ) )
                if not hasattr(form,'check'):
                    warnings.append( warning.NotImplemented(form.name, form ) )
                else:
                    distributionErrors += form.check( info )
                if distributionErrors:
                    warnings.append( warning.context("Distribution %s %s - %s:"
                        % (component.name, uncorrelatedComponent, form.name), distributionErrors ) )

            if isinstance(component, distributions.uncorrelated.component):
                # uncorrelated: same frame required for energy & angle distributions:
                frames = []
                for subcomponent in ('angularComponent','energyComponent'):
                    form = getattr(component, subcomponent).forms[ getattr(component, subcomponent).nativeData]
                    frames.append( form.getProductFrame() )
                    checkForm( form, subcomponent )
                if (len(frames)==2) and (frames[0] != frames[1]):
                    warnings.append( warning.uncorrelatedFramesMismatch( frames[0], frames[1], component ) )

            else:
                for form in component.forms.values():
                    checkForm( form )

        if self.decayChannel is not None:
            for decayProduct in self.decayChannel:
                decayWarnings = decayProduct.check( info )
                if decayWarnings:
                    warnings.append( warning.context("Decay product: %s" % decayProduct.getLabel(), decayWarnings) )

        return warnings

    def toXMLList( self, flags, indent = '' ) :

        indent2 = indent + '  '
        attributeString = ''
        for q in sorted(self.attributes) : attributeString += ' %s="%s"' % ( q, self.attributes[q] )
        xmlString = [ '%s<product name="%s" label="%s" multiplicity="%s"%s>' % \
            ( indent, self.getName( ), self.getLabel( ), self.multiplicity.getXMLAttribute( ), attributeString ) ]
        xmlString += self.distributions.toXMLList( indent2 )
        for genre in self.data : xmlString += self.data[genre].toXMLList( indent = indent2 )
        xmlString += self.multiplicity.toXMLList( indent = indent2 )
        if( not ( self.decayChannel is None ) ) :
            xmlString += self.decayChannel.toXMLList( flags, indent = indent2 )
        xmlString[-1] += '</product>'
        return( xmlString )

    def toENDF6( self, MT, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

        def getPromptOrTotalNubar( self ) :

            nativeDataMultiplicity = self.multiplicity.getNativeDataToken( )
            if( nativeDataMultiplicity == tokens.pointwiseFormToken ) :
                return( self.multiplicity.getFormByToken( tokens.pointwiseFormToken ) )
            elif( nativeDataMultiplicity == tokens.polynomialFormToken ) :
                return( self.multiplicity.getFormByToken( tokens.polynomialFormToken ) )
            else :
                raise Exception( 'Unsupported nubar form = "%s"' % nativeDataMultiplicity )

        targetInfo['product'] = self
        targetInfo['delayedNubarWeight'] = None
        if( 'emissionMode' in self.attributes ) :
            if( self.getAttribute( 'emissionMode' ) == 'delayed' ) :
                MT = 455
                if( MT not in endfMFList[5] ) : endfMFList[5][MT] = [ ]
                targetInfo['delayedNubarWeight'] = self.ENDF6_delayedNubarWeights
            elif( self.getAttribute( 'emissionMode' ) == 'prompt' ) :
                targetInfo['promptNubar'] = getPromptOrTotalNubar( self )
            elif( self.getAttribute( 'emissionMode' ) == 'total' ) :
                    targetInfo['totalNubar'] = getPromptOrTotalNubar( self )

        if( flags['verbosity'] >= 10 ) : print '%s%s: label = %s: to ENDF6:' % ( verbosityIndent, self.getName( ), self.label )
        priorMF6flag = targetInfo['doMF4AsMF6']
        if self.attributes.get('ENDFconversionFlag') == 'MF6':
            targetInfo['doMF4AsMF6'] = True # flag was set in reaction.py, but may need to be overwritten
        if( self.distributions.components ) :
            targetInfo['zapID'] = self.particle.getName( )
            if hasattr(self.particle,'groundState'):
                # ENDF wants ground state mass:
                targetInfo['particleMass'] = self.particle.groundState.getMass( 'eV/c**2' )
            else:
                targetInfo['particleMass'] = self.getMass( 'eV/c**2' )
            targetInfo['multiplicity'] = self.multiplicity
            if 'primary' in self.attributes:
                targetInfo['primaryGammaEnergy'] = self.attributes['primary'].getValueAs('eV')
            self.distributions.toENDF6( MT, endfMFList, flags, targetInfo )
        if( not( self.decayChannel is None ) ) :
            priorIndex, priorToken, priorLabel = targetInfo['productIndex'], targetInfo['productToken'], targetInfo['productLabel']
            for index, product in enumerate( self.decayChannel ) :
                if( product.getName( ) == 'gamma' ) :
                    targetInfo['gammas'].append( product )
                    continue
                targetInfo['productIndex'] = "%s.%s" % ( priorIndex, index )
                targetInfo['productToken'] = product.getToken( )
                targetInfo['productLabel'] = product.getLabel( )
                product.toENDF6( MT, endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent + '    ' )
            targetInfo['productIndex'], targetInfo['productToken'], targetInfo['productLabel'] = priorIndex, priorToken, priorLabel
        targetInfo['doMF4AsMF6'] = priorMF6flag

    def toString( self, simpleString = False, exposeGammaMultiplicity = False ) :
        """Returns a string representation of self. If simpleString is True, the string contains only the final 
    particles, and not any intermediate particles."""

        if( simpleString == True ) :
            s = self.getName( )
            multiplicity = self.multiplicity.getXMLAttribute( )
            if( ( s != 'gamma' ) or exposeGammaMultiplicity ) :
                if( type( multiplicity ) == type( 1 ) ) :
                    if( multiplicity != 1 ) : s = "%s[multiplicity:'%s']" % ( s, multiplicity )
                else :
                    s = "%s[multiplicity:'%s']" % ( s, multiplicity )
        else :
            s = self.getName( )
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
            if( self.decayChannel != None ) : s = '(%s -> %s)' % ( s, self.decayChannel )
        return( s )

def parseXMLNode( productElement, xPath=[], linkData={} ):
    """Translate a <product> element from xml."""

    xPath.append( '%s[@label="%s"]' % (productElement.tag, productElement.get('label')) )
    attrs = dict( productElement.items() )
    particle = linkData['particles'].getParticle( attrs.pop('name') )
    mult = attrs.pop('multiplicity')
    if mult == tokens.unknownFormToken:
        mult = productData.multiplicity.unknown()
    elif mult == productData.multiplicity.energyDependent:
        mult = productData.multiplicity.parseXMLNode( productElement.find('multiplicity'), xPath, linkData )
    elif mult == tokens.partialProductionFormToken:
        mult = productData.multiplicity.partialProduction()
    else:
        mult = float( mult )
        if mult.is_integer(): mult = int(mult)
    prod = product( particle, productElement.tag, label=attrs.pop('label'), multiplicity=mult)
    prod.distributions = productData.distributions.parseXMLNode( productElement.find('distributions'),
            xPath, linkData )
    prod.distributions.parent = prod
    decayChannel = productElement.find( channels.decayChannelToken )
    if decayChannel:
        prod.decayChannel = channels.parseXMLNode( decayChannel, xPath, linkData )
        prod.decayChannel.parent = prod
    for attr in ('decayRate','primary','discrete'):
        if attr in attrs:
            from pqu import physicalQuantityWithUncertainty
            attrs[attr] = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(attrs[attr])
    prod.attributes = attrs

    # extra data may be present in derived files:
    depositionEnergyToken = productData.base.energyDepositionToken
    if productElement.find(depositionEnergyToken) is not None:
        prod.data[depositionEnergyToken] = productData.energyDeposition.parseXMLNode(
                productElement.find(depositionEnergyToken), xPath )
    depositionMomentumToken = productData.base.momentumDepositionToken
    if productElement.find(depositionMomentumToken) is not None:
        prod.data[depositionMomentumToken] = productData.momentumDeposition.parseXMLNode(
                productElement.find(depositionMomentumToken), xPath )
    xPath.pop()
    return prod
