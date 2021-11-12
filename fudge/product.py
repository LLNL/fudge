# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the product class.

GNDS classes store reaction data in a hierarchical format. At the top is the product class. Next is the
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
    Contains a cross section (defined in the gnds.reactionData.crossSection module), along with
    an outgoing channel which must be a outputChannel instance (e.g., 2n + Pu239).
    
:reactionSuite:
    Contains all reactions sharing the same input channel (e.g. n + Pu239)
    
        (e.g., n + Pu239 --> n + Pu239, n + Pu239 --> n + Pu239_e1,
        n + Pu239 --> n + Pu239_e2, n + Pu239 --> 2n + Pu239 ... )
"""

from PoPs import IDs as IDsPoPsModule
from PoPs import misc as miscPoPsModule
from PoPs.decays import misc as miscDecaysPoPsModule
from PoPs.families import particle as PoPsParticleModule

from fudge.core.utilities import brb

from xData import ancestry as ancestryModule
from xData import standards as standardsModule

from fudge import reactionProducts as reactionProductsModule
from fudge import outputChannel as outputChannelModule
from fudge import suites as suitesModule

from .productData.distributions import distribution as distributionModule
from .productData.distributions import unspecified as unspecifiedModule
from .productData.distributions import miscellaneous as miscellaneousModule
from .productData import multiplicity as multiplicityModule
from .productData import energyDeposition as energyDepositionModule
from .productData import momentumDeposition as momentumDepositionModule

__metaclass__ = type

class product( ancestryModule.ancestry ) :
    """
    This is the class for a gnds product. If the product can decay (e.g. for breakup reactions),
    the resulting decay information is defined in the product outputChannel
    """

    moniker = 'product'
    ancestryMembers = ( 'multiplicity', 'distribution', 'outputChannel', 'energyDeposition', 'momentumDeposition' )

    def __init__( self, id, label = None, outputChannel = None ) :
        """Creates a new product object."""

        ancestryModule.ancestry.__init__( self )

        self.__id = id
        self.__label = label
        self.__particle = None  # lazy evaluation

        self.__outputChannel = None
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
        if( self.__outputChannel < other.outputChannel ) : return( -1 )
        if( self.__outputChannel > other.outputChannel ) : return(  1 )
        return( 0 )

    def  __str__( self ) :
        """Converts product object to a string representation."""

        return( self.toString( simpleString = False ) )

    @property
    def id( self ) :

        return( self.__id )

    @property
    def pid( self ) :

        return( self.id )

    @property
    def outputChannel( self ) :

        return( self.__outputChannel )
    
    @property
    def particle( self ) :

        if( self.__particle is None ) :
            pops = self.findAttributeInAncestry( 'PoPs' )
            try :
                self.__particle = pops[self.id]
            except :
                baseName, anti, qualifier = miscPoPsModule.baseAntiQualifierFromID( self.id, qualifierAllowed = True )
                name = baseName + anti
                self.__particle = pops.chemicalElements.getSymbol( name )
        return( self.__particle )

    @property
    def label( self ) :
        """Returns self's label."""

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = value

    @property
    def parentProduct( self ) :
        """Returns the parent product of self or None, if self does not have a parent product."""

        try :
            return( self.ancestor.findClassInAncestry( self.__class__ ) )
        except :
            return( None )

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

        if( isinstance( outputChannel, outputChannelModule.outputChannel ) ) :
            self.__outputChannel = outputChannel
            self.__outputChannel.setAncestor( self )
        else :
            raise TypeError( 'Invalid decay channel = %s' % brb.getType( outputChannel ) )

    def amendForPatch( self, fromLabel, toLabel ) :

        self.__multiplicity.amendForPatch( fromLabel, toLabel )
        self.__energyDeposition.amendForPatch( fromLabel, toLabel )
        self.__momentumDeposition.amendForPatch( fromLabel, toLabel )
        self.__distribution.amendForPatch( fromLabel, toLabel )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.multiplicity.convertUnits( unitMap )
        self.energyDeposition.convertUnits( unitMap )
        self.momentumDeposition.convertUnits( unitMap )
        self.distribution.convertUnits( unitMap )
        if( self.__outputChannel is not None ) : self.__outputChannel.convertUnits( unitMap )

    def checkProductFrame( self ) :
        """
        Calls checkProductFrame for self's distributions and if present for its outputChannel.
        """

        self.distribution.checkProductFrame( )
        if( self.__outputChannel is not None ) : self.__outputChannel.checkProductFrame( )

    def cullStyles( self, styleList ) :

        self.__multiplicity.cullStyles( styleList )
        self.__distribution.cullStyles( styleList )
        self.__energyDeposition.cullStyles( styleList )
        self.__momentumDeposition.cullStyles( styleList )

    @property
    def domainMin( self ) :

        return( self.multiplicity.domainMin )

    @property
    def domainMax( self ) :

        return( self.multiplicity.domainMax )

    @property
    def domainUnit( self ) :

        return( self.multiplicity.domainUnit )

    def diff( self, other, diffResults ) :

        self.distribution.diff( other.distribution, diffResults )

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers :
            if( getattr( self, ancestryMember ) is None ) : continue
            getattr( self, ancestryMember ).findLinks( links )

    def thresholdQAs( self, unit, final = True ) :

        if( self.__outputChannel is not None ) : return( self.__outputChannel.thresholdQAs( unit, final = final ) )
        return( 0. )

    def getLevelAsFloat( self, unit, default = 0. ) :

        if( hasattr( self.particle, 'getLevelAsFloat' ) ) : return( self.particle.getLevelAsFloat( unit, default = default ) )
        return( default )

    def getMass( self, unit ) :
        """Returns the mass of the particle if possible, otherwise None is returned."""

        particle = self.particle
        PoPs = self.findAttributeInAncestry('PoPs')
        if particle.id in PoPs.aliases:
            particle = PoPs[particle.pid]
        return( particle.getMass( unit ) )

    def calculateAverageProductData( self, style, indent, **kwargs ) :

        verbosity = kwargs['verbosity']
        indent2 = indent  + kwargs['incrementalIndent']
        indent3 = indent2 + kwargs['incrementalIndent']

        reactionSuite = kwargs['reactionSuite']
        energyUnit = kwargs['incidentEnergyUnit']
        momentumUnit = kwargs['momentumUnit']
        massUnit = kwargs['massUnit']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']

        kwargs['product'] = self
        try :
            kwargs['productMass'] = reactionSuite.PoPs[self.id].getMass( massUnit )
        except :                        # Can happend when mass is not needed for evaluation and hence not stored.
            kwargs['productMass'] = None

        if( verbosity > 1 ) :
            print('%s%s: label = %s: calculating average product data' % (indent, self.id, self.label))
        if( len( self.distribution ) > 0 ) :
            multiplicity = style.findFormMatchingDerivedStyle( self.multiplicity )
            if( not( isinstance( multiplicity, (multiplicityModule.branching1d, multiplicityModule.unspecified) ) ) ) :
                kwargs['multiplicity'] = multiplicity.toPointwise_withLinearXYs( accuracy = 1e-5, lowerEps = 1e-6, upperEps = 1e-6 )

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
                axes = momentumDepositionModule.defaultAxes( energyUnit = energyUnit, momentumUnit = momentumUnit )
                if( len( momentumData ) == 1 ) :
                    averageMomentum = momentumDepositionModule.XYs1d( data = momentumData[0], axes = axes, label = style.label ) 
                else :
                    averageMomentum = momentumDepositionModule.regions1d( axes = axes, label = style.label ) 
                    for momentumDataRegion in momentumData :
                        averageMomentumRegion = momentumDepositionModule.XYs1d( data = momentumDataRegion, axes = axes ) 
                        averageMomentum.append( averageMomentumRegion )
                self.momentumDeposition.add( averageMomentum )

        if( self.__outputChannel is not None ) :
            self.__outputChannel.calculateAverageProductData( style, indent = indent3, **kwargs )

    def partialProductionIntegral( self, reaction_suite, productID, energyIn, energyOut = None, muOut = None, phiOut = None, 
            frame = standardsModule.frames.productToken, LegendreOrder = 0, **kwargs ) :

        def branchingGammas( initialState, photonBranchingData, probability, LegendreOrder ) :

            partialProductionIntegralSum = 0.0
            if( LegendreOrder != 0 ) : return( partialProductionIntegralSum )           # Currently, assume isotropic scattering in lab frame.
            if( initialState in photonBranchingData ) :
                gammas = photonBranchingData[initialState]['photons']
                for branchingRatio, gammaEnergy, finalState in gammas :
                    gammaEnergy = float( gammaEnergy )
                    energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, gammaEnergy, gammaEnergy )
                    if( energyOutMin <= gammaEnergy <= energyOutMax ) : partialProductionIntegralSum += branchingRatio * probability
                    partialProductionIntegralSum += branchingGammas( finalState, photonBranchingData, branchingRatio * probability, LegendreOrder )

            return( partialProductionIntegralSum )

        partialProductionIntegralSum = 0.0

        if( self.id == productID ) :
            if( self.distribution.hasData( ) ) :
                multiplicity = self.multiplicity[0]
                if( isinstance( multiplicity, multiplicityModule.branching1d ) ) :
                    photonBranchingData = miscDecaysPoPsModule.photonBranchingData( reaction_suite.PoPs, self.parentProduct.id )
                    factor = miscellaneousModule.muPhiEvaluate( muOut, phiOut )
                    partialProductionIntegralSum += factor * branchingGammas( multiplicity.product.parentProduct.id, photonBranchingData, 1.0, LegendreOrder )
                else :
                    domainMin = self.multiplicity[0].domainMin
                    domainMax = self.multiplicity[0].domainMax
                    if( domainMin <= energyIn <= domainMax ) :
                        multiplicity = self.multiplicity.evaluate( energyIn )
                        if( multiplicity != 0.0 ) :
                            partialProductionIntegralSum = self.distribution.integrate( reaction_suite, energyIn, energyOut = energyOut, 
                                    muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder )
                            partialProductionIntegralSum *= multiplicity

        if( self.__outputChannel is not None ) :
            partialProductionIntegralSum += self.__outputChannel.partialProductionIntegral( reaction_suite, productID, energyIn, energyOut = energyOut, 
                    muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder )

        return( partialProductionIntegralSum )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']
        if( verbosity > 1 ) : print('%s%s: label = %s: MonteCarlo_cdf processing' % ( indent, self.id, self.label ))

        self.distribution.processMC_cdf( style, tempInfo, indent )
        if( self.__outputChannel is not None ) : self.__outputChannel.processMC_cdf( style, tempInfo, indent2 )

    def processMultiGroup( self, style, tempInfo, indent ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        tempInfo['workFile'].append( self.label )

        doIt = not( isinstance( self.distribution[0], unspecifiedModule.form ) )
        if( doIt and ( self.id in style.transportables ) ) :
            if( verbosity > 1 ) :
                print('%s%s: label = %s: multiGroup processing' % ( indent, self.id, self.label ))

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

        if( self.__outputChannel is not None ) : self.__outputChannel.processMultiGroup( style, tempInfo, indent2 )

        del tempInfo['workFile'][-1]

    def check( self, info ):
        """ check product multiplicity, distribution and breakup products (if applicable) """
        from fudge import warning
        warnings = []

        multWarnings = self.multiplicity.check( info )
        if multWarnings:
            warnings.append( warning.context("Multiplicity:", multWarnings) )

        if( ( self.label in info['transportables'] ) and ( not self.distribution.hasData( ) ) ) :
            warnings.append( warning.missingDistribution( self.label, self ) )

        distributionWarnings = self.distribution.check( info )
        if distributionWarnings:
            warnings.append( warning.context("Distribution:", distributionWarnings) )

        if self.__outputChannel is not None:
            parentIsTwoBody = info['isTwoBody']
            info['isTwoBody'] = self.__outputChannel.genre == outputChannelModule.Genre.twoBody
            for decayProduct in self.__outputChannel:
                decayWarnings = decayProduct.check( info )
                if decayWarnings:
                    warnings.append( warning.context("Decay product: %s" % decayProduct.label, decayWarnings) )
            info['isTwoBody'] = parentIsTwoBody # reset to parent channel

        return warnings

    def reactionProducts( self, _reactionProducts ) :

        category = reactionProductsModule.Category.particle
        if( self.__outputChannel is not None ) : category = reactionProductsModule.Category.intermediate

        multiplicity1 = self.multiplicity.evaluated
        value = 0
        if( isinstance( multiplicity1, multiplicityModule.constant1d ) ) : value = multiplicity1.value
        if( int( value ) == value ) : value = int( value )

        pid = self.pid
        pops = self.findAttributeInAncestry( 'PoPs' )
        particle = pops[pid]
        if( isinstance( particle, PoPsParticleModule.alias ) ) : pid = pops[particle.pid].id

        if( pid not in _reactionProducts ) : _reactionProducts[pid] = reactionProductsModule.ReactionProduct( category, 0 )
        _reactionProducts[pid].value += value

        if( self.__outputChannel is not None ) : self.__outputChannel.reactionProducts( _reactionProducts )

    def removeStyles( self, styleLabels ) :

        self.__multiplicity.removeStyles( styleLabels )
        self.__distribution.removeStyles( styleLabels )
        self.__energyDeposition.removeStyles( styleLabels )
        self.__momentumDeposition.removeStyles( styleLabels )
        if( self.__outputChannel is not None ) : self.__outputChannel.removeStyles( styleLabels )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s pid="%s" label="%s">' % ( indent, self.moniker, self.id, self.label ) ]

        xmlString += self.multiplicity.toXMLList( indent2, **kwargs )
        xmlString += self.distribution.toXMLList( indent2, **kwargs )
        xmlString += self.energyDeposition.toXMLList( indent2, **kwargs )
        xmlString += self.momentumDeposition.toXMLList( indent2, **kwargs )

        if( self.__outputChannel is not None ) :
            xmlString += self.__outputChannel.toXMLList( indent2, **kwargs )

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
                    iValue = int( _multiplicity.value )
                    if( iValue == _multiplicity.value ) :
                        if( iValue != 1 ) : multiplicity = '%s' % iValue
                    else :
                        multiplicity = '%s ' % _multiplicity.value
                else :
                    multiplicity = 'm(E)*'
            return( multiplicity )

        multiplicity = multiplicityString( )
        if( simpleString == True ) :
            productName = '%s%s' % ( multiplicity, self.id )
        else :
            productName = self.id
            productName = '%s%s' % ( multiplicity, productName )
            if( self.__outputChannel is not None ) : productName = '(%s -> %s)' % ( productName, self.__outputChannel )
        return( productName )

    @staticmethod
    def parseXMLNode( productElement, xPath, linkData ):
        """Translate a <product> element from xml."""

        def parseChildNode( moniker, parser ) :

            child = productElement.find( moniker )
            if( child is not None ) : parser( child, xPath, linkData )

        xPath.append( '%s[@label="%s"]' % ( productElement.tag, productElement.get( 'label' ) ) )

        attrs = dict( productElement.items( ) )
        prod = product( id = attrs.pop( 'pid' ), label = attrs.pop( 'label' ) )
        parseChildNode( multiplicityModule.component.moniker, prod.multiplicity.parseXMLNode )
        parseChildNode( distributionModule.component.moniker, prod.distribution.parseXMLNode )

        outputChannel = productElement.find( outputChannelModule.outputChannel.moniker )
        if( outputChannel is not None ) :
            prod.addOutputChannel( outputChannelModule.parseXMLNode( outputChannel, xPath, linkData ) )

        parseChildNode( energyDepositionModule.component.moniker, prod.energyDeposition.parseXMLNode )
        parseChildNode( momentumDepositionModule.component.moniker, prod.momentumDeposition.parseXMLNode )

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
