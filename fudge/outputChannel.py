# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
In GNDS, an output channel (dubbed outputChannel) contains all outgoing products from either a reaction (also called an interaction)
or since a product may subsequently decay, the decay products for a for. This means that a reaction node has an outputChannel and a product may
have an outputChannel.

For example, consider the following nuclear reaction:

        n + O16 --> n + ( O16_e6 -> He4 + C12 )

This interaction represents a neutron hitting an Oxygen-16 target producing a n and Oxygen-16 in the 6th excited. The
excited Oxygen-16 subsequently decays into He4 and C12. The output channel for this interaction is "n + O16_e6" and
the "O16_e6" contains the decay channel "He4 + C12".

Output and decay channels can be further divided into twoBody or NBody channels.

    :twoBody:       The channel initially has two products whose angular and energy correlates are known. Only the 
                    center-of-mass angular data for one of the products is needed. The products' outgoing energies can 
                    be calculated from their center-of-mass angular information and other data like the masses of the 
                    products involved in a reaction or decay.

    :NBody:         The channel has two or more products whose angular and energy correlates are not known. The full double
                    differential data for each product must be given.
"""

__metaclass__ = type

from PoPs import IDs as IDsPoPsModule

from xData import ancestry as ancestryModule
from xData import standards as standardsModule

from fudge import reactionProducts as reactionProductsModule
from .channelData import Q as QModule

class Genre :

    twoBody = 'twoBody'
    NBody = 'NBody'
    production = 'production'
    sumOfRemainingOutputChannels = 'sumOfRemainingOutputChannels'

    allowed = ( twoBody, NBody, production, sumOfRemainingOutputChannels )

class processes :

    continuum = 'continuum'

class outputChannel( ancestryModule.ancestry ) :
    """This is the GNDS outputChannel which is a Q value and a list of products (i.e., product objects)."""

    moniker = 'outputChannel'
    ancestryMembers = ( 'Q', 'products', 'fissionFragmentData' )

    def __init__( self, genre, process = None ) :
        """Constructor for outputChannel."""

        from . import product as productModule
        from .channelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule

        ancestryModule.ancestry.__init__( self )

        self.__genre = genre                            # Need to check for valid genre.
        self.__process = process

        self.__Q = QModule.component( )
        self.__Q.setAncestor( self )

        self.__products = productModule.products( )
        self.__products.setAncestor( self )

        self.__fissionFragmentData = fissionFragmentDataModule.fissionFragmentData( )
        self.__fissionFragmentData.setAncestor( self )

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers :
            if( ancestryMember in ( 'fissionFragmentData' ) ) : continue
            getattr( self, ancestryMember ).findLinks( links )

    def __getitem__( self, i ) :
        """Returns the i^th product in the channel."""

        return( self.__products[i] )

    def __len__( self ) :
        """Returns the number of products for this channel."""

        return( len( self.__products ) )

    def __str__( self ) :
        """Converts channel object to a string representation."""

        return( self.toString( simpleString = False ) )

    def __cmp__( self, other ) :
        """Compares self to other."""

        if( other is None ) : return( -1 )
        n1 = len( self )
        n2 = len( other )
        n = min( n1, n2 )
        for i1 in range( n ) :
            if( self[i1] < other[i1] ) : return( -1 )
            if( self[i1] > other[i1] ) : return(  1 )
        if( n1 < n2 ) : return( -1 )
        if( n1 > n2 ) : return(  1 )
        return( 0 )

    @property
    def Q( self ) :

        return( self.__Q )

    @property
    def products( self ) :

        return( self.__products )

    @property
    def fissionFragmentData( self ) :

        return( self.__fissionFragmentData )

    @property
    def genre( self ) :

        return( self.__genre )

    @genre.setter
    def genre( self, value ) :

        self.__genre = value

    @property
    def process( self ) :

        return( self.__process )

    @process.setter
    def process( self, value ) :

        self.__process = value

    def amendForPatch( self, fromLabel, toLabel ) :

        self.__Q.amendForPatch( fromLabel, toLabel )
        self.__products.amendForPatch( fromLabel, toLabel )
        self.__fissionFragmentData.amendForPatch( fromLabel, toLabel )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__Q.convertUnits( unitMap )
        self.__products.convertUnits( unitMap )
        self.__fissionFragmentData.convertUnits( unitMap )

    def checkProductFrame( self ) :
        """Calls checkProductFrame for self's products."""

        for product in self : product.checkProductFrame( )

    def diff( self, other, diffResults ) :

        self.Q.diff( other.Q, diffResults )

        aliases = self.rootAncestor.PoPs.aliases
        productIDs = []
        productIDs1 = {}
        for product in self.products :
            pid = product.id
            if( pid in aliases ) : pid = aliases[pid].pid             # Handle old GNDS files which used meta-stable ids instead of nuclear level ids.
            if( pid not in productIDs ) : productIDs.append( pid )
            if( pid not in productIDs1 ) : productIDs1[pid] = []
            productIDs1[pid].append( product )

        aliases = other.rootAncestor.PoPs.aliases
        productIDs2 = {}
        for product in other.products : 
            pid = product.id
            if( pid in aliases ) : pid = aliases[pid].pid             # Handle old GNDS files which used meta-stable ids instead of nuclear level ids.
            if( pid not in productIDs ) : productIDs.append( pid )
            if( pid not in productIDs2 ) : productIDs2[pid] = []
            productIDs2[pid].append( product )

        for productID in productIDs :
            if( productID not in productIDs1 ) :        # This should only happen for photons and for reaction 'sumOfRemainingOutputChannels'.
                diffResults.append( 'Product missing from channel - 1', 'product id = "%s"' % productID, self.toXLink( ), other.toXLink( ) )
            elif( productID not in productIDs2 ) :      # This should only happen for photons and for reaction 'sumOfRemainingOutputChannels'.
                diffResults.append( 'Product missing from channel - 2', 'product id = "%s"' % productID, self.toXLink( ), other.toXLink( ) )
            else :
                if( len( productIDs1[productID] ) != len( productIDs2[productID] ) ) :
                    diffResults.append( 'More than one product with same id', 'product id = "%s"' % productID, self.toXLink( ), other.toXLink( ) )
                else :
                    for i1, product1 in enumerate( productIDs1[productID] ) :
                        product2 = productIDs2[productID][i1]
                        product1.diff( product2, diffResults )

    def getFinalProductList( self ) :
        """Returns a list of [ multiplicity, name ] in order for the final products of this channel."""

        mP = []
        for product in self :
            if( product.outputChannel is None ) :
                try :
                    multiplicity = product.multiplicity.getConstant( )
                except :
                    multiplicity = "(?)"                        # ????? This needs work.
                mP += [ [ multiplicity, product.id ] ]
            else :
                mP += product.outputChannel.getFinalProductList( )
        return( mP )

    def getProductsWithName( self, name ) :
        """Returns a list of all the channel's products with given name ('gamma','n', etc)."""

        return [ prod for prod in self if prod.id == name ]

    def getProductWithName( self, name ) :
        """
        Return the product with given name. If no such product is found, or if more than one are found, raise an Exception.
        """

        prods = self.getProductsWithName( name )
        if( len( prods ) != 1 ) : raise Exception( "Unique product named '%s' not found in %s" % ( name, self ) )
        return( prods[0] )

    def cullStyles( self, styleList ) :

        self.__Q.cullStyles( styleList )
        self.__products.cullStyles( styleList )
        self.__fissionFragmentData.cullStyles( styleList )

    def reactionProducts( self, _reactionProducts ) :

        for product in self.products : product.reactionProducts( _reactionProducts )
        if( self.process is not None ) : _reactionProducts[self.process] = reactionProductsModule.ReactionProduct( reactionProductsModule.Category.process, 0 )

        return( _reactionProducts )

    def removeStyles( self, styleLabels ) :

        self.__Q.removeStyles( styleLabels )
        self.__products.removeStyles( styleLabels )
        self.__fissionFragmentData.removeStyles( styleLabels )

    def thresholdQAs( self, unit, final = True ) :

        if( len( self.__Q ) > 0 ) :
            Q = self.__Q.thresholdQAs( unit )
        else :                      # Calculate from particle masses.
            massUnit = unit + '/c**2'
            reactionSuite = self.getRootAncestor( )
            projectile = reactionSuite.PoPs[reactionSuite.projectile]
            targetID = reactionSuite.target
            if( targetID in reactionSuite.PoPs.aliases ) : targetID = reactionSuite.PoPs[targetID].pid
            target = reactionSuite.PoPs[targetID]
            Q = target.mass[0].float( massUnit ) + projectile.mass[0].float( massUnit )
            for product in self:
                try :   
                    Q -= product.getMass(massUnit) * product.multiplicity.getConstant()
                except :
                    if( product.id == IDsPoPsModule.photon ) : continue
                    raise ValueError( "Non-constant Q-value must be explicitly listed in GNDS!" )
        if( final ) :
            for product in self.__products : Q += product.thresholdQAs( unit, final = final )
        return( Q )  

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        kwargs['outputChannel'] = self

        for productIndex, product in enumerate( self.__products ) :
            kwargs['productIndex'] = str( productIndex )
            product.calculateAverageProductData( style, indent = indent, **kwargs )

        self.__fissionFragmentData.calculateAverageProductData( style, indent, **kwargs )

    def partialProductionIntegral( self, reaction_suite, productID, energyIn, energyOut = None, muOut = None, phiOut = None, 
                frame = standardsModule.frames.labToken, LegendreOrder = 0, **kwargs ) :

        partialProductionIntegralSum = 0.0
        for product in self.__products :
            partialProductionIntegralSum += product.partialProductionIntegral( reaction_suite, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                    phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )

        return( partialProductionIntegralSum )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']

        for productIndex, product in enumerate( self ) :
            tempInfo['productIndex'] = str( productIndex )
            product.processMC_cdf( style, tempInfo, indent2 )

        self.__fissionFragmentData.processMC_cdf( style, tempInfo, indent2 )

    def processMultiGroup( self, style, tempInfo, indent ) :

        indent2 = indent + tempInfo['incrementalIndent']

        self.__Q.processMultiGroup( style, tempInfo, indent )

        for productIndex, product in enumerate( self ) :
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productName'] = product.id
            tempInfo['productLabel'] = product.label
            product.processMultiGroup( style, tempInfo, indent2 )

        self.__fissionFragmentData.processMultiGroup( style, tempInfo, indent2 )

    def QToPointwiseLinear( self, final = True, **kwargs ) :

        linearQ = self.__Q.toPointwise_withLinearXYs( **kwargs )
        if( final ) :
            for product in self.__products :
                if( product.outputChannel is not None ) : linearQ += product.outputChannel.Q.toPointwise_withLinearXYs( final = final, **kwargs )
        return( linearQ )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        qualifiers = ''
        if( self.process is not None ) : qualifiers += ' process="%s"' % str( self.process )
        xmlStringList = [ '%s<%s genre="%s"%s>' % ( indent, self.moniker, self.genre, qualifiers ) ]
        xmlStringList += self.__Q.toXMLList( indent2, **kwargs )
        xmlStringList += self.__products.toXMLList( indent2, **kwargs )
        xmlStringList += self.__fissionFragmentData.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def channelParticlesToString( self, prefix = "", suffix = "", MT = 0 ) :

        if( self.genre == Genre.sumOfRemainingOutputChannels ) : return( 'sumOfRemainingOutputChannels' )

        return( prefix + self.toString( simpleString = True ) + suffix )

    def toStandardString( self ) :
        """
        Returns a more standard nuclear, but less informative, string representation for the channel.
        For output channel it is list of "products)Residual" (e.g., "2n)Pu_238")."""

        mP = self.getFinalProductList( )
        s = ''
        for m, p in mP :
            if( m == 1 ) :
                sp = p
            else :
                sm = m
                if( type( m ) != type( '' ) ) : sm = str( m )
                sp = '%s%s' % ( sm, p )
            if( s == '' ) :
                s = sp
            else :
                s = s + ' + %s' % sp
        return( s )

    def toString( self, indent = '', simpleString = False, exposeGammaMultiplicity = False, MT = 0 ) :
        """Returns a string representation of self. If simpleString is True, the string contains only the initial
        products, and not any decay products."""

        if( self.genre == Genre.sumOfRemainingOutputChannels ) :
            if( self.process is None ) : return( 'sumOfRemainingOutputChannels' )

        s = ''
        p = ''
        gammaString = None
        for product in self.__products :
            if( product.id == IDsPoPsModule.photon ) :
                gammaStringP = product.toString( simpleString = True, exposeGammaMultiplicity = exposeGammaMultiplicity )
                if( gammaString is None ) : gammaString = gammaStringP
                if( 'energyDependent' in gammaStringP ) : gammaString = gammaStringP
            else :
                s += '%s%s' % ( p, product.toString( simpleString = simpleString, exposeGammaMultiplicity = exposeGammaMultiplicity ) )
                p = ' + '
        if( gammaString is not None ) : s += '%s%s' % ( p, gammaString )
        if( self.process is not None ) : s += " [%s]" % self.process
        return( indent + s )

    def parseXMLNode( self, element, xPath, linkData ) :

        from . import product as productModule
        from .channelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule

        xPath.append( element.tag )

        self.genre = element.get( 'genre' )
        self.process = element.get( 'process' )

        for child in element :
            if( child.tag == QModule.component.moniker ) :
                Q = QModule.parseXMLNode( child, xPath, linkData )
                for QForm in Q : self.__Q.add( QForm )
            elif( child.tag == productModule.products.moniker ) :
                self.__products.parseXMLNode( child, xPath, linkData )
            elif( child.tag == fissionFragmentDataModule.fissionFragmentData.moniker ) :
                self.fissionFragmentData.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( "Parsing %s not yet supported" % child.tag )

        xPath.pop( )

def parseXMLNode( element, xPath, linkData ) :
    """Translate '<outputChannel>' from xml."""

    xPath.append( element.tag )

    _outputChannel = outputChannel( Genre.NBody )

    xPath.pop( )

    _outputChannel.parseXMLNode( element, xPath, linkData )

    return( _outputChannel )
