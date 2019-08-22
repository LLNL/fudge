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
In GND, a channel class contains all outgoing products from either a reaction (also called an interaction), 
from a decay or from a collection of interactions. Since products may subsequently decay, a decay channel may 
be nested within a product of another channel. Channels have the following genre:

    :output:        The output channel resulting from an interaction. All products, except gammas,
                    must be present. The ZA of the output channel must match the input channel's ZA.

    :decay:         The output channel resulting from the decay of a particle. All decay products, except gammas,
                    must be present. The ZA of the output channel must match the parent's ZA.

    :fission:     A pseudo output channel resulting from fission. Only prompt neutrons must be present.

    :sumOfRemainingOutputChannels: A pseudo output channel giving product data for all output channels not explicitly given.

    :production:  A pseudo output channel representing a product not fully covered in other output channels.

For example, consider the following nuclear interaction:

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

This module contains the following channel classes::

    channel     : the base class for all other classes.
        outputChannel                   : the base class for the following outout channels,
            twoBodyOutputChannel            : twoBody output channel.
            NBodyOutputChannel              : NBody output channel.
        fissionChannel                  :
        sumOfRemainingOutputChannels    :
        productionChannel               :
"""

__metaclass__ = type

twoBodyGenre = 'twoBody'
NBodyGenre = 'NBody'
productionGenre = 'production'
sumOfRemainingOutputChannelsGenre = 'sumOfRemainingOutputChannels'

outputChannelToken = "outputChannel"

fissionGenreTotal = 'total'
fissionGenreFirstChance = 'firstChance'
fissionGenreSecondChance = 'secondChance'
fissionGenreThirdChance = 'thirdChance'
fissionGenreFourthChance = 'fourthChance'

manyToken = 'many'

from PoPs import IDs as IDsPoPsModule

import xData.ancestry as ancestryModule

from . import tokens as tokensModule

from .channelData import Q as QModule

class channel( ancestryModule.ancestry ) :
    """This is the class for a gnd channel which is a list of products (i.e., product class objects)."""

    moniker = outputChannelToken
    ancestryMembers = ( 'Q', 'products' )

    def __init__( self, genre, numberOfAllowedParticles = manyToken, fissionGenre = None, process = None ) :
        """Creates a new channel object."""

        from . import product as productModule

        ancestryModule.ancestry.__init__( self )
        self.genre = genre
        self.numberOfAllowedParticles = numberOfAllowedParticles
        self.fissionGenre = fissionGenre
        self.process = process

        self.Q = QModule.component( )
        self.Q.setAncestor( self )

        self.products = productModule.products( )
        self.products.setAncestor( self )

    def __cmp__( self, other ) :
        """Compares self to other."""

        if( other is None ) : return( -1 )
        n1 = len( self )
        n2 = len( other )
        n = min( n1, n2 )
        for i1 in xrange( n ) :
            if( self[i1] < other[i1] ) : return( -1 )
            if( self[i1] > other[i1] ) : return(  1 )
        if( n1 < n2 ) : return( -1 )
        if( n1 > n2 ) : return(  1 )
        return( 0 )

    def __getitem__( self, i ) :
        """Returns the i^th product in the channel."""

        return( self.products[i] )

    def __len__( self ) :
        """Returns the number of products for this channel."""

        return( len( self.products ) )

    def  __str__( self ) :
        """Converts channel object to a string representation."""

        return( self.toString( simpleString = False ) )

    def removeProductAtIndex( self, index ) :

        del self.products[index]

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.Q.convertUnits( unitMap )
        self.products.convertUnits( unitMap )

    def checkProductFrame( self ) :
        """Calls checkProductFrame for self's products."""

        for product in self : product.checkProductFrame( )

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
        """Return the product with given name. If no such product is found, or if more than one are found,
        raise an Exception."""

        prods = self.getProductsWithName( name )
        if len(prods)==1: return prods[0]
        else: raise Exception("Unique product named '%s' not found in %s" % (name, self))

    def getFissionGenre( self ) :
        """Returns the channel's fission genre."""

        return( self.fissionGenre )

    def getGenre( self ) :
        """Returns the genre for this channel."""

        return( self.genre )

    def getConstantQAs( self, unit, final = True ) :

        if self.Q is not None:
            Q = self.Q.getConstantAs( unit )
        else:  # calculate from particle masses
            massUnit = unit + '/c**2'
            reactionSuite = self.getRootAncestor( )
            projectile = reactionSuite.PoPs[reactionSuite.projectile]
            target = reactionSuite.PoPs[reactionSuite.target]
            Q = target.mass[0].float( massUnit ) + projectile.mass[0].float( massUnit )
            for product in self:
                try: Q -= product.getMass(massUnit) * product.multiplicity.getConstant()
                except:
                    if( product.id == IDsPoPsModule.photon ) : continue
                    raise ValueError( "Non-constant Q-value must be explicitly listed in GND!" )
        # Q????? need recursion for final = True
        return( Q )  

    def isFission( self ) :

        return( not( self.fissionGenre is None ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        kwargs['outputChannel'] = self

        for productIndex, product in enumerate( self ) :
            kwargs['productIndex'] = str( productIndex )
            product.calculateAverageProductData( style, indent = indent, **kwargs )

    def processMC( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']

        for productIndex, product in enumerate( self ) :
            tempInfo['productIndex'] = str( productIndex )
            product.processMC( style, tempInfo, indent2 )

    def processMultiGroup( self, style, tempInfo, indent ) :

        indent2 = indent + tempInfo['incrementalIndent']

        self.Q.processMultiGroup( style, tempInfo, indent )

        tempInfo['transferMatrixComment'] = tempInfo['reactionSuite'].inputParticlesToReactionString( suffix = " --> " ) +  \
                self.toString( simpleString = True )
# BRBBRB
        for productIndex, product in enumerate( self ) :
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productName'] = product.id
            tempInfo['productLabel'] = product.label
            product.processMultiGroup( style, tempInfo, indent2 )

    def QToPointwiseLinear( self, final = True, **kwargs ) :

        linearQ = self.Q.toPointwise_withLinearXYs( **kwargs )
        if( final ) :
            for product in self.products :
                if( product.outputChannel is not None ) : linearQ += product.outputChannel.Q.toPointwise_withLinearXYs( final = final, **kwargs )
        return( linearQ )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        qualifiers = ''
        if( self.process is not None ) : qualifiers += ' process="%s"' % str( self.process )
        xmlStringList = [ '%s<%s genre="%s"%s>' % ( indent, self.moniker, self.genre, qualifiers ) ]
        xmlStringList += self.Q.toXMLList( indent2, **kwargs )
        xmlStringList += self.products.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def setFissionGenre( self, fissionGenre ) :
        """Sets the channel's fission genre."""

        self.fissionGenre = fissionGenre

    def channelParticlesToString( self, prefix = "", suffix = "", specialChannelsReturnGenre = True ) :

        if( specialChannelsReturnGenre and ( self.genre in [ sumOfRemainingOutputChannelsGenre ] ) ) : return( self.genre )
        return( prefix + self.toString( simpleString = True ) + suffix )

    def toStandardString( self ) :
        """Returns a more standard nuclear, but less informative, string representation for the channel.
        For output channel it is "list of products)Residual" (e.g., "2n)Pu_238")."""

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

    def toString( self, simpleString = False, exposeGammaMultiplicity = False ) :
        """Returns a string representation of self. If simpleString is True, the string contains only the initial
        products, and not any decay products."""

        if( self.genre in [ sumOfRemainingOutputChannelsGenre ] ) : return( self.genre )
        s, p, gammaString, delayNeutrons = '', '', None, 0
        for product in self.products :
            addParticle = True
            if( self.isFission( ) and ( product.id == IDsPoPsModule.neutron ) ) :
# BRB6 'emissionMode'?
                if( product.getAttribute( 'emissionMode' ) == tokensModule.delayedToken ) :
                    addParticle = False
                    delayNeutrons += 1
            if( addParticle ) :
                if( product.id == IDsPoPsModule.photon ) :
                    gammaStringP = product.toString( simpleString = True, exposeGammaMultiplicity = exposeGammaMultiplicity )
                    if( gammaString is None ) : gammaString = gammaStringP
# BRB6 'energyDependent'?
                    if( 'energyDependent' in gammaStringP ) : gammaString = gammaStringP
                else :
                    s += '%s%s' % ( p, product.toString( simpleString = simpleString, exposeGammaMultiplicity = exposeGammaMultiplicity ) )
                    p = ' + '
        if( delayNeutrons > 0 ) :
# BRB6 'n[emissionMode:'%s delayed']'?
            s += '%s%s' % ( p, "n[emissionMode:'%s delayed']" % delayNeutrons )
            p = ' + '
        if( not( gammaString is None ) ) : s += '%s%s' % ( p, gammaString )
        process = []
        if( self.isFission( ) ) :
# BRB6 '%s fission'?
            process.append( "%s fission" % self.fissionGenre )
        elif( self.process is not None ) :
            process.append( str( self.process ) )
        if( len( process ) > 0 ) : s += " [%s]" % ":".join( process )
        return( s )

class outputChannel( channel ) :

    def __init__( self, genre, numberOfAllowedParticles = manyToken, process = None ) :

        channel.__init__( self, genre, numberOfAllowedParticles = numberOfAllowedParticles, process = process )

class twoBodyOutputChannel( outputChannel ) :
    """This is a two-body class version of channel. In particular, two-body is an output channel where only the
    angular information is needed to calculate the outgoing angle/energy for the two-bodies. Either of the bodies may 
    subsequently decay/breakup."""

    moniker = outputChannelToken

    def __init__( self, process = None ) :
        """This is the constructor for the twoBodyOutputChannel class. See class channel for more information."""

        outputChannel.__init__( self, twoBodyGenre, numberOfAllowedParticles = 2, process = process )

class NBodyOutputChannel( outputChannel ) :
    """This is a N-body class version of channel. It is for all output channels that do not fit under the 
    twoBodyOutputChannel class. In particular, the going energy/angular information for each particle is 
    not correlated with the other products. Any of the bodies may decay/breakup."""

    def __init__( self, process = None ) :
        """This is the constructor for the NBodyOutputChannel class. See this class's documentation and
        the class channel for more information."""

        outputChannel.__init__( self, NBodyGenre, process = process )

class fissionChannel( channel ) :

    def __init__( self, fissionGenre ) :

        channel.__init__( self, NBodyGenre, numberOfAllowedParticles = manyToken, fissionGenre = fissionGenre )

class sumOfRemainingOutputChannels( channel ) :
    """This channel is a catchall for all output channels not described for a given input channel."""

    def __init__( self ) :
        """This is the constructor for the sumOfRemainingOutputChannels class. See this class's documentation and
        the class channel for more information."""

        channel.__init__( self, sumOfRemainingOutputChannelsGenre )

class productionChannel( channel ) :
    """
    This is a production channel version of the channel class. A production channel is one that contains
    information about the production of outgoing products for various outgoing channels (that is, independent 
    of the actual outgoing channel).
    """

    def __init__( self ) :
        """
        This is the constructor for the productionChannel class. See this class's documentation and
        the class channel for more information.
        """

        channel.__init__( self, productionGenre, numberOfAllowedParticles = manyToken )

def parseXMLNode( channelElement, xPath, linkData ) :
    """Translate '<outputChannel>' from xml."""

    xPath.append( channelElement.tag )
    outputChannelClass = {
            outputChannelToken : {
                twoBodyGenre : twoBodyOutputChannel, NBodyGenre : NBodyOutputChannel,
                sumOfRemainingOutputChannelsGenre : sumOfRemainingOutputChannels,
                productionGenre : productionChannel }
            }[channelElement.tag][channelElement.get( 'genre' )]
    outputChannel = outputChannelClass( )
    if channelElement.get('process') is not None:
        outputChannel.process = channelElement.get('process')

# BRB6 'Q' and 'products' need to come from class monikers.
    for dat in channelElement :
        if( dat.tag == 'Q' ) :
            Q = QModule.parseXMLNode( dat, xPath, linkData )
            for QForm in Q : outputChannel.Q.add( QForm )
        elif( dat.tag == 'products' ) :
            outputChannel.products.parseXMLNode( dat, xPath, linkData )
        else :
            raise ValueError( "Parsing %s not yet supported" % dat.tag )
    xPath.pop( )
    return( outputChannel )
