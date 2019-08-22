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

This intereaction represent a neutron hitting an Oxygen-16 target producing a n and Oxygen-16 in the 6th excited. The
excited Oxygen-16 subsequently decays into He4 and C12. The output channel for this interaction is "n + O16_e6" and
the "O16_e6" contains the decay channel "He4 + C12".

Output and decay channels can be further divided into twoBody or NBody channels.

    :twoBody:       The channel initially has two products whose angular and energy correlates are known. Only the 
                    center-of-mass angular data for one of the products is needed. The products' outgoing energies can 
                    be calculated from their center-of-mass angular information and other data like the masses of the 
                    particles involved in a reaction or decay.

    :NBody:         The channel has two or more products whose angular and energy correlates are not known. The full double
                    differential data for each product must be given.

This module contains the following channel classes::

    channel     : the base class for all other classes.
        outputChannel                   : the base class for the following outout channels,
            twoBodyOutputChannel            : twoBody output channel.
            NBodyOutputChannel              : NBody output channel.
        decayChannel                    : the base class for the following decay channels,
            twoBodyDecayChannel             : twoBody decay channel.
            NBodyDecayChannel               : NBody decay channel.
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
decayChannelToken = "decayChannel"

fissionGenreTotal = 'total'
fissionGenreFirstChance = 'firstChance'
fissionGenreSecondChance = 'secondChance'
fissionGenreThirdChance = 'thirdChance'
fissionGenreFourthChance = 'fourthChance'

manyToken = 'many'

from fudge.core.utilities import fudgeExceptions, brb
from fudge.core.ancestry import ancestry
import channelData
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty

class channel( ancestry ) :
    """This is the class for a gnd channel which is a list of particles (i.e., product class objects)."""

    moniker = outputChannelToken

    def __init__( self, genre, Q = None, numberOfAllowedParticles = manyToken, fissionGenre = None ) :
        """Creates a new channel object."""

        ancestry.__init__( self, self.moniker, None )
        self.genre = genre
        self.numberOfAllowedParticles = numberOfAllowedParticles
        self.fissionGenre = fissionGenre
        self.fissionEnergyReleased = None
        self.setQ( Q )
        self.particles = []
        self._labels = []

    def __cmp__( self, other ) :
        """Compares self to other."""

        if( other == None ) : return( -1 )
        n1 = len( self )
        n2 = len( other )
        n = min( n1, n2 )
        for i in xrange( n ) :
            if( self[i] < other[i] ) : return( -1 )
            if( self[i] > other[i] ) : return(  1 )
        if( n1 < n2 ) : return( -1 )
        if( n1 > n2 ) : return(  1 )
        return( 0 )

    def __getitem__( self, i ) :
        """Returns the i^th product in the channel."""

        return( self.particles[i] )

    def __len__( self ) :
        """Returns the number of products for this channel."""

        return( len( self.particles ) )

    def  __str__( self ) :
        """Converts channel object to a string representation."""

        return( self.toString( simpleString = False ) )

    def addProduct( self, product_ ) :
        """Adds particle to the end of this channel's particle list. Also generates a unique label if not supplied."""

        import product

        if( not( isinstance( product_, product.product ) ) ) :
            raise Exception( 'Only products (not %s) can be added to a channel!' % type( product_ ) )
        if( ( self.numberOfAllowedParticles != manyToken ) and ( len( self.particles ) >= self.numberOfAllowedParticles ) ) :
            raise fudgeExceptions.FUDGE_Exception( 'Adding more particles than channel allows (=%s)' % self.numberOfAllowedParticles )
        if not product_.getLabel() or product_.getLabel() in self._labels:
            name = product_.getName()
            index = len([prod for prod in self if prod.getName() == name])
            if index>0:
                import string
                suffixCharacters = string.ascii_lowercase
                suffix = ''
                while index>0:
                    index, index1 = divmod( index-1, 26 )
                    suffix = suffixCharacters[index1] + suffix
                product_.label = name+'__'+suffix
            else:
                product_.label = name
        product_.setParent( self )
        self.particles.append( product_ )
        self._labels.append( product_.getLabel() )

    def addFissionEnergyReleased( self, fissionEnergyReleased ) :

        self.fissionEnergyReleased = fissionEnergyReleased

    def removeProductAtIndex( self, index ) :

        del self.particles[index]

    def checkProductFrame( self ) :
        """Calls checkProductFrame for self's products."""

        for product in self : product.checkProductFrame( )

    def getFinalProductList( self ) :
        """Returns a list of [ multiplicity, name ] in order for the final products of this channel."""

        mP = []
        for particle in self :
            if( particle.decayChannel == None ) :
                try :
                    multiplicity = particle.multiplicity.getConstant( )
                except :
                    multiplicity = "(?)"                        # ????? This needs work.
                mP += [ [ multiplicity, particle.particle.getToken( )] ]
            else :
                mP += particle.decayChannel.getFinalProductList( )
        return( mP )

    def getProductsWithName( self, name ) :
        """Returns a list of all the channel's products with given name ('gamma','n', etc)."""

        return [prod for prod in self if prod.getName() == name]

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
            Q = self.getRootParent().target.getMass(massUnit) + self.getRootParent().projectile.getMass(massUnit)
            for product in self:
                try: Q -= product.getMass(massUnit) * product.multiplicity.getConstant()
                except:
                    if product.getName() == 'gamma': continue
                    raise ValueError( "Non-constant Q-value must be explicitly listed in GND!" )
        # Q????? need recursion for final = True
        return( Q )  

    def isFission( self ) :

        return( not( self.fissionGenre is None ) )

    def setQ( self, Q ) :
        """
        Sets channel's Q. Q can be None, a string, a PhysicalQuantityWithUncertainty, a Q form or a Q component. 
        Both string and PhysicalQuantityWithUncertainty must have units of energy.
        """

        if( Q is None ) :
            self.Q = Q
        else :
            Q_ = Q
            if( ( type( Q_ ) == str ) or ( isinstance( Q_, PhysicalQuantityWithUncertainty ) ) ) : Q_ = channelData.Q.constant( Q_ )
            if( not ( isinstance( Q_, channelData.Q.component ) ) ) : Q_ = channelData.Q.component( Q_ )
            self.Q = Q_
            self.Q.setParent( self )

    def toXMLList( self, flags, indent = '' ) :

        indent2 = indent + '  '
        Qstr = ''
        if self.Q is not None: Qstr = ' Q="%s"' % self.Q.getXMLAttribute()
        xmlString = [ '%s<%s genre="%s"%s>' % ( indent, self.moniker, self.genre, Qstr ) ]
        if( not( self.fissionEnergyReleased is None ) ) : xmlString += self.fissionEnergyReleased.toXMLList( indent = indent2 )
        if self.Q is not None: xmlString += self.Q.toXMLList( indent = indent2 )
        for particle in self : xmlString += particle.toXMLList( flags, indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def setFissionGenre( self ) :
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
        particles, and not any decay particles."""

        if( self.genre in [ sumOfRemainingOutputChannelsGenre ] ) : return( self.genre )
        s, p, gammaString, delayNeutrons = '', '', None, 0
        for particle in self.particles :
            addParticle = True
            if( self.isFission( ) and ( particle.getName( ) == 'n' ) ) :
                if( particle.getAttribute( 'emissionMode' ) == 'delayed' ) :
                    addParticle = False
                    delayNeutrons += 1
            if( addParticle ) :
                if( particle.getName( ) == 'gamma' ) :
                    gammaStringP = particle.toString( simpleString = True, exposeGammaMultiplicity = exposeGammaMultiplicity )
                    if( gammaString is None ) : gammaString = gammaStringP
                    if( 'energyDependent' in gammaStringP ) : gammaString = gammaStringP
                else :
                    s += '%s%s' % ( p, particle.toString( simpleString = simpleString, exposeGammaMultiplicity = exposeGammaMultiplicity ) )
                    p = ' + '
        if( delayNeutrons > 0 ) :
            s += '%s%s' % ( p, "n[emissionMode:'%s delayed']" % delayNeutrons )
            p = ' + '
        if( not( gammaString is None ) ) : s += '%s%s' % ( p, gammaString )
        process = []
        if( self.isFission( ) ) :
            process.append( "%s fission" % self.fissionGenre )
        if( len( process ) > 0 ) : s += " [%s]" % ":".join( process )
        return( s )

class outputChannel( channel ) :

    def __init__( self, genre, Q = None, numberOfAllowedParticles = manyToken ) :

        channel.__init__( self, genre, Q, numberOfAllowedParticles = numberOfAllowedParticles )

class twoBodyOutputChannel( outputChannel ) :
    """This is a two-body class version of channel. In particular, two-body is an output channel where only the
    angular information is needed to calculate the outgoing angle/energy for the two-bodies. Either of the bodies may 
    subsequently decay/breakup."""

    moniker = outputChannelToken

    def __init__( self, Q = None ) :
        """This is the constructor for the twoBodyOutputChannel class. See class channel for more information."""

        outputChannel.__init__( self, twoBodyGenre, Q, numberOfAllowedParticles = 2 )

class NBodyOutputChannel( outputChannel ) :
    """This is a N-body class version of channel. It is for all output channels that do not fit under the 
    twoBodyOutputChannel class. In particular, the going energy/angular information for each particle is 
    not correlated with the other particles. Any of the bodies may decay/breakup."""

    def __init__( self, Q = None ) :
        """This is the constructor for the NBodyOutputChannel class. See this class's documentation and
        the class channel for more information."""

        outputChannel.__init__( self, NBodyGenre, Q )

class decayChannel( channel ) :

    def __init__( self, genre, Q = None, numberOfAllowedParticles = manyToken ) :

        channel.__init__( self, genre, Q, numberOfAllowedParticles = numberOfAllowedParticles )

class twoBodyDecayChannel( decayChannel ) :
    """A decay channel producing exactly two products: gamma or alpha decays."""

    moniker = decayChannelToken

    def __init__( self, Q = None ) :

        decayChannel.__init__( self, twoBodyGenre, Q, numberOfAllowedParticles = 2 )

class NBodyDecayChannel( decayChannel ) :
    """A decay channel producing more than two products: beta decay, spontaneous fission """

    moniker = decayChannelToken

    def __init__( self, Q = None ) :

        decayChannel.__init__( self, NBodyGenre, Q )

class fissionChannel( channel ) :

    def __init__( self, Q, fissionGenre ) :

        channel.__init__( self, NBodyGenre, Q, numberOfAllowedParticles = manyToken, fissionGenre = fissionGenre )

class sumOfRemainingOutputChannels( channel ) :
    """This channel is a catchall for all output channels not described for a given input channel."""

    def __init__( self, Q = None ) :
        """This is the constructor for the sumOfRemainingOutputChannels class. See this class's documentation and
        the class channel for more information."""

        channel.__init__( self, sumOfRemainingOutputChannelsGenre, Q )

class productionChannel( channel ) :
    """
    This is a production channel version of the channel class. A production channel is one that contains
    information about the production of outgoing particles for various outgoing channels (that is, independent 
    of the actual outgoing channel).
    """

    def __init__( self, Q = None ) :
        """
        This is the constructor for the productionChannel class. See this class's documentation and
        the class channel for more information.
        """

        if( Q is None ) : Q = channelData.Q.component( channelData.Q.notApplicable( ) )
        channel.__init__( self, productionGenre, Q, numberOfAllowedParticles = manyToken )

def parseXMLNode( channelElement, xPath=[], linkData={} ):
    """Translate '<outputChannel>' or '<decayChannel>' from xml. """
    xPath.append( channelElement.tag )
    from fudge import gnd
    Q_str = channelElement.get( channelData.base.QToken )
    if Q_str is not None:
        try:
            Q = channelData.Q.constant( PhysicalQuantityWithUncertainty( Q_str ) )
        except:
            if Q_str == 'N/A': Q = channelData.Q.notApplicable()
            else: raise ValueError("Can't extract a Q-value from string '%s'!" % Q_str)
        Q = channelData.Q.component( Q )
    else: Q = None
    outputChannelClass = {
            outputChannelToken: { 
                twoBodyGenre: twoBodyOutputChannel, NBodyGenre: NBodyOutputChannel,
                sumOfRemainingOutputChannelsGenre: sumOfRemainingOutputChannels,
                productionGenre: productionChannel },
            decayChannelToken: {
                twoBodyGenre: twoBodyDecayChannel, NBodyGenre: NBodyDecayChannel }
            } [ channelElement.tag ] [ channelElement.get('genre') ]
    outputChannel = outputChannelClass( Q )

    for dat in channelElement:
        if dat.tag==channelData.base.fissionEnergyReleasedToken:
            outputChannel.addFissionEnergyReleased( 
                    gnd.channelData.fissionEnergyReleased.component.parseXMLNode( dat, xPath, linkData ) )
        elif dat.tag=='product':
            outputChannel.addProduct( gnd.product.parseXMLNode( dat, xPath, linkData ) )
        else: raise ValueError("Parsing %s not yet supported" % dat.tag)
    xPath.pop()
    return outputChannel
