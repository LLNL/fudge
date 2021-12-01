# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines classes for storing particle decay information.
Decay data are organized as follows::

   - particle
      - decayData
         - averageEnergies
         - decayModes
            - decayMode
               - probability
               - Q
               - spectra
               - decayPath
                  - decay

All decay info is stored inside the <decayData> section.

<averageEnergies> stores a summary of average outgoing energy for various types of decay radiation,
summed over all possible decay modes.

<decayModes> contains a list of <decayMode>, each with a corresponding probability.
A <decayMode> may also contain a <Q>  (i.e. the decay Q-value), a <spectra> (describing the outgoing energy spectrum
for each types of decay product), and a <decayPath> (listing each step of the decay along with outgoing products).

Ideally the outgoing spectra would be associated with specific decay products under the <decayPath>, but
often the spectra are summed over multiple products.  For example, in B-nn decay (beta-delayed 2-neutron emission),
the spectra of the 2 outgoing neutrons cannot be resolved, so the evaluation usually stores only the summed spectrum
for the two neutrons.
"""

import abc  # FIXME remove unused import?

from xData import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import misc as miscModule

from ..decays import probability as probabilityModule
from ..decays import Q as QModule
from ..decays import product as productModule
from ..decays import spectrum as spectrumModule
from ..decays import averageEnergy as averageEnergyModule

decayModesIT = "isomeric transition"
decayModesSF = "spontaneous fission"
decayModesParticle = ""

decayModeTypes = [decayModesIT, decayModesSF, decayModesParticle]

class decay( miscModule.classWithIndexKey ) :
    """
    Describes one step in a multi-step decay. For example, in beta-delayed 2-n emission
    the first <decay> is for the beta-decay, and the 2nd and 3rd <decay> elements describe the two neutrons
    being emitted sequentially.
    """

    moniker = 'decay'

    def __init__( self, index, type, complete = True ) :
        """
        :param index: int, identifies the position of this decay within the decayPath. 0 = first, 1 = second, etc.
        :param type: string identifying a decay type, such as "isomeric transition"
        :param complete: boolean, whether all outgoing products are listed
        """

        miscModule.classWithIndexKey.__init__( self, index )

        if( not( isinstance( type, str ) ) ) : raise TypeError( 'decay mode type not an instance of str' )
        if( type not in decayModeTypes ) : raise ValueError( 'decay mode type invalid' )
        self.__type = type

        if( not( isinstance( complete, bool ) ) ) : raise TypeError( 'complete not an instance of bool' )
        self.__complete = complete

        self.__products = productModule.suite( )
        self.__products.setAncestor( self )

    @property
    def complete( self ) :

        return( self.__complete )

    @property
    def products( self ) :

        return( self.__products )

    @property
    def type( self ) :

        return( self.__type )

    def convertUnits( self, unitMap ) :
        """See documentation in PoPs.database.convertUnits"""

        self.__products.convertUnits( unitMap )

    def copy( self ) :

        _decay = self.__class__( self.index, self.type, self.complete )
        for item in self.__products : _decay.products.add( item.copy( ) )
        return( _decay )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributes = ''
        if( self.type != '' ) : attributes = ' type="%s"' % self.type
        if( not( self.complete ) ) : attributes = ' complete="false"'
        XMLStringList = [ '%s<%s index="%s"%s>' % ( indent, self.moniker, self.index, attributes ) ]
        XMLStringList += self.products.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        self.products.parseXMLNode( element.find( self.products.moniker ), xPath, linkData )

        xPath.pop()
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        complete = True
        if element.get('complete') == 'false': complete = False
        self = cls( element.get('index'), element.get( 'type', '' ), complete )
        xPath.pop( )
        self.parseXMLNode( element, xPath, linkData )

        return( self )

class decayPath( suiteModule.suite ) :
    """
    Stores all <decay> steps involved in a decayMode
    """

    moniker = 'decayPath'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( decay, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            self.add( decay.parseXMLNodeAsClass( child, xPath, linkData ) )

        xPath.pop( )
        return( self )

class decayMode( miscModule.classWithLabelKey ) :
    """
    Describes one way a particle can decay. The description includes the probability, Q-value, outgoing products
    and spectra. For electro-magnetic decays, may also include a set of internal conversion coefficients
    and/or photon emission probabilities.
    """

    moniker = 'decayMode'

    def __init__( self, label, mode ) :
        """
        :param label: string, must be unique within this <decayModes> section
        :param mode: string, describes the decay mode. Examples: "beta+ or e.c.", "beta-,n", "alpha", etc.
        """

        miscModule.classWithLabelKey.__init__( self, label )

        if( not( isinstance( mode, str ) ) ) : raise TypeError( 'mode not an instance of str' )
        self.__mode = mode

        self.__probability = probabilityModule.suite( )
        self.__probability.setAncestor( self )

        self.__internalConversionCoefficients = spectrumModule.internalConversionCoefficients()
        self.__internalConversionCoefficients.setAncestor( self )

        self.__photonEmissionProbabilities = spectrumModule.photonEmissionProbabilities()
        self.__photonEmissionProbabilities.setAncestor( self )

        self.__Q = QModule.suite( )
        self.__Q.setAncestor( self )

        self.__decayPath = decayPath( )
        self.__decayPath.setAncestor( self )

        self.__spectra = spectrumModule.spectra( )
        self.__spectra.setAncestor( self )

    @property
    def mode( self ) :

        return( self.__mode )

    @property
    def probability( self ) :

        return( self.__probability )

    @property
    def internalConversionCoefficients( self ):

        return( self.__internalConversionCoefficients )

    @property
    def photonEmissionProbabilities( self ):

        return( self.__photonEmissionProbabilities )

    @property
    def decayPath( self ) :

        return( self.__decayPath )

    @property
    def Q( self ) :

        return( self.__Q )

    @property
    def spectra( self ) :

        return( self.__spectra )

    def convertUnits( self, unitMap ) :
        """See convertUnits documentation in PoPs.database"""

        self.__probability.convertUnits( unitMap )
        self.__Q.convertUnits( unitMap )
        self.__decayPath.convertUnits( unitMap )
        self.__spectra.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        other = self.__class__( self.label, self.mode )
        for item in self.__probability : other.probability.add( item.copy( ) )
        for item in self.__Q : other.Q.add( item.copy( ) )
        for item in self.__decayPath : other.decayPath.add( item.copy( ) )
        for item in self.__spectra : other.spectra.add( item.copy( ) )
        return( other )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s" mode="%s">' % ( indent, self.moniker, self.label, self.mode ) ]
        XMLStringList += self.probability.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.internalConversionCoefficients.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.photonEmissionProbabilities.toXMLList(indent=indent2, **kwargs)
        XMLStringList += self.Q.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.decayPath.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.spectra.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( "%s[@label='%s']" % (element.tag, element.get('label')) )

        for child in element :
            if( child.tag == probabilityModule.suite.moniker ) :
                self.probability.parseXMLNode( child, xPath, linkData )
            elif( child.tag == spectrumModule.internalConversionCoefficients.moniker ) :
                self.internalConversionCoefficients.parseXMLNode( child, xPath, linkData )
            elif( child.tag == spectrumModule.photonEmissionProbabilities.moniker ) :
                self.photonEmissionProbabilities.parseXMLNode( child, xPath, linkData )
            elif( child.tag == QModule.suite.moniker ) :
                self.Q.parseXMLNode( child, xPath, linkData )
            elif( child.tag == decayPath.moniker ) :
                self.decayPath.parseXMLNode( child, xPath, linkData )
            elif( child.tag == spectrumModule.spectra.moniker ) :
                self.spectra.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( 'Invalid tag = "%s"' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( element.get('label'), element.get('mode') )
        self.parseXMLNode( element, xPath, linkData )

        xPath.pop( )
        return( self )

class decayModes( suiteModule.suite ) :
    """
    Contains a list of decayMode instances
    """

    moniker = 'decayModes'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( decayMode, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            self.add( decayMode.parseXMLNodeAsClass( child, xPath, linkData ) )

        xPath.pop( )
        return( self )

class decayData( ancestryModule.ancestry ) :
    """
    Contains all decay information for a particle, including average energies for decay products
    and a list of decay modes.
    """

    moniker = 'decayData'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__decayModes = decayModes( )
        self.__decayModes.setAncestor( self )

        self.__averageEnergies = averageEnergyModule.averageEnergies( )
        self.__averageEnergies.setAncestor( self )

    @property
    def decayModes( self ) :

        return( self.__decayModes )

    @property
    def averageEnergies( self ) :

        return( self.__averageEnergies )

    def convertUnits( self, unitMap ) :
        """See convertUnits documentation in PoPs.database"""

        self.__decayModes.convertUnits( unitMap )
        self.__averageEnergies.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _decayData = decayData( )
        self.copyItems( _decayData )
        return( _decayData )

    def copyItems( self, other ) :
        """
        Copy all items in self to other
        :param other: decayData instance where contents of self will be copied
        """

        for item in self.__decayModes : other.decayModes.add( item.copy( ) )
        for item in self.__averageEnergies : other.averageEnergies.add( item.copy( ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList_subs  = self.__decayModes.toXMLList( indent = indent2, **kwargs )
        XMLStringList_subs += self.averageEnergies.toXMLList( indent = indent2, **kwargs )
        if( len( XMLStringList_subs ) == 0 ) : return( [] )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ] + XMLStringList_subs
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        children = { 'decayModes' : self.decayModes, 'averageEnergies' : self.averageEnergies }

        for child in element :
            if( child.tag in children ) :
                children[child.tag].parseXMLNode( child, xPath, linkData )
            else:
                if( not( self.parseExtraXMLElement( child, xPath, linkData ) ) ) :
                    raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )
