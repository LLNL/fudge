# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the class for the suite chemicalElement.
"""

from .. import misc as miscModule

from ..atomic import atomic as atomicDataModule

from . import misc as chemicalElementMiscModule
from . import isotope as isotopeModule

class chemicalElement( miscModule.classWithSymbolKey ) :
    """
    This class represents a chemical element (i.e., Oyxgen, Iron, Uranium). A chemical element
    stores all isotopes that contain the same number of protons (i.e., O16, O17, and O18). 
    The members of a chemicalElement are:

    .. table:: Chemical element members
       :widths: auto

       =========== =========== =========================================================================
       member      type        description
       =========== =========== =========================================================================
       symbol      string      This is the symbol for the chemical element (e.g., 'O' for oxygen).
       Z           int         Number of protons for this chemical element (e.g., 8 for 'O').
       name        string      The name of the chemical element (i.e., 'Oxygen')
       atomicData  atomicData  Container for atomic data, including ionized states and atomic relaxation
       isotopes    suite       The list of isotopes
       =========== =========== =========================================================================

    """

    moniker = 'chemicalElement'

    def __init__( self, symbol, Z, name ) :
        """
        :param symbol: Chemical symbol for the element, e.g. 'Fe' (string)
        :param Z: Atomic number, e.g. 26 (int)
        :param name: Full name of the element, e.g. 'Iron'
        """

        miscModule.classWithSymbolKey.__init__( self, symbol )

        base, anti = miscModule.baseAntiFromID( symbol )

        self.__anti = anti == miscModule.antiSuffix

        chemicalElementMiscModule.checkZ( Z )
        self.__Z = Z
        if( chemicalElementMiscModule.symbolFromZ[Z] != base ) : raise ValueError( 'Z = %s and chemical element symbol = "%s" are not consistent' %
                ( Z, miscModule.toLimitedString( symbol ) ) )

        if( not( isinstance( name, str ) ) ) : TypeError( 'name must be a string' )
        if( chemicalElementMiscModule.nameFromZ[Z] != name ) : 
            raise ValueError( 'Z = "%s" and chemical element name = "%s" are not consistent, should be "%s".' %
                    ( Z, miscModule.toLimitedString( name ), chemicalElementMiscModule.nameFromZ[Z] ) )
        self.__name = name

        self.__atomicData = None

        self.__isotopes = isotopeModule.isotopes( )
        self.__isotopes.setAncestor( self ) 

    def __contains__( self, key ) : 

        return( key in self.__isotopes )

    def __getitem__( self, key ) :

        return( self.__isotopes[key] )

    def __iter__( self ) :
        """ Iterates over isotopes of this element (but not over the excited states of those isotopes) """

        for isotope in self.__isotopes : yield isotope

    @property
    def Z( self ) :

        return( self.__Z )

    @property
    def isAnti( self ) :

        return( self.__anti )

    @property
    def name( self ) :

        return( self.__name )

    @property
    def atomicData( self ) :

        return( self.__atomicData )

    @atomicData.setter
    def atomicData( self, value ) :

        if( not( isinstance( value, atomicDataModule.atomic ) ) ) : raise TypeError( 'value must be an instance of atomicData.' )
        self.__atomicData = value
        self.__atomicData.setAncestor( self )

    @property
    def isotopes( self ) :

        return( self.__isotopes )

    def add( self, item ) :
        """
        :param item: isotope instance, will be added to the list of isotopes
        """

        self.__isotopes.add( item )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        if( self.__atomicData is not None ) : self.__atomicData.convertUnits( unitMap )
        self.__isotopes.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _chemicalElement = chemicalElement( self.symbol, self.Z, self.name )
        if( self.__atomicData is not None ) : _chemicalElement.atomicData = self.__atomicData.copy( )
        for item in self.__isotopes : _chemicalElement.add( item.copy( ) )
        return( _chemicalElement )

    def sortCompare( self, other ) :

        if( not( isinstance( other, chemicalElement ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.Z - other.Z )

    def check(self, info):
        """ See check documentation in PoPs.database """

        from .. import warning as warningModule
        warnings = []

        for isotope in self:
            isotopeWarnings = isotope.check(info)
            if isotopeWarnings:
                warnings.append(warningModule.context('Isotope %s' % isotope.symbol, isotopeWarnings))

        return warnings

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s symbol="%s" Z="%s" name="%s">' % ( indent, self.moniker, self.symbol, self.Z, self.name ) ]
        if( self.__atomicData is not None ) : XMLStringList += self.__atomicData.toXMLList( indent2, **kwargs )
        XMLStringList += self.__isotopes.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        for child in element :
            if( child.tag == isotopeModule.isotopes.moniker ) :
                self.__isotopes.parseXMLNode( child, xPath, linkData )
            elif( child.tag == atomicDataModule.atomic.moniker ) :
                self.atomicData = atomicDataModule.atomic.parseXMLNodeAsClass( child, xPath, linkData )
            else :
                raise ValueError( 'Invalid child = "%s" for %s' % ( child.tag, self.moniker ) )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        self = cls( element.get('symbol'), int( element.get('Z') ), element.get('name') )
        xPath.pop( )

        self.parseXMLNode( element, xPath, linkData )

        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )

        for attributeName in element.keys() :
            if( attributeName == 'symbol' ) :
                symbol = element.get(attributeName)
            elif( attributeName == 'Z' ) :
                Z = int( element.get(attributeName) )
            elif( attributeName == 'name' ) :
                name = element.get(attributeName)
            else :
                raise ValueError( 'Unknown attribute = "%s"' % attributeName )

        return( cls( symbol, Z, name ).parseXMLNode( element, [], [] ) )
