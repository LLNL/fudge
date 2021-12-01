# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the 'atomic' class, used to store properties of atoms.
The primary role of atomic data in PoPs is to store excited and partially-ionized atomic configurations,
along with their binding energy and decay radiation.

Each chemicalElement in a PoPs database may contain an <atomic> section.
"""

from xData import ancestry as ancestryModule

from . import configuration as configurationModule

class atomic( ancestryModule.ancestry ) :

    moniker = 'atomic'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__configurations = configurationModule.configurations( )
        self.__configurations.setAncestor( self )

    @property
    def configurations( self ) :

        return( self.__configurations )

    def convertUnits( self, unitMap ) :
        """See documentation in PoPs.database.convertUnits"""

        self.__configurations.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _atomic = atomic( )
        for configuration in self__.configurations : _atomic.configurations.add( configuration.copy( ) )
        return( _atomic )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__configurations.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            if( child.tag == configurationModule.configurations.moniker ) :
                self.configurations.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( 'Invalid child = "%s" for %s' % ( child.tag, self.moniker ) )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( )
        self.parseXMLNode( element, xPath, linkData )

        xPath.pop( )
        return( self )
