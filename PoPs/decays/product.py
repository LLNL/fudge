# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the decay product classes.
"""

from .. import misc as miscModule
from .. import suite as suiteModule

class product( miscModule.classWithLabelKey ) :

    moniker = 'product'

    def __init__( self, label, pid ) :
        """
        :param label: unique label. Usually equal to the pid unless more than one is emitted
        :param pid: PoPs id for the emitted particle
        """

        miscModule.classWithLabelKey.__init__( self, label )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label not str' )
        self.__label = label

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid not str' )
        self.__pid = pid

    @property
    def pid( self ) :

        return( self.__pid )

    def convertUnits( self, unitMap ) :

        pass

    def copy( self ) :

        return( self.__class__( self.label, self.pid ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        XMLStringList = [ '%s<%s label="%s" pid="%s"/>' % ( indent, self.moniker, self.label, self.pid ) ]
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )
        self = cls( element.get('label'), element.get('pid') )
        xPath.pop( )

        self.parseXMLNode( element, xPath, linkData )
        return( self )

class suite( suiteModule.suite ) :

    moniker = 'products'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( product, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        if( element is None ) : return
        xPath.append( element.tag )

        for child in element :
            self.add( product.parseXMLNodeAsClass( child, xPath, linkData ) )

        xPath.pop( )
        return( self )
