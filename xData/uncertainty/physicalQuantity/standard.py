# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the standard uncertainty class.
"""

from . import uncertainty as uncertaintyModule

class standard( uncertaintyModule.base ) :

    moniker = "standard"

    def __init__( self, value ) :

        uncertaintyModule.base.__init__( self )

        if( not( isinstance( value, uncertaintyModule.quantity ) ) ) : raise TypeError( 'Invalid quantity for %s uncertainty' % self.moniker )
        self.__value = value
        self.__value.setAncestor( self )

    @property
    def value( self ) :

        return( self.__value )

    def copy( self ) :

        return( self.__class__( self.__value.copy( ) ) )

    def parentConvertingUnits( self, factors ) :

        self.__value.parentConvertingUnits( factors )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.value.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += "</%s>" % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        child = element[0]
        if( child.tag == uncertaintyModule.double.moniker ) :
            uncertainty = uncertaintyModule.double.parseXMLNodeAsClass( child, xPath, linkData )
        else :
            raise Exception( 'Invalid child element = "%s".' % child.tag )

        _standard = standard( uncertainty )

        xPath.pop( )
        return( _standard )
