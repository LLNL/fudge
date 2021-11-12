# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule

from . import multiplicity as multiplicityModule
from . import decayConstant as decayConstantModule

class delayedNeutronData( ancestryModule.ancestry ) :

    moniker = 'delayedNeutronData'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__multiplicity = multiplicityModule.suite( )
        self.__multiplicity.setAncestor( self )

        self.__decayConstants = decayConstantModule.suite( )
        self.__decayConstants.setAncestor( self )

    @property
    def multiplicity( self ) :

        return( self.__multiplicity )

    @property
    def decayConstants( self ) :

        return( self.__decayConstants )

    def convertUnits( self, unitMap ) :

        self.__multiplicity.convertUnits( unitMap )
        self.__decayConstants.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__multiplicity.toXMLList( indent2, **kwargs )
        XMLStringList += self.__decayConstants.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 1 ) : return( [] )
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ):

        xPath.append(element.tag)
        for child in element:
            if child.tag == multiplicityModule.suite.moniker:
                self.multiplicity.parseXMLNode(child, xPath, linkData)
            elif child.tag == decayConstantModule.suite.moniker:
                self.decayConstants.parseXMLNode(child, xPath, linkData)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return self

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        self = cls()
        xPath.pop()
        self.parseXMLNode(element, xPath, linkData)

        return self
