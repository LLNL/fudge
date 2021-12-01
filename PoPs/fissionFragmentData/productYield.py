# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .. import misc as miscModule
from .. import suite as suiteModule
from . import nuclides as nuclidesModule
from . import duration as durationModule

class productYield( miscModule.classWithLabelKey ) :

    moniker = 'productYield'

    def __init__( self, label, nuclides ) :

        miscModule.classWithLabelKey.__init__( self, '' )

        self.__nuclides = nuclidesModule.nuclides( nuclides )
        self.__nuclides.setAncestor( self )

        self.__durations = durationModule.suite( )
        self.__durations .setAncestor( self )

    @property
    def nuclides( self ) :

        return( self.__nuclides )

    @property
    def durations( self ) :

        return( self.__durations )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s> ' % ( indent, self.moniker ) ]
        XMLStringList += self.__nuclides.toXMLList( indent2, **kwargs )
        XMLStringList += self.__durations.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ):

        xPath.append(element.tag)
        for child in element:
            if child.tag == nuclidesModule.nuclides.moniker:
                self.nuclides.parseXMLNode(child, xPath, linkData)
            elif child.tag == durationModule.suite.moniker:
                self.durations.parseXMLNode(child, xPath, linkData)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return self

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        self = cls(label='', nuclides=[])
        xPath.pop()
        self.parseXMLNode(element, xPath, linkData)

        return self

class suite( suiteModule.suite ) :

    moniker = 'productYields'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( productYield, ) )
