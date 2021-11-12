# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule

from fudge import suites as suitesModule

from . import duration as durationModule

class productYield( ancestryModule.ancestry ) :

    moniker = 'productYield'

    def __init__( self, label ) :

        ancestryModule.ancestry.__init__( self )

        self.__label = label

        self.__durations = durationModule.suite( )
        self.__durations .setAncestor( self )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def durations( self ) :

        return( self.__durations )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s"> ' % ( indent, self.moniker, self.__label ) ]
        XMLStringList += self.__durations.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        PY = cls(element.get('label'))
        for child in element:
            if child.tag == durationModule.suite.moniker:
                PY.durations.parseXMLNode(child, xPath, linkData)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return PY

class suite( suitesModule.suite ) :

    moniker = 'productYields'

    def __init__( self ) :

        suitesModule.suite.__init__( self, ( productYield, ) )
