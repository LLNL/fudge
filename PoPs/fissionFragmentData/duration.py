# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the duration and durations classes,
used to organize independent and cumulative fission product yields by the time since fission.

The duration consists of a time constant (i.e. elapsed time since fission) and a set of product yields.
It may also contain one or more Q-matrices, which can be used to transform fission product yields to
a later time constant.
"""

from .. import misc as miscModule
from .. import suite as suiteModule
from . import time as timeModule
from . import yields as yieldsModule
from . import QMatrix as QMatrixModule

class duration( miscModule.classWithLabelKey ) :

    moniker = 'duration'

    def __init__( self ) :

        miscModule.classWithLabelKey.__init__( self, '' )

        self.__time = timeModule.suite( )
        self.__time.setAncestor( self )

        self.__yields = yieldsModule.yields( )
        self.__yields.setAncestor( self )

        self.__QMatrix = QMatrixModule.suite( )
        self.__QMatrix.setAncestor( self )

    @property
    def time( self ) :

        return( self.__time )

    @property
    def yields( self ) :

        return( self.__yields )

    @property
    def QMatrix( self ):

        return( self.__QMatrix )

    def sortCompare( self, other ) :

        return( 1 )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.time.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.yields.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ):

        xPath.append(element.tag)
        for child in element:
            if child.tag == timeModule.suite.moniker:
                self.time.parseXMLNode(child, xPath, linkData)
            elif child.tag == yieldsModule.yields.moniker:
                self.yields.parseXMLNode(child, xPath, linkData)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return (self)

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        self = cls()
        xPath.pop()
        self.parseXMLNode(element, xPath, linkData)

        return self


class suite( suiteModule.sortedSuite ) :

    moniker = 'durations'

    def __init__( self ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( duration, ) )
